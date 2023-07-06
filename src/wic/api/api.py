"""CLT utilities."""
import shutil
import subprocess
from dataclasses import dataclass, field
from functools import singledispatch
from pathlib import Path
from typing import Any, List, Optional, TypeVar, Union

import cwl_utils.parser as cu_parser
import typeguard
import yaml
from cwl_utils.parser import CommandLineTool as CWLCommandLineTool
from cwl_utils.parser import DockerRequirement, load_document_by_uri
from pydantic import (  # pylint: disable=E0611
    BaseModel,
    Extra,
    Field,
    ValidationError,
    validator,
)
from pydantic.error_wrappers import ErrorWrapper  # pylint: disable=E0611

from wic.api._types import CWL_TYPES_DICT


class InvalidInputValue(Exception):
    pass


class MissingRequiredValue(Exception):
    pass


class InvalidStep(Exception):
    pass


class InvalidLink(Exception):
    pass


class InvalidCLT(Exception):
    pass


StrPath = TypeVar("StrPath", str, Path)
CWLInputParameter = Union[
    cu_parser.cwl_v1_0.CommandInputParameter,
    cu_parser.cwl_v1_1.CommandInputParameter,
    cu_parser.cwl_v1_2.CommandInputParameter,
]  # cwl_utils does not have CommandInputParameter union


def _is_docker_requirement(req: Any) -> bool:
    """Check if requirement is DockerRequirement."""
    try:
        typeguard.check_type(req, DockerRequirement)  # pylint: disable=E1120
    except BaseException:  # pylint: disable=broad-exception-caught
        return False
    return True


class Tool(BaseModel):  # pylint: disable=too-few-public-methods
    class Config:  # pylint: disable=too-few-public-methods:
        arbitrary_types_allowed = True
        extra = Extra.allow

    cwl_path: Path
    cwlVersion: str
    dockerContainer: str
    cwl: CWLCommandLineTool

    @validator("dockerContainer", pre=True, check_fields=False)
    # type: ignore
    def validate_docker(cls, docker):  # pylint: disable=no-self-argument
        """Validate that there is one docker requirement specified."""
        if not len(docker) > 0:
            raise InvalidCLT("no docker requirement specified")
        # schema salad already checks that len(docker) < 2
        docker = docker[0].dockerPull
        return docker


class CLTInput(BaseModel):  # pylint: disable=too-few-public-methods
    """Input of CLT."""

    inp_type: object
    name: str
    value: Any  # validation happens at assignment
    required: bool = True
    linked: bool = False

    def __init__(self, cwl_inp: CWLInputParameter) -> None:
        inp_type = cwl_inp.type
        if isinstance(inp_type, list) and "null" in inp_type:
            required = False
        else:
            required = True
        super().__init__(inp_type=inp_type, name=cwl_inp.id, required=required)

    @validator("name")
    # type: ignore
    def get_name_from_id(cls, cwl_id) -> Any:  # pylint: disable=no-self-argument
        """Return name of input from InputParameter.id."""
        return cwl_id.split("#")[-1]

    @validator("inp_type")
    # type: ignore
    def set_inp_type(cls, inp) -> object:  # pylint: disable=no-self-argument
        """Return inp_type."""
        if isinstance(inp, list):  # optional inps
            inp = inp[1]
        return CWL_TYPES_DICT[inp]

    def _set_value(
        self, __value: Any, check: bool = True, linked: bool = False
    ) -> None:
        """Set input value."""
        if check:
            if not isinstance(__value, self.inp_type):  # type: ignore
                raise TypeError(
                    f"invalid attribute type for {self.name}: got {__value.__class__.__name__}, expected {self.inp_type.__name__}"  # type: ignore # noqa: E266
                )
            self.value = __value
            return
        self.value = __value  # to be used for linking inputs */&
        if linked:
            self.linked = True


def _default_dict() -> dict:
    return {}


def _get_value_from_cfg(value: Any) -> Any:  # validation happens in Step
    if isinstance(value, dict):
        if "Directory" in value.values():
            try:
                value_ = Path(value["location"])
            except BaseException as exc:
                raise InvalidInputValue() from exc
            if not value_.is_dir():
                raise InvalidInputValue(f"{str(value_)} is not a directory")
            return value_
    else:
        return value


@singledispatch
def _path(path: Any) -> Path:
    raise TypeError


@_path.register
def _(path: str) -> Path:
    return Path(path)


@_path.register
def _(path: Path) -> Path:
    return path


def _is_link(s: str) -> bool:
    """Return True if s starts with & or *."""
    if s.startswith("&") or s.startswith("*"):
        return True
    return False


class Step(Tool):  # pylint: disable=too-few-public-methods
    """Base class for configured CLTs."""

    cwl_name: str
    inputs: List[CLTInput]
    yaml: dict
    cfg_yaml: dict = Field(default_factory=_default_dict)
    _input_names: list

    def __init__(self, cwl_path: StrPath, config_path: Optional[StrPath] = None):
        # validate using cwl.utils
        try:
            cwl = load_document_by_uri(cwl_path)
        except BaseException as s_e:
            e_w = ErrorWrapper(s_e, "invalid cwl file")
            raise ValidationError([e_w], Step)  # pylint: disable=raise-missing-from
        docker = [x for x in cwl.requirements if _is_docker_requirement(x)]
        with open(cwl_path, "r", encoding="utf-8") as file:
            yaml_file = yaml.safe_load(file)
        if config_path:
            with open(config_path, "r", encoding="utf-8") as file:
                cfg_yaml = yaml.safe_load(file)
        else:
            cfg_yaml = _default_dict()  # redundant, to avoid it being unbound
        cwl_name = _path(cwl_path).stem
        input_names = [inp.id.split("#")[-1] for inp in cwl.inputs]
        super().__init__(
            cwl_path=cwl_path,
            cwlVersion=cwl.cwlVersion,
            dockerContainer=docker,
            cwl=cwl,
            inputs=cwl.inputs,
            yaml=yaml_file,
            cfg_yaml=cfg_yaml,
            cwl_name=cwl_name,
            _input_names=input_names,
        )
        if config_path:
            self._set_from_io_cfg()

    @validator("inputs", pre=True)
    # type: ignore
    def cast_to_clt_input_model(
        cls, cwl_inps: List[CWLInputParameter]
    ):  # pylint: disable=no-self-argument
        """Populate inputs from cwl.inputs."""
        return [CLTInput(x) for x in cwl_inps]

    def __setattr__(self, __name: str, __value: Any) -> Any:  # pylint: disable=R1710
        if __name in self._input_names:
            index = self._input_names.index(__name)
            if isinstance(__value, CLTInput):
                if not __value.linked:
                    try:
                        local_input = self.inputs[index]
                        if not local_input.inp_type == __value.inp_type:
                            raise InvalidLink(
                                f"links must have the same input type. cannot link {local_input.name} to {__value.name}"
                            )
                        tmp = f"{__name}{self.cwl_name}"
                        local_input._set_value(f"*{tmp}", check=False, linked=True)
                        __value._set_value(f"&{tmp}", check=False, linked=True)
                    except BaseException as exc:
                        raise exc
                else:  # value is already linked to another inp
                    try:
                        local_input = self.inputs[index]
                        if not local_input.inp_type == __value.inp_type:
                            raise InvalidLink(
                                f"links must have the same input type. cannot link {local_input.name} to {__value.name}"
                            )
                        # tmp = f"{__name}{self.cwl_name}"
                        current_value = __value.value
                        local_input._set_value(
                            current_value.replace("&", "*"), check=False, linked=True
                        )
                    except BaseException as exc:
                        raise exc

            else:
                if isinstance(__value, str) and _is_link(__value):
                    self.inputs[index]._set_value(__value, check=False)
                else:
                    self.inputs[index]._set_value(__value)
        else:
            return super().__setattr__(__name, __value)

    def __getattribute__(self, __name: str) -> Any:
        if __name != "_input_names" and hasattr(self, "_input_names"):
            if __name in self._input_names:
                index = self._input_names.index(__name)
                return self.inputs[index]
        return super().__getattribute__(__name)

    def _set_from_io_cfg(self) -> None:
        for name, value in self.cfg_yaml.items():
            value_ = _get_value_from_cfg(value)
            setattr(self, name, value_)

    def _validate(self) -> None:
        for inp in self.inputs:
            if inp.required and inp.value is None:
                raise MissingRequiredValue(f"{inp.name} is required")

    @property
    def _yml(self) -> dict:
        d = {self.cwl_name: {"in": {inp.name: str(inp.value) for inp in self.inputs}}}
        return d

    def _save_cwl_to_wic_dir(self, wic_dir: StrPath) -> None:
        with open(
            _path(wic_dir).joinpath("cwl_adapters", f"{self.cwl_name}.cwl"),
            "w",
            encoding="utf-8",
        ) as file:
            file.write(yaml.dump(self.yaml))


@dataclass
class Workflow:
    steps: list
    name: str
    wic_dir: Optional[StrPath] = None  # type: ignore
    yml_path: Path = field(init=False)

    def step(self, step_: Step) -> None:
        """Append step to Workflow."""
        self.steps.append(step_)

    def __post_init__(self) -> None:
        for s in self.steps:
            try:
                s._validate()  # pylint: disable=W0212
            except BaseException as exc:
                raise InvalidStep(f"{s.cwl_name} is missing required inputs") from exc

    @property
    def yaml(self) -> Any:
        """WIC YML representation."""
        d = {"steps": [step._yml for step in self.steps]}
        return d

    def _save_yaml(self, wic_dir: StrPath) -> None:
        self.yml_path = _path(wic_dir).joinpath(f"{self.name}.yml")
        with open(self.yml_path, "w", encoding="utf-8") as file:
            file.write(yaml.dump(self.yaml))

    def _save_all_cwl(self, wic_dir: StrPath) -> None:
        for s in self.steps:
            try:
                s._save_cwl_to_wic_dir(wic_dir)
            except BaseException as exc:
                raise exc

    def compile(
        self, wic_dir: Optional[StrPath] = None, out_dir: Optional[StrPath] = None
    ) -> None:
        """Compile Workflow using WIC"""
        # if not cwl_dir:
        #     if not self.cwl_dir:
        #         raise RuntimeError("a CWL directory must be specified")
        #     cwl_dir = self.cwl_dir
        if not wic_dir:
            if not self.wic_dir:
                raise RuntimeError("WIC's directory must be specified")
            wic_dir = self.wic_dir

        self._save_all_cwl(wic_dir)
        self._save_yaml(wic_dir)
        print(f"Compiling {self.name}")
        proc = subprocess.run(
            args=["wic", "--yaml", self.yml_path],
            capture_output=True,
            cwd=wic_dir,
            check=True,
            text=True,
            universal_newlines=True,
        )
        print(proc.stdout)
        if out_dir:
            autogenerated_path = _path(wic_dir).joinpath(
                "autogenerated", f"{self.name}.cwl"
            )
            output = _path(out_dir).joinpath(f"{self.name}.cwl")
            shutil.copy(autogenerated_path, output)
            print(f"Successfully compiled to {str(output)}")
