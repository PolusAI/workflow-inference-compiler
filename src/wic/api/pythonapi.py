# pylint: disable=W1203
"""CLT utilities."""
import logging
import subprocess
from dataclasses import field
from functools import singledispatch
from pathlib import Path
from typing import Any, ClassVar, Optional, TypeVar, Union

import cwl_utils.parser as cu_parser
import yaml
from cwl_utils.parser import CommandLineTool as CWLCommandLineTool
from cwl_utils.parser import load_document_by_uri
from pydantic import BaseModel, Field, PrivateAttr  # pylint: disable=E0611
from pydantic.dataclasses import dataclass as pydantic_dataclass

from wic.api._compat import PYDANTIC_V2
from wic.api._types import CWL_TYPES_DICT

if PYDANTIC_V2:
    from pydantic import ConfigDict, field_validator  # pylint: disable=E0611, C0412
else:
    from pydantic import validator

logger = logging.getLogger("WIC Python API")

_WIC_PATH = Path(__file__).parent.parent.parent.parent  # WIC dir


class InvalidInputValueError(Exception):
    pass


class MissingRequiredValueError(Exception):
    pass


class InvalidStepError(Exception):
    pass


class InvalidLinkError(Exception):
    pass


class InvalidCLTError(ValueError):
    pass


CWLInputParameter = Union[
    cu_parser.cwl_v1_0.CommandInputParameter,
    cu_parser.cwl_v1_1.CommandInputParameter,
    cu_parser.cwl_v1_2.CommandInputParameter,
]  # cwl_utils does not have CommandInputParameter union

StrPath = TypeVar("StrPath", str, Path)


class CLTInput(BaseModel):  # pylint: disable=too-few-public-methods
    """Input of CLT."""

    inp_type: object
    name: str
    value: Any = Field(default=None)  # validation happens at assignment
    required: bool = True
    linked: bool = False

    def __init__(self, cwl_inp: CWLInputParameter) -> None:
        # temporary fix for different versions of cwl_utils
        # where it changed from `type_` to `type`
        if hasattr(cwl_inp, "type_"):
            inp_type = cwl_inp.type_
        elif hasattr(cwl_inp, "type"):
            inp_type = cwl_inp.type
        else:
            raise AttributeError("CWLInputParameter has no attribute type or type_")
        if isinstance(inp_type, list) and "null" in inp_type:
            required = False
        else:
            required = True
        super().__init__(inp_type=inp_type, name=cwl_inp.id, required=required)

    if PYDANTIC_V2:
        @field_validator("name", mode="before")
        # type: ignore
        @classmethod
        def get_name_from_id(cls, cwl_id) -> Any:  # pylint: disable=no-self-argument
            """Return name of input from InputParameter.id."""
            return cwl_id.split("#")[-1]

        @field_validator("inp_type", mode="before")
        # type: ignore
        @classmethod
        def set_inp_type(cls, inp) -> object:  # pylint: disable=no-self-argument
            """Return inp_type."""
            if isinstance(inp, list):  # optional inps
                inp = inp[1]
            return CWL_TYPES_DICT[inp]
    else:
        @validator("name")
        # type: ignore
        @classmethod
        def get_name_from_id(cls, cwl_id) -> Any:  # pylint: disable=no-self-argument
            """Return name of input from InputParameter.id."""
            return cwl_id.split("#")[-1]

        @validator("inp_type")
        # type: ignore
        @classmethod
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
                    f"invalid attribute type for {self.name}: "
                    f"got {__value.__class__.__name__}, "
                    f"expected {self.inp_type.__name__}"  # type: ignore # noqa: E266
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
                raise InvalidInputValueError() from exc
            if not value_.is_dir():
                raise InvalidInputValueError(f"{str(value_)} is not a directory")
            return value_
        return value
    return value


def _is_link(s: str) -> bool:
    """Return True if s starts with & or *."""
    if s.startswith("&") or s.startswith("*"):
        return True
    return False


@singledispatch
def _value_str(val: Any) -> Union[str, bool]:
    """Return value of input as str."""
    return str(val)


@_value_str.register
def _(val: bool) -> bool:
    return val


class Step(BaseModel):  # pylint: disable=too-few-public-methods
    """Base class for Step of Workflow."""
    if PYDANTIC_V2:
        model_config: ClassVar[ConfigDict] = ConfigDict(arbitrary_types_allowed=True)
    else:
        class Config:  # pylint: disable=too-few-public-methods
            arbitrary_types_allowed = True

    clt: CWLCommandLineTool
    clt_path: Path
    clt_name: str
    cwl_version: str
    inputs: list[CLTInput]
    yaml: dict[str, Any]
    cfg_yaml: dict = Field(default_factory=_default_dict)
    _input_names: list[str] = PrivateAttr(default_factory=list)

    def __init__(self, clt_path: StrPath, config_path: Optional[StrPath] = None):
        # validate using cwl.utils
        if not isinstance(clt_path, (Path, str)):
            raise TypeError("cwl_path must be a Path or str")
        clt_path_ = Path(clt_path) if isinstance(clt_path, str) else clt_path
        try:
            clt = load_document_by_uri(clt_path_)
        except Exception as exc:
            raise InvalidCLTError(f"invalid cwl file: {clt_path_}") from exc
        with clt_path_.open("r", encoding="utf-8") as file:
            yaml_file = yaml.safe_load(file)
        if config_path:
            cfg_path_ = Path(config_path) if isinstance(config_path, str) else config_path
            with cfg_path_.open("r", encoding="utf-8") as file:
                cfg_yaml = yaml.safe_load(file)
        else:
            cfg_yaml = _default_dict()  # redundant, to avoid it being unbound
        clt_name = clt_path_.stem
        data = {
            "clt": clt,
            "clt_path": clt_path_,
            "cwl_version": clt.cwlVersion,
            "clt_name": clt_name,
            "inputs": clt.inputs,
            "yaml": yaml_file,
            "cfg_yaml": cfg_yaml,
        }
        super().__init__(**data)
        self._input_names = [inp.id.split("#")[-1] for inp in clt.inputs]
        if config_path:
            self._set_from_io_cfg()

    if PYDANTIC_V2:
        @field_validator("inputs", mode="before")
        @classmethod
        def cast_to_clt_input_model(
            cls, cwl_inps: list[CWLInputParameter]
        ) -> list[CLTInput]:  # pylint: disable=no-self-argument
            """Populate inputs from cwl.inputs."""
            return [CLTInput(x) for x in cwl_inps]
    else:
        @validator("inputs", pre=True)
        # type: ignore
        @classmethod
        def cast_to_clt_input_model(
            cls, cwl_inps: list[CWLInputParameter]
        ):  # pylint: disable=no-self-argument
            """Populate inputs from cwl.inputs."""
            return [CLTInput(x) for x in cwl_inps]

    def __repr__(self) -> str:
        repr_ = f"Step(clt_path={self.clt_path.__repr__()})"
        return repr_

    def __setattr__(self, __name: str, __value: Any) -> Any:  # pylint: disable=R1710
        if __name in ["inputs", "yaml", "cfg_yaml", "clt_name", "_input_names",
                      "__private_attributes__", "__pydantic_private__"]:
            return super().__setattr__(__name, __value)
        if hasattr(self, "_input_names") and __name in self._input_names:
            index = self._input_names.index(__name)
            if isinstance(__value, CLTInput):
                if not __value.linked:
                    try:
                        local_input = self.inputs[index]
                        if not local_input.inp_type == __value.inp_type:
                            raise InvalidLinkError(
                                f"links must have the same input type. "
                                f"cannot link {local_input.name} to {__value.name}"
                            )
                        # Use the current value so we can exactly reproduce hand-crafted yml files.
                        # (Very useful for regression testing!)
                        tmp = __value.value if __value.value else f"{__name}{self.clt_name}"
                        local_input._set_value(f"*{tmp}", check=False, linked=True)
                        __value._set_value(f"&{tmp}", check=False, linked=True)
                    except BaseException as exc:
                        raise exc
                else:  # value is already linked to another inp
                    try:
                        local_input = self.inputs[index]
                        if not local_input.inp_type == __value.inp_type:
                            raise InvalidLinkError(
                                f"links must have the same input type. "
                                f"cannot link {local_input.name} to {__value.name}"
                            )
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

    if PYDANTIC_V2:
        def __getattr__(self, __name: str) -> Any:
            if __name in ["__pydantic_private__", "__class__", "__private_attributes__"]:
                return super().__getattribute__(__name)
            if __name != "_input_names" and __name in self._input_names:
                return self.inputs[self._input_names.index(__name)]
            # pydantic has BaseModel.__getattr__ in a
            # non-TYPE_CHECKING block so mypy doesn't see it
            return super().__getattr__(__name)  # type: ignore
    else:
        def __getattribute__(self, __name: str) -> Any:
            if __name != "_input_names" and hasattr(self, "_input_names"):
                if __name in self._input_names:
                    return self.inputs[self._input_names.index(__name)]
            return super().__getattribute__(__name)

    def _set_from_io_cfg(self) -> None:
        for name, value in self.cfg_yaml.items():
            value_ = _get_value_from_cfg(value)
            setattr(self, name, value_)

    def _validate(self) -> None:
        for inp in self.inputs:
            if inp.required and inp.value is None:
                raise MissingRequiredValueError(f"{inp.name} is required")

    @property
    def _yml(self) -> dict:
        d = {
            self.clt_name: {
                "in": {
                    inp.name: _value_str(inp.value)
                    for inp in self.inputs
                    if inp.value is not None
                }
            }
        }
        return d

    def _save_cwl(self, path: Path) -> None:
        cwl_adapters = path.joinpath("cwl_adapters")
        cwl_adapters.mkdir(exist_ok=True, parents=True)
        with open(
            cwl_adapters.joinpath(f"{self.clt_name}.cwl"),
            "w",
            encoding="utf-8",
        ) as file:
            file.write(yaml.dump(self.yaml))


if PYDANTIC_V2:
    DATACLASS_CONFIG = ConfigDict(validate_assignment=True)
else:
    # mypy marks this incorrect if Pydantic V2
    DATACLASS_CONFIG = ConfigDict(validate_on_init=True, validate_assignment=True)  # type: ignore


@pydantic_dataclass(config=DATACLASS_CONFIG)
class Workflow:
    steps: list[Step]
    name: str
    yml_path: Optional[Path] = field(default=None, init=False, repr=False)  # pylint: disable=E3701

    def __post_init__(self) -> None:
        for s in self.steps:
            try:
                s._validate()  # pylint: disable=W0212
            except BaseException as exc:
                raise InvalidStepError(
                    f"{s.clt_name} is missing required inputs"
                ) from exc

    def append(self, step_: Step) -> None:
        """Append step to Workflow."""
        if not isinstance(step_, Step):
            raise TypeError("step must be a Step")
        self.steps.append(step_)

    @property
    def yaml(self) -> dict[str, Any]:
        """WIC YML representation."""
        d = {"steps": [step._yml for step in self.steps]}  # pylint: disable=W0212
        return d

    def _save_yaml(self) -> None:
        _WIC_PATH.mkdir(parents=True, exist_ok=True)
        self.yml_path = _WIC_PATH.joinpath(f"{self.name}.yml")
        with open(self.yml_path, "w", encoding="utf-8") as file:
            file.write(yaml.dump(self.yaml))

    def _save_all_cwl(self) -> None:
        """Save CWL files to cwl_adapters.

        This is necessary for WIC to compile the workflow.
        """
        _WIC_PATH.mkdir(parents=True, exist_ok=True)
        for s in self.steps:
            try:
                s._save_cwl(_WIC_PATH)  # pylint: disable=W0212
            except BaseException as exc:
                raise exc

    def compile(self, run_local: bool = False) -> Path:
        """Compile Workflow using WIC.

        Returns path to compiled CWL Workflow.
        """

        self._save_all_cwl()
        self._save_yaml()
        logger.info(f"Compiling {self.name}")
        args = ["wic", "--yaml", f"{self.name}.yml"]
        if run_local:
            args.append('--run_local')
        try:
            proc = subprocess.run(
                args=args,
                capture_output=True,
                cwd=_WIC_PATH,
                check=True,
                text=True,
                universal_newlines=True,
            )
        except subprocess.CalledProcessError as exc:
            logger.error(exc.stderr)
            logger.info(exc.stdout)
            raise exc
        logger.info(proc.stdout)
        return _WIC_PATH.joinpath("autogenerated", f"{self.name}.cwl")

    def run(self) -> None:
        """Run compiled workflow."""
        logger.info(f"Running {self.name}")
        self.compile(run_local=True)
