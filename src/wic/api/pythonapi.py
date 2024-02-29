# pylint: disable=W1203
"""CLT utilities."""
import json
import logging
import subprocess
from dataclasses import field
from functools import singledispatch
from pathlib import Path
from typing import Any, ClassVar, Generic, Optional, TypeVar, Union

import cwl_utils.parser as cu_parser
import yaml
from cwl_utils.parser import CommandLineTool as CWLCommandLineTool
from cwl_utils.parser import load_document_by_uri
from pydantic import (  # pylint: disable=E0611
    BaseModel,
    ConfigDict,
    Field,
    PrivateAttr,
    ValidationError,
)
from pydantic.dataclasses import dataclass as pydantic_dataclass

from wic.api._types import CWL_TYPES_DICT


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

    inp_type: Any
    name: str
    value: Any = Field(default=None)  # validation happens at assignment
    required: bool = True
    linked: bool = False

    def __init__(self, cwl_inp: CWLInputParameter) -> None:
        inp_type = cwl_inp.type_  # NOTE: .type in version 0.31 only!
        if isinstance(inp_type, list) and "null" in inp_type:
            required = False
        else:
            required = True
        super().__init__(inp_type=inp_type, name=cwl_inp.id, required=required)

    @field_validator("name", mode="before")
    @classmethod
    def get_name_from_id(cls, cwl_id: str) -> Any:
        """Return name of input from InputParameter.id."""
        return cwl_id.split("#")[-1]

    @field_validator("inp_type", mode="before")
    @classmethod
    def set_inp_type(cls, inp: Any) -> Any:
        """Return inp_type."""
        if isinstance(inp, list):  # optional inps
            inp = inp[1]
        return CWL_TYPES_DICT[inp]

    def _set_value(
        self, __value: Any, check: bool = True, linked: bool = False
    ) -> None:
        """Set input value."""
        if check:
            if not isinstance(__value, self.inp_type):
                raise TypeError(
                    f"invalid attribute type for {self.name}: "
                    f"got {__value.__class__.__name__}, "
                    f"expected {self.inp_type.__name__}"
                )
            self.value = __value
            return
        self.value = __value  # to be used for linking inputs */&
        if linked:
            self.linked = True


InputType = TypeVar("InputType")


class _ConfiguredStepInput(BaseModel, Generic[InputType]):
    """Input of Step used for configuration from JSON.

    This Generic Pydantic Model is used to configure inputs of
    Steps when using a JSON configuration file.
    """
    value: InputType


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
def _yaml_value(val: Any) -> Union[str, bool, int, float]:
    """Convert value to YAML compatible value."""
    return str(val)


@_yaml_value.register
def _(val: int) -> int:
    return val


@_yaml_value.register
def _(val: float) -> float:
    return val


@_yaml_value.register
def _(val: bool) -> bool:
    return val


class Step(BaseModel):  # pylint: disable=too-few-public-methods
    """Base class for Step of Workflow."""
    model_config: ClassVar[ConfigDict] = ConfigDict(arbitrary_types_allowed=True)

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

    @field_validator("inputs", mode="before")
    @classmethod
    def cast_to_clt_input_model(
        cls, cwl_inps: list[CWLInputParameter]
    ) -> list[CLTInput]:
        """Populate inputs from cwl.inputs."""
        return [CLTInput(x) for x in cwl_inps]

    def __repr__(self) -> str:
        repr_ = f"Step(clt_path={self.clt_path.__repr__()})"
        return repr_

    def __setattr__(self, __name: str, __value: Any) -> Any:
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

            elif isinstance(__value, _ConfiguredStepInput):  # when configuring from JSON
                self.inputs[index]._set_value(__value.value)
            else:
                if isinstance(__value, str) and _is_link(__value):
                    self.inputs[index]._set_value(__value, check=False)
                else:
                    self.inputs[index]._set_value(__value)
        else:
            return super().__setattr__(__name, __value)

    def __getattr__(self, __name: str) -> Any:
        if __name in ["__pydantic_private__", "__class__", "__private_attributes__"]:
            return super().__getattribute__(__name)
        if __name != "_input_names" and __name in self._input_names:  # and hasattr(self, "_input_names") ?
            return self.inputs[self._input_names.index(__name)]
        # https://github.com/pydantic/pydantic/blob/812516d71a8696d5e29c5bdab40336d82ccde412/pydantic/main.py#L743-744
        # "We put `__getattr__` in a non-TYPE_CHECKING block because otherwise, mypy allows arbitrary attribute access
        # The same goes for __setattr__ and __delattr__, see: https://github.com/pydantic/pydantic/issues/8643"
        return super().__getattr__(__name)  # type: ignore

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
                    inp.name: _yaml_value(inp.value)
                    for inp in self.inputs
                    if inp.value is not None
                }
            }
        }
        return d

    def _save_cwl(self, path: Path, overwrite: bool) -> None:
        cwl_adapters = path.joinpath("cwl_adapters")
        cwl_adapters.mkdir(exist_ok=True, parents=True)
        if not overwrite:
            if cwl_adapters.joinpath(f"{self.clt_name}.cwl").exists():
                return
        # does not exist:
        with open(
            cwl_adapters.joinpath(f"{self.clt_name}.cwl"),
            "w",
            encoding="utf-8",
        ) as file:
            file.write(yaml.dump(self.yaml))


def _load_config(config_path: StrPath) -> Any:
    """Load configuration JSON file.

    Returns `dict[str, Any]`.
    Type annotation is `Any` because `json.load` type annotation is `Any`.
    """
    with open(config_path, "r", encoding="utf-8") as file:
        cfg = json.load(file)
    return cfg


def configure_step(clt_path: StrPath, config_path: StrPath) -> Step:
    """Configure Step from JSON configuration file.

    This function takes a path to a CLT (Command Line Tool) file and a path to a JSON configuration file,
    and configures a `Step` object based on the provided configuration.

    Args:
        clt_path (StrPath): The path to the CLT file.
        config_path (StrPath): The path to the JSON configuration file.


    Returns:
        Step: The configured `Step` object.

    """
    clt = Step(clt_path)
    input_names = clt._input_names  # pylint: disable=W0212
    config = config_path if isinstance(config_path, dict) else _load_config(config_path)
    # dict with only the inputs that are in the clt
    inp_dict = dict(filter(lambda x: x[0] in input_names, config.items()))
    for inp_name, inp_value in inp_dict.items():
        if isinstance(inp_value, str) and _is_link(inp_value):  # linked * or & inp
            # do not perform validation for linked inputs
            setattr(clt, inp_name, inp_value)
        else:  # not linked inp
            # generic_input performs validation with Pydantic
            inp_type = clt.inputs[input_names.index(inp_name)].inp_type
            try:
                # inp_type is a variable, mypy complains
                generic_input = _ConfiguredStepInput[inp_type](value=inp_value)  # type: ignore
            except ValidationError as exc:
                raise InvalidInputValueError(f"invalid value for {inp_name}") from exc
            setattr(clt, inp_name, generic_input)

    return clt


DATACLASS_CONFIG = ConfigDict(validate_assignment=True)


@pydantic_dataclass(config=DATACLASS_CONFIG)
class Workflow:
    steps: list[Step]
    name: str
    path: Path = field(default=Path.cwd())
    yml_path: Optional[Path] = field(default=None, init=False, repr=False)  # pylint: disable=E3701

    def __post_init__(self) -> None:
        for s in self.steps:
            try:
                s._validate()  # pylint: disable=W0212
            except BaseException as exc:
                raise InvalidStepError(
                    f"{s.clt_name} is missing required inputs"
                ) from exc
        self.path.joinpath(self.name).mkdir(parents=True, exist_ok=True)

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

       # copy to wf path, only informational, not used by WIC
        with open(self.path.joinpath(self.name, f"{self.name}.yml"), "w", encoding="utf-8") as file:
            file.write(yaml.dump(self.yaml))

    def _save_all_cwl(self, overwrite: bool) -> None:
        """Save CWL files to cwl_adapters.

        This is necessary for WIC to compile the workflow.
        """
        _WIC_PATH.mkdir(parents=True, exist_ok=True)
        for s in self.steps:
            try:
                s._save_cwl(_WIC_PATH, overwrite)  # pylint: disable=W0212
            except BaseException as exc:
                raise exc

       # copy to wf path, only informational, not used by WIC
        for s in self.steps:
            try:
                s._save_cwl(self.path.joinpath(self.name), overwrite)  # pylint: disable=W0212
            except BaseException as exc:
                raise exc

    def compile(self, run_local: bool = False, overwrite: bool = False) -> Path:
        """Compile Workflow using WIC.

        Returns path to compiled CWL Workflow.
        """

        self._save_all_cwl(overwrite)
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
            logger.error(exc.stdout)
            return None
        logger.info(proc.stdout)
        # copy files to wf path
        compiled_cwl_path = _WIC_PATH.joinpath("autogenerated", f"{self.name}.cwl")
        self.path.joinpath(self.name, f"{self.name}.cwl").symlink_to(compiled_cwl_path)
        return compiled_cwl_path

    def run(self, overwrite: bool = False) -> None:
        """Run compiled workflow."""
        logger.info(f"Running {self.name}")
        self.compile(run_local=True, overwrite=overwrite)


def configure_workflow(  # pylint: disable=R0913
        clt_list: list[StrPath],
        config_list: list[StrPath],
        name: str,
        path: Union[str, Path] = Path.cwd(),
        compile_workflow: bool = True,
        run: bool = False,
        overwrite: bool = False) -> Workflow:
    """Configure Workflow from list of CLT and configuration files.

    This function takes a list of CLT (Command Line Tool) files and a list of configuration files
    and creates a Workflow object. Each CLT file is paired with a corresponding configuration file
    based on their order in the lists. The Workflow object is then returned.

    Args:
        clt_list (list[StrPath]): A list of paths to CLT files.
        config_list (list[StrPath]): A list of paths to configuration files.
        name (str): The name of the Workflow.
        path (Optional[StrPath]): The path where the Workflow will be created. Defaults to `pathlib.Path.cwd()`.
        compile_workflow (bool): Flag indicating whether to compile the Workflow using WIC. Defaults to True.
        run (bool): Flag indicating whether to run the Workflow after configuration using WIC. Defaults to False.
        overwrite (bool): Flag indicating whether to overwrite existing CLT files in cwl_adapters. Defaults to False.

    Returns:
        Workflow: The configured Workflow object.

    """
    steps = [configure_step(clt, config) for clt, config in zip(clt_list, config_list)]
    path_ = Path(path)
    # mypy incorrectly marks this init
    # of Workflow (pydantic dataclass) as invalid
    # having 'too many arguments'
    workflow = Workflow(steps, name, path_)  # type: ignore
    if run:
        workflow.run(overwrite=overwrite)
    elif compile_workflow:
        workflow.compile(overwrite=overwrite)
    return workflow
