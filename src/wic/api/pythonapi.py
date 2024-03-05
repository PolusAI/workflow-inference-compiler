# pylint: disable=W1203
"""CLT utilities."""
import logging
from dataclasses import field
from functools import singledispatch
from pathlib import Path
import subprocess as sub
from typing import Any, ClassVar, Optional, TypeVar, Union

import cwl_utils.parser as cu_parser
import yaml
from cwl_utils.parser import CommandLineTool as CWLCommandLineTool
from cwl_utils.parser import load_document_by_uri
from pydantic import BaseModel, ConfigDict, Field, PrivateAttr, field_validator
from pydantic.dataclasses import dataclass as pydantic_dataclass

from wic.api._types import CWL_TYPES_DICT
from wic import compiler, input_output, plugins
from wic import run_local as run_local_module
from wic.schemas.wic_schema import get_args
from wic.utils_graphs import get_graph_reps
from wic.wic_types import CompilerInfo, RoseTree, StepId, YamlTree


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

    def _save_cwl(self, path: Path) -> None:
        cwl_adapters = path.joinpath("cwl_adapters")
        cwl_adapters.mkdir(exist_ok=True, parents=True)
        with open(
            cwl_adapters.joinpath(f"{self.clt_name}.cwl"),
            "w",
            encoding="utf-8",
        ) as file:
            file.write(yaml.dump(self.yaml))


DATACLASS_CONFIG = ConfigDict(validate_assignment=True)


@pydantic_dataclass(config=DATACLASS_CONFIG)
class Workflow:
    steps: list[Step]
    name: str
    yml_path: Optional[Path] = field(default=None, init=False, repr=False)

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

    def compile(self, write_to_disk: bool = False) -> CompilerInfo:
        """Compile Workflow using WIC.

        Args:
            write_to_disk (bool, optional): Save the compiled CWL Workflow to disk. Defaults to False.

        Returns:
            CompilerInfo: Contains the data associated with compiled subworkflows\n
            (in the Rose Tree) together with mutable cumulative environment\n
            information which needs to be passed through the recursion.
        """
        args = get_args(self.name)  # Use mock CLI args

        # We do NOT need to save anything to disk
        tools_cwl = plugins.get_tools_cwl(args.homedir, quiet=args.quiet)

        graph = get_graph_reps(self.name)
        # NOTE: Once the Workflow class supports subworkflows, we will need to
        # recursively convert sub-Workflows to YamlTrees here.
        yaml_tree = YamlTree(StepId(self.name, 'global'), self.yaml)

        # The compile_workflow function is 100% in-memory
        compiler_info = compiler.compile_workflow(yaml_tree, args, [], [graph], {}, {}, {}, {},
                                                  tools_cwl, True, relative_run_path=True, testing=False)

        if write_to_disk:
            # Now we can choose whether to write_to_disk or not
            rose_tree: RoseTree = compiler_info.rose
            input_output.write_to_disk(rose_tree, Path('autogenerated/'), relative_run_path=True)

        return compiler_info

    def run(self) -> None:
        """Run compiled workflow."""
        logger.info(f"Running {self.name}")
        compiler_info = self.compile(write_to_disk=True)

        args = get_args(self.name)  # Use mock CLI args
        rose_tree: RoseTree = compiler_info.rose

        # If you don't like it, you can programmatically overwrite anything in args
        # args.no_docker_remove_entrypoints = True
        if not args.no_docker_remove_entrypoints:
            # Requires root, so guard behind CLI option
            root_project_dir = Path(__file__).parent.parent.parent.parent
            cmd = ['python', str(root_project_dir / 'docker_remove_entrypoints.py')]
            sub.run(cmd, check=True, capture_output=True)

            rose_tree = plugins.dockerPull_append_noentrypoint_rosetree(rose_tree)
            input_output.write_to_disk(rose_tree, Path('autogenerated/'), relative_run_path=True)

        # Do NOT capture stdout and/or stderr and pipe warnings and errors into a black hole.
        retval = run_local_module.run_local(args, rose_tree, None, 'cwltool', True)
