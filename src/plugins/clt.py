"""CLT utilities."""
from dataclasses import dataclass
from functools import singledispatch
from pathlib import Path
from typing import Any, Generic, List, Optional, TypeVar, Union

import cwl_utils.parser as cu_parser
import typeguard
import yaml

# from ._types import CWLTypesEnum
from _types import CWL_TYPES_DICT, CWLTypesEnum
from cwl_utils.parser import CommandLineTool as CWLCommandLineTool
from cwl_utils.parser import DockerRequirement, load_document_by_uri
from pydantic import (  # pylint: disable=E0611
    BaseModel,
    Extra,
    Field,
    PrivateAttr,
    ValidationError,
    validator,
)
from pydantic.dataclasses import dataclass as pydanticdataclass
from pydantic.error_wrappers import ErrorWrapper  # pylint: disable=E0611

p = Path(__file__).absolute().with_name("ome.cwl")
cfg = Path(__file__).absolute().with_name("inp.yml")
p2 = Path(__file__).absolute().with_name("montage.cwl")
p3 = Path(__file__).absolute().with_name("filerenaming.cwl")


class InvalidInputValue(Exception):
    pass


class MissingRequiredValue(Exception):
    pass


class InvalidStep(Exception):
    pass


StrPath = TypeVar("StrPath", str, Path)
CWLInputParameter = Union[
    cu_parser.cwl_v1_0.CommandInputParameter,
    cu_parser.cwl_v1_1.CommandInputParameter,
    cu_parser.cwl_v1_2.CommandInputParameter,
]  # cwl_utils does not CommandInputParameter union


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
        assert len(docker) > 0, "no docker requirement specified"
        # schema salad already checks that len(docker) < 2
        docker = docker[0].dockerPull
        return docker


# TODO

T = TypeVar("T")  # Generic


class CLTInput(BaseModel):  # pylint: disable=too-few-public-methods
    """Input of CLT."""

    inp_type: object
    name: str
    value: Any  # validation happens at assignment
    required: bool = True
    # TODO required check yes

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
    def set_inp_type(cls, inp) -> object:
        """Return inp_type."""
        if isinstance(inp, list):  # optional inps
            inp = inp[1]
        return CWL_TYPES_DICT[inp]

    # def __setattr__(self, __name: str, __value: Any) -> None:
    #     if __name == "value":
    #         if not isinstance(__value, self.inp_type):
    #             raise TypeError(
    #                 f"invalid attribute type for {self.name}: got {__value.__class__.__name__}, expected {self.inp_type.__name__}"
    #             )
    #         return super().__setattr__(__name, __value)
    #     else:
    #         return super().__setattr__(__name, __value)
    def _set_value(self, __value: Any, check: bool = True) -> None:
        """Set input value."""
        if check:
            if not isinstance(__value, self.inp_type):  # type: ignore
                raise TypeError(
                    f"invalid attribute type for {self.name}: got {__value.__class__.__name__}, expected {self.inp_type.__name__}"  # type: ignore # noqa: E266
                )
            self.value = __value
            return
        self.value = __value  # to be used for linking inputs */&

    # def __eq__(self, other: str) -> bool:
    #     if not isinstance(other, str):
    #         raise NotImplementedError()
    #     return self.name == other


def _default_dict() -> dict:
    return {}


def _get_value_cfg(value: Any) -> Any:  # validation happens in Step
    if isinstance(value, dict):
        if "Directory" in value.values():
            try:
                value_ = Path(value["path"])
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
        print(type(__value), __name)
        print(isinstance(__value, CLTInput))
        print(f"{__name} {__name in self._input_names}")
        if __name in self._input_names:
            index = self._input_names.index(__name)
            if isinstance(__value, CLTInput):
                print("sdfdfa")
                try:
                    print("HOLA")
                    print(f"{type(__value)}")
                    tmp = f"{__name}{self.cwl_name}"
                    print(tmp)
                    self.inputs[index]._set_value(f"*{tmp}", check=False)
                    __value.value = f"&{tmp}"
                except BaseException as exc:
                    raise exc
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
            value_ = _get_value_cfg(value)
            setattr(self, name, value_)

    def _validate(self) -> None:
        for inp in self.inputs:
            if inp.required and inp.value is None:
                raise MissingRequiredValue(f"{inp.name} is required")

    @property
    def _step_yml(self) -> dict:
        d = {self.cwl_name: {"in": {inp.name: str(inp.value) for inp in self.inputs}}}
        return d


# t = Tool(p)
# c = Step(p, cfg)
step1 = Step(p)
io1 = step1.inpDir
step2 = Step(p3)

step1.inpDir = Path("/Users/camilovelezr/axle/workflow")
step1.filePattern = ".*tif"
step1.fileExtension = ".ome.tif"

step2.inpDir = step1.outDir
step2.outDir = Path("/Users/camilovelezr/axle")
step2.filePattern = ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
step2.outFilePattern = "test"


@dataclass
class Workflow:
    steps: list

    def step(self, step: Step) -> None:
        """Append step to Workflow."""
        self.steps.append(step)

    def __post_init__(self) -> None:
        for step in self.steps:
            try:
                step._validate()  # pylint: disable=W0212
            except BaseException as exc:
                raise InvalidStep(
                    f"{step.cwl_name} is missing required inputs"
                ) from exc

    # def __init__(self, step_list: list) -> None:
    #     for step in step_list:
    #         print(step.cwl_name)
    #         try:
    #             step._validate()  # pylint: disable=W0212
    #         except BaseException as exc:
    #             raise InvalidStep(f"{step.cwl_name} is invalid") from exc
    #         super().__init__()


w1 = Workflow([step1, step2])

t1 = step1._step_yml
t2 = yaml.dump(t1)
t3 = {"steps": [t1]}
t4 = yaml.dump(t3, indent=2)
step1._validate()
step2._validate()
2
