"""CLT utilities."""
from dataclasses import field
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
    ValidationError,
    validator,
)
from pydantic.dataclasses import dataclass as pydanticdataclass
from pydantic.error_wrappers import ErrorWrapper  # pylint: disable=E0611

p = Path(__file__).absolute().with_name("ome.cwl")
cfg = Path(__file__).absolute().with_name("inp.yml")
# p = Path(__name__).absolute().parent.joinpath("src", "plugins", "ome.cwl")
print(p)

# c = load_document_by_uri(p)

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

    # def __init__(self, path: StrPath):
    #     # validate using cwl.utils
    #     try:
    #         cwl = load_document_by_uri(path)
    #     except BaseException as s_e:
    #         e_w = ErrorWrapper(s_e, "invalid cwl file")  # type: ignore
    #         raise ValidationError([e_w], Tool)  # pylint: disable=raise-missing-from
    #     docker = [x for x in cwl.requirements if _is_docker_requirement(x)]
    #     super().__init__(
    #         path=path,
    #         cwlVersion=cwl.cwlVersion,
    #         dockerContainer=docker,
    #         cwl=cwl,
    #     )

    @validator("dockerContainer", pre=True, check_fields=False)
    # type: ignore
    def validate_docker(cls, docker):  # pylint: disable=no-self-argument
        """Validate that there is one docker requirement specified."""
        assert len(docker) > 0, "no docker requirement specified"
        # schema salad already checks that len(docker) < 2
        docker = docker[0].dockerPull
        return docker


# class CLTInputType(BaseModel):  # pylint: disable=R0903
#     """Input type of CWL input."""

#     obj_type: CWLTypesEnum


# TODO

T = TypeVar("T")  # Generic


class CLTInput(BaseModel):  # pylint: disable=too-few-public-methods
    """Input of CLT."""

    inp_type: object
    name: str
    value: Any  # validation happens at assignemnt

    def __init__(self, cwl_inp: CWLInputParameter) -> None:
        super().__init__(inp_type=cwl_inp.type, name=cwl_inp.id)

    @validator("name")
    # type: ignore
    def get_name_from_id(cls, cwl_id) -> str:  # pylint: disable=no-self-argument
        """Return name of input from InputParameter.id."""
        return cwl_id.split("#")[-1]

    @validator("inp_type", pre=True)
    # type: ignore
    def set_inp_type(cls, inp) -> object:
        """Return inp_type."""
        return CWL_TYPES_DICT[inp]

    def __setattr__(self, __name: str, __value: Any) -> None:
        if __name == "value":
            if not isinstance(__value, self.inp_type):
                raise TypeError(
                    f"invalid attribute type for {self.name}: got {__value.__class__.__name__}, expected {self.inp_type.__name__}"
                )
            return super().__setattr__(__name, __value)
        else:
            return super().__setattr__(__name, __value)

    def __eq__(self, other: str) -> bool:
        if not isinstance(other, str):
            raise NotImplementedError()
        return self.name == other


def _default_dict():
    return {}


class Step(Tool):  # pylint: disable=too-few-public-methods
    """Base class for configured CLTs."""

    inputs: List[CLTInput]
    yaml: dict
    cfg_yaml: dict = Field(default_factory=_default_dict)

    def __init__(self, cwl_path: StrPath, config_path: StrPath):
        # validate using cwl.utils
        try:
            cwl = load_document_by_uri(cwl_path)
        except BaseException as s_e:
            e_w = ErrorWrapper(s_e, "invalid cwl file")  # type: ignore
            raise ValidationError([e_w], Step)  # pylint: disable=raise-missing-from
        docker = [x for x in cwl.requirements if _is_docker_requirement(x)]
        with open(cwl_path, "r", encoding="utf-8") as file:
            yaml_file = yaml.safe_load(file)
        with open(config_path, "r", encoding="utf-8") as file:
            cfg_yaml = yaml.safe_load(file)
        super().__init__(  # type: ignore
            cwl_path=cwl_path,  # type: ignore
            cwlVersion=cwl.cwlVersion,
            dockerContainer=docker,  # type: ignore
            cwl=cwl,
            inputs=cwl.inputs,  # type: ignore
            yaml=yaml_file,
            cfg_yaml=cfg_yaml,
        )

    @validator("inputs", pre=True)
    # type: ignore
    def cast_to_clt_input_model(
        cls, cwl_inps: List[CWLInputParameter]
    ):  # pylint: disable=no-self-argument
        """Populate inputs from cwl.inputs."""
        return [CLTInput(x) for x in cwl_inps]

    def __setattr__(self, __name: str, __value: Any) -> Any:  # pylint: disable=R1710
        if __name in self.inputs:
            index = self.inputs.index(__name)  # type: ignore
            self.inputs[index].value = __value
        else:
            return super().__setattr__(__name, __value)


class Workflow(list):
    pass


# t = Tool(p)
c = Step(p, cfg)

s = c.inputs[0]

s.value = "yes"

r = s == "fileExtension"

c.fileExtension = ".ome.zar"
c.inpDir = Path("/Users/camilovelezr/axle/workflow")

y = c.cwl.__dict__
s

2
