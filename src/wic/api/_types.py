"""CWL Types."""
from enum import Enum
from pathlib import Path

from pydantic import BaseModel, Field  # pylint: disable=E0611

# class CWLFileDir(BaseModel):  # pylint: disable=R0903
#     class_: str
#     path: Path


# class CWLFile(CWLFileDir):  # pylint: disable=R0903
#     class_ = Field(default="File", const=True)


# class CWLDir(CWLFileDir):  # pylint: disable=R0903
#     class_ = Field(default="Directory", const=True)


class CWLTypesEnum(str, Enum):
    NULL = "null"
    BOOLEAN = "boolean"
    INT = "int"
    LONG = "long"
    FLOAT = "float"
    DOUBLE = "double"
    STRING = "string"
    FILE = "File"
    DIRECTORY = "Directory"


# NOTE: use strict types

CWL_TYPES_DICT: dict[str, object] = {
    "null": None,
    "boolean": bool,
    "int": int,
    "long": int,
    "float": float,
    "double": float,
    "string": str,
    # "File": CWLFile,
    # "Directory": CWLDir,
    "File": Path,
    "Directory": Path,
}

# cf = CWLFile(path="/Users/camilovelezr/polus-plugins/n.txt")
# cd = CWLDir(path="/Users/camilovelezr/polus-plugins")
2
