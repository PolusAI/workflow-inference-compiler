"""CWL Types."""
from enum import Enum
from pathlib import Path


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


CWL_TYPES_DICT: dict[str, object] = {
    "null": None,
    "boolean": bool,
    "int": int,
    "long": int,
    "float": float,
    "double": float,
    "string": str,
    "File": Path,
    "Directory": Path,
}
