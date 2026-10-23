"""What a name written in a document may be spelled with.

A step id, a workflow input or output and a workflow's own name become the
ids of the compiled CWL, part of the names Sophios makes from them (joined
with `~`, numbered with `@`) and part of file names and `run:` text. So a
written name is letters, digits, `_`, `.` and `-`, and does not start with
`.`: `#`, `?` and `%` change what a `run:` reference means, `/` and `\\` make
a directory, `..` names the parent one, and `: * " < > |` are not file names
on Windows. A tool's own port names are not written names and are not checked:
plain CWL passes through.
"""
import re
from typing import Final

#: A name a document may write.
WRITTEN_NAME: Final = re.compile(r'[A-Za-z0-9_-][A-Za-z0-9_.-]*')


def check_written_name(name: str) -> str | None:
    """Why `name` cannot be written, or None when it can.

    Args:
        name (str): A step id, a workflow input or output name, or a workflow's name.

    Returns:
        str | None: The reason, as a clause that follows the name in a message.
    """
    if WRITTEN_NAME.fullmatch(name):
        return None
    return ("is not a name Sophios can write: use letters, digits, '_', '.' and '-', "
            "and do not start with '.'")
