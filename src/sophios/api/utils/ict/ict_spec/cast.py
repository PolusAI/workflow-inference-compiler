import json
from pathlib import Path
from typing import Union

from yaml import safe_load

from sophios.api.utils.ict.ict_spec.model import ICT


def remove_invalid_keys(ict_data: dict) -> None:
    """remove invalid keys from "ui"""

    if ("ui" in ict_data):

        ui_list = ict_data["ui"]

        # ensure ui is correct data type according to ICT schema
        if (not isinstance(ui_list, list)):
            raise RuntimeError("Error: ui must be a list in ICT data.")

        # remove incorrect keys
        for ui in ui_list:
            ui.pop("format", None)
            ui.pop("required", None)


def cast_to_ict(ict: Union[Path, str, dict]) -> ICT:

    if isinstance(ict, str):
        ict = Path(ict)

    if isinstance(ict, Path):

        if str(ict).endswith(".yaml") or str(ict).endswith(".yml"):
            with open(ict, "r", encoding="utf-8") as f_o:
                data = safe_load(f_o)
        elif str(ict).endswith(".json"):
            with open(ict, "r", encoding="utf-8") as f_o:
                data = json.load(f_o)
        else:
            raise ValueError(f"File extension not supported: {ict}")

        remove_invalid_keys(data)

        return ICT(**data)

    remove_invalid_keys(ict)

    return ICT(**ict)
