import sys
import copy
from typing import Any, Dict, List

from jsonschema import Draft202012Validator

from sophios import cli
from sophios.wic_types import Json

SCHEMA: Json = {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "definitions": {
        "Link": {
            "properties": {
                "id": {
                    "type": "number"
                },
                "inletIndex": {
                    "type": "number"
                },
                "outletIndex": {
                    "type": "number"
                },
                "sourceId": {
                    "type": "number"
                },
                "targetId": {
                    "type": "number"
                },
                "x1": {
                    "type": "number"
                },
                "x2": {
                    "type": "number"
                },
                "y1": {
                    "type": "number"
                },
                "y2": {
                    "type": "number"
                }
            },
            "type": "object",
            "required": ["id", "sourceId", "targetId"]
        },
        "NodeSettings": {
            "properties": {
                "inputs": {
                    "additionalProperties": {
                        "$ref": "#/definitions/T"
                    },
                    "type": "object"
                },
                "outputs": {
                    "additionalProperties": {
                        "$ref": "#/definitions/T"
                    },
                    "type": "object"
                }
            },
            "type": "object"
        },
        "NodeX": {
            "properties": {
                "expanded": {
                    "type": "boolean"
                },
                "height": {
                    "type": "number"
                },
                "id": {
                    "type": "number"
                },
                "internal": {
                    "type": "boolean"
                },
                "name": {
                    "type": "string"
                },
                "pluginId": {
                    "type": "string"
                },
                "settings": {
                    "$ref": "#/definitions/NodeSettings"
                },
                "width": {
                    "type": "number"
                },
                "x": {
                    "type": "number"
                },
                "y": {
                    "type": "number"
                },
                "z": {
                    "type": "number"
                },
            },
            "type": "object",
            "required": ["id", "name", "pluginId", "settings", "internal"]
        },
        "T": {
            "type": "object"
        }
    },
    "properties": {
        "links": {
            "items": {
                "$ref": "#/definitions/Link"
            },
            "type": "array"
        },
        "nodes": {
            "items": {
                "$ref": "#/definitions/NodeX"
            },
            "type": "array"
        },
        "selection": {
            "items": {
                "type": "number"
            },
            "type": "array"
        }
    },
    "type": "object",
    "required": ["links", "nodes"]
}


def del_irrelevant_keys(ldict: List[Dict[Any, Any]], relevant_keys: List[Any]) -> None:
    """deletes irrelevant keys from every dict in the list of dicts"""
    for elem in ldict:
        ekeys = list(elem.keys())
        for ek in ekeys:
            if ek not in relevant_keys:
                # delete the key if it exists
                elem.pop(ek, None)


def validate_schema_and_object(schema: Json, jobj: Json) -> None:
    Draft202012Validator.check_schema(schema)
    df2012 = Draft202012Validator(schema)
    df2012.is_valid(jobj)


def raw_wfb_to_lean_wfb(inp: Json) -> Json:
    """drop all the unnecessary info from incoming wfb object"""
    inp_restrict = copy.deepcopy(inp)
    keys = list(inp.keys())
    # To avoid deserialization
    # required attributes from schema
    prop_req = SCHEMA['required']
    nodes_req = SCHEMA['definitions']['NodeX']['required']
    links_req = SCHEMA['definitions']['Link']['required']
    do_not_rem_nodes_prop = ['cwlScript']
    do_not_rem_links_prop: list = []

    for k in keys:
        if k not in prop_req:
            del inp_restrict[k]
        elif k == 'links':
            lems = inp_restrict[k]
            rel_links_keys = links_req + do_not_rem_links_prop
            del_irrelevant_keys(lems, rel_links_keys)
        elif k == 'nodes':
            nems = inp_restrict[k]
            rel_nodes_keys = nodes_req + do_not_rem_nodes_prop
            del_irrelevant_keys(nems, rel_nodes_keys)
        else:
            pass

    validate_schema_and_object(SCHEMA, inp_restrict)
    return inp_restrict


def wfb_to_wic(request: Json) -> Json:
    """Convert the json object from http request object to a json object that can be used as input to wic compiler.

    Args:
        request (Json): json object from http request

    Returns:
        converted_json (Json): json object that can be used as input to wic compiler"""

    converted_steps: list[Any] = []

    for step in request['steps']:
        step_template = step['template']
        arguments = step['arguments']
        # Get the template name from the step template
        template_name = next((tval['name']
                             for tname, tval in request['templates'].items() if step_template == tname), None)
        # template_name = None
        # for tname, tval in request['templates'].items():
        #     if tname == step_template and tval['name']:
        #         template_name = tval['name']
        #         break
        #     elif tname == step_template and not tval['name']:
        #         break
        #     else:
        #         pass

        converted_step: Json = {}
        if template_name:
            converted_step[template_name] = {
                "in": {}
            }

        for key, value in arguments.items():
            # be aware of usage of periods in the values as delimiters, this may cause an issue when storing in MongoDB
            if value.startswith("steps."):
                parts = value.split('.')
                src_step_idx = int(parts[1][len("step"):])
                src_output_key = parts[3]

                # Get the src template name from the src step name
                src_template = next((step.get("template")
                                    for step in request['steps'] if step.get("name") == parts[1]), None)
                src_template_name = next((stval['name'] for stname, stval in request['templates'].items()
                                          if src_template == stname), None)
                src_converted_step = next((step.get(src_template_name)
                                          for step in converted_steps if step.get(src_template_name)), None)
                if src_converted_step:
                    src_converted_step["in"][src_output_key] = f"&{src_template_name}.{src_output_key}.{src_step_idx}"
                    converted_step[template_name]["in"][key] = f"*{src_template_name}.{src_output_key}.{src_step_idx}"
            else:
                converted_step[template_name]["in"][key] = value
        converted_steps.append(converted_step)

    converted_json: Json = {"steps": converted_steps}
    return converted_json
