import json
# import subprocess as sub
from pathlib import Path
# import signal
# import sys
# from typing import List
# import argparse
import asyncio
from jsonschema import Draft202012Validator

from fastapi import Request

import pytest
from sophios.wic_types import Json


from sophios.api.http import restapi

SCHEMA = {
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
            "required": ["id", "inletIndex", "outletIndex", "sourceId", "targetId"]
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


@pytest.mark.fast
def test_rest_core_single_node() -> None:
    """A simple single node 'hello world' test"""
    # validate schema
    Draft202012Validator.check_schema(SCHEMA)
    df2012 = Draft202012Validator(SCHEMA)
    inp_file = "single_node_helloworld.json"
    inp: Json = {}
    yaml_path = "workflow.json"
    inp_path = Path(__file__).with_name(inp_file)
    with open(inp_path, 'r', encoding='utf-8') as f:
        inp = json.load(f)
    # check if object is conformant with our schema
    df2012.is_valid(inp)
    print('----------- from rest api ----------- \n\n')
    scope = {}
    scope['type'] = 'http'

    async def receive() -> Json:
        inp_byte = json.dumps(inp).encode('utf-8')
        return {"type": "http.request", "body": inp_byte}

    # create a request object and pack it with our json payload
    req: Request = Request(scope)
    req._receive = receive
    res: Json = asyncio.run(restapi.compile_wf(req))  # call to rest api
    assert int(res['retval']) == 0
