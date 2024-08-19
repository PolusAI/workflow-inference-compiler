import copy
from typing import Any, Dict, List
import yaml
from jsonschema import Draft202012Validator
from sophios.utils_yaml import wic_loader

from sophios.wic_types import Json, Cwl

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
    """Validate schema object"""
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
    do_not_rem_nodes_prop = ['cwlScript', 'run']
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


def wfb_to_wic(inp: Json) -> Cwl:
    """convert lean wfb json to compliant wic"""
    # non-schema preserving changes
    inp_restrict = copy.deepcopy(inp)

    for node in inp_restrict['nodes']:
        if node.get('settings'):
            node['in'] = node['settings'].get('inputs')
            if node['settings'].get('outputs'):
                node['out'] = list({k: yaml.load('!& ' + v, Loader=wic_loader())} for k, v in node['settings']
                                   ['outputs'].items())  # outputs always have to be list
            # remove these (now) superfluous keys
            node.pop('settings', None)
            node.pop('pluginId', None)
            node.pop('internal', None)

    # setting the inputs of the non-sink nodes i.e. whose input doesn't depend on any other node's output
    # first get all target node ids
    target_node_ids = []
    for edg in inp_restrict['links']:
        target_node_ids.append(edg['targetId'])
    # now set inputs on non-sink nodes as inline input '!ii '
    # if inputs exist
    non_sink_nodes = [node for node in inp_restrict['nodes'] if node['id'] not in target_node_ids]
    for node in non_sink_nodes:
        if node.get('in'):
            for nkey in node['in']:
                node['in'][nkey] = yaml.load('!ii ' + node['in'][nkey], Loader=wic_loader())

    # After outs are set
    for edg in inp_restrict['links']:
        # links = edge. nodes and edges is the correct terminology!
        src_id = edg['sourceId']
        tgt_id = edg['targetId']
        src_node = next((node for node in inp_restrict['nodes'] if node['id'] == src_id), None)
        tgt_node = next((node for node in inp_restrict['nodes'] if node['id'] == tgt_id), None)
        assert src_node, f'output(s) of source node of edge{edg} must exist!'
        assert tgt_node, f'input(s) of target node of edge{edg} must exist!'
        # flattened list of keys
        if src_node.get('out') and tgt_node.get('in'):
            src_out_keys = [sk for sout in src_node['out'] for sk in sout.keys()]
            tgt_in_keys = tgt_node['in'].keys()
            # we match the source output tag type to target input tag type
            # and connect them through '!* ' for input, all outputs are '!& ' before this
            for sk in src_out_keys:
                tgt_node['in'][sk] = yaml.load('!* ' + tgt_node['in'][sk], Loader=wic_loader())
            # the inputs which aren't dependent on previous/other steps
            # they are by default inline input
            diff_keys = set(tgt_in_keys) - set(src_out_keys)
            for dfk in diff_keys:
                tgt_node['in'][dfk] = yaml.load('!ii ' + tgt_node['in'][dfk], Loader=wic_loader())

    for node in inp_restrict['nodes']:
        node['id'] = node['name']  # just reuse name as node's id, wic id is same as wfb name
        node.pop('name', None)

    workflow_temp: Cwl = {}
    if inp_restrict["links"] != []:
        workflow_temp["steps"] = []
        for node in inp_restrict["nodes"]:
            workflow_temp["steps"].append(node)  # node["cwlScript"]  # Assume dict form
    else:  # A single node workflow
        node = inp_restrict["nodes"][0]
        workflow_temp = node["cwlScript"]
    return workflow_temp
