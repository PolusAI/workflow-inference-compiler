import copy
import json
from pathlib import Path
import asyncio

from fastapi import Request

import pytest
from sophios.wic_types import Json


from sophios.api.http import restapi


@pytest.mark.fast
def test_rest_core_single_node() -> None:
    """A simple single node sophios/restapi test"""
    inp_file = "single_node.json"
    inp: Json = {}
    inp_path = Path(__file__).parent / 'rest_wfb_objects' / inp_file
    with open(inp_path, 'r', encoding='utf-8') as f:
        inp = json.load(f)
    print('----------- from rest api ----------- \n\n')
    scope = {}
    scope['type'] = 'http'

    async def receive() -> Json:
        inp_req: Json = {}
        inp_req['payload'] = copy.deepcopy(inp)
        inp_req['run'] = 'run'
        inp_byte = json.dumps(inp_req).encode('utf-8')
        return {"type": "http.request", "body": inp_byte}

    # create a request object and pack it with our json payload
    req: Request = Request(scope)
    req._receive = receive
    res: Json = asyncio.run(restapi.compile_wf(req))  # call to rest api
    assert int(res['retval']) == 0


@pytest.mark.fast
def test_rest_core_multi_node() -> None:
    """A simple multi node sophios/restapi test"""
    inp_file = "multi_node.json"
    inp: Json = {}
    inp_path = Path(__file__).parent / 'rest_wfb_objects' / inp_file
    with open(inp_path, 'r', encoding='utf-8') as f:
        inp = json.load(f)
    print('----------- from rest api ----------- \n\n')
    scope = {}
    scope['type'] = 'http'

    async def receive() -> Json:
        inp_req: Json = {}
        inp_req['payload'] = copy.deepcopy(inp)
        inp_req['run'] = 'run'
        inp_byte = json.dumps(inp_req).encode('utf-8')
        return {"type": "http.request", "body": inp_byte}

    # create a request object and pack it with our json payload
    req: Request = Request(scope)
    req._receive = receive
    res: Json = asyncio.run(restapi.compile_wf(req))  # call to rest api
    assert int(res['retval']) == 0


@pytest.mark.fast
def test_rest_core_multi_node_inline_cwl() -> None:
    """A simple multi node (inline cwl) sophios/restapi test"""
    inp_file = "multi_node_inline_cwl.json"
    inp: Json = {}
    inp_path = Path(__file__).parent / 'rest_wfb_objects' / inp_file
    with open(inp_path, 'r', encoding='utf-8') as f:
        inp = json.load(f)
    print('----------- from rest api ----------- \n\n')
    scope = {}
    scope['type'] = 'http'

    async def receive() -> Json:
        inp_req: Json = {}
        inp_req['payload'] = copy.deepcopy(inp)
        inp_req['run'] = 'run'
        inp_byte = json.dumps(inp_req).encode('utf-8')
        return {"type": "http.request", "body": inp_byte}

    # create a request object and pack it with our json payload
    req: Request = Request(scope)
    req._receive = receive
    res: Json = asyncio.run(restapi.compile_wf(req))  # call to rest api
    assert int(res['retval']) == 0
