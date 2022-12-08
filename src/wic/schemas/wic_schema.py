import argparse
import json
from pathlib import Path
import random
from unittest.mock import patch
import sys
from typing import Any, Dict, List

import networkx as nx
import graphviz
from jsonschema import RefResolver, Draft202012Validator
import yaml

import wic
from wic import ast, cli, compiler, utils_cwl
from wic.wic_types import GraphData, GraphReps, NodeData, StepId, Yaml, YamlTree
from ..wic_types import Json, Tools
from .biobb import config_schemas


def default_schema(url: bool = False) -> Json:
    """A basic default schema (to avoid copy & paste).

    Args:
        url (bool, optional): Determines whether to include the $schema url. Defaults to False.

    Returns:
        Json: A basic default schema
    """
    schema: Json = {}
    schema['type'] = 'object'
    schema['additionalProperties'] = False
    if url:
        schema['$schema'] = 'https://json-schema.org/draft/2020-12/schema'
    return schema


def named_empty_schema(name: str) -> Json:
    """Creates a schema which starts with name, but is otherwise an empty wildcard

    Args:
        name (str): The identifier of the string

    Returns:
        Json: A schema which matches anything starting with name
    """
    schema = default_schema()
    schema['properties'] = {name: {}} # NOTE: {} is essentially a wildcard
    return schema


def named_null_schema(name: str) -> Json:
    """Creates a schema which starts with name and contains nothing else

    Args:
        name (str): The identifier of the string

    Returns:
        Json: A schema which matches name and nothing else
    """
    # NOTE: Use this together with anyOf to allow no explicit arguments
    schema = default_schema()
    schema['properties'] = {name: {'type': 'null'}}
    return schema


def cwl_type_to_jsonschema_type_schema(type_obj: Json) -> Json:
    """Converts a canonicalized CWL type into the equivalent jsonschema type schema, if possible.

    Args:
        type_obj (Json): A canonical CWL type object

    Returns:
        Json: A JSON type schema corresponding to type_obj if valid else None
    """
    jsontype = cwl_type_to_jsonschema_type(type_obj)
    if jsontype is not None:
        if isinstance(jsontype, str) or isinstance(jsontype, List):
            return {'type': jsontype}
    return jsontype


def cwl_type_to_jsonschema_type(type_obj: Json) -> Json:
    """Converts a canonicalized CWL type into the equivalent jsonschema type schema, if possible.

    Args:
        type_obj (Json): A canonical CWL type object

    Returns:
        Json: A JSON type schema corresponding to type_obj if valid else None
    """
    if isinstance(type_obj, str):
        # Obviously ignore empty types
        if type_obj == '':
            return None
        # The 'null' type is valid; it does not indicate an error.
        if type_obj == 'null':
            return 'null'
        # There are only two numeric types in json
        if type_obj == 'int' or type_obj == 'long':
            return 'integer'
        if type_obj == 'float' or type_obj == 'double':
            return 'number'
        # Rename common abbreviations
        if type_obj == 'bool':
            return 'boolean'
        if type_obj == 'str':
            return 'string'
        # CWL suports an Any type; map this to the empty wildcard {} schema
        if type_obj == 'Any':
            return {}
        # json does not have any File or Directory types
        if type_obj == 'File':
            return None
        if type_obj == 'Directory':
            return None
        if type_obj == 'WritableFile':
            return None
        if type_obj == 'WritableDirectory':
            return None

    if isinstance(type_obj, Dict):
        if type_obj.get('type') == 'array' and 'items' in type_obj:
            items = cwl_type_to_jsonschema_type(type_obj['items'])
            if items is None:
                return None # Propagate any type failures
            if isinstance(type_obj['items'], str):
                # Wrap primitive strings in {'type': ...}
                return {**type_obj, 'items': {'type': items}}
            return {**type_obj, 'items': items}
        # TODO: Other cases?

    if isinstance(type_obj, List):
        items = [cwl_type_to_jsonschema_type(item) for item in type_obj]
        if any([item is None for item in items]):
            return None # Propagate any type failures
        # See https://cswr.github.io/JsonSchema/spec/multiple_types/
        # In a list, if some of the types are themselves arrays or objects,
        # we need to replace them with "array" and "object". This loses
        # information, but that's the specification, so...
        items = ['array' if isinstance(item, Dict) and item.get('type', '') == 'array' else item for item in items]
        items = ['object' if isinstance(item, Dict) and item.get('type', '') == 'object' else item for item in items]
        return items

    # TODO: Support CWL records
    return type_obj


def cwl_schema(name: str, cwl: Json, id_prefix: str) -> Json:
    """Generates a schema (including documentation) based on the inputs of a CWL CommandLineTool or Workflow.

    Args:
        name (str): The name of the CWL CommandLineTool or Workflow
        cwl (Json): The CWL CommandLineTool or Workflow
        id_prefix (str): Either the string 'tools' or 'workflows'

    Returns:
        Json: An autogenerated, documented schema based on the inputs and outputs of a CWL CommandLineTool or Workflow.
    """
    inputs_props: Json = {}
    #required = []
    for key, val in cwl['inputs'].items():
        metadata = {'title': val.get('label', ''), 'description': val.get('doc', '')}

        # Determine required keys
        #if key == 'config' or not ('?' in cwltype or 'default' in val):
        #    required.append(key)

        # Handle special cases
        if key == 'config' and name in config_schemas:
            inputs_props[key] = {'anyOf': [{'type': 'string', **metadata},
                                           {**config_schemas[name], **metadata}]}
            continue

        # Add type information, with exceptions
        cwltype = utils_cwl.canonicalize_type(val.get('type', ''))
        jsontype = cwl_type_to_jsonschema_type_schema(cwltype)
        if jsontype:
            inputs_props[key] = {'anyOf': [{'type': 'string', **metadata}, {**jsontype, **metadata}]}
        else:
            inputs_props[key] = {'type': 'string', **metadata}

    # Do not mark properties which are required for CWL as required for yml,
    # because the whole point of inference is that we shouldn't have to!
    #if not required == []:
    #    inputs_props['required'] = required

    inputs = default_schema()
    inputs['properties'] = inputs_props

    scatter_props = {'type': 'array', 'items': {'anyOf': [{**val, 'const': key} for key, val in inputs_props.items()]}}
    scatterMethod_props: Json = {'type': 'string', 'enum': ['dotproduct', 'flat_crossproduct', 'nested_crossproduct']}

    outputs_props: Json = {}
    for key, val in cwl['outputs'].items():
        metadata = {'title': val.get('label', ''), 'description': val.get('doc', '')}

        # Add type information, with exceptions
        cwltype = utils_cwl.canonicalize_type(val.get('type', ''))
        jsontype = cwl_type_to_jsonschema_type_schema(cwltype)
        if jsontype:
            outputs_props[key] = {'anyOf': [{'type': 'string', **metadata}, {**jsontype, **metadata}]}
        else:
            outputs_props[key] = {'type': 'string', **metadata}

    outputs = default_schema()
    outputs['properties'] = outputs_props

    step_props = default_schema()
    step_props['title'] = cwl.get('label', '')
    step_props['description'] = cwl.get('doc', '')
    step_props['properties'] = {'in': inputs,
                                'out': outputs,
                                'scatter': scatter_props,
                                'scatterMethod': scatterMethod_props}

    schema = default_schema(url=True)
    # NOTE: See comment in get_validator(). Nonetheless, the vscode YAML extension
    # appears to be resolving ids w.r.t. relative local paths. jsonschema
    # (correctly) treats f'tools/{name}.json' as as uninterpreted string,
    # so instead of using name let's just use fake relative paths in ids.
    schema['$id'] = f'{id_prefix}/{name}.json'
    step_name = name + '.yml' if id_prefix == 'workflows' else name
    schema['properties'] = {step_name: step_props}
    return schema


def wic_tag_schema(hypothesis: bool = False) -> Json:
    """The schema of the (recursive) wic: metadata annotation tag.

    Args:
        hypothesis (bool): Determines whether we should restrict the search space.

    Returns:
        Json: The schema of the (recursive) wic: metadata annotation tag.
    """
    # NOTE: This schema needs to be recursive. Use dynamic anchors / references.
    # See https://json-schema.org/draft/2020-12/json-schema-core.html#dynamic-ref
    # and https://stackoverflow.com/questions/69728686/explanation-of-dynamicref-dynamicanchor-in-json-schema-as-opposed-to-ref-and

    graphviz_props: Json = {}
    graphviz_props['label'] = {'type': 'string'}
    graphviz_props['style'] = {'type': 'string'}
    graphviz_props['ranksame'] = {'type': 'array'}
    graphviz_props['ranksame']['items'] = {'type': 'string'}

    graphviz = default_schema()
    graphviz['properties'] = graphviz_props

    # Call recursive reference
    recursive_ref = {'$dynamicRef': '#wic'}
    in_props: Json = {} # TODO: Add yml specific properties

    scatter_props: Json = {} # TODO: Add yml specific properties
    scatterMethod_props: Json = {'type': 'string', 'enum': ['dotproduct', 'flat_crossproduct', 'nested_crossproduct']}

    choices_props = {'wic': recursive_ref, 'scatterMethod': scatterMethod_props}
    if not hypothesis:
        # Empty wildcard {} schemas can cause problems with hypothesis.
        choices_props['in'] = in_props
        choices_props['scatter'] = scatter_props
    choices = default_schema()
    choices['properties'] = choices_props

    # See https://json-schema.org/understanding-json-schema/reference/object.html#patternproperties
    # NOTE: This recursive schema is correct, as determined by jsonschema.validate()
    # However, it seems that the vscode YAML extension does not support recursive
    # schema. (IntelliSense works fine until the first instance of recursion.)
    # TODO: A workaround would be to autogenerate a specific schema for each
    # yml file. We should probably do this anyway for the in: tag.
    steps = default_schema()
    # additionalProperties = False still works with patternProperties FYI
    steps['patternProperties'] = {"\\([0-9]+, [A-Za-z0-9_\\.]+\\)": choices}

    #backends = default_schema()
    backends: Dict[Any, Any] = {}
    backends['type'] = 'object'
    backends['additionalProperties'] = True
    # TODO: Restrict the backend properties and make default_backend an enum

    namespace: Dict[Any, Any] = {}
    namespace['type'] = 'string'
    # namespace['enum'] = ...
    # TODO: Restrict the namespace properties to only those in yml_paths.txt

    backend = {'type': 'string'}
    default_backend = {'type': 'string'}
    inlineable = {'type': 'boolean'}

    environment_props: Json = {}
    environment_props['action'] = {'type': 'string'}
    environment_props['save_defs'] = {'type': 'array'}
    environment_props['save_defs']['items'] = {'type': 'string'}

    environment = default_schema()
    environment['properties'] = environment_props

    schema = default_schema(url=True)
    schema['$id'] = 'wic_tag'
    # Create recursive anchor
    schema['$dynamicAnchor'] = 'wic'
    schema['title'] = 'Metadata annotations'
    schema['description'] = 'Use steps: to recursively overload / pass parameters.\nUse graphviz: to modify the DAGs.'

    schema_props = {'graphviz': graphviz, 'steps': steps, 'backend': backend,
                    'default_backend': default_backend,
                    'namespace': namespace, 'inlineable': inlineable, 'environment': environment}
    if not hypothesis:
        # {'additionalProperties': True} can cause problems with hypothesis.
        schema_props['backends'] = backends
    schema['properties'] = schema_props
    return schema


def wic_main_schema(tools_cwl: Tools, yml_stems: List[str], schema_store: Dict[str, Json], hypothesis: bool = False) -> Json:
    """The main schema which is used to validate yml files.

    Args:
        tools_cwl (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()
        yml_stems (List[str]): The names of the yml workflow definitions found using get_yml_paths()
        schema_store (Dict[str, Json]): A global mapping between ids and schemas
        hypothesis (bool): Determines whether we should restrict the search space.

    Returns:
        Json: The main schema which is used to validate yml files.
    """
    # NOTE: As mentioned below, using $ref's with external schema files
    # (coincidentally?) works with the VSCode YAML extension, and for the
    # jsonschema library we can supply an explicit schemastore. The API of the
    # hypothesis-jsonschema library, however, only takes a schema. So we either
    # need to bundle the external file contents into wic.json (using $def's),
    # or (since there is only one call site per file) simply inline the contents.
    tools_schemas: List[Json] = [{'anyOf': [schema_store.get(f'tools/{step_id.stem}.json', {'$ref': f'tools/{step_id.stem}.json'}),
                                            named_null_schema(step_id.stem)]} for step_id in tools_cwl if not step_id.stem == 'python_script']
#    tools_schemas: List[Json] = [{'anyOf': [{'$ref': f'tools/{step_id.stem}.json'},
#                                            named_null_schema(step_id.stem)]} for step_id in tools_cwl]
    # NOTE: See comment in get_validator(). Nonetheless, the vscode YAML extension
    # appears to be resolving ids w.r.t. relative local paths. jsonschema
    # (correctly) treats f'tools/{name}.json' as an uninterpreted string,
    # so instead of using stem let's just use fake relative paths in ids.

    # NOTE: We could/should re-validate after every AST modification. This will
    # require substantial code changes, so let's not worry about it for now.
    yml_schemas: List[Json] = [{'anyOf': [schema_store.get(f'workflows/{yml_stem}.json', {'$ref': f'workflows/{yml_stem}.json'}),
                                          named_null_schema(f'{yml_stem}.yml')]} for yml_stem in yml_stems]
#    yml_schemas: List[Json] = [{'anyOf': [{'$ref': f'workflows/{yml_stem}.json'},
#                                          named_null_schema(f'{yml_stem}.yml')]} for yml_stem in yml_stems]

    steps: Json = {}
    steps['type'] = 'array'
    steps['description'] = 'A list of workflow steps'

    if hypothesis:
        # For performance reasons, limit the size of the schema. The first time
        # you call .example(), hypothesis will compile the schema and cache the
        # results. Subsequent .example() calls are nearly instantaneous.
        # The time increases fairly rapidly with k, i.e.
        k = 3 # 1-5 minutes...
        # Choose a random subset so we're not testing the same files
        tools_schemas = random.choices(tools_schemas, k=k)
        yml_schemas = random.choices(yml_schemas, k=k)

    steps_schemas = tools_schemas + yml_schemas
    if not hypothesis:
        steps_schemas += [named_empty_schema('python_script')]

    steps['items'] = {'anyOf': steps_schemas, 'minItems': 1, 'title': 'Valid workflow steps'}

    # TODO: Use the real CWL inputs schema
    inputs: Dict[Any, Any] = {}
    inputs['type'] = 'object'
    inputs['additionalProperties'] = True

    # TODO: Use the real CWL outputs schema
    outputs: Dict[Any, Any] = {}
    outputs['type'] = 'object'
    outputs['additionalProperties'] = True

    schema = default_schema(url=True)
    schema['$id'] = 'wic_main'
    schema['title'] = 'Validating against the Workflow Interence Compiler schema'
    #schema['description'] = ''
    #schema['required'] = ['steps'] # steps are not required, e.g. npt.yml
    schema_props = {'steps': steps, 'label': {'type': 'string'}, 'doc': {'type': 'string'}}
    #schema_props['wic'] = wic_tag_schema() # NOTE: This technicaly 'works'
    # with hypothesis, but the wic_tag_schema still needs some work.
    if not hypothesis:
        schema_props['wic'] = wic_tag_schema()
        # {'additionalProperties': True} can cause problems with hypothesis.
        schema_props['inputs'] = inputs
        schema_props['outputs'] = outputs
    schema['properties'] = schema_props

    # https://json-schema.org/understanding-json-schema/structuring.html#bundling
    #import copy
    #schema['$defs'] = copy.deepcopy(schema_store)
    # Without deepcopy, "ValueError: Circular reference detected"
    # "f.write(json.dumps(schema, indent=2))"

    return schema


def get_args(yml_path: str = '') -> argparse.Namespace:
    """This is used to get mock command line arguments.

    Returns:
        argparse.Namespace: The mocked command line arguments
    """
    testargs = ['wic', '--yaml', yml_path, '--cwl_output_intermediate_files', 'True']  # ignore --yaml
    # For now, we need to enable --cwl_output_intermediate_files. See comment in compiler.py
    with patch.object(sys, 'argv', testargs):
        args: argparse.Namespace = wic.cli.parser.parse_args()
    return args


def compile_workflow_generate_schema(yml_path_str: str, yml_path: Path,
                                     tools_cwl: Tools,
                                     yml_paths: Dict[str, Dict[str, Path]],
                                     validator: Draft202012Validator) -> Json:
    """Compiles a workflow and generates a schema which (recursively) includes the inputs/outputs from subworkflows.

    Args:
        yml_path_str (str): The stem of the path to the yml file
        yml_path (Path): The path to the yml file
        tools_cwl (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()
        yml_paths (Dict[str, Dict[str, Path]]): The yml workflow definitions found using get_yml_paths()
        validator (Draft202012Validator): Used to validate the yml files against the autogenerated schema.

    Returns:
        Json: An autogenerated, documented schema based on the inputs and outputs of the Workflow.
    """
    # First compile the workflow.
    # Load the high-level yaml workflow file.
    with open(yml_path, mode='r', encoding='utf-8') as y:
        root_yaml_tree: Yaml = yaml.safe_load(y.read())
    Path('autogenerated/').mkdir(parents=True, exist_ok=True)
    wic_tag = {'wic': root_yaml_tree.get('wic', {})}
    plugin_ns = wic_tag['wic'].get('namespace', 'global')
    step_id = StepId(yml_path_str, plugin_ns)
    y_t = YamlTree(step_id, root_yaml_tree)
    yaml_tree_raw = wic.ast.read_ast_from_disk(y_t, yml_paths, tools_cwl, validator)
    #with open(f'autogenerated/{Path(yml_path).stem}_tree_raw.yml', mode='w', encoding='utf-8') as f:
    #    f.write(yaml.dump(yaml_tree_raw.yml))
    yaml_tree = wic.ast.merge_yml_trees(yaml_tree_raw, {}, tools_cwl)
    #with open(f'autogenerated/{Path(yml_path).stem}_tree_merged.yml', mode='w', encoding='utf-8') as f:
    #    f.write(yaml.dump(yaml_tree.yml))

    graph_gv = graphviz.Digraph(name=f'cluster_{yml_path}')
    graph_gv.attr(newrank='True')
    graph_nx = nx.DiGraph()
    graphdata = GraphData(str(yml_path))
    graph = GraphReps(graph_gv, graph_nx, graphdata)
    compiler_info = wic.compiler.compile_workflow(yaml_tree, get_args(str(yml_path)), [], [graph], {}, {}, {}, {},
                                                    tools_cwl, True, relative_run_path=True, testing=True)
    rose_tree = compiler_info.rose
    sub_node_data: NodeData = rose_tree.data

    #wic.utils.write_to_disk(rose_tree, Path('autogenerated/'), relative_run_path=True)
    schema = cwl_schema(step_id.stem, sub_node_data.compiled_cwl, 'workflows')

    #with open(f'autogenerated/schemas/workflows/{step_id.stem}.json', mode='w', encoding='utf-8') as f:
    #    f.write(json.dumps(schema, indent=2))

    return schema


def get_validator(tools_cwl: Tools, yml_stems: List[str], schema_store: Dict[str, Json] = {},
                  write_to_disk: bool = False, hypothesis: bool = False) -> Draft202012Validator:
    """Generates the main schema used to check the yml files for correctness and returns a validator.

    Args:
        tools_cwl (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()
        yml_stems (List[str]): The names of the yml workflow definitions found using get_yml_paths()
        schema_store (Dict[str, Json]): A global mapping between ids and schemas
        write_to_disk (bool): Controls whether to write the schemas to disk.
        hypothesis (bool): Determines whether we should restrict the search space.

    Returns:
        Draft202012Validator: A validator which is used to check the yml files for correctness.
    """
    for step_id, tool in tools_cwl.items():
        schema_tool = cwl_schema(step_id.stem, tool.cwl, 'tools')
        schema_store[schema_tool['$id']] = schema_tool
        #if write_to_disk:
        #    with open(f'autogenerated/schemas/tools/{step_id.stem}.json', mode='w', encoding='utf-8') as f:
        #        f.write(json.dumps(schema_tool, indent=2))

    # Add temporary placeholders to the schema_store so we don't get
    # "jsonschema.exceptions.RefResolutionError: unknown url type: 'workflows/*.json'"
    for yml_stem in yml_stems:
        if f'workflows/{yml_stem}.json' not in schema_store:
            schema_store[f'workflows/{yml_stem}.json'] = {}

    schema = wic_main_schema(tools_cwl, yml_stems, schema_store, hypothesis)
    schema_store[schema['$id']] = schema
    schema_store['wic_tag'] = wic_tag_schema(hypothesis)
    if write_to_disk:
        with open('autogenerated/schemas/wic.json', mode='w', encoding='utf-8') as f:
            f.write(json.dumps(schema, indent=2))

    # See https://stackoverflow.com/questions/53968770/how-to-set-up-local-file-references-in-python-jsonschema-document
    # The $ref tag refers to URIs defined in $id tags, NOT relative paths on
    # the local filesystem! We need to create a global mapping between ids and schemas
    # i.e. schema_store.
    resolver = RefResolver.from_schema(schema, store=schema_store)
    """ Use check_schema to 'first verify that the provided schema is
    itself valid, since not doing so can lead to less obvious error
    messages and fail in less obvious or consistent ways.'
    """
    # i.e. This should match 'https://json-schema.org/draft/2020-12/schema'
    # NOTE: If you get nasty errors while developing the schema such as:
    # "jsonschema.exceptions.SchemaError: ... is not valid under any of the given schemas"
    # try temporarily commmenting this line out to generate the schema anyway.
    # Then, in any yml file, the very first line should show a "schema stack trace"
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema, resolver=resolver)
    return validator
