from typing import Any, Dict

from .wic_types import Json, Tool, Tools


def cwl_tool_schema(name: str, tool: Tool) -> Json:
    inputs_props: Dict[Any, Any] = {}
    #required = []
    for key, val in tool.cwl['inputs'].items():
        inputs_props[key] = {}
        # Initialize special cases
        if name == 'mdrun' and key == 'config':
            inputs_props[key] = gromacs_mdp_schema()

        inputs_props[key]['title'] = val.get('label', '')
        inputs_props[key]['description'] = val.get('doc', '')

        valtype = val.get('type', '')
        # Determine required keys
        #if key == 'config' or not ('?' in valtype or 'default' in val):
        #    required.append(key)

        # Add type information, with exceptions
        valtype = valtype.replace('?', '')
        if (not (valtype == '' or valtype == 'File') and # Json does not have a File type
            not (name == 'mdrun' and key == 'config')):  # Exclude mdp
            inputs_props[key]['type'] = valtype

    # Do not mark properties which are required for CWL as required for yml,
    # because the whole point of inference is that we shouldn't have to!
    #if not required == []:
    #    inputs_props['required'] = required

    inputs: Dict[Any, Any] = {}
    inputs['type'] = 'object'
    inputs['additionalProperties'] = False
    inputs['properties'] = inputs_props

    step_name_props = {'in': inputs}

    step_name: Dict[Any, Any] = {}
    step_name['type'] = 'object'
    step_name['additionalProperties'] = False
    step_name['properties'] = step_name_props

    name_tag = {name: step_name}

    schema: Dict[Any, Any] = {}
    schema['$schema'] = 'https://json-schema.org/draft/2020-12/schema'
    schema['$id'] = name
    schema['type'] = 'object'
    schema['additionalProperties'] = False
    schema['title'] = tool.cwl.get('label', '')
    schema['description'] = tool.cwl.get('doc', '')
    schema['properties'] = name_tag
    return schema


def wic_tag_schema() -> Json:
    # NOTE: This schema needs to be recursive. Use dynamic anchors / references.
    # See https://json-schema.org/draft/2020-12/json-schema-core.html#dynamic-ref
    # and https://stackoverflow.com/questions/69728686/explanation-of-dynamicref-dynamicanchor-in-json-schema-as-opposed-to-ref-and

    graphviz_props: Dict[Any, Any] = {}
    graphviz_props['label'] = {'type': 'string'}
    graphviz_props['style'] = {'type': 'string'}
    graphviz_props['ranksame'] = {'type': 'array'}
    graphviz_props['ranksame']['items'] = {'type': 'string'}

    graphviz: Dict[Any, Any] = {}
    graphviz['type'] = 'object'
    graphviz['additionalProperties'] = False
    graphviz['properties'] = graphviz_props

    steps: Dict[Any, Any] = {}
    steps['type'] = 'object'
    steps['additionalProperties'] = False # This still works with patternProperties FYI

    # Call recursive reference
    recursive_ref = {"$dynamicRef": "#wic"}
    in_props: Dict[Any, Any] = {} # TODO: Add yml specific properties

    choices: Dict[Any, Any] = {}
    choices['type'] = 'object'
    choices['additionalProperties'] = False
    choices['properties'] = {'in': in_props, 'wic': recursive_ref}

    # See https://json-schema.org/understanding-json-schema/reference/object.html#patternproperties
    # NOTE: This recursive schema is correct, as determined by jsonschema.validate()
    # However, it seems that the vscode YAML extension does not support recursive
    # schema. (IntelliSense works fine until the first instance of recursion.)
    # TODO: A workaround would be to autogenerate a specific schema for each
    # yml file. We should probably do this anyway for the in: tag.
    steps['patternProperties'] = {"\\([0-9]+, [A-Za-z_\\.]+\\)": choices}

    schema_props: Dict[Any, Any] = {}
    schema_props['type'] = 'object'
    schema_props['additionalProperties'] = False
    schema_props['properties'] = {'graphviz': graphviz, 'steps': steps}

    schema: Dict[Any, Any] = {}
    schema['$schema'] = 'https://json-schema.org/draft/2020-12/schema'
    schema['$id'] = 'wic_tag'
    # Create recursive anchor
    schema['$dynamicAnchor'] = 'wic'
    schema['type'] = 'object'
    schema['additionalProperties'] = False
    schema['title'] = 'Metadata annotations'
    schema['description'] = 'Use steps: to recursively overload / pass parameters.\nUse graphviz: to modify the DAGs.'
    schema['properties'] = {'graphviz': graphviz, 'steps': steps}
    return schema


def wic_main_schema(tools_cwl: Tools) -> Json:
    tools_refs = [{'$ref': f'tools/{tool}.json'} for tool in tools_cwl]

    any_tool: Dict[Any, Any] = {}
    any_tool['type'] = 'object'
    any_tool['additionalProperties'] = False
    any_tool['properties'] = {'anyOf': tools_refs, 'title': 'Valid workflow steps'}

    steps: Dict[Any, Any] = {}
    steps['type'] = 'array'
    steps['description'] = 'A list of workflow steps'
    steps['items'] = {'anyOf': tools_refs, 'title': 'Valid workflow steps'}

    schema: Dict[Any, Any] = {}
    schema['$schema'] = 'https://json-schema.org/draft/2020-12/schema'
    schema['$id'] = 'wic_main'
    schema['type'] = 'object'
    schema['additionalProperties'] = False
    schema['title'] = 'Validating against the Workflow Interence Compiler schema'
    #schema['description'] = ''
    #schema['required'] = ['steps']
    schema['properties'] = {'wic': wic_tag_schema(), 'steps': steps} # 'required': ['steps']

    return schema


def gromacs_mdp_schema() -> Json:
    # See https://manual.gromacs.org/documentation/current/user-guide/mdp-options.html
    # TODO: Determine if there is a formal specification of the mdp options.
    # If not, enter the information by hand from the manual.
    mdpoptions = {}
    mdpoptions['integrator'] = {'type': 'string'} # TODO: 'oneOf': ['sd', ...]
    mdpoptions['rvdw'] = {'type': 'number'}

    mdp: Dict[Any, Any] = {}
    mdp['type'] = 'object'
    mdp['additionalProperties'] = False
    mdp['properties'] = mdpoptions

    schema_props = {'mdp': mdp, 'maxwarn': {'type': 'number', "minimum": 0}}

    schema: Dict[Any, Any] = {}
    schema['$schema'] = 'https://json-schema.org/draft/2020-12/schema'
    schema['$id'] = 'gromacs_mdp'
    schema['type'] = 'object'
    schema['additionalProperties'] = False
    schema['properties'] = schema_props
    return schema