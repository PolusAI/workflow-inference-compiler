import copy
from pathlib import Path
from typing import Dict, List

from mergedeep import merge, Strategy
import yaml

from . import utils
from .wic_types import Yaml, Tools, YamlTree, YamlForest

# TODO: Check for inline-ing subworkflows more than once and, if there are not
# any modifications from any parent dsl args, use yaml anchors and aliases.
# That way, we should be able to serialize back to disk without duplication.
def read_yml_forest(yaml_tree_tuple: YamlTree,
                    yaml_dsl_args: Yaml,
                    yml_paths: Dict[str, Path],
                    tools: Tools) -> YamlForest:
    (yaml_name, yaml_tree) = yaml_tree_tuple

    # TODO: Check for top-level yml dsl args
    # yaml_dsl_args_self = yaml_tree.get('wic', {})
    
    if 'backends' in yaml_tree:
        # Recursively expand each backend, but do NOT choose a specific backend.
        # Require back_name to be .yml?
        backends_forest_dict = dict([(back_name, read_yml_forest((back_name, back), {}, yml_paths, tools)) for back_name, back in yaml_tree['backends'].items()])
        backends = dict([f[0] for b, f in backends_forest_dict.items()])
        yaml_tree['backends'] = backends
        return ((yaml_name, yaml_tree), backends_forest_dict)

    steps: List[Yaml] = yaml_tree['steps']

    steps_keys = utils.get_steps_keys(steps)

    subkeys = [key for key in steps_keys if key not in tools]

    yaml_forest_dict = {}

    for i, step_key in enumerate(steps_keys):
        stem = Path(step_key).stem

        args_provided_dict_parent = {}
        
        # Recursively read subworkflows, adding yml file contents
        # (with any modifications due to parameter passing) to yaml_forest_dict
        if step_key in subkeys:
            yaml_path = yml_paths[stem]
            if not (yaml_path.exists() and yaml_path.suffix == '.yml'):
                raise Exception(f'Error! {yaml_path} does not exist or is not a .yml file.')

            # Load the high-level yaml sub workflow file.
            with open(yaml_path, 'r') as y:
                sub_yaml_tree: Yaml = yaml.safe_load(y.read())
                # TODO: Once we have defined a yml DSL schema,
                # check that the file contents actually satisfies the schema.

            # Extract provided yaml args, if any, and (recursively) merge them with
            # provided yaml_dsl_args passed in from the parent, if any.
            yaml_dsl_args_parent = {}
            if f'({i+1}, {step_key})' in yaml_dsl_args:
                yaml_dsl_args_parent = yaml_dsl_args[f'({i+1}, {step_key})']
            yaml_dsl_args_child = {}
            if steps[i][step_key]:
                yaml_dsl_args_child = copy.deepcopy(steps[i][step_key])
                # Delete child yml DSL parameters, if any, so that we are left
                # with only CWL parameters below.
                steps[i].update({step_key: {}}) # delete yml subtree
            # NOTE: To support overloading, the parent args must overwrite the child args!
            sub_yaml_dsl_args = merge(yaml_dsl_args_child, yaml_dsl_args_parent, strategy=Strategy.TYPESAFE_REPLACE) # TYPESAFE_ADDITIVE ? 
            #print('sub_yaml_dsl_args', sub_yaml_dsl_args)

            sub_yml_forest = read_yml_forest((step_key, sub_yaml_tree), sub_yaml_dsl_args, yml_paths, tools)
            (sub_yml_tree_name, sub_yml) = sub_yml_forest[0]
            yaml_forest_dict[sub_yml_tree_name] = sub_yml_forest
            args_provided_dict_parent = sub_yml

        # Extract provided CWL args, if any, and (recursively) merge them with
        # provided CWL args passed in from the parent, if any.
        # (At this point, any DSL args provided from the parent(s) should have
        # all of the initial yml tags removed, leaving only CWL tags remaining.)
        if (f'({i+1}, {step_key})' in yaml_dsl_args
            and (not step_key in subkeys)): # Do NOT add yml tags to the raw CWL!
            args_provided_dict_parent = yaml_dsl_args[f'({i+1}, {step_key})']
        args_provided_dict_child = {}
        if steps[i][step_key]:
            args_provided_dict_child = steps[i][step_key]
        # NOTE: To support overloading, the parent args must overwrite the child args!
        args_provided_dict = merge(args_provided_dict_child, args_provided_dict_parent,
                                   strategy=Strategy.REPLACE) # TYPESAFE_ADDITIVE ?
        # Now mutably overwrite the child args with the merged args
        steps[i][step_key] = args_provided_dict

    # Just return a single yaml tree?
    return ((yaml_name, yaml_tree), yaml_forest_dict)