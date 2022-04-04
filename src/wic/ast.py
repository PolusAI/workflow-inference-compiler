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
def read_AST_from_disk(yaml_tree_tuple: YamlTree,
                       yml_paths: Dict[str, Path],
                       tools: Tools) -> YamlTree:
    (yaml_name, yaml_tree) = yaml_tree_tuple
    
    if 'backends' in yaml_tree:
        # Recursively expand each backend, but do NOT choose a specific backend.
        # Require back_name to be .yml? For now, yes.
        backends_trees_dict = dict([read_AST_from_disk((back_name, back), yml_paths, tools) for back_name, back in yaml_tree['backends'].items()])
        yaml_tree['backends'] = backends_trees_dict
        return (yaml_name, yaml_tree)

    steps: List[Yaml] = yaml_tree['steps']
    steps_keys = utils.get_steps_keys(steps)
    subkeys = [key for key in steps_keys if key not in tools]

    for i, step_key in enumerate(steps_keys):
        stem = Path(step_key).stem

        # Recursively read subworkflows, adding yml file contents
        if step_key in subkeys:
            yaml_path = yml_paths[stem]
            if not (yaml_path.exists() and yaml_path.suffix == '.yml'):
                raise Exception(f'Error! {yaml_path} does not exist or is not a .yml file.')

            # Load the high-level yaml sub workflow file.
            with open(yaml_path, 'r') as y:
                sub_yaml_tree_raw: Yaml = yaml.safe_load(y.read())
                # TODO: Once we have defined a yml DSL schema,
                # check that the file contents actually satisfies the schema.

            (step_key_, sub_yml_tree) = read_AST_from_disk((step_key, sub_yaml_tree_raw), yml_paths, tools)

            # inine sub_yml_tree; Since we are using a top-level wic: key,
            # there shouldn't be any key collisions, but we should check.
            step_i_dict = {} if steps[i][step_key] is None else steps[i][step_key]
            steps[i][step_key] = {**sub_yml_tree, **step_i_dict}

    return (yaml_name, yaml_tree)


def merge_yml_trees(yaml_tree_tuple: YamlTree,
                    wic_parent: Yaml,
                    tools: Tools) -> YamlTree:
    (yaml_name, yaml_tree) = yaml_tree_tuple

    # Check for top-level yml dsl args
    wic_self = {'wic': yaml_tree.get('wic', {})}
    wic = merge(wic_self, wic_parent, strategy=Strategy.TYPESAFE_REPLACE) # TYPESAFE_ADDITIVE ? 
    # Here we want to ADD wic: as a top-level yaml tag.
    # In the compilation phase, we want to remove it.
    yaml_tree['wic'] = wic['wic']
    wic_steps = wic['wic'].get('steps', {})
    
    if 'backends' in yaml_tree:
        # Recursively expand each backend, but do NOT choose a specific backend.
        # Require back_name to be .yml? For now, yes.
        # Pass wic_parent through unmodified (i.e. For now, simply skip backend: steps.)
        backends_trees_dict = dict([merge_yml_trees((back_name, back), wic_parent, tools) for back_name, back in yaml_tree['backends'].items()])
        yaml_tree['backends'] = backends_trees_dict
        return (yaml_name, yaml_tree)

    steps: List[Yaml] = yaml_tree['steps']
    steps_keys = utils.get_steps_keys(steps)
    subkeys = [key for key in steps_keys if key not in tools]

    for i, step_key in enumerate(steps_keys):
        # Recursively merge subworkflows, to implement parameter passing.
        if step_key in subkeys:
            # Extract the sub yaml file that we pre-loaded from disk.
            sub_yml_tree_initial = steps[i][step_key]
            sub_wic = wic_steps.get(f'({i+1}, {step_key})', {})

            (step_key_, sub_yml_tree) = merge_yml_trees((step_key, sub_yml_tree_initial), sub_wic, tools)
            # Now mutably overwrite the self args with the merged args
            steps[i][step_key] = sub_yml_tree

        # Extract provided CWL args, if any, and (recursively) merge them with
        # provided CWL args passed in from the parent, if any.
        # (At this point, any DSL args provided from the parent(s) should have
        # all of the initial yml tags removed, leaving only CWL tags remaining.)
        if not (step_key in subkeys):
            clt_args = wic_steps.get(f'({i+1}, {step_key})', {})
            if 'wic' in clt_args: # Do NOT add yml tags to the raw CWL!
                raise Exception(f'Error! wic: key found in CommandLineTool args.\n' + yaml.dump(sub_wic))
            sub_yml_tree = clt_args
            args_provided_dict_self = {}
            if steps[i][step_key]:
                args_provided_dict_self = steps[i][step_key]
            # NOTE: To support overloading, the parent args must overwrite the child args!
            args_provided_dict = merge(args_provided_dict_self, sub_yml_tree,
                                    strategy=Strategy.TYPESAFE_REPLACE) # TYPESAFE_ADDITIVE ?
            # Now mutably overwrite the self args with the merged args
            steps[i][step_key] = args_provided_dict

    return (yaml_name, yaml_tree)


def tree_to_forest(yaml_tree_tuple: YamlTree, tools: Tools) -> YamlForest:
    (yaml_name, yaml_tree) = yaml_tree_tuple
    
    if 'backends' in yaml_tree:
        backends_forest_dict = dict([(back_name, tree_to_forest((back_name, back), tools)) for back_name, back in yaml_tree['backends'].items()])
        return ((yaml_name, yaml_tree), backends_forest_dict)

    steps: List[Yaml] = yaml_tree['steps']
    steps_keys = utils.get_steps_keys(steps)
    subkeys = [key for key in steps_keys if key not in tools]

    yaml_forest_dict = {}

    for i, step_key in enumerate(steps_keys):

        if step_key in subkeys:
            sub_yaml_tree = steps[i][step_key]
            sub_yml_forest = tree_to_forest((step_key, sub_yaml_tree), tools)
            (sub_yml_tree_name, sub_yml_tree_) = sub_yml_forest[0]
            yaml_forest_dict[sub_yml_tree_name] = sub_yml_forest

    return ((yaml_name, yaml_tree), yaml_forest_dict)