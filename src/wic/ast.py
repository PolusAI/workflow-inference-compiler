import copy
from pathlib import Path
from typing import Dict, List

from mergedeep import merge, Strategy
import yaml

from . import utils
from .wic_types import Namespaces, Yaml, Tools, YamlTree, YamlForest

# NOTE: AST = Abstract Syntax Tree

# TODO: Check for inline-ing subworkflows more than once and, if there are not
# any modifications from any parent dsl args, use yaml anchors and aliases.
# That way, we should be able to serialize back to disk without duplication.
def read_AST_from_disk(yaml_tree_tuple: YamlTree,
                       yml_paths: Dict[str, Path],
                       tools: Tools) -> YamlTree:
    """Reads the yml workflow definition files from disk (recursively) and inlines them into an AST

    Args:
        yaml_tree_tuple (YamlTree): A tuple of a filepath and its Yaml file contents.
        yml_paths (Dict[str, Path]): The yml workflow definitions found using get_yml_paths()
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()

    Raises:
        Exception: If the yml file(s) do not exist

    Returns:
        YamlTree: A tuple of the root filepath and the associated yml AST
    """
    (yaml_name, yaml_tree) = yaml_tree_tuple

    wic = {'wic': yaml_tree.get('wic', {})}
    if 'backends' in wic['wic']:
        # Recursively expand each backend, but do NOT choose a specific backend.
        # Require back_name to be .yml? For now, yes.
        backends_trees_dict = dict([read_AST_from_disk(YamlTree(back_name, back), yml_paths, tools) for back_name, back in wic['wic']['backends'].items()])
        yaml_tree['wic']['backends'] = backends_trees_dict
        return YamlTree(yaml_name, yaml_tree)

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

            (step_key_, sub_yml_tree) = read_AST_from_disk(YamlTree(step_key, sub_yaml_tree_raw), yml_paths, tools)

            # inine sub_yml_tree; Since we are using a top-level wic: key,
            # there shouldn't be any key collisions, but we should check.
            step_i_dict = {} if steps[i][step_key] is None else steps[i][step_key]
            steps[i][step_key] = {**sub_yml_tree, **step_i_dict}

    return YamlTree(yaml_name, yaml_tree)


def merge_yml_trees(yaml_tree_tuple: YamlTree,
                    wic_parent: Yaml,
                    tools: Tools) -> YamlTree:
    """Implements 'parameter passing' by recursively merging wic: yml tags.
    Values from the parent workflow will overwrite / override subworkflows.
    See examples/gromacs/basic.yml for details

    Args:
        yaml_tree_tuple (YamlTree): A tuple of a name and a yml AST
        wic_parent (Yaml): The wic: yml dict from the parent workflow
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()

    Raises:
        Exception: If a wic: tag is found as an argument to a CWL CommandLineTool

    Returns:
        YamlTree: The yml AST with all wic: tags recursively merged.
    """
    (yaml_name, yaml_tree) = yaml_tree_tuple

    # Check for top-level yml dsl args
    wic_self = {'wic': yaml_tree.get('wic', {})}
    wic = merge(wic_self, wic_parent, strategy=Strategy.TYPESAFE_REPLACE) # TYPESAFE_ADDITIVE ? 
    # Here we want to ADD wic: as a top-level yaml tag.
    # In the compilation phase, we want to remove it.
    yaml_tree['wic'] = wic['wic']
    wic_steps = wic['wic'].get('steps', {})

    if 'backends' in wic['wic']:
        # Recursively expand each backend, but do NOT choose a specific backend.
        # Require back_name to be .yml? For now, yes.
        backends_trees_dict = dict([merge_yml_trees(YamlTree(back_name, back), wic_parent, tools) for back_name, back in wic['wic']['backends'].items()])
        yaml_tree['wic']['backends'] = backends_trees_dict
        return YamlTree(yaml_name, yaml_tree)

    steps: List[Yaml] = yaml_tree['steps']
    steps_keys = utils.get_steps_keys(steps)
    subkeys = [key for key in steps_keys if key not in tools]

    for i, step_key in enumerate(steps_keys):
        # Recursively merge subworkflows, to implement parameter passing.
        if step_key in subkeys:
            # Extract the sub yaml file that we pre-loaded from disk.
            sub_yml_tree_initial = steps[i][step_key]
            sub_wic = wic_steps.get(f'({i+1}, {step_key})', {})

            (step_key_, sub_yml_tree) = merge_yml_trees(YamlTree(step_key, sub_yml_tree_initial), sub_wic, tools)
            # Now mutably overwrite the self args with the merged args
            steps[i][step_key] = sub_yml_tree

        # Extract provided CWL args, if any, and (recursively) merge them with
        # provided CWL args passed in from the parent, if any.
        # (At this point, any DSL args provided from the parent(s) should have
        # all of the initial yml tags removed, leaving only CWL tags remaining.)
        if not (step_key in subkeys):
            clt_args = wic_steps.get(f'({i+1}, {step_key})', {})
            if 'wic' in clt_args:
                # Do NOT add yml tags to the raw CWL!
                # We can simply leave any step-specific wic: tags at top-level.
                # Copy so we only delete from the step, not also the top-level.
                clt_args = copy.deepcopy(clt_args)
                del clt_args['wic']
            sub_yml_tree = clt_args
            args_provided_dict_self = {}
            if steps[i][step_key]:
                args_provided_dict_self = steps[i][step_key]
            # NOTE: To support overloading, the parent args must overwrite the child args!
            args_provided_dict = merge(args_provided_dict_self, sub_yml_tree,
                                    strategy=Strategy.TYPESAFE_REPLACE) # TYPESAFE_ADDITIVE ?
            # Now mutably overwrite the self args with the merged args
            steps[i][step_key] = args_provided_dict

    return YamlTree(yaml_name, yaml_tree)


def tree_to_forest(yaml_tree_tuple: YamlTree, tools: Tools) -> YamlForest:
    """The purpose of this function is to abstract away the process of traversing an AST.

    Args:
        yaml_tree_tuple (YamlTree): A tuple of name and yml AST
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()

    Returns:
        YamlForest: A recursive data structure containing all sub-trees encountered while traversing the yml AST.
    """
    (yaml_name, yaml_tree) = yaml_tree_tuple

    wic = {'wic': yaml_tree.get('wic', {})}
    if 'backends' in wic['wic']:
        backends_forest_dict = dict([(back_name, tree_to_forest(YamlTree(back_name, back), tools)) for back_name, back in wic['wic']['backends'].items()])
        return YamlForest(YamlTree(yaml_name, yaml_tree), backends_forest_dict)

    steps: List[Yaml] = yaml_tree['steps']
    steps_keys = utils.get_steps_keys(steps)
    subkeys = [key for key in steps_keys if key not in tools]

    yaml_forest_dict = {}

    for i, step_key in enumerate(steps_keys):

        if step_key in subkeys:
            sub_yaml_tree = steps[i][step_key]
            sub_yml_forest = tree_to_forest(YamlTree(step_key, sub_yaml_tree), tools)
            (sub_yml_tree_name, sub_yml_tree_) = sub_yml_forest.yaml_tree
            yaml_forest_dict[sub_yml_tree_name] = sub_yml_forest

    return YamlForest(YamlTree(yaml_name, yaml_tree), yaml_forest_dict)


def get_inlineable_subworkflows(yaml_tree_tuple: YamlTree,
                                tools: Tools,
                                namespaces_init: Namespaces = []) -> List[Namespaces]:
    """Traverses a yml AST and finds all subworkflows which can be inlined into their parent workflow.

    Args:
        yaml_tree_tuple (YamlTree): A tuple of name and yml AST
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()
        namespaces_init (Namespaces): The initial subworkflow to start the traversal ([] == root)

    Returns:
        List[Namespaces]: The subworkflows which can be inlined into their parent workflows.
    """
    (yaml_name, yaml_tree) = yaml_tree_tuple

    # Check for top-level yml dsl args
    wic = {'wic': yaml_tree.get('wic', {})}

    if 'backends' in wic['wic']:
        # Use yaml_name (instead of back_name) and do not append to namespace_init.
        sub_namespaces_list = [get_inlineable_subworkflows(YamlTree(yaml_name, back), tools, namespaces_init) for back_name, back in wic['wic']['backends'].items()]
        return utils.flatten(sub_namespaces_list)

    steps: List[Yaml] = yaml_tree['steps']
    steps_keys = utils.get_steps_keys(steps)
    subkeys = [key for key in steps_keys if key not in tools]

    # All subworkflows except backends are inlineable.
    inlineable = wic.get('inlineable', True)
    namespaces = [namespaces_init] if inlineable and not (namespaces_init == []) else []

    for i, step_key in enumerate(steps_keys):
        yaml_stem = Path(yaml_name).stem
        step_name_i = utils.step_name_str(yaml_stem, i, step_key)
        if step_key in subkeys:
            sub_yml_tree = steps[i][step_key]

            sub_namespaces = get_inlineable_subworkflows(YamlTree(step_key, sub_yml_tree), tools, namespaces_init + [step_name_i])
            namespaces += sub_namespaces

    return namespaces


def inline_subworkflow(yaml_tree_tuple: YamlTree, tools: Tools, namespaces: Namespaces) -> YamlTree:
    """Inlines the given subworkflow into its immediate parent workflow.

    Args:
        yaml_tree_tuple (YamlTree): A tuple of name and yml AST
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()
        namespaces (Namespaces): Specifies the path in the yml AST to the subworkflow to be inlined.

    Returns:
        YamlTree: The updated root workflow with the given subworkflow inlined into its immediate parent workflow.
    """
    if namespaces == []:
        return yaml_tree_tuple

    (yaml_name, yaml_tree) = yaml_tree_tuple

    wic = {'wic': yaml_tree.get('wic', {})}
    if 'backends' in wic['wic']:
        # Pass namespaces through unmodified
        backends_trees_dict = dict([inline_subworkflow(YamlTree(back_name, back), tools, namespaces) for back_name, back in wic['wic']['backends'].items()])
        yaml_tree['wic']['backends'] = backends_trees_dict
        return YamlTree(yaml_name, yaml_tree)

    steps: List[Yaml] = yaml_tree['steps']
    steps_keys = utils.get_steps_keys(steps)
    subkeys = [key for key in steps_keys if key not in tools]

    for i, step_key in enumerate(steps_keys):
        yaml_stem = Path(yaml_name).stem
        step_name_i = utils.step_name_str(yaml_stem, i, step_key)
        if step_key in subkeys:
            sub_yml_tree = steps[i][step_key]

            if namespaces[0] == step_name_i:
                if len(namespaces) == 1:
                    steps_inits = steps[:i] # Exclude step i
                    steps_tails = steps[i+1:] # Exclude step i
                    # Inline sub-steps.
                    sub_steps: List[Yaml] = sub_yml_tree['steps']
                    yaml_tree['steps'] = steps_inits + sub_steps + steps_tails
                    # Need to re-index both the sub-step numbers as well as the
                    # subsequent steps in this workflow? No, except for wic: steps:
                else:
                    # Strip off one initial namespace
                    (step_key_, sub_yml_tree) = inline_subworkflow(YamlTree(step_key, sub_yml_tree), tools, namespaces[1:])
                    # TODO: re-index wic: steps: ? We probably should, although
                    # inlineing after merging should not affect CWL args.
                    # Re-indexing could be tricky w.r.t. overloading.
                    # TODO: maintain inference boundaries (once feature is added)
                    steps[i][step_key] = sub_yml_tree

    return YamlTree(yaml_name, yaml_tree)