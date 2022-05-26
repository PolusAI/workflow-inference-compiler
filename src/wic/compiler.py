import argparse
import copy
from pathlib import Path
import subprocess as sub
from typing import Dict, List

import json
import graphviz
import networkx  as nx

from . import inference
from . import utils
from .wic_types import Namespaces, Yaml, Tool, Tools, ExplicitEdgeDefs, ExplicitEdgeCalls, CompilerInfo, WorkflowInputsFile, WorkflowOutputs, InternalOutputs, GraphReps, YamlTree, RoseTree, NodeData, EnvData, GraphData

# Use white for dark backgrounds, black for light backgrounds
font_edge_color = 'white'

# NOTE: This must be initialized in main.py and/or cwl_watcher.py
inference_rules: Dict[str, str] = None # type: ignore


def compile_workflow(yaml_tree_ast: YamlTree,
                     args: argparse.Namespace,
                     namespaces: Namespaces,
                     subgraphs_: List[GraphReps],
                     explicit_edge_defs: ExplicitEdgeDefs,
                     explicit_edge_calls: ExplicitEdgeCalls,
                     tools: Tools,
                     is_root: bool,
                     relative_run_path: bool) -> CompilerInfo:
    """fixed-point wrapper around compile_workflow_once\n
    See https://en.wikipedia.org/wiki/Fixed_point_(mathematics)

    Args:
        yaml_tree_ast (YamlTree): A tuple of name and yml AST
        args (Any): all of the other positional arguments for compile_workflow_once
        kwargs (Any): all of the other keyword arguments for compile_workflow_once

    Returns:
        CompilerInfo: Contains the data associated with compiled subworkflows (in the Rose Tree) together with mutable cumulative environment information which needs to be passed through the recursion.
    """
    ast_modified = True
    yaml_tree = yaml_tree_ast
    # There ought to be at most one file format conversion between each step.
    # If everything is working correctly, we should thus reach the fixed point
    # in at most n-1 iterations. However, due to the possibility of bugs in the
    # implementation and/or spurious inputs, we should guarantee termination.
    max_iters = 100 # 100 ought to be plenty. TODO: calculate n-1 from steps:
    iter = 0
    while ast_modified and iter < max_iters:
        subgraphs = copy.deepcopy(subgraphs_) # See comment below!
        compiler_info = compile_workflow_once(yaml_tree, args, namespaces, subgraphs, explicit_edge_defs, explicit_edge_calls, tools, is_root, relative_run_path)
        node_data: NodeData = compiler_info.rose.data
        ast_modified = not (yaml_tree.yml == node_data.yml)
        if ast_modified:
            #import yaml
            #print(yaml.dump(node_data.yml))
            #print()
            yaml_tree = YamlTree(yaml_tree_ast.name, node_data.yml)
        iter += 1

    # Overwrite subgraphs_ element-wise
    # This is a terrible hack due to the fact that the graphviz library API
    # only allows appending to the body. This introduces mutable state, so each
    # time we speculatively compile we accumulate duplicate nodes and edges.
    # The 'correct' solution is to store the nodes and edges that you want to
    # add in a separate data structure, return them from compile_workflow_once,
    # and only once we have reached the fixed point then add them here. Due to
    # labeling and styling and other graphviz metadata that is not trivial, so
    # instead we simply deepcopy and overwite the bodies here.
    # (NOTE: We now have a separate GraphData structure, so the graphviz and
    # networkx representations can probably be removed from the recursion.)
    # (Also note that you have to do this element-wise; you cannot simply write
    # subgraphs_ = subgraphs because that will only overwrite the local binding
    # and thus it will not affect the call site of compile_workflow!)
    # TODO: overwrite the networkx subgraphs. For now this is okay because
    # we are only using the networkx graphs to do an isomorphism check in the
    # regression tests, in which case identical duplication will not matter.
    for i in range(len(subgraphs_)):
        subgraphs_[i].graphviz.body = subgraphs[i].graphviz.body
        subgraphs_[i].graphdata.name = subgraphs[i].graphdata.name
        subgraphs_[i].graphdata.nodes = subgraphs[i].graphdata.nodes
        subgraphs_[i].graphdata.edges = subgraphs[i].graphdata.edges
        subgraphs_[i].graphdata.subgraphs = subgraphs[i].graphdata.subgraphs

    if iter == max_iters:
        import yaml
        print(yaml.dump(node_data.yml))
        raise Exception(f'Error! Maximum number of iterations ({max_iters}) reached in compile_workflow!')
    return compiler_info


def compile_workflow_once(yaml_tree_ast: YamlTree,
                     args: argparse.Namespace,
                     namespaces: Namespaces,
                     subgraphs: List[GraphReps],
                     explicit_edge_defs: ExplicitEdgeDefs,
                     explicit_edge_calls: ExplicitEdgeCalls,
                     tools: Tools,
                     is_root: bool,
                     relative_run_path: bool) -> CompilerInfo:
    """STOP: Have you read the Developer's Guide?? docs/devguide.md\n
    Recursively compiles yml workflow definition ASTs to CWL file contents

    Args:
        yaml_tree_ast (YamlTree): A tuple of name and yml AST
        args (argparse.Namespace): The command line arguments
        namespaces (Namespaces): Specifies the path in the yml AST to the current subworkflow
        subgraphs (List[Graph]): The graphs associated with the parent workflows of the current subworkflow
        explicit_edge_defs (ExplicitEdgeDefs): Stores the (path, value) of the explicit edge definition sites
        explicit_edge_calls (ExplicitEdgeCalls): Stores the (path, value) of the explicit edge call sites
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl(). yml files that have been compiled to CWL SubWorkflows are also added during compilation.
        is_root (bool): True if this is the root workflow
        relative_run_path (bool): Controls whether to use subdirectories or just one directory when writing the compiled CWL files to disk

    Raises:
        Exception: If any errors occur

    Returns:
        CompilerInfo: Contains the data associated with compiled subworkflows (in the Rose Tree) together with mutable cumulative environment information which needs to be passed through the recursion.
    """
    # NOTE: Use deepcopy so that when we delete wic: we don't modify any call sites
    (yaml_path, yaml_tree) = copy.deepcopy(yaml_tree_ast)
    # We also want another copy of the original AST so that if we need to modify it,
    # we can return the modified AST to the call site and re-compile.
    (yaml_path_orig, yaml_tree_orig) = copy.deepcopy(yaml_tree_ast)

    print(' starting', ('  ' * len(namespaces)) + yaml_path)

    # Check for top-level yml dsl args
    wic = {'wic': yaml_tree.get('wic', {})}
    #import yaml; print(yaml.dump(wic))
    if 'wic' in yaml_tree:
        del yaml_tree['wic']
    wic_steps = wic['wic'].get('steps', {})
    
    yaml_stem = Path(yaml_path).stem

    (back_name_, yaml_tree) = utils.extract_backend(yaml_tree, wic['wic'], Path(yaml_path))
    steps: List[Yaml] = yaml_tree['steps']

    steps_keys = utils.get_steps_keys(steps)

    subkeys = [key for key in steps_keys if key not in tools]

    # Add headers
    yaml_tree['cwlVersion'] = 'v1.2' # Use 1.2 to support conditional workflows
    yaml_tree['class'] = 'Workflow'
    yaml_tree['$namespaces'] = {'edam': 'https://edamontology.org/'}
    yaml_tree['$schemas'] = ['https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl']

    # NOTE: currently mutates yaml_tree (maybe)
    maybe_add_subworkflow_requirement(yaml_tree, tools, steps_keys, subkeys)

    # Collect workflow input parameters
    inputs_workflow = {}
    inputs_file_workflow = {}
    
    # Collect the internal workflow output variables
    outputs_workflow = []
    vars_workflow_output_internal = []

    # Copy recursive explicit edge variable definitions and call sites.
    explicit_edge_defs_copy = copy.deepcopy(explicit_edge_defs)
    explicit_edge_calls_copy = copy.deepcopy(explicit_edge_calls)
    # Unlike the first copies which are mutably updated, these are returned
    # unmodified so that we can test that compilation is embedding independent.
    explicit_edge_defs_copy2 = copy.deepcopy(explicit_edge_defs)
    explicit_edge_calls_copy2 = copy.deepcopy(explicit_edge_calls)

    # Collect recursive subworkflow data
    step_1_names = []
    sibling_subgraphs = []
    
    rose_tree_list = []

    graph = subgraphs[-1] # Get the current graph
    graph_gv = graph.graphviz
    graph_nx = graph.networkx
    graphdata = graph.graphdata

    for i, step_key in enumerate(steps_keys):
        step_name_i = utils.step_name_str(yaml_stem, i, step_key)
        stem = Path(step_key).stem
        wic_step_i = wic_steps.get(f'({i+1}, {step_key})', {})

        # Recursively compile subworkflows, adding compiled cwl file contents to tools
        ast_modified = False
        if step_key in subkeys:
            # Extract the sub yaml file that we pre-loaded from disk.
            sub_yaml_tree = YamlTree(step_key, steps[i][step_key])

            # get the label (if any) from the subworkflow
            step_i_wic_graphviz = sub_yaml_tree.yml.get('wic', {}).get('graphviz', {})
            label = step_i_wic_graphviz.get('label', step_key)

            subgraph_gv = graphviz.Digraph(name=f'cluster_{step_key}')
            subgraph_gv.attr(label=label) # str(path)
            subgraph_gv.attr(color='lightblue')  # color of outline
            subgraph_nx = nx.DiGraph()
            graphdata = GraphData(step_key)
            subgraph = GraphReps(subgraph_gv, subgraph_nx, graphdata)

            steps[i].update({step_key: {}}) # delete yml subtree
            sub_compiler_info = compile_workflow(sub_yaml_tree, args, namespaces + [step_name_i], subgraphs + [subgraph], explicit_edge_defs_copy, explicit_edge_calls_copy, tools, False, relative_run_path)

            sub_rose_tree = sub_compiler_info.rose
            rose_tree_list.append(sub_rose_tree)

            sub_node_data: NodeData = sub_rose_tree.data
            sub_env_data = sub_compiler_info.env

            ast_modified = not (sub_yaml_tree.yml == sub_node_data.yml)
            if ast_modified:
                # Propagate the updated yaml_tree (and wic: tags) upwards.
                # Since we already called ast.merge_yml_trees() before initally
                # compiling, the only way the child wic: tags can differ from
                # the parent is if there were modifications during compilation.
                # In other words, it should be safe to simply replace the
                # relevant yaml_tree and wic: tags in the parent with the child
                # values.
                print('AST modified', step_key)
                wic_steps[f'({i+1}, {step_key})'] = {'wic': sub_node_data.yml.get('wic', {})}
                wic_step_i = wic_steps.get(f'({i+1}, {step_key})', {})
                #import yaml
                #print(yaml.dump(wic_steps))

            sibling_subgraphs.append(sub_node_data.graph) # TODO: Just subgraph?
            step_1_names.append(sub_node_data.step_name_1)
            # Add compiled cwl file contents to tools
            # NOTE: We need to consider what to do when stem is not unique.
            # This can happen when two yml files in different subdirectories
            # have the same stem, and it can also happen when we reuse a
            # subworkflow but due to parameter passing the compiled CWL is different.
            if stem in tools:
                # NOTE: Due to speculative compilation, this is probably fine.
                print(f'Overwriting {stem}')
            tools[stem] = Tool(stem + '.cwl', sub_node_data.compiled_cwl)

            # Initialize the above from recursive values.
            # Do not initialize inputs_workflow. See comment below.
            # inputs_workflow.update(sub_node_data.inputs_workflow)
            inputs_namespaced = dict([(f'{step_name_i}___{k}', val) for k, val in sub_env_data.inputs_file_workflow.items()]) # _{step_key}_input___{k}
            inputs_file_workflow.update(inputs_namespaced)
            vars_workflow_output_internal += sub_env_data.vars_workflow_output_internal
            explicit_edge_defs_copy.update(sub_env_data.explicit_edge_defs)
            explicit_edge_calls_copy.update(sub_env_data.explicit_edge_calls)

        tool_i = tools[stem].cwl
        utils.make_tool_DAG(stem, tools[stem])

        # Add run tag, using relative or flat-directory paths
        # NOTE: run: path issues were causing test_cwl_embedding_independence()
        # to fail, so I simply ignore the run tag in that test.
        run_path = tools[stem].run_path
        if relative_run_path:
            if step_key in subkeys:
                run_path = step_name_i + '/' + run_path
            else:
                run_path = ('../' * len(namespaces)) + run_path
                run_path = '../' + run_path # Relative to autogenerated/
        else:
            if step_key in subkeys:
                run_path = '___'.join(namespaces + [step_name_i, run_path])

        if steps[i][step_key]:
            if not 'run' in steps[i][step_key]:
                steps[i][step_key].update({'run': run_path})
        else:
            steps[i] = {step_key: {'run': run_path}}

        # TODO: Add label and doc tags

        # Generate intermediate file names between steps.
        # The inputs for a given step need to come from either
        # 1. the input.yaml file for the overall workflow (extract into separate yml file) or
        # 2. the outputs from the previous steps (most recent first).
        # If there isn't an exact match, remove input_* and output_* from the
        # current step and previous steps, respectively, and then check again.
        # If there still isn't an exact match, explicit renaming may be required.

        # cachedir_path needs to be an absolute path, but for reproducibility
        # we don't want users' home directories in the yml files.
        if 'cwl_watcher' in step_key:
            cachedir_path = Path(args.cachedir).absolute()
            #print('setting cachedir_path to', cachedir_path)
            steps[i][step_key]['in']['cachedir_path'] = str(cachedir_path)

        args_provided = []
        if steps[i][step_key] and 'in' in steps[i][step_key]:
            args_provided = list(steps[i][step_key]['in'])
        #print(args_provided)

        in_tool = tool_i['inputs']
        #print(list(in_tool.keys()))
        if tool_i['class'] == 'CommandLineTool':
            args_required = [arg for arg in in_tool if not (in_tool[arg].get('default') or in_tool[arg]['type'][-1] == '?')]
        elif tool_i['class'] == 'Workflow':
            args_required = [arg for arg in in_tool]
            
            # Add the inputs. For now, assume that all sub-workflows have been
            # auto-generated from a previous application of this script, and
            # thus that all of their transitive inputs have been satisfied.
            # (i.e. simply combine the input yml files using the cat command.)
            steps[i][step_key]['in'] = dict([(key, key) for key in args_required])
        else:
            raise Exception(f'Unknown class', tool_i['class'])
        
        # Note: Some config tags are not required in the cwl files, but are in
        # fact required in the python source code! See check_mandatory_property
        # (Solution: refactor all required arguments out of config and list
        # them as explicit inputs in the cwl files, then modify the python
        # files accordingly.)
        #print(args_required)

        sub_args_provided = [arg for arg in args_required if arg in explicit_edge_calls_copy]
        #print(sub_args_provided)

        label = step_key
        if args.graph_label_stepname:
            label = step_name_i
        step_node_name = '___'.join(namespaces + [step_name_i])

        if not tool_i['class'] == 'Workflow':
            wic_graphviz_step_i = wic_step_i.get('wic', {}).get('graphviz', {})
            label = wic_graphviz_step_i.get('label', label)
            default_style = 'rounded, filled'
            style = wic_graphviz_step_i.get('style', '')
            style = default_style if style == '' else default_style + ', ' + style
            attrs = {'label': label, 'shape': 'box', 'style': style, 'fillcolor': 'lightblue'}
            graph_gv.node(step_node_name, **attrs)
            graph_nx.add_node(step_node_name)
            graphdata.nodes.append((step_node_name, attrs))
        elif not (step_key in subkeys and len(namespaces) < args.graph_inline_depth):
            nssnode = namespaces + [step_name_i]
            # Just like in add_graph_edge(), here we can hide all of the details
            # below a given depth by simply truncating the node's namespaces.
            nssnode = nssnode[:(1 + args.graph_inline_depth)]
            step_node_name = '___'.join(nssnode)
            # NOTE: NOT wic_graphviz_step_i
            # get the label (if any) from the subworkflow
            # TODO: This causes test_cwl_embedding_independence to fail.
            #yml = sub_node_data.yml if ast_modified else sub_yaml_tree.yml
            #step_i_wic_graphviz = yml.get('wic', {}).get('graphviz', {})
            # TODO: For file format conversions, figure out why this is using
            # the label from the parent workflow.
            #label = step_i_wic_graphviz.get('label', label)
            default_style = 'rounded, filled'
            style = '' #step_i_wic_graphviz.get('style', '')
            style = default_style if style == '' else default_style + ', ' + style
            attrs = {'label': label, 'shape': 'box', 'style': style, 'fillcolor': 'lightblue'}
            graph_gv.node(step_node_name, **attrs)
            graph_nx.add_node(step_node_name)
            graphdata.nodes.append((step_node_name, attrs))

        # NOTE: sub_args_provided are handled within the args_required loop below
        for arg_key in args_provided:
            # Extract input value into separate yml file
            # Replace it here with a new variable name
            arg_val = steps[i][step_key]['in'][arg_key]
            # Convert native YAML to a JSON-encoded string for specific tags.
            tags = ['config']
            if arg_key in tags and isinstance(arg_val, Dict):
                arg_val = json.dumps(arg_val)
            in_name = f'{step_name_i}___{arg_key}'  # Use triple underscore for namespacing so we can split later # {step_name_i}_input___{arg_key}
            in_type = in_tool[arg_key]['type'].replace('?', '')  # Providing optional arguments makes them required
            in_dict = {'type': in_type}
            if 'format' in in_tool[arg_key]:
                in_format = in_tool[arg_key]['format']
                in_dict['format'] = in_format
            if arg_val[0] == '&':
                arg_val = arg_val[1:]  # Remove &
                #print('arg_key, arg_val', arg_key, arg_val)
                # NOTE: There can only be one definition, but multiple call sites.
                if not explicit_edge_defs_copy.get(arg_val):
                    # if first time encountering arg_val, i.e. if defining
                    inputs_workflow.update({in_name: in_dict})
                    in_dict = {**in_dict, 'value': arg_val}
                    inputs_file_workflow.update({in_name: in_dict})
                    steps[i][step_key]['in'][arg_key] = in_name
                    explicit_edge_defs_copy.update({arg_val: (namespaces + [step_name_i], arg_key)})
                    # Add a 'dummy' value to explicit_edge_calls, because
                    # that determines sub_args_provided when the recursion returns.
                    explicit_edge_calls_copy.update({in_name: (namespaces + [step_name_i], arg_key)})
                    # TODO: Show input node?
                else:
                    raise Exception(f"Error! Multiple definitions of &{arg_val}!")
            elif arg_val[0] == '*' and 'cwl_watcher' not in step_key:
                # NOTE: Exclude cwl_watcher from explicit edge dereferences.
                # Since cwl_watcher requires explicit filenames for globbing,
                # we do not want to replace them with internal CWL dependencies!
                arg_val = arg_val[1:]  # Remove *
                if not explicit_edge_defs_copy.get(arg_val):
                    if is_root:
                        # TODO: Check this comment.
                        # Even if is_root, we don't want to raise an Exception
                        # here because in test_cwl_embedding_independence, we
                        # recompile all subworkflows as if they were root. That
                        # will cause this code path to be taken but it is not
                        # actually an error. Add a CWL input for testing only.
                        print(f"Error! No definition found for &{arg_val}!")
                        print(f"Creating the CWL input {in_name} anyway, but")
                        print("without any corresponding input value this will fail validation!")
                    inputs_workflow.update({in_name: in_dict})
                    steps[i][step_key]['in'][arg_key] = in_name
                    # Add a 'dummy' value to explicit_edge_calls anyway, because
                    # that determines sub_args_provided when the recursion returns.
                    explicit_edge_calls_copy.update({in_name: (namespaces + [step_name_i], arg_key)})
                else:
                    (nss_def_init, var) =  explicit_edge_defs_copy[arg_val]

                    nss_def_embedded = var.split('___')[:-1]
                    nss_call_embedded = arg_key.split('___')[:-1]
                    nss_def = nss_def_init + nss_def_embedded
                    # [step_name_i] is correct; nss_def_init already contains step_name_j from the recursive call
                    nss_call = namespaces + [step_name_i] + nss_call_embedded

                    nss_def_inits, nss_def_tails = utils.partition_by_lowest_common_ancestor(nss_def, nss_call)
                    nss_call_inits, nss_call_tails = utils.partition_by_lowest_common_ancestor(nss_call, nss_def)
                    # nss_def and nss_call are paths into the abstract syntax tree 'call stack'.
                    # This defines the 'common namespace' in the call stack w.r.t. the inits.
                    assert nss_def_inits == nss_call_inits
                    
                    # TODO: Check this comment.
                    # Relative to the common namespace, if the call site of an explicit
                    # edge is at a depth > 1, (i.e. if it is NOT simply of the form
                    # last_namespace/input_variable) then we
                    # need to create inputs in all of the intervening CWL files
                    # so we can pass in the values from the outer scope(s). Here,
                    # we simply need to use in_name and add to inputs_workflow
                    # and explicit_edge_calls. The outer scope(s) are handled by
                    # the sub_args_provided clause below.
                    # Note that the outputs from the definition site are bubbled
                    # up the call stack until they reach the common namespace.
                    if len(nss_call_tails) > 1:
                        inputs_workflow.update({in_name: in_dict})
                        steps[i][step_key]['in'][arg_key] = in_name
                        # Store explicit edge call site info up through the recursion.
                        explicit_edge_calls_copy.update({in_name: explicit_edge_defs_copy[arg_val]}) # {in_name, (namespaces + [step_name_i], var)} ?
                    elif len(nss_call_tails) == 1:
                        # TODO: Check this comment.
                        # The definition site recursion (only, if any) has completed
                        # and we are already in the common namespace, thus
                        # we need to pass in the value from the definition site.
                        # Note that since len(nss_call_tails) == 1,
                        # there will not be any call site recursion in this case.
                        var_slash = nss_def_tails[0] + '/' + '___'.join(nss_def_tails[1:] + [var])
                        steps[i][step_key]['in'][arg_key] = var_slash
                    else:
                        # Since nss_call includes step_name_i, this should never happen...
                        raise Exception("Error! len(nss_call_tails) == 0! Please file a bug report!\n" +
                                        f'nss_def {nss_def}\n nss_call {nss_call}')

                    # Add an edge, but in a carefully chosen subgraph.
                    # If you add an edge whose head/tail is outside of the subgraph,
                    # graphviz may segfault! Moreover, even if graphviz doesn't
                    # segfault, adding an edge in a given subgraph can cause the
                    # nodes themselves to be rendered in that subgraph, even
                    # though the nodes are defined in a different subgraph!
                    # The correct thing to do is to use the graph associated with
                    # the lowest_common_ancestor of the definition and call site.
                    # (This is the only reason we need to pass in all subgraphs.)
                    label = var.split('___')[-1]
                    graph_init = subgraphs[len(nss_def_inits)]
                    # Let's use regular blue for explicit edges.
                    # Use constraint=false ?
                    utils.add_graph_edge(args, graph_init, nss_def, nss_call, label, color='blue')
            else:
                # NOTE: See comment above about excluding cwl_watcher from
                # explicit edge dereferences.
                if arg_val[0] == '*' and 'cwl_watcher' in step_key and 'file_pattern' not in arg_key:
                    arg_val = arg_val[1:] # Remove *, except for file_pattern

                inputs_workflow.update({in_name: in_dict})
                in_dict = {**in_dict, 'value': arg_val}
                inputs_file_workflow.update({in_name: in_dict})
                steps[i][step_key]['in'][arg_key] = in_name
                if args.graph_show_inputs:
                    input_node_name = '___'.join(namespaces + [step_name_i, arg_key])
                    attrs = {'label': arg_key, 'shape': 'box', 'style': 'rounded, filled', 'fillcolor': 'lightgreen'}
                    graph_gv.node(input_node_name, **attrs)
                    graph_gv.edge(input_node_name, step_node_name, color=font_edge_color)
                    graph_nx.add_node(input_node_name)
                    graph_nx.add_edge(input_node_name, step_node_name)
                    graphdata.nodes.append((input_node_name, attrs))
                    graphdata.edges.append((input_node_name, step_node_name, {}))

        for arg_key in args_required:
            #print('arg_key', arg_key)
            in_name = f'{step_name_i}___{arg_key}'
            if arg_key in args_provided:
                continue  # We already covered this case above.
            if in_name in inputs_file_workflow:
                # We provided an explicit argument (but not an edge) in a subworkflow,
                # and now we just need to pass it up to the root workflow.
                #print('passing', in_name)
                in_type = in_tool[arg_key]['type']
                in_dict = {'type': in_type}
                if 'format' in in_tool[arg_key]:
                    in_format = in_tool[arg_key]['format']
                    in_dict['format'] = in_format
                inputs_workflow.update({in_name: in_dict})
                arg_keyval = {arg_key: in_name}
                steps[i] = utils.add_yamldict_keyval_in(steps[i], step_key, arg_keyval)

                # Obviously since we supplied a value, we do NOT want to perform edge inference.
                continue
            if arg_key in sub_args_provided: # Edges have been explicitly provided
                # The definition site recursion (if any) and the call site
                # recursion (yes, see above), have both completed and we are
                # now in the common namespace, thus 
                # we need to pass in the value from the definition site.
                # Extract the stored defs namespaces from explicit_edge_calls.
                # (See massive comment above.)
                (nss_def_init, var) = explicit_edge_calls_copy[arg_key]

                nss_def_embedded = var.split('___')[:-1]
                nss_call_embedded = arg_key.split('___')[:-1]
                nss_def = nss_def_init + nss_def_embedded
                # [step_name_i] is correct; nss_def_init already contains step_name_j from the recursive call
                nss_call = namespaces + [step_name_i] + nss_call_embedded

                nss_def_inits, nss_def_tails = utils.partition_by_lowest_common_ancestor(nss_def, nss_call)
                nss_call_inits, nss_call_tails = utils.partition_by_lowest_common_ancestor(nss_call, nss_def)
                assert nss_def_inits == nss_call_inits

                nss_call_tails_stems = [utils.parse_step_name_str(x)[0] for x in nss_call_tails]
                arg_val = steps[i][step_key]['in'][arg_key]
                if (nss_call_tails_stems == []) or (yaml_stem in nss_call_tails_stems and nss_call_tails_stems.index(yaml_stem) > 0):
                    # i.e. (if 'dummy' value) or (if it is possible to do more recursion)
                    in_type = in_tool[arg_key]['type']
                    in_dict = {'type': in_type}
                    if 'format' in in_tool[arg_key]:
                        in_format = in_tool[arg_key]['format']
                        in_dict['format'] = in_format
                    inputs_workflow.update({in_name: in_dict})
                    steps[i][step_key]['in'][arg_key] = in_name
                    # Store explicit edge call site info up through the recursion.
                    explicit_edge_calls_copy.update({in_name: explicit_edge_calls_copy[arg_key]})
                else:
                    # TODO: Check this comment.
                    # The definition site recursion (only, if any) has completed
                    # and we are already in the common namespace, thus
                    # we need to pass in the value from the definition site.
                    # Note that since len(nss_call_tails) == 1,
                    # there will not be any call site recursion in this case.
                    var_slash = nss_def_tails[0] + '/' + '___'.join(nss_def_tails[1:] + [var])
                    steps[i][step_key]['in'][arg_key] = var_slash

                # NOTE: We already added an edge to the appropriate subgraph above.
                # TODO: vars_workflow_output_internal?
            else:
                conversions: List[str] = []
                in_name_in_inputs_file_workflow: bool = (in_name in inputs_file_workflow)
                steps[i] = inference.perform_edge_inference(args, tools, steps_keys,
                    yaml_stem, i, steps[i], arg_key, graph, is_root, namespaces,
                    vars_workflow_output_internal, inputs_workflow,
                    in_name_in_inputs_file_workflow, conversions, wic_steps)
                # NOTE: For now, perform_edge_inference mutably appends to
                # inputs_workflow and vars_workflow_output_internal.

                # Automatically insert file format conversion
                conversions = list(set(conversions)) # Remove duplicates
                if not len(conversions) == 0:
                    conversion = conversions[0]
                    print('Automaticaly inserting file format conversion', conversion, i)
                    if not len(conversions) == 1:
                        print('Warning! More than one file format conversion! Choosing', conversion)

                    yaml_tree_mod = yaml_tree_orig
                    steps_mod: List[Yaml] = yaml_tree_mod['steps']
                    steps_mod.insert(i, {conversion: None})

                    # Add inference rules annotations (i.e. for file format conversion)
                    conv_tool = tools[conversion]
                    conv_out_tool = conv_tool.cwl['outputs']
                    inference_rules_dict = dict([(out_key, inference_rules.get(out_val['format'], 'default')) for out_key, out_val in conv_out_tool.items() if 'format' in out_val])
                    keystr = f'({i+1}, {conversion})' # The yml file uses 1-based indexing
                    inf_dict = {'wic': {'inference': inference_rules_dict}}

                    if 'wic' in yaml_tree_mod:
                        if 'steps' in yaml_tree_mod['wic']:
                            yaml_tree_mod['wic']['steps'] = reindex_wic_steps(yaml_tree_mod['wic']['steps'], i)
                            yaml_tree_mod['wic']['steps'][keystr] = inf_dict
                        else:
                            yaml_tree_mod['wic'].update({'steps': {keystr: inf_dict}})
                    else:
                        yaml_tree_mod.update({'wic': {'steps': {keystr: inf_dict}}})

                    node_data = NodeData(namespaces, yaml_stem, yaml_tree_mod, yaml_tree, {}, explicit_edge_defs_copy2, explicit_edge_calls_copy2, graph, inputs_workflow, '')
                    rose_tree = RoseTree(node_data, rose_tree_list)
                    env_data = EnvData(inputs_file_workflow, vars_workflow_output_internal, explicit_edge_defs_copy, explicit_edge_calls_copy)
                    compiler_info = CompilerInfo(rose_tree, env_data)
                    #node_data_dummy = NodeData(None, None, yaml_tree_mod, None, None, None, None, None, None, None)
                    #compiler_info_dummy = CompilerInfo(RoseTree(node_data_dummy, None), None)
                    return compiler_info
        
        # Add CommandLineTool outputs tags to workflow out tags.
        # Note: Add all output tags for now, but depending on config options,
        # not all output files will be generated. This may cause an error.
        out_keyvals = {}
        for out_key, out_dict in tool_i['outputs'].items():
            out_keyvals[out_key] = {'type': out_dict['type'], 'format': out_dict['format']}
            #print(out_key, out_keyvals[out_key])
        if out_keyvals == {}: # FYI this should never happen
            print(f'Error! no outputs for step {step_key}')
        outputs_workflow.append(out_keyvals)

        steps[i] = utils.add_yamldict_keyval_out(steps[i], step_key, list(tool_i['outputs']))

        #print()

    # NOTE: add_subgraphs currently mutates graph
    wic_graphviz = wic['wic'].get('graphviz', {})
    ranksame_strs = wic_graphviz.get('ranksame', [])
    ranksame_pairs = [utils.parse_int_string_tuple(x) for x in ranksame_strs]
    steps_ranksame = ['___'.join(namespaces + [utils.step_name_str(yaml_stem, num-1, name)]) for num, name in ranksame_pairs]
    steps_ranksame = [f'"{x}"' for x in steps_ranksame]  # Escape with double quotes.
    add_subgraphs(args, graph, sibling_subgraphs, namespaces, step_1_names, steps_ranksame)
    step_name_1 = get_step_name_1(step_1_names, yaml_stem, namespaces, steps_keys, subkeys)

    # Add the provided inputs of each step to the workflow inputs
    yaml_tree.update({'inputs': inputs_workflow})

    vars_workflow_output_internal = list(set(vars_workflow_output_internal))  # Get uniques
    # (Why are we getting uniques?)
    workflow_outputs = get_workflow_outputs(args, namespaces, is_root, yaml_stem, steps, outputs_workflow, vars_workflow_output_internal, graph, tools, step_node_name)
    yaml_tree.update({'outputs': workflow_outputs})

    # Finally, rename the steps to be unique
    # and convert the list of steps into a dict
    steps_dict = {}
    for i, step_key in enumerate(steps_keys):
        step_name_i = utils.step_name_str(yaml_stem, i, step_key)
        #steps[i] = {step_name_i: steps[i][step_key]}
        steps_dict.update({step_name_i: steps[i][step_key]})
    yaml_tree.update({'steps': steps_dict})

    # Dump the workflow inputs to a separate yml file.
    yaml_inputs: WorkflowInputsFile = {}
    for key, in_dict in inputs_file_workflow.items():
        new_keyval: WorkflowInputsFile = {}
        if 'File' == in_dict['type']:
            in_format = in_dict['format']
            if isinstance(in_format, List):
                in_format = list(set(in_format)) # get uniques
                if len(in_format) > 1:
                    print(f'NOTE: More than one input file format for {key}')
                    print(f'formats: {in_format}')
                    print(f'Choosing {in_format[0]}')
                in_format = in_format[0]
            new_keyval = {key: {'class': 'File', 'path': in_dict['value'], 'format': in_format}}
        elif 'string' == in_dict['type']:
            # We cannot store string values as a dict, so use type: ignore
            new_keyval = {key: in_dict['value']} # type: ignore 
        else:
            raise Exception(f"Error! Unknown type: {in_dict['type']}")
        yaml_inputs.update(new_keyval)

    if args.cwl_validate:
        print(f'Validating {yaml_stem}.cwl ...')
        cmd = ['cwltool', '--validate', f'{yaml_stem}.cwl']
        sub.run(cmd)

    print('finishing', ('  ' * len(namespaces)) + yaml_path)
    # Note: We do not necessarily need to return inputs_workflow.
    # 'Internal' inputs are encoded in yaml_tree. See Comment above.
    node_data = NodeData(namespaces, yaml_stem, yaml_tree_orig, yaml_tree, yaml_inputs, explicit_edge_defs_copy2, explicit_edge_calls_copy2, graph, inputs_workflow, step_name_1)
    rose_tree = RoseTree(node_data, rose_tree_list)
    env_data = EnvData(inputs_file_workflow, vars_workflow_output_internal, explicit_edge_defs_copy, explicit_edge_calls_copy)
    compiler_info = CompilerInfo(rose_tree, env_data)
    return compiler_info


def reindex_wic_steps(wic_steps: Yaml, index: int) -> Yaml:
    """After inserting a step into a workflow, we need to increment the steps in the wic: metadata annotations tag whose original index is >= the given index.

    Args:
        wic_steps (Yaml): The steps: subtag of the wic: metadata annotations tag.
        index (int): The (zero-based) index of the inserted workflow step.

    Returns:
        Yaml: The updated wic: steps: tag, with the appropriate indices incremented.
    """
    wic_steps_reindexed = {}
    for keystr, val in wic_steps.items():
        (i, s) = utils.parse_int_string_tuple(keystr)
        newstr = f'({i+1}, {s})' if i >= index else keystr
        wic_steps_reindexed[newstr] = val
    return wic_steps_reindexed


def maybe_add_subworkflow_requirement(yaml_tree: Yaml, tools: Tools, steps_keys: List[str], subkeys: List[str]) -> None:
    """Adds a SubworkflowFeatureRequirement if there are any subworkflows

    Args:
        yaml_tree (Yaml): A tuple of name and yml AST
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()
        steps_keys (List[str]): The name of each step in the current CWL workflow
        subkeys (List[str]): The keys associated with subworkflows
    """
    # If there is at least one subworkflow, add a SubworkflowFeatureRequirement
    if (not subkeys == []) or any([tools[Path(key).stem].cwl['class'] == 'Workflow' for key in steps_keys if key not in subkeys]):
        subworkreq = 'SubworkflowFeatureRequirement'
        subworkreqdict = {subworkreq: {'class': subworkreq}}
        if 'requirements' in yaml_tree:
            if not subworkreq in yaml_tree['requirements']:
                new_reqs = dict(list(yaml_tree['requirements'].items()) + list(subworkreqdict))
                yaml_tree['requirements'].update(new_reqs)
        else:
            yaml_tree['requirements'] = subworkreqdict


def add_subgraphs(args: argparse.Namespace,
                  graph: GraphReps,
                  sibling_subgraphs: List[GraphReps],
                  namespaces: Namespaces,
                  step_1_names: List[str],
                  steps_ranksame: List[str]) -> None:
    """Add all subgraphs to the current graph, except for GraphViz subgraphs
    below a given depth, which allows us to hide irrelevant details.

    Args:
        args (argparse.Namespace): The command line arguments
        graph (GraphReps): A tuple of a GraphViz DiGraph and a networkx DiGraph
        sibling_subgraphs (List[Graph]): The subgraphs of the immediate children of the current workflow
        namespaces (Namespaces): Specifies the path in the AST of the current subworkflow
        step_1_names (List[str]): The names of the first step
        steps_ranksame (List[str]): Additional node names to be aligned using ranksame
    """
    graph_gv = graph.graphviz
    graph_nx = graph.networkx
    # Add the cluster subgraphs to the main graph, but we need to add them in
    # reverse order to trick the graphviz layout algorithm.
    for sibling in sibling_subgraphs[::-1]: # Reverse!
        (sib_graph_gv, sib_graph_nx, sib_graphdata) = sibling
        if len(namespaces) < args.graph_inline_depth:
            graph_gv.subgraph(sib_graph_gv)
        graph_nx.add_nodes_from(sib_graph_nx.nodes)
        graph_nx.add_edges_from(sib_graph_nx.edges)
    for sibling in sibling_subgraphs:
        graph.graphdata.subgraphs.append(sibling.graphdata)
    # Align the cluster subgraphs using the same rank as the first node of each subgraph.
    # See https://stackoverflow.com/questions/6824431/placing-clusters-on-the-same-rank-in-graphviz
    if len(namespaces) < args.graph_inline_depth:
        step_1_names_display = [name for name in step_1_names if len(name.split('___')) < 2 + args.graph_inline_depth]
        if len(step_1_names_display) > 1:
            nodes_same_rank = '\t{rank=same; ' + '; '.join(step_1_names_display) + '}\n'
            graph_gv.body.append(nodes_same_rank)
            graph.graphdata.ranksame = step_1_names_display
        if len(steps_ranksame) > 1:
            nodes_same_rank = '\t{rank=same; ' + '; '.join(steps_ranksame) + '}\n'
            graph_gv.body.append(nodes_same_rank)
            graph.graphdata.ranksame = steps_ranksame
        


def get_step_name_1(step_1_names: List[str],
                    yaml_stem: str,
                    namespaces: Namespaces,
                    steps_keys: List[str],
                    subkeys: List[str]) -> str:
    """Finds the name of the first step in the current subworkflow. If the first
    step is itself subworkflow, the call site recurses until it finds a node.
    This is necessary because ranksame in GraphViz can only be applied to
    individual nodes, not cluster_subgraphs.

    Args:
        step_1_names (List[str]): The list of potential first node names
        yaml_stem (str): The name of the current subworkflow (stem of the yaml filepath)
        namespaces (Namespaces): Specifies the path in the AST of the current subworkflow
        steps_keys (List[str]): The name of each step in the current CWL workflow
        subkeys (List[str]): The keys associated with subworkflows

    Returns:
        str: _description_
    """
    if steps_keys[0] in subkeys:
        step_name_1 = step_1_names[0]
    else:
        step_name_1 = utils.step_name_str(yaml_stem, 0, steps_keys[0])
        step_name_1 = '___'.join(namespaces + [step_name_1])
    # NOTE: Since the names of subgraphs '*.yml' contain a period, we need to
    # escape them by enclosing the whole name in double quotes. Otherwise:
    # "Error: *.yml.gv: syntax error in line n near '.'"
        step_name_1 = f'"{step_name_1}"'

    return step_name_1


def get_workflow_outputs(args: argparse.Namespace,
                         namespaces: Namespaces,
                         is_root: bool,
                         yaml_stem: str,
                         steps: List[Yaml],
                         outputs_workflow: WorkflowOutputs,
                         vars_workflow_output_internal: InternalOutputs,
                         graph: GraphReps,
                         tools: Tools,
                         step_node_name: str) -> Dict[str, Dict[str, str]]:
    """Chooses a subset of the CWL outputs: to actually output

    Args:
        args (argparse.Namespace): The command line arguments
        namespaces (Namespaces): Specifies the path in the AST of the current subworkflow
        is_root (bool): True if this is the root workflow
        yaml_stem (str): The name of the current subworkflow (stem of the yaml filepath)
        steps (List[Yaml]): The steps: tag of a CWL workflow
        outputs_workflow (WorkflowOutputs): Contains the contents of the out: tags for each step.
        vars_workflow_output_internal (InternalOutputs): Keeps track of output variables which are internal to the root workflow, but not necessarily to subworkflows.
        graph (GraphReps): A tuple of a GraphViz DiGraph and a networkx DiGraph
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()
        step_node_name (str): The namespaced name of the current step

    Returns:
        Dict[str, Dict[str, str]]: The actual outputs to be specified in the generated CWL file
    """
    # Add the outputs of each step to the workflow outputs
    workflow_outputs = {}
    steps_keys = utils.get_steps_keys(steps)
    for i, step_key in enumerate(steps_keys):
        stem = Path(step_key).stem
        tool_i = tools[stem].cwl
        step_name_i = utils.step_name_str(yaml_stem, i, step_key)
        out_keys = steps[i][step_key]['out']
        for out_key in out_keys:
            out_var = f'{step_name_i}/{out_key}'
            # Avoid duplicating intermediate outputs in GraphViz
            out_key_no_namespace = out_key.split('___')[-1]
            if args.graph_show_outputs:
                case1 = (tool_i['class'] == 'Workflow') and (not out_key in [var.replace('/', '___') for var in vars_workflow_output_internal])
                # Avoid duplicating outputs from subgraphs in parent graphs.
                output_node_name = '___'.join(namespaces + [step_name_i, out_key])
                # TODO: check is_root here
                case1 = case1 and not is_root and (len(step_node_name.split('___')) + 1 == len(output_node_name.split('___')))
                case2 = (tool_i['class'] == 'CommandLineTool') and (not out_var in vars_workflow_output_internal)
                if case1 or case2:
                    graph_gv = graph.graphviz
                    graph_nx = graph.networkx
                    graphdata = graph.graphdata
                    attrs = {'label': out_key_no_namespace, 'shape': 'box', 'style': 'rounded, filled', 'fillcolor': 'lightyellow'}
                    graph_gv.node(output_node_name, **attrs)
                    if args.graph_label_edges:
                        graph_gv.edge(step_node_name, output_node_name, color=font_edge_color, label=out_key_no_namespace)  # Is labeling necessary?
                    else:
                        graph_gv.edge(step_node_name, output_node_name, color=font_edge_color)
                    graph_nx.add_node(output_node_name)
                    graph_nx.add_edge(step_node_name, output_node_name)
                    graphdata.nodes.append((output_node_name, attrs))
                    graphdata.edges.append((step_node_name, output_node_name, {}))
            # NOTE: Unless we are in the root workflow, we always need to
            # output everything. This is because while we are within a
            # subworkflow, we do not yet know if a subworkflow output will be used as
            # an input in a parent workflow (either explicitly, or using inference).
            # Thus, we have to output everything.
            # However, once we reach the root workflow, we can choose whether
            # we want to pollute our output directory with intermediate files.
            # (You can always enable --provenance to get intermediate files.)
            # NOTE: Remove is_root for now because in test_cwl_embedding_independence,
            # we recompile all subworkflows as if they were root. Thus, for now
            # we need to enable args.cwl_output_intermediate_files
            # Exclude intermediate 'output' files.
            if out_var in vars_workflow_output_internal and not args.cwl_output_intermediate_files: # and is_root
                continue
            out_name = f'{step_name_i}___{out_key}'  # Use triple underscore for namespacing so we can split later
            #print('out_name', out_name)

        for out_key, out_dict in outputs_workflow[i].items():
            out_name = f'{step_name_i}___{out_key}'  # Use triple underscore for namespacing so we can split later
            out_var = f'{step_name_i}/{out_key}'
            workflow_outputs.update({out_name: {**out_dict, 'outputSource': out_var}})
        #print('workflow_outputs', workflow_outputs)
    return workflow_outputs