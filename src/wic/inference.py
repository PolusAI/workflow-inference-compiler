import argparse
from pathlib import Path
from typing import List, Tuple

from . import utils
from .wic_types import Namespaces, KV, Tools, WorkflowInputs, InternalOutputs, Graph


def perform_edge_inference(args: argparse.Namespace,
                           tools: Tools, steps_keys: List[str],
                           yaml_stem: str,
                           i: int,
                           steps_i: KV,
                           arg_key: str,
                           graph: Graph,
                           is_root: bool,
                           namespaces: Namespaces,
                           vars_workflow_output_internal: InternalOutputs,
                           inputs_workflow: WorkflowInputs,
                           in_name_in_inputs_file_workflow: bool) -> KV:
    """This function implements the core edge inference feature.
    NOTE: steps_i, vars_workflow_output_internal, inputs_workflow are mutably updated.

    Args:
        args (argparse.Namespace): The command line arguments
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()
        steps_keys (List[str]): The name of each step in the current CWL workflow
        yaml_stem (str): The name (filename without extension) of the current CWL workflow
        i (int): The (zero-based) step number w.r.t. the current subworkflow. Since we are trying to infer inputs from previous outputs, this will not perform any inference (again, w.r.t. the current subworkflow) if i == 0.
        steps_i (KV): The i-th component of the steps: tag of the current CWL workflow
        arg_key (str): The name of the CWL input tag that needs a concrete input value inferred
        graph (Graph): A tuple of a GraphViz DiGraph and a networkx DiGraph
        is_root (bool): True if this is the root workflow (for debugging only)
        namespaces (Namespaces): Specifies the path in the AST of the current subworkflow
        vars_workflow_output_internal (InternalOutputs): Keeps track of output variables which are internal to the root workflow, but not necessarily to subworkflows.
        inputs_workflow (WorkflowInputs): Keeps track of CWL inputs: variables for the current workflow.
        in_name_in_inputs_file_workflow (bool): Used to determine whether failure to find a match should be considered an error.

    Returns:
        KV: steps_i with the input tag arg_key updated with an inferred input value.
    """
    # Use in_name_in_inputs_file_workflow: bool so that at the call site, we don't
    # have to question whether or not this function modifies inputs_file_workflow
    # gmx_trjconv_str assumes gro, but other tools assume tpr; For now,
    # instead of figuring out how to replace('structure_', 'gro_') for
    # gmx_trjconv_str only, just require the user to specify filename.
    step_key = steps_keys[i]
    tool_i = tools[Path(step_key).stem][1]
    step_name_i = utils.step_name_str(yaml_stem, i, step_key)

    in_tool = tool_i['inputs']
    in_type = in_tool[arg_key]['type']
    in_type = in_type.replace('?', '')  # Providing optional arguments makes them required
    in_dict = {'type': in_type}
    in_formats = []
    if 'format' in in_tool[arg_key]:
        in_formats = in_tool[arg_key]['format']
        in_dict['format'] = in_formats
    format_matches_all = []
    attempted_matches_all = []
    for j in range(0, i)[::-1]:  # Reverse order!
        tool_j = tools[Path(steps_keys[j]).stem][1]
        out_tool = tool_j['outputs']
        out_keys = list(tool_j['outputs'])[::-1] # Reverse order!
        format_matches = []
        attempted_matches = []
        for out_key in out_keys:
            out_type = out_tool[out_key]['type']
            out_dict = {'type': out_type}
            out_format = ''
            if 'format' in out_tool[out_key]:
                out_format = out_tool[out_key]['format']
                out_dict['format'] = out_format
            if out_format == '':
                print(f'Warning! No output format! Cannot possibly match!')
                print(f'out_key {out_key}')
            attempted_matches.append((out_key, out_format))
            # Great! We found an 'exact' type and format match.
            if (out_type == in_type) and out_format in in_formats:
                format_matches.append((out_key, out_format))
            # TODO: elif there exists a file format conversion...
        format_matches_all.append(format_matches)
        attempted_matches_all.append(attempted_matches)
        # Most log files just have format_2330 "Textual format", but this can
        # conflict with other structured text files that also have format_2330.
        # Unfortunately, the edam formats are not very finely curated, so they
        # are not very 'exact'. For now, we can perform additional matching
        # based on naming conventions and/or we can simply exclude log files
        # (which are not usually parsed or otherwise used as inputs).
        # Eventually, we will want to improve the format curation.
        # NOTE: Use underscores to prevent excluding e.g. 'topology'
        # This isn't great, but works for now (until someone uses '_log_' ...)
        # format_matches = [x for x in format_matches if not '_log_' in x[0]]
        if not len(format_matches) == 0:
            if len(format_matches) == 1:
                # Great! We found a unique format match.
                out_key = format_matches[0][0]
            else:
                name_matches = []
                # NOTE: The biobb CWL files do not use consistent naming
                # conventions, so we need to perform some renamings here.
                # Eventually, the CWL files themselves should be fixed.
                arg_key_no_namespace = arg_key.split('___')[-1]
                arg_key_renamed = arg_key_no_namespace.replace('input_', '')
                renamings = [('energy_', 'edr_'), ('structure_', 'tpr_'), ('traj_', 'trr_')]
                for name1, name2 in renamings:
                    arg_key_renamed = arg_key_renamed.replace(name1, name2)

                for out_key, out_format in format_matches:
                    out_key_no_namespace = out_key.split('___')[-1]
                    out_key_renamed = out_key_no_namespace.replace('output_', '')
                    if arg_key_renamed == out_key_renamed:
                        name_matches.append((out_key, out_format))

                if len(name_matches) == 0:
                    #print(f'Found multiple outputs with compatible types and formats (but no matching names) for input {arg_key}')
                    #for x in format_matches:
                    #    print(x)
                    #print(f'Arbitrarily choosing the first match {format_matches[0][0]}')
                    out_key = format_matches[0][0]
                elif len(name_matches) == 1:
                    # Great! We found a unique match.
                    out_key = name_matches[0][0]
                else:
                    #print(f'Found multiple outputs with compatible types and formats (and multiple matching names) for input {arg_key}')
                    #for x in name_matches:
                    #    print(x)
                    #print(f'Arbitrarily choosing the first match {name_matches[0][0]}')
                    out_key = name_matches[0][0]
            #print('match!', j)  # We found a match!
            # Generate a new namespace for out_key using the step number and add to inputs
            step_name_j = utils.step_name_str(yaml_stem, j, steps_keys[j])
            arg_val = f'{step_name_j}/{out_key}'

            # We also need to keep track of the 'internal' output variables
            if tool_j['class'] == 'Workflow':
                vars_workflow_output_internal.append(out_key)
            else:
                vars_workflow_output_internal.append(arg_val)

            # Determine which head and tail node to use for the new edge
            # First we need to extract the embedded namespaces
            nss_embedded1 = out_key.split('___')[:-1]
            nss_embedded2 = arg_key.split('___')[:-1]
            nss1 = namespaces + [step_name_j] + nss_embedded1
            nss2 = namespaces + [step_name_i] + nss_embedded2
            # TODO: check this
            out_key_no_namespace = out_key.split('___')[-1]
            label = out_key_no_namespace if tool_j['class'] == 'Workflow' else out_key
            utils.add_graph_edge(args, graph, nss1, nss2, label)

            arg_keyval = {arg_key: arg_val}
            steps_i = utils.add_yamldict_keyval_in(steps_i, step_key, arg_keyval)
            #print(f'inference i {i} y arg_key {arg_key}')
            return steps_i  # Short circuit

    match = False
    if not match:
        #print(f'inference i {i} n arg_key {arg_key}')
        in_name = f'{step_name_i}___{arg_key}'  # Use triple underscore for namespacing so we can split later # {step_name_i}_input___{arg_key}

        # This just means we need to defer to the parent workflow.
        # There will actually be an error only if no parent supplies an
        # input value and thus there is no entry in inputs_file_workflow.
        if is_root and not in_name_in_inputs_file_workflow:
            print('Error! No match found for input', i + 1, step_key, arg_key)
            print('number of attempted matches ', len(utils.flatten(attempted_matches_all)))
            for y in attempted_matches_all:
                for x in y:
                    print(x)

        # Add an input name to this subworkflow (only). Do not add to
        # inputs_file_workflow because this may be an internal input,
        # i.e. it may simply be split across subworkflow boundaries and
        # the input may be supplied in one of the parent workflows.
        # We also do not (necessarily) need to explicitly pass this list
        # to the parent workflow since we are compiling the 'in' tag here,
        # which should match in the parent workflow.
        inputs_workflow.update({in_name: in_dict})

        arg_keyval = {arg_key: in_name}
        steps_i = utils.add_yamldict_keyval_in(steps_i, step_key, arg_keyval)
        return steps_i

    # Add an explicit return statement here so mypy doesn't complain.
    # (mypy's static analysis cannot determine that this is dead code.)
    return steps_i