import argparse
from pathlib import Path
from typing import List, Tuple

from . import utils
from .wic_types import KV, Tools, WorkflowInputs, InternalOutputs, Graph


def perform_edge_inference(args: argparse.Namespace,
                           tools: Tools, steps_keys: List[str],
                           yaml_stem: str,
                           i: int,
                           steps_i: KV,
                           arg_key: str,
                           graph: Graph,
                           is_root: bool,
                           namespaces: List[str],
                           vars_workflow_output_internal: InternalOutputs,
                           inputs_workflow: WorkflowInputs,
                           in_name_in_inputs_file_workflow: bool) -> KV:
    # Use in_name_in_inputs_file_workflow: bool so that at the call site, we don't
    # have to question whether or not this function modifies inputs_file_workflow
    # TODO: Figure out something better than replace. Use type and/or edam format info? Yes.
    arg_key_no_namespace = arg_key.split('___')[-1]
    arg_key_noinput = arg_key_no_namespace.replace('input_', '').replace('solute_', '').replace('energy_', 'edr_').replace('structure_', 'tpr_').replace('traj_', 'trr_')
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
    in_format = None
    if 'format' in in_tool[arg_key]:
        in_format = in_tool[arg_key]['format']
        in_dict['format'] = in_format
    for j in range(0, i)[::-1]:  # Reverse order!
        tool_j = tools[Path(steps_keys[j]).stem][1]
        out_tool = tool_j['outputs']
        out_keys = list(tool_j['outputs'])[::-1] # Reverse order!
        for out_key in out_keys:
            out_type = out_tool[out_key]['type']
            out_dict = {'type': out_type}
            out_format = None
            if 'format' in out_tool[out_key]:
                out_format = out_tool[out_key]['format']
                out_dict['format'] = out_format
            out_key_no_namespace = out_key.split('___')[-1]
            if arg_key_noinput == out_key_no_namespace.replace('output_', ''):
                if in_format and out_format:
                    assert out_format in in_format
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
                label = out_key_no_namespace if tool_j['class'] == 'Workflow' else out_key
                utils.add_graph_edge(args, graph, nss1, nss2, label)

                arg_keyval = {arg_key: arg_val}
                steps_i = utils.add_yamldict_keyval_in(steps_i, step_key, arg_keyval)
                return steps_i  # Short circuit

    match = False
    if not match:
        in_name = f'{step_name_i}___{arg_key}'  # Use triple underscore for namespacing so we can split later # {step_name_i}_input___{arg_key}

        # This just means we need to defer to the parent workflow.
        # There will actually be an error only if no parent supplies an
        # input value and thus there is no entry in inputs_file_workflow.
        if is_root and not in_name_in_inputs_file_workflow:
            print('Error! No match found for input', i + 1, step_key, arg_key)

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