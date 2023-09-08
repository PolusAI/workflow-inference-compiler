import subprocess as sub
from pathlib import Path
import signal
import sys
from typing import List

import graphviz
import networkx as nx
import pytest
import yaml
from networkx.algorithms import isomorphism

import wic.cli
import wic.compiler
import wic.run_local
import wic.utils
from wic import auto_gen_header
from wic.wic_types import GraphData, GraphReps, NodeData, StepId, Yaml, YamlTree

from .test_setup import get_args, tools_cwl, yml_paths, validator, yml_paths_tuples


# Due to the computational complexity of the graph isomorphism problem, we
# need to manually exclude large workflows.
# See https://en.wikipedia.org/wiki/Graph_isomorphism_problem
large_workflows = ['dsb', 'dsb1', 'elm', 'vs_demo_2', 'vs_demo_3', 'vs_demo_4']
yml_paths_tuples_not_large = [(s, p) for (s, p) in yml_paths_tuples if s not in large_workflows]

# NOTE: Most of the workflows in this list have free variables because they are subworkflows
# i.e. if you try to run them, you will get "Missing required input parameter"
run_blacklist: List[str] = [
    'assign_partial_charges_batch',
    'convert_ligand_mol2_to_pdbqt_mdanalysis',
    'download_smiles_ligand_db',
    'convert_ligand_mol2_to_pdbqt_obabel',
    'analysis_realtime_ligand',
    'analysis_realtime_complex',
    'analysis_realtime_protein',
    'ligand_modeling_docking',
    'align_protein_CA_pymol',
    'assign_partial_charges',
    'minimize_ligand_only',
    'analysis_final_steps',
    'autodock_vina_rescore',
    'analysis_final',
    'gen_topol_params',
    'analysis_realtime',
    'convert_pdbqt',
    'download_pdb',
    'setup_vac_min',
    'npt_gromacs',
    'setup_pdb',
    'docking_stability',
    'npt_amber',
    'analysis',
    'solv_ion',
    'topology',
    'stability',
    'docking',
    'l-bfgs',
    'basic',
    'equil',
    'setup',
    'steep',
    'modeling',  # called in tutorial
    'nmr',  # Do NOT run nmr, because it makes the CI go from ~30 mins to >6 hours due to:
    # "MDAnalysis/coordinates/XDR.py:202: UserWarning:
    # Cannot write lock/offset file in same location as
    # .../prod.trr. Using slow offset calculation."
    # So let's un-blacklist tutorial
    # 'tutorial',  # called in nmr
    'prod',
    'flc',
    'dsb',
    'npt',
    'nvt',
    'min',
    'cg',
    'yank',
    'fix_protein',
    # These (currently) always return success, so no point in running them.
    'cwl_watcher_analysis',
    'cwl_watcher_complex',
    'cwl_watcher_ligand',
    'cwl_watcher_protein',
    'run_diffdock',
    'diffdock_workflow'
]


yml_paths_tuples_not_blacklist = [(s, p) for (s, p) in yml_paths_tuples if s not in run_blacklist]
# currently [vs_demo_2, vs_demo_3, vs_demo_4, elm, nmr,
#            multistep1, multistep2, multistep3, helloworld, scattering_scaling]


def is_isomorphic_with_timeout(g_m: isomorphism.GraphMatcher, yml_path_str: str) -> None:
    """Calls the .is_isomorphic() method with a timeout of 10 seconds.

    Args:
        gm (isomorphism.GraphMatcher): The graph isomorphism object.
        yml_path_str (str): The root yml workflow file for error reporting.
    """
    def handler(_signum, _frame):  # type: ignore
        line1 = 'Graph isomorphism check timed out after 10 seconds.'
        line2 = f'Consider adding {yml_path_str} to the large_workflows list.'
        raise Exception(f'{line1}\n{line2}')

    # See https://github.com/bokeh/bokeh/issues/11627#issuecomment-921576787
    if sys.platform == 'win32':
        # Windows does not support alarm, so just use a regular call here.
        # Note that there is a backup timeout in the github action workflow
        assert g_m.is_isomorphic()  # See top-level comment above!
    else:
        # See https://docs.python.org/3/library/signal.html#examples
        # NOTE: You CANNOT use `pytest --workers 8 ...` with this. Otherwise:
        # "ValueError: signal only works in main thread of the main interpreter"
        signal.signal(signal.SIGALRM, handler)

        signal.alarm(10)  # timeout after 10 seconds
        assert g_m.is_isomorphic()  # See top-level comment above!
        signal.alarm(0)  # Disable the alarm


def get_graph_reps(name: str) -> GraphReps:
    """Initialize graph representations

    Args:
        name (str): The name of the graph

    Returns:
        GraphReps: A tuple of graph representations
    """
    graph_gv = graphviz.Digraph(name=f'cluster_{name}')
    graph_gv.attr(newrank='True')
    graph_nx = nx.DiGraph()
    graphdata = GraphData(str(name))
    return GraphReps(graph_gv, graph_nx, graphdata)


@pytest.mark.slow
@pytest.mark.parametrize("yml_path_str, yml_path", yml_paths_tuples_not_blacklist)
def test_run_examples(yml_path_str: str, yml_path: Path, cwl_runner: str) -> None:
    """Runs all of the examples in the examples/ directory. Note that some of
    the yml files lack inputs and cannot be run independently, and are excluded.
    """
    if yml_path_str == 'vs_demo_4':
        return None  # Skip so we don't accidentally DOS pdbbind.org.cn

    args = get_args(str(yml_path))

    # First compile the workflow.
    # Load the high-level yaml workflow file.
    with open(yml_path, mode='r', encoding='utf-8') as y:
        root_yaml_tree: Yaml = yaml.safe_load(y.read())
    Path('autogenerated/').mkdir(parents=True, exist_ok=True)
    wic_tag = {'wic': root_yaml_tree.get('wic', {})}
    plugin_ns = wic_tag['wic'].get('namespace', 'global')
    step_id = StepId(yml_path_str, plugin_ns)
    y_t = YamlTree(step_id, root_yaml_tree)
    yaml_tree_raw = wic.ast.read_ast_from_disk(args.homedir, y_t, yml_paths, tools_cwl, validator)
    with open(f'autogenerated/{Path(yml_path).stem}_tree_raw.yml', mode='w', encoding='utf-8') as f:
        f.write(yaml.dump(yaml_tree_raw.yml))
    yaml_tree = wic.ast.merge_yml_trees(yaml_tree_raw, {}, tools_cwl)
    with open(f'autogenerated/{Path(yml_path).stem}_tree_merged.yml', mode='w', encoding='utf-8') as f:
        f.write(yaml.dump(yaml_tree.yml))

    graph = get_graph_reps(str(yml_path))
    compiler_info = wic.compiler.compile_workflow(yaml_tree, args, [], [graph], {}, {}, {}, {},
                                                  tools_cwl, True, relative_run_path=True, testing=True)
    rose_tree = compiler_info.rose
    sub_node_data: NodeData = rose_tree.data
    yaml_stem = sub_node_data.name

    wic.utils.write_to_disk(rose_tree, Path('autogenerated/'), relative_run_path=True)

    # NOTE: Do not use --cachedir; we want to actually test everything.
    retval = wic.run_local.run_local(args, rose_tree, None, cwl_runner, True)
    assert retval == 0
    return


@pytest.mark.fast
@pytest.mark.serial
@pytest.mark.parametrize("yml_path_str, yml_path", yml_paths_tuples_not_large)
def test_cwl_embedding_independence(yml_path_str: str, yml_path: Path) -> None:
    """Tests that compiling a subworkflow is independent of how it is embedded
    into a parent workflow. Specifically, this compiles the root workflow and
    re-compiles every subworkflow (individually) as if it were a root workflow,
    then checks that the CWL for each subworkflow remains identical and checks
    that the embedded subworkflow DAGs and the re-compiled DAGs are isomorphic.
    """
    args = get_args(str(yml_path))

    # Load the high-level yaml workflow file.
    with open(yml_path, mode='r', encoding='utf-8') as y:
        root_yaml_tree: Yaml = yaml.safe_load(y.read())
    # Write the combined workflow (with all subworkflows as children) to disk.
    Path('autogenerated/').mkdir(parents=True, exist_ok=True)
    wic_tag = {'wic': root_yaml_tree.get('wic', {})}
    plugin_ns = wic_tag['wic'].get('namespace', 'global')
    step_id = StepId(yml_path_str + '.yml', plugin_ns)
    y_t = YamlTree(step_id, root_yaml_tree)
    yaml_tree_raw = wic.ast.read_ast_from_disk(args.homedir, y_t, yml_paths, tools_cwl, validator)
    with open(f'autogenerated/{yml_path.stem}_tree_raw.yml', mode='w', encoding='utf-8') as f:
        f.write(yaml.dump(yaml_tree_raw.yml))
    yaml_tree = wic.ast.merge_yml_trees(yaml_tree_raw, {}, tools_cwl)
    with open(f'autogenerated/{yml_path.stem}_tree_merged.yml', mode='w', encoding='utf-8') as f:
        f.write(yaml.dump(yaml_tree.yml))

    # NOTE: The entire purpose of parsing an entire yaml forest is so we
    # can easily access the subtrees here. (i.e. without re-walking the AST)
    yaml_forest = wic.ast.tree_to_forest(yaml_tree, tools_cwl)
    yaml_forest_lst = wic.utils.flatten_forest(yaml_forest)

    graph = get_graph_reps(str(yml_path))
    is_root = True
    compiler_info = wic.compiler.compile_workflow(yaml_tree, args, [], [graph], {}, {}, {}, {},
                                                  tools_cwl, is_root, relative_run_path=False, testing=True)
    rose_tree = compiler_info.rose
    node_data_lst: List[NodeData] = wic.utils.flatten_rose_tree(rose_tree)

    # This test doesn't necessarily need to write to disk, but useful for debugging.
    wic.utils.write_to_disk(rose_tree, Path('autogenerated/'), relative_run_path=False)

    # Now, for each subworkflow of the given root workflow, compile the
    # subworkflow again from scratch, as if it were the root workflow,
    # and check that the generated CWL is identical. In other words,
    # check that the generated CWL of a subworkflow is independent of its
    # embedding into a parent workflow.
    assert len(node_data_lst[1:]) == len(yaml_forest_lst)
    for sub_node_data, sub_yaml_forest in zip(node_data_lst[1:], yaml_forest_lst):
        sub_name = sub_node_data.name
        assert sub_yaml_forest.yaml_tree.step_id.stem == sub_name + '.yml'

        # NOTE: Do we want to also test embedding independence with args.graph_inline_depth?
        # If so, we will need to patch testargs depending on len(sub_node_data.namespaces)
        # (due to the various instances of `if len(namespaces) < args.graph_inline_depth`)

        graph_fakeroot = get_graph_reps(str(sub_name))
        fake_root = True
        compiler_info_fakeroot = wic.compiler.compile_workflow(sub_yaml_forest.yaml_tree, get_args(str(yml_path)),
                                                               [], [graph_fakeroot], {}, {}, {}, {}, tools_cwl, fake_root, relative_run_path=False, testing=True)
        sub_node_data_fakeroot: NodeData = compiler_info_fakeroot.rose.data
        sub_cwl_fakeroot = sub_node_data_fakeroot.compiled_cwl

        # NOTE: Relative run: paths cause this test to fail, so remove them.
        # Using namespaced filenames in a single flat directory also
        # doesn't work because the namespaces will be of different lengths.
        sub_cwl_embedded = wic.utils.recursively_delete_dict_key('run', sub_node_data.compiled_cwl)
        sub_cwl_fakeroot = wic.utils.recursively_delete_dict_key('run', sub_cwl_fakeroot)

        if sub_cwl_embedded != sub_cwl_fakeroot:
            # Before we crash and burn, write out files for debugging.
            with open(f'{sub_name}_forest_embedded.yml', mode='w', encoding='utf-8') as w:
                w.write(yaml.dump(yaml_forest))
            with open(f'{sub_name}_forest_fakeroot.yml', mode='w', encoding='utf-8') as w:
                w.write(yaml.dump(sub_yaml_forest))
            # NOTE: Use _dot_cwl so we don't glob these files in get_tools_cwl()
            yaml_content = yaml.dump(sub_cwl_embedded, sort_keys=False, line_break='\n', indent=2)
            filename_emb = f'{sub_name}_embedded_dot_cwl'
            with open(filename_emb, mode='w', encoding='utf-8') as w:
                w.write('#!/usr/bin/env cwl-runner\n')
                w.write(auto_gen_header)
                w.write(''.join(yaml_content))
            yaml_content = yaml.dump(sub_cwl_fakeroot, sort_keys=False, line_break='\n', indent=2)
            filename_fake = f'{sub_name}_fakeroot_dot_cwl'
            with open(filename_fake, mode='w', encoding='utf-8') as w:
                w.write('#!/usr/bin/env cwl-runner\n')
                w.write(auto_gen_header)
                w.write(''.join(yaml_content))
            cmd = f'diff {filename_emb} {filename_fake} > {sub_name}.diff'
            sub.run(cmd, shell=True, check=False)
            print(f'Error! Check {filename_emb} and {filename_fake} and {sub_name}.diff')
        assert sub_cwl_embedded == sub_cwl_fakeroot

        # Check that the subgraphs are isomorphic.
        sub_graph_nx = sub_node_data.graph.networkx
        sub_graph_fakeroot_nx = sub_node_data_fakeroot.graph.networkx
        # assert isomorphism.faster_could_be_isomorphic(sub_graph_nx, sub_graph_fakeroot_nx)
        g_m = isomorphism.GraphMatcher(sub_graph_nx, sub_graph_fakeroot_nx)
        print('is_isomorphic()?', yml_path_str, sub_name)
        is_isomorphic_with_timeout(g_m, yml_path_str)


@pytest.mark.serial
@pytest.mark.parametrize("yml_path_str, yml_path", yml_paths_tuples_not_large)
def test_inline_subworkflows(yml_path_str: str, yml_path: Path) -> None:
    """Tests that compiling a workflow is independent of how subworkflows are inlined.
    Specifically, this inlines every subworkflow (individually) and checks that
    the original DAG and the inlined DAGs are isomorphic.
    """
    args = get_args(str(yml_path))
    # Load the high-level yaml workflow file.
    with open(yml_path, mode='r', encoding='utf-8') as y:
        root_yaml_tree: Yaml = yaml.safe_load(y.read())
    Path('autogenerated/').mkdir(parents=True, exist_ok=True)
    wic_tag = {'wic': root_yaml_tree.get('wic', {})}
    plugin_ns = wic_tag['wic'].get('namespace', 'global')
    step_id = StepId(yml_path_str, plugin_ns)
    y_t = YamlTree(step_id, root_yaml_tree)
    yaml_tree_raw = wic.ast.read_ast_from_disk(args.homedir, y_t, yml_paths, tools_cwl, validator)
    with open(f'autogenerated/{Path(yml_path).stem}_tree_raw.yml', mode='w', encoding='utf-8') as f:
        f.write(yaml.dump(yaml_tree_raw.yml))
    yaml_tree = wic.ast.merge_yml_trees(yaml_tree_raw, {}, tools_cwl)
    with open(f'autogenerated/{Path(yml_path).stem}_tree_merged.yml', mode='w', encoding='utf-8') as f:
        f.write(yaml.dump(yaml_tree.yml))

    namespaces_list = wic.ast.get_inlineable_subworkflows(yaml_tree, tools_cwl, 'backend' in wic_tag, [])
    if namespaces_list == []:
        assert True  # There's nothing to test

    graph = get_graph_reps(str(yml_path))
    compiler_info = wic.compiler.compile_workflow(yaml_tree, args, [], [graph], {}, {}, {}, {},
                                                  tools_cwl, True, relative_run_path=True, testing=True)
    rose_tree = compiler_info.rose
    sub_node_data: NodeData = rose_tree.data

    wic.utils.write_to_disk(rose_tree, Path('autogenerated/'), relative_run_path=True)

    # Inline each subworkflow individually and check that the graphs are isomorphic.
    for namespaces in namespaces_list:
        inline_yaml_tree = wic.ast.inline_subworkflow(yaml_tree, tools_cwl, namespaces)

        inline_graph = get_graph_reps(str(yml_path))
        inline_compiler_info = wic.compiler.compile_workflow(inline_yaml_tree, args,
                                                             [], [inline_graph], {}, {}, {}, {}, tools_cwl, True, relative_run_path=True, testing=True)
        inline_rose_tree = inline_compiler_info.rose
        inline_sub_node_data: NodeData = inline_rose_tree.data

        # Check that the subgraphs are isomorphic.
        sub_graph_nx = sub_node_data.graph.networkx
        sub_graph_fakeroot_nx = inline_sub_node_data.graph.networkx
        # assert isomorphism.faster_could_be_isomorphic(sub_graph_nx, sub_graph_fakeroot_nx)
        g_m = isomorphism.GraphMatcher(sub_graph_nx, sub_graph_fakeroot_nx)
        print('is_isomorphic()?', yml_path_str, namespaces)
        is_isomorphic_with_timeout(g_m, yml_path_str)
