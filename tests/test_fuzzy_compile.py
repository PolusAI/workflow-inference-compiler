from datetime import timedelta
from pathlib import Path
import unittest

import graphviz
from hypothesis import given, settings, HealthCheck
import networkx as nx
import pytest

import wic
import wic.ast
import wic.cli
import wic.main
import wic.schemas
import wic.schemas.wic_schema
import wic.utils
from wic.wic_types import GraphData, GraphReps, Yaml, YamlTree, StepId

from .test_setup import get_args, tools_cwl, yml_paths, validator, wic_strategy


class TestFuzzyCompile(unittest.TestCase):

    @pytest.mark.slow
    @given(wic_strategy)
    @settings(max_examples=100,
              suppress_health_check=[HealthCheck.too_slow,
                                     HealthCheck.filter_too_much],
              deadline=timedelta(milliseconds=5000))
    # TODO: Improve schema so we can remove the health checks
    def test_fuzzy_compile(self, yml: Yaml) -> None:
        """Tests that the compiler doesn't crash when given random allegedly valid input.\n
        Note that the full schema has performance limitations, so a random subset of\n
        wic_main_schema is chosen when hypothesis=True, then random values are generated.

        Args:
            yml (Yaml): Yaml input, randomly generated according to a random subset of wic_main_schema
        """
        plugin_ns = 'global'
        y_t = YamlTree(StepId('random_stepid', plugin_ns), yml)
        yaml_tree_raw = wic.ast.read_ast_from_disk(y_t, yml_paths, tools_cwl, validator)
        yaml_tree = wic.ast.merge_yml_trees(yaml_tree_raw, {}, tools_cwl)

        yml_path = Path('random_stepid')
        graph_gv = graphviz.Digraph(name=f'cluster_{yml_path}')
        graph_gv.attr(newrank='True')
        graph_nx = nx.DiGraph()
        graphdata = GraphData(str(yml_path))
        graph = GraphReps(graph_gv, graph_nx, graphdata)
        try:
            compiler_info = wic.compiler.compile_workflow(yaml_tree, get_args(str(yml_path)), [], [graph], {}, {}, {},
                                                          {}, tools_cwl, True, relative_run_path=True, testing=True)
        except Exception as e:
            multi_def_str = 'Error! Multiple definitions of &'
            if multi_def_str in str(e):
                pass
            else:
                #import yaml
                #print(yaml.dump(yml))
                raise e


if __name__ == '__main__':
    unittest.main()
