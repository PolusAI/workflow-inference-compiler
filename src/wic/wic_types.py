from typing import Any, Dict, List, Tuple, NamedTuple

import graphviz
import networkx as nx

# See https://mypy.readthedocs.io/en/stable/kinds_of_types.html#type-aliases

KV = Dict[str, Any]
Cwl = KV
Json = KV
Yaml = KV

# In python there are unfortunately an enormous number of ways to represent the humble struct.
# See https://stackoverflow.com/questions/53409117/what-are-the-main-differences-of-namedtuple-and-typeddict-in-python-mypy
# I have chosen to use NamedTuple to emphasize the immutability aspect (see below).
# See https://mypy.readthedocs.io/en/stable/kinds_of_types.html#named-tuples

class Tool(NamedTuple):
    run_path: str
    cwl: Cwl
Tools = Dict[str, Tool]

# NOTE: Please read the Namespacing section of docs/devguide.md !!!
Namespace = str
Namespaces = List[Namespace]

WorkflowInputs = Dict[str, Dict[str, str]]
WorkflowInputsFile = Dict[str, Dict[str, str]]
WorkflowOutputs = List[Yaml]
InternalOutputs = List[str]
ExplicitEdgeDef = Tuple[Namespaces, str]
ExplicitEdgeDefs = Dict[str, ExplicitEdgeDef]
ExplicitEdgeCalls = Dict[str, ExplicitEdgeDef]
PluginID = int
StepName1 = str
DiGraph = Any # graphviz.DiGraph
# This groups together the classes which represent our graph.
# Excluding --graph_inline_depth related code, all graph
# operations should be performed on all representations.
class GraphReps(NamedTuple):
    graphviz: DiGraph
    networkx: nx.DiGraph
YamlDSLArgs = Yaml

# Since we cannot store extra tags in CWL files, we need a data structure
# to store temporary compiler info that gets passed through the recursion.
# The number of subworkflows is arbitrary (zero or more), so what we want is a
# Rose Tree https://en.wikipedia.org/wiki/Rose_tree
# Unfortunately, since mypy does not support Algebraic Data Types (ADTs)
# we have to break the recursion by replacing the recursive instance of RoseTree with Any :(
DataType = Any
class RoseTree(NamedTuple):
    data: DataType
    sub_trees: List[Any] # Any = RoseTree
# Note that instead of DataType we could provide a specific type, but remember that
# a Rose Tree is defined by its structure, not by the specific type of data it contains.
# We can simply cast to a specific type at each call site, i.e.
# data: SpecificType = rose_tree.data

# Now we can define a specific data type for a single node in our Abstract Syntax Tree.
class NodeData(NamedTuple):
    namespaces: Namespaces
    name: str
    compiled_cwl: Cwl
    workflow_inputs_file: WorkflowInputsFile
    explicit_edge_defs: ExplicitEdgeDefs
    explicit_edge_calls: ExplicitEdgeCalls
    graph: GraphReps
    inputs_workflow: WorkflowInputs
    step_name_1: StepName1

class EnvData(NamedTuple):
    inputs_file_workflow: WorkflowInputsFile
    vars_workflow_output_internal: InternalOutputs
    explicit_edge_defs: ExplicitEdgeDefs
    explicit_edge_calls: ExplicitEdgeCalls

class CompilerInfo(NamedTuple):
    rose: RoseTree
    env: EnvData
# Note that while Tuples and NamedTuples are immutable, they can contain mutable entries!
# Let's partition the data into entries which are immutable / fixed
# (well, at least after compilation of the subworkflow is complete) and entries
# which are mutably updated throughout the recursion. The latter are essentially
# environment variables of sorts, and they do not need to be stored in the Rose Tree.

# Create a type for our Abstract Syntax Tree (AST).
# We can probably use Dict here if str is step_name_i not just yaml_stem.
# If we need to insert steps (i.e. file format conversions), that will happen
# after edge inference, so there should not be a uniqueness issue w.r.t. step
# number re-indexing.
# NOTE: This Dict is currently performing double-duty: if there is only one
# backend, it contains the subworkflows. If there are backends, it contains
# each backend. For now, these are mutually exclusive, so it should be okay.
class YamlTree(NamedTuple):
    name: str
    yml: Yaml
class YamlForest(NamedTuple):
    yaml_tree: YamlTree
    sub_forests: Dict[str, Any] # Any = YamlForest