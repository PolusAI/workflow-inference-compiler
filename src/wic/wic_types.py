from typing import Any, Dict, List, Tuple

import graphviz
import networkx as nx

KV = Dict[str, Any]
Cwl = KV
Json = KV
Yaml = KV

Tool = Tuple[str, Cwl]
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
Graph = Tuple[DiGraph, nx.DiGraph]
YamlDSLArgs = Yaml

# Since we cannot store extra tags in CWL files, we need a data structure
# to store temporary compiler info that gets passed through the recursion.
# Unfortunately, since mypy does not support Algebraic Data Types (ADTs)
# we have to weaken the type of RecursiveData to Any :(
NodeData = Tuple[Namespaces, str, Cwl, WorkflowInputsFile, ExplicitEdgeDefs, ExplicitEdgeCalls, Graph]
# RecursiveData = Tuple[NodeData, List[RecursiveData]]
RecursiveData = Any
CompilerInfo = Tuple[RecursiveData, WorkflowInputs, WorkflowInputsFile, InternalOutputs, ExplicitEdgeDefs, ExplicitEdgeCalls, StepName1]

# Create a type for our Abstract Syntax Tree (AST).
# We can probably use Dict here if str is step_name_i not just yaml_stem.
# If we need to insert steps (i.e. file format conversions), that will happen
# after edge inference, so there should not be a uniqueness issue w.r.t. step
# number re-indexing.
# NOTE: This Dict is currently performing double-duty: if there is only one
# backend, it contains the subworkflows. If there are backends, it contains
# each backend. For now, these are mutually exclusive, so it should be okay.
YamlTree = Tuple[str, Yaml]
# YamlForest = Tuple(YamlTree, Dict[str, YamlForest])
YamlForest = Any