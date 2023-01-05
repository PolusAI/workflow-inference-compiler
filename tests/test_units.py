import os # pylint: disable=too-many-lines
from pathlib import Path
import random
import copy
from string import printable
from typing import Dict,List,Any,Tuple, Generator
import graphviz
import networkx as nx
import unittest




from hypothesis import given, settings,HealthCheck
from hypothesis.strategies import text, integers, booleans, floats, none
import hypothesis.strategies as s
from hypothesis.strategies._internal.strategies import SearchStrategy
import hypothesis_jsonschema as hj


import wic
import wic.cli
import wic.main
import wic.schemas
import wic.schemas.wic_schema
from wic import utils, utils_cwl, labshare, compiler, ast
from wic.wic_types import (RoseTree, StepId, Yaml, YamlForest, YamlTree, GraphData, GraphReps)
from .test_setup import wic_strategy, get_args, tools_cwl


dict_or_list=hj.from_schema({}).filter(lambda x: isinstance(x, Dict) or isinstance(x,List))
dict_or_str=hj.from_schema({}).filter(lambda x: isinstance(x, Dict) or isinstance(x,str))
filtered_pretty_text: SearchStrategy[str] = text(printable,min_size=1).filter(lambda x: '_' not in x)
any_types: Any = none() | booleans() | integers() | floats() | filtered_pretty_text
most_types: Any = booleans() | integers() | floats() | filtered_pretty_text
nested_lists: Any = s.recursive(s.lists(any_types,min_size=1),\
     lambda children: s.lists(children,min_size=1))
# _ reserved for flattening dictionary keys
filtered_text: SearchStrategy[str] = s.text(min_size=1).filter(lambda x: '_' not in x) 
recursive_fn: Any = lambda children: s.dictionaries(text(min_size=1),children,min_size=1)
nested_dict: Any = hj.from_schema({}).filter(lambda x: isinstance(x, Dict))
# nested_dict: Any = s.recursive(s.dictionaries(filtered_text,any_types).filter(lambda x: len(x)!=0),recursive_fn )
str_list_nonempty = s.lists(s.text(min_size=2),min_size=1,max_size=100)
str_list = s.lists(s.text(),max_size=100)
wic_step = s.tuples(s.integers(min_value=0,max_value=100000),\
            s.text(min_size=1).filter(lambda x: ',' not in x))
wic_steps = s.lists(wic_step,min_size=1,max_size=10)
wic_step_with_name = s.tuples(s.text(min_size=1).filter(lambda x: ',' not in x),\
    s.integers(min_value=0,max_value=100000), s.text(min_size=1).filter(lambda x: ',' not in x))
wic_steps_with_name = s.lists(wic_step_with_name,min_size=1,max_size=10)

def create_compiled_rose_tree_from_nested_dict(test_dict: Dict) -> RoseTree:
    """Generate RoseTree with fake compiled work flow structure

    Args:
        test_dict (Dict): Nested Dictionary

    Returns:
        RoseTree: Rose Tree associated with nested dictionary
    """
    yaml_tree = create_yaml_tree (test_dict)
    sub_name = 'Fake Root'
    yml_path = ''
    graph_fakeroot_gv = graphviz.Digraph(name=f'cluster_{sub_name}')
    graph_fakeroot_gv.attr(newrank='True')
    graph_fakeroot_nx = nx.DiGraph()
    graphdata_fakeroot = GraphData(str(sub_name))
    graph_fakeroot = GraphReps(graph_fakeroot_gv, graph_fakeroot_nx, graphdata_fakeroot)
    fake_root = True
    compiler_info_fakeroot = wic.compiler.compile_workflow(yaml_tree, get_args(str(yml_path)),
        [], [graph_fakeroot], {}, {}, {}, {}, tools_cwl, fake_root, relative_run_path=False, testing=True)
    return compiler_info_fakeroot.rose

def check_if_duplicate_keys(first_dict: Dict[str,Any],second_dict: Dict[str,Any]) -> bool:
    """Check for duplicate keys between two dictionaries

    Args:
        first_dict (Dict[str,Any]): Dictionary of keys
        second_dict (Dict[str,Any]): Another dictionary of keys

    Returns:
        bool: If duplicate key found
    """
    has_dup=False
    for key,value in first_dict.items():
        if key in second_dict.keys():
            has_dup=True
    return has_dup


def add_wic_steps(wic_dict: Dict,step_tuple_list: List[Tuple[int,str]], \
    in_steps_dict: Dict = {}, out_steps: List[str] = []) -> Dict:
    """Add steps into nested dictionary to simulate run time
    dictionary if wic_strategy doesnt have in generated dictionary
        need this because we chose not to include the property wic in our schema

    Args:
        wic_dict (Dict): Nested dictionary

    Returns:
        Dict: Modified dictionary
    """
    plugin_ns = wic_dict.get('namespace', 'global')
    if 'steps' not in wic_dict.keys():
        steps_list = wic_dict['steps']=[]
        for i, item in enumerate(step_tuple_list):
            backend=item[1]
            wic_dict=add_backend_steps(wic_dict,backend,item[0])
            stepid = StepId(backend, plugin_ns)
            steps_list.append({stepid:None})
        wic_dict['steps']=steps_list

    if 'wic' not in wic_dict.keys():
        wic_dict['wic']={}
    inner_wic_dict=wic_dict['wic']
    if 'steps' not in inner_wic_dict.keys():
        inner_wic_dict['steps']={}
    the_wic_steps = inner_wic_dict['steps']
    if len(the_wic_steps.keys())==0:
        for item in step_tuple_list:
            if len(in_steps_dict)!=0:
                the_wic_steps[str(item)]={}
                the_wic_steps[str(item)]['in']=in_steps_dict
            if len(out_steps)!=0:
                the_wic_steps[str(item)]={}
                the_wic_steps[str(item)]['out']=out_steps
            elif len(in_steps_dict)==0 and len(out_steps)==0:
                the_wic_steps[str(item)] = ""
    inner_wic_dict['steps']=the_wic_steps
    wic_dict['wic']=inner_wic_dict
    return wic_dict


def check_if_steps_in_dict(wic_dict: Dict) -> bool:
    """Check if steps are in dictionary

    Args:
        wic_dict (Dict): Nested dictionary

    Returns:
        bool: True if found the appropriate keys
    """
    check=False
    if 'steps' in wic_dict.keys():
        if 'wic' in wic_dict.keys():
            check=True

    return check

def flatten_list(test_list: List) -> List:
    """Flatten a nested list

    Args:
        S (List): Nested list

    Returns:
        List: Flattened list
    """
    if test_list == []:
        return test_list
    if isinstance(test_list[0], list):
        return flatten_list(test_list[0]) + flatten_list(test_list[1:])
    return test_list[:1] + flatten_list(test_list[1:])

def grab_nested_dicts(d : Dict) -> Generator:

    """Recursively grab nested dicts from a dict

    Args:
        d (Dict): Recursive nested dictionary

    Yields:
        Generator[Dict]: All the sub nested dictionaries in generator object
    """
    for value in d.values():
        if isinstance(value, Dict):
            yield value
            yield from grab_nested_dicts(value)


def flatten_nested_dict(d : Dict) -> Dict:
    """Flatten a dict

    Args:
        d (Dict): Nested dictionary

    Returns:
        Dict: All items in dictionary are flattened into 1 dimension of keys
    """
    def items()-> Generator:
        for key, value in d.items():
            if isinstance(value, Dict):
                for subkey, subvalue in flatten_nested_dict(value).items():
                    yield key +'_'+subkey, subvalue
                if len(value) == 0:
                    yield key, value
            else:
                yield key, value

    return dict(items())


def flatten_nested_keys(d : Any) -> List[str]:
    """Flatten dictionary keys or list keys

    Args:
        d (Any): Nested dictionary

    Returns:
        List: All keys in dictionary are flattened into 1 dimension of keys
    """
    if isinstance(d,Dict):
        for key, value in d.items():
            if isinstance(value, Dict):
                return [key]+flatten_nested_keys(value)
            else:
                return [key]
    elif isinstance(d,List):
        return utils.flatten([flatten_nested_keys(x) for x in d])

    return []

def convert_dict_keys(d : Dict) -> Dict:
    """Convert dict keys from string to tuples

    Args:
        d (Dict): Used in conjunction with flatten_nested_dict,
        this converts keys that are stored as key1+_key2+_ .. into a tuple (key1,key2,..)

    Returns:
        Dict: Modified dictionary
    """
    new_dict = {}
    for key, value in d.items():
        if isinstance(key, str):
            key = tuple(key.split('_'))

        new_dict[tuple(key)] = value
    return new_dict




def add_backend_steps(wic_dict : Dict, backend : str, steps : int) -> Dict:
    """Add backend and steps to a wic dict

    Args:
        wic_dict (Dict): Nested dictionary, that does not have desired keys from schema
        backend (str): A string attempting to mimic how dictionary has backend key in run time
        steps (int): Step number attempting to mimic how dictionaries have step number like in run time

    Returns:
        Dict: Modified dictionary
    """
    wic_dict['backend']=backend
    wic_dict['backends']={}
    plugin_ns = wic_dict.get('namespace', 'global')
    stepid = StepId(backend, plugin_ns)
    wic_dict['backends'][stepid]={}
    wic_dict['backends'][stepid]['steps']=steps
    return wic_dict


def randomly_add_contents_lists_to_list(test : List , index : int) -> List:
    """Randomly add contents of a list to another list

    Args:
        test (List): A list of text
        index (int): Random index of list test

    Returns:
        List: Modified list with content added
    """
    newls=[]
    contentlist=test[index]
    for i,item in enumerate(test):
        if i !=index:
            templs=item.copy()
            item=random.choice(templs)
            contentlist.append(item)
            newls.append(templs)
        else:
            newls.append(contentlist)

    return newls



def recursive_tree_items(dictionary : Dict) -> Generator:
    """Recursively grab nested dicts from a dict and generate RoseTrees

    Args:
        dictionary (Dict): Nested dictionary to allow generation of RoseTree structure

    Yields:
        Generator: All RoseTrees from nested dictionaries in original dictionary
    """
    for key, value in dictionary.items():
        if isinstance(value,dict):
            yield RoseTree(data=value,sub_trees=[])
            yield from recursive_tree_items(value)

def recursive_forest_items(yaml_tree : YamlTree) -> Generator:
    """Recursively grab nested dicts from a dict and generate YamlForests

    Args:
        yaml_tree (YamlTree): YamlTree generated from nested dictionary

    Yields:
        Generator: Yields all the sub forests generated from nested dictionaries in original YamlTree
    """
    stored_dicts=[]
    for key, value in yaml_tree.yml.items():
        if isinstance(value,dict):
            yml_tree=create_yaml_tree(value)
            forest=YamlForest(yaml_tree=yml_tree,sub_forests=[])
            backend=''
            testdict=yml_tree.yml
            if 'default_backend' in testdict:
                backend = testdict['default_backend']
            if 'backend' in testdict:
                backend = testdict['backend']
            plugin_ns = testdict.get('namespace', 'global')
            stepid = StepId(backend, plugin_ns)
            tup=(stepid,forest)
            if testdict not in stored_dicts and len(testdict)!=0:
                yield tup
                stored_dicts.append(testdict) # dont want duplicates
            yield from recursive_forest_items(yml_tree)

def create_rose_tree_from_nested_dict(test_dict : Dict) ->RoseTree:
    """Create a RoseTree from a nested dict

    Args:
        test_dict (Dict): Nested dictionary that will be used to create a RoseTree

    Returns:
        RoseTree: Generated from nested dictionary
    """
    subdicts=recursive_tree_items(test_dict)
    rosetree=RoseTree(data=test_dict,sub_trees=list(subdicts))
    return rosetree

def create_yaml_tree(test_dict: Dict) ->YamlTree:
    """Create a YamlTree from a nested dict

    Args:
        test_dict (Dict): Nested dictionary that will be used to create a YamlTree

    Returns:
        YamlTree: Generated from nested dictionary
    """
    backend=''
    if 'default_backend' in test_dict:
        backend = test_dict['default_backend']
    if 'backend' in test_dict:
        backend = test_dict['backend']
    plugin_ns = test_dict.get('namespace', 'global')
    stepid = StepId(backend, plugin_ns)
    yaml_tree=YamlTree(step_id=stepid,yml=test_dict)
    return yaml_tree

def create_yaml_forest(test_yaml_tree: YamlTree) -> YamlForest:
    """Create a YamlForest from a YamlTree

    Args:
        test_yaml_tree (YamlTree): YamlTree (originally made from
        nested dictionary) that will be used to create YamlForest

    Returns:
        YamlForest: Generated from YamlTree
    """
    subfors=list(recursive_forest_items(test_yaml_tree))
    yaml_forest=YamlForest(yaml_tree=test_yaml_tree,sub_forests=subfors)
    return yaml_forest




class TestUnits(unittest.TestCase):

    def test_read_lines_pairs(self)-> None: # generalize to making files on the fly
        """Test read_lines_pairs
        """
        testpath=os.path.abspath(os.path.join(__file__,os.pardir))
        wicpath=os.path.abspath(os.path.join(testpath,os.pardir))
        srcpath=os.path.join(wicpath,'src')
        wicsrcpath=os.path.join(srcpath,'wic')
        cwl_dirs_path=Path(os.path.join(wicsrcpath,'cwl_dirs.txt'))
        self.assertEqual([('global', 'biobb/'), ('global', 'cwl_adapters/')],utils.read_lines_pairs(cwl_dirs_path))

    @given(item1=filtered_text,item2=s.integers(min_value=0,max_value=10 ** 6),item3=filtered_text)
    def test_parse_step_name_str(self,item1 : str, item2 : int, item3 : str)-> None:
        """Test parse_step_name_str

        Args:
            item1 (str): Random string without _
            item2 (int): Random integer >= 0
            item3 (str): Random string without _
        """
        output=utils.step_name_str(item1, item2 ,item3)
        self.assertEqual((item1, item2, item3),utils.parse_step_name_str(output))

    @given(item1=filtered_text,item2=filtered_text)
    def test_shorten_restore_namespaced_output_name(self,item1 : str,item2 : str)->None:
        """Test shorten_restore_namespaced_output_name

        Args:
            item1 (str): Random string without _
            item2 (str): Random string without _
        """
        thesep='__'
        namespaced_output_name=item1+thesep+item2
        tup=utils.shorten_namespaced_output_name(namespaced_output_name,sep=thesep)
        restored_name=utils.restore_namespaced_output_name(tup[0],tup[1],sep=thesep)
        self.assertEqual(restored_name,namespaced_output_name)

    

    @given(common = str_list, child_1 = str_list, child_2 = str_list)
    def test_partition_by_lowest_common_ancestor(self, common: List[str], child_1: List[str], child_2: List[str]) -> None:
        path_1 = common.copy() + child_1
        path_2 = common.copy() + child_2
        init_1, tails_1 = utils.partition_by_lowest_common_ancestor(path_1, path_2)
        init_2, tails_2 = utils.partition_by_lowest_common_ancestor(path_2, path_1)
        # Make sure that common is an initial segment of init_1 and init_2
        assert all([x==y for x,y in zip(common,init_1)])
        assert all([x==y for x,y in zip(common,init_2)])
        # Make sure that tails is an initial segment of the children.
        # We want to pair the indices starting from the end, hence reverse
        assert all([x==y for x,y in zip(tails_1[::-1],child_1[::-1])])
        assert all([x==y for x,y in zip(tails_2[::-1],child_2[::-1])])



    @settings(
    suppress_health_check=[HealthCheck.too_slow,HealthCheck.filter_too_much,HealthCheck.data_too_large],deadline=1000
    )
    @given(dict_list=s.lists(nested_dict,min_size=1,max_size=10))
    def test_get_steps_keys(self,dict_list : List[Dict])->None:
        """Test get_steps_keys

        Args:
            dict_list (List[Dict]): Listed of nested dictionaries
        """
        firstkeylist=[]
        for unflattendict in dict_list:
            firstkeylist.extend(list(unflattendict.keys()))
        self.assertEqual(firstkeylist,utils.get_steps_keys(dict_list))

    @settings(
    suppress_health_check=[HealthCheck.too_slow,HealthCheck.filter_too_much,HealthCheck.data_too_large],deadline=1000
    )
    @given(yaml_tree_tuple=\
        wic_strategy.map(create_yaml_tree),stringlist=str_list_nonempty,astringlist=str_list_nonempty) # type: ignore [arg-type]
    def test_get_subkeys(self,yaml_tree_tuple : YamlTree,stringlist : List[str], astringlist : List[str])->None:

        """Test get_subkeys

        Args:
            yaml_tree_tuple (YamlTree): YamlTree generated from nested dictionary
            stringlist (List[str]): List of random strings
            astringlist (List[str]): List of random strings
        """
        args = get_args()
        (step_id, yaml_tree) = yaml_tree_tuple
        if check_if_steps_in_dict(yaml_tree) is True:
            steps: List[Yaml] = yaml_tree['steps']
            the_wic_steps = yaml_tree['wic'].get('steps', {})
            steps_keys = utils.get_steps_keys(steps)
            tools = wic.main.get_tools_cwl(args.cwl_dirs_file)
            tools_stems = [stepid.stem for stepid in tools]

        else:
            output=randomly_add_contents_lists_to_list([stringlist,astringlist],0)
            steps_keys=output[0]
            tools_stems=output[1]




        results=utils.get_subkeys(steps_keys,tools_stems)
        expected_results=[i for i in steps_keys if i not in tools_stems]
        self.assertEqual(expected_results,results)

    @settings(
    suppress_health_check=[HealthCheck.too_slow,HealthCheck.filter_too_much,HealthCheck.data_too_large],deadline=1000
    )
    @given(
        yaml_tree=wic_strategy,
        wic_dict=wic_strategy,
        yaml_path=s.text(),
        backend=s.text(min_size=1),
        steps=s.integers(min_value=0)
    )
    def test_extract_backend(self, yaml_tree : Dict, wic_dict : Dict, yaml_path : str ,backend : str,steps : int)->None: # pylint: disable=too-many-arguments
        """Test extract_backend

        Args:
            yaml_tree (Dict): Random dictionary with wic schema
            wic_dict (Dict): Random dictionary with wic schema
            yaml_path (str): Random string
            backend (str): Random string
            steps (int): Random integer >= 0
        """
        yaml_tree_copy=copy.deepcopy(yaml_tree)
        if 'backend' in wic_dict.keys():
            backend=wic_dict['backend']
        else:
            wic_dict=add_backend_steps(wic_dict,backend,steps)
            yaml_tree_copy.update({'steps': steps})
        results=utils.extract_backend(yaml_tree,wic_dict,Path(yaml_path))
        self.assertEqual((backend,yaml_tree_copy),results)

    @given(masterlist=nested_lists)
    def test_flatten(self, masterlist : List[List])->None:
        """Test flatten

        Args:
            masterlist (List[List]): Nested lists of random data types
        """
        if any(isinstance(el, list) for el in masterlist):
            expected_output=[x for lst in masterlist for x in lst]
            self.assertEqual(expected_output,utils.flatten(masterlist))

    @settings(
    suppress_health_check=[HealthCheck.too_slow,HealthCheck.filter_too_much,HealthCheck.data_too_large],deadline=1000
    )
    @given(rosetree=nested_dict.map(create_rose_tree_from_nested_dict))
    def test_flatten_rose_tree(self, rosetree : RoseTree) -> None:
        """Test flatten_rose_tree

        Args:
            rosetree (RoseTree): RoseTree generated from random nested dictionary
        """
        def flatten_rose_tree(rosetree : RoseTree) -> List:
            data=rosetree.data
            flattened_dict=list(grab_nested_dicts(data))
            return [data]+flattened_dict

        expected_results=flatten_rose_tree(rosetree)
        results=utils.flatten_rose_tree(rosetree)
        self.assertEqual(expected_results,results)

    @given(test_obj = dict_or_list)
    def test_recursively_dict_key(self,test_obj: Any) -> None:
        """

        Args:
            test_obj (Any): _description_
        """
        flatten_keylist=flatten_nested_keys(test_obj)
        if len(flatten_keylist)!=0:
            key=random.choice(flatten_keylist) # choose random key to delete
            results = utils.recursively_delete_dict_key(key, test_obj)
            assert not utils.recursively_contains_dict_key(key,results)


    @settings(
    suppress_health_check=[HealthCheck.too_slow,HealthCheck.filter_too_much,HealthCheck.data_too_large],deadline=1000
    )
    @given(
        wic_dict=wic_strategy,
        step_tuple_list=wic_steps,
    )
    def test_reindex_wic_steps(self,wic_dict: Dict,step_tuple_list: List[Tuple[int,str]]) -> None:
        """Test reindex_wic_steps

        Args:a
            wic_dict (Dict): Dictionary with wic_strategy
        """
        has_steps=check_if_steps_in_dict(wic_dict)
        if has_steps is False:
            wic_dict = add_wic_steps(wic_dict,step_tuple_list)
        the_wic_steps = wic_dict['wic'].get('steps', {})
        random_step=random.choice(list(the_wic_steps.keys()))
        random_step=utils.parse_int_string_tuple(random_step)
        random_index=random_step[0]
        results=utils.reindex_wic_steps(the_wic_steps,random_index)
        new_wic_steps={}
        for key,value in the_wic_steps.items():
            index,the_string=utils.parse_int_string_tuple(key)
            newstr = f'({index+1}, {the_string})' if index >= random_index else key
            new_wic_steps[newstr]=value
        expected_results=new_wic_steps
        self.assertEqual(expected_results,results)


    @given(common = str_list, child_1 = str_list, child_2 = str_list)
    def test_recursively_insert_into_dict_tree(self, common: List[str], child_1: List[str], child_2: List[str]) -> None:
        """ Test recursively_insert_into_dict_tree

        Args:
            common (List[str]): Common root
            child_1 (List[str]): Path after common root
            child_2 (List[str]): Path after common root
        """
        tails_12 = utils.partition_by_lowest_common_ancestor(child_1,child_2)[1]
        tails_21 = utils.partition_by_lowest_common_ancestor(child_2,child_1)[1]
        disjoint_leaves = len(tails_12) == 0 and len(tails_21) == 0
        if disjoint_leaves:
            path_1 = common.copy() + child_1
            path_2 = common.copy() + child_2
            dummy_value_1 = 0
            dummy_value_2 = 1

            tree_1: Dict = {}
            tree_1 = utils.recursively_insert_into_dict_tree(tree_1, path_1, dummy_value_1)
            tree_1 = utils.recursively_insert_into_dict_tree(tree_1, path_2, dummy_value_2)

            tree_2: Dict = {}
            tree_2 = utils.recursively_insert_into_dict_tree(tree_2, path_2, dummy_value_2)
            tree_2 = utils.recursively_insert_into_dict_tree(tree_2, path_1, dummy_value_1)

            if path_1 != path_2:
                self.assertEqual(tree_1, tree_2)


    @given(input_mapping = s.dictionaries(text(min_size=1),\
        s.lists(text(min_size=1),min_size=1),min_size=1).filter(lambda x: len(x)!=0))
    def test_get_input_mappings(self,input_mapping: Dict[str,List[str]]) -> None: # may need to change inputs for more general testing
        """Test get_input_mappings

        Args:
            input_mapping (Dict): Dict[str,List[str]])
        """
        found_dup=False
        all_values=[]
        for key,value_list in input_mapping.items():
            all_values.extend(value_list)
            if key in value_list:
                found_dup=True
        if found_dup is False:
            key_list=list(input_mapping.keys())
            k=random.randint(1,len(key_list))
            random_subset=random.sample(key_list,k)
            arg_keys=random_subset
            results=utils.get_input_mappings(input_mapping,arg_keys,False)
            for item in results:
                assert item in all_values



    @given(output_mapping = s.dictionaries(text(min_size=1),text(min_size=1),min_size=1))
    def test_get_output_mappings(self,output_mapping: Dict[str,str]) -> None: # may need to change inputs for more general testing
        """Test get_output_mappings

        Args:
            output_mapping (Dict): Dict[str,str])
        """
        is_dup=False
        for key,value in output_mapping.items():
            if key==value:
                is_dup=True
        if is_dup is False:
            out_key=random.choice(list(output_mapping.keys()))
            results=utils.get_output_mapping(output_mapping,out_key)
            all_vals=list(output_mapping.values())
            assert results in all_vals

    @given(wic_dict = wic_strategy, keyval=nested_dict, inkeyval=nested_dict,
    step_tuple_list=wic_steps)
    @settings(
    suppress_health_check=[HealthCheck.too_slow,HealthCheck.filter_too_much,HealthCheck.data_too_large],deadline=1000
    )
    def test_add_yamldict_keyval_in(self,wic_dict: Dict, keyval: Dict, inkeyval: Dict, \
        step_tuple_list: List[Tuple[int,str]]) -> None:
        """Test add_yamldict_keyval_in

        Args:
            wic_dict (Dict): wic_strategy that has steps from work flow
            keyval (Dict): Extra steps added into workflow
            inkeyval (Dict): Steps added into workflow for inputs
        """
        in_steps=inkeyval
        if check_if_duplicate_keys(in_steps,keyval) is False:
            wic_dict = add_wic_steps(wic_dict,step_tuple_list,in_steps_dict=in_steps)
            steps_i = wic_dict['wic'].get('steps', {})
            key_list=list(steps_i.keys())
            step_key=random.choice(key_list)
            results=utils_cwl.add_yamldict_keyval_in(steps_i,step_key,keyval)
            result_step_key=results[step_key]
            result_step_key_in=result_step_key['in']
            result_step_key_in_flat=flatten_nested_dict(result_step_key_in)
            result_step_key_in_flat_values=list(result_step_key_in_flat.values())
            keyval_flat=flatten_nested_dict(keyval)
            keyval_flat_values=list(keyval_flat.values())
            for value in keyval_flat_values:
                assert value in result_step_key_in_flat_values



    @settings(
    suppress_health_check=[HealthCheck.too_slow,HealthCheck.filter_too_much,HealthCheck.data_too_large],deadline=1000
    )
    @given(wic_dict = wic_strategy, keyval=s.lists(s.text(min_size=2)), outvals=s.lists(s.text(min_size=1)),
    step_tuple_list=wic_steps)
    def test_add_yamldict_keyval_out(self,wic_dict: Dict, keyval: List[str], outvals: List[str], \
        step_tuple_list: List[Tuple[int,str]]) -> None:
        """Test add_yamldict_keyval_out

        Args:
            wic_dict (Dict): wic_strategy that has steps from work flow
            keyval (List[str]): Extra steps added into workflow for outputs
            outvals (List[str]): Steps added into workflow for outputs
        """

        the_out_steps=outvals
        wic_dict = add_wic_steps(wic_dict,step_tuple_list,out_steps=the_out_steps)
        steps_i = wic_dict['wic'].get('steps', {})
        key_list=list(steps_i.keys())
        step_key=random.choice(key_list)
        results=utils_cwl.add_yamldict_keyval_out(steps_i,step_key,keyval)
        result_step_key=results[step_key]
        result_step_key_out=result_step_key['out']
        for value in keyval:
            assert value in result_step_key_out

    @settings(
    suppress_health_check=[HealthCheck.too_slow,HealthCheck.filter_too_much,HealthCheck.data_too_large],deadline=1000
    )
    @given(test_obj = dict_or_str, items_dict= nested_dict)
    def test_canonicalize_type(self,test_obj: Any,  items_dict: Dict) -> None:
        """Test canonicalize_type

        Args:
            test_obj (Any): Nested dictionary or string
        """
        if isinstance(test_obj,Dict):
            test_obj.update({'type':'array'})
            test_obj.update({'items':items_dict})
            results=utils_cwl.canonicalize_type(test_obj)
            flat_results_dict=flatten_nested_dict(results)
            flat_test_dict=flatten_nested_dict(test_obj)
            self.assertEqual(flat_results_dict,flat_test_dict)
        elif isinstance(test_obj,str):
            test_list=['?','[]']
            item=random.choice(test_list)
            test_obj+=item
            results=utils_cwl.canonicalize_type(test_obj)
            if item=='?':
                the_item=results[1]
                assert the_item in test_obj
            elif item=='[]':
                the_item=results['items']
                assert the_item in test_obj

    @settings(
    suppress_health_check=[HealthCheck.too_slow,HealthCheck.filter_too_much,HealthCheck.data_too_large],deadline=1000
    )
    @given(list_nested_dicts = s.lists(nested_dict,min_size=3,max_size=3))
    def test_copy_cwl_input_output_dict(self,list_nested_dicts: List[Dict]) -> None:
        """Test copy_cwl_input_output_dict

        Args:
            list_nested_dics (List[Dict]): List of nested dictionaries
        """
        keylist=['format', 'label', 'doc']
        io_dict=dict(zip(keylist,list_nested_dicts))
        io_dict.update({'type':{}})
        results=utils_cwl.copy_cwl_input_output_dict(io_dict)
        for key,value in io_dict.items():
            if key in keylist:
                self.assertEqual(value,results[key])

    @settings(
    suppress_health_check=[HealthCheck.too_slow,HealthCheck.filter_too_much,HealthCheck.data_too_large],deadline=1000
    )
    @given(test_dict = nested_dict)
    def test_remove_dot_dollar(self,test_dict: Dict) -> None: # Claimed may not be correct (but passes test)
        """Test remove_dot_dollar

        Args:
            test_dict (Dict): Nested dictionary
        """
        key_list=['$namespaces','$schemas','.yml']
        for key in key_list:
            test_dict.update({key:''})

        results=labshare.remove_dot_dollar(test_dict)
        for key in key_list:
            assert key not in results.keys()

    @settings(
    suppress_health_check=[HealthCheck.too_slow,HealthCheck.filter_too_much,HealthCheck.data_too_large]
    )
    @given(wic_dict = wic_strategy, step_tuple_list=wic_steps, single_step = wic_step)
    def test_insert_step_into_workflow(self, wic_dict: Dict, step_tuple_list: \
        List[Tuple[int,str]],single_step: Tuple[int,str]) -> None:
        """Test insert_step_into_workflow

        Args:
            wic_dict (Dict): wic_strategy
            step_tuple_list (List[Tuple[int,str]]): List of steps to initialize wic_dict
            single_step (Tuple[int,str]): Extra step to add for testing function
        """
        has_steps=check_if_steps_in_dict(wic_dict)
        if has_steps is False:
            wic_dict = add_wic_steps(wic_dict,step_tuple_list)
        plugin_ns = wic_dict.get('namespace', 'global')
        i=single_step[0]
        backend=single_step[1]
        stepid = StepId(backend, plugin_ns)
        first_tool_key = list(tools_cwl.keys())[0]
        first_tool_value = tools_cwl[first_tool_key]
        tools = {stepid : first_tool_value} # make a fake tool
        results=compiler.insert_step_into_workflow(wic_dict,stepid,tools,i)
        results_wic=results['wic']
        results_wic_steps=results_wic['steps']
        found=False
        for key,value in results_wic_steps.items():
            step=utils.parse_int_string_tuple(key)
            index=step[0]
            if index==i+1:
                found=True
        self.assertEqual(found,True)

    @settings(
    suppress_health_check=[HealthCheck.too_slow,HealthCheck.filter_too_much,HealthCheck.data_too_large],deadline=1000
    )
    @given(yaml_tree_tuple=wic_strategy.map(create_yaml_tree)) # type: ignore [arg-type]
    def test_tree_to_forest(self,yaml_tree_tuple: YamlTree) -> None:
        """Test tree_to_forest

        Args:
            yaml_tree_tuple (Dict): _description_
        """
        input_yml=yaml_tree_tuple.yml
        results=ast.tree_to_forest(yaml_tree_tuple,tools_cwl)
        results_tree=results.yaml_tree
        results_tree_yml=results_tree.yml
        assert input_yml.items() <= results_tree_yml.items()
