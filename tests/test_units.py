import unittest
import pytest
import sys
import os
import yaml
import string
import random
import copy
import argparse
from unittest.mock import patch
from string import ascii_letters, printable
from typing import Dict,List,Any,Tuple, Generator, Union
import datetime
from pathlib import Path

from hypothesis import given, settings,HealthCheck
from hypothesis.strategies import text, integers, booleans, floats, dictionaries, recursive, lists, none
import hypothesis.strategies as s
from hypothesis_jsonschema import from_schema
from hypothesis.strategies._internal.strategies import SearchStrategy


import wic
import wic.cli
import wic.main
import wic.schemas
import wic.schemas.wic_schema
import wic.utils as utils
from wic.wic_types import (RoseTree, StepId, Yaml, YamlForest, YamlTree, Json)
from .test_setup import wic_strategy, get_args


def check_nested_dict_for_key(test_dict: Dict,key_item: str) -> bool:
    """Check to see if substring is in any nested dictionary keys

    Args:
        test_dict (Dict): Nested dictionary

    Returns:
        bool: True if found subkey
    """
    cont=True
    dict_list=grab_nested_dicts(test_dict)
    
    for key,value in test_dict.items():
        if key_item in key: # using this as delimiter for flattening dictionary
            cont=False
    for otest_dict in dict_list:
        for key,value in otest_dict.items():
            if key_item in key:
                cont=False
    return cont


def add_wic_steps(wic_dict: Dict,step_tuple_list: List[Tuple[int,str]]) -> Dict:
    """Add steps into nested dictionary to simulate run time dictionary if wic_strategy doesnt have in generated dictionary
        need this because we chose not to include the property wic in our schema

    Args:
        wic_dict (Dict): Nested dictionary

    Returns:
        Dict: Modified dictionary
    """
    if 'wic' not in wic_dict.keys():
        wic_dict['wic']={}
    inner_wic_dict=wic_dict['wic']
    if 'steps' not in inner_wic_dict.keys():
        inner_wic_dict['steps']={}
    wic_steps = inner_wic_dict['steps']
    if len(wic_steps.keys())==0:
        for item in step_tuple_list:
            wic_steps[str(item)] = ""
    inner_wic_dict['steps']=wic_steps
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

def flatten_list(S: List) -> List:
    """Flatten a nested list

    Args:
        S (List): Nested list

    Returns:
        List: Flattened list
    """
    if S == []:
        return S
    if isinstance(S[0], list):
        return flatten_list(S[0]) + flatten_list(S[1:])
    return S[:1] + flatten_list(S[1:])

def grab_nested_dicts(d : Dict) -> Generator:
    
    """Recursively grab nested dicts from a dict

    Args:
        d (Dict): Recursive nested dictionary

    Yields:
        Generator[Dict]: All the sub nested dictionaries in generator object
    """
    for value in d.values():
        if isinstance(value, Dict):
            yield(value)
            yield from grab_nested_dicts(value)


def unflatten_dict(test_dict: Dict) -> Dict:
    """Unflatten nested dictionary that was flattend with flatten_nested_dict

    Args:
        test_dict (_type_): Flattened dictionary

    Returns:
        _type_: Unflattened dictionary
    """
    resultDict : Dict = dict()
    for key, value in test_dict.items():
        parts = list(key)
        d = resultDict
        for part in parts[:-1]:
            if part not in d:
                d[part] = dict()
            d = d[part]
        d[parts[-1]] = value
        
    return resultDict


def flatten_nested_dict(d : Dict) -> Dict:
    """Flatten a dict

    Args:
        d (Dict): Nested dictionary

    Returns:
        Dict: All items in dictionary are flattened into 1 dimension of keys
    """
    def items()-> Generator:
        for key, value in d.items():
            if isinstance(value, dict):
                for subkey, subvalue in flatten_nested_dict(value).items():
                    yield key +'_'+subkey, subvalue
                if len(value) == 0:
                    yield key, value
            else:
                yield key, value

    return dict(items())

def convert_dict_keys(d : Dict) -> Dict:
    """Convert dict keys from string to tuples

    Args:
        d (Dict): Used in conjunction with flatten_nested_dict, this converts keys that are stored as key1+_key2+_ .. into a tuple (key1,key2,..)

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
        test_yaml_tree (YamlTree): YamlTree (originally made from nested dictionary) that will be used to create YamlForest

    Returns:
        YamlForest: Generated from YamlTree 
    """
    subfors=list(recursive_forest_items(test_yaml_tree))
    if len(subfors)!=0:
        yaml_forest=YamlForest(yaml_tree=test_yaml_tree,sub_forests=[subfors]) 
    else:
        yaml_forest=YamlForest(yaml_tree=test_yaml_tree,sub_forests=[tuple([StepId('', ''),\
            YamlForest(yaml_tree=create_yaml_tree({}),sub_forests=[])])])
    return yaml_forest


filtered_pretty_text: SearchStrategy[str] = text(printable,min_size=1).filter(lambda x: '_' not in x)
any_types: Any = none() | booleans() | integers() | floats() | filtered_pretty_text
nested_lists: Any = s.recursive(s.lists(any_types,min_size=1),\
     lambda children: s.lists(children,min_size=1))
filtered_text: SearchStrategy[str] = s.text(min_size=1).filter(lambda x: '_' not in x)
recursive_fn: Any = lambda children: s.dictionaries(text(min_size=1),children,min_size=1)
nested_dict: Any = s.recursive(s.dictionaries(filtered_text,any_types).filter(lambda x: len(x)!=0),recursive_fn )

class TestUnits(unittest.TestCase):
    @pytest.mark.fast
    def test_read_lines_pairs(self)-> None: # generalize to making files on the fly
        """Test read_lines_pairs
        """
        testpath=os.path.abspath(os.path.join(__file__,os.pardir))
        wicpath=os.path.abspath(os.path.join(testpath,os.pardir))
        srcpath=os.path.join(wicpath,'src')
        wicsrcpath=os.path.join(srcpath,'wic')
        cwl_dirs_path=Path(os.path.join(wicsrcpath,'cwl_dirs.txt'))
        self.assertEqual([('global', 'biobb/'), ('global', 'cwl_adapters/')],utils.read_lines_pairs(cwl_dirs_path))

    @pytest.mark.fast
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

    @pytest.mark.fast 
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

    @pytest.mark.fast
    @given(stringlist=s.lists(s.text(min_size=2),min_size=1,max_size=100),astringlist=s.lists(s.text(min_size=2),min_size=1,max_size=50))
    def test_partition_by_lowest_common_ancestor(self,stringlist : List[str], astringlist : List[str])->None: 
        """Test partition_by_lowest_common_ancestor

        Args:
            stringlist (List[str]): List of random strings
            astringlist (List[str]): List of random strings
        """
        testlist=astringlist.copy()
        for item in astringlist:
            templist=[item,None]
            citem : Any =random.choice(templist)
            if citem is not None:
                testlist.insert(random.randrange(len(testlist)+1), citem) # randomly add different things to dupicate list
                
 
        results=utils.partition_by_lowest_common_ancestor(stringlist,testlist)
        found=False
        for i,test_string in enumerate(stringlist):
            same=True
            teststring=testlist[i]
            if i<=len(testlist)-1:
                if test_string!=teststring:
                    same=False
                if i==0 and same is False:
                    expected_results : Tuple[List,List] =([],stringlist)
                    found=True
                    break
                else:
                    if same is False:
                        lastsameindex=i-1
                        expected_results=(stringlist[:lastsameindex+1],stringlist[lastsameindex+1:])
                        found=True
                        break
        if found is False:
            expected_results=(stringlist,[])
        self.assertEqual(expected_results,results)

    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much,HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
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

    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much,HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
    )  
    @given(yaml_tree_tuple=wic_strategy.map(create_yaml_tree),stringlist=s.lists(s.text(min_size=2),min_size=1,max_size=100),astringlist=s.lists(s.text(min_size=2),min_size=1,max_size=50))
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
                    wic_steps = yaml_tree['wic'].get('steps', {})
                    steps_keys = utils.get_steps_keys(steps)
                    tools = wic.main.get_tools_cwl(args.cwl_dirs_file)
                    tools_stems = [stepid.stem for stepid in tools]
                    
            else:
                the_ls=randomly_add_contents_lists_to_list([stringlist,astringlist],0) # account for case where wic_strategy doesnt give correct data type
                steps_keys,tools_stems=the_ls[:]

                    
            results=utils.get_subkeys(steps_keys,tools_stems)
            expected_results=[i for i in steps_keys if i not in tools_stems]
            self.assertEqual(expected_results,results)

    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much,HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
    )  
    @given(
        yaml_tree=wic_strategy,
        wic_dict=wic_strategy,
        yaml_path=s.text(),
        backend=s.text(min_size=1),
        steps=s.integers(min_value=0)
    )
    def test_extract_backend(self, yaml_tree : Dict, wic_dict : Dict, yaml_path : str ,backend : str,steps : int)->None:
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
            wic_dict=add_backend_steps(wic_dict,backend,steps) # in case wic_strategy doesnt give correct data structure
            yaml_tree_copy.update({'steps': steps})
        results=utils.extract_backend(yaml_tree,wic_dict,Path(yaml_path))
        self.assertEqual((backend,yaml_tree_copy),results)
    
    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much, HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
    )  
    @given(masterlist=nested_lists)
    def test_flatten(self, masterlist : List[List])->None:
        """Test flatten

        Args:
            masterlist (List[List]): Nested lists of random data types
        """
        if any(isinstance(el, list) for el in masterlist): # recursive funciton sometimes gives list with no nested lists
            expected_output=[x for lst in masterlist for x in lst]
            self.assertEqual(expected_output,utils.flatten(masterlist))

    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much, HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
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

    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much, HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
    )
    @given(test_dict = nested_dict)
    def test_recursively_delete_dict_key(self, test_dict : Dict) -> None:
        """Test recursively_delete_dict_key

        Args:
            test_dict (Dict): Random nested dictionary
        """
        if len(test_dict)!=0:
            cont=check_nested_dict_for_key(test_dict,'_')
            flatten_dict=convert_dict_keys(flatten_nested_dict(test_dict)) # convert to format to call unflatten function later
            if cont is True:
                fullkeylist=list(flatten_dict.keys())
                if len(fullkeylist)>1:
                    keylist=random.choice(fullkeylist)
                elif len(fullkeylist)==1:
                    keylist=fullkeylist[0]
                else:
                    cont=False
                if cont is True:
                    flatten_keylist=[item for sublist in keylist for item in sublist]
                    if len(flatten_keylist)!=0:
                        key=random.choice(flatten_keylist) # choose random key to delete
                        results = utils.recursively_delete_dict_key(key, test_dict)
                        keylists_to_delete=[]
                        for key_list in fullkeylist:
                            if key in key_list:
                                keylists_to_delete.append(key_list)
                        
                        new_keys=[]
                        for key_list in keylists_to_delete:
                            del flatten_dict[key_list] 
                            if len(key_list)>1 and key!=key_list[0]:
                                new_key_list=[]
                                for i in key_list:
                                    if i!=key:
                                        new_key_list.append(i)
                                    else:
                                        break
                                
                                new_keys.append(tuple(new_key_list))
                                flatten_dict[tuple(new_key_list)]={}
                        dups_to_remove=[]
                        for key_list in flatten_dict.keys():
                            if key_list in new_keys:
                                for okey_list in flatten_dict.keys():
                                    if okey_list!=key_list and len(key_list)<len(okey_list):
                                        is_dup=True
                                        for i,ovalue in enumerate(okey_list):
                                            if i<=len(key_list)-1:
                                                value=key_list[i]
                                                if value!=ovalue:
                                                    is_dup=False
                                        if is_dup is True:
                                            if key_list not in dups_to_remove:
                                                dups_to_remove.append(key_list)            
                                
                        for key_list in dups_to_remove:
                            del flatten_dict[key_list] 
                        
                        expected_results=unflatten_dict(flatten_dict)
                        self.assertEqual(expected_results,results)


    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much, HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
    )
    @given(test_dict = nested_dict)
    def test_recursively_contains_dict_key(self,test_dict: Dict) -> None: # a bit trivial, flatten keys first
        dict_list=grab_nested_dicts(test_dict)
        keylist=[list(i.keys()) for i in dict_list]
        keylist=list(flatten_list(keylist))
        key=random.choice(keylist)
        results=utils.recursively_contains_dict_key(key,test_dict)
        self.assertEqual(True,results)
        

        
    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much, HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
    )
    @given(
        wic_dict=wic_strategy,
        step_tuple_list=s.lists(s.tuples(s.integers(min_value=0),s.text(min_size=1).filter(lambda x: ',' not in x)),min_size=1,max_size=10),
    )
    def test_reindex_wic_steps(self,wic_dict: Dict,step_tuple_list: List[Tuple[int,str]]) -> None: # a bit trivial, no alternative solution, replace with a weaker test?
        """Test reindex_wic_steps

        Args:
            wic_dict (Dict): Dictionary with wic_strategy
        """
        has_steps=check_if_steps_in_dict(wic_dict)
        if has_steps is False:
            wic_dict = add_wic_steps(wic_dict,step_tuple_list)
        wic_steps = wic_dict['wic'].get('steps', {})
        random_step=random.choice(list(wic_steps.keys()))
        random_step=utils.parse_int_string_tuple(random_step)
        random_index=random_step[0]
        results=utils.reindex_wic_steps(wic_steps,random_index)
        new_wic_steps={}
        for key,value in wic_steps.items():
            tupkey=list(utils.parse_int_string_tuple(key))
            index = tupkey[0]
            newstr = f'({index+1}, {tupkey[1]})' if index >= random_index else key
            new_wic_steps[newstr]=value
        expected_results=new_wic_steps
        self.assertEqual(expected_results,results)


    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much, HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
    )
    @given(test_dict = nested_dict, key=filtered_pretty_text,val=any_types)
    def test_recursively_insert_into_dict_tree(self,test_dict: Dict, key: str, val: Any) -> None:
        """_summary_

        Args:
            test_dict (Dict): Nested dictionary
            key (str): Random string
            val (Any): Random type
        """
        
        dict_list=list(grab_nested_dicts(test_dict))
        cont=check_nested_dict_for_key(test_dict,'_')
        if len(dict_list)>1 and cont is True:
            flatten_dict=convert_dict_keys(flatten_nested_dict(test_dict))
            fullkeylist=list(flatten_dict.keys())
            key_list=random.choice(fullkeylist)
            if type(flatten_dict[key_list]) is dict:
                del flatten_dict[key_list]
                key_list=list(key_list)
                key_list.append(key)
                results=utils.recursively_insert_into_dict_tree(test_dict,key_list,val)
                key_list=tuple(key_list)
                flatten_dict[key_list]=[val]
            else: 
                results=test_dict
            expected_results=unflatten_dict(flatten_dict)
            self.assertEqual(expected_results,results)

    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much, HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
    )
    @given(input_mapping = s.dictionaries(text(min_size=1),s.lists(text(min_size=1),min_size=1),min_size=1).filter(lambda x: len(x)!=0))
    def test_get_input_mappings(self,input_mapping: Dict[str,List[str]]) -> None:
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



    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much, HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
    )
    @given(output_mapping = s.dictionaries(text(min_size=1),text(min_size=1),min_size=1))
    def test_get_output_mappings(self,output_mapping: Dict[str,str]) -> None:
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
            



    



        
        
                    
                    


   
 

    


