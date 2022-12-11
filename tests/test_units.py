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
from typing import Dict,List,Any,Tuple, Generator
import datetime
from pathlib import Path
from flatten_dict import flatten,unflatten

from hypothesis import given, settings,HealthCheck
from hypothesis.strategies import text, integers, booleans, emails, floats, dictionaries, recursive, lists, \
    composite, datetimes, none, fixed_dictionaries
import hypothesis.strategies as s
from hypothesis import given
from hypothesis_jsonschema import from_schema
from hypothesis.strategies._internal.strategies import SearchStrategy


import wic
import wic.cli
import wic.main
import wic.schemas
import wic.schemas.wic_schema
import wic.utils as utils
from wic.wic_types import (ExplicitEdgeCalls, Namespaces, NodeData, RoseTree, StepId, Yaml, YamlForest, YamlTree, Json)



def grab_nested_dicts(d : Dict) -> Generator:
    """Recursively grab nested dicts from a dict"""
    for value in d.values():
        if isinstance(value, Dict):
            yield(value)
            yield from grab_nested_dicts(value)



def flatten_nested_dict(d : Dict) -> Dict:
    """Flatten a dict"""
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
    """Convert dict keys from string to tuples"""
    new_dict = {}
    for key, value in d.items():
        if isinstance(key, str):
            key = tuple(key.split('_'))
  
        new_dict[tuple(key)] = value
    return new_dict
    
    


def add_backend_steps(wic_dict : Dict, backend : str, steps : int) -> Dict:
    """Add backend and steps to a wic dict"""
    wic_dict['backend']=backend
    wic_dict['backends']={}
    plugin_ns = wic_dict.get('namespace', 'global')
    stepid = StepId(backend, plugin_ns)
    wic_dict['backends'][stepid]={}
    wic_dict['backends'][stepid]['steps']=steps
    return wic_dict


def randomly_add_contents_lists_to_list(test : List , index : int) -> List:
    """Randomly add contents of a list to another list"""
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


def get_args(yml_path: str = '') -> argparse.Namespace:
    """Get args for wic"""
    testargs = ['wic', '--yaml', yml_path, '--cwl_output_intermediate_files', 'True']  # ignore --yaml
    with patch.object(sys, 'argv', testargs):
        args: argparse.Namespace = wic.cli.parser.parse_args()
    return args

def grab_wic_schemas() -> Tuple[Any,Dict]:
    """Grab wic schemas"""
    args = get_args()
    tools_cwl = wic.main.get_tools_cwl(args.cwl_dirs_file)
    yml_paths = wic.main.get_yml_paths(args.yml_dirs_file)
    yaml_stems = wic.utils.flatten([list(p) for p in yml_paths.values()])
    schema_store : Dict[str, Json]= {} 
    validator = wic.schemas.wic_schema.get_validator(tools_cwl, yaml_stems, schema_store, write_to_disk=True,hypothesis=True)
    if args.generate_schemas_only:
        yml_paths_tuples = [(yml_path_str, yml_path)
            for yml_namespace, yml_paths_Dict in yml_paths.items()
            for yml_path_str, yml_path in yml_paths_Dict.items()]
        for yml_path_str, yml_path in yml_paths_tuples:
            schema = wic.schemas.wic_schema.compile_workflow_generate_schema(yml_path_str, yml_path,tools_cwl, yml_paths, validator)
            # overwrite placeholders in schema_store. See comment in get_validator()
            schema_store[schema['$id']] = schema
    wic_main_schema_ins = wic.schemas.wic_schema.wic_main_schema(tools_cwl,yaml_stems,schema_store,hypothesis=True)
    wic_strategy_ins = from_schema(wic_main_schema_ins)
    
    return wic_strategy_ins,wic_main_schema_ins


def recursive_tree_items(dictionary : Dict) -> Generator:
    """Recursively grab nested dicts from a dict and generate RoseTrees"""
    for key, value in dictionary.items():
        if isinstance(value,dict):
            yield RoseTree(data=value,sub_trees=[])
            yield from recursive_tree_items(value)

def recursive_forest_items(yaml_tree : YamlTree) -> Generator:
    """Recursively grab nested dicts from a dict and generate YamlForests"""
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
            tup=tuple([stepid,forest])
            if testdict not in stored_dicts and len(testdict)!=0:
                yield tup
                stored_dicts.append(testdict) # dont want duplicates 
            yield from recursive_forest_items(yml_tree)

def create_rose_tree_from_nested_dict(testdict : Dict) ->RoseTree:
    """Create a RoseTree from a nested dict"""
    subdicts=recursive_tree_items(testdict)
    rosetree=RoseTree(data=testdict,sub_trees=list(subdicts))
    return rosetree

def create_yaml_tree(test_dict : Dict) ->YamlTree:
    """Create a YamlTree from a nested dict"""
    backend=''
    if 'default_backend' in test_dict:
        backend = test_dict['default_backend']
    if 'backend' in test_dict:
        backend = test_dict['backend']
    plugin_ns = test_dict.get('namespace', 'global')
    stepid = StepId(backend, plugin_ns)
    yaml_tree=YamlTree(step_id=stepid,yml=test_dict)
    return yaml_tree

def create_yaml_forest(test_yaml_tree : YamlTree) -> YamlForest:
    """Create a YamlForest from a nested dict"""
    subfors=list(recursive_forest_items(test_yaml_tree))
    if len(subfors)!=0:
        yaml_forest=YamlForest(yaml_tree=test_yaml_tree,sub_forests=[subfors]) 
    else:
        yaml_forest=YamlForest(yaml_tree=test_yaml_tree,sub_forests=[tuple([StepId('', ''),\
            YamlForest(yaml_tree=create_yaml_tree({}),sub_forests=[])])])
    return yaml_forest

nested_lists : Any = s.recursive(s.lists(none() | booleans() | integers() | floats() | text(printable),min_size=1),\
     lambda children: s.lists(children,min_size=1))
nested_dict : Any = s.recursive(s.dictionaries(s.text(min_size=1).filter(lambda x: '_' not in x),none() | booleans() | integers() | floats() | text(printable).filter(lambda x: '_' not in x)), \
    lambda children: s.dictionaries(text(min_size=1),children,min_size=1))
wic_strategy,wic_main_schema=grab_wic_schemas()

class TestUnits(unittest.TestCase):
    @pytest.mark.fast
    def test_read_lines_pairs(self)-> None: # generalize to making files on the fly
        """Test read_lines_pairs"""
        testpath=os.path.abspath(os.path.join(__file__,os.pardir))
        wicpath=os.path.abspath(os.path.join(testpath,os.pardir))
        srcpath=os.path.join(wicpath,'src')
        wicsrcpath=os.path.join(srcpath,'wic')
        cwl_dirs_path=Path(os.path.join(wicsrcpath,'cwl_dirs.txt'))
        self.assertEqual([('global', 'biobb/'), ('global', 'cwl_adapters/')],utils.read_lines_pairs(cwl_dirs_path))

    @pytest.mark.fast
    @given(item1=s.text().filter(lambda x: '_' not in x),item2=s.integers(min_value=0,max_value=10 ** 6),item3=s.text().filter(lambda x: '_' not in x))
    def test_parse_step_name_str(self,item1 : str, item2 : int, item3 : str)-> None:
        """Test parse_step_name_str"""
        output=utils.step_name_str(item1, item2 ,item3)
        self.assertEqual((item1, item2, item3),utils.parse_step_name_str(output))

    @pytest.mark.fast 
    @given(item1=s.text(min_size=1).filter(lambda x: '_' not in x),item2=s.text(min_size=1).filter(lambda x: '_' not in x))
    def test_shorten_restore_namespaced_output_name(self,item1 : str,item2 : str)->None:
        """Test shorten_restore_namespaced_output_name"""
        thesep='__'
        namespaced_output_name=item1+thesep+item2
        tup=utils.shorten_namespaced_output_name(namespaced_output_name,sep=thesep)
        restored_name=utils.restore_namespaced_output_name(tup[0],tup[1],sep=thesep)
        self.assertEqual(restored_name,namespaced_output_name)

    @pytest.mark.fast
    @given(stringlist=s.lists(s.text(min_size=2),min_size=1,max_size=100),astringlist=s.lists(s.text(min_size=2),min_size=1,max_size=50))
    def test_partition_by_lowest_common_ancestor(self,stringlist : List[str], astringlist : List[str])->None: 
        """Test partition_by_lowest_common_ancestor"""
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
    @given(dictlist=s.lists(nested_dict,min_size=1,max_size=10))
    def test_get_steps_keys(self,dictlist : List[Dict])->None: 
        """Test get_steps_keys"""
        firstkeylist=[]
        for unflattendict in dictlist:
            firstkeylist.extend(list(unflattendict.keys()))
        self.assertEqual(firstkeylist,utils.get_steps_keys(dictlist))

    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much,HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
    )  
    @given(yaml_tree_tuple=wic_strategy.map(create_yaml_tree),stringlist=s.lists(s.text(min_size=2),min_size=1,max_size=100),astringlist=s.lists(s.text(min_size=2),min_size=1,max_size=50))
    def test_get_subkeys(self,yaml_tree_tuple : YamlTree,stringlist : List[str], astringlist : List[str])->None:
            """Test get_subkeys"""
            args = get_args()
            found=False
            (step_id, yaml_tree) = yaml_tree_tuple
            if 'steps' in yaml_tree.keys():
                steps: List[Yaml] = yaml_tree['steps']
                if 'wic' in yaml_tree.keys():
                    wic_steps = yaml_tree['wic'].get('steps', {})
                    steps_keys = utils.get_steps_keys(steps)
                    tools = wic.main.get_tools_cwl(args.cwl_dirs_file)
                    tools_stems = [stepid.stem for stepid in tools]
                    found=True
            if found is False:
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
        """Test extract_backend"""
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
        """Test flatten"""
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
        """Test flatten_rose_tree"""
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
        """Test recursively_delete_dict_key"""
        if len(test_dict)!=0:
            cont=True
            dict_list=grab_nested_dicts(test_dict)
            
            for key,value in test_dict.items():
                if '_' in key: # using this as delimiter for flattening dictionary
                    cont=False
            for otest_dict in dict_list:
                for key,value in otest_dict.items():
                    if '_' in key:
                        cont=False
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
                        key=random.choice(flatten_keylist)
                        results = utils.recursively_delete_dict_key(key, test_dict)
                        keylists_to_delete=[]
                        for key_list in fullkeylist:
                            if key in key_list:
                                keylists_to_delete.append(key_list)
                        for key_list in keylists_to_delete:
                            del flatten_dict[key_list] 
                            if len(key_list)>1 and key!=key_list[0]:
                                new_key_list=[]
                                for i in key_list:
                                    if i!=key:
                                        new_key_list.append(i)
                                    else:
                                        break
                                
                                isduplicate=False
                                for okey_list in fullkeylist:
                                    if okey_list[0]==new_key_list[0] and okey_list!=key_list:
                                        isduplicate=True
                                if isduplicate is False:
                                    flatten_dict[tuple(new_key_list)]={}
                       
                        expected_results=unflatten(flatten_dict)
                        self.assertEqual(expected_results,results)
                    
                    


   
 

    


