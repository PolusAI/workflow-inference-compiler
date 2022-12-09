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
  for v in d.values():
    if isinstance(v, Dict):
        yield(v)
        yield from grab_nested_dicts(v)


def add_backend_steps(wic : Dict, backend : str, steps : int) -> Dict:
    wic['backend']=backend
    wic['backends']={}
    plugin_ns = wic.get('namespace', 'global')
    stepid = StepId(backend, plugin_ns)
    wic['backends'][stepid]={}
    wic['backends'][stepid]['steps']=steps
    return wic


def randomly_add_contents_lists_to_list(ls : List , index : int) -> List:
    newls=[]
    contentlist=ls[index]
    for i in range(len(ls)):
        if i !=index:
            templs=ls[i].copy()
            item=random.choice(templs)
            contentlist.append(item)
            newls.append(templs)
        else:
            newls.append(contentlist)

    return newls


def get_args(yml_path: str = '') -> argparse.Namespace:
    testargs = ['wic', '--yaml', yml_path, '--cwl_output_intermediate_files', 'True']  # ignore --yaml
    with patch.object(sys, 'argv', testargs):
        args: argparse.Namespace = wic.cli.parser.parse_args()
    return args

def grab_wic_schemas() -> Tuple[Any,Dict]:
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
    wic_main_schema = wic.schemas.wic_schema.wic_main_schema(tools_cwl,yaml_stems,schema_store,hypothesis=True)
    wic_strategy = from_schema(wic_main_schema)
    
    return wic_strategy,wic_main_schema


def recursive_tree_items(dictionary : Dict) -> Generator:
    for key, value in dictionary.items():
        if type(value) is dict:
            yield RoseTree(data=value,sub_trees=[])
            yield from recursive_tree_items(value)

def recursive_forest_items(yaml_tree : YamlTree) -> Generator:
    storedDicts=[]
    for key, value in yaml_tree.yml.items():
        if type(value) is dict:
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
            if testdict not in storedDicts and len(testdict)!=0:
                yield tup
                storedDicts.append(testdict) # dont want duplicates 
            yield from recursive_forest_items(yml_tree)

def create_rose_tree_from_nested_dict(testdict : Dict) ->RoseTree:
    subdicts=recursive_tree_items(testdict)
    rosetree=RoseTree(data=testdict,sub_trees=list(subdicts)) 
    return rosetree

def create_yaml_tree(dict : Dict) ->YamlTree:
    backend=''
    if 'default_backend' in dict:
        backend = dict['default_backend']
    if 'backend' in dict:
        backend = dict['backend']
    plugin_ns = dict.get('namespace', 'global')
    stepid = StepId(backend, plugin_ns)
    yaml_tree=YamlTree(step_id=stepid,yml=dict)
    return yaml_tree

def create_yaml_forest(test_yaml_tree : YamlTree) -> YamlForest:
    subfors=list(recursive_forest_items(test_yaml_tree))
    if len(subfors)!=0:
        yaml_forest=YamlForest(yaml_tree=test_yaml_tree,sub_forests=[subfors]) 
    else:
        yaml_forest=YamlForest(yaml_tree=test_yaml_tree,sub_forests=[tuple([StepId('', ''),YamlForest(yaml_tree=create_yaml_tree({}),sub_forests=[])])]) 
    return yaml_forest

nested_lists : Any = s.recursive(s.lists(none() | booleans() | integers() | floats() | text(printable),min_size=1), lambda children: s.lists(children,min_size=1))
nested_dict : Any = s.recursive(s.dictionaries(s.text(),none() | booleans() | integers() | floats() | text(printable)), lambda children: s.dictionaries(text(),children,min_size=1))
wic_strategy,wic_main_schema=grab_wic_schemas()

class TestUnits(unittest.TestCase):
    @pytest.mark.fast
    def test_read_lines_pairs(self)-> None: # generalize to making files on the fly
        testpath=os.path.abspath(os.path.join(__file__,os.pardir))
        wicpath=os.path.abspath(os.path.join(testpath,os.pardir))
        srcpath=os.path.join(wicpath,'src')
        wicsrcpath=os.path.join(srcpath,'wic')
        cwl_dirs_path=Path(os.path.join(wicsrcpath,'cwl_dirs.txt'))
        self.assertEqual([('global', 'biobb/'), ('global', 'cwl_adapters/')],utils.read_lines_pairs(cwl_dirs_path))

    @pytest.mark.fast
    @given(item1=s.text().filter(lambda x: '_' not in x),item2=s.integers(min_value=0,max_value=10 ** 6),item3=s.text().filter(lambda x: '_' not in x))
    def test_parse_step_name_str(self,item1 : str, item2 : int, item3 : str)-> None:
        output=utils.step_name_str(item1, item2 ,item3)
        self.assertEqual((item1, item2, item3),utils.parse_step_name_str(output))

    @pytest.mark.fast 
    @given(item1=s.text(min_size=1).filter(lambda x: '_' not in x),item2=s.text(min_size=1).filter(lambda x: '_' not in x))
    def test_shorten_restore_namespaced_output_name(self,item1 : str,item2 : str)->None:
        thesep='__'
        namespaced_output_name=item1+thesep+item2
        tup=utils.shorten_namespaced_output_name(namespaced_output_name,sep=thesep)
        restored_name=utils.restore_namespaced_output_name(tup[0],tup[1],sep=thesep)
        self.assertEqual(restored_name,namespaced_output_name)

    @pytest.mark.fast
    @given(stringlist=s.lists(s.text(min_size=2),min_size=1,max_size=100),astringlist=s.lists(s.text(min_size=2),min_size=1,max_size=50))
    def test_partition_by_lowest_common_ancestor(self,stringlist : List[str], astringlist : List[str])->None: 
        testlist=astringlist.copy()
        for item in astringlist:
            templist=[item,None]
            citem : Any =random.choice(templist)
            if citem!=None:
                testlist.insert(random.randrange(len(testlist)+1), citem) # randomly add different things to dupicate list
                
 
        results=utils.partition_by_lowest_common_ancestor(stringlist,testlist)
        found=False
        for i in range(len(stringlist)):
            same=True
            string=stringlist[i]
            teststring=testlist[i]
            if i<=len(testlist)-1:
                if string!=teststring:
                    same=False
                if i==0 and same==False:
                    expected_results : Tuple[List,List] =([],stringlist)
                    found=True
                    break
                else:
                    if same==False:
                        lastsameindex=i-1
                        expected_results=(stringlist[:lastsameindex+1],stringlist[lastsameindex+1:])
                        found=True
                        break
        if found==False:
            expected_results=(stringlist,[])
        self.assertEqual(expected_results,results)

    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much,HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
    )  
    @given(dictlist=s.lists(nested_dict,min_size=1,max_size=10))
    def test_get_steps_keys(self,dictlist : List[Dict])->None: # trivial test case?
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
            if found==False:
                ls=randomly_add_contents_lists_to_list([stringlist,astringlist],0) # account for case where wic_strategy doesnt give correct data type
                steps_keys,tools_stems=ls[:]

                    
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
        wic=wic_strategy,
        yaml_path=s.text(),
        backend=s.text(min_size=1),
        steps=s.integers(min_value=0)
    )
    def test_extract_backend(self, yaml_tree : Dict, wic : Dict, yaml_path : str ,backend : str,steps : int)->None:
        yaml_tree_copy=copy.deepcopy(yaml_tree)
        if 'backend' in wic.keys():
            backend=wic['backend']
        else:
            wic=add_backend_steps(wic,backend,steps) # in case wic_strategy doesnt give correct data structure
            yaml_tree_copy.update({'steps': steps})
        results=utils.extract_backend(yaml_tree,wic,Path(yaml_path))
        self.assertEqual((backend,yaml_tree_copy),results)
    
    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much, HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
    )  
    @given(masterlist=nested_lists)
    def test_flatten(self, masterlist : List[List])->None: # trivial unit test? doesnt handle nested lists only
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
        def flatten_rose_tree(rosetree : RoseTree) -> List:
            data=rosetree.data
            flattened_dict=list(grab_nested_dicts(data))
            return [data]+[i for i in flattened_dict]
            
        expected_results=flatten_rose_tree(rosetree)
        results=utils.flatten_rose_tree(rosetree)
        self.assertEqual(expected_results,results)

    @pytest.mark.fast
    @settings(
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much, HealthCheck.data_too_large],
    deadline=datetime.timedelta(milliseconds=1000),
    )  
    @given(forest=nested_dict.map(create_yaml_tree).map(create_yaml_forest))
    def test_flatten_forest(self,forest : YamlForest) -> None: # FIX ME, 
        yaml_tree=forest.yaml_tree
        dict=yaml_tree.yml
        if len(dict)==0: # trivial case, seems to have issues too i.e if testdict=={}
            pass
        else:
            def flatten_forest(forest : YamlForest) -> List:
                yaml_tree=forest.yaml_tree
                testdict=yaml_tree.yml
                flattened_dict=list(grab_nested_dicts(testdict))
                forests=[create_yaml_forest(create_yaml_tree(d)) for d in flattened_dict]
                return forests
            results=utils.flatten_forest(forest)
            expected_results=flatten_forest(forest)
            self.assertEqual(expected_results,results)

 

    


