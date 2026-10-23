"""Behavioural coverage for `sophios.input_output` config handling, and for the paths a written job names."""
import json
from pathlib import Path, PureWindowsPath

import pytest

from sophios import input_output as io
from sophios.ir.locations import absolute_paths


@pytest.mark.fast
@pytest.mark.parametrize('tag', ['search_paths_cwl', 'search_paths_wic'])
def test_a_search_path_given_as_a_string_is_a_type_error(tmp_path: Path, tag: str) -> None:
    """A string is not iterated into characters; the error names the config key."""
    config = tmp_path / 'c.json'
    config.write_text(json.dumps({'search_paths_cwl': {'global': []}, 'search_paths_wic': {'global': []},
                                  tag: {'global': '/whole/disk'}}), encoding='utf-8')
    with pytest.raises(TypeError, match=f'{tag}.global'):
        io.read_config_from_disk(config)


@pytest.mark.fast
def test_a_relative_job_location_is_written_as_a_file_uri() -> None:
    """A `location` is a URI wherever the job is written, a Windows drive too; a `path` stays a path."""
    values = {'f': {'class': 'File', 'location': 'data/a.txt',
                    'secondaryFiles': [{'class': 'File', 'path': 'data/a.idx'}]}}
    written = absolute_paths(values, PureWindowsPath('C:/project'))
    assert written['f']['location'] == 'file:///C:/project/data/a.txt'
    assert written['f']['secondaryFiles'][0]['path'] == str(PureWindowsPath('C:/project/data/a.idx'))
    assert values['f']['location'] == 'data/a.txt'
