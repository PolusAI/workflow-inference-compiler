# See https://docs.pytest.org/en/7.1.x/example/simple.html
import pytest


def pytest_addoption(parser) -> None:
    parser.addoption('--cwl_runner', type=str, required=False, default='cwltool', choices=['cwltool', 'toil-cwl-runner'],
                     help='The CWL runner to use for running workflows locally.')


@pytest.fixture
def cwl_runner(request):
    return request.config.getoption("--cwl_runner")
