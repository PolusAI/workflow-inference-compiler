# pylint: disable=W0621, W0212
from pathlib import Path
from typing import Any

import cwl_utils.parser as cu_parser
import hypothesis
import pytest
import yaml

from wic.api import Step, Workflow
from wic.api._types import CWL_TYPES_DICT
from wic.api.pythonapi import (
    _WIC_PATH,
    CLTInput,
    InvalidStepError,
    MissingRequiredValueError,
)

RSRC_DIR = Path(__file__).parent / "resources"


@pytest.fixture
def assembler_step() -> Step:
    """Assembler Step."""
    assembler_ = RSRC_DIR / "assembler.cwl"
    return Step(assembler_)


@pytest.fixture
def renaming_step() -> Step:
    """Renaming Step."""
    renaming_ = RSRC_DIR / "filerenaming.cwl"
    return Step(renaming_)


@pytest.fixture
def basic_step() -> Step:
    """Basic Step."""
    basic_ = RSRC_DIR / "basicfl.cwl"
    return Step(basic_)


@pytest.fixture
def montage_step() -> Step:
    """Basic Step."""
    montage_ = RSRC_DIR / "montage.cwl"
    return Step(montage_)


@pytest.fixture
def omeconverter_step() -> Step:
    """Basic Step."""
    omeconverter_ = RSRC_DIR / "omeconverter.cwl"
    return Step(omeconverter_)


@pytest.fixture(scope="session")
def basic_step_session() -> Step:
    """Basic Step."""
    basic_ = RSRC_DIR / "basicfl.cwl"
    return Step(basic_)


@pytest.fixture
def steps(assembler_step: Step, renaming_step: Step, basic_step: Step) -> list[Step]:
    """Return list of steps."""
    return [assembler_step, renaming_step, basic_step]


@pytest.fixture
def assembler_cwl() -> Any:
    """Assembler CWL."""
    assembler_ = RSRC_DIR / "assembler.cwl"
    return cu_parser.load_document_by_uri(assembler_)


@pytest.fixture(scope="session")
def cwl_adapters(tmp_path_factory: Path) -> Path:
    """CWL Adapters Tmp Dir."""
    tmp = tmp_path_factory.mktemp("tmp")  # type: ignore
    return tmp  # type: ignore


def test_cltinp_type(assembler_cwl: Any) -> None:
    """Test CLTInput type."""
    for inp in assembler_cwl.inputs:
        clt_inp = CLTInput(inp)
        assert isinstance(clt_inp.inp_type, object)
        # set type_ variable used for comparison
        if hasattr(inp, "type_"):
            if isinstance(inp.type_, list):  # optional (null, boolean)
                type_ = inp.type_[1]
            else:
                type_ = inp.type_
        elif hasattr(inp, "type"):
            if isinstance(inp.type, list):  # optional (null, boolean)
                type_ = inp.type[1]
            else:
                type_ = inp.type
        assert clt_inp.inp_type == CWL_TYPES_DICT[type_]


def test_step1(assembler_step: Step) -> None:
    """Test Step 1."""
    st = assembler_step
    assert isinstance(st, Step)
    assert st.clt_name == "assembler"
    assert all(isinstance(i, CLTInput) for i in st.inputs)
    with RSRC_DIR.joinpath("assembler.cwl").open("r", encoding="utf-8") as file:
        assert st.yaml == yaml.safe_load(file)
    assert st.cwl_version == "v1.2"
    assert st.clt_path == RSRC_DIR / "assembler.cwl"


def test_step1_cwl(assembler_step: Step, assembler_cwl: Any) -> None:
    """Test Step 1."""
    st = assembler_step
    cwl_ = assembler_cwl
    assert st.clt == cwl_


def test_step1_str() -> None:
    """Test Step 1 with str."""
    st = Step(str(RSRC_DIR / "assembler.cwl"))
    assert isinstance(st, Step)
    assert st.clt_name == "assembler"
    assert all(isinstance(i, CLTInput) for i in st.inputs)
    with RSRC_DIR.joinpath("assembler.cwl").open("r", encoding="utf-8") as file:
        assert st.yaml == yaml.safe_load(file)
    assert st.cwl_version == "v1.2"
    assert st.clt_path == RSRC_DIR / "assembler.cwl"


def test_step2(renaming_step: Step) -> None:
    """Test Step 2."""
    st = renaming_step
    assert isinstance(st, Step)
    assert st.clt_name == "filerenaming"
    assert all(isinstance(i, CLTInput) for i in st.inputs)
    with RSRC_DIR.joinpath("filerenaming.cwl").open("r", encoding="utf-8") as file:
        assert st.yaml == yaml.safe_load(file)
    assert st.cwl_version == "v1.1"
    assert st.clt_path == RSRC_DIR / "filerenaming.cwl"


def test_step3(basic_step: Step) -> None:
    """Test Step 3."""
    st = basic_step
    assert isinstance(st, Step)
    assert st.clt_name == "basicfl"
    assert all(isinstance(i, CLTInput) for i in st.inputs)
    with RSRC_DIR.joinpath("basicfl.cwl").open("r", encoding="utf-8") as file:
        assert st.yaml == yaml.safe_load(file)
    assert st.cwl_version == "v1.0"
    assert st.clt_path == RSRC_DIR / "basicfl.cwl"


def test_cfg_yml1(renaming_step: Step) -> None:
    """Test IO cfg yaml 1."""
    st1 = renaming_step
    st2 = Step(RSRC_DIR / "filerenaming.cwl", RSRC_DIR / "filerenaming_io.yml")
    with RSRC_DIR.joinpath("filerenaming_io.yml").open("r", encoding="utf-8") as file:
        assert st2.cfg_yaml == yaml.safe_load(file)
    assert st1.clt_name == st2.clt_name
    assert all(isinstance(i, CLTInput) for i in st2.inputs)
    assert [x for x in st2.inputs if x.name ==
            "outFilePattern"][0].value == "x{row:dd}_y{col:dd}_p{s:dd}_c{channel:d}.ome.tif"
    assert [x for x in st2.inputs if x.name ==
            "filePattern"][0].value == ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
    assert [x for x in st2.inputs if x.name == "mapDirectory"][0].value == "raw"
    assert st2.cwl_version == "v1.1"


def test_cfg_yml1_str(renaming_step: Step) -> None:
    """Test IO cfg yaml 1 with str."""
    st1 = renaming_step
    st2 = Step(RSRC_DIR.joinpath("filerenaming.cwl"), str(RSRC_DIR / "filerenaming_io.yml"))
    with RSRC_DIR.joinpath("filerenaming_io.yml").open("r", encoding="utf-8") as file:
        assert st2.cfg_yaml == yaml.safe_load(file)
    assert st1.clt_name == st2.clt_name
    assert all(isinstance(i, CLTInput) for i in st2.inputs)
    assert [x for x in st2.inputs if x.name ==
            "outFilePattern"][0].value == "x{row:dd}_y{col:dd}_p{s:dd}_c{channel:d}.ome.tif"
    assert [x for x in st2.inputs if x.name ==
            "filePattern"][0].value == ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
    assert [x for x in st2.inputs if x.name == "mapDirectory"][0].value == "raw"
    assert st2.cwl_version == "v1.1"


def test_cfg_yml2() -> None:
    """Test IO cfg yaml 2."""
    st = Step(RSRC_DIR / "filerenaming.cwl", RSRC_DIR / "filerenaming_io.yml")
    assert st.clt_path == RSRC_DIR / "filerenaming.cwl"


def test_set_inp1(assembler_step: Step) -> None:
    """Test set input 1."""
    st = assembler_step
    st.imgPath = RSRC_DIR
    with pytest.raises(TypeError):
        st.outDir = str(RSRC_DIR / "out")
    with pytest.raises(TypeError, match="got str, expected Path"):
        st.outDir = str(RSRC_DIR / "out")
    st.outDir = RSRC_DIR / "out"
    st.preview = True
    st.timesliceNaming = False
    assert st.imgPath.value == RSRC_DIR  # type: ignore
    assert st.outDir.value == RSRC_DIR / "out"  # type: ignore


def test_set_inp2(assembler_step: Step) -> None:
    """Test set input 2."""
    st = assembler_step
    st.imgPath = RSRC_DIR
    st.outDir = RSRC_DIR / "out"
    st.preview = True
    st.timesliceNaming = False
    wic_yaml = {
        "assembler": {
            "in": {
                "imgPath": str(RSRC_DIR),
                "outDir": str(RSRC_DIR / "out"),
                "preview": True,
                "timesliceNaming": False
            }
        }
    }
    assert st._yml == wic_yaml


@hypothesis.given(
    filePattern=hypothesis.strategies.text(min_size=8),
    darkfield=hypothesis.strategies.booleans(),
    groupby=hypothesis.strategies.text(min_size=3)
)
def test_set_inp3(basic_step_session: Step,
                  filePattern: str,
                  darkfield: bool,
                  groupby: str) -> None:
    """Test set input 3."""
    st = basic_step_session
    st.filePattern = filePattern
    st.getDarkfield = darkfield
    st.groupBy = groupby
    assert st.filePattern.value == filePattern  # type: ignore
    assert st.getDarkfield.value == darkfield  # type: ignore
    assert st.groupBy.value == groupby  # type: ignore


def test_link_inp1(assembler_step: Step,
                   renaming_step: Step) -> None:
    """Test link inputs 1."""
    renaming_step.inpDir = assembler_step.outDir
    inpDir = renaming_step.inpDir.value.split("*")[1]
    outDir = assembler_step.outDir.value.split("&")[1]
    assert inpDir == outDir


def test_link_inp2(assembler_step: Step,
                   basic_step: Step) -> None:
    """Test link inputs 2."""
    assembler_step.imgPath = basic_step.outDir
    imgPath = assembler_step.imgPath.value.split("*")[1]
    outDir = basic_step.outDir.value.split("&")[1]
    assert imgPath == outDir


def test_link_inp3(assembler_step: Step,
                   renaming_step: Step,
                   basic_step: Step) -> None:
    """Test link inputs if already linked 1."""
    assembler_step.imgPath = basic_step.outDir
    renaming_step.inpDir = basic_step.outDir
    imgPath = assembler_step.imgPath.value.split("*")[1]
    outDir = basic_step.outDir.value.split("&")[1]
    inpDir = renaming_step.inpDir.value.split("*")[1]
    assert len(set([imgPath, outDir, inpDir])) == 1


def test_link_inp4(assembler_step: Step,
                   renaming_step: Step,
                   basic_step: Step) -> None:
    """Test link inputs if already linked 2."""
    basic_step.inpDir = renaming_step.outDir
    assembler_step.imgPath = basic_step.outDir
    assembler_step.stitchPath = basic_step.outDir
    imgPath = assembler_step.imgPath.value.split("*")[1]
    stitchPath = assembler_step.stitchPath.value.split("*")[1]
    inpDir = basic_step.inpDir.value.split("*")[1]
    outDir_basic = basic_step.outDir.value.split("&")[1]
    outDir_renaming = renaming_step.outDir.value.split("&")[1]
    assert imgPath == outDir_basic
    assert imgPath == stitchPath
    assert inpDir == outDir_renaming


def test_inp_names1(assembler_step: Step) -> None:
    """Test input names 1."""
    st = assembler_step
    assert sorted(st._input_names) == ["imgPath", "outDir", "preview", "stitchPath", "timesliceNaming"]


def test_inp_names2(renaming_step: Step) -> None:
    """Test input names 2."""
    st = renaming_step
    assert sorted(st._input_names) == ["filePattern", "inpDir", "mapDirectory", "outDir", "outFilePattern"]


def test_inp_names3(basic_step: Step) -> None:
    """Test input names 3."""
    st = basic_step
    assert sorted(st._input_names) == ["filePattern", "getDarkfield", "groupBy", "inpDir",  "outDir"]


@pytest.mark.parametrize("step_index", range(3))
def test_save_cwl(steps: list[Step], step_index: int, cwl_adapters: Path) -> None:
    """Test save yaml 1."""
    st = steps[step_index]
    if st.clt_name == "assembler":
        st.imgPath = RSRC_DIR
        st.preview = True
        st.timesliceNaming = False
    if st.clt_name == "filerenaming":
        st.filePattern = ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
        st.outFilePattern = "x{row:dd}_y{col:dd}_p{s:dd}_c{channel:d}.ome.tif"
        st.mapDirectory = "raw"
        st.inpDir = RSRC_DIR
    if st.clt_name == "basicfl":
        st.filePattern = ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
        st.getDarkfield = False
        st.groupBy = "channel"
        st.inpDir = RSRC_DIR
    st.outDir = RSRC_DIR / "out"
    with RSRC_DIR.joinpath("assembler.cwl").open("r", encoding="utf-8") as file:
        cwl_ = yaml.safe_load(file)
    st._save_cwl(cwl_adapters)
    with (cwl_adapters / "cwl_adapters" / "assembler.cwl").open("r", encoding="utf-8") as file:
        assert yaml.safe_load(file) == cwl_


def test_validate1(assembler_step: Step, tmp_path: Path) -> None:
    """Test validate 1."""
    st = assembler_step
    st.imgPath = tmp_path
    st.outDir = tmp_path / "out"
    st.stitchPath = tmp_path / "stitch"
    st._validate()


def test_validate2(assembler_step: Step, tmp_path: Path) -> None:
    """Test validate 2."""
    st = assembler_step
    st.imgPath = tmp_path
    st.stitchPath = tmp_path / "stitch"
    with pytest.raises(MissingRequiredValueError):
        st._validate()


def test_wf1(basic_step: Step,
             assembler_step: Step,
             renaming_step: Step,
             ) -> None:
    """Test WF1"""
    with pytest.raises(InvalidStepError):
        Workflow([basic_step, assembler_step, renaming_step],
                 name="wf1")


def test_wf2(basic_step: Step,
             assembler_step: Step,
             renaming_step: Step,
             tmp_path: Path) -> None:
    """Test WF2"""
    basic_step.inpDir = tmp_path / "in"
    basic_step.filePattern = ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
    basic_step.getDarkfield = False
    assembler_step.imgPath = basic_step.outDir
    assembler_step.stitchPath = tmp_path / "stitch"
    renaming_step.inpDir = assembler_step.outDir
    with pytest.raises(InvalidStepError):
        Workflow([basic_step, assembler_step, renaming_step],
                 name="wf1")


def test_wf3(basic_step: Step,
             assembler_step: Step,
             renaming_step: Step,
             tmp_path: Path) -> None:
    """Test WF3"""
    basic_step.inpDir = tmp_path / "in"
    basic_step.filePattern = ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
    basic_step.getDarkfield = False
    assembler_step.imgPath = basic_step.outDir
    assembler_step.stitchPath = tmp_path / "stitch"
    renaming_step.inpDir = assembler_step.outDir
    renaming_step.filePattern = ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
    renaming_step.outFilePattern = "x{row:dd}_y{col:dd}_p{s:dd}_c{channel:d}.ome.tif"
    renaming_step.outDir = tmp_path / "out"
    wf = Workflow([basic_step, assembler_step, renaming_step],
                  name="wf1")
    assert wf.name == "wf1"


def test_wf4(basic_step: Step,
             assembler_step: Step,
             renaming_step: Step,
             tmp_path: Path) -> None:
    """Test WF4"""
    basic_step.inpDir = tmp_path / "in"
    basic_step.filePattern = ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
    basic_step.getDarkfield = False
    assembler_step.imgPath = basic_step.outDir
    # assembler_step.stitchPath = tmp_path / "stitch"
    renaming_step.inpDir = assembler_step.outDir
    renaming_step.filePattern = ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
    renaming_step.outFilePattern = "x{row:dd}_y{col:dd}_p{s:dd}_c{channel:d}.ome.tif"
    renaming_step.outDir = tmp_path / "out"
    with pytest.raises(InvalidStepError):
        Workflow([basic_step, assembler_step, renaming_step],
                 name="wf1")


def test_compile_1(basic_step: Step,
                   assembler_step: Step,
                   renaming_step: Step,
                   tmp_path: Path) -> None:
    """Test Compile 1."""
    basic_step.inpDir = tmp_path / "in"
    basic_step.filePattern = ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
    basic_step.getDarkfield = False
    assembler_step.imgPath = basic_step.outDir
    assembler_step.stitchPath = tmp_path / "stitch"
    renaming_step.inpDir = assembler_step.outDir
    renaming_step.filePattern = ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
    renaming_step.outFilePattern = "x{row:dd}_y{col:dd}_p{s:dd}_c{channel:d}.ome.tif"
    renaming_step.outDir = tmp_path / "out"
    wf = Workflow([basic_step, assembler_step, renaming_step],
                  name="wf1")

    compiled = wf.compile()
    assert wf.yml_path == _WIC_PATH.joinpath("wf1.yml")
    assert _WIC_PATH.joinpath("autogenerated/wf1.cwl") == compiled
    assert compiled.exists()


def test_compile_2(omeconverter_step: Step,
                   montage_step: Step,
                   tmp_path: Path) -> None:
    """Test Compile 2 with new Steps."""
    ome_converter = omeconverter_step
    montage = montage_step
    ome_converter.inpDir = tmp_path / "in"
    ome_converter.filePattern = ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
    ome_converter.fileExtension = "default"
    montage.inpDir = ome_converter.outDir
    montage.filePattern = ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
    montage.outDir = tmp_path / "out"
    wf = Workflow([ome_converter, montage],
                  name="wf2")

    compiled = wf.compile()
    assert wf.yml_path == _WIC_PATH.joinpath("wf2.yml")
    assert _WIC_PATH.joinpath("autogenerated/wf2.cwl") == compiled
    assert compiled.exists()
