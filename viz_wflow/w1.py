"""Visualization Workflow."""
from pathlib import Path

from wic.api import Step, Workflow

PATH = Path(__file__).absolute()

assemble = Step(PATH.with_name("assemble.cwl"))
filerenaming = Step(PATH.with_name("filerenaming.cwl"))
montage = Step(PATH.with_name("montage.cwl"))
omeconverter = Step(PATH.with_name("omeconverter.cwl"))
precomputeslide = Step(PATH.with_name("precomputeslide.cwl"))

omeconverter.inpDir = Path("/Users/camilovelezr/axle/workflow/hamdah/417744")
omeconverter.filePattern = ".*.tif"
omeconverter.fileExtension = ".ome.tif"
# omeconverter.outDir = "&converted"
filerenaming.inpDir = omeconverter.outDir
filerenaming.filePattern = ".*_{row:c}{col:dd}_s{s:d}_w{channel:d}.*.ome.tif"
filerenaming.outFilePattern = "x{row:dd}_y{col:dd}_p{s:dd}_c{channel:d}.ome.tif"

montage.inpDir = filerenaming.outDir
montage.filePattern = "x{x:d+}_y{y:d+}_p{p:d+}_c{c:d}.ome.tif"
montage.layout = "p,yx"
montage.imageSpacing = "1"
montage.gridSpacing = "20"

assemble.stitchPath = montage.outDir
assemble.imgPath = filerenaming.outDir
assemble.timesliceNaming = False

precomputeslide.inpDir = assemble.outDir
precomputeslide.pyramidType = "Zarr"
precomputeslide.imageType = "image"
precomputeslide.filePattern = "x(00-15)_y(01-24)_p01_c{c}.ome.tif"
precomputeslide.outDir = Path("/Users/camilovelezr/axle/workflow/hamdah")

wf = Workflow(
    [omeconverter, filerenaming, montage, assemble, precomputeslide],
    "firstwfapi",
    wic_dir=Path("/Users/camilovelezr/workflow-inference-compiler"),
)
wf.compile(out_dir=Path("/Users/camilovelezr/axle/workflow/hamdah"))
