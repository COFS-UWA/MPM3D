from abaqus import *
from abaqusConstants import *
from caeModules import *

blk_len = 1.0

md = mdb.models['Model-1']

spudcan_skt = md.ConstrainedSketch(name='__profile__', sheetSize=10.0)
spudcan_skt.Line(point1=(0.0, 0.0), point2=(blk_len*3.0, 0.0))
spudcan_skt.Line(point1=(blk_len*3.0, 0.0), point2=(blk_len*3.0, blk_len*2.0))
spudcan_skt.Line(point1=(blk_len*3.0, blk_len*2.0), point2=(blk_len*4.0, blk_len*2.0))
spudcan_skt.Line(point1=(blk_len*4.0, blk_len*2.0), point2=(blk_len*4.0, 0.0))
spudcan_skt.Line(point1=(blk_len*4.0, 0.0), point2=(blk_len*5.0, 0.0))
spudcan_skt.Line(point1=(blk_len*5.0, 0.0), point2=(blk_len*5.0, blk_len*3.0))
spudcan_skt.Line(point1=(blk_len*5.0, blk_len*3.0), point2=(blk_len*2.0, blk_len*3.0))
spudcan_skt.Line(point1=(blk_len*2.0, blk_len*3.0), point2=(blk_len*2.0, blk_len))
spudcan_skt.Line(point1=(blk_len*2.0, blk_len), point2=(blk_len, blk_len))
spudcan_skt.Line(point1=(blk_len, blk_len), point2=(blk_len, blk_len*3.0))
spudcan_skt.Line(point1=(blk_len, blk_len*3.0), point2=(0.0, blk_len*3.0))
spudcan_skt.Line(point1=(0.0, blk_len*3.0), point2=(0.0, 0.0))

spudcan_prt = md.Part(name = 'Part-1', dimensionality = THREE_D, type = DEFORMABLE_BODY)
spudcan_prt.BaseSolidExtrude(sketch = spudcan_skt, depth=blk_len*3.0)
del spudcan_skt

spudcan_prt.DatumAxisByPrincipalAxis(principalAxis = XAXIS) # 1
spudcan_prt.DatumPointByCoordinate(coords = (blk_len, 0.0, 0.0)) # 2
spudcan_prt.DatumPointByCoordinate(coords = (blk_len*2.0, 0.0, 0.0))
spudcan_prt.DatumPointByCoordinate(coords = (blk_len*3.0, 0.0, 0.0))
spudcan_prt.DatumPointByCoordinate(coords = (blk_len*4.0, 0.0, 0.0))
datum_keys = spudcan_prt.datums.keys()
x_axis = spudcan_prt.datums[datum_keys[0]]
cpt1 = spudcan_prt.datums[datum_keys[1]]
cpt2 = spudcan_prt.datums[datum_keys[2]]
cpt3 = spudcan_prt.datums[datum_keys[3]]
cpt4 = spudcan_prt.datums[datum_keys[4]]

# cut1
cl_tmp = spudcan_prt.cells.getByBoundingBox( \
    xMin = -blk_len*0.01, xMax = blk_len*5.01, \
    yMin = -blk_len*0.01, yMax = blk_len*3.01, \
    zMin = -blk_len*0.01, zMax = blk_len*3.01)
spudcan_prt.PartitionCellByPlanePointNormal(point = cpt1, normal = x_axis, cells = cl_tmp)
# cut2
cl_tmp = spudcan_prt.cells.getByBoundingBox( \
    xMin = blk_len*(1.0-0.01), xMax = blk_len*5.01, \
    yMin = -blk_len*0.01, yMax = blk_len*3.01, \
    zMin = -blk_len*0.01, zMax = blk_len*3.01)
spudcan_prt.PartitionCellByPlanePointNormal(point = cpt2, normal = x_axis, cells = cl_tmp)
# cut3
cl_tmp = spudcan_prt.cells.getByBoundingBox( \
    xMin = blk_len*(2.0-0.01), xMax = blk_len*5.01, \
    yMin = -blk_len*0.01, yMax = blk_len*3.01, \
    zMin = -blk_len*0.01, zMax = blk_len*3.01)
spudcan_prt.PartitionCellByPlanePointNormal(point = cpt3, normal = x_axis, cells = cl_tmp)
# cut4
cl_tmp = spudcan_prt.cells.getByBoundingBox( \
    xMin = blk_len*(3.0-0.01), xMax = blk_len*5.01, \
    yMin = -blk_len*0.01, yMax = blk_len*3.01, \
    zMin = -blk_len*0.01, zMax = blk_len*3.01)
spudcan_prt.PartitionCellByPlanePointNormal(point = cpt4, normal = x_axis, cells = cl_tmp)

spudcan_prt.seedPart(size = blk_len*0.5, deviationFactor = 0.1, minSizeFactor = 0.1)
cl_tmp = spudcan_prt.cells.getByBoundingBox( \
    xMin = -blk_len*0.01, xMax = blk_len*5.01, \
    yMin = -blk_len*0.01, yMax = blk_len*3.01, \
    zMin = -blk_len*0.01, zMax = blk_len*3.01)
spudcan_prt.setMeshControls(regions = cl_tmp, elemShape = TET, technique = FREE)
elemType1 = mesh.ElemType(elemCode = C3D8R, elemLibrary = STANDARD)
elemType2 = mesh.ElemType(elemCode = C3D6, elemLibrary = STANDARD)
elemType3 = mesh.ElemType(elemCode = C3D4, elemLibrary = STANDARD, 
    secondOrderAccuracy = OFF, distortionControl = DEFAULT)
spudcan_prt.setElementType(regions = (cl_tmp, ), \
    elemTypes=(elemType3, ))

spudcan_prt.generateMesh()

assembly = md.rootAssembly
assembly.DatumCsysByDefault(CARTESIAN)
assembly.Instance(name = 'Inst-1', part = spudcan_prt, dependent = ON)

job = mdb.Job(name='weird_block', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, 
    numDomains=1, activateLoadBalancing=False, multiprocessingMode=DEFAULT, 
    numCpus=1)
job.writeInput(consistencyChecking=OFF)

mdb.saveAs(pathName='./weird_block')