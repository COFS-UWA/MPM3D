from abaqus import *
from abaqusConstants import *
from caeModules import *

# model parameters
# diameter = 3 m
pt1 = (0.0, 0.0)
pt2 = (0.3375, 0.43)
pt3 = (1.5, 0.7)
pt4 = (1.5, 0.775)
pt5 = (0.365, 1.045)
pt6 = (0.365, 4.5)
pt7 = (0.0, 4.5)
mesh_size = 0.1

# Build model
# top
h1 = pt2[1]
h2 = pt3[1]
h3 = pt4[1]
h4 = pt5[1]
h5 = pt6[1]
r1 = abs(pt2[0])
r2 = abs(pt3[0])
r3 = abs(pt5[0])

# create spudcan model
md = mdb.models['Model-1']
spudcan_skt = md.ConstrainedSketch(name = '__profile__', sheetSize = 10.0)
spudcan_skt.ConstructionLine(point1=(0.0, -5.0), point2=(0.0, 5.0))
spudcan_skt.Line(point1 = pt1, point2 = pt2)
spudcan_skt.Line(point1 = pt2, point2 = pt3)
spudcan_skt.Line(point1 = pt3, point2 = pt4)
spudcan_skt.Line(point1 = pt4, point2 = pt5)
spudcan_skt.Line(point1 = pt5, point2 = pt6)
spudcan_skt.Line(point1 = pt6, point2 = pt7)
spudcan_skt.Line(point1 = pt7, point2 = pt1)
spudcan_prt = md.Part(name = 'Part-1', dimensionality = THREE_D, type = DEFORMABLE_BODY)
spudcan_prt.BaseSolidRevolve(sketch = spudcan_skt, angle = 360.0, flipRevolveDirection = OFF)
del spudcan_skt

spudcan_prt.DatumAxisByPrincipalAxis(principalAxis = YAXIS)
spudcan_prt.DatumPointByCoordinate(coords = (0.0, pt2[1], 0.0))
spudcan_prt.DatumPointByCoordinate(coords = (0.0, pt3[1], 0.0))
spudcan_prt.DatumPointByCoordinate(coords = (0.0, pt4[1], 0.0))
spudcan_prt.DatumPointByCoordinate(coords = (0.0, pt5[1], 0.0))
datum_keys = spudcan_prt.datums.keys()
y_axis = spudcan_prt.datums[datum_keys[1]]
cpt1 = spudcan_prt.datums[datum_keys[2]]
cpt2 = spudcan_prt.datums[datum_keys[3]]
cpt3 = spudcan_prt.datums[datum_keys[4]]
cpt4 = spudcan_prt.datums[datum_keys[5]]
# cut
cl_tmp = spudcan_prt.cells.getByBoundingBox( \
    xMin = -r2*1.01, xMax = r2*1.01, \
    yMin = -h5*0.01, yMax = h5*1.01, \
    zMin = -r2*1.01, zMax = r2*1.01)
spudcan_prt.PartitionCellByPlanePointNormal(point = cpt1, normal = y_axis, cells = cl_tmp)
#
cl_tmp = spudcan_prt.cells.getByBoundingBox( \
    xMin = -r2*1.01, xMax = r2*1.01, \
    yMin = h1 - h5*0.01, yMax = h5*1.01, \
    zMin = -r2*1.01, zMax = r2*1.01)
spudcan_prt.PartitionCellByPlanePointNormal(point = cpt2, normal = y_axis, cells = cl_tmp)
#
cl_tmp = spudcan_prt.cells.getByBoundingBox( \
    xMin = -r2*1.01, xMax = r2*1.01, \
    yMin = h2 - h5*0.01, yMax = h5*1.01, \
    zMin = -r2*1.01, zMax = r2*1.01)
spudcan_prt.PartitionCellByPlanePointNormal(point = cpt3, normal = y_axis, cells = cl_tmp)
#
cl_tmp = spudcan_prt.cells.getByBoundingBox( \
    xMin = -r2*1.01, xMax = r2*1.01, \
    yMin = h3 - h5*0.01, yMax = h5*1.01, \
    zMin = -r2*1.01, zMax = r2*1.01)
spudcan_prt.PartitionCellByPlanePointNormal(point = cpt4, normal = y_axis, cells = cl_tmp)

spudcan_prt.seedPart(size = mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
cl_tmp = spudcan_prt.cells.getByBoundingBox( \
    xMin = -r2*1.01, xMax = r2*1.01, \
    yMin = -h5*0.01, yMax = h5*1.01, \
    zMin = -r2*1.01, zMax = r2*1.01)
spudcan_prt.setMeshControls(regions = cl_tmp, elemShape = TET, technique = FREE)
elemType1 = mesh.ElemType(elemCode = C3D8R, elemLibrary = STANDARD)
elemType2 = mesh.ElemType(elemCode = C3D6, elemLibrary = STANDARD)
elemType3 = mesh.ElemType(elemCode = C3D4, elemLibrary = STANDARD, 
    secondOrderAccuracy = OFF, distortionControl = DEFAULT)
spudcan_prt.setElementType(regions = (cl_tmp, ), \
    elemTypes=(elemType1, elemType2, elemType3))

spudcan_prt.generateMesh()

assembly = md.rootAssembly
assembly.DatumCsysByDefault(CARTESIAN)
assembly.Instance(name = 'Inst-1', part = spudcan_prt, dependent = ON)

job = mdb.Job(name='spudcan_model_Hossain_2006', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, 
    numDomains=1, activateLoadBalancing=False, multiprocessingMode=DEFAULT, 
    numCpus=1)
job.writeInput(consistencyChecking=OFF)

mdb.saveAs(pathName='./spudcan_model_Hossain_2006')