from abaqus import *
from abaqusConstants import *
from caeModules import *

# model parameters
diameter = 2.0
height = 2.0
mesh_size = diameter * 0.1
md_name = 'cylinder_2x2_model'

#
radius = diameter * 0.5
# create cylinder model
md = mdb.models['Model-1']
cylinder_skt = md.ConstrainedSketch(name = '__profile__', sheetSize = 20.0)
cylinder_skt_geo = cylinder_skt.geometry
cylinder_skt_vrt = cylinder_skt.vertices
cylinder_skt.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
#s1.assignCenterline(line=g[3])
cylinder_skt.FixedConstraint(entity=cylinder_skt_geo[2])
cylinder_skt.rectangle(point1=(0.0, 0.0), point2=(diameter*0.5, height))
cylinder_skt.CoincidentConstraint(entity1=cylinder_skt_vrt[0], entity2=cylinder_skt_geo[2], addUndoState=False)
cylinder_prt = md.Part(name = 'Part-1', dimensionality = THREE_D, type = DEFORMABLE_BODY)
cylinder_prt.BaseSolidRevolve(sketch = cylinder_skt, angle = 360.0, flipRevolveDirection=OFF)
del cylinder_skt

cylinder_prt.DatumAxisByPrincipalAxis(principalAxis = ZAXIS)
cylinder_prt.DatumPointByCoordinate(coords = (0.0, 0.0, 0.0))
datum_keys = cylinder_prt.datums.keys()
z_axis = cylinder_prt.datums[datum_keys[1]] # z axis
xpt1 = cylinder_prt.datums[datum_keys[2]] # x point 1
# cut
cl_tmp = cylinder_prt.cells.getByBoundingBox( \
    xMin = -radius*1.01, xMax = radius*1.01, \
    yMin = -height*0.01, yMax = height*1.01, \
    zMin = -radius*1.01, zMax = radius*1.01)
cylinder_prt.PartitionCellByPlanePointNormal(point = xpt1, normal = z_axis, cells = cl_tmp)

cylinder_prt.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
cl_tmp = cylinder_prt.cells.getByBoundingBox( \
    xMin = -radius*1.01, xMax = radius*1.01, \
    yMin = -height*0.01, yMax = height*1.01, \
    zMin = -radius*1.01, zMax = radius*1.01)
cylinder_prt.setMeshControls(regions = cl_tmp, elemShape = TET, technique = FREE)
elemType1 = mesh.ElemType(elemCode = C3D8R, elemLibrary = STANDARD)
elemType2 = mesh.ElemType(elemCode = C3D6, elemLibrary = STANDARD)
elemType3 = mesh.ElemType(elemCode = C3D4, elemLibrary = STANDARD, 
    secondOrderAccuracy = OFF, distortionControl = DEFAULT)
cylinder_prt.setElementType(regions = (cl_tmp, ), \
    elemTypes=(elemType1, elemType2, elemType3))

cylinder_prt.generateMesh()

assembly = md.rootAssembly
assembly.DatumCsysByDefault(CARTESIAN)
assembly.Instance(name = 'Inst-1', part = cylinder_prt, dependent = ON)
#assembly.rotate(instanceList=('Inst-1', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0), angle=90.0)

job = mdb.Job(name=md_name, model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, 
    numDomains=1, activateLoadBalancing=False, multiprocessingMode=DEFAULT, 
    numCpus=1)
job.writeInput(consistencyChecking=OFF)

mdb.saveAs(pathName='./'+md_name)