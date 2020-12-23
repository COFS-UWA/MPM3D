from abaqus import *
from abaqusConstants import *
from caeModules import *

# model parameters
h1 = 1.44
h2 = 0.46
h3 = 2.92
h4 = 4.00
r1 = 5.375
r2 = 0.351

# create spudcan model
md = mdb.models['Model-1']
spudcan_skt = md.ConstrainedSketch(name = '__profile__', sheetSize = 20.0)
spudcan_skt_geo = spudcan_skt.geometry
spudcan_skt_vrt = spudcan_skt.vertices
spudcan_skt.ConstructionLine(point1 = (0.0, -100.0), point2 = (0.0, 100.0))
spudcan_skt.FixedConstraint(entity = spudcan_skt_geo[2])
spudcan_skt.Line(point1 = (0.0, 0.0), point2 = (r1, h1))
spudcan_skt.CoincidentConstraint(entity1 = spudcan_skt_vrt[0], \
    entity2 = spudcan_skt_geo[2], addUndoState = False)
spudcan_skt.Line(point1 = (r1, h1), point2 = (r1, h1+h2))
spudcan_skt.VerticalConstraint(entity = spudcan_skt_geo[4], addUndoState = False)
spudcan_skt.Line(point1 = (r1, h1+h2), point2 = (r2, h1+h2+h3))
spudcan_skt.Line(point1 = (r2, h1+h2+h3), point2 = (r2, h1+h2+h3+h4))
spudcan_skt.VerticalConstraint(entity = spudcan_skt_geo[6], addUndoState = False)
spudcan_skt.Line(point1 = (r2, h1+h2+h3+h4), point2 = (0.0, h1+h2+h3+h4))
spudcan_skt.HorizontalConstraint(entity = spudcan_skt_geo[7], addUndoState = False)
spudcan_skt.PerpendicularConstraint(entity1 = spudcan_skt_geo[6], \
    entity2 = spudcan_skt_geo[7], addUndoState = False)
spudcan_skt.CoincidentConstraint(entity1 = spudcan_skt_vrt[5], \
    entity2 = spudcan_skt_geo[2], addUndoState = False)
spudcan_skt.Line(point1 = (0.0, h1+h2+h3+h4), point2 = (0.0, 0.0))
spudcan_skt.VerticalConstraint(entity = spudcan_skt_geo[8], addUndoState = False)
spudcan_skt.PerpendicularConstraint(entity1 = spudcan_skt_geo[7], \
    entity2 = spudcan_skt_geo[8], addUndoState = False)
spudcan_prt = md.Part(name = 'Part-1', dimensionality = THREE_D, type = DEFORMABLE_BODY)
spudcan_prt.BaseSolidRevolve(sketch = spudcan_skt, angle = 360.0, flipRevolveDirection=OFF)
del spudcan_skt

spudcan_prt.DatumAxisByPrincipalAxis(principalAxis = ZAXIS)
spudcan_prt.DatumPointByCoordinate(coords = (0.0, 0.0, 0.0))
datum_keys = spudcan_prt.datums.keys()
z_axis = spudcan_prt.datums[datum_keys[1]] # z axis
xpt1 = spudcan_prt.datums[datum_keys[2]] # x point 1
# cut
cl_tmp = spudcan_prt.cells.getByBoundingBox( \
    xMin = -r1*1.01, xMax = r1*1.01, \
    yMin = -(h1+h2+h3+h4)*0.01, yMax = (h1+h2+h3+h4)*1.01, \
    zMin = -r1*1.01, zMax = r1*1.01)
spudcan_prt.PartitionCellByPlanePointNormal(point = xpt1, normal = z_axis, cells = cl_tmp)

spudcan_prt.seedPart(size=0.5, deviationFactor=0.1, minSizeFactor=0.1)
cl_tmp = spudcan_prt.cells.getByBoundingBox( \
    xMin = -r1*1.01, xMax = r1*1.01, \
    yMin = -(h1+h2+h3+h4)*0.01, yMax = (h1+h2+h3+h4)*1.01, \
    zMin = -r1*1.01, zMax = r1*1.01)
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

job = mdb.Job(name='spudcan_model', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, 
    numDomains=1, activateLoadBalancing=False, multiprocessingMode=DEFAULT, 
    numCpus=1)
job.writeInput(consistencyChecking=OFF)

mdb.saveAs(pathName='./spudcan_model')