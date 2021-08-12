from abaqus import *
from abaqusConstants import *
from caeModules import *

# model parameters
diameter = 2.0
height = 2.0
mesh_size = diameter * 0.1
md_name = 'cylinder_2x2_quad_model'

#
diameter2 = diameter
radius2 = diameter2 * 0.5
diameter += mesh_size * 2.0
radius = diameter * 0.5
# create cylinder model
md = mdb.models['Model-1']
cylinder_skt = md.ConstrainedSketch(name = '__profile__', sheetSize = 20.0)
cylinder_skt_geo = cylinder_skt.geometry
cylinder_skt_vrt = cylinder_skt.vertices
cylinder_skt.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
#s1.assignCenterline(line=g[3])
cylinder_skt.FixedConstraint(entity=cylinder_skt_geo[2])
cylinder_skt.rectangle(point1=(0.0, 0.0), point2=(radius, height))
cylinder_skt.CoincidentConstraint(entity1=cylinder_skt_vrt[0], entity2=cylinder_skt_geo[2], addUndoState=False)
cylinder_prt = md.Part(name = 'Part-1', dimensionality = THREE_D, type = DEFORMABLE_BODY)
cylinder_prt.BaseSolidRevolve(sketch = cylinder_skt, angle = 90.0, flipRevolveDirection=OFF)
del cylinder_skt

# cut
p1 = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['Part-1']
f1, e1, d1 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f1[0], sketchUpEdge=e1[0], 
    sketchPlaneSide=SIDE1, origin=(0.0, height, 0.0))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=4.89, 
    gridSpacing=0.12, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
#s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-1']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
s.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, radius2), 
    point2=(-radius2, 0.0), direction=COUNTERCLOCKWISE)
s.CoincidentConstraint(entity1=v[3], entity2=g[2], addUndoState=False)
s.CoincidentConstraint(entity1=v[4], entity2=g[4], addUndoState=False)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
e, d2 = p.edges, p.datums
p.PartitionFaceBySketch(sketchUpEdge=e[0], faces=pickedFaces, sketch=s)
#s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-1']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
e1, d1 = p.edges, p.datums
pickedEdges =(e1[0], )
p.PartitionCellByExtrudeEdge(line=d1[1], cells=pickedCells, edges=pickedEdges, 
    sense=FORWARD)

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