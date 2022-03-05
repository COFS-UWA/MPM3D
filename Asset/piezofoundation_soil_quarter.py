# soil block of spudcan
# only one quarter
from abaqus import *
from abaqusConstants import *
from caeModules import *

spudcan_radius = 1.0

top = 0.5 * spudcan_radius
depth = 7.0 * spudcan_radius
dense_depth = 3.5 * spudcan_radius
coarse_depth = 3.8 * spudcan_radius

width = 5.5 * spudcan_radius
dense_width = 3.0 * spudcan_radius
coarse_width = 3.3 * spudcan_radius

dense_elem_size = 0.125 * spudcan_radius
coarse_elem_size = 0.25 * spudcan_radius

md = mdb.models['Model-1']

# Extrusion 1
rect_skt = md.ConstrainedSketch(name = 'rect_extrude', sheetSize = 10.0)
rect_skt.rectangle(point1 = (0.0, 0.0), point2 = (width, width))
prt1 = md.Part(name = 'Part-1', dimensionality = THREE_D, type = DEFORMABLE_BODY)
prt1.BaseSolidExtrude(sketch = rect_skt, depth = top)
del rect_skt

# Extrusion 2
prt1_faces = prt1.faces
prt1_edges = prt1.edges
prt1_tf = prt1.MakeSketchTransform(sketchPlane = prt1_faces[5], \
    sketchUpEdge = prt1_edges[2], sketchPlaneSide = SIDE1, \
    sketchOrientation = RIGHT, origin = (0.0, 0.0, 0.0))
rect2_skt = md.ConstrainedSketch(name = 'rect2_extrude', sheetSize = 10.0, transform = prt1_tf)
prt1.projectReferencesOntoSketch(sketch = rect2_skt, filter = COPLANAR_EDGES)
rect2_skt.rectangle(point1 = (0.0, 0.0), point2 = (-width, width))
prt1.SolidExtrude(sketchPlane = prt1_faces[5], sketchUpEdge = prt1_edges[2], \
    sketchPlaneSide = SIDE1, sketchOrientation = RIGHT,
    sketch = rect2_skt, depth = depth, flipExtrudeDirection = OFF)
del rect2_skt

# Datum point
prt1.DatumAxisByPrincipalAxis(principalAxis = ZAXIS)
prt1.DatumPointByCoordinate(coords = (0.0, 0.0, 0.0))
prt1.DatumPointByCoordinate(coords = (0.0, 0.0, -dense_depth))
prt1.DatumPointByCoordinate(coords = (0.0, 0.0, -coarse_depth))
prt1.DatumAxisByPrincipalAxis(principalAxis = XAXIS)
prt1.DatumPointByCoordinate(coords = (dense_width, 0.0, 0.0))
prt1.DatumPointByCoordinate(coords = (coarse_width, 0.0, 0.0))
prt1.DatumAxisByPrincipalAxis(principalAxis = YAXIS)
prt1.DatumPointByCoordinate(coords = (0.0, dense_width, 0.0))
prt1.DatumPointByCoordinate(coords = (0.0, coarse_width, 0.0))

datum_keys = prt1.datums.keys()
z_axis = prt1.datums[datum_keys[0]] # z axis
zpt1 = prt1.datums[datum_keys[1]] # z point 1
zpt2 = prt1.datums[datum_keys[2]] # z point 2
zpt3 = prt1.datums[datum_keys[3]] # z point 3
x_axis = prt1.datums[datum_keys[4]] # x axis
xpt1 = prt1.datums[datum_keys[5]] # x point 1
xpt2 = prt1.datums[datum_keys[6]] # x point 2
y_axis = prt1.datums[datum_keys[7]] # y axis
ypt1 = prt1.datums[datum_keys[8]] # y point 1
ypt2 = prt1.datums[datum_keys[9]] # y point 2

# z cut 1
cl_tmp = prt1.cells.getByBoundingBox(xMin = -width*0.01, xMax = width*1.01, \
    yMin = -width*0.01, yMax = width*1.01, zMin = -depth*1.01, zMax = top*1.01)
prt1.PartitionCellByPlanePointNormal(point = zpt1, normal = z_axis, cells = cl_tmp)
# z cut 2
cl_tmp = prt1.cells.getByBoundingBox(xMin = -width*0.01, xMax = width*1.01, \
    yMin = -width*0.01, yMax = width*1.01, zMin = -depth*1.01, zMax = top*0.01)
prt1.PartitionCellByPlanePointNormal(point = zpt2, normal = z_axis, cells = cl_tmp)
# z cut 3
cl_tmp = prt1.cells.getByBoundingBox(xMin = -width*0.01, xMax = width*1.01, \
    yMin = -width*0.01, yMax = width*1.01, zMin = -depth*1.01, zMax = -dense_depth*0.99)
prt1.PartitionCellByPlanePointNormal(point = zpt3, normal = z_axis, cells = cl_tmp)
# x cut 1
cl_tmp = prt1.cells.getByBoundingBox(xMin = -width*0.01, xMax = width*1.01, \
    yMin = -width*0.01, yMax = width*1.01, zMin = -depth*1.01, zMax = top*1.01)
prt1.PartitionCellByPlanePointNormal(point = xpt1, normal = x_axis, cells = cl_tmp)
# x cut 2
cl_tmp = prt1.cells.getByBoundingBox(xMin = -width*0.01, xMax = width*1.01, \
    yMin = -width*0.01, yMax = width*1.01, zMin = -depth*1.01, zMax = top*1.01)
prt1.PartitionCellByPlanePointNormal(point = xpt2, normal = x_axis, cells = cl_tmp)
# y cut 1
cl_tmp = prt1.cells.getByBoundingBox(xMin = -width*0.01, xMax = width*1.01, \
    yMin = -width*0.01, yMax = width*1.01, zMin = -depth*1.01, zMax = top*1.01)
prt1.PartitionCellByPlanePointNormal(point = ypt1, normal = y_axis, cells = cl_tmp)
# y cut 2
cl_tmp = prt1.cells.getByBoundingBox(xMin = -width*0.01, xMax = width*1.01, \
    yMin = -width*0.01, yMax = width*1.01, zMin = -depth*1.01, zMax = top*1.01)
prt1.PartitionCellByPlanePointNormal(point = ypt2, normal = y_axis, cells = cl_tmp)

# seed globally
prt1.seedPart(size = coarse_elem_size, deviationFactor = 0.1, minSizeFactor = 0.1)
all_cells = prt1.cells.getByBoundingBox(xMin = -width*0.01, xMax = width*1.01, \
    yMin = -width*0.01, yMax = width*1.01, zMin = -depth*1.01, zMax = top*1.01)
prt1.setMeshControls(regions = all_cells, elemShape = TET, technique = FREE, allowMapped = False)
#elem_type1 = mesh.ElemType(elemCode = C3D8R, elemLibrary = STANDARD)
#elem_type2 = mesh.ElemType(elemCode = C3D6, elemLibrary = STANDARD)
elem_type3 = mesh.ElemType(elemCode = C3D4, elemLibrary = STANDARD, \
    secondOrderAccuracy = OFF, distortionControl = DEFAULT)
ac_set = prt1.Set(name = 'AllCells', cells = all_cells)
prt1.setElementType(regions = ac_set, elemTypes = (elem_type3, ))
# seed edge
seed_edg = prt1.edges.findAt(
    # horizontal lines
    ((0.0, dense_width*0.5, top), ),
    ((dense_width, dense_width*0.5, top), ),
    ((dense_width*0.5, 0.0, top), ),
    ((dense_width*0.5, dense_width, top), ),
    ((0.0, dense_width*0.5, 0.0), ),
    ((dense_width, dense_width*0.5, 0.0), ),
    ((dense_width*0.5, 0.0, 0.0), ),
    ((dense_width*0.5, dense_width, 0.0), ),
    ((0.0, dense_width*0.5, -dense_depth), ),
    ((dense_width, dense_width*0.5, -dense_depth), ),
    ((dense_width*0.5, 0.0, -dense_depth), ),
    ((dense_width*0.5, dense_width, -dense_depth), ),
    # vertical lines
    ((dense_width, dense_width, top*0.5), ),
    ((dense_width, 0.0, top*0.5), ),
    ((0.0, 0.0, top*0.5), ),
    ((0.0, dense_width, top*0.5), ),
    ((dense_width, dense_width, -dense_depth*0.5), ),
    ((dense_width, 0.0, -dense_depth*0.5), ),
    ((0.0, 0.0, -dense_depth*0.5), ),
    ((0.0, dense_width, -dense_depth*0.5), )
    )
seed_edge_set = prt1.Set(name = 'SeedEdges', edges = seed_edg)
prt1.seedEdgeBySize(edges = seed_edg, size = dense_elem_size, \
    deviationFactor = 0.1, minSizeFactor = 0.1, constraint = FIXED)

prt1.generateMesh()

assembly = md.rootAssembly
assembly.DatumCsysByDefault(CARTESIAN)
assembly.Instance(name = 'Inst-1', part = prt1, dependent = ON)

job = mdb.Job(name = 'piezofoundation_soil_quarter', model = 'Model-1',
    description = '', type = ANALYSIS, atTime = None, waitMinutes = 0,
    waitHours = 0, queue = None, memory = 90, memoryUnits = PERCENTAGE,
    getMemoryFromAnalysis = True, explicitPrecision = SINGLE,
    nodalOutputPrecision = SINGLE, echoPrint = OFF, modelPrint = OFF,
    contactPrint = OFF, historyPrint=OFF, userSubroutine = '', scratch = '',
    resultsFormat = ODB, parallelizationMethodExplicit = DOMAIN, numDomains = 1,
    activateLoadBalancing = False, multiprocessingMode = DEFAULT, numCpus = 1)
job.writeInput(consistencyChecking = OFF)

mdb.saveAs(pathName = './piezofoundation_soil_quarter.cae')
