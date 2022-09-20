import math
from abaqus import *
from abaqusConstants import *
from caeModules import *

# input parameters
footing_radius = 1.5
cae_name = 'spudcan_soil_half_cylinder_8R'

cy_radius = 8.0 * footing_radius
cy_dense_radius = 3.0 * footing_radius
cy_coarse_radius = 3.5 * footing_radius
cy_top = 0.5 * footing_radius
cy_depth = 8.0 * footing_radius
cy_dense_depth = 3.5 * footing_radius
cy_coarse_depth = 4.0 * footing_radius

dense_elem_size = 0.125 * footing_radius
coarse_elem_size = 0.25 * footing_radius

# build cae model
half_sqrt_2 = 0.5*math.sqrt(2.0)
cy_len = cy_top + cy_depth

cy_md = mdb.models['Model-1']

# Form cylinder by extrusion and create part
cy_skt = cy_md.ConstrainedSketch(name='__profile__', sheetSize=100.0)
cy_geo = cy_skt.geometry
cy_skt.ArcByCenterEnds(center=(0.0, 0.0), point1=(cy_radius, 0.0), point2=(-cy_radius, 0.0), direction=COUNTERCLOCKWISE)
cy_skt.Line(point1=(-cy_radius, 0.0), point2=(cy_radius, 0.0))
cy_prt = cy_md.Part(name='Part-1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
cy_prt.BaseSolidExtrude(sketch=cy_skt, depth=cy_len)

del cy_skt

# Draw two arcs at top face
cy_prt_faces, cy_prt_edges = cy_prt.faces, cy_prt.edges
cut_transform = cy_prt.MakeSketchTransform(sketchPlane=cy_prt_faces[2], sketchUpEdge=cy_prt_edges[4], 
    sketchPlaneSide=SIDE1, origin=(0.0, 0.0, cy_len))
cut_skt = cy_md.ConstrainedSketch(name='__profile__', sheetSize=50.0, gridSpacing=5.0, transform=cut_transform)
cy_prt.projectReferencesOntoSketch(sketch=cut_skt, filter=COPLANAR_EDGES)
cut_skt.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, -cy_coarse_radius), point2=(0.0, cy_coarse_radius), direction=COUNTERCLOCKWISE)
cut_skt.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, -cy_dense_radius), point2=(0.0, cy_dense_radius), direction=COUNTERCLOCKWISE)
#cut_skt.setPrimaryObject(option=SUPERIMPOSE)
cy_prt.PartitionFaceBySketch(sketchUpEdge=cy_prt_edges[4], faces=cy_prt_faces[2], sketch=cut_skt)
del cut_skt

cy_prt_cells = cy_prt.cells
# circular partition 1
part_sweep_dir1 = cy_prt_edges.findAt(((cy_radius, 0.0, cy_len*0.5), ), )[0]
part_sweep_edg1 = cy_prt_edges.findAt(((cy_coarse_radius*half_sqrt_2, cy_coarse_radius*half_sqrt_2, cy_len), ), )
part_cell = cy_prt_cells.findAt(((0.0, cy_radius*0.5, cy_depth*0.5), ), )
cy_prt.PartitionCellByExtrudeEdge(line=part_sweep_dir1, cells=part_cell, edges=part_sweep_edg1, sense=REVERSE)
# circular partition 2
part_sweep_dir2 = cy_prt_edges.findAt(((cy_radius, 0.0, cy_len*0.5), ), )[0]
part_sweep_edg2 = cy_prt_edges.findAt(((cy_dense_radius*half_sqrt_2, cy_dense_radius*half_sqrt_2, cy_len), ), )
part_cell = cy_prt_cells.findAt(((0.0, cy_coarse_radius*0.5, cy_depth*0.5), ), )
cy_prt.PartitionCellByExtrudeEdge(line=part_sweep_dir2, cells=part_cell, edges=part_sweep_edg2, sense=REVERSE)

# z partitions
cy_prt.DatumAxisByPrincipalAxis(principalAxis = ZAXIS)
cy_prt.DatumPointByCoordinate(coords = (0.0, 0.0, cy_depth))
cy_prt.DatumPointByCoordinate(coords = (0.0, 0.0, cy_depth - cy_dense_depth))
cy_prt.DatumPointByCoordinate(coords = (0.0, 0.0, cy_depth - cy_coarse_depth))
datum_keys = cy_prt.datums.keys()
z_axis = cy_prt.datums[datum_keys[0]] # z axis
zpt1 = cy_prt.datums[datum_keys[1]] # z point 1
zpt2 = cy_prt.datums[datum_keys[2]] # z point 2
zpt3 = cy_prt.datums[datum_keys[3]] # z point 3

# z cut 1
cl_tmp = cy_prt_cells.getByBoundingBox(xMin = -cy_radius*1.01, xMax = cy_radius*1.01, \
    yMin = -cy_radius*0.01, yMax = cy_radius*1.01, zMin = -cy_len*0.01, zMax = cy_len*1.01)
cy_prt.PartitionCellByPlanePointNormal(point = zpt1, normal = z_axis, cells = cl_tmp)
# z cut 2
cl_tmp = cy_prt_cells.getByBoundingBox(xMin = -cy_radius*1.01, xMax = cy_radius*1.01, \
    yMin = -cy_radius*0.01, yMax = cy_radius*1.01, zMin = -cy_len*0.01, zMax = cy_len*1.01)
cy_prt.PartitionCellByPlanePointNormal(point = zpt2, normal = z_axis, cells = cl_tmp)
# z cut 3
cl_tmp = cy_prt_cells.getByBoundingBox(xMin = -cy_radius*1.01, xMax = cy_radius*1.01, \
    yMin = -cy_radius*0.01, yMax = cy_radius*1.01, zMin = -cy_len*0.01, zMax = cy_len*1.01)
cy_prt.PartitionCellByPlanePointNormal(point = zpt3, normal = z_axis, cells = cl_tmp)

# mesh the part
# seed globally
cy_prt.seedPart(size = coarse_elem_size, deviationFactor = 0.1, minSizeFactor = 0.1)
all_cells = cy_prt_cells.getByBoundingBox(xMin = -cy_radius*1.01, xMax = cy_radius*1.01, \
    yMin = -cy_radius*0.01, yMax = cy_radius*1.01, zMin = -cy_len*0.01, zMax = cy_len*1.01)
cy_prt.setMeshControls(regions = all_cells, elemShape = TET, technique = FREE, allowMapped = False)
elem_type3 = mesh.ElemType(elemCode = C3D4, elemLibrary = STANDARD, secondOrderAccuracy = OFF, distortionControl = DEFAULT)
ac_set = cy_prt.Set(name = 'AllCells', cells = all_cells)
cy_prt.setElementType(regions = ac_set, elemTypes = (elem_type3, ))
# seed edge
seed_edg = cy_prt_edges.findAt(
    # arcs
    ((0.0, cy_dense_radius, cy_len), ),
    ((0.0, cy_dense_radius, cy_depth), ),
    ((0.0, cy_dense_radius, cy_depth - cy_dense_depth), ),
    # horizontal lines
    ((0.0, 0.0, cy_len), ),
    ((0.0, 0.0, cy_depth), ),
    ((0.0, 0.0, cy_depth - cy_dense_depth), ),
    # vertical lines
    (( cy_dense_radius, 0.0, cy_depth + 0.5*cy_top), ),
    ((-cy_dense_radius, 0.0, cy_depth + 0.5*cy_top), ),
    (( cy_dense_radius, 0.0, cy_depth - 0.5*cy_dense_depth), ),
    ((-cy_dense_radius, 0.0, cy_depth - 0.5*cy_dense_depth), ),
    )
#seed_edge_set = cy_prt.Set(name = 'SeedEdges', edges = seed_edg)
cy_prt.seedEdgeBySize(edges = seed_edg, size = dense_elem_size, \
    deviationFactor = 0.1, minSizeFactor = 0.1, constraint = FIXED)
# generate mesh
cy_prt.generateMesh()

assembly = cy_md.rootAssembly
assembly.DatumCsysByDefault(CARTESIAN)
assembly.Instance(name = 'Inst-1', part = cy_prt, dependent = ON)

job = mdb.Job(name = cae_name, model = 'Model-1',
    description = '', type = ANALYSIS, atTime = None, waitMinutes = 0,
    waitHours = 0, queue = None, memory = 90, memoryUnits = PERCENTAGE,
    getMemoryFromAnalysis = True, explicitPrecision = SINGLE,
    nodalOutputPrecision = SINGLE, echoPrint = OFF, modelPrint = OFF,
    contactPrint = OFF, historyPrint=OFF, userSubroutine = '', scratch = '',
    resultsFormat = ODB, parallelizationMethodExplicit = DOMAIN, numDomains = 1,
    activateLoadBalancing = False, multiprocessingMode = DEFAULT, numCpus = 1)
job.writeInput(consistencyChecking = OFF)

mdb.saveAs(pathName = './' + cae_name + '.cae')
