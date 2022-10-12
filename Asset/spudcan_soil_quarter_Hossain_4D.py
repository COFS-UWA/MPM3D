import math
from abaqus import *
from abaqusConstants import *
from caeModules import *

# input parameters
footing_radius = 1.5
cae_name = 'spudcan_soil_quarter_Hossain_4D'

cy_radius = 8.0 * footing_radius
cy_dense_radius = 3.0 * footing_radius
cy_coarse_radius = 3.5 * footing_radius
cy_top = 1.0 * footing_radius
cy_depth = 10.0 * footing_radius
cy_dense_depth = 5.0 * footing_radius
cy_coarse_depth = 5.5 * footing_radius

dense_elem_size = 0.16 * footing_radius
coarse_elem_size = 0.3 * footing_radius

# build cae model
half_sqrt_2 = 0.5*math.sqrt(2.0)
cy_len = cy_top + cy_depth

ft_md = mdb.models['Model-1']

# Form cylinder by extrusion and create part
cy_skt = ft_md.ConstrainedSketch(name='__profile__', sheetSize=100.0)
cy_geo = cy_skt.geometry
cy_skt.ArcByCenterEnds(center=(0.0, 0.0), point1=(cy_radius, 0.0), point2=(0.0, cy_radius), direction=COUNTERCLOCKWISE)
cy_skt.Line(point1=(0.0, cy_radius), point2=(0.0, 0.0))
cy_skt.VerticalConstraint(entity=cy_geo[3], addUndoState=False)
cy_skt.PerpendicularConstraint(entity1=cy_geo[2], entity2=cy_geo[3], addUndoState=False)
cy_skt.Line(point1=(0.0, 0.0), point2=(cy_radius, 0.0))
cy_skt.HorizontalConstraint(entity=cy_geo[4], addUndoState=False)
cy_skt.PerpendicularConstraint(entity1=cy_geo[3], entity2=cy_geo[4], addUndoState=False)

ft_prt = ft_md.Part(name='Part-1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
ft_prt.BaseSolidExtrude(sketch=cy_skt, depth=cy_len)

del cy_skt

# Draw two arcs at top face
prt_faces = ft_prt.faces
prt_edges = ft_prt.edges
trans = ft_prt.MakeSketchTransform(sketchPlane=prt_faces[3], sketchUpEdge=prt_edges[7], \
    sketchPlaneSide=SIDE1, origin=(0.0, 0.0, cy_len))
cut_arcs_skt = ft_md.ConstrainedSketch(name='__profile__', transform=trans, sheetSize=100.0)
cl_geo = cut_arcs_skt.geometry
cl_vrt = cut_arcs_skt.vertices
ft_prt.projectReferencesOntoSketch(sketch=cut_arcs_skt, filter=COPLANAR_EDGES)
cut_arcs_skt.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, -cy_dense_radius), point2=(cy_dense_radius, 0.0), direction=COUNTERCLOCKWISE)
cut_arcs_skt.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, -cy_coarse_radius), point2=(cy_coarse_radius, 0.0), direction=COUNTERCLOCKWISE)

partition_faces = prt_faces.findAt(((cy_radius*0.5, cy_radius*0.5, cy_len), ), )
ft_prt.PartitionFaceBySketch(sketchUpEdge=prt_edges[7], faces=partition_faces, sketch=cut_arcs_skt)

del cut_arcs_skt

prt_cells = ft_prt.cells
# circular partition 1
part_sweep_dir1 = prt_edges.findAt(((0.0, 0.0, cy_depth*0.5), ), )[0]
part_sweep_edg1 = prt_edges.findAt(((cy_coarse_radius*half_sqrt_2, cy_coarse_radius*half_sqrt_2, cy_len), ), )
part_cell = prt_cells.findAt(((cy_radius*0.5, cy_radius*0.5, cy_depth*0.5), ), )
ft_prt.PartitionCellByExtrudeEdge(line=part_sweep_dir1, cells=part_cell, edges=part_sweep_edg1, sense=REVERSE)
# circular partition 2
part_sweep_dir2 = prt_edges.findAt(((0.0, 0.0, cy_depth*0.5), ), )[0]
part_sweep_edg2 = prt_edges.findAt(((cy_dense_radius*half_sqrt_2, cy_dense_radius*half_sqrt_2, cy_len), ), )
part_cell = prt_cells.findAt(((cy_coarse_radius*0.5, cy_coarse_radius*0.5, cy_depth*0.5), ), )
ft_prt.PartitionCellByExtrudeEdge(line=part_sweep_dir2, cells=part_cell, edges=part_sweep_edg2, sense=REVERSE)

# z partitions
ft_prt.DatumAxisByPrincipalAxis(principalAxis = ZAXIS)
ft_prt.DatumPointByCoordinate(coords = (0.0, 0.0, cy_depth))
ft_prt.DatumPointByCoordinate(coords = (0.0, 0.0, cy_depth - cy_dense_depth))
ft_prt.DatumPointByCoordinate(coords = (0.0, 0.0, cy_depth - cy_coarse_depth))
datum_keys = ft_prt.datums.keys()
z_axis = ft_prt.datums[datum_keys[0]] # z axis
zpt1 = ft_prt.datums[datum_keys[1]] # z point 1
zpt2 = ft_prt.datums[datum_keys[2]] # z point 2
zpt3 = ft_prt.datums[datum_keys[3]] # z point 3

# z cut 1
cl_tmp = prt_cells.getByBoundingBox(xMin = -cy_radius*0.01, xMax = cy_radius*1.01, \
    yMin = -cy_radius*0.01, yMax = cy_radius*1.01, zMin = -cy_len*0.01, zMax = cy_len*1.01)
ft_prt.PartitionCellByPlanePointNormal(point = zpt1, normal = z_axis, cells = cl_tmp)
# z cut 2
cl_tmp = prt_cells.getByBoundingBox(xMin = -cy_radius*0.01, xMax = cy_radius*1.01, \
    yMin = -cy_radius*0.01, yMax = cy_radius*1.01, zMin = -cy_len*0.01, zMax = cy_len*1.01)
ft_prt.PartitionCellByPlanePointNormal(point = zpt2, normal = z_axis, cells = cl_tmp)
# z cut 3
cl_tmp = prt_cells.getByBoundingBox(xMin = -cy_radius*0.01, xMax = cy_radius*1.01, \
    yMin = -cy_radius*0.01, yMax = cy_radius*1.01, zMin = -cy_len*0.01, zMax = cy_len*1.01)
ft_prt.PartitionCellByPlanePointNormal(point = zpt3, normal = z_axis, cells = cl_tmp)

# mesh the part
# seed globally
ft_prt.seedPart(size = coarse_elem_size, deviationFactor = 0.1, minSizeFactor = 0.1)
all_cells = prt_cells.getByBoundingBox(xMin = -cy_radius*0.01, xMax = cy_radius*1.01, \
    yMin = -cy_radius*0.01, yMax = cy_radius*1.01, zMin = -cy_len*0.01, zMax = cy_len*1.01)
ft_prt.setMeshControls(regions = all_cells, elemShape = TET, technique = FREE, allowMapped = False)
elem_type3 = mesh.ElemType(elemCode = C3D4, elemLibrary = STANDARD, secondOrderAccuracy = OFF, distortionControl = DEFAULT)
ac_set = ft_prt.Set(name = 'AllCells', cells = all_cells)
ft_prt.setElementType(regions = ac_set, elemTypes = (elem_type3, ))
# seed edge
seed_edg = prt_edges.findAt(
    # arcs
    ((cy_dense_radius*half_sqrt_2, cy_dense_radius*half_sqrt_2, cy_len), ),
    ((cy_dense_radius*half_sqrt_2, cy_dense_radius*half_sqrt_2, cy_depth), ),
    ((cy_dense_radius*half_sqrt_2, cy_dense_radius*half_sqrt_2, cy_depth - cy_dense_depth), ),
    # horizontal lines
    ((cy_dense_radius*0.5, 0.0, cy_len), ),
    ((cy_dense_radius*0.5, 0.0, cy_depth), ),
    ((cy_dense_radius*0.5, 0.0, cy_depth - cy_dense_depth), ),
    ((0.0, cy_dense_radius*0.5, cy_len), ),
    ((0.0, cy_dense_radius*0.5, cy_depth), ),
    ((0.0, cy_dense_radius*0.5, cy_depth - cy_dense_depth), ),
    # vertical lines
    ((cy_dense_radius, 0.0, cy_depth + 0.5*cy_top), ),
    ((0.0, cy_dense_radius, cy_depth + 0.5*cy_top), ),
    ((0.0, 0.0, cy_depth + 0.5*cy_top), ),
    ((cy_dense_radius, 0.0, cy_depth - 0.5*cy_dense_depth), ),
    ((0.0, cy_dense_radius, cy_depth - 0.5*cy_dense_depth), ),
    ((0.0, 0.0, cy_depth - 0.5*cy_dense_depth), ),
    )
seed_edge_set = ft_prt.Set(name = 'SeedEdges', edges = seed_edg)
ft_prt.seedEdgeBySize(edges = seed_edg, size = dense_elem_size, \
    deviationFactor = 0.1, minSizeFactor = 0.1, constraint = FIXED)
# generate mesh
ft_prt.generateMesh()

assembly = ft_md.rootAssembly
assembly.DatumCsysByDefault(CARTESIAN)
assembly.Instance(name = 'Inst-1', part = ft_prt, dependent = ON)

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
