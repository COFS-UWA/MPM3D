import math
from abaqus import *
from abaqusConstants import *
from caeModules import *

# input parameters
footing_radius = 1.5
cae_name = 'spudcan_soil_half_cylinder_3_5R'

cy_radius = 3.5 * footing_radius
cy_top = 0.5 * footing_radius
cy_depth = 4.0 * footing_radius

elem_size = 0.25 * footing_radius

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

cy_prt_cells = cy_prt.cells

# z partitions
cy_prt.DatumAxisByPrincipalAxis(principalAxis = ZAXIS)
cy_prt.DatumPointByCoordinate(coords = (0.0, 0.0, cy_depth))
datum_keys = cy_prt.datums.keys()
z_axis = cy_prt.datums[datum_keys[0]] # z axis
zpt1 = cy_prt.datums[datum_keys[1]] # z point 1
# z cut 1
cl_tmp = cy_prt_cells.getByBoundingBox(xMin = -cy_radius*1.01, xMax = cy_radius*1.01, \
    yMin = -cy_radius*0.01, yMax = cy_radius*1.01, zMin = -cy_len*0.01, zMax = cy_len*1.01)
cy_prt.PartitionCellByPlanePointNormal(point = zpt1, normal = z_axis, cells = cl_tmp)

# mesh the part
cy_prt.seedPart(size = elem_size, deviationFactor = 0.1, minSizeFactor = 0.1)
all_cells = cy_prt_cells.getByBoundingBox(xMin = -cy_radius*1.01, xMax = cy_radius*1.01, \
    yMin = -cy_radius*0.01, yMax = cy_radius*1.01, zMin = -cy_len*0.01, zMax = cy_len*1.01)
cy_prt.setMeshControls(regions = all_cells, elemShape = TET, technique = FREE, allowMapped = False)
elem_type3 = mesh.ElemType(elemCode = C3D4, elemLibrary = STANDARD, secondOrderAccuracy = OFF, distortionControl = DEFAULT)
ac_set = cy_prt.Set(name = 'AllCells', cells = all_cells)
cy_prt.setElementType(regions = ac_set, elemTypes = (elem_type3, ))
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
