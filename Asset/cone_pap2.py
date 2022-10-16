import math
from abaqus import *
from abaqusConstants import *
from caeModules import *

cae_name = "cone_pap2"
cone_radius = 0.2
mh_size_ratio = 0.25

md = mdb.models['Model-1']
cone_prt = md.Part(name='Part-1', dimensionality=THREE_D, type=DEFORMABLE_BODY)

cone_skt = md.ConstrainedSketch(name='__profile__', sheetSize=10.0)
cone_skt.ConstructionLine(point1=(0.0, -10.0), point2=(0.0, 10.0))
cone_skt.Line(point1=(0.0, 0.0), point2=(cone_radius, 0.0))
cone_skt.Line(point1=(cone_radius, 0.0), point2=(0.0, cone_radius))
cone_skt.Line(point1=(0.0, cone_radius), point2=(0.0, 0.0))
#cone_skt.setPrimaryObject(option=STANDALONE)
cone_prt.BaseSolidRevolve(sketch=cone_skt, angle=360.0, flipRevolveDirection=OFF)
del cone_skt

cone_prt.seedPart(size = cone_radius*mh_size_ratio, deviationFactor = 0.1, minSizeFactor = 0.1)
cone_prt_cells = cone_prt.cells
all_cells = cone_prt_cells.getByBoundingBox( \
    xMin = -cone_radius*1.01, xMax = cone_radius*1.01, \
    yMin = -cone_radius*0.01, yMax = cone_radius*1.01, \
    zMin = -cone_radius*1.01, zMax = cone_radius*1.01)
cone_prt.setMeshControls(regions = all_cells, elemShape = TET, technique = FREE, allowMapped = False)
ac_set = cone_prt.Set(name = 'AllCells', cells = all_cells)
elem_type3 = mesh.ElemType(elemCode = C3D4, elemLibrary = STANDARD, secondOrderAccuracy = OFF, distortionControl = DEFAULT)
cone_prt.setElementType(regions = ac_set, elemTypes = (elem_type3, ))
cone_prt.generateMesh()

assembly = md.rootAssembly
assembly.DatumCsysByDefault(CARTESIAN)
assembly.Instance(name = 'Inst-1', part = cone_prt, dependent = ON)

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
