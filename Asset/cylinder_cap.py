from abaqus import *
from abaqusConstants import *
from caeModules import *

radius = 0.15
height = 0.05

md = mdb.models['Model-1']

cy_skt = md.ConstrainedSketch(name='__profile__', sheetSize=10.0)
cy_skt.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, radius))

cy_prt = md.Part(name = 'Part-1', dimensionality = THREE_D, type = DEFORMABLE_BODY)
cy_prt.BaseSolidExtrude(sketch=cy_skt, depth=height)
del cy_skt

cy_prt.seedPart(size = radius*0.2, deviationFactor = 0.1, minSizeFactor = 0.1)
cl_tmp = cy_prt.cells.getByBoundingBox( \
    xMin = -radius*1.01, xMax = radius*1.01, \
    yMin = -radius*1.01, yMax = radius*1.01, \
    zMin = -height*0.01, zMax = height*1.01)
cy_prt.setMeshControls(regions = cl_tmp, elemShape = TET, technique = FREE)
elemType1 = mesh.ElemType(elemCode = C3D8R, elemLibrary = STANDARD)
elemType2 = mesh.ElemType(elemCode = C3D6, elemLibrary = STANDARD)
elemType3 = mesh.ElemType(elemCode = C3D4, elemLibrary = STANDARD, 
    secondOrderAccuracy = OFF, distortionControl = DEFAULT)
cy_prt.setElementType(regions = (cl_tmp, ), elemTypes=(elemType3, ))

cy_prt.generateMesh()

assembly = md.rootAssembly
assembly.DatumCsysByDefault(CARTESIAN)
assembly.Instance(name = 'Inst-1', part = cy_prt, dependent = ON)

job = mdb.Job(name='cylinder_cap', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, 
    numDomains=1, activateLoadBalancing=False, multiprocessingMode=DEFAULT, 
    numCpus=1)
job.writeInput(consistencyChecking=OFF)

mdb.saveAs(pathName='./cylinder')