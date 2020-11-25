#////////////////////////////////////////////////////////////////////////////
# Explicit FE Model - Sinusoidal Input for modelling IBUS 
# Author: Xavie REGAL 
# Adapted from python code developed by LLoyd Fletcher/Alex Guigue
# Date: 28/7/2017
#////////////////////////////////////////////////////////////////////////////

#----------------------------------------------------------------------------
# Import preamble to access abaqus commands
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup() 
import visualization 
import displayGroupOdbToolset as dgo 
import regionToolset
import os 
import section
import step
import job
import material
import sys
import time
import math

#del mymodel

#----------------------------------------------------------------------------
# Define file paths and working directory
directory='C:/Abaqus_DATA/'
extension='Beam_isotrope/'
jobFileName='strainrate1'
SubroutinesPath = 'C:\\Users\\radg0\\Documents\\abaqus_data\\VUMAT\\VUMAT_isotrope_strainrate_dependent_VUAMP.for'

# If the directory does not exist, create it
if not os.path.exists(directory+extension):
    os.makedirs(directory+extension)#Make the directory if need be

# Set the working directory
os.chdir(directory+extension)

#----------------------------------------------------------------------------
# Model Parameters

# SIMULATION FLAGS
submitJob = 1
extractData = 1
extractDispEnerOnly = 0

# SIMULATION PARAMETERS
blankFrames = 5
frameRate = 5e5
outputStep = 2e-6
specOutputFrames = 128
solveStep = 0.01e-6
outputFreq = int(outputStep/solveStep)
betaDamping = 10e-7
bVisLDamp = 0.0
bVisQDamp = 0.0

# LOADING FREQUENCY DEFINITION
freq0     = 20000.0
omega0  = 2*math.pi*freq0
disp_ampl = 3e-5 

# SPECIMEN
specElemSize = 0.50e-3
specElemSize = 0.50e-3
specLength = 41e-3
specHeight = 12e-3
specThick = 2e-3
# PVC
specRho = 1290.0
specExx1 = 2.0e9
specExx2 = 0.5e9
specNuxy = 0.41
Cv = 1450 
specWaveSpeed = math.sqrt(specExx1/specRho)
beta = 0 
delta1 = 0.1
delta2 = 0.3
lambda1 = 0.2


# Total Simulation Time Calculation
guideOutputFrames = 0	# only used for projectile impact model
totalOutputSteps = int(specOutputFrames+blankFrames)
totalSimTime = totalOutputSteps*outputStep
totalSimTime = 0.002; 
#----------------------------------------------------------------------------
# Model and viewport definitions
mymodel=mdb.Model(name=jobFileName)

myviewport=session.Viewport(name='IBUS',origin=(0.0, 0.0))
myviewport.makeCurrent()
myviewport.maximize()

#----------------------------------------------------------------------------
# Part Definitions 
import part
myviewport.partDisplay.geometryOptions.setValues(
        referenceRepresentation=ON)

# SPECIMEN
s1 = mymodel.ConstrainedSketch(name='profile', 
        sheetSize=200.0)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=STANDALONE)
s1.rectangle(point1=(0.0, 0.0), point2=(specLength, specHeight))

part_specimen = mymodel.Part(name='Specimen', dimensionality=TWO_D_PLANAR, 
        type=DEFORMABLE_BODY)
part_specimen.BaseShell(sketch=s1)
s1.unsetPrimaryObject()
myviewport.setValues(displayedObject=part_specimen)
del mymodel.sketches['profile']

#----------------------------------------------------------------------------
# Assembly Definition
import assembly
myassembly = mymodel.rootAssembly

# Create an instance of the specimen
specimen_inst=myassembly.Instance(name='SpecimenInstance', part=part_specimen, dependent=ON)

#----------------------------------------------------------------------------
# Material Definition

specMat=mymodel.Material(name='specMat')
#General 
specMat.Density(table=((specRho, ), ))
#Mechanical 
specMat.Depvar(n=3)
specMat.UserMaterial(mechanicalConstants=(specExx1,specNuxy, delta1))

#Thermal 
specMat.SpecificHeat( table=((Cv, ), ))
specMat.InelasticHeatFraction(fraction = (0.9))
specMat.Conductivity( table=((lambda1, ), ))
	
part_specimen.DatumCsysByThreePoints(name='Datum csys-1', coordSysType=CARTESIAN, origin=(
	0.0, 0.0, 0.0), line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0))

#----------------------------------------------------------------------------
# Section Definition			
myassembly.regenerate()    
section_specimen=mymodel.HomogeneousSolidSection(name='Section-Specimen', 
	material='specMat', thickness=specThick)
	
# Regions
region_specimen=(part_specimen.faces,)

# Section Assignment	
part_specimen.SectionAssignment(region=region_specimen, sectionName='Section-Specimen', offset=0.0, 
	offsetType=MIDDLE_SURFACE, offsetField='', 
	thicknessAssignment=FROM_SECTION)

#----------------------------------------------------------------------------
# Loading Definition
import load
	
	
mymodel.PeriodicAmplitude(name='Sinwave', timeSpan=STEP, 
    frequency=omega0, start=0.0, a_0=0.0, data=((0.0, 1.0), ))	
	
mymodel.UserAmplitude(name='SinwaveU', numVariables=0)	
	
#----------------------------------------------------------------------------
# Meshing
# Seed parts with element edge size
part_specimen.seedPart(size=specElemSize, deviationFactor=0.1, minSizeFactor=0.1)
	
# Element type definition
elementType=mesh.ElemType(elemCode=CPS4RT, elemLibrary=EXPLICIT,secondOrderAccuracy=OFF, 
hourglassControl=DEFAULT,  distortionControl=DEFAULT)

# Assign element type to parts
part_specimen.setElementType(regions=(region_specimen,),elemTypes=(elementType,))

# Generate the mesh
part_specimen.generateMesh()
myviewport.partDisplay.setValues(mesh=ON)
myviewport.assemblyDisplay.setValues(mesh=ON)
myviewport.assemblyDisplay.meshOptions.setValues(meshTechnique=ON)

# Define element and node sets for output
specimen_el=part_specimen.elements
specimen_nd=part_specimen.nodes

tol = specElemSize/2
specimen_elements=specimen_el.getByBoundingBox(0,-specHeight/2-tol,-tol,specLength+tol,specHeight/2+tol,tol)
part_specimen.Set(name='Set_elements_specimen', elements=specimen_elements)

specimen_nodes=specimen_nd.getByBoundingBox(0,-specHeight/2-tol,-tol,specLength+tol,specHeight/2+tol,tol)
part_specimen.Set(name='Set_node_specimen', nodes=specimen_nodes)

# Show the mesh
myviewport.assemblyDisplay.setValues(mesh=ON)
myviewport.assemblyDisplay.meshOptions.setValues(meshTechnique=ON)

#----------------------------------------------------------------------------
# Explicit Dynamics - Load Step Definition
name_step='Sonotrode'
#ExplicitDynamicsStep TempDisplacementDynamicsStep
mymodel.TempDisplacementDynamicsStep(name=name_step, previous='Initial', 
	timePeriod=totalSimTime, description='Sonotrode sollicitation',
	timeIncrementationMethod=FIXED_USER_DEFINED_INC, userDefinedInc=solveStep)
#mymodel.steps['Sonotrode'].setValues(
#    timeIncrementationMethod=AUTOMATIC_GLOBAL, scaleFactor=1.0, 
#    maxIncrement=solveStep)
step_impact=mymodel.steps[name_step]
step_impact.setValues(nlgeom=OFF)

# Specify bulk viscosity
step_impact.setValues(linearBulkViscosity=bVisLDamp, 
	quadBulkViscosity=bVisQDamp)

myviewport.assemblyDisplay.setValues(loads=ON, bcs=ON, 
	predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)

#----------------------------------------------------------------------------
# Explicit Dynamics - Displacement Pulse

# Edges of each part
edge_specimen = specimen_inst.edges 

# Specimen Contact Face
X1 = specLength-tol
Y1 = -tol
Z1 = -tol
X2 = specLength+tol
Y2 = specHeight+tol
Z2 = tol
ed_specimen=edge_specimen.getByBoundingBox(X1,Y1,Z1,X2,Y2,Z2)

# Create the surfaces
surf_specimen=myassembly.Surface(name='Surf-Specimen',side1Edges=ed_specimen)
region1 = regionToolset.Region(edges=ed_specimen)

# Define the imposed displacement
mymodel.DisplacementBC(name='Imposed_disp', 
    createStepName='Sonotrode', region=region1, u1=disp_ampl, u2=UNSET, ur3=UNSET, 
    amplitude='SinwaveU', fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)	
	
# Define the initial temperature 
#region_specimen=(part_specimen.faces,)

face_temp = specimen_inst.faces ;
region2 = regionToolset.Region(faces=face_temp)

mymodel.Temperature(name='Initial-Temp', createStepName='Initial', region=region2, 
   distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(20.0, ))

a = mymodel.rootAssembly
#region3 = a.surfaces['Surf-Specimen']
#mymodel.SurfaceHeatFlux(name='Load-1', createStepName='Sonotrode', region=region3, magnitude=20.0)
	
	
myviewport.assemblyDisplay.setValues(interactions=OFF, 
	constraints=OFF, connectors=OFF, engineeringFeatures=OFF)

#----------------------------------------------------------------------------
# Explicit Dynamics - Job and Output Definition

# Field Outputs
if extractDispEnerOnly:
	mymodel.fieldOutputRequests['F-Output-1'].setValues(variables=('UT','S','COORD'))
else:
	mymodel.fieldOutputRequests['F-Output-1'].setValues(variables=(
		'S', 'E', 'ER', 'UT', 'AT','COORD','TEMP','ENER'))
mymodel.fieldOutputRequests['F-Output-1'].setValues(timeInterval=outputStep,timeMarks=ON)

# History Outputs - Energy
mymodel.HistoryOutputRequest(name='H-Output-1', 
    createStepName='Sonotrode', variables=PRESELECT)
mymodel.historyOutputRequests['H-Output-1'].setValues(timeInterval=outputStep)

# Define the Job
Sonotrode=mdb.Job(name=jobFileName, model=jobFileName, description='', type=ANALYSIS, 
	atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
	memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE, 
	nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF, 
	contactPrint=OFF, historyPrint=OFF,
	userSubroutine=SubroutinesPath, 
	scratch='', 
	parallelizationMethodExplicit=DOMAIN, numDomains=1, 
	activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=1)
print 'Job created. \n'

	
if not os.path.exists(directory+extension+jobFileName):
	os.makedirs(directory+extension+jobFileName)
# Set as working directory
os.chdir(directory+extension+jobFileName)

# Submit the job
Sonotrode.submit(consistencyChecking=OFF)
print'Job '+jobFileName+' submitted'+'\n'+''

# Wait for job completion
print'Waiting for '+jobFileName+' to finish'+'\n'+'   '
Sonotrode.waitForCompletion()

print jobFileName+' FINISHED'+'\n'+'   '

#////////////////////////////////////////////////////////////////////////////
# Data Extraction
#////////////////////////////////////////////////////////////////////////////


from abaqusConstants import *
from odbAccess import *
import visualization
# Point at the directory where the odb file is
os.chdir(directory+extension+jobFileName)
	
# Open the database
o1 = session.openOdb(name = jobFileName, path = jobFileName + '.odb', readOnly=True) # Open the output database :o1
myviewport.setValues(displayedObject=o1) # Put the output in the viewport

if extractData == 1:



	# Create folders for extracted values and later for images
	extractFolder='/Extracted_Data/'
	if not os.path.exists(directory+extension+jobFileName+extractFolder):
		os.makedirs(directory+extension+jobFileName+extractFolder)	
	
	# Open the odb file
	odb=session.odbs[jobFileName +'.odb']

	nsteps=len(odb.steps['Sonotrode'].frames)-1

	for i in range(nsteps-130,nsteps):
	#for i in range(1,nsteps):
		fileName_E  = directory+extension+jobFileName+extractFolder+'E_'+str(i)+'.rpt'
		fileName_ER = directory+extension+jobFileName+extractFolder+'ER_'+str(i)+'.rpt'
		fileName_S  = directory+extension+jobFileName+extractFolder+'S_'+str(i)+'.rpt'
		fileName_U  = directory+extension+jobFileName+extractFolder+'U_'+str(i)+'.rpt'   
		fileName_A  = directory+extension+jobFileName+extractFolder+'A_'+str(i)+'.rpt'  
		fileName_T  = directory+extension+jobFileName+extractFolder+'T_'+str(i)+'.rpt'  
		fileName_DE = directory+extension+jobFileName+extractFolder+'DE_'+str(i)+'.rpt' 
		session.writeFieldReport(fileName=fileName_E,  append=OFF, sortItem='Element Label', odb=odb, step=0, frame=i, outputPosition=NODAL, variable=(('E', INTEGRATION_POINT, ((COMPONENT, 'E11'), (COMPONENT, 'E22'), (COMPONENT, 'E12'), )), ))
		session.writeFieldReport(fileName=fileName_ER, append=OFF, sortItem='Element Label', odb=odb, step=0, frame=i, outputPosition=NODAL, variable=(('ER', INTEGRATION_POINT, ((COMPONENT, 'ER11'), (COMPONENT, 'ER22'), (COMPONENT, 'ER12'), )), ))
		session.writeFieldReport(fileName=fileName_S,  append=OFF, sortItem='Element Label', odb=odb, step=0, frame=i, outputPosition=NODAL,  variable=(('S', INTEGRATION_POINT, ((COMPONENT, 'S11'), (COMPONENT, 'S22'), (COMPONENT, 'S12'), )), ))
		session.writeFieldReport(fileName=fileName_U,  append=OFF, sortItem='Node Label', odb=odb, step=0, frame=i, outputPosition=NODAL, variable=(('COORD', NODAL, ((COMPONENT, 'COOR1'), (COMPONENT, 'COOR2'), )),('UT', NODAL, ((COMPONENT, 'U1'), (COMPONENT, 'U2'), )), ))
		session.writeFieldReport(fileName=fileName_A,  append=OFF, sortItem='Node Label', odb=odb, step=0, frame=i, outputPosition=NODAL, variable=(('AT', NODAL, ((COMPONENT, 'A1'), (COMPONENT, 'A2'), )), ))
		session.writeFieldReport(fileName=fileName_T,  append=OFF, sortItem='Node Label', odb=odb, step=0, frame=i, outputPosition=NODAL, variable=(('TEMP',INTEGRATION_POINT, ), ))
		session.writeFieldReport(fileName=fileName_DE, append=OFF, sortItem='Node Label', odb=odb, step=0, frame=i, outputPosition=NODAL, variable=(('PENER',INTEGRATION_POINT, ), ))


	os.chdir(directory+extension+jobFileName)

	# Open the database
	o1 = session.openOdb(name = jobFileName, path = jobFileName + '.odb', readOnly=True) # Open the output database :o1
	myviewport.setValues(displayedObject=o1) # Put the output in the viewport

	# Create folders for extracted values and later for images
	extractFolder='/Extracted_Data/'
	if not os.path.exists(directory+extension+jobFileName+extractFolder):
		os.makedirs(directory+extension+jobFileName+extractFolder)
	
	# Open the odb file
	odb=session.odbs[jobFileName +'.odb']
 
	
	f = open(directory+extension+jobFileName+extractFolder+'nTimeSteps.txt', 'w')
	f.write(str(nsteps));
	f.close()

	f = open(directory+extension+jobFileName+extractFolder+'param.txt', 'w')
	f.write(str(nsteps)+'\n');
	f.write(str(frameRate)+'\n');
	f.write(str(specRho)+'\n');
	f.write(str(specExx1)+'\n');
	f.write(str(specNuxy)+'\n');
	f.write(str(beta)+'\n');
	f.write(str(freq0)+'\n');
	f.write(str(disp_ampl)+'\n');
	f.write(str(specLength)+'\n');
	f.write(str(specHeight)+'\n');
	f.write(str(specThick)+'\n');
	f.write(str(delta1)+'\n');
	f.close()
	