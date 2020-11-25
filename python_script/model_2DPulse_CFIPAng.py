#////////////////////////////////////////////////////////////////////////////
# Explicit FE Model - Pressure Pulse on CFRP Plate
# Author: Lloyd Fletcher
# Adapted from python code developed by Alex Guigue
# Date: 22/1/2019
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
directory='C:/Abaqus_WorkingDirectory/'
extension='Model_2D_CFRP_OA/'
jobFileName='2DPulse_CFRP_Ang'

# If the directory does not exist, createn it
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
outputGlobalCoords = 1

# SIMULATION PARAMETERS
blankFrames = 3
frameRate = 2e6
outputStep = 1/frameRate
specOutputFrames = 128
solveStep = 0.050e-6
outputFreq = int(outputStep/solveStep)
betaDamping = 1e-8
bVisLDamp = 0.0
bVisQDamp = 0.0

# SPECIMEN
specElemSize = 0.5e-3
specLength = 70e-3
specHeight = 44e-3
specThick = 1e-3
# CFRP - ARL Plate
specAngle = 90
specRho = 1600
specE11 = 124e9
specE22 = 8.3e9
specNu12= 0.3
specG12 = 3.7e9

# LOADING PULSE DEFINITION
pulseStart = blankFrames*outputStep
pulseTriangle = 1
pulseDuration = 18e-6
pulseRise = 9e-6	
pulseFall = 9e-6
pulsePressure = 300e6	# Pa

# Total Simulation Time Calculation
guideOutputFrames = 0
totalOutputSteps = int(specOutputFrames)
totalSimTime = totalOutputSteps*outputStep

#----------------------------------------------------------------------------
# Model and viewport definitions
mymodel=mdb.Model(name=jobFileName)

myviewport=session.Viewport(name='Pulse Model 2D',origin=(0.0, 0.0))
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
s1.rectangle(point1=(0.0, -specHeight/2), point2=(specLength, specHeight/2))

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
specMat.Density(table=((specRho, ), ))
specMat.Elastic(type=LAMINA, table=((specE11, specE22, specNu12, specG12, 0.0, 0.0,), ))
specMat.Damping(beta=betaDamping)

part_specimen.DatumCsysByThreePoints(name='Datum csys-1', coordSysType=CARTESIAN, origin=(
	0.0, 0.0, 0.0), line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0))

#----------------------------------------------------------------------------
# Section Definition			
myassembly.regenerate()    
section_specimen=mymodel.HomogeneousSolidSection(name='Section-Specimen', 
	material='specMat', thickness=specThick)
	
# Regions
region_specimen=(part_specimen.faces,)

# Co-ords for orthotropic material
orientation = part_specimen.datums[2]

part_specimen.MaterialOrientation(region=region_specimen, 
	orientationType=SYSTEM, axis=AXIS_3, localCsys=orientation, fieldName='', 
	additionalRotationType=ROTATION_ANGLE, angle=specAngle, 
	additionalRotationField='', stackDirection=STACK_3) 

# Section Assignment	
part_specimen.SectionAssignment(region=region_specimen, sectionName='Section-Specimen', offset=0.0, 
	offsetType=MIDDLE_SURFACE, offsetField='', 
	thicknessAssignment=FROM_SECTION)

#----------------------------------------------------------------------------
# Loading Definition
import load
if pulseTriangle:
	t1 = pulseStart
	t2 = t1+pulseRise
	t3 = t1+pulseDuration
	mymodel.TabularAmplitude(name='Amp-1', timeSpan=STEP,smooth=SOLVER_DEFAULT, 
		data=((0.0, 0.0), (t1, 0.0), (t2, 1.0), (t3, 0.0)))
else:
	t1 = pulseStart
	t2 = t1+pulseRise
	t3 = t1+pulseDuration-pulseFall
	t4 = t1+pulseDuration
	mymodel.TabularAmplitude(name='Amp-1', timeSpan=STEP,smooth=SOLVER_DEFAULT, 
		data=((0.0, 0.0), (t1, 0.0), (t2, 1.0), (t3, 1.0), (t4, 0.0)))
	
#----------------------------------------------------------------------------
# Meshing
# Seed parts with element edge size
part_specimen.seedPart(size=specElemSize, deviationFactor=0.1, minSizeFactor=0.1)
	
# Element type definition
elementType=mesh.ElemType(elemCode=CPS4R, elemLibrary=EXPLICIT,secondOrderAccuracy=OFF, 
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
name_step='Impact'
mymodel.ExplicitDynamicsStep(name=name_step, previous='Initial', 
	timePeriod=totalSimTime, description='Inertial Impact',
	timeIncrementationMethod=FIXED_USER_DEFINED_INC, userDefinedInc=solveStep)
mymodel.steps['Impact'].setValues(
    timeIncrementationMethod=AUTOMATIC_GLOBAL, scaleFactor=1.0, 
    maxIncrement=solveStep)
step_impact=mymodel.steps[name_step]
step_impact.setValues(nlgeom=OFF)

# Specify bulk viscosity
step_impact.setValues(linearBulkViscosity=bVisLDamp, 
	quadBulkViscosity=bVisQDamp)

myviewport.assemblyDisplay.setValues(loads=ON, bcs=ON, 
	predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)

#----------------------------------------------------------------------------
# Explicit Dynamics - Pressure Pulse

# Edges of each part
edge_specimen = specimen_inst.edges 

# Specimen Contact Face
X1 = specLength-tol
Y1 = -specHeight/2
Z1 = -tol
X2 = specLength+tol
Y2 = specHeight/2
Z2 = tol
ed_specimen=edge_specimen.getByBoundingBox(X1,Y1,Z1,X2,Y2,Z2)

# Create the surfaces
surf_specimen=myassembly.Surface(name='Surf-Specimen',side1Edges=ed_specimen)

# Define the pressure pulse on the specimen edge
mymodel.Pressure(name='Load-1', createStepName='Impact', 
    region=surf_specimen, distributionType=UNIFORM, field='', magnitude=pulsePressure, 
    amplitude='Amp-1')

myviewport.assemblyDisplay.setValues(interactions=OFF, 
	constraints=OFF, connectors=OFF, engineeringFeatures=OFF)

#----------------------------------------------------------------------------
# Explicit Dynamics - Job and Output Definition

# Field Outputs
if extractDispEnerOnly:
	mymodel.fieldOutputRequests['F-Output-1'].setValues(variables=('UT','S','COORD'))
else:
	mymodel.fieldOutputRequests['F-Output-1'].setValues(variables=(
		'S', 'E', 'UT', 'AT','COORD'))
mymodel.fieldOutputRequests['F-Output-1'].setValues(timeInterval=outputStep,timeMarks=ON)

# History Outputs - Energy
mymodel.HistoryOutputRequest(name='H-Output-1', 
    createStepName='Impact', variables=('ALLAE', 'ALLDC', 'ALLDMD', 
    'ALLIE', 'ALLKE', 'ALLSE', 'ALLVD', 'ALLWK', 'ALLCW', 'ALLPW', 'ETOTAL'))
mymodel.historyOutputRequests['H-Output-1'].setValues(timeInterval=outputStep)

# Define the Job
Impact=mdb.Job(name=jobFileName, model=jobFileName, description='', type=ANALYSIS, 
	atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
	memoryUnits=PERCENTAGE, explicitPrecision=SINGLE, 
	nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
	contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
	parallelizationMethodExplicit=DOMAIN, numDomains=1, 
	activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=1)
print 'Job created. \n'

if not os.path.exists(directory+extension+jobFileName):
	os.makedirs(directory+extension+jobFileName)
# Set as working directory
os.chdir(directory+extension+jobFileName)

# Submit the job
Impact.submit(consistencyChecking=OFF)
print'Job '+jobFileName+' submitted'+'\n'+''

# Wait for job completion
print'Waiting for '+jobFileName+' to finish'+'\n'+'   '
Impact.waitForCompletion()

print jobFileName+' FINISHED'+'\n'+'   '

#////////////////////////////////////////////////////////////////////////////
# Data Extraction
#////////////////////////////////////////////////////////////////////////////
# Point at the directory where the odb file is
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

# Create a global co-ord system
odb.rootAssembly.DatumCsysByThreePoints(name='CSYS-1', 
    coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), 
    point2=(0.0, 1.0, 0.0))

if outputGlobalCoords:	
	dtm = odb.rootAssembly.datumCsyses['CSYS-1']
	session.viewports['Pulse Model 2D'].odbDisplay.basicOptions.setValues(
		transformationType=USER_SPECIFIED, datumCsys=dtm)
	
# Check to see if the simulation was completed
nsteps = len(odb.steps['Impact'].frames) #Count the number of frames
print 'For '+jobFileName+ ' the number of steps is '+str(nsteps)+'\n'+'   '	

#--------------------------------------------------------------------
# Write the simulation configuration to file
configFilePath = directory+extension+jobFileName+extractFolder+jobFileName+'_SimConfig.txt'
simConfigFile = open(configFilePath,'w')

simConfigFile.write('SIMULATION PARAMETERS\n')
simConfigFile.write('{0:10} {1:10} {2:10}'.format('Ttot','Tsol','Tout')+'\n')
simConfigFile.write('{0:<10.3e} {1:<10.3e} {2:<10.3e}'.format(
	totalSimTime,solveStep,outputStep)+'\n')
simConfigFile.write('{0:10} {1:10} {2:10} {3:10}'.format('Ftot','Ftst','Fwvg','Fblank')+'\n')
simConfigFile.write('{0:<10.3e} {1:<10.3e} {2:<10.3e} {3:<10.3e}'.format(
	totalOutputSteps,specOutputFrames,guideOutputFrames,blankFrames)+'\n')
simConfigFile.write('{0:10} {1:10} {2:10} '.format('BetaDamp','BVisLDmp','BVisQDmp')+'\n')
simConfigFile.write('{0:<10.3e} {1:<10.3e} {2:<10.3e}'.format(
	betaDamping,bVisLDamp,bVisQDamp)+'\n')
simConfigFile.write('\n')

simConfigFile.write('SPECIMEN CONFIGURATION\n')
simConfigFile.write('{0:10} {1:10} {2:10} {3:10}'.format(
	'S_elem','S_leng','S_height','S_thick')+'\n')
simConfigFile.write('{0:<10.3e} {1:<10.3e} {2:<10.3e} {3:<10.3e}'.format(
	specElemSize,specLength,specHeight,specThick)+'\n')
simConfigFile.write('{0:10} {1:10} {2:10} {3:10} {4:10} {5:10}'.format(
	'S_Rho','S_E11','S_N12','S_E22','S_G12','S_Ang')+'\n')
simConfigFile.write('{0:<10.3e} {1:<10.3e} {2:<10.3e} {3:<10.3e} {4:<10.3e} {5:<10.3e}'.format(
	specRho,specE11,specNu12,specE22,specG12,specAngle)+'\n')
simConfigFile.write('\n')

simConfigFile.close()
#--------------------------------------------------------------------
if extractData:
	if nsteps < totalOutputSteps:
		print 'For '+jobFileName+' solveStep is too high'+'\n'+'Extraction not performed'+'\n'+'  '
		failure+=1
	else:
		# Define the file extensione for the reports
		rptExt = '.rpt'	
		
		# Create a display group that only includes the specimen nodes
		leaf_dg = dgo.LeafFromNodeSets(nodeSets='SPECIMENINSTANCE.SET_NODE_SPECIMEN',)
		spec_dg_nodes = session.DisplayGroup(leaf=leaf_dg,name='DisGroup-SpecimenNodes')
		# Create a display group that only includes the specimen elements
		leaf_dg2 = dgo.LeafFromElementSets(elementSets='SPECIMENINSTANCE.SET_ELEMENTS_SPECIMEN',)
		spec_dg_elems = session.DisplayGroup(leaf=leaf_dg2, name='DisGroup-SpecimenElements')

		# Write Abaqus reports
		session.fieldReportOptions.setValues(printMinMax=OFF, printTotal=OFF)

		# Output the node table
		rptName= directory+extension+jobFileName+extractFolder+jobFileName+'_NodeLocs'+rptExt
		if os.path.exists(rptName):
			os.remove(rptName)
		session.writeFieldReport(fileName=rptName, append=OFF, sortItem='Node Label', 
			odb=odb, step=0, frame=1, outputPosition=NODAL,displayGroup=spec_dg_nodes,
			variable=(('COORD', NODAL, ((COMPONENT, 'COOR1'), (COMPONENT, 'COOR2'), )),))
			
		# Create file names for ouput variables	
		rptName_Disp = directory+extension+jobFileName+extractFolder+jobFileName+'_Disp'
		rptName_CZDMG = directory+extension+jobFileName+extractFolder+jobFileName+'_CZDMG'
		rptName_TotEnergy = directory+extension+jobFileName+extractFolder +jobFileName+'_TotEnergy'
		if not extractDispEnerOnly:
			rptName_Accel = directory+extension+jobFileName+extractFolder +jobFileName+'_Accel'
			rptName_Stress = directory+extension+jobFileName+extractFolder +jobFileName+'_Stress'
			rptName_Strain = directory+extension+jobFileName+extractFolder +jobFileName+'_Strain'
		
		#--------------------------------------------------------------------
		# Extract the energy time history for the model
		# If the file exists clear it first so they are not appended to
		rptName = rptName_TotEnergy + rptExt
		if os.path.exists(rptName):
			os.remove(rptName)
		# Export energy time history as XY data:
		session.XYDataFromHistory(name='ENERGY-1', odb=odb, 
			outputVariableName='Artificial strain energy: ALLAE for Whole Model', )
		session.XYDataFromHistory(name='ENERGY-2', odb=odb, 
			outputVariableName='Damage dissipation energy: ALLDMD for Whole Model', )
		session.XYDataFromHistory(name='ENERGY-3', odb=odb, 
			outputVariableName='External work: ALLWK for Whole Model', )
		session.XYDataFromHistory(name='ENERGY-4', odb=odb, 
			outputVariableName='Internal energy: ALLIE for Whole Model', )
		session.XYDataFromHistory(name='ENERGY-5', odb=odb, 
			outputVariableName='Internal work by constraint penalty: ALLCW for Whole Model', )
		session.XYDataFromHistory(name='ENERGY-6', odb=odb, 
			outputVariableName='Internal work by penalty contact: ALLPW for Whole Model', )
		session.XYDataFromHistory(name='ENERGY-7', odb=odb, 
			outputVariableName='Kinetic energy: ALLKE for Whole Model', )
		session.XYDataFromHistory(name='ENERGY-8', odb=odb, 
			outputVariableName='Strain energy: ALLSE for Whole Model', )
		session.XYDataFromHistory(name='ENERGY-9', odb=odb, 
			outputVariableName='Total energy of the output set: ETOTAL for Whole Model', )
		session.XYDataFromHistory(name='ENERGY-10', odb=odb, 
			outputVariableName='Viscous dissipation: ALLVD for Whole Model', )
		x0 = session.xyDataObjects['ENERGY-1']
		x1 = session.xyDataObjects['ENERGY-2']
		x2 = session.xyDataObjects['ENERGY-3']
		x3 = session.xyDataObjects['ENERGY-4']
		x4 = session.xyDataObjects['ENERGY-5']
		x5 = session.xyDataObjects['ENERGY-6']
		x6 = session.xyDataObjects['ENERGY-7']
		x7 = session.xyDataObjects['ENERGY-8']
		x8 = session.xyDataObjects['ENERGY-9']
		x9 = session.xyDataObjects['ENERGY-10']
		session.writeXYReport(fileName=rptName, xyData=(x0, x1, x2, x3, x4, x5, 
			x6, x7, x8, x9 ))
		#--------------------------------------------------------------------

		print '\n'
		for i in range (1, nsteps): # Loop over result steps
			print str(jobFileName)+': extracting frame '+str(i)+'\n'
			
			# Write the nodal displacement data to file
			rptName = rptName_Disp + str(i) + rptExt
			session.writeFieldReport(fileName=rptName, append=OFF, sortItem='Node Label', 
				odb=odb, step=0, frame=i, outputPosition=NODAL, displayGroup=spec_dg_nodes,
				variable=(('UT', NODAL, ((COMPONENT, 'U1'), (COMPONENT, 'U2'), )), ))
				
			if not extractDispEnerOnly:
				# Write the nodal acceleration data to file
				rptName = rptName_Accel + str(i) + rptExt
				session.writeFieldReport(fileName=rptName, append=OFF, sortItem='Node Label', 
					odb=odb, step=0, frame=i, outputPosition=NODAL, displayGroup=spec_dg_nodes,
					variable=(('AT', NODAL, ((COMPONENT, 'A1'), (COMPONENT, 'A2'), )), ))
					
				# Write strains to file, interpolate to the nodes
				rptName = rptName_Strain + str(i) + rptExt
				session.writeFieldReport(fileName=rptName, append=OFF, sortItem='Element Label', odb=odb, 
					step=0, frame=i, outputPosition=NODAL, displayGroup=spec_dg_elems,
					variable=(('E', INTEGRATION_POINT, ((COMPONENT, 'E11'), (COMPONENT, 'E22'),
					(COMPONENT, 'E12'), )), ))
						
				# Write stress to file, interpolate to the nodes
				rptName = rptName_Stress + str(i) + rptExt
				session.writeFieldReport(fileName=rptName, append=OFF, sortItem='Element Label', odb=odb, 
					step=0, frame=i, outputPosition=NODAL, displayGroup=spec_dg_elems,
					variable=(('S', INTEGRATION_POINT, ((COMPONENT, 'S11'), (COMPONENT, 'S22'),
					(COMPONENT, 'S12'), )), ))
							
		print '\n'
		print '\n'+'Data extraction complete.' +'\n'+'   '
			
		# Close the database when finished
		session.odbs[jobFileName + '.odb'].close()
		
		
	

