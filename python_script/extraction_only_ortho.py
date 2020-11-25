#Just for extraction

# Create folders for extracted values and later for images
extractFolder='/Extracted_Data/'
if not os.path.exists(directory+extension+jobFileName+extractFolder):
	os.makedirs(directory+extension+jobFileName+extractFolder)	

# Open the odb file
odb=session.odbs[jobFileName +'.odb']

# Create a global co-ord system
odb.rootAssembly.DatumCsysByThreePoints(name='CSYS-1',coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), point2=(0.0, 1.0, 0.0))
	
dtm = odb.rootAssembly.datumCsyses['CSYS-1']
session.viewports['ExpDynSonotrode2D'].odbDisplay.basicOptions.setValues(transformationType=USER_SPECIFIED, datumCsys=dtm)

nsteps=len(odb.steps['Sonotrode'].frames)-1

for i in range(nsteps-130,nsteps):
#for i in range(1,nsteps):
	fileName_E   =directory+extension+jobFileName+extractFolder+'E_'+str(i)+'.rpt'
	fileName_ER  =directory+extension+jobFileName+extractFolder+'ER_'+str(i)+'.rpt'
	fileName_S   =directory+extension+jobFileName+extractFolder+'S_'+str(i)+'.rpt'
	fileName_U   =directory+extension+jobFileName+extractFolder+'U_'+str(i)+'.rpt'   
	fileName_A   =directory+extension+jobFileName+extractFolder+'A_'+str(i)+'.rpt'  
	session.writeFieldReport(fileName=fileName_E, append=OFF, sortItem='Element Label', odb=odb, step=0, frame=i, outputPosition=NODAL, variable=(('E', INTEGRATION_POINT, ((COMPONENT, 'E11'), (COMPONENT, 'E22'), (COMPONENT, 'E12'), )), ))
	session.writeFieldReport(fileName=fileName_ER, append=OFF, sortItem='Element Label', odb=odb, step=0, frame=i, outputPosition=NODAL, variable=(('ER', INTEGRATION_POINT, ((COMPONENT, 'ER11'), (COMPONENT, 'ER22'), (COMPONENT, 'ER12'), )), ))
	session.writeFieldReport(fileName=fileName_S, append=OFF, sortItem='Element Label', odb=odb, step=0, frame=i, outputPosition=NODAL,  variable=(('S', INTEGRATION_POINT, ((COMPONENT, 'S11'), (COMPONENT, 'S22'), (COMPONENT, 'S12'), )), ))
	session.writeFieldReport(fileName=fileName_U, append=OFF, sortItem='Node Label', odb=odb, step=0, frame=i, outputPosition=NODAL, variable=(('COORD', NODAL, ((COMPONENT, 'COOR1'), (COMPONENT, 'COOR2'), )),('UT', NODAL, ((COMPONENT, 'U1'), (COMPONENT, 'U2'), )), ))
	session.writeFieldReport(fileName=fileName_A, append=OFF, sortItem='Node Label', odb=odb, step=0, frame=i, outputPosition=NODAL, variable=(('AT', NODAL, ((COMPONENT, 'A1'), (COMPONENT, 'A2'), )), ))

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
f.write(str(specE11)+'\n');
f.write(str(specE22)+'\n');
f.write(str(specNu12)+'\n');
f.write(str(beta)+'\n');
f.write(str(freq0)+'\n');
f.write(str(disp_ampl)+'\n');
f.write(str(specLength)+'\n');
f.write(str(specHeight)+'\n');
f.write(str(specThick)+'\n');
f.write(str(betaDamping)+'\n');
f.write(str(specAngle)+'\n');
f.close()