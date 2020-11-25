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
	rptName_TotEnergy = directory+extension+jobFileName+extractFolder +jobFileName+'_TotEnergy'
	if not extractDispEnerOnly:
		rptName_Vel = directory+extension+jobFileName+extractFolder+jobFileName+'_Vel'
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
			# Write the nodal velocity data to file
			rptName = rptName_Vel + str(i) + rptExt
			session.writeFieldReport(fileName=rptName, append=OFF, sortItem='Node Label', 
				odb=odb, step=0, frame=i, outputPosition=NODAL, displayGroup=spec_dg_nodes,
				variable=(('VT', NODAL, ((COMPONENT, 'V1'), (COMPONENT, 'V2'), )), ))
			
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
	
	
	

