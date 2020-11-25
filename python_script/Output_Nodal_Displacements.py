from abaqusConstants import *
from odbAccess import *
import visualization
odb_name = 'Job-1.odb'
session.openOdb(odb_name)
odb=session.odbs[odb_name]
nsteps=len(odb.steps['Step-1'].frames)

for i in range(1,nsteps):
    fileName_Ex='exportfolder/'+'Ex_'+str(i)+'.rpt'
    fileName_Ey='exportfolder/'+'Ey_'+str(i)+'.rpt'
    fileName_Exy='exportfolder/'+'Exy_'+str(i)+'.rpt'
    fileName_area = 'exportfolder/'+'area_'+str(i)+'.rpt'
    fileName_sx='exportfolder/'+'sx_'+str(i)+'.rpt'
    fileName_sy='exportfolder/'+'sy_'+str(i)+'.rpt'
    fileName_sxy='exportfolder/'+'sxy_'+str(i)+'.rpt'
    fileName_svm='exportfolder/'+'svm_'+str(i)+'.rpt' 
    fileName_peeq='exportfolder/'+'peeq_'+str(i)+'.rpt'


    
    session.writeFieldReport(fileName=fileName_Ex, append=OFF, sortItem='Element Label', odb=odb, step=0, frame=i, outputPosition=INTEGRATION_POINT, variable=(('E', INTEGRATION_POINT, ((COMPONENT, 'E11'), )), ))
    session.writeFieldReport(fileName=fileName_Ey, append=OFF, sortItem='Element Label', odb=odb, step=0, frame=i, outputPosition=INTEGRATION_POINT, variable=(('E', INTEGRATION_POINT, ((COMPONENT, 'E22'), )), ))
    session.writeFieldReport(fileName=fileName_Exy, append=OFF, sortItem='Element Label', odb=odb, step=0, frame=i, outputPosition=INTEGRATION_POINT, variable=(('E', INTEGRATION_POINT, ((COMPONENT, 'E12'), )), ))
#    session.writeFieldReport(fileName=fileName_sx, append=OFF, sortItem='Element Label', odb=odb, step=0, frame=i, outputPosition=ELEMENT_CENTROID, variable=(('S', INTEGRATION_POINT, ((COMPONENT, 'S11'), )), ))
#    session.writeFieldReport(fileName=fileName_sy, append=OFF, sortItem='Element Label', odb=odb, step=0, frame=i, outputPosition=ELEMENT_CENTROID, variable=(('S', INTEGRATION_POINT, ((COMPONENT, 'S22'), )), ))
#    session.writeFieldReport(fileName=fileName_sxy, append=OFF, sortItem='Element Label', odb=odb, step=0, frame=i, outputPosition=ELEMENT_CENTROID, variable=(('S', INTEGRATION_POINT, ((COMPONENT, 'S12'), )), ))
#    session.writeFieldReport(fileName=fileName_peeq, append=OFF, sortItem='Element Label', odb=odb, step=0, frame=i, outputPosition=ELEMENT_CENTROID, variable=(('PEEQ', INTEGRATION_POINT, ), ))


f = open('exportfolder/nTimeSteps.txt', 'w')
f.write(str(nsteps));
f.close()
