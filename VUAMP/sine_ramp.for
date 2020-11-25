!*********************************************************************
!    VUAMP FOR STRAIN RATE DEPENDENT 
!    WRITTEN BY XAVIER REGAL
!*********************************************************************
      Subroutine VUAMP(
!          passed in for information and state variables
     *     ampName, time, ampValueOld, dt, nprops, props, nSvars, 
     *     svars, lFlagsInfo, nSensor, sensorValues, sensorNames, 
     *     jSensorLookUpTable, 
!          to be defined
     *     ampValueNew, 
     *     lFlagsDefine,
     *     AmpDerivative, AmpSecDerivative, AmpIncIntegral)
      
      include 'vaba_param.inc'

!     svars - additional state variables, similar to (V)UEL
      dimension sensorValues(nSensor), props(nprops), 
     *     svars(nSvars)
      character*80 sensorNames(nSensor)
      character*80 ampName

!     time indices
      parameter( iStepTime        = 1,
     *           iTotalTime       = 2,
     *           nTime            = 2)
!     flags passed in for information
      parameter( iInitialization   = 1,
     *           iRegularInc       = 2,
     *           ikStep            = 3,
     *           nFlagsInfo        = 3)
!     optional flags to be defined
      parameter( iComputeDeriv     = 1,
     *           iComputeSecDeriv  = 2,
     *           iComputeInteg     = 3,
     *           iStopAnalysis     = 4,
     *           iConcludeStep     = 5,
     *           nFlagsDefine      = 5)

      parameter( freq0=20000.0, TempS = 0.001) 

      dimension time(nTime), lFlagsInfo(nFlagsInfo),
     *          lFlagsDefine(nFlagsDefine)
      dimension jSensorLookUpTable(*)


	  
      lFlagsDefine(iComputeDeriv)    = 1
      lFlagsDefine(iComputeSecDeriv) = 0
	  
!     Local variables 	  
!*********************************************************************
!     REAL f1, f2, f3, g1, g2, g3
!     REAL omega0, tim, period  
!*********************************************************************
!      freq0 = props(1)
!      TempS = props(2)

	  
      pi = 4.D0*DATAN(1.D0) 	  
      omega0 = 2.0*pi*freq0
      tim = time(iStepTime)
      period = 1.0/freq0
      	
      	
      IF (tim.LT.(0.25*period)) THEN 
        f1 =  0.5*(COS(2.0*omega0*tim-3.14159)+1.0)
        f2 = -1.0*omega0*(SIN(2.0*omega0*tim-3.14159))
        f3 = -2.0*(omega0*omega0)*(COS(2.0*omega0*tim-3.14159))
      ELSE
        f1 = SIN(omega0*tim) 
        f2 = omega0*COS(omega0*tim) 	
        f3 = -1.0*(omega0*omega0)*SIN(omega0*tim) 
      END IF
      	
      IF (tim.LT.TempS) THEN
        g1 = 1.0/TempS*tim  
        g2 = 1.0/TempS
        g3 = 0.0
      ELSE
        g1 = 1.0
        g2 = 0.0
        g3 = 0.0
      END IF
      
      ampValueNew = g1*f1 
      AmpDerivative	= (g1*f2)+(g2*f1)

      RETURN
      END
