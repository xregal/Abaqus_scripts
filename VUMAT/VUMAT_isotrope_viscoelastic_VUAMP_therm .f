!*********************************************************************
!    VUAMP FOR STRAIN RATE DEPENDENT                                 * 
!    WRITTEN BY XAVIER REGAL                                         *  
!*********************************************************************
!    Generate an amplitude as a progressive sine wave, it by part to *
!    be C1 since the beginin                                         *
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

!     freq0: frequency ; TempS time of establishment 
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

	  
      omega0 = 2.0*3.14159*freq0
      tim = time(iStepTime)
      period = 1.0/freq0
      	
!*********************************************************************
!         Amplitude is define by ampValue = f(t)*g(t)
!         f(t) is the sin wave, by part to have a slope nul a t=0
!         g(t) is the ramp, by part establish when reach TempS
!          
!         The 1st derivative is given 
!         The 2nd derivative is not given because is not continue
!*********************************************************************   	
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
      
      ampValueNew   = g1*f1 
      AmpDerivative	= (g1*f2)+(g2*f1)

      RETURN
      END


!*********************************************************************
!    VUMAT FOR STRAIN RATE DEPENDENT 
!    WRITTEN BY XAVIER REGAL
!*********************************************************************
      SUBROUTINE VUMAT(
! READ ONLY (UNMODIFIABLE)VARIABLES -
     1  NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL,
     2  STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH,
     3  PROPS, DENSITY, STRAININC, RELSPININC,
     4  TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,
     5  STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD,
     6  TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW,
! WRITE ONLY (MODIFIABLE) VARIABLES -
     7  STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW )
!
      INCLUDE 'VABA_PARAM.INC'

!      IMPLICIT NONE
!
      DIMENSION PROPS(NPROPS), DENSITY(NBLOCK), COORDMP(NBLOCK,*),
     1  CHARLENGTH(NBLOCK), STRAININC(NBLOCK,NDIR+NSHR),
     2  RELSPININC(NBLOCK,NSHR), TEMPOLD(NBLOCK),
     3  STRETCHOLD(NBLOCK,NDIR+NSHR),
     4  DEFGRADOLD(NBLOCK,NDIR+NSHR+NSHR),
     5  FIELDOLD(NBLOCK,NFIELDV), STRESSOLD(NBLOCK,NDIR+NSHR),
     6  STATEOLD(NBLOCK,NSTATEV), ENERINTERNOLD(NBLOCK),
     7  ENERINELASOLD(NBLOCK), TEMPNEW(NBLOCK),
     8  STRETCHNEW(NBLOCK,NDIR+NSHR),
     8  DEFGRADNEW(NBLOCK,NDIR+NSHR+NSHR),
     9  FIELDNEW(NBLOCK,NFIELDV),
     1  STRESSNEW(NBLOCK,NDIR+NSHR), STATENEW(NBLOCK,NSTATEV),
     2  ENERINTERNNEW(NBLOCK), ENERINELASNEW(NBLOCK)
!
      CHARACTER*80 CMNAME

!     Local variables
!*********************************************************************
      REAL E1, E2, E, nu, BETA, DELTA1, DELTA2, DELTA , Q11, Q12, Q33 
      REAL ES1, ES2, ES3, ES4, SR1, SR2, SR3, SR4, T1, T2 
      REAL sig1, sig2, sig3, sig4, deltaIntE, deltaIneE
!**********************************************************************  

      E1     = PROPS(1) 
      E2     = PROPS(1) 
      nu     = PROPS(2)  
	  DELTA1 = PROPS(3)
	  DELTA2 = PROPS(3)
	  OMEGA  = PROPS(4)
	  
	  T1 = 20.0
      T2 = 50.0

!-----Loop over the nodes----------------------------------------------
      DO NP = 1, NBLOCK
	  
        E = (E2-E1/(T2-T1)*TEMPOLD(NP)) 
     1    + (E1 - (E2-E21/(T2-T1)*T1) 
      
        DELTA = (DELTA1-DELTA2/(T2-T1)*TEMPOLD(NP))
     1        + (DELTA1 - (DELTA1-DELTA2/(T2-T1)*T1) 
	  
      BETA = TAN(DELTA)/OMEGA  
      Q11 = E/(1.-(nu**2.0))
      Q12 = Q11*nu 
      Q33 = Q11*0.5*(1.-nu)
	  
	  
! Copy the values of strains into temp variables

        ES1 = STRAININC(NP,1)
        ES2 = STRAININC(NP,2)
        ES4 = STRAININC(NP,4)
		
!     STRAIN RATE
        SR1=1.*(strainInc(NP,1)/DT)
        SR2=1.*(strainInc(NP,2)/DT)
        SR4=1.*(strainInc(NP,4)/DT)
        
        STATENEW(NP,1) = BETA*(Q11*SR1+Q12*SR2)
        STATENEW(NP,2) = BETA*(Q12*SR1+Q11*SR2)
        STATENEW(NP,4) = BETA*(Q33*SR4*2)
		
!**********************************************************************  		
        sig1 =STRESSOLD(NP,1)+STATENEW(NP,1)-STATEOLD(NP,1)
        sig2 =STRESSOLD(NP,2)+STATENEW(NP,2)-STATEOLD(NP,2)		
        sig4 =STRESSOLD(NP,4)+STATENEW(NP,4)-STATEOLD(NP,4)	  
!**********************************************************************
        STRESSNEW(NP,1) = sig1+Q11*ES1+Q12*ES2
        STRESSNEW(NP,2) = sig2+Q12*ES1+Q11*ES2
        STRESSNEW(NP,3) = 0
        STRESSNEW(NP,4) = sig4+Q33*ES4*2.
		
		
!**********************************************************************
!        Energy
!**********************************************************************	

!  Update the specific internal energy
        deltaIntE = 0.5*(
     1     ES1*(STRESSOLD(NP,1)+STRESSNEW(NP,1)) 
     2  +  ES2*(STRESSOLD(NP,2)+STRESSNEW(NP,2))
     3  +2*ES4*(STRESSOLD(NP,4)+STRESSNEW(NP,4)))

        ENERINTERNNEW(NP) = ENERINTERNOLD(NP) + deltaIntE/DENSITY(NP)
		
!  Update the dissipated energy
		
        deltaIneE = 0.5*(
     1     (STATENEW(NP,1)+STATEOLD(NP,1))*ES1
     2  +  (STATENEW(NP,2)+STATEOLD(NP,2))*ES2
     3  +2*(STATENEW(NP,4)+STATEOLD(NP,4))*ES4)

        ENERINELASNEW(NP) = ENERINELASOLD(NP) + deltaIneE/DENSITY(NP)

	
		
      END DO
     
      RETURN
      END