!**********************************************************
!    VUMAT FOR STRAIN RATE DEPENDENT 
!    WRITTEN BY XAVIER REGAL
!**********************************************************
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
      REAL E, nu, BETA, BETADAMPING, Q11, Q12, Q33
      REAL ES1, ES2, ES3, ES4, SR1, SR2, SR3, SR4
      REAL sig1, sig2, sig3, sig4, deltaEE, deltaIE
!**********************************************************************  

      E     = PROPS(1) 
      nu    = PROPS(2)  
	  DELTA  = PROPS(3)
	  OMEGA  = PROPS(4)
	  
      BETA = TAN(DELTA)/OMEGA  
      Q11 = E/(1.-(nu**2.0))
      Q12 = Q11*nu 
      Q33 = Q11*0.5*(1.-nu)


      DO NP = 1, NBLOCK
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
        stressNew(NP,1) = sig1+Q11*ES1+Q12*ES2
        stressNew(NP,2) = sig2+Q12*ES1+Q11*ES2
        stressNew(NP,3) = 0
        stressNew(NP,4) = sig4+Q33*ES4*2.
		
		
!**********************************************************************
!        Energy
!**********************************************************************	
!  Update the specific internal energy
        deltaIntE = 0.5*(
     1       ES1*(STRESSOLD(NP,1)+STRESSNEW(NP,1)) 
     2  +    ES2*(STRESSOLD(NP,2)+STRESSNEW(NP,2))
     3  +2.0*ES4*(STRESSOLD(NP,4)+STRESSNEW(NP,4)))

        ENERINTERNNEW(NP) = ENERINTERNOLD(NP) + deltaIntE/DENSITY(NP)
		
!  Update the dissipated energy
		
        deltaIneE = 0.5*(
     1       (STATENEW(NP,1)+STATEOLD(NP,1))*ES1
     2  +    (STATENEW(NP,2)+STATEOLD(NP,2))*ES2
     3  +2.0*(STATENEW(NP,4)+STATEOLD(NP,4))*ES4)

        ENERINELASNEW(NP) = ENERINELASOLD(NP) + deltaIneE/DENSITY(NP)
		
      END DO
     
      RETURN
      END