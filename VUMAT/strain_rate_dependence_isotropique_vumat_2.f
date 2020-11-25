!**********************************************************************
!    VUMAT FOR STRAIN RATE DEPENDENT STIFFNESS MODEL,                 *
!    Isotropic material                                               *
!    WRITTEN BY Xavier RÉGAL after HAIBIN ZHU      2019               *
!**********************************************************************
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
!     Local variables
!**********************************************************************
      INTEGER NP
      REAL E_ANTE,nu,BETA,ALPHADAMPING,BETADAMPING,Q11,Q12,Q33
      REAL ES1, ES2, ES3, ES4, SR1, SR2, SR3, SR4
      REAL sig1, sig2, sig3, sig4

!     Reading the material parameters
      E_ANTE       = PROPS(1)
      nu           = PROPS(2)
      BETA         = PROPS(3)
      BETADAMPING  = PROPS(4)
!	  ALPHADAMPING = PROPS(5)
!**********************************************************************     

      DO NP = 1, NBLOCK
! Copy the values of strains increment into temporary variables
        ES1 = STRAININC(NP,1)
        ES2 = STRAININC(NP,2)
        ES3 = STRAININC(NP,3)
        ES4 = STRAININC(NP,4)
		
! Determine the all strains 		

        STATENEW(NP,5) = STATEOLD(NP,5) + ES1
        STATENEW(NP,6) = STATEOLD(NP,6) + ES2
        STATENEW(NP,7) = STATEOLD(NP,7) + ES3
        STATENEW(NP,8) = STATEOLD(NP,8) + ES4

        
! Copy the values of strainrates into temporary variables
        SR1=(STRAININC(NP,1)/DT)
        SR2=(STRAININC(NP,2)/DT)
        SR3=(STRAININC(NP,3)/DT)
        SR4=(STRAININC(NP,4)/DT)
		
!     Defining the material behaviour

        Q11 = (E_ANTE+(BETA*ALOG(1.0*ABS(SR1)+1.0)))/(1.-(nu**2.0))
        Q12 = Q11*nu 
        Q33 = Q11*0.5*(1.-nu)
		
        STATENEW(NP,1) = BETADAMPING*(Q11*SR1+Q12*SR2)
        STATENEW(NP,2) = BETADAMPING*(Q12*SR1+Q11*SR2)
        STATENEW(NP,3) = 0
        STATENEW(NP,4) = BETADAMPING*(Q33*SR4)
		
!**********************************************************************  		
!        sig1 =STRESSOLD(NP,1)+STATENEW(NP,1)-STATEOLD(NP,1)
!        sig2 =STRESSOLD(NP,2)+STATENEW(NP,2)-STATEOLD(NP,2)
!        sig3 =STRESSOLD(NP,3)+STATENEW(NP,3)-STATEOLD(NP,3)
!        sig4 =STRESSOLD(NP,4)+STATENEW(NP,4)-STATEOLD(NP,4)
        sig1 = STATENEW(NP,1)-STATEOLD(NP,1)
        sig2 = STATENEW(NP,2)-STATEOLD(NP,2)
        sig3 = STATENEW(NP,3)-STATEOLD(NP,3)
        sig4 = STATENEW(NP,4)-STATEOLD(NP,4)
		
		
        STRESSNEW(NP,1) = sig1+Q11*STATENEW(NP,5)+Q12*STATENEW(NP,6)
        STRESSNEW(NP,2) = sig2+Q12*STATENEW(NP,5)+Q11*STATENEW(NP,6)
        STRESSNEW(NP,3) = 0
        STRESSNEW(NP,4) = sig4+Q33*STATENEW(NP,8)
        
!     CHECK THE TIME INCREMENT    
        write(6,*) DT
!**********************************************************************  

      END DO
     
      RETURN
      END