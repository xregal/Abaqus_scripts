!**********************************************************
!    VUMAT FOR STRAIN RATE DEPENDENT STIFFNESS MODEL, ORTHOTROPIC MATERIELS
!    WRITTEN BY HAIBIN ZHU
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
!**********************************************************************
      CHARACTER*80 CMNAME
!**********************************************************************
      INTEGER NP
      REAL BETADAMPING,Q11,Q12,Q13, Q22, Q33, Q23, Q66,A2, A6, B2,B6
      REAL ES1, ES2, ES3, ES4, SR1, SR2, SR3, SR4
      REAL sig1, sig2, sig3, sig4
!**********************************************************************  
      PARAMETER( ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, THREE = 3.D0,
     1  THIRD = ONE/THREE, HALF = .5D0, TWOTHIRDS = TWO/THREE,
     2  THREEHALFS = 1.5D0 )
	 
	 
      Q11          = PROPS(1)
      Q12          = PROPS(2)
	  Q22          = PROPS(3)
      Q13          = PROPS(4)
      Q23          = PROPS(5)
      Q33          = PROPS(6)   
      Q66          = PROPS(7)	  
      BETADAMPING  = PROPS(8)   
	  
      A2  = Q22
      B2  = 0.0
!	  6E8
      A6  = Q66
      B6  = 0.0
	  
!**********************************************************************
      DO NP = 1, NBLOCK
! Copy the values of strains into temp variables
!**********************************************************************
        ES1 = strainInc(NP,1)
        ES2 = strainInc(NP,2)
        ES3 = strainInc(NP,3)
        ES4 = strainInc(NP,4)     

!**********************************************************************  		
!       COMPUTE THE STRAIN RATE DEPEDENT STIFFNESS 
!       THE SHEAR STRAIN COMPONENT IS TENSORIAL STRAIN RATHER THAN ENGINEERING SHEAR STRAIN

!       Fibre STRAIN RATE
        SR1=1.*ABS(strainInc(NP,1)/DT)
!       TRANSVERSE STRAIN RATE
        SR2=1.*ABS(strainInc(NP,2)/DT)
        
!       SHEAR STRAIN RATE
        SR4=1.*ABS(strainInc(NP,4)/DT)
        
!       STRAIN RATE EFFECT ON TRANSVERSE COMPONENT
        Q22 = A2+B2*LOG(SR2+1)  
        Q66 = A6+B6*LOG(SR4+1)  
!        Q66=A6
      
	  
        STATENEW(NP,1) = BETADAMPING*(Q11*SR1+Q12*SR2)
        STATENEW(NP,2) = BETADAMPING*(Q12*SR1+Q22*SR2)
        STATENEW(NP,3) = 0
        STATENEW(NP,4) = BETADAMPING*(Q66*SR4*2)
		
!**********************************************************************  		
        sig1 =STRESSOLD(NP,1)+STATENEW(NP,1)-STATEOLD(NP,1)
        sig2 =STRESSOLD(NP,2)+STATENEW(NP,2)-STATEOLD(NP,2)
        sig3 =STRESSOLD(NP,3)+STATENEW(NP,3)-STATEOLD(NP,3)
        sig4 =STRESSOLD(NP,4)+STATENEW(NP,4)-STATEOLD(NP,4)	  
!**********************************************************************
        stressNew(NP,1) = sig1+Q11*ES1+Q12*ES2+Q13*ES3
        stressNew(NP,2) = sig2+Q12*ES1+Q22*ES2+Q23*ES3
        stressNew(NP,3) = sig3+Q13*ES1+Q23*ES2+Q33*ES3
        stressNew(NP,4) = sig4+Q66*ES4*2.

!    CHECK THE TIME INCREMENT    
      write(6,*) DT

      END DO
     
      RETURN
      END