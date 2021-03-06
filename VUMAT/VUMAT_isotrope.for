C**********************************************************
C    VUMAT FOR STRAIN RATE DEPENDENT STIFFNESS MODEL, ORTHOTROPIC MATERIELS
C    WRITTEN BY HAIBIN ZHU
C**********************************************************
      SUBROUTINE VUMAT(
C READ ONLY (UNMODIFIABLE)VARIABLES -
     1  NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL,
     2  STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH,
     3  PROPS, DENSITY, STRAININC, RELSPININC,
     4  TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,
     5  STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD,
     6  TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW,
C WRITE ONLY (MODIFIABLE) VARIABLES -
     7  STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW )
C
      INCLUDE 'VABA_PARAM.INC'

C      IMPLICIT NONE
C
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
C
      CHARACTER*80 CMNAME

C     Local variables

      REAL E, nu, BETA, BETADAMPING, Q11, Q12, Q33, Q11T
      REAL ES1, ES2, ES3, ES4, SR1, SR2, SR3, SR4



      E            = PROPS(1) 
      nu           = PROPS(2)  

      Q11 = E/(1.-(nu**2.0))
      Q12 = Q11*nu 
      Q33 = Q11*0.5*(1.-nu)


      DO NP = 1, NBLOCK
C Copy the values of strains into temp variables

        ES1 = STRAININC(NP,1)
        ES2 = STRAININC(NP,2)
        ES3 = STRAININC(NP,3)
        ES4 = STRAININC(NP,4)
  	
        STRESSNEW(NP,1) = STRESSOLD(NP,1)+(Q11*ES1)+(Q12*ES2)                          
        STRESSNEW(NP,2) = STRESSOLD(NP,2)+(Q12*ES1)+(Q11*ES2)   
        STRESSNEW(NP,3) = 0 		
        STRESSNEW(NP,4) = STRESSOLD(NP,4)+Q33*(ES4) 
        
C     CHECK THE TIME INCREMENT    
        write(6,*) DT

      END DO
     
      RETURN
      END