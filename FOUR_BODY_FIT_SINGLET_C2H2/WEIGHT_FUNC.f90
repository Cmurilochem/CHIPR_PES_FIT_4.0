!###############################################################################
! THIS IS A USER SUPPLIED SUBROUTINE THAT RELATES THE WEIGHTS OF EACH AB INITIO
! POINT WITH THE INPUT WW(I) VALUES
! NOTE THAT THE SIMPLEST FORM WOULD BE WEIGHT=WW
!###############################################################################

      DOUBLE PRECISION FUNCTION WEIGHT(R,Y,WW)
      USE COMMON_VAR, ONLY : NX,NP,PI,DEGS
      INTEGER :: I
      DOUBLE PRECISION,DIMENSION(NX) :: R
      DOUBLE PRECISION :: Y,WW,EMIN

      EMIN=-0.645D+00

      IF (Y.LE.-0.59D+00) THEN 
        WEIGHT=1.85D+00
      ELSE IF (Y.GT.-0.59D+00 .AND. Y.LE.-0.55D+00) THEN  
        WEIGHT=1.5D+00
      ELSE IF (Y.GT.-0.55D+00 .AND. Y.LE.-0.30D+00)  THEN
        WEIGHT=1.25D+00
      ELSE 
        WEIGHT=0.8D+00!(1.0D+00-TANH(1.0*(Y-EMIN)))
      END IF
      
      WRITE(1001,*) R,Y,WEIGHT

!      
!     THE SIMPLEST FORM WOULD BE
!
!     WEIGHT=WW
      END FUNCTION WEIGHT
