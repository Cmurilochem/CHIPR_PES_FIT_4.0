!###############################################################################
! THIS IS A USER SUPPLIED SUBROUTINE THAT RELATES THE WEIGHTS OF EACH AB INITIO
! POINT WITH THE INPUT WW(I) VALUES
! NOTE THAT THE SIMPLEST FORM WOULD BE WEIGHT=WW
!###############################################################################

      DOUBLE PRECISION FUNCTION WEIGHT(R,Y,WW)
      USE COMMON_VAR, ONLY : NX,NP,PI,DEGS
      INTEGER :: I
      DOUBLE PRECISION,DIMENSION(NX) :: R
      DOUBLE PRECISION :: Y,WW
!      
!     THE SIMPLEST FORM WOULD BE
!
      IF (Y.LE.-0.15) THEN
        WEIGHT=1.5D+00
      ELSE IF (Y.GE.0.1) THEN  
        WEIGHT=WW!0.25D+00 
      ELSE
        WEIGHT=WW
      END IF
      END FUNCTION WEIGHT
