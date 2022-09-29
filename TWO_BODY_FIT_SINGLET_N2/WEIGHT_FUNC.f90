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

      EMIN=-0.36D+00

      IF (Y.LE.-0.35D+00) THEN 
        WEIGHT=1.5D+00
      ELSE 
        WEIGHT=1.0D+00!(1.0D+00-TANH(0.5*(Y-EMIN)))
      END IF

!      
!     THE SIMPLEST FORM WOULD BE
!
      WEIGHT=WW
      END FUNCTION WEIGHT
