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
      WEIGHT=WW
      END FUNCTION WEIGHT
