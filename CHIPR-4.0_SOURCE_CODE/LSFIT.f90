!###############################################################################
! THIS IS THE MAIN ROUTINE TO PERFORM LEAST SQUARES FITTING
!###############################################################################

      SUBROUTINE LSFIT(NP,NC,C,FVEC,DFVEC)
      USE COMMON_VAR, ONLY : NX,M,L,X,Y,INFO 
      USE COMMON_VAR, ONLY : W,CONT,EH2CM
      IMPLICIT NONE
      INTEGER :: I,NP,NC
      DOUBLE PRECISION, DIMENSION(NC) :: C
      DOUBLE PRECISION, DIMENSION(NP) :: FVEC
      DOUBLE PRECISION, DIMENSION(NP,NC) :: DFVEC
      INTEGER :: LWA
      DOUBLE PRECISION :: TOL,DPMPAR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IWA
      DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:) :: WA

      EXTERNAL TOFIT
      
!
!     SET TOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
!     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
!     THIS IS THE RECOMMENDED SETTING.
!      

      TOL=DSQRT(DPMPAR(1))
!
!     CALLING MINPACK
!
      LWA=NP*NC+5*NC+NP

      ALLOCATE(IWA(NC),WA(LWA))

      CALL LMDIF1(TOFIT,NP,NC,C,FVEC,TOL,INFO,IWA,WA,LWA)

      DEALLOCATE(IWA,WA)

      END SUBROUTINE LSFIT
