!###############################################################################
! THIS IS A USER SUPPLIED SUBROUTINE THAT RELATES THE WEIGHTS OF EACH AB INITIO
! POINT WITH THE INPUT WW(I) VALUES
! NOTE THAT THE SIMPLEST FORM WOULD BE WEIGHT=WW
!###############################################################################

      DOUBLE PRECISION FUNCTION WEIGHT(R,Y,WW)
      USE COMMON_VAR, ONLY : NX,NP,PI,DEGS
      DOUBLE PRECISION,DIMENSION(NX) :: R
      DOUBLE PRECISION :: Y,WW,EMIN
      DOUBLE PRECISION :: COSTHETA,THETA

      COSTHETA=(R(1)**2-R(2)**2-R(3)**2)/(-2.0D+00*R(2)*R(3))
      IF (COSTHETA.LT.(-1.00D+00)) THEN 
        COSTHETA=-1.00D+00
      ELSE IF (COSTHETA.GT.(1.00D+00)) THEN
        COSTHETA=1.00D+00
      END IF
      THETA=ACOS(COSTHETA)*DEGS/PI

      WEIGHT=1.0D+00

      EMIN=-0.5D+00
      
      IF (Y.LE.-0.40D+00) THEN 
        WEIGHT=6.5D+00
      ELSE IF (Y.GT.-0.40D+00 .AND. Y.LE.-0.26D+00) THEN
        WEIGHT=4.0D+00
      ELSE IF (Y.GT.-0.26D+00 .AND. Y.LE.0.0D+00) THEN
        WEIGHT=1.0D+00!(1.0D+00-TANH(0.5*(Y-EMIN))) 
      ELSE IF (Y.GT.(0.00D+00)) THEN 
        WEIGHT=1.0D+00
      END IF

!      IF (Y.GT.(0.00D+00)) THEN
!         WEIGHT=0.5D+00
!      END IF

      !WEIGHT=(1.0D+00+EXP(-1.5*(Y-EMIN)**2)) !(1.0D+00-TANH(0.15*(Y-EMIN)))

      WRITE(1000,*) R(1),R(2),R(3),WEIGHT,THETA,Y
!
!
! 
!     THE SIMPLEST FORM WOULD BE
!      WEIGHT=WW
      END FUNCTION WEIGHT
