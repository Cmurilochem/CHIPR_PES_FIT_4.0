!###############################################################################
! THIS IS A USER SUPPLIED SUBROUTINE THAT RELATES THE WEIGHTS OF EACH AB INITIO
! POINT WITH THE INPUT WW(I) VALUES
! NOTE THAT THE SIMPLEST FORM WOULD BE WEIGHT=WW
!###############################################################################

      DOUBLE PRECISION FUNCTION WEIGHT(R,Y,WW)
      USE COMMON_VAR, ONLY : NX,NP,PI,DEGS
      DOUBLE PRECISION,DIMENSION(NX) :: R
      DOUBLE PRECISION :: Y,WW
      DOUBLE PRECISION :: COSTHETA,THETA
      COSTHETA=(R(3)**2-R(1)**2-R(2)**2)/(-2.0D+00*R(1)*R(2))
      IF (COSTHETA.LT.(-1.00D+00)) THEN 
        COSTHETA=-1.00D+00
      ELSE IF (COSTHETA.GT.(1.00D+00)) THEN
        COSTHETA=1.00D+00
      END IF
      THETA=ACOS(COSTHETA)*DEGS/PI
      WEIGHT=1.0D+00
!
!
!
      IF (THETA.GE.(59.0D+00) .AND. THETA.LE.(71.00D+00)) THEN
        IF (R(1).GE.(2.2D+00) .AND. R(1).LE.(3.00D+00)) THEN
          IF (R(2).GE.(2.2D+00) .AND. R(2).LE.(3.00D+00)) THEN
            WEIGHT=2.0D+00
          END IF
        END IF 
      END IF
!
!
!
      IF (THETA.GE.(99.0D+00) .AND. THETA.LE.(121.00D+00)) THEN
        IF (R(1).GE.(2.2D+00) .AND. R(1).LE.(3.00D+00)) THEN
          IF (R(2).GE.(2.2D+00) .AND. R(2).LE.(3.00D+00)) THEN
            WEIGHT=2.0D+00
          END IF
        END IF 
      END IF
!
!
!
      IF (THETA.GE.(169.0D+00) .AND. THETA.LE.(181.00D+00)) THEN
        IF (R(1).GE.(2.2D+00) .AND. R(1).LE.(3.00D+00)) THEN
          IF (R(2).GE.(2.2D+00) .AND. R(2).LE.(3.00D+00)) THEN
            WEIGHT=2.0D+00
          END IF
        END IF 
      END IF
!
!
! 
      IF (Y.GT.(0.00D+00)) THEN
         WEIGHT=0.5D+00
      END IF

!      WRITE(1000,*) R(1),R(2),R(3),WEIGHT,THETA,Y
!
!
! 
!     THE SIMPLEST FORM WOULD BE
!      WEIGHT=WW
      END FUNCTION WEIGHT
