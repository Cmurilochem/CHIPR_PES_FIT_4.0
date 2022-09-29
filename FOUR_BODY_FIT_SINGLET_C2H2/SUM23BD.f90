!###############################################################################
! THIS IS A USER SUPPLIED SUBROUTINE THAT CALCULATES THE SUM OF TWO- AND THREE-BOBY ENERGIES
! FOR A GIVEN TETRATOMIC MOLECULE
!###############################################################################

      SUBROUTINE SUM23BD(R1,R2,R3,R4,R5,R6,SUMV23)
      IMPLICIT NONE
      INTEGER :: I
      DOUBLE PRECISION :: R1,R2,R3,R4,R5,R6,SUMV23
      DOUBLE PRECISION, DIMENSION(6) :: W
!
!     R1 IS ALWAYS H-H DISTANCES
!     R3-R5 ARE H-C DISTANCES
!     R6 IS ALWAYS C-C DISTANCES
!
      W(1)=R1
      W(2)=R2
      W(3)=R3
      W(4)=R4
      W(5)=R5
      W(6)=R6
!
      CALL POT23C2H2(W,SUMV23)
      RETURN 
      END
