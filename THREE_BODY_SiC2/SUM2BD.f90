!###############################################################################
! THIS IS A USER SUPPLIED SUBROUTINE THAT CALCULATES THE SUM OF TWO-DOBY ENERGIES
! FOR A GIVEN TRIATOMIC MOLECULE
!###############################################################################

      SUBROUTINE SUM2BD(R1,R2,R3,SUMV2)
      IMPLICIT NONE
      DOUBLE PRECISION :: R1,R2,R3,SUMV2
      DOUBLE PRECISION :: V1,V2,V3
      DOUBLE PRECISION :: DIAT1,DIAT2,DIAT3
      V1=DIAT1(R1)
      V2=DIAT2(R2)
      V3=DIAT3(R3)
      SUMV2=V1+V2+V3
      RETURN 
      END

!###############################################################################
! THE USER MUST DEFINE THESE FUNCTIONS:
! DIAT1(R,V21): THE DIATOMIC CURVE FOR R1 COORDINATE
! DIAT2(R,V21): THE DIATOMIC CURVE FOR R2 COORDINATE
! DIAT3(R,V21): THE DIATOMIC CURVE FOR R3 COORDINATE
!
! NOTE THAT, EVEN FOR AN A3 MOLECULE, THE DIATOMIC CURVES SHOULD BE DEFINED
! IN THIS CASE, ALL THREE FUNCTIONS ARE DEFINED BY THE SAME SUBROUTINE
!
! FOR AN AB2, R1 IS DEFINED BY ONE DIFFERENT SUBROUTINE, AND TWO OTHER EQUAL FOR 
! R2 AND R3
!
! FOR AN ABC, ALL SUBROUTINES ARE DIFFERENT
!
!###############################################################################

      DOUBLE PRECISION FUNCTION DIAT1(R1)
      IMPLICIT NONE
      DOUBLE PRECISION :: R1,V21,DVDR
!
!     DEFINE HERE YOUR SUBROUTINE FOR DIATOMIC R1 
!
      CALL CHIPR_C2(R1,V21,.FALSE.,DVDR)

      DIAT1=V21

      RETURN 
      END

!###############################################################################

      DOUBLE PRECISION FUNCTION DIAT2(R2)
      IMPLICIT NONE
      DOUBLE PRECISION :: R2,V22,DVDR
!
!     DEFINE HERE YOUR SUBROUTINE FOR DIATOMIC R2
!
      CALL CHIPR_SiC(R2,V22,.FALSE.,DVDR)

      DIAT2=V22

      RETURN 
      END

!###############################################################################

      DOUBLE PRECISION FUNCTION DIAT3(R3)
      IMPLICIT NONE
      DOUBLE PRECISION :: R3,V23,DVDR
!
!     DEFINE HERE YOUR SUBROUTINE FOR DIATOMIC R3 
!
      CALL CHIPR_SiC(R3,V23,.FALSE.,DVDR)

      DIAT3=V23

      RETURN 
      END