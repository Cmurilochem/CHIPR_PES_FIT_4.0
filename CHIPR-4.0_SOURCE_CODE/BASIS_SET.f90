      SUBROUTINE BASIS_CONTRACT(DEG,M,C,R,YVAL)
!######################################################################################################
! DEG: NUMBER OF DIFFERENT SYMMETRY UNRELATED DEGREES OF FREEDON: 
! DEG=1 FOR AN A3, DEG=2 FOR AN AB2, AND DEG=3 FOR AN ABC-TYPE MOLECULE
! M: DEFINES THE LENGHT OF CONTRACTION OF THE BASIS FUNCTIONS
! YVAL: IS THE VALUE OF THE COORDINATE IN THE CONTRACTION
!######################################################################################################
      USE COMMON_VAR, ONLY : JOBTYP,OPTTYP,NX,EVALID
      USE COMMON_VAR, ONLY : GAMAINIT,RREF0INIT,ZETAINIT
      USE COMMON_VAR, ONLY : CBSINITPOL,GAMAINITPOL
      USE COMMON_VAR, ONLY : RREFINITPOL,ZETAINITPOL
      IMPLICIT NONE
      INTEGER :: I,J,DEG
! M: DEFINES THE LENGHT OF CONTRACTION OF THE BASIS FUNCTIONS
      INTEGER :: M
      DOUBLE PRECISION :: RREF0,ZETA
      DOUBLE PRECISION, DIMENSION(2*M+2) :: C
      DOUBLE PRECISION, DIMENSION(M) :: GAMA,VAL
      DOUBLE PRECISION :: R,YVAL
!
!     STRATEGY TO ALLOW ONE TO FIX NON-LINEAR PARAMETERS
!
      IF (NX.EQ.1 .OR. NX.EQ.3 .OR. NX.EQ.6) THEN
!
!     IN CASE OF BASIS SET OPTIMIZATION
! 
        IF (JOBTYP.EQ."OPTBASIS") THEN
          IF (EVALID.EQ.0) THEN 
            IF (OPTTYP.EQ."FIXNONLIN") THEN
              C(2*M+1)=RREF0INIT
              C(2*M+2)=ZETAINIT
              DO I=1,M
                C(M+I)=GAMAINIT(I)
              END DO
            ELSE IF (OPTTYP.EQ."OPTNONLIN") THEN 
              C(2*M+1)=RREF0INIT
              C(2*M+2)=ZETAINIT
            ELSE IF (OPTTYP.EQ."OPTZETARREF") THEN
              DO I=1,M
                C(M+I)=GAMAINIT(I)
              END DO
            ELSE IF (OPTTYP.EQ."OPTALL") THEN
              CONTINUE
            END IF
          ELSE IF (EVALID.EQ.1) THEN
            CONTINUE
          ELSE IF (EVALID.EQ.2) THEN
!
!         FOR THE CASE OF EVALUATING THE DIATOMIC PART
! 
            CONTINUE            
          END IF
        END IF
!
!     IN CASE OF POLYNOMIAL OPTIMIZATION
!
        IF (JOBTYP.EQ."FITPOL" .OR. JOBTYP.EQ."DIRECTFIT") THEN
          IF (EVALID.EQ.0) THEN 
            IF (OPTTYP.EQ."FIXBAS") THEN
              C(2*M+1)=RREFINITPOL(1,DEG)
              C(2*M+2)=ZETAINITPOL(1,DEG)
              DO I=1,M
                C(I)=CBSINITPOL(I,DEG)
                C(M+I)=GAMAINITPOL(I,DEG)
              END DO
            ELSE IF (OPTTYP.EQ."FIXBASNONLIN") THEN 
              C(2*M+1)=RREFINITPOL(1,DEG)
              C(2*M+2)=ZETAINITPOL(1,DEG)
              DO I=1,M
                C(M+I)=GAMAINITPOL(I,DEG)
              END DO
            ELSE IF (OPTTYP.EQ."OPTBASNONLIN") THEN 
              C(2*M+1)=RREFINITPOL(1,DEG)
              C(2*M+2)=ZETAINITPOL(1,DEG)
            ELSE IF (OPTTYP.EQ."OPTBASZETARREF") THEN
              DO I=1,M
                C(M+I)=GAMAINITPOL(I,DEG)
              END DO
            ELSE IF (OPTTYP.EQ."OPTALL") THEN
              CONTINUE
            END IF
          ELSE IF (EVALID.EQ.1) THEN
            CONTINUE            
          ELSE IF (EVALID.EQ.2) THEN
!
!         FOR THE CASE OF EVALUATING THE DIATOMIC PART
! 
            CONTINUE            
          END IF
        END IF
!
!     END OF CASES
! 
      END IF

      DO I=1,M
        GAMA(I)=C(M+I)
      END DO

      RREF0=C(2*M+1)
      ZETA=C(2*M+2)

! SETTING UP BASIS FUNCTIONS

! SETTING SECH BASIS 
      DO I=1,M-1
        CALL PHISECBASIS(DEG,I,RREF0,ZETA,GAMA(I),1,R,VAL(I))
      END DO

! SETTING CSECH BASIS FOR LONG RANGE BEHAVIOUR 
! IN THIS VERSION, CSECH BASIS IS ALWAYS THE LAST ONE IN THE CONTRACTION
      DO I=M,M
        CALL PHICSECBASIS(DEG,I,RREF0,ZETA,GAMA(I),1,6,R,VAL(I))
      END DO

      YVAL=0.00D+00

! SETTING UP THE CONSTRACTED COORDINATES
      DO J=1,M
        YVAL=YVAL+C(J)*VAL(J)
      END DO

      RETURN
      END

      SUBROUTINE PHISECBASIS(DEG,IND,RREF0,ZETA,GAMA,ETA,R,VAL)
!######################################################################################################
! DEG: DEFINES THE DEGREE OF FREEDOM: 1,2,3...N(N-1)/2
! IND: DEFINES THE RANK NUMBER OF THE BASIS IN THE TOTAL CONTRACTION OF LENGHT M
! RREFIND: DEFINES THE DISTRIBUTED ORIGIN OF EACH BASIS ACCORDING TO RREFIND=ZETA*(RREF0)**(IND-1)
! GAMA: DEFINES THE NON-LINEAR PARAMETER DETERMINING THE DECAY OF THE SECH FUNCTION
! VAL: IS THE VALUE OF THE FUNCTION
!######################################################################################################
      IMPLICIT NONE
      INTEGER :: DEG, IND, ETA
      DOUBLE PRECISION :: RREF0, ZETA, RREFIND 
      DOUBLE PRECISION :: GAMA, R, VAL, RHO
      DOUBLE PRECISION :: SECH, ORIG
! DEFINING THE ORIGIN OF THE BASIS SET
      RREFIND=ORIG(IND,ZETA,RREF0)
! DEFINING THE DISPLACEMENT FROM THE ORIGIN
      RHO=(R-RREFIND)
! THE VALUE OF THE BASIS
      VAL=(SECH(GAMA*RHO))**(DBLE(ETA))
      RETURN
      END

      SUBROUTINE PHICSECBASIS(DEG,IND,RREF0,ZETA,GAMA,ETA,LR,R,VAL)
!######################################################################################################
! DEG: DEFINES THE DEGREE OF FREEDOM: 1,2,3...N(N-1)/2
! IND: DEFINES THE RANK NUMBER OF THE BASIS IN THE TOTAL CONTRACTION OF LENGHT M
! RREFIND: DEFINES THE DISTRIBUTED ORIGIN OF EACH BASIS ACCORDING TO RREFIND=ZETA*(RREF0)**(IND-1)
! GAMA: DEFINES THE NON-LINEAR PARAMETER DETERMINING THE DECAY OF THE SECH FUNCTION
! VAL: IS THE VALUE OF THE FUNCTION
! LR: IS THE EXPONENT OF THE TANH PART
!######################################################################################################
      IMPLICIT NONE
      INTEGER :: DEG, IND, LR, ETA
      DOUBLE PRECISION :: RREF0, ZETA, RREFIND 
      DOUBLE PRECISION :: GAMA, R, VAL, RHO
      DOUBLE PRECISION :: SECH, BETA, FAC, ORIG
      BETA=1.00D+00/5.0D+00
! DEFINING THE ORIGIN OF THE BASIS SET
      RREFIND=ORIG(IND,ZETA,RREF0)
! DEFINING THE DISPLACEMENT FROM THE ORIGIN
      RHO=(R-RREFIND)
! DEFINING THE TANH FACTOR
      FAC=(TANH(BETA*R)/R)**(DBLE(LR))
! THE VALUE OF THE BASIS
      VAL=FAC*(SECH(GAMA*RHO))**(DBLE(ETA))
      RETURN
      END

!######################################################################################################
! DEFINING THE ORIGIN OF EACH BASIS FUNCTION
!######################################################################################################

      DOUBLE PRECISION FUNCTION ORIG(IND,ZETA,RREF0)
      IMPLICIT NONE
      INTEGER :: IND
      DOUBLE PRECISION :: ZETA,RREF0
      ORIG=ZETA*(RREF0)**(DBLE(IND)-1.0D+00)
      RETURN
      END

!######################################################################################################
! DEFINING THE HYPERBOLIC SECANT FUNCTION
!###################################################################################################### 

      DOUBLE PRECISION FUNCTION SECH(X)
      IMPLICIT NONE
      DOUBLE PRECISION :: X
      SECH=1.00D+00/(COSH(X))      
      RETURN
      END
