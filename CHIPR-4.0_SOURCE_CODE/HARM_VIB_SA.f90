!######################################################################################################
! THIS SUBROUTINE TRIES TO FIND THE GLOBAL MINIMUM OF THE FITTED SURFACE 
! USING A SIMULATED ANNELING CODE. THE OPTIMUM STRUCTURE SO FOUND ARE SUBJECT TO 
! HARMONIC VIBRATIONAL ANALYSIS USING WILSON'S FG FORMALISM
!######################################################################################################
      SUBROUTINE GLOBALOPT
      USE COMMON_VAR, ONLY : NATOM,NATOMMAX,CARTX
      USE COMMON_VAR, ONLY : NX,ATOM,AMU,JOBTYP
      IMPLICIT NONE
      INTEGER :: I,J,ERR
      DOUBLE PRECISION :: POT,GRAD,GRADMAX
      LOGICAL :: CHECK
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: CARTXOPT,Q
      DOUBLE PRECISION :: GMENERGY
      INTEGER :: IERGM
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: PDOT,EIG
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX,3*NATOMMAX) :: FA,A,VEC
      DOUBLE PRECISION, DIMENSION(NATOMMAX) :: W
      CHARACTER(2), DIMENSION(NATOMMAX) :: TMPATOM
      
      IF (JOBTYP.EQ."OPTBASIS") THEN

!     DO NOTHING
      
      ELSE IF (JOBTYP.EQ."FITPOL" .OR. JOBTYP.EQ."DIRECTFIT") THEN  
      
        WRITE(6,100)
        WRITE(6,10)
!
!     OPENING THE FILE FOR SA OUTPUT
!
        OPEN(UNIT=123,FILE='sa.res',IOSTAT=ERR,&
     &   FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
        IF (ERR /= 0) THEN
        PRINT*,"PROBLEMS IN OPENING THE FILE SA.RES"
        END IF      
!
!     TESTING THE FUNCTION
!
        CALL FCNSA(3*NATOM,CARTX,POT)        
!
!     CALLING SIMULATED ANNELING  
!
        DO I=1,3*NATOMMAX
          CARTXOPT(I)=0.0D+00
        END DO
        
        CALL SIMANN(3*NATOM,CARTX,CARTXOPT,GMENERGY,IERGM)
        
        IF (IERGM .EQ. 0) THEN
!      
! WRITTING THE FINAL STRUCTURE TO A FILE
!     
          OPEN(UNIT=124,FILE='gmin.res',IOSTAT=ERR,&
     &     FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
            IF (ERR /= 0) THEN
              PRINT*, "PROBLEMS IN OPENING THE FILE GMIN.RES"
            END IF
!
! WRITTING THE NORMAL MODE ANALYSIS TO A FILE
!     
          OPEN(UNIT=125,FILE='normalmode.res',IOSTAT=ERR,&
     &     FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
            IF (ERR /= 0) THEN
              PRINT*, "PROBLEMS IN OPENING THE FILE NORMALMODE.RES"
            END IF              
!
! USING RICHARDSON'S EXTRAPOLATION FORMULA TO CALCULATE THE GRADIENT
!
          DO I=1,3*NATOMMAX
            PDOT(I)=0.0D+00
            EIG(I)=0.0D+00
            Q(I)=0.0D+00
            DO J=1,3*NATOMMAX
              FA(I,J)=0.0D+00
              A(I,J)=0.0D+00
              VEC(I,J)=0.0D+00
            END DO
          END DO
!
! TRANFERING ATOMIC MASSES TO W VECTOR AND ATOM CHARACTERS TO TMPATOM
!
          DO I=1,NATOM
            W(I)=AMU(I)
            TMPATOM=ATOM(I)
          END DO
          
          CALL DER(CARTXOPT,PDOT,3*NATOM)
!        
! CALCULATING THE MODULUS OF THE GRADIENT VECTOR AND THE MAXIMUM COMPONENT
!
          GRAD=0.0D+00
          DO I=1,3*NATOM
            GRAD=GRAD+(PDOT(I)**2)            
          END DO
          GRAD=SQRT(GRAD)
          GRADMAX=MAXVAL(PDOT)
          
          WRITE(6,11)
          WRITE(6,12) GMENERGY
          WRITE(6,110) GRAD
          WRITE(6,210) GRADMAX
          WRITE(6,113)
!        
! CALCULATING THE FORCE CONSTANT MATRIX (IN A.U)
!
          CALL FMTRX(CARTXOPT,PDOT,FA,A,3*NATOM)
!
! CONSTRUCTING AN EFFECTIVE (MASS-WEIGHTED) WILSON -FG- MATRIX
!
          CALL FGMTRX(FA,A,W,NATOM,EIG,VEC)
!
! WRITING THE OPTIMIZED STRUCTURE TO A GMIN.RES FILE
!
          CALL WRITMOLDEN(1,CARTXOPT,EIG,VEC,NATOM,ATOM)
!
! CONVERTING CARTESIAN TO INTERNAL AGAIN
!
          CALL C2IN(CARTXOPT,Q,3*NATOM)
        
          WRITE(6,13)
          DO I=1,NX
            WRITE(6,14) I, Q(I)
          END DO
          WRITE(6,113)
          WRITE(6,111)
          DO I=1,3*NATOM
            WRITE(6,112) I, EIG(I)
          END DO        
          WRITE(6,113) 
          WRITE(6,114)
          WRITE(6,115)        
                
        ELSE
      
          WRITE(6,100)
          WRITE(6,116)
        
        END IF
        
      END IF  

   10 FORMAT ("ATTEMPTING TO FIND THE GLOBAL MINIMUM",/)    
   11 FORMAT ("GLOBAL MINIMUM FOUND",/)
   12 FORMAT ("ENERGY (IN HARTREE):",12X,F12.5)
   13 FORMAT ("MINIMUM GEOMETRY (IN BOHR):",/)
   14 FORMAT (15X,"R(",I1,")=",1X,F15.5)
  110 FORMAT ("GRADIENT NORM IS (IN HARTREE/BOHR):",6E12.4,/)
  210 FORMAT ("GRAD. MAX. (IN HARTREE/BOHR):",6X,6E12.4)
  111 FORMAT ("HARMONIC FREQUENCIES ARE (IN CM-1):",/)
  112 FORMAT (15X,I3,F15.2)
  113 FORMAT (A)
  114 FORMAT ("NOTE: THERE ARE 6 TRANSLATIONS AND ROTATIONS.")
  115 FORMAT ("FOR LINEAR MOLECULES ONLY 2 ROTATIONS ARE DEFINED")
  116 FORMAT("ATTEMPT TO FIND GLOBAL MINIMUM FAILED: SEE SA.RES FILE",/)  
  100  FORMAT (/,80("#"),/,80("#"),/) 
      
      CLOSE(UNIT=123)
      INQUIRE(UNIT=124,OPENED=CHECK)
      IF (CHECK) THEN 
        CLOSE(UNIT=124)
      ELSE
      END IF
      INQUIRE(UNIT=125,OPENED=CHECK)
      IF (CHECK) THEN 
        CLOSE(UNIT=125)
      ELSE
      END IF      
      
      RETURN
      END

!###############################################################################
! THIS SUBROUTINE IS THE MAIN INTERFACE TO THE SIMULATED ANNELING PROGRAM
!###############################################################################            
      SUBROUTINE FCNSA(N,X,POT)
!###############################################################################      
! N IS 3X THE NUMBER OF ATOMS
! X IS THE INITIAL CARTESIAN COORDINATES
! POT IS THE POTENTIAL
!###############################################################################!
      USE COMMON_VAR, ONLY : NATOMMAX,NX,NC,MBS,L,C
      IMPLICIT NONE
      INTEGER :: I, N
      DOUBLE PRECISION, DIMENSION(NX) :: W
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: X, R
      DOUBLE PRECISION, DIMENSION(NC) :: DC
      DOUBLE PRECISION :: TWOBD,THREBD,TWOTHREBD,FOURBD,POT
      CALL C2IN(X,R,N)
      DO I=1,NX
        W(I)=R(I)
      END DO
      
      IF (NX.EQ.1) THEN
      
        CALL CHIPR_DIAT(MBS,L,C,W,POT,DC)
       
      ELSE IF (NX.EQ.3) THEN 
      
        CALL SUM2BD(W(1),W(2),W(3),TWOBD)
        CALL CHIPR_TRIAT(MBS,L,C,W,THREBD,DC)
        POT=TWOBD+THREBD
      
      ELSE IF (NX.EQ.6) THEN
      
        CALL SUM23BD(W(1),W(2),W(3),&
     &         W(4),W(5),W(6),TWOTHREBD)
        CALL CHIPR_TETRA(MBS,L,C,W,FOURBD,DC)
        POT=TWOTHREBD+FOURBD
      
      END IF
        
      RETURN
      END
      
!####################################################################################
! CONVERT CARTESIAN TO INTERPARTICLE DISTANCES (A.U)
!####################################################################################

      SUBROUTINE C2IN(X,R,NTA)
      IMPLICIT NONE
      INTEGER, PARAMETER :: NA=100,NT=3*NA
      INTEGER :: NTA,I,K,J
      DOUBLE PRECISION, DIMENSION (NT):: X,R
      K=0
      DO I=1,NTA,3
       DO J=I,NTA,3
        IF(I.NE.J)THEN
         K=K+1
         R(K)=SQRT((X(I)-X(J))**2+(X(I+1)-X(J+1))**2+&
     &   (X(I+2)-X(J+2))**2)
        ENDIF
       ENDDO
      ENDDO
      RETURN
      END      
      
!###############################################################################
! CONSTRUCT AN EFFECTIVE WILSON -FG- MATRIX BY INTRODUCING MASS 
! DEPENDENCE INTO THE POTENTIAL FORCE CONSTANT MATRIX. THE NORMAL
! MODES AND THE SPECTROSCOPIC FREQUENCIES ARE THEN EVALUATED.          
!###############################################################################

      SUBROUTINE FGMTRX(FA,A,W,NATOM,EIG,VEC)
      IMPLICIT NONE
      INTEGER :: I, J, K
      INTEGER, PARAMETER :: NATOMMAX=100
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX,3*NATOMMAX) :: FA, A 
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX,3*NATOMMAX) :: VEC, DB
      DOUBLE PRECISION, DIMENSION(NATOMMAX) :: W
      INTEGER :: NATOM, NCARTCOORD
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: DIM, DA, EIG
      DOUBLE PRECISION :: RHO, JH, B2A, A2M, AMU2KG, PI, C, TFACT
      DOUBLE PRECISION :: DUM, SGN
  101 FORMAT(A)
  102 FORMAT(/,10X,'EIGENVECTORS OF CARTESIAN FORCE CONSTANT MATRIX'/)
  103 FORMAT(/,18X,' NORMAL MODES AND VIBRATIONAL FREQUENCIES (CM-1)')
  104 FORMAT(/)

      NCARTCOORD=3*NATOM
      RHO=1.0D-08
! CONVERTION FROM HARTREES TO JOULES
      JH=4.3597439D-18
! CONVERTION FROM BOHR TO ANGS
      B2A=0.5291772086D+00
! CONVERTION FROM ANGS TO METER
      A2M=1.0D-10
! CONVERTION FROM A.M.U TO KG
      AMU2KG=1.660538921D-27
! SPEED OF LIGHT IN CM.S-1
      C=2.99792458D+10
      PI=3.14159265358979D+00
! FACTOR TO MULTIPLY THE EIGENVALUES
! FACTOR IS IN CM-2
! EIGENVALUES ARE IN Eh.BOHR-2
      TFACT=(JH)/((B2A*A2M)**2*AMU2KG*4.0D+0*(PI**2)*(C**2))

! CALCULATE ARRAY DIM USED FOR MASS-WEIGHTING

      K=0
      DO 51 I=1,NATOM
         DO 50 J=1,3
            K=K+1
            DIM(K)=1.D0/DSQRT(W(I))
   50    CONTINUE
   51 CONTINUE

! MOVE DIAGONAL OF -A- TO -DA- AND SAVE

      DO 10 I=1,NCARTCOORD
         DA(I)=A(I,I)
  10  CONTINUE
      DO 12 I=1,NCARTCOORD
         DO 11 J=1,I
            A(J,I)=A(I,J)
  11     CONTINUE
         A(I,I)=DA(I)
  12  CONTINUE

! NOW DIAGONALIZE THE SYMETRIZED CARTESIAN FORCE CONSTANT MATRIX

      CALL EIGN(A,EIG,VEC,NCARTCOORD,RHO)
      WRITE(125,101)
      WRITE(125,102)
      CALL EIGOUT(VEC,NCARTCOORD,125,EIG)

! CREATE FULL FORCE MATRIX FOR MASS-WEIGHTED TRANSFORMATION

      DO 20 I=1,NCARTCOORD
         DO 15 J=1,I
            A(I,J)=A(J,I)
  15     CONTINUE
         A(I,I)=DA(I)
  20  CONTINUE

! CONSTRUCT MASS WEIGHTED MATRIX

      DO 30 I=1,NCARTCOORD
         DO 29 J=1,NCARTCOORD
            FA(I,J)=DIM(I)*FA(I,J)*DIM(J)
!           A(I,J)=DIM(I)*A(I,J)*DIM(J)
  29     CONTINUE
  30  CONTINUE

! DIAGONALIZE MASS-WEIGHTED MATRIX AND CONVERT FREQUENCIES
! TO WAVE NUMBERS
      CALL EIGN(FA,EIG,VEC,NCARTCOORD,RHO)
!      CALL EIGN(A,EIG,VEC,NCARTCOORD,RHO)
! THE EIGENVALUES ARE GIVEN BY EIG=4.0*(PI**2)*(W**2)
! WHERE W ARE THE HARMONIC VIBRATIONAL FREQUENCIES
      DO 35 I=1,NCARTCOORD
! FINDING THE HARMONIC VIBRATIONAL FREQUENCIES IN CM-1
        IF (EIG(I).LT.0) THEN 
          EIG(I)=-DSQRT(DABS(TFACT*EIG(I)))
        ELSE 
          EIG(I)=DSQRT(DABS(TFACT*EIG(I)))
        END IF
  35  CONTINUE

! NORMALIZE THE EIGENVECTORS

      DO 39 I=1,NCARTCOORD
         DO 38 J=1,NCARTCOORD
            A(I,J)=VEC(I,J)
  38     CONTINUE
  39  CONTINUE
      DO 42 I=1,NCARTCOORD
         DUM=0.0D0
         DO 40 J=1,NCARTCOORD
            DUM=DUM+A(I,J)*A(I,J)
  40     CONTINUE
         DUM=DSQRT(DUM)
         DO 41 J=1,NCARTCOORD
            A(I,J)=A(I,J)/DUM
  41     CONTINUE
  42  CONTINUE

! SET PHASE OF EIGENVECTORS

      DO 45 I=1,NCARTCOORD
         SGN=1.D0
         DO J=1,NCARTCOORD
            IF (A(I,J).GE.1.D-3) THEN
               SGN=1.D0
               GOTO 43
            ELSEIF (A(I,J).LE.-1.D-3) THEN
               SGN=-1.D0
               GOTO 43
            ENDIF
         ENDDO
   43    CONTINUE
         DO J=1,NCARTCOORD
            A(I,J)=SGN*A(I,J)
         ENDDO
   45 CONTINUE

! MASS-WEIGHT THE EIGENVECTORS SO THAT FURTHER NORMAL MODE
! TRANSFORMATION IS PERFORMED FROM THE CARTESIAN COORDINATES

      DO 47 I=1,NCARTCOORD
         DO 46 J=1,NCARTCOORD
            A(I,J)=DIM(J)*A(I,J)
            VEC(I,J)=A(I,J)
   46    CONTINUE
   47 CONTINUE

      WRITE(125,103)
      CALL EIGOUT(A,NCARTCOORD,125,EIG)      
      WRITE(125,104)

      END

!###############################################################################
! EVALUATE THE FORCE CONSTANT MATRIX BY DIFFERENCING
! THE GRADIENT OF THE POTENTIAL ENERGY FUNCTION 
!###############################################################################

      SUBROUTINE FMTRX(Q,PDOT,FA,A,NCARTCOORD)
      IMPLICIT NONE
      INTEGER :: I, J, K, IND
      INTEGER, PARAMETER :: NATOMMAX=100, NPTS=2
      INTEGER :: NCARTCOORD
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: Q, PDOT, GZ, EIG
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX,NPTS) :: DG
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX,3*NATOMMAX) :: FA, A
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX,3*NATOMMAX) :: P, AUX1
      DOUBLE PRECISION, DIMENSION(NPTS) :: DIST
      DOUBLE PRECISION, PARAMETER :: HINC=0.0001D+00
      DOUBLE PRECISION :: DU, DL
!HINC=0.0005D+00
  100 FORMAT(10X,5(1H*),' CALCULATION OF FORCE CONSTANT  ',5(1H*))
  102 FORMAT(A)
  103 FORMAT(15X,' FORCE CONSTANT MATRIX', /)
  104 FORMAT(/,10X,'SYMMETRIZED CARTESIAN FORCE CONSTANT MATRIX'/)

      DIST(1)=HINC
      DIST(2)=-HINC

! INITIALIZING VECTORS

      WRITE(125,100)
      DO I=1,NCARTCOORD
        EIG(I)=0.0D0
        DO J=1,NCARTCOORD
          FA(I,J)=0.0D+00
        END DO 
      END DO 

      DO I=1,NCARTCOORD
        GZ(I)=PDOT(I)
      END DO

      DO I=1,NCARTCOORD
        DO J=1,NPTS
          Q(I)=Q(I)+DIST(J)

! DISPLACE COORDINATE BY +- HINC AND CALCULATE GRADIENT INTO DG

! USING CENTRAL FINITE DIFFERENCE 
!        CALL DVDQ(Q,GZ,NCARTCOORD)

! USING RICHARDSON'S EXTRAPOLATION FORMULA
          CALL DER(Q,GZ,NCARTCOORD)
          Q(I)=Q(I)-DIST(J)

          DO K=1,NCARTCOORD
            DG(K,J)=GZ(K)
          END DO

        END DO

! DEFINE EACH LINE OF THE CARTESIAN FORCE CONSTANT MATRIX

        DO K=1,NCARTCOORD 
          FA(K,I)=(DG(K,1)-DG(K,2))/(2.0D+00*HINC) 
        END DO
      END DO

      WRITE(125,102)
      WRITE(125,103)

      CALL EIGOUT(FA,NCARTCOORD,125,EIG)

      DO 25 I=1,NCARTCOORD
         DO 24 J=1,I
            DU=0.5D+00*(FA(I,J)+FA(J,I))
            DL=0.5D+00*(FA(I,J)-FA(J,I))
            A(I,J)=DU
            A(J,I)=DL
   24    CONTINUE
   25 CONTINUE
      DO 27 I=1,NCARTCOORD
         A(I,I)=FA(I,I)
   27 CONTINUE

! A IS THE SYMETRIZED CARTESIAN FORCE CONSTANT MATRIX

      WRITE(125,102)
      WRITE(125,104)
      CALL EIGOUT(A,NCARTCOORD,125,EIG)

! PROJECTION TEST BEFORE ENDING

! INICIALIZING VECTOR

      DO I=1,3*NATOMMAX
       DO J=1,3*NATOMMAX
         P(I,J)=0.0D+00 
         AUX1(I,J)=0.0D+00
       END DO 
      END DO

      IND=0

      CALL PROJECT(Q,P,IND,NCARTCOORD)

      CALL DDER(Q,FA,NCARTCOORD)

      CALL MATP(P,FA,AUX1,NCARTCOORD,NCARTCOORD,NCARTCOORD,3*NATOMMAX)

      CALL MATP(AUX1,P,FA,NCARTCOORD,NCARTCOORD,NCARTCOORD,3*NATOMMAX)

      END
      
!###############################################################################
! EVALUATE THE FORCE CONSTANT MATRIX BY DIFFERENCING
! THE GRADIENT OF THE POTENTIAL ENERGY FUNCTION 
! THIS APPROACH USES RICHARDSON'S EXTRAPOLATION FORMULA
!###############################################################################

      SUBROUTINE DDER(X,VAL,NTA)
      IMPLICIT NONE
      INTEGER :: I,J,NTA
      INTEGER, PARAMETER :: NA=100, NT=3*NA
      DOUBLE PRECISION, PARAMETER :: PREC=1.0D-4
      DOUBLE PRECISION, DIMENSION(NT) :: X
      DOUBLE PRECISION, DIMENSION(NT,NT) :: VAL
      DOUBLE PRECISION :: D2VL
      DO I=1,NTA
      DO J=I,NTA
      VAL(I,J)=D2VL(X,PREC,I,J,NTA)
      ENDDO
      ENDDO
      DO I=1,NTA
      DO J=I,NTA
      VAL(J,I)=VAL(I,J)
      ENDDO
      ENDDO
      RETURN
      END

      FUNCTION D2VL(X,STEP,I,J,NTA)
      IMPLICIT NONE
      INTEGER :: I,J,IN,NTA
      INTEGER, PARAMETER :: NA=100, NT=3*NA
      DOUBLE PRECISION, DIMENSION(NT) :: X,XD
      DOUBLE PRECISION :: STEP,FMIN1,FMAX1
      DOUBLE PRECISION :: FMIN2,FMAX2,DVL,D2VL      
      DO IN=1,NTA
      XD(IN)=X(IN)
      ENDDO
      XD(J)=X(J)-STEP
      FMIN1=DVL(XD,STEP,I,NTA)
      XD(J)=X(J)+STEP
      FMAX1=DVL(XD,STEP,I,NTA)
      XD(J)=X(J)-2*STEP
      FMIN2=DVL(XD,STEP,I,NTA)
      XD(J)=X(J)+2*STEP
      FMAX2=DVL(XD,STEP,I,NTA)
      D2VL=(FMIN2-8.0*FMIN1+8.0*FMAX1-FMAX2)/(12.D0*STEP)
      RETURN
      END

!###############################################################################
! NUMERICAL GRADIENT OF THE POTENTIAL ENERGY FUNCTION IN CARTESIAN COORDINATES
! THIS APPROACH USES CENTRAL FINITE DIFFERENCE
!###############################################################################

      SUBROUTINE DVDQ(Q,PDOT,NCARTCOORD)
      IMPLICIT NONE
      INTEGER :: I, J, K
      INTEGER, PARAMETER :: NATOMMAX=100, NPTS=2
      INTEGER :: NCARTCOORD
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: Q, PDOT
      DOUBLE PRECISION, DIMENSION(NPTS) :: DIST, E
      DOUBLE PRECISION, PARAMETER :: HINC=0.0001D+00
!HINC=0.0005D+00

      DIST(1)=HINC
      DIST(2)=-HINC

! INITIALIZING VECTORS

      DO I=1,NCARTCOORD
        PDOT(I)=0.0D+00
      END DO

      DO I=1,NCARTCOORD
        DO J=1,NPTS
          Q(I)=Q(I)+DIST(J)

! DISPLACE COORDINATE BY +- HINC AND CALCULATE THE ENERGY

          CALL FCNSA(NCARTCOORD,Q,E(J))
          Q(I)=Q(I)-DIST(J)
        END DO

! USING TWO POINT CENTRAL DIFFERENCE FORMULA

        PDOT(I)=(E(1)-E(2))/(2.0D+00*HINC)
      END DO

      END      
            
!###############################################################################
! NUMERICAL GRADIENT OF THE POTENTIAL ENERGY FUNCTION IN CARTESIAN COORDINATES
! THIS APPROACH USES RICHARDSON'S EXTRAPOLATION FORMULA
!###############################################################################

      SUBROUTINE DER(X,VAL,NTA)
      IMPLICIT NONE
      INTEGER :: NTA,I
      INTEGER, PARAMETER :: NA=100, NT=3*NA
      DOUBLE PRECISION :: PREC=1.0D-4
      DOUBLE PRECISION, DIMENSION(NT) :: X,VAL
      DOUBLE PRECISION :: DVL
      DO I=1,NTA
      VAL(I)=DVL(X,PREC,I,NTA)
      ENDDO
      RETURN
      END

      FUNCTION DVL(X,STEP,I,NTA)
      IMPLICIT NONE
      INTEGER :: NTA,I,IN
      INTEGER, PARAMETER :: NA=100, NT=3*NA
      DOUBLE PRECISION, DIMENSION(NT) :: X,XD
      DOUBLE PRECISION :: STEP,FMIN1,FMAX1
      DOUBLE PRECISION :: FMIN2,FMAX2,VL,DVL
      DO IN=1,NTA
      XD(IN)=X(IN)
      ENDDO
      XD(I)=X(I)-STEP
      FMIN1=VL(XD,NTA)
      XD(I)=X(I)+STEP
      FMAX1=VL(XD,NTA)
      XD(I)=X(I)-2*STEP
      FMIN2=VL(XD,NTA)
      XD(I)=X(I)+2*STEP
      FMAX2=VL(XD,NTA)
      DVL=(FMIN2-8.0*FMIN1+8.0*FMAX1-FMAX2)/(12.D0*STEP)
      RETURN
      END

!###############################################################################
! FUNCTION EVALUATION
!###############################################################################

      FUNCTION VL(X,NTA)
      IMPLICIT NONE
      INTEGER :: NTA,NR
      INTEGER, PARAMETER :: NA=100, NT=3*NA
      DOUBLE PRECISION, DIMENSION(NT) :: X,R
      DOUBLE PRECISION :: VL
      NR=NTA*(NTA/3-1)/6
      CALL FCNSA(NTA,X,VL)
      RETURN
      END

!###############################################################################
! WRITE EIGENVALUES AND EIGENVECTORS IN UNIT IP
!###############################################################################

      SUBROUTINE EIGOUT(A,N,IP,EIG)
      IMPLICIT NONE
      INTEGER :: I,J,N,IP,L,K
      INTEGER, PARAMETER :: NATOMMAX=100
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: EIG
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX,3*NATOMMAX) :: A
      DOUBLE PRECISION :: D
  101 FORMAT(//,2X,6(9X,I3))
  102 FORMAT(5X,6E12.4)
  103 FORMAT(I3,1X,6F12.6)
  104 FORMAT(/)
      DO 2 I=1,N
         DO 1 J=1,I
            D=A(I,J)
            A(I,J)=A(J,I)
            A(J,I)=D
    1    CONTINUE
    2 CONTINUE
      K=0
    3 L=K+1
      K=K+6
      IF (N.LT.K) K=N
      WRITE(IP,101)(I,I=L,K)
      WRITE(IP,104)
      WRITE(IP,102)(EIG(I),I=L,K)
      WRITE(IP,104)
      DO 4 I=1,N
         WRITE(IP,103)I,(A(I,J),J=L,K)
    4 CONTINUE
      IF(K.LT.N) GOTO 3
      DO 7 I=1,N
         DO 6 J=1,I
            D=A(I,J)
            A(I,J)=A(J,I)
            A(J,I)=D
    6    CONTINUE
    7 CONTINUE
      RETURN
      END

!###############################################################################
! DIAGONALIZE A MATRIX A, OF WHICH ONLY LOWER TRIANGLE IS USED
! AND DESTROYED, USING THE GIVENS-HOUSHOLDER ALGORITHM.
! EIGENVALUES ARE RETURNED IN ALGEBRAIC ASCENDING ORDER IN ARRAY
! EIG THE EIGENVECTORS ARE RETURNED IN VEC.   
!
! PARAMETERS PASSED                   
! RHO IS THE UPPER LIMIT FOR OFF-DIAGONAL     
! NN IS THE SIZE OF THE MATRIX              
!###############################################################################

      SUBROUTINE EIGN(A,EIG,VEC,NN,RHO)
      IMPLICIT NONE
      INTEGER :: I,J,K,NN,N,N1,N2,NR,I1,NPAS,L,M,NRR,NT,NP,NV,LV,ITEMP
      INTEGER, PARAMETER :: NATOMMAX=100
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: EIG
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: W,BETASQ
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: GAMMA,BETA
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX,3*NATOMMAX) :: A,VEC
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: P, Q
      INTEGER, DIMENSION(3*NATOMMAX) :: IPOSV, IVPOS, IORD
      DOUBLE PRECISION :: RHOSQ,RHO,SHIFT,S,SGN,D,B,SQRTS,TEMP,WTAW
      DOUBLE PRECISION :: QJ,WJ,COSA,PP,PPBS,PPBR,COSAP,SINA,SINA2,DIA
      DOUBLE PRECISION :: R1,R2,DIF,A2,G,R12,SUM,U
      RHOSQ=RHO*RHO
      N=NN
      IF (N.EQ.0) RETURN
      SHIFT = 0.0D0
      N1=N-1
      N2=N-2
      GAMMA(1)=A(1,1)
      IF (N2.NE.0) THEN 
         IF (N2.LT.0) GOTO 280
         DO 260 NR=1,N2
            B=A(NR+1,NR)
            S=0.0D0
            DO 130 I=NR,N2
               S=S+A(I+2,NR)**2
  130       CONTINUE
!
!         PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION
!
            A(NR+1,NR)=0.0D0
            IF (S.GT.0) THEN
               S=S+B*B
               SGN=+1.0D0
               IF (B.LT.0) SGN = -1.D0
               SQRTS=DSQRT(S)
               D=SGN/(SQRTS+SQRTS)
               TEMP=DSQRT(0.5D0+B*D)
               W(NR)=TEMP
               A(NR+1,NR)=TEMP
               D=D/TEMP
               B=-SGN*SQRTS
!
!         D IS FACTOR OF PROPORTIONALITY. COMPUTE AND SAVE W VECTOR
!
               DO 170 I=NR,N2
                  TEMP=D*A(I+2,NR)
                  W(I+1)=TEMP
                  A(I+2,NR)=TEMP
  170          CONTINUE
!
!         PREMULTIPLY VECTOR W BY MATRIX A TO OBTAIN P VECTOR.
!         SIMULTANEOUSLY ACCUMULATE DOT PRODUCT WP,(THE SCALAR K).
!
               WTAW = 0.0D0
               DO 210 I=NR,N1
                  SUM = 0.0D0
                  DO 180 J=NR,I
                     SUM=SUM+A(I+1,J+1)*W(J)
  180             CONTINUE
                  I1=I+1
                  IF ((N1-I1).GE.0) THEN
                     DO 200 J=I1,N1
                        SUM=SUM+A(J+1,I+1)*W(J)
  200                CONTINUE
                  ENDIF
                  P(I)=SUM
                  WTAW=WTAW+SUM*W(I)
  210          CONTINUE
!
!         P VECTOR AND SCALAR K NOW STORED. NEXT COMPUTE Q VECTOR
!         AND FORM PAP MATRIX.
!
               DO 220 I=NR,N1
                  Q(I)=P(I)-WTAW*W(I)
  220          CONTINUE
               DO 240 J=NR,N1
                  QJ=Q(J)
                  WJ=W(J)
                  DO 230 I=J,N1
                     A(I+1,J+1)=A(I+1,J+1)-2.*(W(I)*QJ+WJ*Q(I))
  230             CONTINUE
  240          CONTINUE
  250          BETA(NR)=B
            ENDIF
            BETASQ(NR)=B*B
            GAMMA(NR+1)=A(NR+1,NR+1)
  260    CONTINUE
      ENDIF
      B=A(N,N-1)
      BETA(N-1)=B
      BETASQ(N-1)=B*B
      GAMMA(N)=A(N,N)
  280 BETASQ(N)=0.
!
!         ADJOIN AN IDENTITY MATRIX TO BE POSTMULTIPLIED BY ROTATIONS
!
      DO 300 I=1,N
         DO 290 J=1,N
            VEC(I,J)=0.0D0
  290    CONTINUE
         VEC(I,I)=1.0D0
  300 CONTINUE
      M=N
      SUM=0.0D0
      NPAS=1
      GOTO 400
  310 SUM=SUM+SHIFT
      COSA=1.
      G=GAMMA(1)-SHIFT
      PP=G
      PPBS=PP*PP+BETASQ(1)
      PPBR=DSQRT(PPBS)
      DO 330 J=1,M
         COSAP=COSA
         IF (PPBS.EQ.0) THEN
            SINA = 0.0D0
            SINA2=0.0D0
            COSA=1.0D0
         ELSE
            SINA=BETA(J)/PPBR
            SINA2=BETASQ(J)/PPBS
            COSA=PP/PPBR
!
!         POSTMULTIPLY IDENTITY BY P-TRANSPOSE
!
            NT=J+NPAS
            IF (NT.GT.N) NT=N
            DO 320 I=1,NT
               TEMP=COSA*VEC(J,I)+SINA*VEC(J+1,I)
               VEC(J+1,I)=-SINA*VEC(J,I)+COSA*VEC(J+1,I)
               VEC(J,I)=TEMP
  320       CONTINUE
         ENDIF
         DIA=GAMMA(J+1)-SHIFT
         U=SINA2*(G+DIA)
         GAMMA(J)=G+U
         G=DIA-U
         PP=DIA*COSA-SINA*COSAP*BETA(J)
         IF (J.EQ.M) THEN
            BETA(J)=SINA*PP
            BETASQ(J)=SINA2*PP*PP
            GOTO 340
         ENDIF
         PPBS=PP*PP+BETASQ(J+1)
         PPBR=DSQRT(PPBS)
         BETA(J)=SINA*PPBR
         BETASQ(J)=SINA2*PPBS
  330 CONTINUE
  340 GAMMA(M+1)=G
!
!         TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT
!
      NPAS=NPAS+1
      IF (BETASQ(M).GT.RHOSQ) GOTO 410
  390 EIG(M+1)=GAMMA(M+1)+SUM
  400 BETA(M)=0.0D0
      BETASQ(M)=0.0D0
      M=M-1
      IF (M.EQ.0) GOTO 430
      IF (BETASQ(M).LE.RHOSQ) GOTO 390
  410 CONTINUE 
!
!         TAKE ROOT OF CORNER 2 BY 2 NEAREST TO LOWER DIAGONAL IN
!         VALUE AS ESTIMATE OF EIGENVALUE TO USE FOR SHIFT
!
      A2=GAMMA(M+1)
      R2=0.5D0*A2
      R1=0.5D0*GAMMA(M)
      R12=R1+R2
      DIF=R1-R2
      TEMP=DSQRT(DIF*DIF+BETASQ(M))
      R1=R12+TEMP
      R2=R12-TEMP
      DIF=DABS(A2-R1)-DABS(A2-R2)
      IF (DIF.GE.0) THEN
         SHIFT=R2
      ELSE
         SHIFT=R1
      ENDIF
      GOTO 310
  430 EIG(1)=GAMMA(1)+SUM
!
!         INITIALIZE AUXILIARY TABLES REQUIRED FOR
!         REARANGING THE VECTORS
!
      DO 440 J=1,N
         IPOSV(J)=J
         IVPOS(J)=J
         IORD(J) = J
  440 CONTINUE
!
!         USE A TRANSPOSITON SORT TO ORDER THE EIGENVALUES
!
      M=N
      GOTO 470
  450 DO 460 J=1,M
         IF (EIG(J).GT.EIG(J+1)) THEN
            TEMP=EIG(J)
            EIG(J)=EIG(J+1)
            EIG(J+1)=TEMP
            ITEMP=IORD(J)
            IORD(J)=IORD(J+1)
            IORD(J+1)=ITEMP
         ENDIF
  460 CONTINUE
  470 M=M-1
      IF (M.NE.0) GOTO 450
      IF (N1.NE.0) THEN
         DO 490 L=1,N1
            NV=IORD(L)
            NP=IPOSV(NV)
            IF (NP.NE.L) THEN
               LV=IVPOS(L)
               IVPOS(NP)=LV
               IPOSV(LV)=NP
               DO 480 I=1,N
                  TEMP=VEC(L,I)
                  VEC(L,I)=VEC(NP,I)
                  VEC(NP,I) = TEMP
  480          CONTINUE
            ENDIF
  490    CONTINUE
      ENDIF
!
!         BACK TRANSFORM THE VECTORS OF THE TRIPLE DIAGONAL MATRIX
!
      DO 550 NRR=1,N
         K=N1
  510    K=K-1
         IF (K.GT.0) THEN
            SUM = 0.0
            DO 520 I=K,N1
               SUM=SUM+VEC(NRR,I+1)*A(I+1,K)
  520       CONTINUE
            SUM=SUM+SUM
            DO 530 I=K,N1
               VEC(NRR,I+1)=VEC(NRR,I+1)-SUM*A(I+1,K)
  530       CONTINUE
            GOTO 510
         ENDIF
  550 CONTINUE
      RETURN
      END

!###############################################################################
! WRITE STRUCTURE TO A MOLDEN READABLE FILE
!###############################################################################

      SUBROUTINE WRITMOLDEN(NINDEX,X,EIG,VEC,NATOM,ATOM)
      IMPLICIT NONE
      INTEGER :: I, J, K
      INTEGER, PARAMETER :: NATOMMAX=100
      INTEGER :: NINDEX, NATOM, NCARTCOORD
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: X, EIG
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX,3*NATOMMAX) :: VEC
      CHARACTER(2), DIMENSION(NATOMMAX) :: ATOM

      NCARTCOORD=3*NATOM

      WRITE(124,*)'[Molden Format]                 ',NINDEX
      WRITE(124,*)'               '
      WRITE(124,*)'[FREQ]'

      DO K=1,NCARTCOORD
        WRITE(124,*)EIG(K)
      ENDDO

      WRITE(124,*)'               '
      WRITE(124,*)'[FR-COORD]'
      DO I=1,NATOM
        WRITE(124,*) ATOM(I),(X(3*(I-1)+J),J=1,3) 
      END DO

      WRITE(124,*)'               '
      WRITE(124,*)'[FR-NORM-COORD]'

      DO K=1,NCARTCOORD
        WRITE(124,*)'vibration',K
        WRITE(124,'(3F20.7)')(VEC(K,J),J=1,NCARTCOORD)
      ENDDO
      WRITE(124,*)'               '
      WRITE(124,*)'               '
      WRITE(124,*)'               '
      WRITE(124,*)'               '

      END

!###############################################################################
! MATRIX MULTIPLICATION 
!###############################################################################

      SUBROUTINE MATP(A,B,C,N,L,M,ND)
      IMPLICIT NONE
      INTEGER :: ND,N,L,M,I,J,K
      DOUBLE PRECISION, DIMENSION(ND,L) :: A
      DOUBLE PRECISION, DIMENSION(ND,M) :: B,C
      DO I=1,N
         DO J=1,M
            C(I,J)=0.0
            DO K=1,L
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END

!###############################################################################
! FROM THE NR PROGRAM  
!###############################################################################

      SUBROUTINE PROJECT(X,P,IND,NTA)
      IMPLICIT NONE
      INTEGER :: IND,NTA,I,J
      REAL :: RNT
      INTEGER, PARAMETER :: NA=100, NT=3*NA
      DOUBLE PRECISION, DIMENSION(NT) :: X
      DOUBLE PRECISION, DIMENSION(NT,NT) :: B,S,SINV,AUX1,BTR,P,U
      DO I=1,NTA
      DO J=1,6
      B(I,J)=0.0
      ENDDO
      ENDDO
      DO I=1,NTA
      DO J=1,NTA
      U(I,J)=0.0
      ENDDO
      U(I,I)=1.0
      ENDDO
      DO I=1,NTA,3
      B(I,  1)=1.0
      B(I+1,2)=1.0
      B(I+2,3)=1.0
      B(I+1,4)= X(I+2)
      B(I+2,4)=-X(I+1)
      B(I,  5)=-X(I+2)
      B(I+2,5)= X(I)      
      B(I,  6)= X(I+1)
      B(I+1,6)=-X(I)
      ENDDO
      RNT=REAL(NTA/3)
      DO I=1,6
      DO J=1,6
      S(I,J)=0.0
      ENDDO
      ENDDO
      S(1,1)=RNT
      S(2,2)=RNT
      S(3,3)=RNT
      DO I=1,NTA,3
      S(1,5)=S(1,5)-X(I+2)
      S(1,6)=S(1,6)+X(I+1)
      S(2,4)=S(2,4)+X(I+2)
      S(2,6)=S(2,6)-X(I)
      S(3,4)=S(3,4)-X(I+1)
      S(3,5)=S(3,5)+X(I)
      S(4,5)=S(4,5)-X(I)*X(I+1)
      S(4,6)=S(4,6)-X(I)*X(I+2)
      S(5,6)=S(5,6)-X(I+1)*X(I+2)
      S(4,4)=S(4,4)+X(I+1)**2+X(I+2)**2
      S(5,5)=S(5,5)+X(I)**2  +X(I+2)**2
      S(6,6)=S(6,6)+X(I)**2  +X(I+1)**2
      ENDDO
      S(6,1)=S(1,6)
      S(6,2)=S(2,6)
      S(6,3)=S(3,6)
      S(6,4)=S(4,6)
      S(6,5)=S(5,6)
      S(5,1)=S(1,5)
      S(5,2)=S(2,5)
      S(5,3)=S(3,5)
      S(5,4)=S(4,5)
      S(4,1)=S(1,4)
      S(4,2)=S(2,4)
      S(4,3)=S(3,4)
      S(3,1)=S(1,3)
      S(3,2)=S(2,3)
      S(2,1)=S(1,2)
      CALL MINVER(S,SINV,IND)
      IF(IND.EQ.100)RETURN
      CALL MATP(B,SINV,AUX1,NTA,6,6,NT)
      DO I=1,6
        DO J=1,NTA
          BTR(I,J)=B(J,I)
        ENDDO
      ENDDO
      CALL MATP(AUX1,BTR,P,NTA,6,NTA,NT)
      DO I=1,NTA
      DO J=1,NTA
      P(I,J)=U(I,J)-P(I,J)
      ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE MINVER(A,AINV,IND)
      IMPLICIT NONE
      INTEGER :: IND,I,J
      INTEGER, PARAMETER :: NA=100, NT=3*NA, ND=6
      DOUBLE PRECISION, DIMENSION(NT,NT+1) :: RMAT
      DOUBLE PRECISION, DIMENSION(NT) :: X
      DOUBLE PRECISION, DIMENSION(NT,NT) :: A,AINV
      DO I=1,ND
        DO J=1,ND
          RMAT(I,J)=A(I,J)
        ENDDO
      ENDDO
      DO I=1,ND
      DO J=1,ND
        RMAT(J,ND+1)=0.0
        ENDDO
        RMAT(I,ND+1)=1.0
        CALL SVD(RMAT,ND,X,IND)
        IF(IND.EQ.100)RETURN
        DO J=1,ND
        AINV(J,I)=X(J)
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE SVD(RMAT,N,A,IND)
      IMPLICIT NONE
      INTEGER :: N,I,J,IND
      INTEGER, PARAMETER :: NA=100, NT=3*NA
      INTEGER, PARAMETER :: NMAX=1000, MMAX=100
      DOUBLE PRECISION, PARAMETER :: TOL=1.0D-9
      DOUBLE PRECISION, DIMENSION(NT,NT+1) :: RMAT
      DOUBLE PRECISION, DIMENSION(NT) :: A,W
      DOUBLE PRECISION, DIMENSION(NT,NT) :: V,U
      DOUBLE PRECISION, DIMENSION(NMAX) :: B
      DOUBLE PRECISION :: THRESH,WMAX
      DO 12 I=1,N
        DO 11 J=1,N
          U(I,J)=RMAT(I,J)
11      CONTINUE
        B(I)=RMAT(I,N+1)
12    CONTINUE
      CALL SVDCMP(U,N,N,W,V,IND)
      IF(IND.EQ.100)RETURN
      WMAX=0.0
      DO 13 J=1,N
        IF(W(J).GT.WMAX)WMAX=W(J)
13    CONTINUE
      THRESH=TOL*WMAX
      DO 14 J=1,N
        IF(W(J).LT.THRESH)W(J)=0.0
14    CONTINUE
      CALL SVBKSB(U,W,V,N,N,B,A)
      RETURN
      END

      SUBROUTINE SVDCMP(A,M,N,W,V,IND)
      IMPLICIT NONE
      INTEGER :: N,M,IND,I,J,K,L,NM,ITS
      INTEGER, PARAMETER :: NA=100, NT=3*NA
      INTEGER, PARAMETER :: NMAX=1000
      DOUBLE PRECISION, DIMENSION(NT) :: W
      DOUBLE PRECISION, DIMENSION(NT,NT) :: A,V
      DOUBLE PRECISION, DIMENSION(NMAX) :: RV1
      DOUBLE PRECISION :: G,SCALE,ANORM,S,F,H,C,Y,Z,X
      G=0.0
      SCALE=0.0
      ANORM=0.0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0
        S=0.0
        SCALE=0.0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+ABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.0
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0
        S=0.0
        SCALE=0.0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.0
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.0) THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=0.0
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0
            V(J,I)=0.0
31        CONTINUE
        ENDIF
        V(I,I)=1.0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=0.0
33        CONTINUE
        ENDIF
        IF (G.NE.0.0) THEN
          G=1.0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.0
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              F=(S/A(I,I))*G
              DO 35 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=0.0
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+1.0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=0.0
          S=1.0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((ABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=SQRT(F*F+G*G)
              W(I)=H
              H=1.0/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.30)THEN
          IND=100
          RETURN
          ENDIF
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
          G=SQRT(F*F+1.0)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C=1.0
          S=1.0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=SQRT(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=SQRT(F*F+H*H)
            W(J)=Z
            IF (Z.NE.0.0) THEN
              Z=1.0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 NM=1,M
              Y=A(NM,J)
              Z=A(NM,I)
              A(NM,J)= (Y*C)+(Z*S)
              A(NM,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END

      SUBROUTINE SVBKSB(U,W,V,M,N,B,X)
      IMPLICIT NONE
      INTEGER :: I,J,N,M,JJ
      INTEGER, PARAMETER :: NMAX=1000
      INTEGER, PARAMETER :: NA=100,NT=3*NA
      DOUBLE PRECISION, DIMENSION(NT) :: W,X
      DOUBLE PRECISION, DIMENSION(NT,NT) :: U,V
      DOUBLE PRECISION, DIMENSION(NMAX) :: B,TMP
      DOUBLE PRECISION :: S
      DO 12 J=1,N
        S=0.
        IF(W(J).NE.0.)THEN
          DO 11 I=1,M
            S=S+U(I,J)*B(I)
11        CONTINUE
          S=S/W(J)
        ENDIF
        TMP(J)=S
12    CONTINUE
      DO 14 J=1,N
        S=0.
        DO 13 JJ=1,N
          S=S+V(J,JJ)*TMP(JJ)
13      CONTINUE
        X(J)=S
14    CONTINUE
      RETURN
      END

!###############################################################################      