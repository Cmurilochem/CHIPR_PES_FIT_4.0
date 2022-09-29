C###############################################################################
C C2H2 GROUND-STATE POTENTIAL IN CARTESIAN COORDINATES (A.U)
C
C C. M. R. ROCHA AND A. J. C. VARANDAS (TO BE SUBMITTED) 2020
C
C N IS 3X THE NUMBER OF ATOMS, IN THIS CASE 12
C
C X IS THE INITIAL CARTESIAN COORDINATES IN BOHR
C 
C THE CONVENSION IS
C
C X(1) - X(3)  : X, Y, Z FOR H1
C X(4) - X(6)  : X, Y, Z FOR H2
C X(7) - X(9)  : X, Y, Z FOR C1
C X(10)- X(12) : X, Y, Z FOR C2
C
C POT IS THE POTENTIAL IN HARTREES
C
C###############################################################################

      SUBROUTINE FCNC2H2(N,X,POT)
      IMPLICIT NONE
      INTEGER, PARAMETER :: NATOM=4
      INTEGER, PARAMETER :: INTDIST=(NATOM*(NATOM-1))/2
      INTEGER, PARAMETER :: NATOMMAX=100
      DOUBLE PRECISION, DIMENSION(INTDIST) :: W,DVDR
      INTEGER :: I, N
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: X, R
      DOUBLE PRECISION :: POT
      LOGICAL :: DER
      DER=.FALSE.
      CALL CART2INTER(X,R,N)
      DO I=1,INTDIST
        W(I)=R(I)
      END DO
      CALL POTC2H2(W,POT,DER,DVDR)
      RETURN
      END

C###############################################################################
C C2H2 GROUND-STATE POTENTIAL IN CARTESIAN COORDINATES (A.U) AND FIRST DERIVATIVES
C
C C. M. R. ROCHA AND A. J. C. VARANDAS (TO BE SUBMITTED) 2020
C
C N IS 3X THE NUMBER OF ATOMS, IN THIS CASE 12
C
C X IS THE INITIAL CARTESIAN COORDINATES IN BOHR
C 
C THE CONVENSION IS
C
C X(1) - X(3)  : X, Y, Z FOR H1
C X(4) - X(6)  : X, Y, Z FOR H2
C X(7) - X(9)  : X, Y, Z FOR C1
C X(10)- X(12) : X, Y, Z FOR C2
C
C F IS THE POTENTIAL IN HARTREES
C
C G IS THE GRADIENT
C
C FIRST DERIVATIVES ARE CALCULATED SEMI-NUMERICALLY 
C NUMERICALLY FOR 2-PLUS-3-BODIES AND
C ANALITICALLY FOR 4-BODIES
C 
C###############################################################################

      SUBROUTINE FCNC2H2G(N,X,F,G)
      IMPLICIT NONE
      INTEGER :: I, N
      INTEGER, PARAMETER :: NATOMMAX=100
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: X, G
      DOUBLE PRECISION :: POT, F
      CALL FCNC2H2(N,X,POT)
      F=POT
      DO I=1,N
       G(I)=0.0D+00
      END DO
      CALL CARTDERPES(X,G)
      RETURN
      END

C###############################################################################
C C2H2 GROUND-STATE POTENTIAL ENERGY SURFACE IN A.U
C
C C. M. R. ROCHA AND A. J. C. VARANDAS (TO BE SUBMITTED) 2020
C
C DOUBLET C2H POTENTIAL FROM JOSEPH AND VARANDAS, J. PHYS. CHEM. A 114, 2655-2664 (2010)
C
C SINGLET CH2 POTENTIAL FROM JOSEPH AND VARANDAS, J. PHYS. CHEM. A 113, 4175-4183 (2009)
C
C FOUR-BODY TERMS USING CHIPR AND CALIBRATED FROM 42538 CCSD(T)/CBS(t,q) ENERGIES
C
C R1 DEFINES H-H BOND DISTANCE
C R2, R3, R4, R5 DEFINE C-H
C R6 DEFINES C-C BOND DISTANCE
C
C R IS THE 6D VECTOR OF BOND DISTANCES IN BOHR
C POT IS THE POTENTIAL IN HARTREES
C
C SET DER=.TRUE. TO CALCULATE THE GRADIENT IN INTERNAL COORDINATES
C THE DERIVATIVES ARE KEPT INTO THE 6D VECTOR DVDR
C
C##############################################################################

      SUBROUTINE POTC2H2(R,POT,DER,DVDR)
      IMPLICIT NONE
      INTEGER :: I
      INTEGER, PARAMETER :: NX=4
      INTEGER, PARAMETER :: INTDIST=NX*(NX-1)/2      
      DOUBLE PRECISION, DIMENSION(INTDIST) :: R
      DOUBLE PRECISION :: POT,POT23,POT4
      LOGICAL :: DER
      DOUBLE PRECISION, DIMENSION(INTDIST) :: DVDR23BD
      DOUBLE PRECISION, DIMENSION(INTDIST) :: DVDR4BD
      DOUBLE PRECISION, DIMENSION(INTDIST) :: DVDR      
      
      CALL POT23C2H2(R,POT23)
C     CALL CHIPR_TETRA_C2H2(R,POT4,DER,DVDR4BD)
      
      POT=POT23
C+POT4
      
      IF (DER) THEN
        CALL INTDERIVE(R,DVDR23BD,INTDIST)
        DO I=1,INTDIST
          DVDR(I)=DVDR23BD(I)+DVDR4BD(I)
        END DO        
      ELSE        
        DO I=1,INTDIST
          DVDR(I)=1000
        END DO        
      END IF   
      
      RETURN
      END
      
C####################################################################################
C CALCULATES DERIVATIVES OF THE POTENTIAL WITH RESPECT TO CARTESIAN COORDINATES
C THIS USES THE CHAIN RULE
C####################################################################################

      SUBROUTINE CARTDERPES(X,DVDX)
      IMPLICIT NONE 
      INTEGER :: I
      INTEGER, PARAMETER :: NATOM=4     
      DOUBLE PRECISION, DIMENSION(3*NATOM) :: X
      DOUBLE PRECISION, DIMENSION(3*NATOM) :: DVDX,O
      DOUBLE PRECISION, DIMENSION(3*NATOM-6) :: Y,R,DVDR
      DOUBLE PRECISION :: V
C
C     This subprogram uses the chain rule to calculate the derivatives of the
C     energy with respect to the cartesian coordinates from the derivatives
C     with respect to the internal coordinates for a three-body system.
C     The convention assumed in this subprogram is as follows:
C
C     R(1) : R(A-B)
C     R(2) : R(A-C)
C     R(3) : R(A-D)
C     R(4) : R(B-C)
C     R(5) : R(B-D)
C     R(6) : R(C-D)
C     X(1)  - X(3)  : X, Y, Z for atom A
C     X(4)  - X(6)  : X, Y, Z for atom B
C     X(7)  - X(9)  : X, Y, Z for atom C
C     X(10) - X(12) : X, Y, Z for atom D
C
C
      CALL C2I(X,O,3*NATOM)

      R(1)=O(1)
      R(2)=O(2)
      R(3)=O(3)
      R(4)=O(4)
      R(5)=O(5)
      R(6)=O(6)      

C     CALL POTC2H2(R,V,.TRUE.,DVDR)

      DO I=1,6
        Y(I)=DVDR(I)/R(I)
      END DO

      DVDX( 1)=(X(1)-X(4))*Y(1)+(X(1)-X(7))*Y(2)+(X(1)-X(10))*Y(3) 
      DVDX( 2)=(X(2)-X(5))*Y(1)+(X(2)-X(8))*Y(2)+(X(2)-X(11))*Y(3) 
      DVDX( 3)=(X(3)-X(6))*Y(1)+(X(3)-X(9))*Y(2)+(X(3)-X(12))*Y(3) 
      DVDX( 4)=(X(4)-X(1))*Y(1)+(X(4)-X(7))*Y(4)+(X(4)-X(10))*Y(5) 
      DVDX( 5)=(X(5)-X(2))*Y(1)+(X(5)-X(8))*Y(4)+(X(5)-X(11))*Y(5)
      DVDX( 6)=(X(6)-X(3))*Y(1)+(X(6)-X(9))*Y(4)+(X(6)-X(12))*Y(5)
      DVDX( 7)=(X(7)-X(1))*Y(2)+(X(7)-X(4))*Y(4)+(X(7)-X(10))*Y(6)
      DVDX( 8)=(X(8)-X(2))*Y(2)+(X(8)-X(5))*Y(4)+(X(8)-X(11))*Y(6)
      DVDX( 9)=(X(9)-X(3))*Y(2)+(X(9)-X(6))*Y(4)+(X(9)-X(12))*Y(6)
      DVDX(10)=(X(10)-X(1))*Y(3)+(X(10)-X(4))*Y(5)+(X(10)-X(7))*Y(6)
      DVDX(11)=(X(11)-X(2))*Y(3)+(X(11)-X(5))*Y(5)+(X(11)-X(8))*Y(6)
      DVDX(12)=(X(12)-X(3))*Y(3)+(X(12)-X(6))*Y(5)+(X(12)-X(9))*Y(6)      
   
      RETURN
      END 
      
C####################################################################################
C CONVERT CARTEZIAN TO INTERPARTICLE DISTANCES (A.U)
C####################################################################################

      SUBROUTINE C2I(X,R,NTA)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NA=4,NT=3*NA)
      DIMENSION X(NT),R(NT)
      K=0
      DO I=1,NTA,3
       DO J=I,NTA,3
        IF(I.NE.J)THEN
         K=K+1
         R(K)=SQRT((X(I)-X(J))**2+(X(I+1)-X(J+1))**2+
     *   (X(I+2)-X(J+2))**2)
        ENDIF
       ENDDO
      ENDDO
      RETURN
      END      

C###############################################################################
C SUM OF TWO- PLUS THREE-BODY TERMS FOR GROUND-STATE C2H2
C###############################################################################      
      
      SUBROUTINE POT23C2H2(R,POT)
      IMPLICIT NONE
      INTEGER :: I
      INTEGER, PARAMETER :: NATOM=4
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0000000D+00
      INTEGER, PARAMETER :: INTDIST=NATOM*(NATOM-1)/2
      DOUBLE PRECISION, DIMENSION(INTDIST) :: R, W
      INTEGER :: NTWOBD, NTHREEBD, NFOURBD
      REAL :: BICO, RECU1, RECU2
      DOUBLE PRECISION :: POT2BDC2H21, POT2BDC2H22 
      DOUBLE PRECISION :: POT3BDC2H21, POT3BDC2H22
      DOUBLE PRECISION :: POT, POT4BDC2H2
      INTEGER, PARAMETER :: TWOBD=2 
      INTEGER, PARAMETER :: THREEBD=3
      INTEGER, PARAMETER :: FOURBD=4
      LOGICAL :: MTC
      INTEGER :: DI, TI, QI
      INTEGER, DIMENSION(TWOBD) :: DCOMB
      INTEGER, DIMENSION(THREEBD) :: TCOMB
      INTEGER, DIMENSION(FOURBD) :: QCOMB
      INTEGER, PARAMETER :: DIMAX=2415
      INTEGER, PARAMETER :: TIMAX=54740
      INTEGER, PARAMETER :: QIMAX=916895
      INTEGER, DIMENSION(DIMAX) :: DA1, DA2
      DOUBLE PRECISION, DIMENSION(DIMAX,DIMAX) :: IPD
      DOUBLE PRECISION, DIMENSION(DIMAX) :: WIPD
      DOUBLE PRECISION, DIMENSION(DIMAX) :: POT2BD1, POT2BD2
      DOUBLE PRECISION, DIMENSION(DIMAX) :: VEHF, VDC
      INTEGER, DIMENSION(TIMAX) :: TA1, TA2, TA3
      DOUBLE PRECISION, DIMENSION(3) :: WIPT, RPROV
      DOUBLE PRECISION, DIMENSION(TIMAX) :: POT3BD1, POT3BD2
      INTEGER, DIMENSION(QIMAX) :: QA1, QA2, QA3, QA4
      DOUBLE PRECISION, DIMENSION(6) :: WIPQ   
      DOUBLE PRECISION, DIMENSION(QIMAX) :: POT4BD
      INTEGER, PARAMETER :: NDIM=2
      DOUBLE PRECISION, DIMENSION(NDIM,NDIM) :: V, VEC
      DOUBLE PRECISION, DIMENSION(NDIM) :: EIGVAL
      DOUBLE PRECISION :: RHO,RANG23BD,X0,GAM,MARK
      DOUBLE PRECISION :: EMIN,EREF,DEREF,DE,ES

C NUMBER OF TWO-BODY TERMS
      NTWOBD=BICO(NATOM,TWOBD)
C NUMBER OF THREE-BODY TERMS  USING RECURRENCE RELATION
      RECU1=((NATOM-TWOBD)/(TWOBD+1.00))
      NTHREEBD=NINT(RECU1*NTWOBD)
C NUMBER OF FOUR-BODY TERMS USING RECURRENCE RELATION
      RECU2=((NATOM-THREEBD)/(THREEBD+1.00))
      NFOURBD=NINT(RECU2*NTHREEBD)

C      WRITE(*,*) NTWOBD, NTHREEBD, NFOURBD

C INICIALIZING VARIABLES

      POT=0.000D+00
      POT2BDC2H21=0.000D+00
      POT2BDC2H22=0.000D+00
      POT3BDC2H21=0.000D+00
      POT3BDC2H22=0.000D+00
      POT4BDC2H2=0.000D+00

C###############################################################################
C COMBINATORIAL ANALYSIS OF TWO-BODY TERMS
C###############################################################################

      MTC=.FALSE.

      DI=0

      DO I=1,TWOBD
        DCOMB(I)=0
      END DO

  100 CALL NEXKSB(NATOM,TWOBD,DCOMB,MTC)

      DI=DI+1

C EACH DIATOMIC INDEX IS SAVED IN THE VECTORS TWOBDATOM1 AND TWOBDATOM2

      DA1(DI)=DCOMB(1) 

      DA2(DI)=DCOMB(2)

      IF (MTC) GOTO 100

C      DO I=1,NTWOBD
C        WRITE(*,*) DA1(I), DA2(I)
C      END DO

C DEFINING EACH INTERPARTICLE DISTANCE
      DO I=1,NTWOBD 
        IPD(DA1(I),DA2(I))=R(I)
C JUST IN CASE
        IPD(DA2(I),DA1(I))=IPD(DA1(I),DA2(I))  
C        WRITE(*,*) IPD(DA1(I),DA2(I)), DA1(I), DA2(I)       
      END DO

C CALCULATING TWO-DODY ENERGY
      DO I=1,NTWOBD 
        WIPD(I)=IPD(DA1(I),DA2(I))
C       WRITE(*,*) WIPD(I), DA1(I), DA2(I)
C JUST IN CASE
        IF (WIPD(I).EQ.ZERO) THEN 
          PRINT*, "ERROR IN THE TWO-BODY INDEXES"
          STOP
        ELSE
          IF (DA1(I).EQ.1 .AND. DA2(I).EQ.2) THEN
C TWO-BODY TERMS FOR HH          
            CALL H2SIG(WIPD(I),POT2BD1(I))
C           WRITE(*,*) WIPD(I),POT2BD1(I)
C            
          ELSE IF (DA1(I).EQ.1 .AND. DA2(I).EQ.3  
     1    .OR. DA1(I).EQ.1 .AND. DA2(I).EQ.4  
     2    .OR. DA1(I).EQ.2 .AND. DA2(I).EQ.3
     3    .OR. DA1(I).EQ.2 .AND. DA2(I).EQ.4) THEN
C TWO-BODY TERMS FOR CH 
            CALL CHPI(WIPD(I),POT2BD1(I))
C           WRITE(*,*) WIPD(I),POT2BD1(I)
C            
          ELSE IF (DA1(I).EQ.3 .AND. DA2(I).EQ.4) THEN 
C TWO-BODY TERMS FOR CC  
            CALL C2SIG(WIPD(I),POT2BD1(I))
C           WRITE(*,*) WIPD(I),POT2BD1(I)
C            
          ELSE
            PRINT*, "ERROR IN THE TWO-BODY INDEXES"
            STOP
          END IF
C TOTAL TWO-BODY ENERGY 
C         WRITE(*,*) WIPD(I), DA1(I), DA2(I), POT2BD1(I)
          POT2BDC2H21=POT2BDC2H21+POT2BD1(I)
        END IF
      END DO

C###############################################################################
C COMBINATORIAL ANALYSIS OF THREE-BODY TERMS
C###############################################################################

      MTC=.FALSE.

      TI=0

      DO I=1,THREEBD
        TCOMB(I)=0
      END DO 

  200 CALL NEXKSB(NATOM,THREEBD,TCOMB,MTC)

      TI=TI+1

C EACH TRIATOMIC INDEX IS SAVED IN THE VECTORS TA1, TA2 AND TA3

      TA1(TI)=TCOMB(1) 

      TA2(TI)=TCOMB(2)

      TA3(TI)=TCOMB(3)

      IF (MTC) GOTO 200

C      DO I=1,NTHREEBD
C        WRITE(*,*) TA1(I), TA2(I), TA3(I) 
C      END DO

C CALCULATING THREE-DODY ENERGY
      DO I=1,NTHREEBD
        WIPT(1)=IPD(TA1(I),TA2(I))
        WIPT(2)=IPD(TA1(I),TA3(I))
        WIPT(3)=IPD(TA2(I),TA3(I))
C       WRITE(*,*) TA1(I), TA2(I), TA3(I)
C       WRITE(*,*) WIPT(1), WIPT(2), WIPT(3)
C JUST IN CASE
        IF  (WIPT(1).EQ.ZERO 
     1  .OR. WIPT(2).EQ.ZERO 
     2  .OR. WIPT(3).EQ.ZERO) THEN 
          PRINT*, "ERROR IN THE THREE-BODY INDEXES"
          STOP
        ELSE
          IF ( TA1(I).EQ.1 .AND. 
     1         TA2(I).EQ.2 .AND. 
     3         TA3(I).EQ.3
     4    .OR. TA1(I).EQ.1 .AND. 
     5         TA2(I).EQ.2 .AND. 
     6         TA3(I).EQ.4) THEN
C THREE-BODY TERMS FOR CH2
C ADJUSTING THE CORRECT ORDER OF THE DISTANCES INSIDE "THRBODYCH2"
            CALL THRBODYCH2(WIPT,POT3BD1(I))
C            WRITE(*,*) WIPT,POT3BD1(I)
C
          ELSE IF ( TA1(I).EQ.1 .AND. 
     1         TA2(I).EQ.3 .AND. 
     3         TA3(I).EQ.4
     4    .OR. TA1(I).EQ.2 .AND. 
     5         TA2(I).EQ.3 .AND. 
     6         TA3(I).EQ.4) THEN 
C THREE-BODY TERMS FOR C2H
C THE ADJUSTMENT IS MADE INSIDE "THRBODYC2H" 
            CALL THRBODYC2H(WIPT,POT3BD1(I))
C           WRITE(*,*) WIPT,POT3BD1(I)
C        
          ELSE 
            PRINT*, "ERROR IN THE THREE-BODY INDEXES"
            STOP
          END IF
C         WRITE(*,*) WIPT(1),WIPT(2),WIPT(3),POT3BD1(I)
          POT3BDC2H21=POT3BDC2H21+POT3BD1(I)
        END IF
      END DO

C###############################################################################
C COMBINATORIAL ANALYSIS OF FOUR-BODY TERMS 
C###############################################################################
C
C      MTC=.FALSE.
C
C      QI=0
C
C      DO I=1,FOURBD
C        QCOMB(I)=0
C      END DO 
C
C  300 CALL NEXKSB(NATOM,FOURBD,QCOMB,MTC)
C
C      QI=QI+1
C
CC EACH TRETRATOMIC INDEX IS SAVED IN THE VECTORS QA1, QA2, QA3, QA4
C
C      QA1(QI)=QCOMB(1) 
C
C      QA2(QI)=QCOMB(2)
C
C      QA3(QI)=QCOMB(3)
C
C      QA4(QI)=QCOMB(4)
C
C      IF (MTC) GOTO 300
C
CC CALCULATING FOUR-DODY ENERGY
C      DO I=1,NFOURBD
C        WIPQ(1)=IPD(QA1(I),QA2(I))
C        WIPQ(2)=IPD(QA1(I),QA3(I))
C        WIPQ(3)=IPD(QA1(I),QA4(I))
C        WIPQ(4)=IPD(QA2(I),QA3(I))
C        WIPQ(5)=IPD(QA2(I),QA4(I))
C        WIPQ(6)=IPD(QA3(I),QA4(I))
C        IF (WIPQ(1).EQ.ZERO 
C     1  .OR. WIPQ(2).EQ.ZERO 
C     2  .OR. WIPQ(3).EQ.ZERO 
C     3  .OR. WIPQ(4).EQ.ZERO
C     4  .OR. WIPQ(5).EQ.ZERO
C     5  .OR. WIPQ(6).EQ.ZERO) THEN 
C          PRINT*, "ERROR IN THE FOUR-BODY INDEXES"
C          STOP
C        ELSE
C           CALL CHIPR_TETRA_C2H2(WIPQ,POT4BD(I),DER,DVDR3BD)
C           POT4BDC2H2=POT4BDC2H2+POT4BD(I)
C        END IF 
C      END DO

C###############################################################################
C FINAL POTENTIAL
C###############################################################################

      POT=POT2BDC2H21+POT3BDC2H21

      RETURN
      END
      
C###############################################################################
C NUMERICAL GRADIENT OF THE SUM OF TWO PLUS THREE BODIES IN INTERNAL COORDINATES
C THIS APPROACH USES CENTRAL DIFFERENCE
C###############################################################################

      SUBROUTINE INTDERIVE(X,VAL,NTA)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (PREC=1.0D-4)
      DIMENSION X(NTA),VAL(NTA)
      DO I=1,NTA
        VAL(I)=DVDR(X,PREC,I,NTA)
      END DO
      RETURN
      END

      FUNCTION DVDR(X,STEP,I,NTA)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NTA),XD(NTA)
      DO J=1,NTA
        XD(J)=X(J)
      END DO
      XD(I)=X(I)-STEP
      FMIN1=VNUM(XD,NTA)
      XD(I)=X(I)+STEP
      FMAX1=VNUM(XD,NTA)
C      XD(I)=X(I)-2*STEP
C      FMIN2=VNUM(XD,NTA)
C      XD(I)=X(I)+2*STEP
C      FMAX2=VNUM(XD,NTA)
C      DVDR=(FMIN2-8.0*FMIN1+8.0*FMAX1-FMAX2)/(12.D0*STEP)
      DVDR=(FMAX1-FMIN1)/(2.D0*STEP)
      RETURN
      END

      FUNCTION VNUM(X,NTA)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NTA)
      CALL POT23C2H2(X,POT)
      VNUM=POT
      RETURN
      END      
           
C###############################################################################
C FOURD-BODY TERM FOR GROUND-STATE C2H2
C###############################################################################
   

C###############################################################################

      FUNCTION BICO(N,K)
      INTEGER k,n
      REAL bico
      REAL factln
      bico=nint(exp(factln(n)-factln(k)-factln(n-k)))
      END

C###############################################################################

      FUNCTION FACTLN(N)
      INTEGER n
      REAL factln
      REAL a(100),gammln
      SAVE a
      DATA a/100*-1./
C      if (n.lt.0) pause 'negative factorial in factln'
      if (n.le.99) then 
      if (a(n+1).lt.0.) a(n+1)=gammln(n+1.)
      factln=a(n+1)
      else
      factln=gammln(n+1.) 
      endif
      return
      END

C###############################################################################

      FUNCTION GAMMLN(XX)
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     1 24.01409824083091d0,-1.231739572450155d0,
     2 0.1208650973866179d-2,-0.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
      y=y+1.d0
      ser=ser+cof(j)/y
      enddo
      gammln=tmp+log(stp*ser/x)
      return
      END

C###############################################################################

      SUBROUTINE NEXKSB(N,K,A,MTC)
      INTEGER :: A(K),H 
      LOGICAL :: MTC
      IF (K .gt.0 ) GOTO 30
      DO 1 I=1,N
 1    A(I)=0
      MTC=.FALSE.
      RETURN
   30 IF (MTC) GOTO 40
   20 M2=0
      H=K
      GOTO 50
   40 IF (M2 .lt. N-H) H=0
      H=H+1
      M2=A(K+1-H)
   50 DO 51 J=1,H
   51 A(K+J-H)=M2+J
      MTC=A(1) .ne. N-K+1
      RETURN
      END SUBROUTINE NEXKSB
      
C###############################################################################
C FINDING THREE BODY ENERGIES FOR CH2 SINGLET 
C###############################################################################

      SUBROUTINE THRBODYCH2(W,THREEBDCH2)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3) :: R, W
      DOUBLE PRECISION :: POTTOT, THREEBDCH2
      DOUBLE PRECISION :: IFN, V2EHFCH2, CM, FR

      R(1)=W(1)
      R(2)=W(2)
      R(3)=W(3)

      CALL CH2JV(R,POTTOT)
            
      CM=SQRT(2.0D0*R(3)**2+2.0D0*R(2)**2-R(1)**2)*1.0D0/2.0D0

      IFN=50.000D+00

C JUST IN CASE, SET BY BRUTE FORE THE THREE-BODY AS ZERO BEYOND 50 BORHS
      IF (R(1).GE.IFN .OR. R(2).GE.IFN .OR. R(3).GE.IFN) THEN 
        THREEBDCH2=0.0D+00
      ELSE
        THREEBDCH2=POTTOT-(FR(R(1),CM)+V2EHFCH2(R))
C       WRITE(*,*) POTTOT,FR(R(1),CM),THREEBDCH2,V2EHFCH2(R)
      END IF

      END SUBROUTINE THRBODYCH2      

C###############################################################################
C FINDING THREE BODY ENERGIES FOR C2H SIGMA PLUS DOUBLET 
C###############################################################################

      SUBROUTINE THRBODYC2H(W,THREEBDC2H)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3) :: R, W
      DOUBLE PRECISION :: POTTOT, THREEBDC2H
      DOUBLE PRECISION :: IFN, V2EHF

      R(1)=W(3)
      R(2)=W(2)
      R(3)=W(1)

      CALL C2HJV(R,POTTOT)

      IFN=50.000D+00

C JUST IN CASE, SET BY BRUTE FORE THE THREE-BODY AS ZERO BEYOND 50 BORHS
      IF (R(1).GE.IFN .OR. R(2).GE.IFN .OR. R(3).GE.IFN) THEN 
        THREEBDC2H=0.0D+00
      ELSE
        THREEBDC2H=POTTOT-V2EHF(R)
      END IF

      END SUBROUTINE THRBODYC2H

C###############################################################################
C CH PI DOUBLET POTENTIAL ENERGY CURVE - JOSEPH AND VARANDAS, J. PHYS. CHEM. A 114, 2655-2664 (2010)
C###############################################################################

      SUBROUTINE CHPI(R,POT)
      IMPLICIT NONE
      DOUBLE PRECISION :: R, POT
      DOUBLE PRECISION :: VEHF_CH, VDC2
      POT=VEHF_CH(R,2)+VDC2(R,2)
      RETURN
      END 
      
C###############################################################################
C C2 SIGMA SINGLET POTENTIAL ENERGY CURVE - JOSEPH AND VARANDAS, J. PHYS. CHEM. A 114, 2655-2664 (2010)
C###############################################################################

      SUBROUTINE C2SIG(R,POT)
      IMPLICIT NONE
      DOUBLE PRECISION :: R, POT
      DOUBLE PRECISION :: VEHF_CC, VDC2
      POT=VEHF_CC(R,1)+VDC2(R,1)
      RETURN
      END
      
C###############################################################################
C H2 SIGMA SINGLET POTENTIAL ENERGY CURVE - JOSEPH AND VARANDAS, J. PHYS. CHEM. A 113, 4175-4183 (2009)
C############################################################################### 

      SUBROUTINE H2SIG(R,POT)
      IMPLICIT NONE
      DOUBLE PRECISION :: R, POT
      DOUBLE PRECISION :: VEHF_HH, VDC2CH2
      POT=VEHF_HH(R,1)+VDC2CH2(R,1)
      RETURN
      END

C####################################################################################
C C2H SIGMA PLUS DOUBLET POTENTIAL ENERGY SURFACE
C JOSEPH AND VARANDAS, J. PHYS. CHEM. A 114, 2655-2664 (2010)
C####################################################################################
C R1=CC IN A.U.
C R2=CH IN A.U.
C R3=CH IN A.U.
C POTENTIAL IN A.U. 
C####################################################################################

      SUBROUTINE C2HJV(R,POT)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3) :: R
      DOUBLE PRECISION :: POT, V2EHF, EHF, V3DC
      POT=V2EHF(R)+EHF(R)+V3DC(R)
      RETURN 
      END

C     SINGLE SHEETED POTENTIAL ENERGY SURFACE IS WRITTEN AS,

C     V(R)=V2EHF(2,BODY)+VDC(2BODY)+VEHF(3,BODY)+VDC(3,BODY)  

C####################################################################################
C CH2 SINGLET POTENTIAL ENERGY SURFACE
C JOSEPH AND VARANDAS, J. PHYS. CHEM. A 113, 4175-4183 (2009)
C####################################################################################
C R1=HH IN A.U.
C R2=CH IN A.U.
C R3=CH IN A.U.
C POTENTIAL IN A.U.  
C####################################################################################

      SUBROUTINE CH2JV(R,POT)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3) :: R
      DOUBLE PRECISION :: POT, V2EHFCH2, FR
      DOUBLE PRECISION :: EHFCH2, V3DCCH2, CM
      CM=SQRT(2.0D0*R(3)**2+2.0D0*R(2)**2-R(1)**2)*1.0D0/2.0D0
      POT=FR(R(1),CM)+V2EHFCH2(R)+EHFCH2(R)+V3DCCH2(R)
      RETURN 
      END 
      
C     SINGLE SHEETED POTENTIAL ENERGY SURFACE IS WRITTEN AS

C     V(R)=V(1D,1 BODY)*F(R)+VEHF(2,BODY)+VDC(2BODY)+VEHF(3,BODY)+VDC(3,BODY) 
C####################################################################################

       DOUBLE PRECISION FUNCTION V2EHF(R)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION R(3)

c       OPEN (UNIT=40,FILE='two_body_ehf.dat')
C **  two body potential (EHFACE2U model used)

       V2EHF=VEHF_CC(R(1),1)+VDC2(R(1),1)+
     & VEHF_CH(R(2),2)+VDC2(R(2),2)+
     & VEHF_CH(R(3),3)+VDC2(R(3),3)


       RETURN
       END

c------------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION VEHF_CH(R,I)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A1(3),A2(3),A3(3),DD(3),RE(3),G0(3),G1(3),G2(3)
      DIMENSION A4(3),A5(3),A6(3),A7(3),A8(3),A9(3)

      DATA DD(2),A1(2),A2(2),A3(2),A4(2),A5(2),A6(2),A7(2),A8(2),
     &   A9(2)/0.2569653D0,1.69707421D0,0.964370989D0,1.31064307D0,
     &   0.299246996D0,0.434916463D0,0.0860477559D0,
     &   0.0352759271D0,-0.0147471197D0,0.00412075906D0/
      
      DATA DD(3),A1(3),A2(3),A3(3),A4(3),A5(3),A6(3),A7(3),A8(3),
     &   A9(3)/0.2569653D0,1.69707421D0,0.964370989D0,1.31064307D0,
     &   0.299246996D0,0.434916463D0,0.0860477559D0,
     &   0.0352759271D0,-0.0147471197D0,0.00412075906D0/

c**********************************************************************

c     &   A9(3)/0.2569653D0,1.69707421D0,0.964370989D0,1.31064307D0,
c     &   0.299246996D0,0.434916463D0,0.0860477559D0,
c     &   0.0352759271D0,-0.0147471197D0,0.00412075906D0/


c cn using 1.54 for c2
c   DD=      -0.2569653
c  1.1630215  1.69366893  0.483900724  0.964370989  1.31064307  0.299246996
c  0.434916463  0.0860477559  0.0352759271 -0.0147471197  0.00412075906
c A1=  1.69707421

c    1.1630215  1.69366893  0.483900724
c***********************************************************************
      G1(2)= 1.69366893D0/1.1630215D0
      DATA G0(2),G2(2)/1.1630215D0,0.483900724D0/

      G1(3)= 1.69366893D0/1.1630215D0
      DATA G0(3),G2(3)/1.1630215D0,0.483900724D0/
   
      DATA RE(1),RE(2),RE(3)/2.358D0,2.120D0,2.120D0/
CC R(1) DONE

      RS=R-RE(I)

      GAMA=G0(I)*(1.0D0+G1(I)*TANH(G2(I)*RS))

      VEHF_CH =-DD(I)*(1.0D0+A1(I)*RS+A2(I)*RS**2+A3(I)*RS**3+
     &          A4(I)*RS**4+A5(I)*RS**5+A6(I)*RS**6+A7(I)*RS**7+
     &          A8(I)*RS**8+A9(I)*RS**9)*EXP(-GAMA*RS)/R

    
       RETURN
       END
C---------------------------------------------------------------------------
C VEHF_CC DONE
      DOUBLE PRECISION FUNCTION VEHF_CC(R,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A1(3),A2(3),A3(3),DD(3),RE(3),G0(3),G1(3),G2(3)
      DIMENSION A4(3),A5(3),A6(3),A7(3),A8(3),A9(3),A10(3),A11(3)

      DATA DD(1),A1(1),A2(1),A3(1),A4(1),A5(1),A6(1),A7(1),A8(1),A9(1)/
     & 0.5028416d0,2.14475014d0,2.70507527d0,4.88078109d0, 
     & 1.51711031d0,0.0317173349d0,0.211157188d0,  
     & 0.0551671144d0,-0.0664174164d0,0.013815779d0/

c cn using 1.54 for c2


c**********************************************************************
c  DD=      -0.5028416
c  1.65664499  1.16195315  1.80986503  2.70507527  4.88078109  1.51711031
c  0.0317173349  0.211157188  0.0551671144 -0.0664174164  0.013815779
c  A1=  2.14475014

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c      1.65664499  1.16195315  1.80986503 

	DATA G0(1),G2(1)/1.65664499D0,1.80986503D0/
	G1(1)= 1.16195315D0/1.65664499D0

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      DATA RE(1),RE(2),RE(3)/2.3580d0,2.120D0,2.120D0/
      RS=R-RE(I)

      GAMA=G0(I)*(1.0D0+G1(I)*TANH(G2(I)*RS))

      VEHF_CC =-DD(I)*(1.0D0+A1(I)*RS+A2(I)*RS**2+A3(I)*RS**3+
     &          A4(I)*RS**4+A5(I)*RS**5+A6(I)*RS**6+A7(I)*RS**7+
     &          A8(I)*RS**8+A9(I)*RS**9)*EXP(-GAMA*RS)/R
c       WRITE(10,*)R,VEHF_CC
       RETURN
       END
c----------------------------------------------
C VDC2 DONE
      FUNCTION VDC2(R,I)                                 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CA(9,9)
C ** dynamical correlation, diatom  calculation
c     C **  atom-atom dispersion coefficients for HH
C      C-C 40.0D0,1023.8095D0,29611.2625D0
      DATA CA(1,1),CA(1,2),CA(1,3)/40.0D0,
     &    962.303D0,26160.246D0/
c   cn using 1.54 for c2

C **  atom-atom dispersion coefficients for CH
      DATA CA(2,1),CA(2,2),CA(2,3)/16.1D0,351.81D0,8687.1D0/
C **  atom-atom dispersion coefficients for CH
      DATA CA(3,1),CA(3,2),CA(3,3)/16.1D0,351.81D0,8687.1D0/
      RIN1 =1.0D0/R
      RIN2 =RIN1*RIN1
      RIN4 =RIN2*RIN2
      RIN5 =RIN4*RIN1
      RIN6 =RIN2*RIN2*RIN2
      RIN8 =RIN6*RIN2
      RIN10=RIN8*RIN2

      VDC2=-CA(I,1)*CHI(6,R,I)*RIN6-
     &      CA(I,2)*CHI(8,R,I)*RIN8-CA(I,3)*CHI(10,R,I)*RIN10

c      WRITE(12,*)R,VDC2
      RETURN
      END
c-------------------------------------------------------------------
      FUNCTION CHI(N,R,I)
CC CHI DONE
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A0(16),B0(16),R0(3)
C **  damping functions
      DATA (A0(J), J=1,16)/16.3661D0,10.062652D0,7.5709297D0,
     &      6.1869941D0,5.29027135D0,4.6549655D0,
     &      4.17772554D0,3.8040565D0,3.50229909D0,3.2527089D0,
     &      3.04213207D0,2.86194363D0,2.70562629,2.56852131D0,
     &      2.44713191D0,2.33877835D0/

      DATA (B0(J), J=1,16)/15.623647D0,14.19721174D0,12.9010098D0,
     &      11.723150819D0,10.65283005D0,9.6802293D0,
     &      8.79642676D0,7.9933152D0,7.26352748D0,6.6003693D0,
     &      5.99775016D0,5.45015706D0,4.95255907D0,4.50039165D0,
     &      4.089507D0,3.71613603D0/

      DATA R0(1),R0(2),R0(3)/7.887088D0,7.4096D0,7.4096D0/


      RHO=(11.0D0+2.5D0*R0(I))/2.0D0
      RR =R/RHO
      AUX=(1.0D0-EXP(-A0(N)*RR-B0(N)*RR*RR))**N
      CHI=AUX


      RETURN
      END
c------------------------------------------------------------
      DOUBLE PRECISION FUNCTION V3DC(R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(3),RS(3),CT(3)
      s1=R(1)
      s2=R(2)
      s3=R(3)

      f1=0.5D0*(1.D0-DTANH(1.D0*(6.D0*s1-s2-s3)))
      f2=0.5D0*(1.D0-DTANH(1.D0*(6.D0*s2-s3-s1)))
      f3=0.5D0*(1.D0-DTANH(1.D0*(6.D0*s3-s1-s2)))

      CALL INT_JAC(R,RS,CT)

 
      v3dc=f1 * (CHC2(6,R(1),CT(1)) * CHID(6,RS(1),3) / RS(1)**6
     & +CHC2(8,R(1),CT(1)) * CHID(8,RS(1),3) / RS(1)**8
     & +CHC2(10,R(1),CT(1)) * CHID(10,RS(1),3) / RS(1)**10)
     & +f2 * (CCCH(6,R(2),CT(2)) * CHID(6,RS(2),3) / RS(2)**6
     & +CCCH(8,R(2),CT(2)) * CHID(8,RS(2),3) / RS(2)**8
     & +CCCH(10,R(2),CT(2)) * CHID(10,RS(2),2) / RS(2)**10)
     & +f3 * (CCCH(6,R(3),CT(3)) * CHID(6,RS(3),3) / RS(3)**6
     & +CCCH(8,R(3),CT(3)) * CHID(8,RS(3),3) / RS(3)**8
     & +CCCH(10,R(3),CT(3)) * CHID(10,RS(3),3) / RS(3)**10)

c      v3dc=1.0d0*v3dc
      v3dc=-1.d0*v3dc

      RETURN
      END      
c--------------------------------------------------------------------
       DOUBLE PRECISION FUNCTION CHC2(N,R,CT)
CC CHC2  DONE
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION D(6),RM(6),C1(6),C2(6),C3(6),C4(6),C5(6),CINF(6)
c---------------------------------------------------------------
      CND(X,I)=D(I)*(1.0D0+C1(I)*(X-RM(I))+C2(I)*(X-RM(I))*(X-RM(I))
     &         +C3(I)*(X-RM(I))*(X-RM(I))*(X-RM(I)))*EXP(-C1(I)
     &         *(X-RM(I))-C4(I)*(X-RM(I))*(X-RM(I))
     &         -C5(I)*(X-RM(I))*(X-RM(I))*(X-RM(I)))+CINF(I)
c--------------------------------------------------------------------
c COEFFICIENTS ARE NUMBERED IN THE FOLLOWING SEQUENCE 
C    C60,C62,C80,C82,C10,C84(I=1,2,3,4,5,6 RESPECTIVELY) IMPORTANT !!!!!!!!!
C RM AND CINF ALSO FOLLOW THE SAME ORDER----------------------------
c------------------------------------------------------------------------
                          
       DATA (RM(I),I=1,6)/4.30D0,4.30D0,4.2843D0,4.2920D0,4.2946D0,
     &       4.2745d0/ 
       DATA (D(I),I=1,6)/ 8.074D0,5.0004D0,252.8254D0,450.7644D0,
     &      3635.7732D0,34.0875D0/ 
			    	
       DATA (CINF(I),I=1,6)/32.2D0,0.0D0,23.844D0,0.0D0,20124.0689D0,
     &      0.00D0/

        DATA (C1(I),I=1,6)/0.964455687D0,0.516503154D0,0.94120995D0,
     &      0.652702369D0,0.857938221D0,0.922252894D0/ 

       DATA (C2(I),I=1,6)/0.210044883D0,0.0354992997D0,0.207425096D0,
     &      0.103480321D0,0.222173958D0,0.200867481D0/

       DATA (C3(I),I=1,6)/-6.88549889D-05,-0.0156719496D0,
     &     0.000658242731D0,-0.00681056517D0,0.00210617391D0,
     &     1.44541987D-05/

      DATA (C4(I),I=1,6)/0.152402039D0,0.25091714D0,0.158689556D0,
     &     0.275239256D0,0.623506258D0,0.173158778D0/

      DATA (C5(I),I=1,6)/0.00626892328,1.2779794D-08,0.01136466D0,
     &     6.34100044D-09,4.33722898D-09,0.0177756391D0/
c-----------------------------------------------------------------
      IF(N.EQ.6)THEN

      CHC2=0.0d0 
      CHC2=CND(R,1)+CND(R,2)*POLEG(CT,2)

      END IF

      IF(N.EQ.8)THEN
      CHC2=0.0d0
      CHC2=CND(R,3)+CND(R,4)*POLEG(CT,2)+CND(R,6)*POLEG(CT,3)
      

      END IF

      IF(N.EQ.10)THEN
      CHC2=0.0d0
      CHC2=CND(R,5)
      END IF

      RETURN
      END
       
c--------------------------------------------------------------------------

     
      DOUBLE PRECISION FUNCTION CCCH(N,R,CT)
CC CCCH DONE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D(6),RM(6),C1(6),C2(6),C3(6),C4(6),C5(6),CINF(6)

      CND(X,I)=D(I)*(1.0D0+C1(I)*(X-RM(I))+C2(I)*(X-RM(I))*(X-RM(I))
     &         +C3(I)*(X-RM(I))*(X-RM(I))*(X-RM(I)))*EXP(-C1(I)
     &         *(X-RM(I))-C4(I)*(X-RM(I))*(X-RM(I))
     &         -C5(I)*(X-RM(I))*(X-RM(I))*(X-RM(I)))+CINF(I)

c--------------------------------------------------------------------
c COEFFICIENTS ARE NUMBERED IN THE FOLLOWING SEQUENCE
C    C60,C62,C80,C82,C10,C84(I=1,2,3,4,5,6 RESPECTIVELY)
C RM AND CINF ALSO FOLLOW THE SAME ORDER----------------------------
c------------------------------------------------------------------------

                       
      DATA (D(I),I=1,6)/9.6384D0,9.7588D0,615.8112,1108.5820D0,
     &      34402.7812D0,100.7716D0/
                        
      DATA (RM(I),I=1,6)/3.9560D0,3.7384D0,3.7769D0,3.6684D0,
     &      3.6503D0,3.6596D0/

      DATA (CINF(I),I=1,6)/56.10D0,0.0D0,1313.967D0,0.0D0,
     &      40389.7743D0,0.0D0/

      DATA (C1(I),I=1,6)/1.40818493D0,0.499552981D0,1.13491655D0,
     &      0.880042107D0,1.08965922D0,0.646217246D0/ 
 
      DATA (C2(I),I=1,6)/0.571748431D0,0.00991842255D0,0.437419296D0,
     &      0.208765467D0,0.34084621D0,0.172728605D0/

      DATA (C3(I),I=1,6)/0.0811681724D0,0.000591839355D0,0.0572007246D0,
     &      0.0118822525D0,0.0212747502D0,0.0138379383D0/

      DATA (C4(I),I=1,6)/0.321210295D0,0.203513464D0,0.181240363D0,
     &      0.136378753D0,0.408891165D0,0.22475278D0/ 

      DATA (C5(I),I=1,6)/0.0147024311D0,1.25589347D-08,4.63612401D-09,
     &      2.30879391D-09,1.91742759D-08,1.35341175D-09/

c----------------------------------------------------------------------

      IF(N.EQ.6)THEN
      CCCH=0.0d0 
      CCCH=CND(R,1)+CND(R,2)*POLEG(CT,2)
      C60=CND(R,1)

      END IF

      IF(N.EQ.8)THEN
      CCCH=0.0d0
      CCCH=CND(R,3)+CND(R,4)*POLEG(CT,2)+CND(R,6)*POLEG(CT,3)
      C80=CND(R,3)

      END IF

      IF(N.EQ.10)THEN
      CCCH=0.0d0
      CCCH=CND(R,5)
      END IF

      RETURN
      END

c---------------------------------------------------------------------------
c---------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION CHID(N,RS,I)
CC CHID DONE
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A0(16),B0(16),R0(3)
C ** R0 A-BC  taken as (R0 AB + R0 AC)/2
      DATA R0(1),R0(2),R0(3)/7.887088D0,7.4096D0,7.4096D0/
      DATA (A0(J), J=1,16)/16.3661D0,10.062652D0,7.5709297D0,
     &      6.1869941D0,5.29027135D0,4.6549655D0,
     &      4.17772554D0,3.8040565D0,3.50229909D0,3.2527089D0,
     &      3.04213207D0,2.86194363D0,2.70562629,2.56852131D0,
     &      2.44713191D0,2.33877835D0/

      DATA (B0(J), J=1,16)/15.623647D0,14.19721174D0,12.9010098D0,
     &      11.723150819D0,10.65283005D0,9.6802293D0,
     &      8.79642676D0,7.9933152D0,7.26352748D0,6.6003693D0,
     &      5.99775016D0,5.45015706D0,4.95255907D0,4.50039165D0,
     &      4.089507D0,3.71613603D0/

      RHO=(11.0D0+2.5D0*R0(I))/2.0D0
      RR =RS/RHO
      AUX1=(1.0D0-EXP(-A0(N)*RR-B0(N)*RR*RR))**N
      CHID=AUX1

      RETURN
      END

c--------------------------------------------------------------------------

       DOUBLE PRECISION FUNCTION POLEG(RAD,N)
       
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION POL(3) 

       POL(1)=1.0d0
       POL(2)=0.5d0*(3.d0*RAD*RAD-1.0d0)
c       write(*,*)POL(2) 
       POL(3)=(1.0d0/8.0d0*(35.0d0*(RAD)**4
     & -30.0D0*(RAD)**2+3.0d0))

       POLEG=POL(N)

       RETURN
       END

C*********************************************************
C   ROUTINE FOR INTENAL-JACOBI TRANSFORMATION
C   W1,W2,W3 ARE THE ATOMIC MASSES
C*********************************************************
      SUBROUTINE INT_JAC(R,RS,CCITA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=3)
      DIMENSION R(N),W(N),RS(N),CCITA(N),RPR(N)
      EXTERNAL RADIO,ANGLE
c      DATA (W(I),I=1,N)/1.00783D0,1.00783D0,12.010D0/
      DATA (W(I),I=1,N)/12.01D0,12.01D0,1.0783D0/
         RPR(1)=R(1)
         RPR(2)=R(2)
         RPR(3)=R(3)
      CM=W(1)/(W(1)+W(2))
      RS(1)=RADIO(CM,RPR)
      CCITA(1)=ANGLE(CM,RPR)/(RS(1)*R(1))
      IF(CCITA(1).GT.1.0D0) CCITA(1)=1.0D0
      IF(CCITA(1).LT.-1.0D0) CCITA(1)=-1.0D0

      CALL PERMUTA(RPR)
      CM=W(2)/(W(2)+W(3))
      RS(2)=RADIO(CM,RPR)
      CCITA(2)=ANGLE(CM,RPR)/(RS(2)*R(2))
      IF(CCITA(2).GT.1.0D0) CCITA(2)=1.0D0
      IF(CCITA(2).LT.-1.0D0) CCITA(2)=-1.0D0

      CALL PERMUTA(RPR)
      CM=W(3)/(W(3)+W(1))
      RS(3)=RADIO(CM,RPR)
      CCITA(3)=ANGLE(CM,RPR)/(RS(3)*R(3))
      IF(CCITA(3).GT.1.0D0) CCITA(3)=1.0D0
      IF(CCITA(3).LT.-1.0D0) CCITA(3)=-1.0D0
      RETURN
      END

      DOUBLE PRECISION FUNCTION RADIO(CM,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=3,ZERO=1.D-12)
      DIMENSION R(N)
      RADIO=ZERO
      ARG=(CM*CM-CM)*R(1)*R(1)+
     &    (1.0D0-CM)*R(2)*R(2)+CM*R(3)*R(3)
      IF(ARG.GT.ZERO) RADIO=DSQRT(ARG)
      RETURN
      END

      DOUBLE PRECISION FUNCTION ANGLE(CM,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=3)
      DIMENSION R(N)
      ANGLE=(CM-0.5D0)*R(1)*R(1)-0.5D0*R(2)*R(2)+
     &          0.5D0 *R(3)*R(3)
      RETURN
      END

      SUBROUTINE PERMUTA(R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=3)
      DIMENSION R(N)
      RAUX=R(1)
      R(1)=R(2)
      R(2)=R(3)
      R(3)=RAUX
      RETURN
      END


c**********************************************************************************

       DOUBLE PRECISION FUNCTION EHF(R)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION R(3)
     
       EHF=V3EHF(R)
c       write(*,*)R(1),R(2),R(3),EHF
       RETURN
       END


C************************************************************************

c      subroutine polycal(ap,ma,ndat,g1,g2,g3,R1_REF,R2_REF,R3_REF,x,bb,c
c      * ,y,afunc,np)

 
      DOUBLE PRECISION FUNCTION V3EHF(R)
      implicit none
      INTEGER ma,j,i,np,k,kk,jj,ndat,NCOF
      DOUBLE PRECISION x,bb,c,g1(6),g2(6),g3(6)

      DOUBLE PRECISION R1_REF(6),R2_REF(6),R3_REF(6),mon(183)

      DOUBLE PRECISION R1,R2,R3,temp(6),yy,kk1,R(3),CM,a(183)

      DOUBLE PRECISION A11,A12,A13,A14,Q1(6),Q2(6),Q3(6)
      DOUBLE PRECISION S2A_SQ(6),S2B_SQ(6),S3_CU(6)

      DOUBLE PRECISION  afunc,afunc1,afunc2,afunc3,afunc4,afunc5,afunc6
      DOUBLE PRECISION  afunc7
c      DOUBLE PRECISION temp1,g11,R11_REF,g22,R22_REF,g33,R33_REF

c      DOUBLE PRECISION Q11,Q22,Q33,S2A_SQt,S2B_SQt,S3_CUt
c-------------------------------------------------------------------


      ma=183
      np=6

      R1=R(1)
      R2=R(2)
      R3=R(3)


c**************************************************************************

      R1_REF(1)=1.5d0
      R2_REF(1)=4.0d0
      R3_REF(1)=4.0d0

      g1(1)=0.90d0  
      g2(1)=1.25d0
      g3(1)=1.25d0 

c-------------------------------------------
      R1_REF(2)=2.0d0
      R2_REF(2)=2.5d0
      R3_REF(2)=2.5d0

      g1(2)=1.25d0    
      g2(2)=0.85d0 
      g3(2)=0.85d0 

c------------------------------------------------
      R1_REF(3)=2.5d0
      R2_REF(3)=3.0d0
      R3_REF(3)=3.0d0

      g1(3)=0.50d0 
      g2(3)=0.60d0
      g3(3)=0.60d0 

c---------------------------------------------

      R1_REF(4)=3.0d0
      R2_REF(4)=2.5d0
      R3_REF(4)=2.5d0

      g1(4)=1.00d0
      g2(4)=0.60d0
      g3(4)=0.60d0 
c---------------------------------------------      
       
      R1_REF(5)=3.5d0
      R2_REF(5)=2.0d0
      R3_REF(5)=2.0d0

      g1(5)=1.00d0
      g2(5)=0.75d0 
      g3(5)=0.75d0 


c---------------------------------------------

      R1_REF(6)=2.0d0
      R2_REF(6)=4.5d0
      R3_REF(6)=4.5d0

      g1(6)=0.85d0   
      g2(6)=0.65d0
      g3(6)=0.65d0   

c-----------------------------------------------------------------


      do kk=1,ma
      a(kk)=0.0d0
      mon(kk)=0.0d0
      end do


      do k=1,np
      temp(k)=0.0d0

      Q1(k)=0.0d0
      Q2(k)=0.0d0
      Q3(k)=0.0d0
      S2A_SQ(k)=0.0d0
      S2B_SQ(k)=0.0d0
      S3_CU(k)=0.0d0

      end do

	 a(1)  =  16.6864496d0
	 a(2)  =  1.03049245d0
	 a(3)  = -5.4035779d0
	 a(4)  =  7.10900933d0
	 a(5)  =  7.21784199d0
	 a(6)  = -9.81460949d0
	 a(7)  =  3.27795053d0
	 a(8)  =  0.270209816d0
	 a(9)  =  1.66480651d0
	 a(10)  =  0.56933257d0
	 a(11)  = -3.03336156d0
	 a(12)  = -1.19491786d0
	 a(13)  = -1.6929539d0
	 a(14)  =  0.760150263d0
	 a(15)  =  0.186926635d0
	 a(16)  =  0.584617124d0
	 a(17)  = -0.974403053d0
	 a(18)  = -1.88951351d0
	 a(19)  = -2.44390804d0
	 a(20)  =  0.0500109437d0
	 a(21)  =  0.650709688d0
	 a(22)  =  1.0790761d0
	 a(23)  =  0.0618708853d0
	 a(24)  =  0.308258412d0
	 a(25)  =  0.0293130911d0
	 a(26)  = -0.402147465d0
	 a(27)  =  0.0651701558d0
	 a(28)  = -0.311597479d0
	 a(29)  = -0.376327803d0
	 a(30)  = -0.376190065d0
	 a(31)  =  0.307382831d0
	 a(32)  =  0.0502479089d0
	 a(33)  = -0.00580226945d0
	 a(34)  =  0.0981776188d0
	 a(35)  =  49.6773241d0
	 a(36)  =  27.5239869d0
	 a(37)  =  40.8921062d0
	 a(38)  =  12.5091524d0
	 a(39)  =  11.5502905d0
	 a(40)  =  24.6050314d0
	 a(41)  = -6.26253413d0
	 a(42)  =  4.30810551d0
	 a(43)  =  14.6023677d0
	 a(44)  =  3.03484874d0
	 a(45)  =  18.0233737d0
	 a(46)  = -8.12013158d0
	 a(47)  =  15.6482259d0
	 a(48)  =  0.538975643d0
	 a(49)  =  0.928167123d0
	 a(50)  =  0.2415591d0
	 a(51)  = -0.468667624d0
	 a(52)  =  0.153578339d0
	 a(53)  = -0.571476775d0
	 a(54)  =  1.06253969d0
	 a(55)  = -0.224403657d0
	 a(56)  = -0.530005417d0
	 a(57)  = -0.212499916d0
	 a(58)  = -1.14618372d0
	 a(59)  =  0.965687855d0
	 a(60)  = -1.76982611d0
	 a(61)  =  0.677881851d0
	 a(62)  = -0.703118845d0
	 a(63)  =  1.80516261d0
	 a(64)  = -1.15587852d0
	 a(65)  = -1.08937122d0
	 a(66)  = -0.978501756d0
	 a(67)  =  1.62080348d0
	 a(68)  = -0.110600716d0
	 a(69)  =  47.3833717d0
	 a(70)  = -28.5591456d0
	 a(71)  =  41.3691335d0
	 a(72)  =  9.2141863d0
	 a(73)  = -6.07989523d0
	 a(74)  = -19.1107885d0
	 a(75)  =  11.9401871d0
	 a(76)  = -1.24864143d0
	 a(77)  =  3.04173829d0
	 a(78)  =  0.248387303d0
	 a(79)  =  6.49991418d0
	 a(80)  = -2.22509392d0
	 a(81)  =  1.14525585d0
	 a(82)  =  0.231150805d0
	 a(83)  = -0.501595507d0
	 a(84)  = -0.000324742263d0
	 a(85)  =  0.385466925d0
	 a(86)  = -1.3119784d0
	 a(87)  =  0.695059788d0
	 a(88)  = -0.00803483093d0
	 a(89)  =  0.0955011224d0
	 a(90)  =  0.211871637d0
	 a(91)  = -0.0246001451d0
	 a(92)  =  0.0930468099d0
	 a(93)  = -0.0423302442d0
	 a(94)  = -0.0515798751d0
	 a(95)  = -0.000621421818d0
	 a(96)  =  0.0449171635d0
	 a(97)  = -0.123763784d0
	 a(98)  = -0.028657903d0
	 a(99)  = -0.0311010212d0
	 a(100)  = -0.0148710084d0
	 a(101)  =  0.0240814673d0
	 a(102)  =  0.00506141999d0
	 a(103)  = -180.977706d0
	 a(104)  =  37.3464498d0
	 a(105)  = -29.5512823d0
	 a(106)  = -29.6522417d0
	 a(107)  = -21.9564201d0
	 a(108)  = -54.1506702d0
	 a(109)  =  7.34409057d0
	 a(110)  =  1.70095192d0
	 a(111)  = -2.38385847d0
	 a(112)  = -1.21644833d0
	 a(113)  =  5.22014777d0
	 a(114)  = -1.18242336d0
	 a(115)  = -5.04015721d0
	 a(116)  = -0.169088162d0
	 a(117)  =  2.12317895d0
	 a(118)  =  1.41924159d0
	 a(119)  =  4.03443932d0
	 a(120)  =  0.234047303d0
	 a(121)  = -4.26183016d0
	 a(122)  =  6.96901925d0
	 a(123)  =  2.32886194d0
	 a(124)  = -1.71991013d0
	 a(125)  = -0.551012829d0
	 a(126)  = -1.86779998d0
	 a(127)  = -0.405932061d0
	 a(128)  = -0.179005099d0
	 a(129)  =  0.0714017909d0
	 a(130)  = -1.77340354d0
	 a(131)  =  1.30228377d0
	 a(132)  = -2.21605717d0
	 a(133)  =  0.257181499d0
	 a(134)  =  0.490239434d0
	 a(135)  = -0.125214527d0
	 a(136)  = -0.0493525608d0
	 a(137)  =  39.1412756d0
	 a(138)  =  0.838360148d0
	 a(139)  = -21.3253729d0
	 a(140)  =  8.04310997d0
	 a(141)  =  8.00336394d0
	 a(142)  =  6.66654828d0
	 a(143)  = -1.76903742d0
	 a(144)  =  2.56213951d0
	 a(145)  =  0.369909687d0
	 a(146)  =  0.148542497d0
	 a(147)  = -2.53867335d0
	 a(148)  =  4.17834717d0
	 a(149)  = -3.45615688d0
	 a(150)  = -0.161041445d0
	 a(151)  = -0.649358861d0
	 a(152)  =  0.255526393d0
	 a(153)  = -1.19260804d0
	 a(154)  = -1.09004702d0
	 a(155)  =  1.68390774d0
	 a(156)  = -1.10471914d0
	 a(157)  = -0.413778686d0
	 a(158)  = -0.104693826d0
	 a(159)  =  0.43549487d0
	 a(160)  =  0.892597814d0
	 a(161)  =  0.541298161d0
	 a(162)  = -0.224682773d0
	 a(163)  =  0.19937264d0
	 a(164)  =  0.881267899d0
	 a(165)  = -0.234665553d0
	 a(166)  =  1.12702073d0
	 a(167)  =  0.416016149d0
	 a(168)  = -0.550898813d0
	 a(169)  = -0.0282879081d0
	 a(170)  = -0.136164453d0
	 a(171)  =  1.04180918d0
	 a(172)  =  1.54261556d0
	 a(173)  =  1.09517131d0
	 a(174)  = -0.456960478d0
	 a(175)  =  0.268479403d0
	 a(176)  = -0.52597555d0
	 a(177)  =  0.112851431d0
	 a(178)  =  0.324481057d0
	 a(179)  =  0.156415487d0
	 a(180)  = -0.0141176961d0
	 a(181)  =  0.490278034d0
	 a(182)  = -0.08397785d0
	 a(183)  =  0.102457566d0

      temp(1)=temp(1)+(1.0d0-tanh(g1(1)*(R1-R1_REF(1))))*
     *  (1.0d0-tanh(g2(1)*(R2-R2_REF(1))))*
     *  (1.0d0-tanh(g3(1)*(R3-R3_REF(1))))
c-----------------------------------------------------------------
      A11=SQRT(1.0D0/3.0D0)
      A12=SQRT(1.0D0/2.0D0)
      A13=SQRT(2.0D0/3.0D0)
      A14=SQRT(1.0D0/6.0D0)

      Q1(1)=(A11*(R1-R1_REF(1)))+(A11*(R2-R2_REF(1)))+
     * (A11*(R3-R3_REF(1)))
      Q2(1)=(A12*(R2-R2_REF(1)))-(A12*(R3-R3_REF(1)))      
      Q3(1)=(A13*(R1-R1_REF(1)))-(A14*(R2-R2_REF(1)))-
     * (A14*(R3-R3_REF(1)))

      S2A_SQ(1)=(Q2(1)*Q2(1))+(Q3(1)*Q3(1))
      S2B_SQ(1)=(Q2(1)*Q2(1))-(Q3(1)*Q3(1))
      S3_CU(1)=(Q3(1)*Q3(1)*Q3(1))-(3.0D0*Q2(1)*Q2(1)*Q3(1))


c-----------------------------------------------------------------

ctemp1,g11,R11_REF,g22,R22_REF,g33,R33_REF,Q11,Q22,Q33,S2A_SQt,S2B_SQt,S3_CUt

c---------------------------------------------------------------------
     
      mon(1)=1.0d0
      mon(2)=Q1(1)
      mon(3)=Q3(1)
      mon(4)=Q1(1)*Q1(1) 
      mon(5)=S2A_SQ(1)
      mon(6)=Q1(1)*Q3(1)  
      mon(7)=S2B_SQ(1)
      mon(8)=Q1(1)*Q1(1)*Q1(1)
      mon(9)=Q1(1)*S2A_SQ(1)
      mon(10)=S3_CU(1)
      mon(11)=Q1(1)*Q1(1)*Q3(1)
      mon(12)=Q1(1)*S2B_SQ(1)
      mon(13)=Q3(1)*S2A_SQ(1)
      mon(14)=Q1(1)*Q1(1)*Q1(1)*Q1(1)
      mon(15)=Q1(1)*Q1(1)*S2A_SQ(1)
      mon(16)=S2A_SQ(1)*S2A_SQ(1)
      mon(17)=Q1(1)*S3_CU(1)
      mon(18)=Q1(1)*Q1(1)*Q1(1)*Q3(1)
      mon(19)=Q1(1)*Q1(1)*S2B_SQ(1)
      mon(20)=Q1(1)*Q3(1)*S2A_SQ(1)
      mon(21)=Q3(1)*S3_CU(1)
      mon(22)=S2A_SQ(1)*S2B_SQ(1)
      mon(23)=Q1(1)**5
      mon(24)=Q1(1)*Q1(1)*Q1(1)*S2A_SQ(1)
      mon(25)=Q1(1)*S2A_SQ(1)*S2A_SQ(1)
      mon(26)=Q1(1)*Q1(1)*S3_CU(1)
      mon(27)=S2A_SQ(1)*S3_CU(1)
      mon(28)=Q1(1)*Q1(1)*Q1(1)*Q1(1)*Q3(1)
      mon(29)=Q1(1)*Q1(1)*Q1(1)*S2B_SQ(1)
      mon(30)=Q1(1)*Q1(1)*Q3(1)*S2A_SQ(1)
      mon(31)=Q1(1)*Q3(1)*S3_CU(1)
      mon(32)=Q1(1)*S2A_SQ(1)*S2B_SQ(1)
      mon(33)=Q3(1)*S2A_SQ(1)*S2A_SQ(1)
      mon(34)=S2B_SQ(1)*S3_CU(1)
c-----------------------------------------------------
      afunc1=0.0d0
      DO j=1,34

      afunc1=afunc1+a(j)*mon(j)

      END DO
      afunc1=afunc1*temp(1)
c     afunc=afunc1	

c**************************************************************************

      temp(2)=temp(2)+(1.0d0-tanh(g1(2)*(R1-R1_REF(2))))*
     *  (1.0d0-tanh(g2(2)*(R2-R2_REF(2))))*
     *  (1.0d0-tanh(g3(2)*(R3-R3_REF(2))))
c-----------------------------------------------------------------
      A11=SQRT(1.0D0/3.0D0)
      A12=SQRT(1.0D0/2.0D0)
      A13=SQRT(2.0D0/3.0D0)
      A14=SQRT(1.0D0/6.0D0)

      Q1(2)=(A11*(R1-R1_REF(2)))+(A11*(R2-R2_REF(2)))+
     * (A11*(R3-R3_REF(2)))
      Q2(2)=(A12*(R2-R2_REF(2)))-(A12*(R3-R3_REF(2)))      
      Q3(2)=(A13*(R1-R1_REF(2)))-(A14*(R2-R2_REF(2)))-
     * (A14*(R3-R3_REF(2)))

      S2A_SQ(2)=(Q2(2)*Q2(2))+(Q3(2)*Q3(2))
      S2B_SQ(2)=(Q2(2)*Q2(2))-(Q3(2)*Q3(2))
      S3_CU(2)=(Q3(2)*Q3(2)*Q3(2))-(3.0D0*Q2(2)*Q2(2)*Q3(2))

c---------------------------------------------------------------------
     
      mon(35)=1.0d0
      mon(36)=Q1(2)
      mon(37)=Q3(2)
      mon(38)=Q1(2)*Q1(2) 
      mon(39)=S2A_SQ(2)
      mon(40)=Q1(2)*Q3(2)  
      mon(41)=S2B_SQ(2)
      mon(42)=Q1(2)*Q1(2)*Q1(2)
      mon(43)=Q1(2)*S2A_SQ(2)
      mon(44)=S3_CU(2)
      mon(45)=Q1(2)*Q1(2)*Q3(2)
      mon(46)=Q1(2)*S2B_SQ(2)
      mon(47)=Q3(2)*S2A_SQ(2)
      mon(48)=Q1(2)*Q1(2)*Q1(2)*Q1(2)
      mon(49)=Q1(2)*Q1(2)*S2A_SQ(2)
      mon(50)=S2A_SQ(2)*S2A_SQ(2)
      mon(51)=Q1(2)*S3_CU(2)
      mon(52)=Q1(2)*Q1(2)*Q1(2)*Q3(2)
      mon(53)=Q1(2)*Q1(2)*S2B_SQ(2)
      mon(54)=Q1(2)*Q3(2)*S2A_SQ(2)
      mon(55)=Q3(2)*S3_CU(2)
      mon(56)=S2A_SQ(2)*S2B_SQ(2)
      mon(57)=Q1(2)**5
      mon(58)=Q1(2)*Q1(2)*Q1(2)*S2A_SQ(2)
      mon(59)=Q1(2)*S2A_SQ(2)*S2A_SQ(2)
      mon(60)=Q1(2)*Q1(2)*S3_CU(2)
      mon(61)=S2A_SQ(2)*S3_CU(2)
      mon(62)=Q1(2)*Q1(2)*Q1(2)*Q1(2)*Q3(2)
      mon(63)=Q1(2)*Q1(2)*Q1(2)*S2B_SQ(2)
      mon(64)=Q1(2)*Q1(2)*Q3(2)*S2A_SQ(2)
      mon(65)=Q1(2)*Q3(2)*S3_CU(2)
      mon(66)=Q1(2)*S2A_SQ(2)*S2B_SQ(2)
      mon(67)=Q3(2)*S2A_SQ(2)*S2A_SQ(2)
      mon(68)=S2B_SQ(2)*S3_CU(2)

      afunc2=0.0d0 
      DO j=35,68

      afunc2=afunc2+a(j)*mon(j)

      END DO
      afunc2=afunc2*temp(2)



      afunc=afunc1+afunc2
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
      temp(3)=temp(3)+(1.0d0-tanh(g1(3)*(R1-R1_REF(3))))*
     *  (1.0d0-tanh(g2(3)*(R2-R2_REF(3))))*
     *  (1.0d0-tanh(g3(3)*(R3-R3_REF(3))))
c-----------------------------------------------------------------
      A11=SQRT(1.0D0/3.0D0)
      A12=SQRT(1.0D0/2.0D0)
      A13=SQRT(2.0D0/3.0D0)
      A14=SQRT(1.0D0/6.0D0)

      Q1(3)=(A11*(R1-R1_REF(3)))+(A11*(R2-R2_REF(3)))+
     * (A11*(R3-R3_REF(3)))
      Q2(3)=(A12*(R2-R2_REF(3)))-(A12*(R3-R3_REF(3)))      
      Q3(3)=(A13*(R1-R1_REF(3)))-(A14*(R2-R2_REF(3)))-
     * (A14*(R3-R3_REF(3)))

      S2A_SQ(3)=(Q2(3)*Q2(3))+(Q3(3)*Q3(3))
      S2B_SQ(3)=(Q2(3)*Q2(3))-(Q3(3)*Q3(3))
      S3_CU(3)=(Q3(3)*Q3(3)*Q3(3))-(3.0D0*Q2(3)*Q2(3)*Q3(3))

c---------------------------------------------------------------------
      mon(69)=1.0d0
      mon(70)=Q1(3)
      mon(71)=Q3(3)
      mon(72)=Q1(3)*Q1(3) 
      mon(73)=S2A_SQ(3)
      mon(74)=Q1(3)*Q3(3)  
      mon(75)=S2B_SQ(3)
      mon(76)=Q1(3)*Q1(3)*Q1(3)
      mon(77)=Q1(3)*S2A_SQ(3)
      mon(78)=S3_CU(3)
      mon(79)=Q1(3)*Q1(3)*Q3(3)
      mon(80)=Q1(3)*S2B_SQ(3)
      mon(81)=Q3(3)*S2A_SQ(3)
      mon(82)=Q1(3)*Q1(3)*Q1(3)*Q1(3)
      mon(83)=Q1(3)*Q1(3)*S2A_SQ(3)
      mon(84)=S2A_SQ(3)*S2A_SQ(3)
      mon(85)=Q1(3)*S3_CU(3)
      mon(86)=Q1(3)*Q1(3)*Q1(3)*Q3(3)
      mon(87)=Q1(3)*Q1(3)*S2B_SQ(3)
      mon(88)=Q1(3)*Q3(3)*S2A_SQ(3)
      mon(89)=Q3(3)*S3_CU(3)
      mon(90)=S2A_SQ(3)*S2B_SQ(3)
      mon(91)=Q1(3)**5
      mon(92)=Q1(3)*Q1(3)*Q1(3)*S2A_SQ(3)
      mon(93)=Q1(3)*S2A_SQ(3)*S2A_SQ(3)
      mon(94)=Q1(3)*Q1(3)*S3_CU(3)
      mon(95)=S2A_SQ(3)*S3_CU(3)
      mon(96)=Q1(3)*Q1(3)*Q1(3)*Q1(3)*Q3(3)
      mon(97)=Q1(3)*Q1(3)*Q1(3)*S2B_SQ(3)
      mon(98)=Q1(3)*Q1(3)*Q3(3)*S2A_SQ(3)
      mon(99)=Q1(3)*Q3(3)*S3_CU(3)
      mon(100)=Q1(3)*S2A_SQ(3)*S2B_SQ(3)
      mon(101)=Q3(3)*S2A_SQ(3)*S2A_SQ(3)
      mon(102)=S2B_SQ(3)*S3_CU(3)
 
  
  
      afunc3=0.0d0 
      DO j=69,102

      afunc3=afunc3+a(j)*mon(j)

      END DO
      afunc3=afunc3*temp(3)


      
      
      temp(4)=temp(4)+(1.0d0-tanh(g1(4)*(R1-R1_REF(4))))*
     *  (1.0d0-tanh(g2(4)*(R2-R2_REF(4))))*
     *  (1.0d0-tanh(g3(4)*(R3-R3_REF(4))))
c-----------------------------------------------------------------
      A11=SQRT(1.0D0/3.0D0)
      A12=SQRT(1.0D0/2.0D0)
      A13=SQRT(2.0D0/3.0D0)
      A14=SQRT(1.0D0/6.0D0)

      Q1(4)=(A11*(R1-R1_REF(4)))+(A11*(R2-R2_REF(4)))+
     * (A11*(R3-R3_REF(4)))
      Q2(4)=(A12*(R2-R2_REF(4)))-(A12*(R3-R3_REF(4)))      
      Q3(4)=(A13*(R1-R1_REF(4)))-(A14*(R2-R2_REF(4)))-
     * (A14*(R3-R3_REF(4)))

      S2A_SQ(4)=(Q2(4)*Q2(4))+(Q3(4)*Q3(4))
      S2B_SQ(4)=(Q2(4)*Q2(4))-(Q3(4)*Q3(4))
      S3_CU(4)=(Q3(4)*Q3(4)*Q3(4))-(3.0D0*Q2(4)*Q2(4)*Q3(4))

c---------------------------------------------------------------------
     
      mon(103)=1.0d0
      mon(104)=Q1(4)
      mon(105)=Q3(4)
      mon(106)=Q1(4)*Q1(4) 
      mon(107)=S2A_SQ(4)
      mon(108)=Q1(4)*Q3(4)  
      mon(109)=S2B_SQ(4)
      mon(110)=Q1(4)*Q1(4)*Q1(4)
      mon(111)=Q1(4)*S2A_SQ(4)
      mon(112)=S3_CU(4)
      mon(113)=Q1(4)*Q1(4)*Q3(4)
      mon(114)=Q1(4)*S2B_SQ(4)
      mon(115)=Q3(4)*S2A_SQ(4)
      mon(116)=Q1(4)*Q1(4)*Q1(4)*Q1(4)
      mon(117)=Q1(4)*Q1(4)*S2A_SQ(4)
      mon(118)=S2A_SQ(4)*S2A_SQ(4)
      mon(119)=Q1(4)*S3_CU(4)
      mon(120)=Q1(4)*Q1(4)*Q1(4)*Q3(4)
      mon(121)=Q1(4)*Q1(4)*S2B_SQ(4)
      mon(122)=Q1(4)*Q3(4)*S2A_SQ(4)
      mon(123)=Q3(4)*S3_CU(4)
      mon(124)=S2A_SQ(4)*S2B_SQ(4)
      mon(125)=Q1(4)**5
      mon(126)=Q1(4)*Q1(4)*Q1(4)*S2A_SQ(4)
      mon(127)=Q1(4)*S2A_SQ(4)*S2A_SQ(4)
      mon(128)=Q1(4)*Q1(4)*S3_CU(4)
      mon(129)=S2A_SQ(4)*S3_CU(4)
      mon(130)=Q1(4)*Q1(4)*Q1(4)*Q1(4)*Q3(4)
      mon(131)=Q1(4)*Q1(4)*Q1(4)*S2B_SQ(4)
      mon(132)=Q1(4)*Q1(4)*Q3(4)*S2A_SQ(4)
      mon(133)=Q1(4)*Q3(4)*S3_CU(4)
      mon(134)=Q1(4)*S2A_SQ(4)*S2B_SQ(4)
      mon(135)=Q3(4)*S2A_SQ(4)*S2A_SQ(4)
      mon(136)=S2B_SQ(4)*S3_CU(4)
 
  
      afunc4=0.0d0 
      DO j=103,136

      afunc4=afunc4+a(j)*mon(j)

      END DO
      afunc4=afunc4*temp(4)

      temp(5)=temp(5)+(1.0d0-tanh(g1(5)*(R1-R1_REF(5))))*
     *  (1.0d0-tanh(g2(5)*(R2-R2_REF(5))))*
     *  (1.0d0-tanh(g3(5)*(R3-R3_REF(5))))
c-----------------------------------------------------------------
      A11=SQRT(1.0D0/3.0D0)
      A12=SQRT(1.0D0/2.0D0)
      A13=SQRT(2.0D0/3.0D0)
      A14=SQRT(1.0D0/6.0D0)

      Q1(5)=(A11*(R1-R1_REF(5)))+(A11*(R2-R2_REF(5)))+
     * (A11*(R3-R3_REF(5)))
      Q2(5)=(A12*(R2-R2_REF(5)))-(A12*(R3-R3_REF(5)))      
      Q3(5)=(A13*(R1-R1_REF(5)))-(A14*(R2-R2_REF(5)))-
     * (A14*(R3-R3_REF(5)))

      S2A_SQ(5)=(Q2(5)*Q2(5))+(Q3(5)*Q3(5))
      S2B_SQ(5)=(Q2(5)*Q2(5))-(Q3(5)*Q3(5))
      S3_CU(5)=(Q3(5)*Q3(5)*Q3(5))-(3.0D0*Q2(5)*Q2(5)*Q3(5))

c---------------------------------------------------------------------
     
      
      mon(137)=1.0d0
      mon(138)=Q1(5)
      mon(139)=Q3(5)
      mon(140)=Q1(5)*Q1(5) 
      mon(141)=S2A_SQ(5)
      mon(142)=Q1(5)*Q3(5)  
      mon(143)=S2B_SQ(5)
      mon(144)=Q1(5)*Q1(5)*Q1(5)
      mon(145)=Q1(5)*S2A_SQ(5)
      mon(146)=S3_CU(5)
      mon(147)=Q1(5)*Q1(5)*Q3(5)
      mon(148)=Q1(5)*S2B_SQ(5)
      mon(149)=Q3(5)*S2A_SQ(5)
      mon(150)=Q1(5)*Q1(5)*Q1(5)*Q1(5)
      mon(151)=Q1(5)*Q1(5)*S2A_SQ(5)
      mon(152)=S2A_SQ(5)*S2A_SQ(5)
      mon(153)=Q1(5)*S3_CU(5)
      mon(154)=Q1(5)*Q1(5)*Q1(5)*Q3(5)
      mon(155)=Q1(5)*Q1(5)*S2B_SQ(5)
      mon(156)=Q1(5)*Q3(5)*S2A_SQ(5)
      mon(157)=Q3(5)*S3_CU(5)
      mon(158)=S2A_SQ(5)*S2B_SQ(5)
      mon(159)=Q1(5)**5
      mon(160)=Q1(5)*Q1(5)*Q1(5)*S2A_SQ(5)
      mon(161)=Q1(5)*S2A_SQ(5)*S2A_SQ(5)
      mon(162)=Q1(5)*Q1(5)*S3_CU(5)
      mon(163)=S2A_SQ(5)*S3_CU(5)
      mon(164)=Q1(5)*Q1(5)*Q1(5)*Q1(5)*Q3(5)
      mon(165)=Q1(5)*Q1(5)*Q1(5)*S2B_SQ(5)
      mon(166)=Q1(5)*Q1(5)*Q3(5)*S2A_SQ(5)
      mon(167)=Q1(5)*Q3(5)*S3_CU(5)
      mon(168)=Q1(5)*S2A_SQ(5)*S2B_SQ(5)
      mon(169)=Q3(5)*S2A_SQ(5)*S2A_SQ(5)
      mon(170)=S2B_SQ(5)*S3_CU(5)

      afunc5=0.0d0 
      DO j=137,170

      afunc5=afunc5+a(j)*mon(j)

      END DO
      afunc5=afunc5*temp(5)

    
      temp(6)=temp(6)+(1.0d0-tanh(g1(6)*(R1-R1_REF(6))))*
     *  (1.0d0-tanh(g2(6)*(R2-R2_REF(6))))*
     *  (1.0d0-tanh(g3(6)*(R3-R3_REF(6))))
c-----------------------------------------------------------------
      A11=SQRT(1.0D0/3.0D0)
      A12=SQRT(1.0D0/2.0D0)
      A13=SQRT(2.0D0/3.0D0)
      A14=SQRT(1.0D0/6.0D0)

      Q1(6)=(A11*(R1-R1_REF(6)))+(A11*(R2-R2_REF(6)))+
     * (A11*(R3-R3_REF(6)))
      Q2(6)=(A12*(R2-R2_REF(6)))-(A12*(R3-R3_REF(6)))      
      Q3(6)=(A13*(R1-R1_REF(6)))-(A14*(R2-R2_REF(6)))-
     * (A14*(R3-R3_REF(6)))

      S2A_SQ(6)=(Q2(6)*Q2(6))+(Q3(6)*Q3(6))
      S2B_SQ(6)=(Q2(6)*Q2(6))-(Q3(6)*Q3(6))
      S3_CU(6)=(Q3(6)*Q3(6)*Q3(6))-(3.0D0*Q2(6)*Q2(6)*Q3(6))

c---------------------------------------------------------------------
     
      mon(171)=1.0d0
      mon(172)=Q1(6)
      mon(173)=Q3(6)
      mon(174)=Q1(6)*Q1(6) 
      mon(175)=S2A_SQ(6)
      mon(176)=Q1(6)*Q3(6)  
      mon(177)=S2B_SQ(6)
      mon(178)=Q1(6)*Q1(6)*Q1(6)
      mon(179)=Q1(6)*S2A_SQ(6)
      mon(180)=S3_CU(6)
      mon(181)=Q1(6)*Q1(6)*Q3(6)
      mon(182)=Q1(6)*S2B_SQ(6)
      mon(183)=Q3(6)*S2A_SQ(6)
  
      afunc6=0.0d0 
      DO j=171,183

      afunc6=afunc6+a(j)*mon(j)

      END DO
      afunc6=afunc6*temp(6)


      afunc=afunc1+afunc2+afunc3+afunc4+afunc5+afunc6
      V3EHF=AFUNC
       
      return 
      end 
       
c-----------------------------------------------------------------------------
       DOUBLE PRECISION FUNCTION FR(R,CM)
       IMPLICIT NONE
       DOUBLE PRECISION ALPHA1,ALPH2,BETA1,BETA2
       DOUBLE PRECISION R110,R120,R111,R121,R1,R(3)
       DOUBLE PRECISION HR1,Gr1,RR1,RR10,delta,f,CM

       DOUBLE PRECISION ALPHA,R10,gr,hr,RR,RHH,RCHH

       alpha1          = 0.57135          
       R110            = 2.31574          
       BETA1           = 1.00176          
       R111            = 4.6892           
       ALPH2           = 0.70561          
       R120            = 3.8463           
       BETA2           = 0.138958         
       R121            = 5.71409  
        
       delta= 0.0466042608d0

       RHH=R(1)
       RCHH=CM

       ALPHA=0.75d0
       R10=5.5

       gr=0.5d0*(1.0d0+tanh(alpha*(RCHH-R10)))

       hr =0.0d0
       hr=delta*(0.25*(1.0d0-tanh(alpha1*(RHH-R110)+BETA1*(RHH-R111)**3)
     &   +(1.0d0-tanh(ALPH2*(RHH-R120)+BETA2*(RHH-R121)**3))))


       FR=GR*HR
c       write(33,*)RHH,RCHH,FR

       RETURN
       END

c------------------------------------------------------------------------------

       DOUBLE PRECISION FUNCTION V2EHFCH2(R)

c       DOUBLE PRECISION FUNCTION V2EHF(R)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION R(3)

c       OPEN (UNIT=40,FILE='two_body_ehf.dat')
C **  two body potential (EHFACE2U model used)

       V2EHFCH2=VEHF_HH(R(1),1)+VDC2CH2(R(1),1)+
     & VEHF_CH_CH2(R(2),2)+VDC2CH2(R(2),2)+
     & VEHF_CH_CH2(R(3),3)+VDC2CH2(R(3),3)

       RETURN
       END

c------------------------------------------------------------------------------

       DOUBLE PRECISION FUNCTION VEHF_CH_CH2(R,I)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A1(3),A2(3),A3(3),DD(3),RE(3),G0(3),G1(3),G2(3)
      DIMENSION A4(3),A5(3),A6(3),A7(3),A8(3),A9(3)

      DATA DD(2),A1(2),A2(2),A3(2),A4(2),A5(2),A6(2),A7(2),A8(2),
     &   A9(2)/0.2564169,1.7111764D0,1.05024831D0,1.43177832D0,
     &    0.366998978D0,0.5157509D0,0.113220497D0,0.0513221284D0,
     &   -0.0159231522D0,0.00502354172D0/
      
      DATA DD(3),A1(3),A2(3),A3(3),A4(3),A5(3),A6(3),A7(3),A8(3),
     &   A9(3)/0.2564169,1.7111764D0,1.05024831D0,1.43177832D0,
     &    0.366998978D0,0.5157509D0,0.113220497D0,0.0513221284D0,
     &   -0.0159231522D0,0.00502354172D0/

c A1=  1.7111764
c  1.05024831  1.43177832  0.366998978
c  0.5157509  0.113220497  0.0513221284 -0.0159231522  0.00502354172


c  1.17591785  1.75894627  0.502592586

C     DATA G0(1),G1(1),G2(1)/1.17591785D0,1.75894627d0,0.502592586d0/
      G1(2)=1.75894627d0/1.17591785D0
      DATA G0(2),G2(2)/1.17591785D0,0.502592586d0/

      G1(3)=1.75894627d0/1.17591785D0
      DATA G0(3),G2(3)/1.17591785D0,0.502592586d0/

   
      DATA RE(1),RE(2),RE(3)/2.35683488D0,2.1176D0,2.1176D0/
CC R(1) DONE

      RS=R-RE(I)

      GAMA=G0(I)*(1.0D0+G1(I)*TANH(G2(I)*RS))

      VEHF_CH_CH2 =-DD(I)*(1.0D0+A1(I)*RS+A2(I)*RS**2+A3(I)*RS**3+
     &          A4(I)*RS**4+A5(I)*RS**5+A6(I)*RS**6+A7(I)*RS**7+
     &          A8(I)*RS**8+A9(I)*RS**9)*EXP(-GAMA*RS)/R
    
       RETURN
       END
c---------------------------------------------------------------------------

      FUNCTION VEHF_HH(R,I)

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DD(1),A1(1),A2(1),A3(1),G0(1),G1(1),G2(1)
      DIMENSION RE(3)
      DATA DD(1),A1(1),A2(1),A3(1)/0.22979439D0,1.82027480D0,
     &     0.52437767D0,0.36999610D0/
      DATA G0(1),G1(1),G2(1)/1.094670D0,1.009737D0,0.235856D0/
      DATA RE(1),RE(2),RE(3)/1.401D0,2.120D0,2.120D0/

      RS=R-RE(I)

      GAMA=G0(I)*(1.0D0+G1(I)*TANH(G2(I)*RS))

      VEHF_HH =-DD(I)*(1.0D0+A1(I)*RS+A2(I)*RS*RS+A3(I)*RS*RS*RS)
     &      *EXP(-GAMA*RS)/R + AUXH(R,I)

      RETURN
      END
c---------------------------------------------------------------------------------

      FUNCTION AUXH(R,I)
      IMPLICIT REAL*8 (A-H,O-Z)

      DATA AP,ALPHAP,A1P,GAMAP/-0.8205D0,2.5D0,0.0d0,2.0D0/
c     CHIEXC=0.0d0

      CHIEXC=1.0d0*CHICH2(6,R,I)

      AUXH = AP*(R**ALPHAP)*(1.0D0+A1P*R)*EXP(-GAMAP*R)
      AUXH = AUXH*CHIEXC
c      write(*,*)AUXH,CHIEXC

      RETURN
      END


c----------------------------------------------

      FUNCTION VDC2CH2(R,I)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CA(9,9)
C ** dinamical correlation, diatom  calculation
c     C **  atom-atom dispersion coefficients for HH
      DATA CA(1,1),CA(1,2),CA(1,3),CA(1,4),CA(1,5),CA(1,6),CA(1,7),
     &  CA(1,8),CA(1,9)/6.499D0,124.4D0,3286.0D0,-3475.0D0,121500.0D0,
     &  -291400.0D0,6061000.0D0,-23050000.0D0,393800000.0D0/

C **  atom-atom dispersion coefficients for CH
      DATA CA(2,1),CA(2,2),CA(2,3)/16.1D0,351.109D0,10030.689D0/
C **  atom-atom dispersion coefficients for CH
      DATA CA(3,1),CA(3,2),CA(3,3)/16.1D0,351.109D0,10030.689D0/
      RIN1 =1.0D0/R

      RIN2 =RIN1*RIN1
      RIN6 =RIN2*RIN2*RIN2
      RIN8 =RIN6*RIN2
      RIN10=RIN8*RIN2
      RIN11=RIN10*RIN1
      RIN12=RIN10*RIN2
      RIN13=RIN12*RIN1
      RIN14=RIN12*RIN2
      RIN15=RIN14*RIN1
      RIN16=RIN14*RIN2
      
      IF(I.EQ.1)THEN
      VDC2CH2=-CA(I,1)*CHICH2(6,R,I)*RIN6-CA(I,2)*CHICH2(8,R,I)*RIN8
     &     -CA(I,3)*CHICH2(10,R,I)*RIN10-CA(I,4)*CHICH2(11,R,I)*RIN11
     &     -CA(I,5)*CHICH2(12,R,I)*RIN12-CA(I,6)*CHICH2(13,R,I)*RIN13
     &     -CA(I,7)*CHICH2(14,R,I)*RIN14-CA(I,8)*CHICH2(15,R,I)*RIN15
     &     -CA(I,9)*CHICH2(16,R,I)*RIN16

      ELSE 
      VDC2CH2=-CA(I,1)*CHICH2(6,R,I)*RIN6-CA(I,2)*CHICH2(8,R,I)*RIN8
     &     -CA(I,3)*CHICH2(10,R,I)*RIN10
      END IF

      RETURN
      END

c-------------------------------------------------------------------

      FUNCTION CHICH2(N,R,I)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A0(16),B0(16),R0(3)
C **  damping functions
      DATA (A0(J), J=1,16)/16.3661D0,10.062652D0,7.5709297D0,
     &      6.1869941D0,5.29027135D0,4.65478733D0,
     &      4.17772554D0,3.80388898D0,3.50229909D0,3.25255118D0,
     &      3.04213207D0,2.86194363D0,2.70562629,2.56852131D0,
     &      2.44713191D0,2.33877835D0/

      DATA (B0(J), J=1,16)/15.623647D0,14.19721174D0,12.9010098D0,
     &      11.723150819D0,10.65283005D0,9.68021804D0,
     &      8.79642676D0,7.99330588D0,7.26352748D0,6.60036154D0,
     &      5.99775016D0,5.45015706D0,4.95255907D0,4.50039165D0,
     &      4.089507D0,3.71613603D0/

      DATA R0(1),R0(2),R0(3)/6.9282D0,7.4096D0,7.4096D0/

      RHO=(11.0D0+2.5D0*R0(I))/2.0D0
      RR =R/RHO
      AUX=(1.0D0-EXP(-A0(N)*RR-B0(N)*RR*RR))**N
      CHICH2=AUX
C      DO IC=2,N
C      CHI=CHI*AUX
C      ENDDO

      RETURN		
      END

c-----------------------------------------------------------------------------
       DOUBLE PRECISION FUNCTION EHFCH2(R)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION R(3)
     
       EHFCH2=V3EHFCH2(R)
c       write(*,*)R(1),R(2),R(3),EHF
       RETURN
       END


C************************************************************************


c      subroutine polycal(ap,ma,ndat,g1,g2,g3,R1_REF,R2_REF,R3_REF,x,bb,c
c      * ,y,afunc,np)

 
      DOUBLE PRECISION FUNCTION V3EHFCH2(R)
      implicit none
      INTEGER ma,j,i,np,k,kk,jj,ndat,NCOF
      DOUBLE PRECISION x,bb,c,g1(7),g2(7),g3(7)

      DOUBLE PRECISION R1_REF(7),R2_REF(7),R3_REF(7),mon(140)

      DOUBLE PRECISION R1,R2,R3,temp(7),yy,kk1,R(3),CM,a(140)

      DOUBLE PRECISION A11,A12,A13,A14,Q1(7),Q2(7),Q3(7)
      DOUBLE PRECISION S2A_SQ(7),S2B_SQ(7),S3_CU(7)

      DOUBLE PRECISION  afunc,afunc1,afunc2,afunc3,afunc4,afunc5,afunc6
      DOUBLE PRECISION  afunc7
c      DOUBLE PRECISION temp1,g11,R11_REF,g22,R22_REF,g33,R33_REF

c      DOUBLE PRECISION Q11,Q22,Q33,S2A_SQt,S2B_SQt,S3_CUt
c-------------------------------------------------------------------


      ma=140
      np=7

      R1=R(1)
      R2=R(2)
      R3=R(3)

      CM=SQRT(2.0D0*R(3)**2+2.0D0*R(2)**2-R(1)**2)*1.0D0/2.0D0



c**************************************************************************

      R1_REF(1)=1.4d0
      R2_REF(1)=4.0d0
      R3_REF(1)=4.0d0

      g1(1)=1.5d0
      g2(1)=0.50d0
      g3(1)=0.50d0 

c-----------------------------------------------------------------

      R1_REF(2)=2.0d0
      R2_REF(2)=2.5d0
      R3_REF(2)=2.5d0

      g1(2)=0.25d0
      g2(2)=0.75d0
      g3(2)=0.75d0
 
c-----------------------------------------------------------------

      R1_REF(3)=2.6d0
      R2_REF(3)=3.0d0
      R3_REF(3)=3.0d0

      g1(3)=0.25d0
      g2(3)=0.85d0
      g3(3)=0.85d0

c--------------------------------------------------------------------

      R1_REF(4)=4.18d0
      R2_REF(4)=2.09d0
      R3_REF(4)=2.09d0

      g1(4)=0.5d0
      g2(4)=0.8d0
      g3(4)=0.8d0

c-----------------------------------------------------------------

      R1_REF(5)=3.5d0
      R2_REF(5)=2.0d0
      R3_REF(5)=2.0d0

      g1(5)=1.00d0
      g2(5)=0.65d0
      g3(5)=0.65d0

c     g1(5)=1.25d0 gives better result
c-----------------------------------------------------------------

      R1_REF(6)=1.8d0
      R2_REF(6)=4.0d0
      R3_REF(6)=4.0d0

      g1(6)=1.25d0
      g2(6)=0.50d0
      g3(6)=0.50d0

c-----------------------------------------------------------------

      R1_REF(7)=1.40d0
      R2_REF(7)=7.5d0
      R3_REF(7)=7.5d0

      g1(7)=1.55d0
      g2(7)=0.50d0
      g3(7)=0.50d0

c-----------------------------------------------------------------


      do kk=1,ma
      a(kk)=0.0d0
      mon(kk)=0.0d0
      end do


      do k=1,np
      temp(k)=0.0d0

      Q1(k)=0.0d0
      Q2(k)=0.0d0
      Q3(k)=0.0d0
      S2A_SQ(k)=0.0d0
      S2B_SQ(k)=0.0d0
      S3_CU(k)=0.0d0

      end do
       a(1)=    -0.2548616756
       a(2)=     1.5147884986
       a(3)=    -0.7647569043
       a(4)=     0.0629362915
       a(5)=    -0.0169770387
       a(6)=     0.0217555338
       a(7)=    -0.0188944626
       a(8)=     0.4580169007
       a(9)=     0.0951527773
       a(10)=    -0.1487483758
       a(11)=     0.9760440271
       a(12)=    -0.0447137723
       a(13)=    -0.3426147134
       a(14)=    10.8354811567
       a(15)=   -10.5865202619
       a(16)=     7.4824695662
       a(17)=     1.0228346686
       a(18)=     0.5622127298
       a(19)=     2.0367011360
       a(20)=     1.9322027027
       a(21)=    -0.2761664053
       a(22)=    -0.1094912406
       a(23)=    -0.1676689851
       a(24)=     0.7277462111
       a(25)=     0.9997555730
       a(26)=     0.1247901442
       a(27)=     0.0882297187
       a(28)=     0.0404032676
       a(29)=     0.0743084519
       a(30)=    -0.0489327878
       a(31)=    -0.1991435898
       a(32)=    -0.2551258851
       a(33)=    -0.1899806346
       a(34)=     0.0526724060
       a(35)=    -0.0812363486
       a(36)=    -2.7327136557
       a(37)=     0.0081689512
       a(38)=    -0.4891841316
       a(39)=    -0.7737525740
       a(40)=    -0.6153943709
       a(41)=     1.2136158104
       a(42)=    -0.2205732853
       a(43)=     0.0157693506
       a(44)=    -0.0450150645
       a(45)=     0.0227956860
       a(46)=     0.0870311840
       a(47)=     0.0021740287
       a(48)=    -0.0467568504
       a(49)=    -0.0513873112
       a(50)=    -0.0618249586
       a(51)=    -0.0127101500
       a(52)=     0.0438626594
       a(53)=     0.1344959804
       a(54)=     0.1265143486
       a(55)=     0.0403433756
       a(56)=    -0.0145197248
       a(57)=    -0.0369158032
       a(58)=   -19.7909022415
       a(59)=    -3.6225414115
       a(60)=    -2.8759754175
       a(61)=    -1.6924170761
       a(62)=    -4.0495070802
       a(63)=     3.4673921891
       a(64)=     0.4873716655
       a(65)=    -0.1862519656
       a(66)=    -0.6670666051
       a(67)=     0.0840928720
       a(68)=     0.2376873666
       a(69)=    -0.5095111076
       a(70)=    -0.1721073218
       a(71)=    19.1311858893
       a(72)=    15.7918273878
       a(73)=    12.2611168574
       a(74)=     7.6687615685
       a(75)=     6.9697824842
       a(76)=    11.2340678037
       a(77)=    -7.1021863474
       a(78)=     3.0767202484
       a(79)=     5.5077236609
       a(80)=     1.2396197211
       a(81)=     8.0099396151
       a(82)=    -5.7672599863
       a(83)=     4.4680493749
       a(84)=     0.7163962805
       a(85)=     1.8200446430
       a(86)=     0.9601146581
       a(87)=     1.2599572054
       a(88)=     2.6003690370
       a(89)=    -1.9843479125
       a(90)=     3.8655789216
       a(91)=     0.6027127596
       a(92)=    -0.8779701920
       a(93)=     0.1127192932
       a(94)=     0.3137806563
       a(95)=     0.5537130011
       a(96)=     0.3031883058
       a(97)=     0.1439066888
       a(98)=     0.4167627410
       a(99)=    -0.3150089962
       a(100)=     0.8147976847
       a(101)=     0.3350065719
       a(102)=    -0.5046906058
       a(103)=     0.2726524075
       a(104)=    -0.0372204848
       a(105)=     0.0074638143
       a(106)=     0.0254808161
       a(107)=     0.0723021983
       a(108)=     0.0072632298
       a(109)=     0.0318593496
       a(110)=     0.0076525518
       a(111)=     0.0010359332
       a(112)=     0.0291959021
       a(113)=    -0.0389215242
       a(114)=     0.0319926866
       a(115)=     0.0369035911
       a(116)=    -0.0594273971
       a(117)=     0.0968703646
       a(118)=    -0.0199404621
       a(119)=     0.0003288123
       a(120)=    -0.0131151715
       a(121)=    -0.5704030626
       a(122)=    -1.2580255815
       a(123)=     0.0613735136
       a(124)=    -0.1033935767
       a(125)=    -0.1294763492
       a(126)=    -0.2718137827
       a(127)=     0.0956992581
       a(128)=    -0.2524117518
       a(129)=    -0.1813929507
       a(130)=     0.0624384194
       a(131)=    -0.6686709847
       a(132)=     0.1468953485
       a(133)=     0.1113853186
       a(134)=    -0.0029512607
       a(135)=    -0.0029562823
       a(136)=    -0.0050546780
       a(137)=    -0.0038821988
       a(138)=    -0.0039487057
       a(139)=    -0.0105361164
       a(140)=     0.0038290866


      temp(1)=temp(1)+(1.0d0-tanh(g1(1)*(R1-R1_REF(1))))*
     *  (1.0d0-tanh(g2(1)*(R2-R2_REF(1))))*
     *  (1.0d0-tanh(g3(1)*(R3-R3_REF(1))))
c-----------------------------------------------------------------
      A11=SQRT(1.0D0/3.0D0)
      A12=SQRT(1.0D0/2.0D0)
      A13=SQRT(2.0D0/3.0D0)
      A14=SQRT(1.0D0/6.0D0)

      Q1(1)=(A11*(R1-R1_REF(1)))+(A11*(R2-R2_REF(1)))+
     * (A11*(R3-R3_REF(1)))
      Q2(1)=(A12*(R2-R2_REF(1)))-(A12*(R3-R3_REF(1)))      
      Q3(1)=(A13*(R1-R1_REF(1)))-(A14*(R2-R2_REF(1)))-
     * (A14*(R3-R3_REF(1)))

      S2A_SQ(1)=(Q2(1)*Q2(1))+(Q3(1)*Q3(1))
      S2B_SQ(1)=(Q2(1)*Q2(1))-(Q3(1)*Q3(1))
      S3_CU(1)=(Q3(1)*Q3(1)*Q3(1))-(3.0D0*Q2(1)*Q2(1)*Q3(1))



c      write(10,*),temp,Q1,Q2,Q3,S2A_SQ,S2B_SQ,S3_CU
c-----------------------------------------------------------------

ctemp1,g11,R11_REF,g22,R22_REF,g33,R33_REF,Q11,Q22,Q33,S2A_SQt,S2B_SQt,S3_CUt

c---------------------------------------------------------------------
     
      mon(1)=1.0d0
      mon(2)=Q1(1)
      mon(3)=Q3(1)
      mon(4)=Q1(1)*Q1(1) 
      mon(5)=S2A_SQ(1)
      mon(6)=Q1(1)*Q3(1)  
      mon(7)=S2B_SQ(1)
      mon(8)=Q1(1)*Q1(1)*Q1(1)
      mon(9)=Q1(1)*S2A_SQ(1)
      mon(10)=S3_CU(1)
      mon(11)=Q1(1)*Q1(1)*Q3(1)
      mon(12)=Q1(1)*S2B_SQ(1)
      mon(13)=Q3(1)*S2A_SQ(1)

c-----------------------------------------------------
      afunc1=0.0d0
      DO j=1,13

      afunc1=afunc1+a(j)*mon(j)

      END DO
      afunc1=afunc1*temp(1)

c*************************************************************************

     
      temp(2)=temp(2)+(1.0d0-tanh(g1(2)*(R1-R1_REF(2))))*
     *  (1.0d0-tanh(g2(2)*(R2-R2_REF(2))))*
     *  (1.0d0-tanh(g3(2)*(R3-R3_REF(2))))
c-----------------------------------------------------------------
      A11=SQRT(1.0D0/3.0D0)
      A12=SQRT(1.0D0/2.0D0)
      A13=SQRT(2.0D0/3.0D0)
      A14=SQRT(1.0D0/6.0D0)

      Q1(2)=(A11*(R1-R1_REF(2)))+(A11*(R2-R2_REF(2)))+
     * (A11*(R3-R3_REF(2)))
      Q2(2)=(A12*(R2-R2_REF(2)))-(A12*(R3-R3_REF(2)))      
      Q3(2)=(A13*(R1-R1_REF(2)))-(A14*(R2-R2_REF(2)))-
     * (A14*(R3-R3_REF(2)))

      S2A_SQ(2)=(Q2(2)*Q2(2))+(Q3(2)*Q3(2))
      S2B_SQ(2)=(Q2(2)*Q2(2))-(Q3(2)*Q3(2))
      S3_CU(2)=(Q3(2)*Q3(2)*Q3(2))-(3.0D0*Q2(2)*Q2(2)*Q3(2))

c---------------------------------------------------------------------

      mon(14)=1.0d0
      mon(15)=Q1(2)
      mon(16)=Q3(2)
      mon(17)=Q1(2)*Q1(2) 
      mon(18)=S2A_SQ(2)
      mon(19)=Q1(2)*Q3(2)  
      mon(20)=S2B_SQ(2)
      mon(21)=Q1(2)*Q1(2)*Q1(2)
      mon(22)=Q1(2)*S2A_SQ(2)
      mon(23)=S3_CU(2)
      mon(24)=Q1(2)*Q1(2)*Q3(2)
      mon(25)=Q1(2)*S2B_SQ(2)
      mon(26)=Q3(2)*S2A_SQ(2)
      mon(27)=Q1(2)*Q1(2)*Q1(2)*Q1(2)
      mon(28)=Q1(2)*Q1(2)*S2A_SQ(2)
      mon(29)=S2A_SQ(2)*S2A_SQ(2)
      mon(30)=Q1(2)*S3_CU(2)
      mon(31)=Q1(2)*Q1(2)*Q1(2)*Q3(2)
      mon(32)=Q1(2)*Q1(2)*S2B_SQ(2)
      mon(33)=Q1(2)*Q3(2)*S2A_SQ(2)
      mon(34)=Q3(2)*S3_CU(2)
      mon(35)=S2A_SQ(2)*S2B_SQ(2)
  
      afunc2=0.0d0 
      DO j=14,35

      afunc2=afunc2+a(j)*mon(j)

      END DO
      afunc2=afunc2*temp(2)




    
      temp(3)=temp(3)+(1.0d0-tanh(g1(3)*(R1-R1_REF(3))))*
     *  (1.0d0-tanh(g2(3)*(R2-R2_REF(3))))*
     *  (1.0d0-tanh(g3(3)*(R3-R3_REF(3))))
c-----------------------------------------------------------------
      A11=SQRT(1.0D0/3.0D0)
      A12=SQRT(1.0D0/2.0D0)
      A13=SQRT(2.0D0/3.0D0)
      A14=SQRT(1.0D0/6.0D0)

      Q1(3)=(A11*(R1-R1_REF(3)))+(A11*(R2-R2_REF(3)))+
     * (A11*(R3-R3_REF(3)))
      Q2(3)=(A12*(R2-R2_REF(3)))-(A12*(R3-R3_REF(3)))      
      Q3(3)=(A13*(R1-R1_REF(3)))-(A14*(R2-R2_REF(3)))-
     * (A14*(R3-R3_REF(3)))

      S2A_SQ(3)=(Q2(3)*Q2(3))+(Q3(3)*Q3(3))
      S2B_SQ(3)=(Q2(3)*Q2(3))-(Q3(3)*Q3(3))
      S3_CU(3)=(Q3(3)*Q3(3)*Q3(3))-(3.0D0*Q2(3)*Q2(3)*Q3(3))

c---------------------------------------------------------------------
     
      mon(36)=1.0d0
      mon(37)=Q1(3)
      mon(38)=Q3(3)
      mon(39)=Q1(3)*Q1(3) 
      mon(40)=S2A_SQ(3)
      mon(41)=Q1(3)*Q3(3)  
      mon(42)=S2B_SQ(3)
      mon(43)=Q1(3)*Q1(3)*Q1(3)
      mon(44)=Q1(3)*S2A_SQ(3)
      mon(45)=S3_CU(3)
      mon(46)=Q1(3)*Q1(3)*Q3(3)
      mon(47)=Q1(3)*S2B_SQ(3)
      mon(48)=Q3(3)*S2A_SQ(3)
      mon(49)=Q1(3)*Q1(3)*Q1(3)*Q1(3)
      mon(50)=Q1(3)*Q1(3)*S2A_SQ(3)
      mon(51)=S2A_SQ(3)*S2A_SQ(3)
      mon(52)=Q1(3)*S3_CU(3)
      mon(53)=Q1(3)*Q1(3)*Q1(3)*Q3(3)
      mon(54)=Q1(3)*Q1(3)*S2B_SQ(3)
      mon(55)=Q1(3)*Q3(3)*S2A_SQ(3)
      mon(56)=Q3(3)*S3_CU(3)
      mon(57)=S2A_SQ(3)*S2B_SQ(3)
 
      afunc3=0.0d0 
      DO j=36,57

      afunc3=afunc3+a(j)*mon(j)

      END DO
      afunc3=afunc3*temp(3)


      temp(4)=temp(4)+(1.0d0-tanh(g1(4)*(R1-R1_REF(4))))*
     *  (1.0d0-tanh(g2(4)*(R2-R2_REF(4))))*
     *  (1.0d0-tanh(g3(4)*(R3-R3_REF(4))))
c-----------------------------------------------------------------
      A11=SQRT(1.0D0/3.0D0)
      A12=SQRT(1.0D0/2.0D0)
      A13=SQRT(2.0D0/3.0D0)
      A14=SQRT(1.0D0/6.0D0)

      Q1(4)=(A11*(R1-R1_REF(4)))+(A11*(R2-R2_REF(4)))+
     * (A11*(R3-R3_REF(4)))
      Q2(4)=(A12*(R2-R2_REF(4)))-(A12*(R3-R3_REF(4)))      
      Q3(4)=(A13*(R1-R1_REF(4)))-(A14*(R2-R2_REF(4)))-
     * (A14*(R3-R3_REF(4)))

      S2A_SQ(4)=(Q2(4)*Q2(4))+(Q3(4)*Q3(4))
      S2B_SQ(4)=(Q2(4)*Q2(4))-(Q3(4)*Q3(4))
      S3_CU(4)=(Q3(4)*Q3(4)*Q3(4))-(3.0D0*Q2(4)*Q2(4)*Q3(4))

c---------------------------------------------------------------------
     
     
      mon(58)=1.0d0
      mon(59)=Q1(4)
      mon(60)=Q3(4)
      mon(61)=Q1(4)*Q1(4) 
      mon(62)=S2A_SQ(4)
      mon(63)=Q1(4)*Q3(4)  
      mon(64)=S2B_SQ(4)
      mon(65)=Q1(4)*Q1(4)*Q1(4)
      mon(66)=Q1(4)*S2A_SQ(4)
      mon(67)=S3_CU(4)
      mon(68)=Q1(4)*Q1(4)*Q3(4)
      mon(69)=Q1(4)*S2B_SQ(4)
      mon(70)=Q3(4)*S2A_SQ(4)
 
  
      afunc4=0.0d0 
      DO j=58,70

      afunc4=afunc4+a(j)*mon(j)

      END DO
      afunc4=afunc4*temp(4)

      temp(5)=temp(5)+(1.0d0-tanh(g1(5)*(R1-R1_REF(5))))*
     *  (1.0d0-tanh(g2(5)*(R2-R2_REF(5))))*
     *  (1.0d0-tanh(g3(5)*(R3-R3_REF(5))))
c-----------------------------------------------------------------
      A11=SQRT(1.0D0/3.0D0)
      A12=SQRT(1.0D0/2.0D0)
      A13=SQRT(2.0D0/3.0D0)
      A14=SQRT(1.0D0/6.0D0)

      Q1(5)=(A11*(R1-R1_REF(5)))+(A11*(R2-R2_REF(5)))+
     * (A11*(R3-R3_REF(5)))
      Q2(5)=(A12*(R2-R2_REF(5)))-(A12*(R3-R3_REF(5)))      
      Q3(5)=(A13*(R1-R1_REF(5)))-(A14*(R2-R2_REF(5)))-
     * (A14*(R3-R3_REF(5)))

      S2A_SQ(5)=(Q2(5)*Q2(5))+(Q3(5)*Q3(5))
      S2B_SQ(5)=(Q2(5)*Q2(5))-(Q3(5)*Q3(5))
      S3_CU(5)=(Q3(5)*Q3(5)*Q3(5))-(3.0D0*Q2(5)*Q2(5)*Q3(5))

c---------------------------------------------------------------------
     
      mon(71)=1.0d0
      mon(72)=Q1(5)
      mon(73)=Q3(5)
      mon(74)=Q1(5)*Q1(5) 
      mon(75)=S2A_SQ(5)
      mon(76)=Q1(5)*Q3(5)  
      mon(77)=S2B_SQ(5)
      mon(78)=Q1(5)*Q1(5)*Q1(5)
      mon(79)=Q1(5)*S2A_SQ(5)
      mon(80)=S3_CU(5)
      mon(81)=Q1(5)*Q1(5)*Q3(5)
      mon(82)=Q1(5)*S2B_SQ(5)
      mon(83)=Q3(5)*S2A_SQ(5)
      mon(84)=Q1(5)*Q1(5)*Q1(5)*Q1(5)
      mon(85)=Q1(5)*Q1(5)*S2A_SQ(5)
      mon(86)=S2A_SQ(5)*S2A_SQ(5)
      mon(87)=Q1(5)*S3_CU(5)
      mon(88)=Q1(5)*Q1(5)*Q1(5)*Q3(5)
      mon(89)=Q1(5)*Q1(5)*S2B_SQ(5)
      mon(90)=Q1(5)*Q3(5)*S2A_SQ(5)
      mon(91)=Q3(5)*S3_CU(5)
      mon(92)=S2A_SQ(5)*S2B_SQ(5)
      mon(93)=Q1(5)**5
      mon(94)=Q1(5)*Q1(5)*Q1(5)*S2A_SQ(5)
      mon(95)=Q1(5)*S2A_SQ(5)*S2A_SQ(5)
      mon(96)=Q1(5)*Q1(5)*S3_CU(5)
      mon(97)=S2A_SQ(5)*S3_CU(5)
      mon(98)=Q1(5)*Q1(5)*Q1(5)*Q1(5)*Q3(5)
      mon(99)=Q1(5)*Q1(5)*Q1(5)*S2B_SQ(5)
      mon(100)=Q1(5)*Q1(5)*Q3(5)*S2A_SQ(5)
      mon(101)=Q1(5)*Q3(5)*S3_CU(5)
      mon(102)=Q1(5)*S2A_SQ(5)*S2B_SQ(5)
      mon(103)=Q3(5)*S2A_SQ(5)*S2A_SQ(5)
      mon(104)=S2B_SQ(5)*S3_CU(5)
      mon(105)=Q1(5)**6
      mon(106)=(Q1(5)**4)*S2A_SQ(5)
      mon(107)=Q1(5)*Q1(5)*S2A_SQ(5)*S2A_SQ(5)
      mon(108)=Q1(5)*Q1(5)*Q1(5)*S3_CU(5)
      mon(109)=Q1(5)*S2A_SQ(5)*S3_CU(5)
      mon(110)=S2A_SQ(5)*S2A_SQ(5)*S2A_SQ(5)
      mon(111)=S3_CU(5)*S3_CU(5)
      mon(112)=(Q1(5)**5)*Q3(5)
      mon(113)=(Q1(5)**4)*S2B_SQ(5)
      mon(114)=Q1(5)*Q1(5)*Q1(5)*Q3(5)*S2A_SQ(5)
      mon(115)=Q1(5)*Q1(5)*Q3(5)*S3_CU(5)
      mon(116)=Q1(5)*Q1(5)*S2A_SQ(5)*S2B_SQ(5)
      mon(117)=Q1(5)*Q3(5)*S2A_SQ(5)*S2A_SQ(5)
      mon(118)=Q1(5)*S2B_SQ(5)*S3_CU(5)
      mon(119)=Q3(5)*S2A_SQ(5)*S3_CU(5)
      mon(120)=S2A_SQ(5)*S2A_SQ(5)*S2B_SQ(5)

  
      afunc5=0.0d0 
      DO j=71,120

      afunc5=afunc5+a(j)*mon(j)

      END DO
      afunc5=afunc5*temp(5)




    
      temp(6)=temp(6)+(1.0d0-tanh(g1(6)*(R1-R1_REF(6))))*
     *  (1.0d0-tanh(g2(6)*(R2-R2_REF(6))))*
     *  (1.0d0-tanh(g3(6)*(R3-R3_REF(6))))
c-----------------------------------------------------------------
      A11=SQRT(1.0D0/3.0D0)
      A12=SQRT(1.0D0/2.0D0)
      A13=SQRT(2.0D0/3.0D0)
      A14=SQRT(1.0D0/6.0D0)

      Q1(6)=(A11*(R1-R1_REF(6)))+(A11*(R2-R2_REF(6)))+
     * (A11*(R3-R3_REF(6)))
      Q2(6)=(A12*(R2-R2_REF(6)))-(A12*(R3-R3_REF(6)))      
      Q3(6)=(A13*(R1-R1_REF(6)))-(A14*(R2-R2_REF(6)))-
     * (A14*(R3-R3_REF(6)))

      S2A_SQ(6)=(Q2(6)*Q2(6))+(Q3(6)*Q3(6))
      S2B_SQ(6)=(Q2(6)*Q2(6))-(Q3(6)*Q3(6))
      S3_CU(6)=(Q3(6)*Q3(6)*Q3(6))-(3.0D0*Q2(6)*Q2(6)*Q3(6))

c---------------------------------------------------------------------
     
      mon(121)=1.0d0
      mon(122)=Q1(6)
      mon(123)=Q3(6)
      mon(124)=Q1(6)*Q1(6) 
      mon(125)=S2A_SQ(6)
      mon(126)=Q1(6)*Q3(6)  
      mon(127)=S2B_SQ(6)
      mon(128)=Q1(6)*Q1(6)*Q1(6)
      mon(129)=Q1(6)*S2A_SQ(6)
      mon(130)=S3_CU(6)
      mon(131)=Q1(6)*Q1(6)*Q3(6)
      mon(132)=Q1(6)*S2B_SQ(6)
      mon(133)=Q3(6)*S2A_SQ(6)
 
 
      afunc6=0.0d0 
      DO j=121,133

      afunc6=afunc6+a(j)*mon(j)

      END DO
      afunc6=afunc6*temp(6)


    
      temp(7)=temp(7)+(1.0d0-tanh(g1(7)*(R1-R1_REF(7))))*
     *  (1.0d0-tanh(g2(7)*(R2-R2_REF(7))))*
     *  (1.0d0-tanh(g3(7)*(R3-R3_REF(7))))
c-----------------------------------------------------------------
      A11=SQRT(1.0D0/3.0D0)
      A12=SQRT(1.0D0/2.0D0)
      A13=SQRT(2.0D0/3.0D0)
      A14=SQRT(1.0D0/6.0D0)

      Q1(7)=(A11*(R1-R1_REF(7)))+(A11*(R2-R2_REF(7)))+
     * (A11*(R3-R3_REF(7)))
      Q2(7)=(A12*(R2-R2_REF(7)))-(A12*(R3-R3_REF(7)))      
      Q3(7)=(A13*(R1-R1_REF(7)))-(A14*(R2-R2_REF(7)))-
     * (A14*(R3-R3_REF(7)))

      S2A_SQ(7)=(Q2(7)*Q2(7))+(Q3(7)*Q3(7))
      S2B_SQ(7)=(Q2(7)*Q2(7))-(Q3(7)*Q3(7))
      S3_CU(7)=(Q3(7)*Q3(7)*Q3(7))-(3.0D0*Q2(7)*Q2(7)*Q3(7))

c---------------------------------------------------------------------
     
      mon(134)=1.0d0
      mon(135)=Q1(7)
      mon(136)=Q3(7)
      mon(137)=Q1(7)*Q1(7) 
      mon(138)=S2A_SQ(7)
      mon(139)=Q1(7)*Q3(7)  
      mon(140)=S2B_SQ(7)
c      mon(123)=Q1(7)*Q1(7)*Q1(7)
c      mon(124)=Q1(7)*S2A_SQ(7)
c      mon(125)=S3_CU(7)
c      mon(126)=Q1(7)*Q1(7)*Q3(7)
c      mon(127)=Q1(7)*S2B_SQ(7)
c      mon(128)=Q3(7)*S2A_SQ(7)

      afunc7=0.0d0 

      DO j=134,140

      afunc7=afunc7+a(j)*mon(j)

      END DO
   
      afunc7=afunc7*temp(7)


      afunc=afunc1+afunc2+afunc3+afunc4+afunc5+afunc6+afunc7
      V3EHFCH2=AFUNC

c      write(*,*)'afunc',afunc,afunc1,afunc2,afunc3,afunc4
      return 
      end 


c------------------------------------------------------------
      DOUBLE PRECISION FUNCTION V3DCCH2(R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(3),RS(3),CT(3)
      s1=R(1)
      s2=R(2)
      s3=R(3)

      f1=0.5D0*(1.D0-DTANH(1.D0*(6.D0*s1-s2-s3)))
      f2=0.5D0*(1.D0-DTANH(1.D0*(6.D0*s2-s3-s1)))
      f3=0.5D0*(1.D0-DTANH(1.D0*(6.D0*s3-s1-s2)))

      CALL INT_JAC_CH2(R,RS,CT)

 
      V3DCCH2=f1 * (CCH2(6,R(1),CT(1)) * CHIDCH2(6,RS(1),3) / RS(1)**6
     & +CCH2(8,R(1),CT(1)) * CHIDCH2(8,RS(1),3) / RS(1)**8
     & +CCH2(10,R(1),CT(1)) * CHIDCH2(10,RS(1),3) / RS(1)**10)
     & +f2 * (CHCH(6,R(2),CT(2)) * CHIDCH2(6,RS(2),3) / RS(2)**6
     & +CHCH(8,R(2),CT(2)) * CHIDCH2(8,RS(2),3) / RS(2)**8
     & +CHCH(10,R(2),CT(2)) * CHIDCH2(10,RS(2),2) / RS(2)**10)
     & +f3 * (CHCH(6,R(3),CT(3)) * CHIDCH2(6,RS(3),3) / RS(3)**6
     & +CHCH(8,R(3),CT(3)) * CHIDCH2(8,RS(3),3) / RS(3)**8
     & +CHCH(10,R(3),CT(3)) * CHIDCH2(10,RS(3),3) / RS(3)**10)

c      v3dc=1.0d0*v3dc
      V3DCCH2=-1.d0*V3DCCH2

      RETURN
      END      
c--------------------------------------------------------------------
       DOUBLE PRECISION FUNCTION CCH2(N,R,CT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION D(6),RM(6),C1(6),C2(6),C3(6),C4(6),C5(6),CINF(6)
c---------------------------------------------------------------
      CND(X,I)=D(I)*(1.0D0+C1(I)*(X-RM(I))+C2(I)*(X-RM(I))*(X-RM(I))
     &         +C3(I)*(X-RM(I))*(X-RM(I))*(X-RM(I)))*EXP(-C1(I)
     &         *(X-RM(I))-C4(I)*(X-RM(I))*(X-RM(I))
     &         -C5(I)*(X-RM(I))*(X-RM(I))*(X-RM(I)))+CINF(I)
c--------------------------------------------------------------------
c COEFFICIENTS ARE NUMBERED IN THE FOLLOWING SEQUENCE 
C    C60,C62,C80,C82,C10,C84(I=1,2,3,4,5,6 RESPECTIVELY)
C RM AND CINF ALSO FOLLOW THE SAME ORDER----------------------------
c------------------------------------------------------------------------
       DATA (RM(I),I=1,6)/3.4158D0,2.1763D0,3.4084D0,2.2491D0,3.4034,
     &       2.1349d0/ 
       DATA (D(I),I=1,6)/8.3452D0,5.4108D0,237.5104D0,372.7948D0,
     &     8471.3D0,36.02d0/
       DATA (CINF(I),I=1,6)/32.2D0,0.0D0,702.216D0,0.0D0,20061.0D0,0.00/
       DATA (C1(I),I=1,6)/1.21159742D0,0.539024102D0,1.26370383D0,
     &     1.20519469D0,0.112850404D0,1.31253801d0/
       DATA (C2(I),I=1,6)/0.38234715D0,-0.0093106543D0,0.43348456D0,
     &     0.38802817D0,-0.361286218D0,0.4293d0/
      DATA (C3(I),I=1,6)/0.0477008753D0,0.01262754D0,0.052455744D0,
     &     0.0201211879D0,0.169153437D0,0.00916632d0/
      DATA (C4(I),I=1,6)/0.21049927D0,0.190215654D0,0.1988556338D0,
     &     -0.0157850995D0,0.282169598D0,0.32045507d0/
      DATA (C5(I),I=1,6)/0.000000001D0,0.00000001D0,0.000000001D0,
     &     0.000000001D0,0.032247798,0.00000001/
c-----------------------------------------------------------------
      IF(N.EQ.6)THEN

      CCH2=0.0d0 
      CCH2=CND(R,1)+CND(R,2)*POLEGCH2(CT,2)
c     write(60,*)R,CND(R,1)
c     write(62,*)R,CND(R,2)*POLEG(CT,2)
      END IF

      IF(N.EQ.8)THEN
      CCH2=0.0d0
      CCH2=CND(R,3)+CND(R,4)*POLEGCH2(CT,2)+CND(R,6)*POLEGCH2(CT,3)
      
c     write(60,*)R,CCH2
      END IF

      IF(N.EQ.10)THEN
      CCH2=0.0d0
      CCH2=CND(R,5)
      END IF

      RETURN
      END
       
c--------------------------------------------------------------------------

     
      DOUBLE PRECISION FUNCTION CHCH(N,R,CT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D(6),RM(6),C1(6),C2(6),C3(6),C4(6),C5(6),CINF(6)

      CND(X,I)=D(I)*(1.0D0+C1(I)*(X-RM(I))+C2(I)*(X-RM(I))*(X-RM(I))
     &         +C3(I)*(X-RM(I))*(X-RM(I))*(X-RM(I)))*EXP(-C1(I)
     &         *(X-RM(I))-C4(I)*(X-RM(I))*(X-RM(I))
     &         -C5(I)*(X-RM(I))*(X-RM(I))*(X-RM(I)))+CINF(I)

c--------------------------------------------------------------------
c COEFFICIENTS ARE NUMBERED IN THE FOLLOWING SEQUENCE
C    C60,C62,C80,C82,C10,C84(I=1,2,3,4,5,6 RESPECTIVELY)
C RM AND CINF ALSO FOLLOW THE SAME ORDER----------------------------
c------------------------------------------------------------------------

      DATA (D(I),I=1,6)/4.8004D0,3.9195D0,112.8635D0,323.1503D0,
     &      4092.2957D0,28.9647D0/

      DATA (RM(I),I=1,6)/3.9560D0,3.7371D0,3.9404D0,3.7537D0,
     &      3.9305D0,3.7142D0/

      DATA (CINF(I),I=1,6)/22.60D0,0.0D0,479.1828D0,0.0D0,13336.0D0,0.0/

      DATA (C1(I),I=1,6)/1.40840499D0,0.49960201D0,1.38918659D0,
     &     0.984737909D0,1.37904956D0,0.883426768/

      DATA (C2(I),I=1,6)/0.57185799D0,0.0098876426D0,0.56988421D0,
     &     0.29958805D0,0.57002801D0,0.152664106/

      DATA (C3(I),I=1,6)/0.0811532658D0,0.000629703D0,0.081266045D0,
     &     0.028894432D0,0.081674607D0,-0.0206191526/

      DATA (C4(I),I=1,6)/0.320208885D0,0.203051132D0,0.3082019D0,
     &     0.11244944D0,0.304710584D0,0.397514965/


      DATA (C5(I),I=1,6)/0.0145981138D0,1.27207289D-08,0.0103661606D0,
     &     4.742902D-09,0.00789332256D0,4.81907112D-09/

c----------------------------------------------------------------------

      IF(N.EQ.6)THEN
      CHCH=0.0d0 
      CHCH=CND(R,1)+CND(R,2)*POLEGCH2(CT,2)
      C60=CND(R,1)
c      write(60,*)r,C60

      END IF

      IF(N.EQ.8)THEN
      CHCH=0.0d0
      CHCH=CND(R,3)+CND(R,4)*POLEGCH2(CT,2)+CND(R,6)*POLEGCH2(CT,3)
      C80=CND(R,3)
c      write(80,*)r,C80
      END IF

      IF(N.EQ.10)THEN
      CHCH=0.0d0
      CHCH=CND(R,5)
      END IF

      RETURN
      END

c---------------------------------------------------------------------------
c---------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION CHIDCH2(N,RS,I)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A0(10),B0(10),R0(3)
C ** R0 A-BC  taken as (R0 AB + R0 AC)/2
      DATA R0(1),R0(2),R0(3)/7.4096D0,7.1688D0,7.1688475D0/
      DATA (A0(J), J=1,10)/16.3661D0,10.062652D0,7.5709297D0,
     1      6.1869941D0,5.29027135D0,4.6549655D0,
     2      4.17772554D0,3.8040565D0,3.50229909D0,3.2527089D0/
      DATA (B0(J), J=1,10)/15.623647D0,14.19721174D0,12.9010098D0,
     1      11.723150819D0,10.65283005D0,9.6802293D0,
     2      8.79642676D0,7.9933152D0,7.26352748D0,6.6003693D0/
      RHO=(11.0D0+2.5D0*R0(I))/2.0D0
      RR =RS/RHO
      AUX1=(1.0D0-EXP(-A0(N)*RR-B0(N)*RR*RR))**N
      CHIDCH2=AUX1
c      CHI3=AUX1
c      DO IC=2,N
c      CHI=CHI3*AUX1
c      ENDDO
c     write(*,*)RS,CHI3
      RETURN
      END

c--------------------------------------------------------------------------

       DOUBLE PRECISION FUNCTION POLEGCH2(RAD,N)
       
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION POL(3) 

       POL(1)=1.0d0
       POL(2)=0.5d0*(3.d0*RAD*RAD-1.0d0)
c       write(*,*)POL(2) 
       POL(3)=(1.0d0/8.0d0*(35.0d0*(RAD)**4
     & -30.0D0*(RAD)**2+3.0d0))

       POLEGCH2=POL(N)

       RETURN
       END

C*********************************************************
C   ROUTINE FOR INTENAL-JACOBI TRANSFORMATION
C   W1,W2,W3 ARE THE ATOMIC MASSES
C*********************************************************
      SUBROUTINE INT_JAC_CH2(R,RS,CCITA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=3)
      DIMENSION R(N),W(N),RS(N),CCITA(N),RPR(N)
      EXTERNAL RADIOCH2,ANGLECH2
      DATA (W(I),I=1,N)/1.00783D0,1.00783D0,12.010D0/
         RPR(1)=R(1)
         RPR(2)=R(2)
         RPR(3)=R(3)
      CM=W(1)/(W(1)+W(2))
      RS(1)=RADIOCH2(CM,RPR)
      CCITA(1)=ANGLECH2(CM,RPR)/(RS(1)*R(1))
      IF(CCITA(1).GT.1.0D0) CCITA(1)=1.0D0
      IF(CCITA(1).LT.-1.0D0) CCITA(1)=-1.0D0

      CALL PERMUTACH2(RPR)
      CM=W(2)/(W(2)+W(3))
      RS(2)=RADIOCH2(CM,RPR)
      CCITA(2)=ANGLECH2(CM,RPR)/(RS(2)*R(2))
      IF(CCITA(2).GT.1.0D0) CCITA(2)=1.0D0
      IF(CCITA(2).LT.-1.0D0) CCITA(2)=-1.0D0

      CALL PERMUTACH2(RPR)
      CM=W(3)/(W(3)+W(1))
      RS(3)=RADIOCH2(CM,RPR)
      CCITA(3)=ANGLECH2(CM,RPR)/(RS(3)*R(3))
      IF(CCITA(3).GT.1.0D0) CCITA(3)=1.0D0
      IF(CCITA(3).LT.-1.0D0) CCITA(3)=-1.0D0
      RETURN
      END

      DOUBLE PRECISION FUNCTION RADIOCH2(CM,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=3,ZERO=1.D-12)
      DIMENSION R(N)
      RADIOCH2=ZERO
      ARG=(CM*CM-CM)*R(1)*R(1)+
     &    (1.0D0-CM)*R(2)*R(2)+CM*R(3)*R(3)
      IF(ARG.GT.ZERO) RADIOCH2=DSQRT(ARG)
      RETURN
      END

      DOUBLE PRECISION FUNCTION ANGLECH2(CM,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=3)
      DIMENSION R(N)
      ANGLECH2=(CM-0.5D0)*R(1)*R(1)-0.5D0*R(2)*R(2)+
     &          0.5D0 *R(3)*R(3)
      RETURN
      END

      SUBROUTINE PERMUTACH2(R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=3)
      DIMENSION R(N)
      RAUX=R(1)
      R(1)=R(2)
      R(2)=R(3)
      R(3)=RAUX
      RETURN
      END

C####################################################################################
C CONVERT CARTEZIAN TO INTERPARTICLE DISTANCES (A.U)
C####################################################################################

      SUBROUTINE CART2INTER(X,R,NTA)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NA=100,NT=3*NA)
      DIMENSION X(NT),R(NT)
      K=0
      DO I=1,NTA,3
       DO J=I,NTA,3
        IF(I.NE.J)THEN
         K=K+1
         R(K)=SQRT((X(I)-X(J))**2+(X(I+1)-X(J+1))**2+
     1   (X(I+2)-X(J+2))**2)
        ENDIF
       ENDDO
      ENDDO
      RETURN
      END

C###############################################################################
C NUMERICAL GRADIENT OF THE POTENTIAL ENERGY FUNCTION IN CARTEZIAN COORDINATES
C THIS APPROACH USES RICHARDSON'S EXTRAPOLATION FORMULA
C###############################################################################

      SUBROUTINE DERIVE(X,VAL,NTA)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=100, NT=3*NA,PREC=1.0D-4)
      DIMENSION X(NT),VAL(NT)
      DO I=1,NTA
      VAL(I)=DVLP(X,PREC,I,NTA)
      ENDDO
      RETURN
      END

      FUNCTION DVLP(X,STEP,I,NTA)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=100, NT=3*NA)
      DIMENSION X(NT),XD(NT)
      DO IN=1,NTA
      XD(IN)=X(IN)
      ENDDO
      XD(I)=X(I)-STEP
      FMIN1=VLP(XD,NTA)
      XD(I)=X(I)+STEP
      FMAX1=VLP(XD,NTA)
      XD(I)=X(I)-2*STEP
      FMIN2=VLP(XD,NTA)
      XD(I)=X(I)+2*STEP
      FMAX2=VLP(XD,NTA)
      DVLP=(FMIN2-8.0*FMIN1+8.0*FMAX1-FMAX2)/(12.D0*STEP)
      RETURN
      END

C###############################################################################
C FUNCTION EVALUATION
C###############################################################################

      FUNCTION VLP(X,NTA)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=100,NT=3*NA)
      DIMENSION X(NT),R(NT)
      NR=NTA*(NTA/3-1)/6
      CALL FCNC2H2(NTA,X,VLP)
      RETURN
      END
            
C####################################################################################
C CH2 TRIPLET POTENTIAL ENERGY SURFACE
C CARTER, MILLS AND MURREL, MOL. PHYS. 39, 455-469 (1980)
C####################################################################################
C R1=CH IN A.U.
C R2=CH IN A.U.
C R3=HH IN A.U.
C POTENTIAL IN A.U. 
C#################################################################################### 

      SUBROUTINE CH2(RAU,POTEH)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3) :: RAU
      DOUBLE PRECISION :: R1, R2, R3
      DOUBLE PRECISION, DIMENSION(3) :: TWOBD
      DOUBLE PRECISION :: POTEH, THREEBD
      CALL CHPICMM(RAU(1),TWOBD(1))
      CALL CHPICMM(RAU(2),TWOBD(2))
      CALL HHCMM(RAU(3),TWOBD(3))
      CALL CH2CMM(RAU,THREEBD)
      POTEH=SUM(TWOBD)+THREEBD
      RETURN
      END
      
C####################################################################################
C CH2 TRIPLET THREE-BODY ENERGY - CARTER, MILLS AND MURREL 1980
C####################################################################################

      SUBROUTINE CH2CMM(RAU,POTEH)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3) :: RAU
      DOUBLE PRECISION :: BOHR2ANGS, EH2EV
      DOUBLE PRECISION :: POT, POTEH, VZERO
      DOUBLE PRECISION, DIMENSION(3) :: CPOL1
      DOUBLE PRECISION, DIMENSION(6) :: CPOL2
      DOUBLE PRECISION :: POL1, POL2, GAMMA, FAC
      DOUBLE PRECISION :: RHO1, RHO2, RHO3
      DOUBLE PRECISION :: R1, R2, R3

      BOHR2ANGS=0.529177208D+00

      EH2EV=27.211385D+00

      R1=RAU(1)*BOHR2ANGS
      R2=RAU(2)*BOHR2ANGS
      R3=RAU(3)*BOHR2ANGS

      VZERO=-0.445D+00

      RHO1=R1-1.088D+00
      RHO2=R2-1.088D+00
      RHO3=R3-2.018D+00

      CPOL1(1)=-0.511D+00
      CPOL1(2)=-0.511D+00
      CPOL1(3)=4.997D+00

      POL1=CPOL1(1)*RHO1+
     1     CPOL1(2)*RHO2+
     2     CPOL1(3)*RHO3

      CPOL2(1)=-22.113D+00
      CPOL2(2)=-22.113D+00
      CPOL2(3)=-32.801D+00
      CPOL2(4)=-42.790D+00
      CPOL2(5)=64.826D+00
      CPOL2(6)=64.826D+00

      POL2=CPOL2(1)*(RHO1**2)+
     1     CPOL2(2)*(RHO2**2)+
     2     CPOL2(3)*(RHO3**2)+
     3     CPOL2(4)*(RHO1*RHO2)+
     4     CPOL2(5)*(RHO1*RHO3)+
     5     CPOL2(6)*(RHO2*RHO3)

      GAMMA=6.00D+00

      FAC=2.00D+00*SQRT(3.0D+00)

      POT=VZERO*(1.00D+00+POL1+POL2)*
     1          (1.00D+00-TANH(GAMMA*(RHO1+RHO2+RHO3)/FAC))

      POTEH=POT/EH2EV

      RETURN
      END      
      
C####################################################################################
C CH PI DOUBLET POTENTIAL ENERGY CURVE - CARTER, MILLS AND MURREL 1980
C####################################################################################

      SUBROUTINE CHPICMM(RAU,POTEH)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(4) :: COEF 
      INTEGER :: I
      DOUBLE PRECISION :: R, RAU, RE, X
      DOUBLE PRECISION :: POT, POTEH, D
      DOUBLE PRECISION :: BOHR2ANGS, EH2EV

      BOHR2ANGS=0.529177208D+00

      EH2EV=27.211385D+00

      R=RAU*BOHR2ANGS

      D=3.648D+00
  
      RE=1.120D+00

      X=R-RE

      COEF(1)=1.000D+00
      COEF(2)=3.726D+00
      COEF(3)=3.086D+00
      COEF(4)=1.877D+00

      POT=-D*(COEF(1)+COEF(2)*X+
     1        COEF(3)*(X**2)+
     2        COEF(4)*(X**3))*EXP(-COEF(2)*X)

      POTEH=POT/EH2EV

      RETURN
      END
      
C####################################################################################
C HH POTENTIAL ENERGY CURVE - CARTER, MILLS AND MURREL 1980
C####################################################################################      

      SUBROUTINE HHCMM(RAU,POTEH)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(4) :: COEF 
      INTEGER :: I
      DOUBLE PRECISION :: R, RAU, RE, X
      DOUBLE PRECISION :: POT, POTEH, D
      DOUBLE PRECISION :: BOHR2ANGS, EH2EV

      BOHR2ANGS=0.529177208D+00

      EH2EV=27.211385D+00

      R=RAU*BOHR2ANGS

      D=4.750D+00
  
      RE=0.742D+00

      X=R-RE

      COEF(1)=1.000D+00
      COEF(2)=3.880D+00
      COEF(3)=3.741D+00
      COEF(4)=3.255D+00

      POT=-D*(COEF(1)+COEF(2)*X+
     1        COEF(3)*(X**2)+
     2        COEF(4)*(X**3))*EXP(-COEF(2)*X)

      POTEH=POT/EH2EV

      RETURN
      END      
