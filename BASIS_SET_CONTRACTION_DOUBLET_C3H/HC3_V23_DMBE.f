C###############################################################################
C HC3 GROUND-STATE POTENTIAL IN CARTESIAN COORDINATES
C###############################################################################

      SUBROUTINE FCN(N,X,POT)
C N IS 3X THE NUMBER OF ATOMS
C X IS THE INITIAL CARTESIAN COORDINATES
C POT IS THE POTENTIAL
      IMPLICIT NONE
      INTEGER, PARAMETER :: NATOM=4
      INTEGER, PARAMETER :: INTDIST=(NATOM*(NATOM-1))/2
      INTEGER, PARAMETER :: NATOMMAX=100
      DOUBLE PRECISION, DIMENSION(INTDIST) :: W
      INTEGER :: I, N
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: X, R
      DOUBLE PRECISION :: POT
      CALL CART2INTER(X,R,N)
      DO I=1,INTDIST
        W(I)=R(I)
      END DO
      CALL POTHC3(W,POT)
      RETURN
      END

C###############################################################################
C HC3 GROUND-STATE POTENTIAL IN CARTESIAN COORDINATES (A.U) AND DERIVATIVES
C###############################################################################

      SUBROUTINE FCNCG(N,X,F,G)
C N IS 3X THE NUMBER OF ATOMS
C X IS THE INITIAL CARTESIAN COORDINATES
C F IS THE POTENTIAL
C G IS THE GRADIENT
      IMPLICIT NONE
      INTEGER :: I, N
      INTEGER, PARAMETER :: NATOMMAX=100
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: X, G
      DOUBLE PRECISION :: POT, F
      CALL FCN(N,X,POT)
      F=POT
      DO I=1,N
       G(I)=0.0D+00
      END DO
      CALL DERIVE(X,G,N)
      RETURN
      END

C###############################################################################
C HC3 GROUND-STATE POTENTIAL WITH ONLY 2- AND 3-BODY CONTRIBUTIONS
C C2H POTENTIAL FROM JOSEPH AND VARANDAS, J. PHYS. CHEM. A 114, 2655-2664 (2010)
C WITH EXPERIMENTAL C2 SIGMA SINGLET DIATOMIC, CMRR, NOT PUBLISHED (2018)
C C3 POTENTIAL FROM ROCHA AND VARANDAS, CHEM. PHYS. LETT. 700, 36-43 (2018)
C WITH EXPERIMENTAL C2 PIU TRIPLET DIATOMIC, CMRR, NOT PUBLISHED (2018)
C AND ALSO WITH THREE-BODY TERMS FOR GROUND-STATE TRIPLET C3
C ROCHA AND VARANDAS 2019
C###############################################################################

      SUBROUTINE POTHC3(R,POTHC3TOT)
      IMPLICIT NONE
      INTEGER :: I
      INTEGER, PARAMETER :: NATOM=4
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0000000D+00
      INTEGER, PARAMETER :: INTDIST=NATOM*(NATOM-1)/2
      DOUBLE PRECISION, DIMENSION(INTDIST) :: R, W
      INTEGER :: NTWOBD, NTHREEBD, NFOURBD
      REAL :: BICO, RECU1, RECU2
      DOUBLE PRECISION :: POT2BDHC31, POT2BDHC32 
      DOUBLE PRECISION :: POT3BDHC31, POT3BDHC32
      DOUBLE PRECISION :: POTHC3TOT, POT4BDHC3 
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
      DOUBLE PRECISION, DIMENSION(3) :: WIPT
      DOUBLE PRECISION, DIMENSION(TIMAX) :: POT3BD1, POT3BD2
      INTEGER, DIMENSION(QIMAX) :: QA1, QA2, QA3, QA4
      DOUBLE PRECISION, DIMENSION(6) :: WIPQ   
      DOUBLE PRECISION, DIMENSION(QIMAX) :: POT4BD
      INTEGER, PARAMETER :: NDIM=2
      DOUBLE PRECISION, DIMENSION(NDIM,NDIM) :: V, VEC
      DOUBLE PRECISION, DIMENSION(NDIM) :: EIGVAL
      DOUBLE PRECISION :: RHO

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

      POTHC3TOT=0.000D+00
      POT2BDHC31=0.000D+00
      POT2BDHC32=0.000D+00
      POT3BDHC31=0.000D+00
      POT3BDHC32=0.000D+00
      POT4BDHC3=0.000D+00

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
C        WRITE(*,*) WIPD(I), DA1(I), DA2(I)
C JUST IN CASE
        IF (WIPD(I).EQ.ZERO) THEN 
          PRINT*, "ERROR IN THE TWO-BODY INDEXES"
          STOP
        ELSE
          IF  (DA1(I).EQ.1 .AND. DA2(I).EQ.2  
     1    .OR. DA1(I).EQ.1 .AND. DA2(I).EQ.3  
     2    .OR. DA1(I).EQ.1 .AND. DA2(I).EQ.4) THEN
C TWO-BODY TERMS FOR CH 
C TWO VECTORS ARE CALCULATED 
            CALL CHPI(WIPD(I),POT2BD1(I))
            CALL CHPI(WIPD(I),POT2BD2(I))
          ELSE 
C ONE VECTOR WITH ENERGIES FROM C2 SIGMA G PLUS SINGLET
C            CALL C2SIGMAAJCV(WIPD(I),POT2BD1(I))
            CALL CHIPR_C2SIGPS_EXP(WIPD(I),POT2BD1(I))
C THE OTHER WITH ENERGIES FROM C2 PI U TRIPLET
            CALL POTC2PIUCMRR(WIPD(I),VEHF(I),VDC(I),POT2BD2(I))
C            CALL CHIPR_C2PIUT_EXP(WIPD(I),POT2BD2(I))
          END IF
          POT2BDHC31=POT2BDHC31+POT2BD1(I)
          POT2BDHC32=POT2BDHC32+POT2BD2(I)
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
C        WRITE(*,*) TA1(I), TA2(I), TA3(I)
C        WRITE(*,*) WIPT(1), WIPT(2), WIPT(3)
C JUST IN CASE
        IF  (WIPT(1).EQ.ZERO 
     1  .OR. WIPT(2).EQ.ZERO 
     2  .OR. WIPT(3).EQ.ZERO) THEN 
          PRINT*, "ERROR IN THE THREE-BODY INDEXES"
          STOP
        ELSE
        IF (TA1(I).EQ.2 .AND. 
     1      TA2(I).EQ.3 .AND. 
     3      TA3(I).EQ.4) THEN 
C THREE-BODY TERMS FOR C3 SIGMA G PLUS SINGLET
          CALL THRBODYSPECC3(WIPT,POT3BD1(I))
C THREE-BODY TERMS FOR C3 TRIPLET GROUND-STATE
          CALL CHIPR_TRIAT_C3T(WIPT,POT3BD2(I))          
        ELSE 
C THREE-BODY TERMS FOR C2H SIGMA PLUS DOUBLET
          CALL THRBODYC2H(WIPT,POT3BD1(I))
          CALL THRBODYC2H(WIPT,POT3BD2(I))
        END IF 
          POT3BDHC31=POT3BDHC31+POT3BD1(I)
          POT3BDHC32=POT3BDHC32+POT3BD2(I)
        END IF
      END DO

C###############################################################################
C COMBINATORIAL ANALYSIS OF FOUR-BODY TERMS 
C###############################################################################

      MTC=.FALSE.

      QI=0

      DO I=1,FOURBD
        QCOMB(I)=0
      END DO 

  300 CALL NEXKSB(NATOM,FOURBD,QCOMB,MTC)

      QI=QI+1

C EACH TRETRATOMIC INDEX IS SAVED IN THE VECTORS QA1, QA2, QA3, QA4

      QA1(QI)=QCOMB(1) 

      QA2(QI)=QCOMB(2)

      QA3(QI)=QCOMB(3)

      QA4(QI)=QCOMB(4)

      IF (MTC) GOTO 300

C      DO I=1,NFOURBD
C        WRITE(*,*) QA1(I), QA2(I), QA3(I), QA4(I) 
C      END DO

      V(1,1)=POT2BDHC31+POT3BDHC32
      V(2,2)=POT2BDHC32+POT3BDHC31
      V(1,2)=0.0025D+00
      V(2,1)=V(1,2)

      RHO=1.0D-08

      CALL DIAG23BD(V,EIGVAL,VEC,2,RHO)

C      WRITE(*,*) (POT3BDHC31), (POT3BDHC32)

      POTHC3TOT=EIGVAL(1)

      RETURN
      END

C###############################################################################
C DIAGONALIZE A MATRIX A, OF WHICH ONLY LOWER TRIANGLE IS USED
C AND DESTROYED, USING THE GIVENS-HOUSHOLDER ALGORITHM.
C EIGENVALUES ARE RETURNED IN ALGEBRAIC ASCENDING ORDER IN ARRAY
C EIG THE EIGENVECTORS ARE RETURNED IN VEC.   
C
C PARAMETERS PASSED                   
C RHO IS THE UPPER LIMIT FOR OFF-DIAGONAL     
C NN IS THE SIZE OF THE MATRIX              
C###############################################################################

      SUBROUTINE DIAG23BD(A,EIG,VEC,NN,RHO)
      PARAMETER(NDIM=2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EIG(NDIM)
      DIMENSION W(NDIM),BETASQ(NDIM)
      DIMENSION GAMMA(NDIM),BETA(NDIM)
      DIMENSION A(NDIM,NDIM),VEC(NDIM,NDIM)
      DIMENSION P(NDIM),Q(NDIM),IPOSV(NDIM)
      DIMENSION IVPOS(NDIM),IORD(NDIM)
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
C
C         PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION
C
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
C
C         D IS FACTOR OF PROPORTIONALITY. COMPUTE AND SAVE W VECTOR
C
               DO 170 I=NR,N2
                  TEMP=D*A(I+2,NR)
                  W(I+1)=TEMP
                  A(I+2,NR)=TEMP
  170          CONTINUE
C
C         PREMULTIPLY VECTOR W BY MATRIX A TO OBTAIN P VECTOR.
C         SIMULTANEOUSLY ACCUMULATE DOT PRODUCT WP,(THE SCALAR K).
C
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
C
C         P VECTOR AND SCALAR K NOW STORED. NEXT COMPUTE Q VECTOR
C         AND FORM PAP MATRIX.
C
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
C
C         ADJOIN AN IDENTITY MATRIX TO BE POSTMULTIPLIED BY ROTATIONS
C
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
C
C         POSTMULTIPLY IDENTITY BY P-TRANSPOSE
C
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
C
C         TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT
C
      NPAS=NPAS+1
      IF (BETASQ(M).GT.RHOSQ) GOTO 410
  390 EIG(M+1)=GAMMA(M+1)+SUM
  400 BETA(M)=0.0D0
      BETASQ(M)=0.0D0
      M=M-1
      IF (M.EQ.0) GOTO 430
      IF (BETASQ(M).LE.RHOSQ) GOTO 390
  410 CONTINUE 
C
C         TAKE ROOT OF CORNER 2 BY 2 NEAREST TO LOWER DIAGONAL IN
C         VALUE AS ESTIMATE OF EIGENVALUE TO USE FOR SHIFT
C
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
C
C         INITIALIZE AUXILIARY TABLES REQUIRED FOR
C         REARANGING THE VECTORS
C
      DO 440 J=1,N
         IPOSV(J)=J
         IVPOS(J)=J
         IORD(J) = J
  440 CONTINUE
C
C         USE A TRANSPOSITON SORT TO ORDER THE EIGENVALUES
C
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
C
C         BACK TRANSFORM THE VECTORS OF THE TRIPLE DIAGONAL MATRIX
C
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

C####################################################################################
C C2H SIGMA PLUS DOUBLET POTENTIAL ENERGY SURFACE
C JOSEPH AND VARANDAS, J. PHYS. CHEM. A 114, 2655-2664 (2010)
C####################################################################################
C R1=CC IN A.U.
C R2=CH IN A.U.
C R3=CH IN A.U.
C POTENTIAL IN A.U. 
C####################################################################################
c  f77  code for C2H Potentials.

c       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c       DIMENSION R(3)

c    2.2851    4.2940    2.0089 !geometry of global minimum, with  V= -0.4210
c    2.4232    2.4085    2.4085 !geometry of T-shaped minimum, with  V=-0.3861

c       R(1)=2.2851d0
c       R(2)=4.2940d0
c       R(3)=2.0089d0

c       V=POTEN(R)
       
c       write(*,6)R(1),R(2),R(3),V
c6      format(3(1x,f13.8),3(2x,f13.8),2x,f13.8,2x,f6.2)

       
c       END
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

C###############################################################################
C C2 SIGMA G PLUS SINGLET POTENTIAL - AJCV CHEM. PHYS. LETT. 471, 315-321 (2009)
C###############################################################################

      SUBROUTINE C2SIGMAAJCV(W,VCC1)
      DOUBLE PRECISION :: V11, V22, V12, R, POT_1, VEHF_1
      DOUBLE PRECISION :: VDC_1, POT_2, VEHF_2
      DOUBLE PRECISION :: VDC_2, V12_A, V12_B, COUPLING, VCC1 
      DOUBLE PRECISION :: AN, BN 
      DOUBLE PRECISION :: W
      DOUBLE PRECISION, PARAMETER :: ANGS2BOHR=0.529177209097955D+00

      R=W*ANGS2BOHR

      CALL DPOCC1(R,VEHF_1,VDC_1,V11)

      POT_1=V11
 
      CALL DPOCC2(R,VEHF_2,VDC_2,V22)

      POT_2=V22

      CALL DCOUPLING(R,V12_A,V12_B,V12)

      COUPLING=V12
                                      
      VCC1=(POT_1+POT_2)/2.00D+00-
     *      SQRT(((POT_1-POT_2)/2.00D+00)**2+COUPLING**2)
  
      END SUBROUTINE C2SIGMAAJCV

C###############################################################################

      SUBROUTINE DPOCC1(R,VEHF_1,VDC_1,V11)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(4) :: A, B, DAMP, CN, DC
      INTEGER, DIMENSION(4) :: ID
      DOUBLE PRECISION, DIMENSION(6) :: COEF, POLI 
      INTEGER :: I
      DOUBLE PRECISION :: R, RE, X, R0, RHO, RHO_LINHA
      DOUBLE PRECISION :: VDC_1, D, GAMMA, VEHF_1, V11
      DOUBLE PRECISION :: AN, BN

      ID=(/5D+00,6D+00,8D+00,10D+00/)

      RE=1.2424D+00

      X=R-RE

      D=0.267457D+00

      GAMMA=3.19584D+00 

      COEF(1)=4.08283D+00
      COEF(2)=1.40751D+00
      COEF(3)=-3.09673D+00
      COEF(4)=-2.90848D+00
      COEF(5)=2.97621D+00
      COEF(6)=-0.60801D+00

      DO I=1,6
        POLI(I)=COEF(I)*X**I	
      END DO
 
      VEHF_1=-(D/R)*(1.0D+00+SUM(POLI))*EXP(-GAMMA*X)

      R0=4.17577394D+00

      RHO=(5.5D+00*0.5291772192D+00)+1.25D+00*R0

      RHO_LINHA=R/RHO
 
      CN(1)=-0.905048D+00
      CN(2)=1.73142D+00
      CN(3)=CN(2)*R0**(1.54)
      CN(4)=CN(2)*1.31*R0**(2*1.54)

      DO I=1,4
        A(I) = AN(ID(I))
      END DO

      DO I=1,4
        B(I) = BN(ID(I))
      END DO

      DO I=1,4
        DAMP(I)=(1.0D+00-EXP(-A(I)*RHO_LINHA-
     *           B(I)*RHO_LINHA**2))**INT(ID(I))
      END DO 

      DO I=1,4
        DC(I)=-CN(I)*DAMP(I)/R**INT(ID(I)) 
      END DO 
 
      VDC_1=SUM(DC)  

      V11=VEHF_1+VDC_1 

      END SUBROUTINE DPOCC1

C###############################################################################

      SUBROUTINE DPOCC2(R,VEHF_2,VDC_2,V22)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(4) :: A, B, DAMP, CN, DC
      INTEGER, DIMENSION(4) :: ID
      DOUBLE PRECISION, DIMENSION(7) :: COEF, POLI 
      INTEGER :: I
      DOUBLE PRECISION :: R, RE, X, R0, RHO, RHO_LINHA
      DOUBLE PRECISION :: VDC_2, D, GAMMA, VEHF_2, V22
      DOUBLE PRECISION :: AN, BN

      ID = (/5D+00,6D+00,8D+00,10D+00/)

      RE=1.3707D+00

      X=R-RE

      D=0.14217D+00

      GAMMA=5.36985D+00 

      COEF(1)=6.93459D+00
      COEF(2)=16.2156D+00
      COEF(3)=23.9559D+00
      COEF(4)=23.2267D+00
      COEF(5)=-2.5254D+00
      COEF(6)=-18.4952D+00
      COEF(7)=5.47728D+00

      DO I=1,7
        POLI(I)=COEF(I)*X**I	
      END DO
 
      VEHF_2 =-(D/R)*(1.0D+00+SUM(POLI))*EXP(-GAMMA*X) 

      R0=4.17577394D+00

      RHO=(5.5D+00*0.5291772192D+00)+1.25D+00*R0

      RHO_LINHA=R/RHO

      CN(1) = 0.0D+00
      CN(2) = 2.18932D+00
      CN(3) = CN(2)*R0**(1.54)
      CN(4) = CN(2)*1.31*R0**(2*1.54)

      DO I=1,4
       A(I) = AN(ID(I))
      END DO

      DO I=1,4
        B(I) = BN(ID(I))
      END DO

      DO I=1,4
        DAMP(I)=(1.0D+00-EXP(-A(I)*RHO_LINHA-
     *           B(I)*RHO_LINHA**2))**INT(ID(I))
      END DO 

      DO I=1,4
        DC(I)=-CN(I)*DAMP(I)/R**INT(ID(I)) 
      END DO 

      VDC_2=SUM(DC)

      V22=VEHF_2+VDC_2 

      END SUBROUTINE DPOCC2

C###############################################################################

      SUBROUTINE DCOUPLING(R,V12_A,V12_B,V12)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(5) :: A_linha, POLI_A
      DOUBLE PRECISION, DIMENSION(7) :: B_linha, POLI_B
      DOUBLE PRECISION :: V12_A, V12_B, V12, Rc, Rm 
      DOUBLE PRECISION :: GAMMA_A, GAMMA_B, S_tilde, R, S_2_tilde
      INTEGER :: I
      INTEGER, DIMENSION(5) :: INDEX_A = (/0,1,2,3,4/)
      INTEGER, DIMENSION(7) :: INDEX_B = (/0,1,2,3,4,5,6/)
      DOUBLE PRECISION :: AN, BN

      A_linha(1)=-0.00307997D+00
      A_linha(2)=0.0231469D+00
      A_linha(3)=0.0204273D+00
      A_linha(4)=-0.0693208D+00
      A_linha(5)=0.330466D+00 

      Rc=1.610505D+00
  
      GAMMA_A=5.36985D+00
 
      S_tilde=(R-Rc)
 
      DO I=1,5
        POLI_A(I) = A_linha(I)*S_tilde**INDEX_A(I)
      END DO

      V12_A=SUM(POLI_A)/COSH(GAMMA_A*S_tilde)

      B_linha(1)=-0.0000446415D+00
      B_linha(2)=0.0675524D+00
      B_linha(3)=0.114319D+00
      B_linha(4)=-0.433484D+00
      B_linha(5)=-0.360928D+00
      B_linha(6)=1.22331D+00
      B_linha(7)=-0.710691D+00

      Rm=1.2424D+00

      GAMMA_B=2.9832D+00

      S_2_tilde=(R-Rm)

      DO I=1,7
        POLI_B(I)=B_linha(I)*S_2_tilde**INDEX_B(I)
      END DO

      V12_B=SUM(POLI_B)*EXP(-GAMMA_B*S_2_tilde**2)
 
      V12=V12_A + V12_B 

      END SUBROUTINE DCOUPLING

C###############################################################################
C FINDING THREE BODY ENERGIES FOR C3 SIGMA G PLUS SINGLET
C###############################################################################

      SUBROUTINE THRBODYSPECC3(R,THREEBDSPEC)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3) :: R
      DOUBLE PRECISION :: SCALEDPOT23TOT
      DOUBLE PRECISION :: POTTOT
      DOUBLE PRECISION :: THREEBDSPEC
      DOUBLE PRECISION :: IFN

      CALL SCALEDDYNPOT23(R(1),R(2),R(3),SCALEDPOT23TOT)

      CALL POTC3SPEC(R,POTTOT)

      IFN=50.000D+00

C JUST IN CASE, SET BY BRUTE FORE THE THREE-BODY AS ZERO BEYOND 50 BORHS
      IF (R(1).GE.IFN .OR. R(2).GE.IFN .OR. R(3).GE.IFN) THEN 
        THREEBDSPEC=0.0D+00
      ELSE
        THREEBDSPEC=POTTOT-SCALEDPOT23TOT
      END IF

      END SUBROUTINE THRBODYSPECC3

C################################################################################
C# ENERGY SWITCHING GROUND-STATE GLOBAL PES OF SINGLET C3 - DMBE/ES/SS BETA=50  #
C################################################################################
C##### C.M.R. Rocha and A.J.C. Varandas, Chem. Phys. Lett. in press (2018) ######
C################################################################################
C                                                                              ##
C USE CALL POTC3SPEC(R,POT)                                                    ##
C                                                                              ##
C R IS AN ARRAY OF INPUT COORDINATES IN BOHR                                   ##
C                                                                              ## 
C OUTPUT ENERGIES IN HARTREE                                                   ##
C AND GIVEN WITH RESPECT TO INFINITELY SEPARATED C(3P)+C2(3PIu) FRAGMENTS      ##
C                                                                              ## 
C 1.) COORDINATES FOR THE DINFH GLOBAL MINIMUM (IN BOHR)                       ##
C                                                                              ##
C R(1)=2.445D+00                                                               ##
C                                                                              ##  
C R(2)=4.890D+00                                                               ##
C                                                                              ##
C R(3)=2.445D+00                                                               ##
C                                                                              ##
C ENERGY (IN HARTRE)=-0.290427D+00                                             ##
C                                                                              ##
C HARMONIC VIBRATIONAL FREQUENCIES (IN CM-1)                                   ##
C                                                                              ## 
C W1=1206.7                                                                    ##
C                                                                              ## 
C W2=42.8                                                                      ##
C                                                                              ##
C W3=2101.3                                                                    ##
C                                                                              ##
C 2.) COORDINATES FOR THE C2V ISOMERIZATION TRANSITION STATE (IN BOHR)      ##                                                               
C                                                                              ##
C R(1)=2.771D+00                                                               ##
C                                                                              ##  
C R(2)=2.401D+00                                                               ##
C                                                                              ##
C R(3)=2.771D+00                                                               ##
C                                                                              ##
C ENERGY (IN HARTRE)=-0.256353D+00                                             ##
C                                                                              ##
C HARMONIC VIBRATIONAL FREQUENCIES (IN CM-1)                                   ##
C                                                                              ## 
C W1=1295.2                                                                    ##
C                                                                              ## 
C W2=1840.5                                                                    ##
C                                                                              ##
C W3=1047.3i                                                                   ##
C                                                                              ## 
C 3.) COORDINATES FOR THE C2V VdW (LONG-RANGE) TRANSITION STATE (IN BOHR)      ##                                                               
C                                                                              ##
C R(1)=7.249D+00                                                               ##
C                                                                              ##  
C R(2)=7.249D+00                                                               ##
C                                                                              ##
C R(3)=2.470D+00                                                               ##
C                                                                              ##
C ENERGY (IN HARTRE)=-0.002558D+00                                             ##
C                                                                              ##
C HARMONIC VIBRATIONAL FREQUENCIES (IN CM-1)                                   ##
C                                                                              ## 
C W1=1618.1                                                                    ##
C                                                                              ## 
C W2=160.0                                                                     ##
C                                                                              ##
C W3=130.0i                                                                    ##
C                                                                              ##
C 4.) COORDINATES FOR THE C2V (LONG-RANGE) 2ND ORDER SADDLE POINT (IN BOHR)    ##                                                               
C                                                                              ##
C R(1)=5.511D+00                                                               ##
C                                                                              ##  
C R(2)=5.511D+00                                                               ##
C                                                                              ##
C R(3)=2.478D+00                                                               ##
C                                                                              ##
C ENERGY (IN HARTRE)=0.012504D+00                                              ##
C                                                                              ##
C HARMONIC VIBRATIONAL FREQUENCIES (IN CM-1)                                   ##
C                                                                              ## 
C W1=1345.3                                                                    ##
C                                                                              ## 
C W2=459.5i                                                                    ##
C                                                                              ##
C W3=546.5i                                                                    ##
C                                                                              ##
C INTEL AND GNU FORTRAN CONPILERS TESTED ON X86-64 FEDORA PLATFORM             ##
C                                                                              ##
C################################################################################

      SUBROUTINE POTC3SPEC(R,POT)
      IMPLICIT NONE
      INTEGER :: N
      DOUBLE PRECISION, PARAMETER :: CM2EH=219474.63137050000000D+00
      DOUBLE PRECISION, DIMENSION(3) :: R
      DOUBLE PRECISION :: V1, V2, POT
      DOUBLE PRECISION :: EMIN, EREF
      DOUBLE PRECISION :: DE, DEREF
      DOUBLE PRECISION :: ES, F
      DOUBLE PRECISION :: BETA
      DOUBLE PRECISION, PARAMETER :: ZERO=1.000D-12
      DOUBLE PRECISION :: POTZERO

      EMIN=-0.523221D+00 
C SUM OF DE OF DIATOMIC PLUS GROUND STATE ENERGY
      EREF=EMIN+(6250.0000D+00/CM2EH)

      CALL POTC3(R,V1)
      CALL SSC3(R,V2)

      IF (V1==EMIN) THEN
        DE=ZERO
      ELSE 
        DE=(V1-EMIN)
      END IF

      DEREF=(EREF-EMIN)

      BETA=10.0000D+00
      N=2

      F=ES(BETA,DE,DEREF,N)

      POTZERO=0.00000000000000000D+00

      POT=(V2+(V1-V2)*F)-POTZERO

      END SUBROUTINE POTC3SPEC

C###############################################################################
C C3 GROUND-STATE LOCAL SURFACE OF B. SCHRODER AND P. SEBALD
C The Journal of Chemical Physics 144, 044307 (2016); doi: 10.1063/1.4940780
C COORDINATES R1, R2 AND R3 IN BOHRS
C COEFFICIENTS ALSO IN A.U
C###############################################################################

      SUBROUTINE SSC3(R,LP)
      IMPLICIT NONE 
      INTEGER :: I, J, K, L, NPARAM
      INTEGER :: ORDERMIN, ORDERMAX
      DOUBLE PRECISION, DIMENSION(3)  :: R                                    
C BOHRS
      DOUBLE PRECISION, DIMENSION(3) :: VECDIST
C      DOUBLE PRECISION :: MINVECDIST, MAXVECDIST
      INTEGER, DIMENSION(1) :: POSITMIN, POSITMAX
      INTEGER :: POSIT
      DOUBLE PRECISION :: F1, F2, F3, THETA, COSTHETA, V, LP
      DOUBLE PRECISION :: DD1, DD2, DTHETA
      DOUBLE PRECISION, PARAMETER :: PI=3.14159265358979D+00
      DOUBLE PRECISION, PARAMETER :: DEG=180.0000000000000000000D+00
      DOUBLE PRECISION, PARAMETER :: ANGS2BOHR=0.529177209097955D+00
      DOUBLE PRECISION, PARAMETER :: RE=1.2939700000000000000000D+00          
C ANGS
      DOUBLE PRECISION, PARAMETER :: VREF=-0.523221D+00          
C SUM OF DE OF DIATOMIC PLUS GROUND STATE ENERGY - A.U
      DOUBLE PRECISION, PARAMETER :: ANGMAX=100.0000000000000000D+00          
C DEGS
      DOUBLE PRECISION, PARAMETER :: WALL=0.00000000000000000000D+00          
C A.U
      DOUBLE PRECISION, DIMENSION(369) :: W

      VECDIST(1)=R(1)
      VECDIST(2)=R(2)
      VECDIST(3)=R(3)
      POSITMIN=MINLOC(VECDIST)
      POSITMAX=MAXLOC(VECDIST)
      POSIT=6-(POSITMIN(1)+POSITMAX(1))
      IF (R(1)==R(2).AND.R(1)==R(3).AND.R(2)==R(3)) THEN
        F3=R(3)
        F1=R(1)
        F2=R(2)
      ELSE
        F3=VECDIST(POSITMAX(1))
        F1=VECDIST(POSITMIN(1))
        F2=VECDIST(POSIT)
      END IF
      COSTHETA=((F3**2-F1**2-F2**2)/(-2.00D+00*F1*F2))                       
C BRUTE-FORCE PERMUTATION - ALWAYS CHANNEL #3

      IF (COSTHETA>1.00D+00) THEN 
        COSTHETA=1.00D+00
      ELSE IF (COSTHETA<-1.00D+00) THEN 
        COSTHETA=-1.00D+00
      END IF 

      THETA=ACOS(COSTHETA)

      DD1=(F1-(RE/ANGS2BOHR))
      DD2=(F2-(RE/ANGS2BOHR))
      DTHETA=PI-THETA

      V=0.000000D+00
      LP=0.000000D+00
      ORDERMIN=2
      ORDERMAX=14
      NPARAM=0

      W(1  )=0.000413870D+00 !002
      W(2  )=0.332406930D+00 !020
      W(3  )=-0.00355710D+00 !110
      W(4  )=0.332406930D+00 !200
      W(5  )=-0.00449159D+00 !012
      W(6  )=-0.35060064D+00 !030
      W(7  )=-0.00449159D+00 !102
      W(8  )=-0.00614366D+00 !120
      W(9  )=-0.00614366D+00 !210
      W(10 )=-0.35060064D+00 !300
      W(11 )=0.000646540D+00 !004
      W(12 )=-0.00103148D+00 !022
      W(13 )=0.226902090D+00 !040
      W(14 )=-0.01222084D+00 !112
      W(15 )=-0.00075277D+00 !130
      W(16 )=-0.00103148D+00 !202
      W(17 )=0.000396800D+00 !220
      W(18 )=-0.00075277D+00 !310
      W(19 )=0.226902090D+00 !400
      W(20 )=-0.00042866D+00 !014
      W(21 )=0.001555570D+00 !032
      W(22 )=-0.11822982D+00 !050
      W(23 )=-0.00042866D+00 !104
      W(24 )=0.014602430D+00 !122
      W(25 )=0.002159950D+00 !140
      W(26 )=0.014602430D+00 !212
      W(27 )=0.005448950D+00 !230
      W(28 )=0.001555570D+00 !302
      W(29 )=0.005448950D+00 !320
      W(30 )=0.002159950D+00 !410
      W(31 )=-0.11822982D+00 !500
      W(32 )=0.000165470D+00 !006
      W(33 )=0.003274710D+00 !024
      W(34 )=0.015707020D+00 !042
      W(35 )=0.053588550D+00 !060
      W(36 )=-0.00263969D+00 !114
      W(37 )=-0.01948533D+00 !132
      W(38 )=-0.02856759D+00 !150
      W(39 )=0.003274710D+00 !204
      W(40 )=0.005767650D+00 !222
      W(41 )=0.001766880D+00 !240
      W(42 )=-0.01948533D+00 !312
      W(43 )=0.013937690D+00 !330
      W(44 )=0.015707020D+00 !402
      W(45 )=0.001766880D+00 !420
      W(46 )=-0.02856759D+00 !510
      W(47 )=0.053588550D+00 !600
      W(48 )=0.000455290D+00 !016
      W(49 )=-0.00321668D+00 !034
      W(50 )=0.000000000D+00 !052
      W(51 )=-0.02849987D+00 !070
      W(52 )=0.000455290D+00 !106
      W(53 )=-0.00004568D+00 !124
      W(54 )=0.000000000D+00 !142
      W(55 )=0.000000000D+00 !160
      W(56 )=-0.00004568D+00 !214
      W(57 )=0.000000000D+00 !232
      W(58 )=0.000000000D+00 !250
      W(59 )=-0.00321668D+00 !304
      W(60 )=0.000000000D+00 !322
      W(61 )=0.000000000D+00 !340
      W(62 )=0.000000000D+00 !412
      W(63 )=0.000000000D+00 !430
      W(64 )=0.000000000D+00 !502
      W(65 )=0.000000000D+00 !520
      W(66 )=0.000000000D+00 !610
      W(67 )=-0.02849987D+00 !700
      W(68 )=-0.00004851D+00 !008
      W(69 )=-0.00074881D+00 !026
      W(70 )=-0.06550030D+00 !044
      W(71 )=0.000000000D+00 !062
      W(72 )=0.020014950D+00 !080
      W(73 )=0.001621080D+00 !116
      W(74 )=0.069681790D+00 !134
      W(75 )=0.000000000D+00 !152
      W(76 )=0.000000000D+00 !170
      W(77 )=-0.00074881D+00 !206
      W(78 )=-0.05757693D+00 !224
      W(79 )=0.000000000D+00 !242
      W(80 )=0.000000000D+00 !260
      W(81 )=0.069681790D+00 !314
      W(82 )=0.000000000D+00 !332
      W(83 )=0.000000000D+00 !350
      W(84 )=-0.06550030D+00 !404
      W(85 )=0.000000000D+00 !422
      W(86 )=0.000000000D+00 !440
      W(87 )=0.000000000D+00 !512
      W(88 )=0.000000000D+00 !530
      W(89 )=0.000000000D+00 !602
      W(90 )=0.000000000D+00 !620
      W(91 )=0.000000000D+00 !710
      W(92 )=0.020014950D+00 !800
      W(93 )=-0.00020893D+00 !018
      W(94 )=0.005371820D+00 !036
      W(95 )=0.000000000D+00 !054 
      W(96 )=0.000000000D+00 !072
      W(97 )=0.000000000D+00 !090
      W(98 )=-0.00020893D+00 !108
      W(99 )=-0.00172969D+00 !126
      W(100)=0.000000000D+00 !144
      W(101)=0.000000000D+00 !162
      W(102)=0.000000000D+00 !180
      W(103)=-0.00172969D+00 !216
      W(104)=0.000000000D+00 !234
      W(105)=0.000000000D+00 !252
      W(106)=0.000000000D+00 !270
      W(107)=0.005371820D+00 !306
      W(108)=0.000000000D+00 !324
      W(109)=0.000000000D+00 !342
      W(110)=0.000000000D+00 !360
      W(111)=0.000000000D+00 !414
      W(112)=0.000000000D+00 !432
      W(113)=0.000000000D+00 !450
      W(114)=0.000000000D+00 !504
      W(115)=0.000000000D+00 !522
      W(116)=0.000000000D+00 !540
      W(117)=0.000000000D+00 !612
      W(118)=0.000000000D+00 !630
      W(119)=0.000000000D+00 !702
      W(120)=0.000000000D+00 !720
      W(121)=0.000000000D+00 !810
      W(122)=0.000000000D+00 !900
      W(123)=0.000106040D+00 !0010
      W(124)=-0.00157956D+00 !028
      W(125)=0.063146770D+00 !046
      W(126)=0.000000000D+00 !064
      W(127)=0.000000000D+00 !082
      W(128)=0.000000000D+00 !0100
      W(129)=0.000683530D+00 !118
      W(130)=-0.06693212D+00 !136
      W(131)=0.000000000D+00 !154
      W(132)=0.000000000D+00 !172 
      W(133)=0.000000000D+00 !190
      W(134)=-0.00157956D+00 !208
      W(135)=0.017050800D+00 !226
      W(136)=0.000000000D+00 !244
      W(137)=0.000000000D+00 !262
      W(138)=0.000000000D+00 !280
      W(139)=-0.06693212D+00 !316
      W(140)=0.000000000D+00 !334
      W(141)=0.000000000D+00 !352
      W(142)=0.000000000D+00 !370
      W(143)=0.063146770D+00 !406
      W(144)=0.000000000D+00 !424
      W(145)=0.000000000D+00 !442
      W(146)=0.000000000D+00 !460
      W(147)=0.000000000D+00 !514
      W(148)=0.000000000D+00 !532
      W(149)=0.000000000D+00 !550
      W(150)=0.000000000D+00 !604
      W(151)=0.000000000D+00 !622
      W(152)=0.000000000D+00 !640
      W(153)=0.000000000D+00 !712
      W(154)=0.000000000D+00 !730
      W(155)=0.000000000D+00 !802
      W(156)=0.000000000D+00 !820
      W(157)=0.000000000D+00 !910
      W(158)=0.000000000D+00 !1000
      W(159)=0.000000000D+00 !0110
      W(160)=-0.00220384D+00 !038
      W(161)=0.000000000D+00 !056
      W(162)=0.000000000D+00 !074
      W(163)=0.000000000D+00 !092
      W(164)=0.000000000D+00 !0110
      W(165)=0.000000000D+00 !1010
      W(166)=-0.00059965D+00 !128
      W(167)=0.000000000D+00 !146
      W(168)=0.000000000D+00 !164
      W(169)=0.000000000D+00 !182
      W(170)=0.000000000D+00 !1100
      W(171)=-0.00059965D+00 !218
      W(172)=0.000000000D+00 !236
      W(173)=0.000000000D+00 !254
      W(174)=0.000000000D+00 !272
      W(175)=0.000000000D+00 !290
      W(176)=-0.00220384D+00 !308
      W(177)=0.000000000D+00 !326
      W(178)=0.000000000D+00 !344
      W(179)=0.000000000D+00 !362
      W(180)=0.000000000D+00 !380
      W(181)=0.000000000D+00 !416
      W(182)=0.000000000D+00 !434
      W(183)=0.000000000D+00 !452
      W(184)=0.000000000D+00 !470
      W(185)=0.000000000D+00 !506
      W(186)=0.000000000D+00 !524
      W(187)=0.000000000D+00 !542
      W(188)=0.000000000D+00 !560
      W(189)=0.000000000D+00 !614
      W(190)=0.000000000D+00 !632
      W(191)=0.000000000D+00 !650
      W(192)=0.000000000D+00 !704
      W(193)=0.000000000D+00 !722
      W(194)=0.000000000D+00 !740
      W(195)=0.000000000D+00 !812
      W(196)=0.000000000D+00 !830
      W(197)=0.000000000D+00 !902
      W(198)=0.000000000D+00 !920
      W(199)=0.000000000D+00 !1010
      W(200)=0.000000000D+00 !1100
      W(201)=-0.00004245D+00 !0012
      W(202)=0.000000000D+00 !0210
      W(203)=-0.00935219D+00 !048
      W(204)=0.000000000D+00 !066
      W(205)=0.000000000D+00 !084
      W(206)=0.000000000D+00 !0102
      W(207)=0.000000000D+00 !0120
      W(208)=0.000000000D+00 !1110
      W(209)=0.001070210D+00 !138
      W(210)=0.000000000D+00 !156
      W(211)=0.000000000D+00 !174
      W(212)=0.000000000D+00 !192
      W(213)=0.000000000D+00 !1110
      W(214)=0.000000000D+00 !2010
      W(215)=0.045224820D+00 !228
      W(216)=0.000000000D+00 !246
      W(217)=0.000000000D+00 !264
      W(218)=0.000000000D+00 !282
      W(219)=0.000000000D+00 !2100
      W(220)=0.001070210D+00 !318
      W(221)=0.000000000D+00 !336
      W(222)=0.000000000D+00 !354
      W(223)=0.000000000D+00 !372
      W(224)=0.000000000D+00 !390
      W(225)=-0.00935219D+00 !408
      W(226)=0.000000000D+00 !426
      W(227)=0.000000000D+00 !444
      W(228)=0.000000000D+00 !462
      W(229)=0.000000000D+00 !480
      W(230)=0.000000000D+00 !516
      W(231)=0.000000000D+00 !534
      W(232)=0.000000000D+00 !552
      W(233)=0.000000000D+00 !570
      W(234)=0.000000000D+00 !606
      W(235)=0.000000000D+00 !624
      W(236)=0.000000000D+00 !642
      W(237)=0.000000000D+00 !660
      W(238)=0.000000000D+00 !714
      W(239)=0.000000000D+00 !732
      W(240)=0.000000000D+00 !750
      W(241)=0.000000000D+00 !804
      W(242)=0.000000000D+00 !822
      W(243)=0.000000000D+00 !840
      W(244)=0.000000000D+00 !912
      W(245)=0.000000000D+00 !930
      W(246)=0.000000000D+00 !1002
      W(247)=0.000000000D+00 !1020
      W(248)=0.000000000D+00 !1110
      W(249)=0.000000000D+00 !1200
      W(250)=0.000000000D+00 !0112
      W(251)=0.000000000D+00 !0310
      W(252)=0.000000000D+00 !058
      W(253)=0.000000000D+00 !076
      W(254)=0.000000000D+00 !094
      W(255)=0.000000000D+00 !0112
      W(256)=0.000000000D+00 !0130
      W(257)=0.000000000D+00 !1012
      W(258)=0.000000000D+00 !1210
      W(259)=0.000000000D+00 !148
      W(260)=0.000000000D+00 !166
      W(261)=0.000000000D+00 !184
      W(262)=0.000000000D+00 !1102
      W(263)=0.000000000D+00 !1120
      W(264)=0.000000000D+00 !2110
      W(265)=0.000000000D+00 !238
      W(266)=0.000000000D+00 !256
      W(267)=0.000000000D+00 !274
      W(268)=0.000000000D+00 !292
      W(269)=0.000000000D+00 !2110
      W(270)=0.000000000D+00 !3010
      W(271)=0.000000000D+00 !328
      W(272)=0.000000000D+00 !346
      W(273)=0.000000000D+00 !364
      W(274)=0.000000000D+00 !382
      W(275)=0.000000000D+00 !3100
      W(276)=0.000000000D+00 !418
      W(277)=0.000000000D+00 !436
      W(278)=0.000000000D+00 !454
      W(279)=0.000000000D+00 !472
      W(280)=0.000000000D+00 !490
      W(281)=0.000000000D+00 !508
      W(282)=0.000000000D+00 !526
      W(283)=0.000000000D+00 !544
      W(284)=0.000000000D+00 !562
      W(285)=0.000000000D+00 !580
      W(286)=0.000000000D+00 !616
      W(287)=0.000000000D+00 !634
      W(288)=0.000000000D+00 !652
      W(289)=0.000000000D+00 !670
      W(290)=0.000000000D+00 !706
      W(291)=0.000000000D+00 !724
      W(292)=0.000000000D+00 !742
      W(293)=0.000000000D+00 !760
      W(294)=0.000000000D+00 !814
      W(295)=0.000000000D+00 !832
      W(296)=0.000000000D+00 !850
      W(297)=0.000000000D+00 !904
      W(298)=0.000000000D+00 !922
      W(299)=0.000000000D+00 !940
      W(300)=0.000000000D+00 !1012
      W(301)=0.000000000D+00 !1030
      W(302)=0.000000000D+00 !1102
      W(303)=0.000000000D+00 !1120
      W(304)=0.000000000D+00 !1210
      W(305)=0.000000000D+00 !1300
      W(306)=0.000006440D+00 !0014
      W(307)=0.000000000D+00 !0212
      W(308)=0.000000000D+00 !0410
      W(309)=0.000000000D+00 !068
      W(310)=0.000000000D+00 !086
      W(311)=0.000000000D+00 !0104
      W(312)=0.000000000D+00 !0122
      W(313)=0.000000000D+00 !0140
      W(314)=0.000000000D+00 !1112
      W(315)=0.000000000D+00 !1310
      W(316)=0.000000000D+00 !158
      W(317)=0.000000000D+00 !176
      W(318)=0.000000000D+00 !194
      W(319)=0.000000000D+00 !1112
      W(320)=0.000000000D+00 !1130
      W(321)=0.000000000D+00 !2012
      W(322)=0.000000000D+00 !2210
      W(323)=0.000000000D+00 !248
      W(324)=0.000000000D+00 !266
      W(325)=0.000000000D+00 !284
      W(326)=0.000000000D+00 !2102
      W(327)=0.000000000D+00 !2120
      W(328)=0.000000000D+00 !3110
      W(329)=0.000000000D+00 !338
      W(330)=0.000000000D+00 !356
      W(331)=0.000000000D+00 !374
      W(332)=0.000000000D+00 !392
      W(333)=0.000000000D+00 !3110
      W(334)=0.000000000D+00 !4010
      W(335)=0.000000000D+00 !428
      W(336)=0.000000000D+00 !446
      W(337)=0.000000000D+00 !464
      W(338)=0.000000000D+00 !482
      W(339)=0.000000000D+00 !4100
      W(340)=0.000000000D+00 !518
      W(341)=0.000000000D+00 !536
      W(342)=0.000000000D+00 !554
      W(343)=0.000000000D+00 !572
      W(344)=0.000000000D+00 !590
      W(345)=0.000000000D+00 !608
      W(346)=0.000000000D+00 !626
      W(347)=0.000000000D+00 !644
      W(348)=0.000000000D+00 !662
      W(349)=0.000000000D+00 !680
      W(350)=0.000000000D+00 !716
      W(351)=0.000000000D+00 !734
      W(352)=0.000000000D+00 !752
      W(353)=0.000000000D+00 !770
      W(354)=0.000000000D+00 !806
      W(355)=0.000000000D+00 !824
      W(356)=0.000000000D+00 !842
      W(357)=0.000000000D+00 !860
      W(358)=0.000000000D+00 !914
      W(359)=0.000000000D+00 !932
      W(360)=0.000000000D+00 !950
      W(361)=0.000000000D+00 !1004
      W(362)=0.000000000D+00 !1022
      W(363)=0.000000000D+00 !1040
      W(364)=0.000000000D+00 !1112
      W(365)=0.000000000D+00 !1130
      W(366)=0.000000000D+00 !1202
      W(367)=0.000000000D+00 !1220
      W(368)=0.000000000D+00 !1310
      W(369)=0.000000000D+00 !1400

      DO L=ORDERMIN,ORDERMAX
        DO I=0,L
          DO J=0,(L-I)
            K=L-I-J
              IF (MOD(K,2)==0) THEN
                NPARAM=NPARAM+1
                V=V+W(NPARAM)*(DD1**(I))*(DD2**(J))*(DTHETA**(K))
              END IF
          END DO
        END DO     
      END DO

      IF (DTHETA>(ANGMAX*PI/DEG)) THEN
        LP=V+VREF                                                 
C LP=WALL FOR DVR VIBRATIONAL CALCULATIONS
      ELSE
        LP=V+VREF
      END IF

      END SUBROUTINE SSC3

C####################################################################################
C C3 GROUND-STATE GLOBAL SURFACE OF C. M. R. ROCHA. AND A. J. C. VARANDAS           #
C COORDINATES R1, R2 AND R3 IN BOHRS                                                #
C ENERGY IN A.U                                                                     #
C####################################################################################

      SUBROUTINE POTC3(R,GP)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3) :: R
      DOUBLE PRECISION :: EHFD3H, EHFC2V, EHFTOT 
      DOUBLE PRECISION :: SCALEDPOT23, GP
      DOUBLE PRECISION, PARAMETER :: DD=0.0000000000000000000000000D+00
C     DOUBLE PRECISION, PARAMETER :: DD=0.2323414039409390507984199D+00
      DOUBLE PRECISION :: ONEBODY

      CALL SCALEDDYNPOT23(R(1),R(2),R(3),SCALEDPOT23)

      CALL THRBODYTOT(R(1),R(2),R(3),EHFD3H,EHFC2V,EHFTOT)

      ONEBODY=DD

      GP=SCALEDPOT23+EHFTOT+ONEBODY

      END SUBROUTINE POTC3

C####################################################################################

      SUBROUTINE THRBODYTOT(D1,D2,D3,EHFD3H,EHFC2V,EHFTOT)
      IMPLICIT NONE
      DOUBLE PRECISION :: D1, D2, D3
      DOUBLE PRECISION, DIMENSION(144) :: C 
      DOUBLE PRECISION :: EHFC2V1, EHFC2V2
      DOUBLE PRECISION :: EHFD3H, EHFC2V, EHFTOT

       C(  1)=  0.939419642867667104D-01
       C(  2)=  0.814094356965121627D-02
       C(  3)= -0.137991765583622972D+00
       C(  4)= -0.120752320819536626D-01
       C(  5)= -0.526762995195253561D-01
       C(  6)=  0.344140515403429492D+00
       C(  7)= -0.144177778289186331D-01
       C(  8)= -0.216556741282060983D+00
       C(  9)= -0.118516399062917804D+00
       C( 10)=  0.194334439638676221D+00
       C( 11)= -0.672855804010018006D-02
       C( 12)=  0.445640961467239749D+00
       C( 13)=  0.210270114680647513D+01
       C( 14)= -0.473015437845383593D-01
       C( 15)=  0.424675612328713759D-01
       C( 16)= -0.142313896757015378D-02
       C( 17)=  0.922426072159836935D-02
       C( 18)=  0.779262046319446761D-01
       C( 19)= -0.430756618506507479D-03
       C( 20)=  0.101106392051225646D+01
       C( 21)=  0.168776997019681810D-01
       C( 22)=  0.310076807579411050D-01
       C( 23)= -0.555669706556694617D-03
       C( 24)=  0.189073348055313062D+00
       C( 25)=  0.531839785809717306D-02
       C( 26)=  0.659753281861308993D+00
       C( 27)= -0.405365137471814449D-01
       C( 28)=  0.187354019661505999D+00
       C( 29)=  0.813305566580927765D-02
       C( 30)=  0.782028657391569890D-02
       C( 31)= -0.102320710237646184D-03
       C( 32)=  0.139389147527451004D-02
       C( 33)=  0.551339890164715286D-02
       C( 34)=  0.871095286063372325D-02
       C( 35)=  0.142558452764864302D-02
       C( 36)=  0.982099724137223407D-01
       C( 37)=  0.766319032385726507D-02
       C( 38)=  0.253278266783274787D-01
       C( 39)=  0.815312358197797056D-03
       C( 40)=  0.920746721066789567D-03
       C( 41)=  0.437428583519401120D-04
       C( 42)= -0.381045127804672845D-04
       C( 43)=  0.175798202807236079D-02
       C( 44)= -0.466517441974516019D-03
       C( 45)=  0.351449350549256755D-02
       C( 46)= -0.235648639269382502D-02
       C( 47)=  0.672025193130964838D-04
       C( 48)= -0.212710137098015605D-02
       C( 49)=  0.215021725683087460D-03
       C( 50)=  0.134146357895935923D-03
       C( 51)=  0.641835269147592122D-04
       C( 52)= -0.113967961247292232D-04
       C( 53)=  0.102605859070874463D-04
       C( 54)=  0.552705799921816565D-01
       C( 55)=  0.867056281944996593D-01
       C( 56)= -0.292061509407157527D+00
       C( 57)=  0.340031286778664996D-01
       C( 58)=  0.127723520264307916D+00
       C( 59)=  0.126702812674413323D+01
       C( 60)= -0.153752772318713638D-01
       C( 61)= -0.253399962637563180D-02
       C( 62)= -0.121847965474670805D+00
       C( 63)=  0.762417017494457006D+00
       C( 64)= -0.102626123630331217D-01
       C( 65)=  0.439870788773157273D+00
       C( 66)=  0.168735339849976040D+01
       C( 67)= -0.593511612131693914D-01
       C( 68)=  0.207750399009005826D+00
       C( 69)= -0.884966701418230949D-03
       C( 70)=  0.722406868036933741D-02
       C( 71)=  0.380312893499337584D-01
       C( 72)=  0.340641951855910896D-01
       C( 73)=  0.535047684436933180D+00
       C( 74)=  0.262522419438343813D-01
       C( 75)=  0.643720125996914139D-01
       C( 76)=  0.101519265922605704D-02
       C( 77)=  0.341658750722979229D-01
       C( 78)=  0.149953120015734655D-02
       C( 79)=  0.108240178322483449D+00
       C( 80)= -0.157367215202812889D-01
       C( 81)=  0.438820536029984748D-01
       C( 82)=  0.583518008253752166D-02
       C( 83)=  0.679892780345060874D-02
       C( 84)=  0.138470992549592043D-03
       C( 85)=  0.326796638933332995D+01
       C( 86)=  0.609999999999999987D+00
       C( 87)= -0.818307840867775037D-03
       C( 88)= -0.763963412105609719D-03
       C( 89)= -0.585126332145601002D-03
       C( 90)= -0.254858222399380373D-02
       C( 91)= -0.507353537340084482D-01
       C( 92)=  0.246576863041835619D-01
       C( 93)= -0.891309408429813476D-03
       C( 94)=  0.850094949972518971D-01
       C( 95)=  0.400507858455017435D-01
       C( 96)=  0.272964768595807573D-01
       C( 97)= -0.901669218460573541D-03
       C( 98)= -0.117434181088866413D+00
       C( 99)= -0.157822962081838825D+00
       C(100)= -0.144728632687294745D+01
       C(101)= -0.354159905365007022D+00
       C(102)= -0.278247474532728756D+01
       C(103)= -0.221459200604135731D+01
       C(104)= -0.107676338335859406D+00
       C(105)= -0.136264888251007368D-01
       C(106)= -0.160230674075789725D-01
       C(107)=  0.789401227026172853D-01
       C(108)= -0.469425271984457021D-01
       C(109)=  0.892272218467323429D-02
       C(110)=  0.354883893697421487D-01
       C(111)= -0.139408714066350751D-01
       C(112)=  0.260960999999999999D+01
       C(113)=  0.261999637999999990D+01
       C(114)=  0.260440127999999982D+01
       C(115)=  0.329751996671108705D+01
       C(116)=  0.323232707276617127D-05
       C(117)=  0.486143040388596260D-04
       C(118)= -0.212510693796359815D-02
       C(119)=  0.213367046321904380D-03
       C(120)= -0.128324510130825532D-02
       C(121)= -0.573522106229466290D-02
       C(122)=  0.357905162584196482D-03
       C(123)=  0.197085359974491720D-02
       C(124)= -0.178597719776261250D-02
       C(125)=  0.215027561383702992D-03
       C(126)=  0.219562598118572434D-03
       C(127)= -0.135147884296237010D-01
       C(128)= -0.541527773071893911D-01
       C(129)=  0.197689938834732071D+00
       C(130)= -0.851429167298263928D-01
       C(131)=  0.323164050737253100D-01
       C(132)=  0.288432405228098654D+00
       C(133)= -0.478502565254514192D-01
       C(134)=  0.512262982168152151D-02
       C(135)=  0.117985605912901449D-01
       C(136)= -0.244890509588525018D-02
       C(137)=  0.745602598096316554D-02
       C(138)=  0.641217471119934646D-04
       C(139)= -0.782019372756585958D-02
       C(140)= -0.534177225704332440D-03
       C(141)=  0.331976404999999986D+01
       C(142)=  0.331120610000000015D+01
       C(143)=  0.332404301999999996D+01
       C(144)=  0.430430404981779002D+01

      CALL THRBODYD3H(D1,D2,D3,C,EHFD3H)

      CALL THRBODYC2V1(D1,D2,D3,C,EHFC2V1)

      CALL THRBODYC2V2(D1,D2,D3,C,EHFC2V2)

      EHFC2V=EHFC2V1+EHFC2V2
 
      EHFTOT=EHFD3H+EHFC2V

      END SUBROUTINE THRBODYTOT

C####################################################################################

      SUBROUTINE THRBODYD3H(D1,D2,D3,C,EHFD3H)
      IMPLICIT NONE
      DOUBLE PRECISION :: D1, D2, D3
      INTEGER :: I, J, K, L, NUM
      DOUBLE PRECISION, DIMENSION(144) :: C
      DOUBLE PRECISION, DIMENSION(3) :: R  
      DOUBLE PRECISION, DIMENSION(3) :: Q 
      DOUBLE PRECISION, DIMENSION(3,3) :: A
      DOUBLE PRECISION :: TAU1, TAU2, TAU3
      DOUBLE PRECISION :: G1P1, G2P1, G3P1
      DOUBLE PRECISION :: G1P2, G2P2, G3P2
      INTEGER, PARAMETER :: ORDER1=9 
      INTEGER, PARAMETER :: ORDER2=7 
      DOUBLE PRECISION :: RDISP 
      DOUBLE PRECISION :: GAM 
      DOUBLE PRECISION :: POL1, POL2, EHFD3H 
      DOUBLE PRECISION :: DISP, RANG
      DOUBLE PRECISION :: LIN, RANGJT

      RDISP=C(85)
      GAM=C(86)
 
      R(1)=D1
      R(2)=D2
      R(3)=D3
      
      A(1,1)=SQRT(1.000D+00/3.00D+00)
      A(1,2)=A(1,1)
      A(1,3)=A(1,1)
      A(2,1)=0.000D+00
      A(2,2)=SQRT(1.000D+00/2.000D+00)
      A(2,3)=-A(2,2)
      A(3,1)=SQRT(2.000D+00/3.000D+00)
      A(3,2)=-SQRT(1.000D+00/6.000D+00)
      A(3,3)=A(3,2)
      
      NUM=0  
      POL1=0.000D+00
      POL2=0.000D+00
      EHFD3H=0.000D+00
        
      DO I=1,3
        Q(I)=0.00D+00
          DO J=1,3
            Q(I)=Q(I)+A(I,J)*DISP(R(J),RDISP)
          END DO
      END DO     
          
      TAU1=Q(1)
      TAU2=Q(2)**2+Q(3)**2
      TAU3=Q(3)**3-3.00D+00*Q(3)*Q(2)**2 
               
      DO L=0,ORDER1
        DO I=0,L
         IF (I==0) THEN
           G1P1=1.0000D+00
         ELSE 
           G1P1=TAU1**I
         END IF
         DO J=0,(L-I),2
           IF (J==0) THEN
             G2P1=1.0000D+00
           ELSE
             G2P1=TAU2**(J/2)
           END IF
           K=L-I-J
           IF (MOD(K,3)==0) THEN
             IF (K==0) THEN
               G3P1=1.0000D+00
             ELSE
               G3P1=TAU3**(K/3)
             END IF
             NUM=NUM+1
             POL1=POL1+C(NUM)*G1P1*G2P1*G3P1 
           END IF
         END DO
        END DO
      END DO

      NUM=53
  
      DO L=0,ORDER2
        DO I=0,L
          IF (I==0) THEN
            G1P2=1.0000D+00
          ELSE 
            G1P2=TAU1**I
          END IF
          DO J=0,(L-I),2
            IF (J==0) THEN
              G2P2=1.0000D+00
            ELSE
              G2P2=TAU2**(J/2)
            END IF
            K=L-I-J
            IF (MOD(K,3)==0) THEN
              IF (K==0) THEN
                G3P2=1.0000D+00
              ELSE
                G3P2=TAU3**(K/3)
              END IF
              NUM=NUM+1
              POL2=POL2+C(NUM)*G1P2*G2P2*G3P2 
            END IF
          END DO
        END DO
      END DO

      LIN=SQRT(TAU2)*RANGJT(D1,D2,D3)

      EHFD3H=(POL1-LIN*POL2)*RANG(D1,D2,D3,RDISP,GAM)
             
      END SUBROUTINE THRBODYD3H

C####################################################################################

      SUBROUTINE THRBODYC2V1(D1,D2,D3,C,EHFC2V1)
      IMPLICIT NONE
      DOUBLE PRECISION :: D1, D2, D3
      INTEGER :: I, J, K, L, NUMC2V
      DOUBLE PRECISION, DIMENSION(144) :: C
      DOUBLE PRECISION, DIMENSION(3) :: RC2V  
      DOUBLE PRECISION, DIMENSION(3) :: QC2V 
      DOUBLE PRECISION, DIMENSION(3,3) :: A
      DOUBLE PRECISION :: TAU1C2V, TAU2C2V, TAU3C2V
      DOUBLE PRECISION :: G1P1C2V, G2P1C2V, G3P1C2V
      DOUBLE PRECISION :: G1P2C2V, G2P2C2V, G3P2C2V
      DOUBLE PRECISION :: G1P3C2V, G2P3C2V, G3P3C2V
      INTEGER, PARAMETER :: ORDER1C2V=4 
      INTEGER, PARAMETER :: ORDER2C2V=3
      INTEGER, PARAMETER :: ORDER3C2V=3
      DOUBLE PRECISION :: POL1C2V, POL2C2V, POL3C2V, EHFC2V1 
      DOUBLE PRECISION :: DISP, RDISPC2V, X1, X2, RANGC2V
      DOUBLE PRECISION :: GAM1
      DOUBLE PRECISION :: S1, S2, S3
      DOUBLE PRECISION :: DELTAS,  Q31, RHO0, Q1ABS, Q2ABS, Q3ABS
      DOUBLE PRECISION, DIMENSION(3) :: DELTA, THETA
      DOUBLE PRECISION, DIMENSION(3) :: Q3SEAM, Q2SEAM
      DOUBLE PRECISION :: PI, ZERO
      INTEGER :: PHASE
      DOUBLE PRECISION :: LIN, RANGJT

      RDISPC2V=C(112)
      X1=C(113)
      X2=C(114)
      GAM1=C(115)

      RC2V(1)=D1
      RC2V(2)=D2
      RC2V(3)=D3
      
      A(1,1)=SQRT(1.000D+00/3.00D+00)
      A(1,2)=A(1,1)
      A(1,3)=A(1,1)
      A(2,1)=0.000D+00
      A(2,2)=SQRT(1.000D+00/2.000D+00)
      A(2,3)=-A(2,2)
      A(3,1)=SQRT(2.000D+00/3.000D+00)
      A(3,2)=-SQRT(1.000D+00/6.000D+00)
      A(3,3)=A(3,2)
      
      NUMC2V=86
      POL1C2V=0.000D+00
      POL2C2V=0.000D+00
      POL3C2V=0.000D+00
      EHFC2V1=0.000D+00
        
      DO I=1,3
        QC2V(I)=0.00D+00
          DO J=1,3
            QC2V(I)=QC2V(I)+A(I,J)*DISP(RC2V(J),RDISPC2V)
          END DO
      END DO   
          
      TAU1C2V=QC2V(1)
      TAU2C2V=QC2V(2)**2+QC2V(3)**2
      TAU3C2V=QC2V(3)**3-3.00D+00*QC2V(3)*QC2V(2)**2 
               
      DO L=0,ORDER1C2V
        DO I=0,L
         IF (I==0) THEN
           G1P1C2V=1.0000D+00
         ELSE 
           G1P1C2V=TAU1C2V**I
         END IF
         DO J=0,(L-I),2
           IF (J==0) THEN
             G2P1C2V=1.0000D+00
           ELSE
             G2P1C2V=TAU2C2V**(J/2)
           END IF
           K=L-I-J
           IF (MOD(K,3)==0) THEN
             IF (K==0) THEN
               G3P1C2V=1.0000D+00
             ELSE
               G3P1C2V=TAU3C2V**(K/3)
             END IF
             NUMC2V=NUMC2V+1
             POL1C2V=POL1C2V+C(NUMC2V)*G1P1C2V*G2P1C2V*G3P1C2V 
           END IF
         END DO
        END DO
      END DO 

      NUMC2V=97
  
      DO L=0,ORDER2C2V
        DO I=0,L
          IF (I==0) THEN
            G1P2C2V=1.0000D+00
          ELSE 
            G1P2C2V=TAU1C2V**I
          END IF
          DO J=0,(L-I),2
            IF (J==0) THEN
              G2P2C2V=1.0000D+00
            ELSE
              G2P2C2V=TAU2C2V**(J/2)
            END IF
            K=L-I-J
            IF (MOD(K,3)==0) THEN
              IF (K==0) THEN
                G3P2C2V=1.0000D+00
              ELSE
                G3P2C2V=TAU3C2V**(K/3)
              END IF
              NUMC2V=NUMC2V+1
              POL2C2V=POL2C2V+C(NUMC2V)*G1P2C2V*G2P2C2V*G3P2C2V 
            END IF
          END DO
        END DO
      END DO

      NUMC2V=104

      DO L=0,ORDER3C2V
        DO I=0,L
          IF (I==0) THEN
            G1P3C2V=1.0000D+00
          ELSE 
            G1P3C2V=TAU1C2V**I
          END IF
          DO J=0,(L-I),2
            IF (J==0) THEN
              G2P3C2V=1.0000D+00
            ELSE
              G2P3C2V=TAU2C2V**(J/2)
            END IF
            K=L-I-J
            IF (MOD(K,3)==0) THEN
              IF (K==0) THEN
                G3P3C2V=1.0000D+00
              ELSE
                G3P3C2V=TAU3C2V**(K/3)
              END IF
              NUMC2V=NUMC2V+1
              POL3C2V=POL3C2V+C(NUMC2V)*G1P3C2V*G2P3C2V*G3P3C2V 
            END IF
          END DO
        END DO
      END DO

      Q1ABS=A(1,1)*(D1+D2+D3)
      Q2ABS=A(2,2)*D2+A(2,3)*D3
      Q3ABS=A(3,1)*D1+A(3,2)*D2+A(3,3)*D3

      PI=3.1415926535897932D+00
 
      ZERO=0.00000000000000D+00

      THETA=(/(PI/2.00D+00),(7.00D+00*PI/6.00D+00),
     &        (11.00D+00*PI/6.00D+00)/)

      RHO0=ABS(Q31(Q1ABS))

      IF (Q31(Q1ABS)>=ZERO) THEN
        PHASE=0
      ELSE IF (Q31(Q1ABS)<ZERO) THEN
        PHASE=1
      ENDIF

      Q3SEAM(1)=RHO0*SIN(THETA(1)+DBLE(PHASE)*PI)
      Q2SEAM(1)=RHO0*COS(THETA(1)+DBLE(PHASE)*PI)
      Q3SEAM(2)=RHO0*SIN(THETA(2)+DBLE(PHASE)*PI)
      Q2SEAM(2)=RHO0*COS(THETA(2)+DBLE(PHASE)*PI)
      Q3SEAM(3)=RHO0*SIN(THETA(3)+DBLE(PHASE)*PI)
      Q2SEAM(3)=RHO0*COS(THETA(3)+DBLE(PHASE)*PI)

      DELTA(1)=SQRT((Q2ABS-Q2SEAM(1))**2+
     &              (Q3ABS-Q3SEAM(1))**2)
      DELTA(2)=SQRT((Q2ABS-Q2SEAM(2))**2+
     &              (Q3ABS-Q3SEAM(2))**2)
      DELTA(3)=SQRT((Q2ABS-Q2SEAM(3))**2+
     &              (Q3ABS-Q3SEAM(3))**2)

      S1=A(1,1)*(DELTA(1)+DELTA(2)+DELTA(3))
      S2=A(2,2)*DELTA(2)+A(2,3)*DELTA(3)
      S3=A(3,1)*DELTA(1)+A(3,2)*DELTA(2)+A(3,3)*DELTA(3)

      DELTAS=SQRT(S2**2+S3**2)

      LIN=SQRT(TAU2C2V)*RANGJT(D1,D2,D3)

      EHFC2V1=(POL1C2V-DELTAS*POL2C2V-
     &         LIN*POL3C2V)*
     &         RANGC2V(D1,D2,D3,X1,X2,GAM1)

      END SUBROUTINE THRBODYC2V1

C####################################################################################

      SUBROUTINE THRBODYC2V2(D1,D2,D3,C,EHFC2V2)
      IMPLICIT NONE
      DOUBLE PRECISION :: D1, D2, D3
      INTEGER :: I, J, K, L, NUMC2V
      DOUBLE PRECISION, DIMENSION(144) :: C
      DOUBLE PRECISION, DIMENSION(3) :: RC2V  
      DOUBLE PRECISION, DIMENSION(3) :: QC2V 
      DOUBLE PRECISION, DIMENSION(3,3) :: A
      DOUBLE PRECISION :: TAU1C2V, TAU2C2V, TAU3C2V
      DOUBLE PRECISION :: G1P1C2V, G2P1C2V, G3P1C2V
      DOUBLE PRECISION :: G1P2C2V, G2P2C2V, G3P2C2V
      DOUBLE PRECISION :: G1P3C2V, G2P3C2V, G3P3C2V
      INTEGER, PARAMETER :: ORDER1C2V=4 
      INTEGER, PARAMETER :: ORDER2C2V=3
      INTEGER, PARAMETER :: ORDER3C2V=3
      DOUBLE PRECISION :: POL1C2V, POL2C2V, POL3C2V, EHFC2V2
      DOUBLE PRECISION :: DISP, RDISPC2V, X1, X2, RANGC2V
      DOUBLE PRECISION :: GAM1
      DOUBLE PRECISION :: S1, S2, S3
      DOUBLE PRECISION :: DELTAS, Q31, RHO0, Q1ABS, Q2ABS, Q3ABS
      DOUBLE PRECISION, DIMENSION(3) :: DELTA, THETA
      DOUBLE PRECISION, DIMENSION(3) :: Q3SEAM, Q2SEAM
      DOUBLE PRECISION :: PI, ZERO
      INTEGER :: PHASE
      DOUBLE PRECISION :: LIN, RANGJT

      RDISPC2V=C(141)
      X1=C(142)
      X2=C(143)
      GAM1=C(144)

      RC2V(1)=D1
      RC2V(2)=D2
      RC2V(3)=D3
      
      A(1,1)=SQRT(1.000D+00/3.00D+00)
      A(1,2)=A(1,1)
      A(1,3)=A(1,1)
      A(2,1)=0.000D+00
      A(2,2)=SQRT(1.000D+00/2.000D+00)
      A(2,3)=-A(2,2)
      A(3,1)=SQRT(2.000D+00/3.000D+00)
      A(3,2)=-SQRT(1.000D+00/6.000D+00)
      A(3,3)=A(3,2)
      
      NUMC2V=115
      POL1C2V=0.000D+00
      POL2C2V=0.000D+00
      POL3C2V=0.000D+00
      EHFC2V2=0.000D+00
        
      DO I=1,3
        QC2V(I)=0.00D+00
          DO J=1,3
            QC2V(I)=QC2V(I)+A(I,J)*DISP(RC2V(J),RDISPC2V)
          END DO
      END DO   
          
      TAU1C2V=QC2V(1)
      TAU2C2V=QC2V(2)**2+QC2V(3)**2
      TAU3C2V=QC2V(3)**3-3.00D+00*QC2V(3)*QC2V(2)**2 
               
      DO L=0,ORDER1C2V
        DO I=0,L
         IF (I==0) THEN
           G1P1C2V=1.0000D+00
         ELSE 
           G1P1C2V=TAU1C2V**I
         END IF
         DO J=0,(L-I),2
           IF (J==0) THEN
             G2P1C2V=1.0000D+00
           ELSE
             G2P1C2V=TAU2C2V**(J/2)
           END IF
           K=L-I-J
           IF (MOD(K,3)==0) THEN
             IF (K==0) THEN
               G3P1C2V=1.0000D+00
             ELSE
               G3P1C2V=TAU3C2V**(K/3)
             END IF
             NUMC2V=NUMC2V+1
             POL1C2V=POL1C2V+C(NUMC2V)*G1P1C2V*G2P1C2V*G3P1C2V 
           END IF
         END DO
        END DO
      END DO 

      NUMC2V=126
  
      DO L=0,ORDER2C2V
        DO I=0,L
          IF (I==0) THEN
            G1P2C2V=1.0000D+00
          ELSE 
            G1P2C2V=TAU1C2V**I
          END IF
          DO J=0,(L-I),2
            IF (J==0) THEN
              G2P2C2V=1.0000D+00
            ELSE
              G2P2C2V=TAU2C2V**(J/2)
            END IF
            K=L-I-J
            IF (MOD(K,3)==0) THEN
              IF (K==0) THEN
                G3P2C2V=1.0000D+00
              ELSE
                G3P2C2V=TAU3C2V**(K/3)
              END IF
              NUMC2V=NUMC2V+1
              POL2C2V=POL2C2V+C(NUMC2V)*G1P2C2V*G2P2C2V*G3P2C2V 
            END IF
          END DO
        END DO
      END DO

      NUMC2V=133

      DO L=0,ORDER3C2V
        DO I=0,L
          IF (I==0) THEN
            G1P3C2V=1.0000D+00
          ELSE 
            G1P3C2V=TAU1C2V**I
          END IF
          DO J=0,(L-I),2
            IF (J==0) THEN
              G2P3C2V=1.0000D+00
            ELSE
              G2P3C2V=TAU2C2V**(J/2)
            END IF
            K=L-I-J
            IF (MOD(K,3)==0) THEN
              IF (K==0) THEN
                G3P3C2V=1.0000D+00
              ELSE
                G3P3C2V=TAU3C2V**(K/3)
              END IF
              NUMC2V=NUMC2V+1
              POL3C2V=POL3C2V+C(NUMC2V)*G1P3C2V*G2P3C2V*G3P3C2V 
            END IF
          END DO
        END DO
      END DO

      Q1ABS=A(1,1)*(D1+D2+D3)
      Q2ABS=A(2,2)*D2+A(2,3)*D3
      Q3ABS=A(3,1)*D1+A(3,2)*D2+A(3,3)*D3

      PI=3.1415926535897932D+00

      ZERO=0.00000000000000D+00

      THETA=(/(PI/2.00D+00),(7.00D+00*PI/6.00D+00),
     &        (11.00D+00*PI/6.00D+00)/)

      RHO0=ABS(Q31(Q1ABS))

      IF (Q31(Q1ABS)>=ZERO) THEN
        PHASE=0
      ELSE IF (Q31(Q1ABS)<ZERO) THEN
        PHASE=1
      ENDIF

      Q3SEAM(1)=RHO0*SIN(THETA(1)+DBLE(PHASE)*PI)
      Q2SEAM(1)=RHO0*COS(THETA(1)+DBLE(PHASE)*PI)
      Q3SEAM(2)=RHO0*SIN(THETA(2)+DBLE(PHASE)*PI)
      Q2SEAM(2)=RHO0*COS(THETA(2)+DBLE(PHASE)*PI)
      Q3SEAM(3)=RHO0*SIN(THETA(3)+DBLE(PHASE)*PI)
      Q2SEAM(3)=RHO0*COS(THETA(3)+DBLE(PHASE)*PI)

      DELTA(1)=SQRT((Q2ABS-Q2SEAM(1))**2+
     &              (Q3ABS-Q3SEAM(1))**2)
      DELTA(2)=SQRT((Q2ABS-Q2SEAM(2))**2+
     &              (Q3ABS-Q3SEAM(2))**2)
      DELTA(3)=SQRT((Q2ABS-Q2SEAM(3))**2+
     &              (Q3ABS-Q3SEAM(3))**2)

      S1=A(1,1)*(DELTA(1)+DELTA(2)+DELTA(3))
      S2=A(2,2)*DELTA(2)+A(2,3)*DELTA(3)
      S3=A(3,1)*DELTA(1)+A(3,2)*DELTA(2)+A(3,3)*DELTA(3)

      DELTAS=SQRT(S2**2+S3**2)

      LIN=SQRT(TAU2C2V)*RANGJT(D1,D2,D3)

      EHFC2V2=(POL1C2V-DELTAS*POL2C2V-
     &         LIN*POL3C2V)*
     &         RANGC2V(D1,D2,D3,X1,X2,GAM1)

      END SUBROUTINE THRBODYC2V2

C####################################################################################

      SUBROUTINE SCALEDDYNPOT23(R1,R2,R3,SCALEDPOT23)
      IMPLICIT NONE
      DOUBLE PRECISION :: R1, R2, R3, SCALEDPOT23 
      DOUBLE PRECISION :: TWOBODYDYN, TWOBODYHF, VDC3
      DOUBLE PRECISION :: DEE_DIATOMIC, SWIT
      DOUBLE PRECISION, DIMENSION(3) :: EHF, DYN, V2, TWOBODY 

      DEE_DIATOMIC=0.000000000D+00

C      CALL POT2BODY(R1,EHF(1),DYN(1),V2(1))
C      CALL POT2BODY(R2,EHF(2),DYN(2),V2(2))
C      CALL POT2BODY(R3,EHF(3),DYN(3),V2(3))

      CALL POTC2PIUCMRR(R1,EHF(1),DYN(1),V2(1))
      CALL POTC2PIUCMRR(R2,EHF(2),DYN(2),V2(2))
      CALL POTC2PIUCMRR(R3,EHF(3),DYN(3),V2(3))
 
      TWOBODY(1)=DYN(1)*(1-SWIT(R2,R3,R1))*(1-SWIT(R3,R1,R2))

      TWOBODY(2)=DYN(2)*(1-SWIT(R1,R2,R3))*(1-SWIT(R3,R1,R2))

      TWOBODY(3)=DYN(3)*(1-SWIT(R1,R2,R3))*(1-SWIT(R2,R3,R1))

      TWOBODYDYN=SUM(TWOBODY)

      TWOBODYHF=SUM(EHF)

      CALL POT3DC(R1,R2,R3,VDC3)
  
      SCALEDPOT23=TWOBODYHF+TWOBODYDYN+3*DEE_DIATOMIC+VDC3
  
      END SUBROUTINE SCALEDDYNPOT23

C####################################################################################

      SUBROUTINE POT3DC(R1,R2,R3,VDC3)
      IMPLICIT NONE
      DOUBLE PRECISION :: R1, R2, R3, DC3, VDC3
      DOUBLE PRECISION, DIMENSION(3) :: RE, rjac, cosgamma, RB, RC
      DOUBLE PRECISION, DIMENSION(10,3) :: DM, RM
      DOUBLE PRECISION, DIMENSION(3,10,3) :: a, b
      DOUBLE PRECISION, DIMENSION(3,3) :: PLegend
      DOUBLE PRECISION, DIMENSION(10) :: C_DISP_DIATOM
      DOUBLE PRECISION, DIMENSION(10,3,3) :: C_RE
      DOUBLE PRECISION, DIMENSION(10,3) :: C_RE_GAMMA
      INTEGER :: lmax
      INTEGER :: n, l, channel, poli
      DOUBLE PRECISION :: r
      DOUBLE PRECISION, DIMENSION(10,3) :: POLI_A, POLI_B
      DOUBLE PRECISION :: R0_ATOM_DIATOM
C LEE ROY RADIUS (R0) OF C-Mg
      DOUBLE PRECISION, PARAMETER :: R0 = 7.89107343775231D+00 
      DOUBLE PRECISION :: JABOBI_RADIUS_EXP, DAMP_PART, SWIT_PART 
      DOUBLE PRECISION :: DAMP, SWIT

      R0_ATOM_DIATOM = 10.983323302739D+00

      CALL TRIANG2JACOBI2(R1,R2,R3,RE,rjac,cosgamma,RB,RC)

      RM(6,1) = 4.5000D+00;       DM(6,1) = 27.5279D+00;
      a(1,6,1) = 0.99499999D+00;  a(2,6,1) = 0.25778783D+00;
      a(3,6,1) = 0.00139584D+00;  b(1,6,1) = a(1,6,1);
      b(2,6,1) = 0.27435311D+00;  b(3,6,1) = 0.02546246D+00;

      RM(6,2) = 4.5000D+00;       DM(6,2) = 28.5233D+00;
      a(1,6,2) = 0.76362610D+00;  a(2,6,2) = 0.21435356D+00;
      a(3,6,2) = -0.00008571D+00; b(1,6,2) = a(1,6,2);
      b(2,6,2) = 0.24930236D+00;  b(3,6,2) = 0.00676850D+00;
 
      RM(8,1) = 4.4711D+00;       DM(8,1) = 1037.6370D+00
      a(1,8,1) = 0.93135973D+00;  a(2,8,1) = 0.23731725D+00;
      a(3,8,1) = -0.00013049D+00; b(1,8,1) = a(1,8,1);
      b(2,8,1) = 0.24399061D+00;  b(3,8,1)= 0.02486858;
 
      RM(8,2) = 4.4823D+00;       DM(8,2) = 2941.9456D+00;
      a(1,8,2) = 0.79819039D+00;  a(2,8,2) = 0.22589941D+00;
      a(3,8,2) = -0.00082011D+00; b(1,8,2) = a(1,8,2);
      b(2,8,2) = 0.29644898D+00;  b(3,8,2) = 0.02109793D+00;
 
      RM(8,3) = 4.4873D+00;       DM(8,3) = 452.6876D+00;
      a(1,8,3) = 1.23539671D+00;  a(2,8,3) = 0.57644149D+00;
      a(3,8,3) = 0.07460865D+00;  b(1,8,3) = a(1,8,3);
      b(2,8,3) = 0.58644002D+00;  b(3,8,3) = 0.07326621D+00;
  
      RM(10,1) = 4.4525D+00;      DM(10,1) = 46263.4796D+00;
      a(1,10,1) = 0.89618614D+00; a(2,10,1) = 0.22325495D+00;
      a(3,10,1)= -0.00182731D+00; b(1,10,1) = a(1,10,1);
      b(2,10,1)= 0.23140031D+00;  b(3,10,1) = 0.02681374D+00;

      DO channel=1,3  !! LOOP FOR EACH CHANNEL

      PLegend(1,channel) = 1.0D+00      

      PLegend(2,channel) = 0.5D+00*(3.0D+00*cosgamma(channel)**2
c23456
     + - 1.0D+00)                                     

      PLegend(3,channel) = (35.0D+00*cosgamma(channel)**4 
c23456
     + - 30.0D+00*cosgamma(channel)**2 + 3.0D+00)/8.0D+00   

      END DO

      C_DISP_DIATOM(6) = 40.9D+00

      C_DISP_DIATOM(8) = C_DISP_DIATOM(6)*R0**(1.54)  

      C_DISP_DIATOM(10) = C_DISP_DIATOM(6)*1.31D+00*R0**(2.0*1.54)

      DO channel=1,3         
      DO n=6,10,2            

      IF (n == 6) THEN  

      lmax = 2           
 
      ELSE IF (n == 8) THEN

      lmax = 3             

      ELSE IF (n == 10) THEN 

      lmax = 1             

      END IF 

      DO l=1,lmax

      r = RE(channel) - RM(n,l)

      POLI_A(n,l) = 0.0D+00

      POLI_B(n,l) = 0.0D+00

      DO poli = 1,3
     
      POLI_A(n,l) = POLI_A(n,l) + a(poli,n,l)*r**(poli)

      POLI_B(n,l) = POLI_B(n,l) - b(poli,n,l)*r**(poli)
           
      END DO

      IF (l == 1) THEN

      C_RE(n,l,channel) = C_DISP_DIATOM(n) + C_DISP_DIATOM(n)
c23456
     1 + DM(n,l)*(1.0D+00 + POLI_A(n,l))*EXP(POLI_B(n,l))		 

      ELSE

      C_RE(n,l,channel) =
c23456
     1  DM(n,l)*(1.0D+00 + POLI_A(n,l))*EXP(POLI_B(n,l)) 
                     
      END IF

      END DO
      END DO                  
      END DO                  

      DO channel=1,3         
      DO n=6,10,2            

      IF (n == 6) THEN  

      lmax = 2            
 
      ELSE IF (n == 8) THEN

      lmax = 3             

      ELSE IF (n == 10) THEN 

      lmax = 1             

      END IF 

      C_RE_GAMMA(n,channel) = 0.0D+00

      DO l=1,lmax

      C_RE_GAMMA(n,channel) = C_RE_GAMMA(n,channel) 
c23456
     1 + C_RE(n,l,channel)*PLegend(l,channel)   
     
      END DO

      END DO                
      END DO                 

      DC3 = 0.0D+00

      JABOBI_RADIUS_EXP = 0.0D+00    

      DAMP_PART = 0.0D+00            

      SWIT_PART = 0.0D+00            

      DO channel = 1,3

      DO n=6,10,2
  
      JABOBI_RADIUS_EXP = rjac(channel)**(-n)

      DAMP_PART = DAMP(rjac(channel),R0_ATOM_DIATOM,n)

      SWIT_PART = SWIT(RE(channel),RB(channel),RC(channel))

      DC3 = DC3 
c23456
     1 + SWIT_PART*C_RE_GAMMA(n,channel)*DAMP_PART*JABOBI_RADIUS_EXP   

      END DO

      END DO

      VDC3 = - DC3

      END SUBROUTINE POT3DC


C####################################################################################

      SUBROUTINE TRIANG2JACOBI2(R1,R2,R3,RE,RJAC,COSGAMMA,RB,RC)
      IMPLICIT NONE
      DOUBLE PRECISION :: R1, R2, R3
      DOUBLE PRECISION, PARAMETER  :: ZERO = 1.0D-12
      DOUBLE PRECISION, DIMENSION(3) :: RE
      DOUBLE PRECISION, DIMENSION(3) :: RJAC
      DOUBLE PRECISION, DIMENSION(3) :: RB, RC
      DOUBLE PRECISION, DIMENSION(3) :: SQUAREDRJAC
      DOUBLE PRECISION, DIMENSION(3) :: COSGAMMA
      DOUBLE PRECISION :: PI
      INTEGER :: I

      PI=3.1415926535897932384626433832795D+00

      DO I=1,3                !! CYCLIC PERMUTATIONS OF 1,2,3

        IF (I==1) THEN          !! CHANNEL 1 
 
          RE(I)=R1
          RB(I)=R2
          RC(I)=R3

        ELSE IF (I==2) THEN     !! CHANNEL 2

          RE(I)=R2
          RB(I)=R3
          RC(I)=R1

        ELSE IF (I==3) THEN     !! CHANNEL 3

          RE(I)=R3
          RB(I)=R1
          RC(I)=R2

        END IF

        RJAC(I)=ZERO

        SQUAREDRJAC(I)=0.5D+00*(RB(I)**2+RC(I)**2-0.5D+00*RE(I)**2)

          IF (SQUAREDRJAC(I) > 0) THEN 

            RJAC(I)=sqrt(SQUAREDRJAC(I)) 

          END IF
    
        COSGAMMA(I)=(0.5D+00*(RC(I)**2-RB(I)**2))/(RJAC(I)*RE(I))

          IF (COSGAMMA(I)>1.0D+00) THEN 

            COSGAMMA(I)=1.0D+00 

          ELSE IF (COSGAMMA(I)<-1.0D+00) THEN 
   
            COSGAMMA(I)=-1.0D+00

          END IF 

      END DO

      END SUBROUTINE TRIANG2JACOBI2

C####################################################################################
C EXPERIMENTALLY-DETERMINED C2 PI U TRIPLET CURVE CMRR (2018)
C####################################################################################

      SUBROUTINE POTC2PIUCMRR(R,VEHF,VDC,VCC)
      IMPLICIT NONE
      INTEGER, PARAMETER :: P=16
      DOUBLE PRECISION, DIMENSION(P) :: C
      DOUBLE PRECISION, DIMENSION(4) :: A, B, QUI, CN, DC
      INTEGER, DIMENSION(4) :: ID = (/5,6,8,10/)
      DOUBLE PRECISION, DIMENSION(8) :: COEF, POLI 
      INTEGER :: I
      DOUBLE PRECISION :: R, RE, X, R0, RHO, RHO_LINHA  
      DOUBLE PRECISION :: VDC, D, GAMMA, VEHF
      DOUBLE PRECISION :: G0, G1, G2, VCC, DEE
      DOUBLE PRECISION :: AN, BN

       C(  1)= 0.24805692D+01
       C(  2)= 0.46782041D+00
       C(  3)= 0.16367516D+01
       C(  4)= 0.11584156D+01
       C(  5)= 0.51780588D+00
       C(  6)= 0.00000000D+00
       C(  7)= 0.14540400D+02
       C(  8)= 0.40900000D+02
       C(  9)= 0.21856288D+01
       C( 10)= 0.16318276D+01
       C( 11)= 0.21122385D+01
       C( 12)= 0.95811379D+00
       C( 13)= 0.78546378D+00
       C( 14)= 0.28131851D+00
       C( 15)= 0.14908734D+00
       C( 16)=-0.16704804D-01

      RE=C(1)
      D=C(2)
      G0=C(3)
      G1=C(4)
      G2=C(5)
      DEE=C(6)

      X=R-RE

      GAMMA=G0*(1.0D+00+G1*TANH(G2*X))

      COEF(1)=C(9)
      COEF(2)=C(10)
      COEF(3)=C(11)
      COEF(4)=C(12)
      COEF(5)=C(13)
      COEF(6)=C(14)
      COEF(7)=C(15)
      COEF(8)=C(16)

      DO I=1,8
        POLI(I)=COEF(I)*X**I	
      END DO

      VEHF=-(D/R)*(1.0D+00+SUM(POLI))*EXP(-GAMMA*X) 

      R0=7.89107343775231D+00
      RHO=5.5D+00+1.25D+00*R0 
      RHO_LINHA=R/RHO

      CN(1)=C(7)
      CN(2)=C(8)
      CN(3)=CN(2)*R0**(1.54)
      CN(4)=CN(2)*1.31D+00*R0**(2.0*1.54)

      DO I=1,4
        A(I)=AN(ID(I))
      END DO

      DO I=1,4
        B(I)=BN(ID(I))
      END DO

      DO I=1,4
        QUI(I)=(1.0D+00-EXP(-A(I)*RHO_LINHA-B(I)*RHO_LINHA**2))**(ID(I))
      END DO 

      DO I=1,4
        DC(I)=-CN(I)*QUI(I)/R**(ID(I)) 
      END DO 
                            
      VDC=SUM(DC)                           
                                          
      VCC=VEHF+VDC+DEE

      END SUBROUTINE POTC2PIUCMRR

C####################################################################################

      SUBROUTINE POT2BODY(R,VEHF,VDC,VCC)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(4) :: A, B, QUI, CN, DC
      INTEGER, DIMENSION(4) :: ID = (/5,6,8,10/)
      DOUBLE PRECISION, DIMENSION(10) :: COEF, POLI 
      INTEGER :: I
      DOUBLE PRECISION :: R, RE, X, R0, RHO, RHO_LINHA  
      DOUBLE PRECISION :: VDC, D, GAMMA, VEHF, G0, G1, G2, VCC, DEE
      DOUBLE PRECISION :: AN, BN
  
      RE=2.47932D+00
      X=R-RE
      D=0.467489D+00
      G0=0.286063D+00
      G1=15.5532D+00
      G2=0.0736335D+00	
      DEE=0.000000000D+00
      GAMMA=G0*(1.0D+00+G1*TANH(G2*X))

      COEF(1)=0.850571D+00
      COEF(2)=-1.07115D+00
      COEF(3)=0.969122D+00
      COEF(4)=-0.612554D+00
      COEF(5)=0.310833D+00
      COEF(6)=-0.0944244D+00
      COEF(7)=-0.00736976D+00
      COEF(8)=0.0137852D+00
      COEF(9)=-0.00333289D+00
      COEF(10)=0.000257898D+00

      DO I=1,10
        POLI(I)=COEF(I)*X**I	
      END DO

      VEHF=-(D/R)*(1.0D+00+SUM(POLI))*EXP(-GAMMA*X) 

      R0=7.89107343775231D+00
      RHO=5.5D+00+1.25D+00*R0 
      RHO_LINHA=R/RHO

      CN(1)=14.5404D+00
      CN(2)=40.9D+00
      CN(3)=CN(2)*R0**(1.54)
      CN(4)=CN(2)*1.31D+00*R0**(2.0*1.54)

      DO I=1,4
        A(I)=AN(ID(I))
      END DO

      DO I=1,4
        B(I)=BN(ID(I))
      END DO

      DO I=1,4
        QUI(I)=(1.0D+00-EXP(-A(I)*RHO_LINHA-B(I)*RHO_LINHA**2))**(ID(I))
      END DO 

      DO I=1,4
        DC(I)=-CN(I)*QUI(I)/R**(ID(I)) 
      END DO 
                            
      VDC=SUM(DC)                           
                                          
      VCC=VEHF+VDC+DEE

      END SUBROUTINE POT2BODY

C###############################################################################

      DOUBLE PRECISION FUNCTION ES(BETA,DEMIN,DEREF,N)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: BETA, DEMIN, DEREF
      INTEGER, INTENT(IN) :: N
      DOUBLE PRECISION :: RATIO
      IF (DEMIN>=DEREF) THEN
        ES=1.000D+00
      ELSE IF (DEMIN<DEREF) THEN
        RATIO=((DEREF/DEMIN)-1.000D+00)
        ES=EXP(-BETA*RATIO**N)
      END IF
      END FUNCTION ES

C####################################################################################

      DOUBLE PRECISION FUNCTION Q31(Q1)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: Q1
      DOUBLE PRECISION :: A, B, C, D, E, F, G, H 
      A=0.00454872897741226D+00
      B=0.02509103210095290D+00
      C=1.05938652067574000D+00
      D=4.78291330706927000D+00
      E=-0.9261734635796440D+00
      F=-0.2164776500990330D+00
      G=0.17221936619583700D+00
      H=0.68281911615943100D+00
      Q31=A-B*TANH(C*(Q1-D)+E*(Q1-D)**2+
     &     F*(Q1-D)**3+G*(Q1-D)**4+H*(Q1-D)**5)
      END FUNCTION Q31

C####################################################################################

      DOUBLE PRECISION FUNCTION RANGJT(R1,R2,R3)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: R1,R2,R3
      DOUBLE PRECISION, PARAMETER :: X0=2.885169410D+00
      DOUBLE PRECISION, PARAMETER :: GAM1=10000.00D+00
      RANGJT=(1.00D+00-EXP(-GAM1*((R1-X0)**2+(R2-X0)**2+(R3-X0)**2)))
      END FUNCTION RANGJT

C####################################################################################

      DOUBLE PRECISION FUNCTION DISP(R,RDISP)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: R,RDISP
      DISP=(R-RDISP)
      END FUNCTION DISP 

C####################################################################################

      DOUBLE PRECISION FUNCTION RANG(R1,R2,R3,X0,GAM1)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: R1,R2,R3,X0,GAM1
      RANG=(1.00D+00-TANH(GAM1*(R1-X0)))*
     &     (1.00D+00-TANH(GAM1*(R2-X0)))*
     &     (1.00D+00-TANH(GAM1*(R3-X0)))
      END FUNCTION RANG

C####################################################################################

      DOUBLE PRECISION FUNCTION RANGC2V(R1,R2,R3,X1,X2,GAM1)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: R1,R2,R3,X1,X2,GAM1
      DOUBLE PRECISION, DIMENSION(3) :: RANG
      RANG(1)=(1.00D+00-TANH(GAM1*(R1-X1)))*
     &        (1.00D+00-TANH(GAM1*(R2-X2)))*
     &        (1.00D+00-TANH(GAM1*(R3-X2)))
      RANG(2)=(1.00D+00-TANH(GAM1*(R3-X1)))*
     &        (1.00D+00-TANH(GAM1*(R1-X2)))*
     &        (1.00D+00-TANH(GAM1*(R2-X2)))
      RANG(3)=(1.00D+00-TANH(GAM1*(R2-X1)))*
     &        (1.00D+00-TANH(GAM1*(R3-X2)))*
     &        (1.00D+00-TANH(GAM1*(R1-X2)))
      RANGC2V=RANG(1)+RANG(2)+RANG(3)
      END FUNCTION RANGC2V

C####################################################################################

      DOUBLE PRECISION FUNCTION DAMP(RJAC,R0,N)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: RJAC, R0
      INTEGER, INTENT(IN):: N
      DOUBLE PRECISION :: RHO, RHO_LINHA, AN, BN
      DOUBLE PRECISION, DIMENSION(10) :: A, B
      RHO = 5.5D+00 + 1.25D+00*R0 
      RHO_LINHA=RJAC/RHO
      A(N)=AN(N)
      B(N)=BN(N)
      DAMP=(1.0D+00-EXP(-A(N)*RHO_LINHA-B(N)*RHO_LINHA**2))**N
      END FUNCTION DAMP

C####################################################################################

      DOUBLE PRECISION FUNCTION AN(N)
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: ALFA0=16.36606000D+00
      DOUBLE PRECISION, PARAMETER :: ALFA1=0.70172000D+00
      INTEGER, INTENT(IN):: N
      AN=ALFA0/(DBLE(N))**ALFA1
      END FUNCTION AN 

C####################################################################################

      DOUBLE PRECISION FUNCTION BN(N)
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: BETA0=17.19338D+00
      DOUBLE PRECISION, PARAMETER :: BETA1=0.09574D+00
      INTEGER, INTENT(IN):: N
      BN=BETA0*DEXP(-BETA1*DBLE(N))
      END FUNCTION BN

C####################################################################################

      DOUBLE PRECISION FUNCTION SWIT(RA,RB,RC)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: RA, RB, RC
      DOUBLE PRECISION, PARAMETER :: ETA=6.0D+00, CSI=1.0D+00
      SWIT=0.5D+00*(1.0D+00-TANH(CSI*(ETA*RA-RB-RC)))
      END FUNCTION SWIT

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
      CALL FCN(NTA,X,VLP)
      RETURN
      END

C###############################################################################
C THREE-BODY TERM FOR GROUND-STATE TRIPLET C3
C C.M.R. ROCHA AND A.J.C. VARANDAS 2019
C###############################################################################

      SUBROUTINE CHIPR_TRIAT_C3T(R,POT)
      IMPLICIT NONE
      CHARACTER(LEN=3), PARAMETER :: MOLTYP="A3"
      INTEGER, PARAMETER :: DEG=1
      INTEGER, PARAMETER :: NC=705
      INTEGER, PARAMETER :: NX=3
      INTEGER :: I,J,K,L,M,S,O,ID
      INTEGER :: POLORDER
      INTEGER, DIMENSION(DEG) :: BSORDER,NCBAS
      DOUBLE PRECISION, DIMENSION(NC) :: C
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: POL,BS
      DOUBLE PRECISION, DIMENSION(3) :: R
      DOUBLE PRECISION :: POT,REPDAMP
      INTEGER :: NCPOL,NCTOTAL,SUMC
      DOUBLE PRECISION, DIMENSION(NX) :: Y
      INTEGER :: TOTNUM,MNUM
      INTEGER, DIMENSION(3,6) :: P

      C(  1)=  0.386462686309568504D+01
      C(  2)= -0.164070340740438724D+03
      C(  3)=  0.136650571156685929D+01
      C(  4)=  0.594066162558292763D+04
      C(  5)= -0.546477967176414313D+04
      C(  6)=  0.183929701733514275D+04
      C(  7)= -0.411484961577390204D+05
      C(  8)=  0.465346682154638984D+05
      C(  9)= -0.201479044882052622D+05
      C( 10)= -0.106266971129748049D+04
      C( 11)=  0.105129333333067843D+06
      C( 12)=  0.287366233979977333D+05
      C( 13)= -0.175843553177008784D+06
      C( 14)=  0.848014350072227535D+05
      C( 15)= -0.652616109502490726D+05
      C( 16)=  0.687610666957922513D+05
      C( 17)= -0.824770458315127325D+05
      C( 18)= -0.228966336982557521D+06
      C( 19)=  0.496001565208565735D+06
      C( 20)= -0.147135788096068572D+06
      C( 21)= -0.196825500374073818D+04
      C( 22)=  0.412547120801116456D+06
      C( 23)= -0.423924571153475728D+06
      C( 24)= -0.289628220453600661D+05
      C( 25)= -0.164119875379085424D+06
      C( 26)=  0.113606806152462712D+06
      C( 27)= -0.503364434600585490D+06
      C( 28)= -0.931645031300408300D+04
      C( 29)=  0.306355800293331966D+06
      C( 30)= -0.781411000648871181D+06
      C( 31)=  0.667840341512434970D+05
      C( 32)=  0.492959098238230392D+06
      C( 33)= -0.140439067081986141D+06
      C( 34)=  0.535783144921971834D+06
      C( 35)=  0.316492429735922837D+06
      C( 36)= -0.263191665884998201D+03
      C( 37)=  0.468330498258402047D+06
      C( 38)= -0.322067137468170622D+05
      C( 39)= -0.109943240767227928D+07
      C( 40)= -0.201933304945216252D+06
      C( 41)=  0.840807848529611132D+06
      C( 42)=  0.914569881201875745D+06
      C( 43)= -0.394969927406983974D+06
      C( 44)=  0.468655104409140244D+06
      C( 45)=  0.565470036196962581D+06
      C( 46)= -0.803256314635783201D+06
      C( 47)=  0.475311931162403489D+06
      C( 48)=  0.877056189519908075D+05
      C( 49)= -0.402104521226408659D+06
      C( 50)= -0.504692505928197352D+06
      C( 51)= -0.878469432875571511D+05
      C( 52)=  0.165611433509597252D+07
      C( 53)= -0.121656931896444899D+07
      C( 54)=  0.884262131974278484D+06
      C( 55)= -0.188322133399405191D+07
      C( 56)= -0.720437908603700227D+06
      C( 57)=  0.251332616845693876D+06
      C( 58)= -0.535172339965660012D+04
      C( 59)= -0.878706024441424915D+03
      C( 60)= -0.249303517025618094D+05
      C( 61)= -0.399360973575325670D+04
      C( 62)= -0.107991348420195375D+07
      C( 63)=  0.480309770961557864D+06
      C( 64)=  0.144359064771516528D+07
      C( 65)=  0.459367991809543967D+07
      C( 66)=  0.110411095646536894D+06
      C( 67)= -0.179148465572928800D+07
      C( 68)= -0.175824683433538041D-01
      C( 69)= -0.903363882408217323D+05
      C( 70)=  0.190726092384676449D+06
      C( 71)= -0.515995806228434049D+04
      C( 72)= -0.201351868607002450D+07
      C( 73)= -0.765268079156843829D+06
      C( 74)= -0.674263074700887315D+07
      C( 75)=  0.381434763608807058D+06
      C( 76)=  0.275952510072562518D+07
      C( 77)= -0.290499452501060069D+07
      C( 78)=  0.146439056385161378D+07
      C( 79)= -0.386353287657811365D+06
      C( 80)= -0.165570288952999702D+07
      C( 81)= -0.821392491639254558D+02
      C( 82)= -0.545179297968710400D+07
      C( 83)=  0.128669768245610595D+07
      C( 84)=  0.392840055047508758D+04
      C( 85)= -0.171244951927851094D+07
      C( 86)= -0.528238854746572815D+04
      C( 87)=  0.164611428449600982D+07
      C( 88)=  0.132126536321979156D+07
      C( 89)=  0.243783386917905556D+07
      C( 90)=  0.189741459249548009D+07
      C( 91)=  0.117812970625130197D+04
      C( 92)=  0.404801736171090696D+07
      C( 93)=  0.793365421890168800D+06
      C( 94)= -0.357791020895039692D+05
      C( 95)= -0.139543231313403671D+01
      C( 96)= -0.419813425719083977D+01
      C( 97)=  0.323444998093351629D+07
      C( 98)=  0.856278029792423826D+06
      C( 99)=  0.142255917834987006D+00
      C(100)= -0.217572292258743383D+07
      C(101)=  0.829937478479702986D+05
      C(102)=  0.209572409792559594D+07
      C(103)=  0.616175026903293110D+04
      C(104)= -0.562007304357786570D+07
      C(105)=  0.283963316780139168D+06
      C(106)=  0.303532069061397607D+02
      C(107)= -0.132924837776597613D+07
      C(108)=  0.537860011326507447D+05
      C(109)= -0.878101016674569459D+06
      C(110)= -0.254903479186414827D+01
      C(111)= -0.159395983478838066D+07
      C(112)=  0.110355129660328142D+08
      C(113)=  0.371800440138146933D+07
      C(114)= -0.338211904683074774D+07
      C(115)= -0.232573766033829749D+07
      C(116)=  0.564600981012214068D+07
      C(117)=  0.216413105607046699D+07
      C(118)=  0.118775121517732204D+07
      C(119)=  0.363776010460503120D+07
      C(120)= -0.340930368700196454D+07
      C(121)= -0.139333710602482641D+07
      C(122)=  0.638155253402769798D+06
      C(123)=  0.413298902613121551D+07
      C(124)= -0.855835380753751285D+07
      C(125)=  0.930903018531226189D+04
      C(126)= -0.427422775554060005D+07
      C(127)=  0.744287432505308464D+07
      C(128)=  0.136816544491493795D+08
      C(129)=  0.127807340726956185D+08
      C(130)= -0.110783887228928152D+08
      C(131)= -0.706831071566507779D+07
      C(132)= -0.108224045292860275D+06
      C(133)= -0.833531024908608058D+06
      C(134)= -0.107882665627681147D+08
      C(135)=  0.116237540576072074D+08
      C(136)= -0.799615053139260598D+07
      C(137)=  0.683923129990253225D+07
      C(138)=  0.152356936959359991D+05
      C(139)= -0.349714145926694013D+07
      C(140)= -0.833698789967783168D+07
      C(141)= -0.749162661306171212D+07
      C(142)=  0.534174839569014031D+07
      C(143)= -0.697272494506542478D+07
      C(144)=  0.975704784493953735D+07
      C(145)= -0.801353356090133358D+07
      C(146)=  0.510442221873073088D+01
      C(147)=  0.892683207053336147D+01
      C(148)= -0.143951275691942200D+08
      C(149)= -0.172080162112242659D+07
      C(150)=  0.146333166476405468D+08
      C(151)=  0.167383831971946899D+08
      C(152)= -0.372886157474681735D+07
      C(153)=  0.226993333195227906D+08
      C(154)= -0.143255642105382867D+08
      C(155)= -0.129893833297079839D+08
      C(156)= -0.308642905840877332D+08
      C(157)= -0.181628906245966516D+08
      C(158)=  0.433491387812812105D+08
      C(159)=  0.676631203361410275D+07
      C(160)= -0.259130313897013020D+00
      C(161)= -0.961687765642873943D+07
      C(162)=  0.517456277323312008D+01
      C(163)= -0.954764685009295680D+07
      C(164)=  0.528658799552335190D+00
      C(165)=  0.171143291579988960D+06
      C(166)=  0.292257712662187172D+07
      C(167)=  0.636572789685283136D+07
      C(168)= -0.224761384416173138D+08
      C(169)= -0.127610350221335851D+08
      C(170)= -0.215765049424951933D+08
      C(171)= -0.185157996439657477D+07
      C(172)=  0.833408767796584126D+07
      C(173)=  0.464798656919078156D+07
      C(174)=  0.180800452303782925D+08
      C(175)= -0.917335487233920954D+07
      C(176)= -0.478531052602523845D+07
      C(177)=  0.144587474624043945D+08
      C(178)=  0.171547678450389877D+08
      C(179)= -0.169633099706774810D+07
      C(180)= -0.275554110375516862D+08
      C(181)=  0.311808379863139689D+08
      C(182)=  0.721591452277420554D+07
      C(183)= -0.126093994285622947D+08
      C(184)=  0.350830001008094405D+05
      C(185)= -0.438553308708212301D+08
      C(186)= -0.333407868351951838D+08
      C(187)=  0.543348450918029845D+08
      C(188)= -0.688191786799972970D+07
      C(189)=  0.263269719674374834D+08
      C(190)= -0.786502292794325228D+02
      C(191)= -0.122184650683076866D+07
      C(192)= -0.778847871111222630D+05
      C(193)=  0.342459288507881714D+07
      C(194)= -0.110116452321333084D+08
      C(195)=  0.201261170329717621D+08
      C(196)=  0.181327144860393517D+08
      C(197)= -0.286799910762527725D+07
      C(198)= -0.108213523826157134D+08
      C(199)=  0.713060724239228172D+00
      C(200)=  0.115025923937183776D+06
      C(201)=  0.281230631887282841D+08
      C(202)= -0.167934266525046514D+05
      C(203)=  0.599405529788874105D+02
      C(204)= -0.269990479973975243D+04
      C(205)=  0.112379469777573599D+07
      C(206)=  0.214728820388520658D+08
      C(207)= -0.165441102660050448D+08
      C(208)=  0.322836733364111818D+08
      C(209)=  0.237017195763823256D+06
      C(210)= -0.429112893354753107D+08
      C(211)=  0.570003063182277977D+08
      C(212)=  0.467064999735803902D+08
      C(213)=  0.339297717310808375D+08
      C(214)= -0.122271582230637092D+08
      C(215)=  0.973739051453130879D+07
      C(216)= -0.593646046950091347D+08
      C(217)= -0.663562296538148075D+08
      C(218)= -0.405976655894962279D+07
      C(219)=  0.786264280664065182D+08
      C(220)= -0.109204908988351468D+08
      C(221)=  0.253097778337490149D+08
      C(222)= -0.247334369080773667D+08
      C(223)= -0.366138953869888978D+07
      C(224)= -0.437955064935936034D+08
      C(225)=  0.361342107202808410D+08
      C(226)=  0.165027251049523689D+08
      C(227)= -0.157805452567070816D+08
      C(228)= -0.595364718303564861D+08
      C(229)= -0.111782915638171147D+03
      C(230)=  0.734902567933269404D+07
      C(231)=  0.552554321561603770D+08
      C(232)=  0.601661529328814223D+08
      C(233)= -0.281908449999750592D+08
      C(234)=  0.545908126279407367D+07
      C(235)= -0.162993096643027868D+08
      C(236)= -0.249384616817251630D+08
      C(237)= -0.178199523329258809D+06
      C(238)=  0.688283439799518287D+08
      C(239)=  0.563785141822805861D+06
      C(240)= -0.419209807013922255D+05
      C(241)= -0.576943260261163712D+08
      C(242)= -0.825466742070768028D+08
      C(243)= -0.663709367215185985D+08
      C(244)= -0.944825212311260253D+08
      C(245)=  0.734274224333852082D+08
      C(246)=  0.356874410632450879D+08
      C(247)=  0.790361466078581810D+08
      C(248)=  0.771682446062758416D+08
      C(249)= -0.585620868054117262D+08
      C(250)= -0.565546977497173026D+08
      C(251)= -0.695553965075203776D+08
      C(252)= -0.764317227010574639D+08
      C(253)= -0.461840975582094193D+08
      C(254)= -0.123288103432185203D+07
      C(255)=  0.177094438884616606D+08
      C(256)= -0.956874010479725075D+04
      C(257)=  0.144862389942891467D+08
      C(258)= -0.586207142941665575D+08
      C(259)= -0.499223183290395811D+08
      C(260)=  0.552429289963045704D+02
      C(261)= -0.163181768882941850D+07
      C(262)=  0.194263156149757393D+08
      C(263)=  0.246903197050104663D+07
      C(264)=  0.234735475801909110D+07
      C(265)= -0.324525889978753626D+08
      C(266)=  0.535706769212359637D+08
      C(267)= -0.325003472399994619D+08
      C(268)=  0.513113711422717348D+08
      C(269)=  0.194237284470249265D+08
      C(270)= -0.251088151150525175D+08
      C(271)=  0.529738283852853552D+08
      C(272)= -0.628214268531709604D+05
      C(273)= -0.515821083421782851D+08
      C(274)= -0.452129332280869484D+08
      C(275)=  0.165402545939466059D+09
      C(276)= -0.414580819556501582D+08
      C(277)= -0.850372207740547210D+08
      C(278)= -0.147049070983322952D+08
      C(279)=  0.535216522197716311D+08
      C(280)= -0.829653346956873089D+08
      C(281)= -0.126580667436876878D+09
      C(282)=  0.345459316658869758D+08
      C(283)= -0.375560440986494571D+08
      C(284)=  0.117897553507444769D+09
      C(285)=  0.144392325537598193D+09
      C(286)=  0.301733682162267491D+08
      C(287)=  0.668367409418370053D+08
      C(288)= -0.479943196287087351D+08
      C(289)= -0.133740970734640596D+01
      C(290)= -0.112249304177569777D+09
      C(291)=  0.615258876777558252D+08
      C(292)=  0.418447685374310254D+03
      C(293)=  0.865704305552286655D+08
      C(294)= -0.286779039760492481D+08
      C(295)=  0.159377750260621272D+08
      C(296)= -0.213998818570050634D+08
      C(297)=  0.116063842862110212D+08
      C(298)=  0.897899218775052077D+03
      C(299)=  0.576030931849952787D+08
      C(300)= -0.591207759814174566D+07
      C(301)= -0.730268818813881427D+08
      C(302)=  0.692903663049481511D+08
      C(303)=  0.499917373817384839D+08
      C(304)=  0.718548493092974126D+08
      C(305)= -0.676727227168784440D+08
      C(306)= -0.150405449068333069D+07
      C(307)=  0.187287789251625501D+08
      C(308)= -0.545304944439832941D+08
      C(309)=  0.916743293225164115D+08
      C(310)= -0.402350106458656569D+03
      C(311)= -0.757941286488061845D+08
      C(312)=  0.135078049848364711D+09
      C(313)=  0.107238095605589882D+09
      C(314)= -0.139065614310157508D+09
      C(315)=  0.174091788284730256D+09
      C(316)= -0.100348152636076454D+08
      C(317)= -0.186046620337421507D+09
      C(318)= -0.975444572069500983D+08
      C(319)=  0.190121080342686772D+08
      C(320)= -0.928488832087469101D+08
      C(321)= -0.162894296617259502D+09
      C(322)=  0.131522307986930587D+01
      C(323)=  0.635902796971568465D+08
      C(324)= -0.484489801807977334D+08
      C(325)=  0.101969308458226845D+09
      C(326)=  0.726833289745701998D+08
      C(327)=  0.215802548604931831D+09
      C(328)= -0.358126988735353723D+08
      C(329)= -0.113436810387704596D+09
      C(330)= -0.576120844455416314D+06
      C(331)=  0.119503034072731286D+09
      C(332)= -0.691297121845252812D+08
      C(333)=  0.152061507246922106D+09
      C(334)=  0.218460145214267015D+09
      C(335)=  0.140941410644769937D+09
      C(336)=  0.449410679032283872D+08
      C(337)= -0.147552630084336251D+09
      C(338)=  0.455198686201960470D+03
      C(339)= -0.157148872470094919D+09
      C(340)=  0.822791517197841927D+03
      C(341)=  0.248890338599490821D+09
      C(342)=  0.931804492411251068D+08
      C(343)=  0.515307969192987382D+08
      C(344)= -0.146571206858828038D+09
      C(345)=  0.552813746598252465D+03
      C(346)= -0.150947042579000397D+07
      C(347)= -0.973106879798005223D+08
      C(348)=  0.309894582702212706D+08
      C(349)=  0.175067248639494985D+09
      C(350)=  0.531544228531821519D+08
      C(351)= -0.101508070493108481D+09
      C(352)= -0.138809936987558365D+09
      C(353)= -0.756371759056946486D+08
      C(354)= -0.239159017919550854D+06
      C(355)=  0.758572091009290218D+08
      C(356)= -0.906302728813334834D+06
      C(357)=  0.209139135526038399D+04
      C(358)= -0.194898004403473675D+09
      C(359)=  0.104779754494096383D+09
      C(360)=  0.155613416403043687D+09
      C(361)= -0.254562881634380817D+09
      C(362)= -0.545874705055054277D+08
      C(363)= -0.525810952752391156D+07
      C(364)= -0.513776806989356577D+08
      C(365)= -0.103854115627336368D+09
      C(366)=  0.173844118808950037D+08
      C(367)=  0.257452683541916519D+09
      C(368)= -0.160941523589171052D+09
      C(369)= -0.904585295022329665D+06
      C(370)=  0.557232834051621333D+08
      C(371)= -0.124754550563966413D+07
      C(372)=  0.262367141073879510D+09
      C(373)= -0.313431115416375405D+02
      C(374)= -0.252016324080120802D+09
      C(375)= -0.128171592253549471D+09
      C(376)= -0.145647865260627389D+09
      C(377)= -0.133262423647790402D+09
      C(378)= -0.154631344724482034D+00
      C(379)=  0.872985539990470558D+08
      C(380)=  0.263473065539835572D+09
      C(381)=  0.219174160845999867D+09
      C(382)=  0.271765044027515233D+09
      C(383)= -0.412567801399783015D+09
      C(384)= -0.543098069261015058D+09
      C(385)=  0.295617012359393500D+08
      C(386)=  0.560522393169209063D+08
      C(387)= -0.119271370910373144D+08
      C(388)=  0.102276399985316500D+09
      C(389)=  0.189149636837033182D+09
      C(390)=  0.589202450211664289D+08
      C(391)= -0.636684260378074124D+08
      C(392)=  0.112785810944169760D+09
      C(393)=  0.191322321804751568D+08
      C(394)= -0.827447816885721236D+08
      C(395)= -0.680053456737457309D+07
      C(396)= -0.203293093964415928D+07
      C(397)= -0.697415291504790187D+08
      C(398)= -0.139415509498597443D+09
      C(399)= -0.544735221539729312D+08
      C(400)=  0.253551927002947092D+09
      C(401)= -0.353851427073399842D+09
      C(402)=  0.450888193387679756D+09
      C(403)=  0.414691802057749033D+09
      C(404)=  0.344683624661482424D+08
      C(405)=  0.100880718967824340D+09
      C(406)= -0.184340419233662218D+09
      C(407)= -0.304667375464704037D+08
      C(408)=  0.197592697815884352D+09
      C(409)= -0.332269876197082460D+09
      C(410)=  0.148708032858467966D+09
      C(411)=  0.369283710519427836D+09
      C(412)=  0.151189874933755606D+09
      C(413)= -0.115383613186886025D+06
      C(414)= -0.137731828802842766D+09
      C(415)= -0.184865443299308708D+03
      C(416)=  0.103172575777985230D+08
      C(417)=  0.975399909062853158D+08
      C(418)=  0.264746480686026454D+09
      C(419)= -0.259932461081043810D+09
      C(420)= -0.885639479687149823D+08
      C(421)=  0.363147510430193841D+08
      C(422)=  0.126817185677753925D+09
      C(423)= -0.510163435198088609D+04
      C(424)= -0.667340406552204013D+09
      C(425)= -0.564450139277918786D+08
      C(426)= -0.455830877860542297D+09
      C(427)= -0.170043434399561942D+09
      C(428)= -0.767732878260688186D+08
      C(429)=  0.108733610810910329D+09
      C(430)=  0.371343610708439350D+08
      C(431)= -0.168553369879585892D+09
      C(432)=  0.480785222433780968D+09
      C(433)=  0.380843714128129244D+09
      C(434)=  0.163505642454762965D+09
      C(435)= -0.464583787840582252D+09
      C(436)= -0.223344128561521694D+08
      C(437)= -0.126549295857765496D+09
      C(438)=  0.224415195791764498D+09
      C(439)= -0.361180738317449629D+09
      C(440)= -0.159867928394387484D+09
      C(441)= -0.244395142787895977D+09
      C(442)= -0.176733907014004320D+09
      C(443)=  0.347954267905110419D+09
      C(444)=  0.729548032385139763D+08
      C(445)= -0.132044448496250272D+09
      C(446)= -0.101664891482720785D+08
      C(447)=  0.319046503327601194D+09
      C(448)=  0.253274818839847982D+09
      C(449)= -0.541578716269258261D+09
      C(450)= -0.397242725268816710D+09
      C(451)=  0.102971175967201829D+09
      C(452)= -0.448269167338930905D+09
      C(453)= -0.274118075098442912D+09
      C(454)= -0.219930325581451327D+09
      C(455)= -0.525512375783562243D+09
      C(456)= -0.561826198425127387D+09
      C(457)=  0.197736421511560827D+09
      C(458)= -0.103107333226767106D+01
      C(459)=  0.170874128622146964D+09
      C(460)=  0.131070455383041531D+09
      C(461)= -0.385987960670464993D+09
      C(462)=  0.386218900928325117D+09
      C(463)=  0.365985496005640090D+09
      C(464)= -0.585751831155496649D+06
      C(465)= -0.379636025104878604D+09
      C(466)=  0.433868608205206871D+09
      C(467)= -0.275170597268814027D+09
      C(468)= -0.340636891501707956D+08
      C(469)=  0.355965221080349445D+09
      C(470)=  0.469506228951948643D+09
      C(471)=  0.110287198305887341D+09
      C(472)= -0.174877302083864152D+09
      C(473)= -0.162606531361586094D+09
      C(474)=  0.106112306809511319D+09
      C(475)=  0.982287691051959157D+09
      C(476)= -0.241492405341778606D+09
      C(477)= -0.335843355455635667D+09
      C(478)= -0.798354559319236279D+09
      C(479)= -0.619977623993471861D+09
      C(480)= -0.734402653322786570D+09
      C(481)= -0.726458006682692766D+09
      C(482)= -0.532722147046656087D+08
      C(483)=  0.423217201609222963D+08
      C(484)=  0.820299090645987093D+08
      C(485)=  0.107954255328184295D+10
      C(486)=  0.405461879367607415D+09
      C(487)=  0.406083037851314306D+09
      C(488)=  0.119253277572059661D+09
      C(489)=  0.318114622228210270D+09
      C(490)=  0.781317937159617897D+06
      C(491)=  0.479923365173723642D+07
      C(492)=  0.327113665767753720D+09
      C(493)=  0.396089098352810323D+09
      C(494)=  0.500864736179708969D+07
      C(495)= -0.581158333819566727D+09
      C(496)= -0.473968549301478863D+09
      C(497)= -0.390429593525136530D+09
      C(498)=  0.562138532518208865D+07
      C(499)= -0.876334910154533863D+09
      C(500)=  0.430317970836516023D+09
      C(501)= -0.134398491073701084D+09
      C(502)=  0.766375019305797666D+07
      C(503)= -0.142553051520190030D+09
      C(504)=  0.228211199085313857D+09
      C(505)=  0.118776378011460099D+08
      C(506)= -0.437117119518903136D+09
      C(507)=  0.161514821328624415D+10
      C(508)= -0.119504690659863204D+09
      C(509)= -0.374646837731896937D+09
      C(510)=  0.136426656389880270D+09
      C(511)=  0.282883618674007714D+09
      C(512)= -0.672721843924823165D+09
      C(513)= -0.848169348835755885D+08
      C(514)=  0.857917520579238772D+09
      C(515)= -0.112110588066887960D+09
      C(516)= -0.170510699957159281D+09
      C(517)=  0.315824797850696206D+09
      C(518)=  0.526705052495776653D+09
      C(519)=  0.109755787667413926D+10
      C(520)=  0.234481909425346941D+09
      C(521)=  0.136981648223553038D+10
      C(522)=  0.425088248874861717D+09
      C(523)=  0.857035754703413010D+09
      C(524)= -0.893088360744966984D+09
      C(525)= -0.104124370768322647D+10
      C(526)= -0.740861325992428780D+09
      C(527)=  0.493505433196177304D+09
      C(528)= -0.527451199624193966D+09
      C(529)=  0.519267540552635193D+08
      C(530)=  0.321317702437673807D+09
      C(531)=  0.947426978938322216D+08
      C(532)=  0.706350896608492136D+09
      C(533)=  0.191246737329526782D+10
      C(534)= -0.297333102545530200D+09
      C(535)= -0.962339741861212850D+09
      C(536)= -0.141572668790918779D+10
      C(537)= -0.385400447775380104D+06
      C(538)= -0.191088384556865543D+09
      C(539)=  0.599000684868682742D+09
      C(540)= -0.402661252166199964D+06
      C(541)= -0.483464966256011307D+09
      C(542)= -0.467146374238409698D+09
      C(543)=  0.167556332489675307D+10
      C(544)= -0.776090337040742159D+09
      C(545)=  0.568368165286387563D+09
      C(546)=  0.988867236447566152D+09
      C(547)= -0.227492514211942144D+08
      C(548)= -0.142619190028003544D+09
      C(549)=  0.201474377923165530D+09
      C(550)= -0.108840370218646812D+10
      C(551)= -0.921227697939270735D+08
      C(552)=  0.809621539461791217D+08
      C(553)= -0.678213474138400197D+09
      C(554)= -0.461026319192961633D+09
      C(555)=  0.525238044709189951D+09
      C(556)=  0.848937081977301359D+09
      C(557)=  0.758803310161771894D+09
      C(558)=  0.995393172631168485D+09
      C(559)=  0.284000289772026360D+09
      C(560)=  0.987510842373637080D+09
      C(561)=  0.315029088409410743D+07
      C(562)= -0.139723850381990014D+05
      C(563)= -0.411835225411671877D+09
      C(564)= -0.247874804899713039D+10
      C(565)=  0.758255862800950766D+09
      C(566)= -0.154985624364118780D+04
      C(567)=  0.221885345511809206D+10
      C(568)= -0.451726153618574515D+08
      C(569)= -0.106553178704491246D+10
      C(570)=  0.125250912700743005D+09
      C(571)= -0.163326975703739977D+10
      C(572)= -0.295057475395786995D+07
      C(573)=  0.489609437947264731D+09
      C(574)= -0.104480353098022994D+08
      C(575)=  0.177438771495141602D+10
      C(576)= -0.145849339846289635D+09
      C(577)=  0.380004235404810369D+09
      C(578)= -0.604682094324273109D+09
      C(579)=  0.361723957230835687D+06
      C(580)=  0.269013417368267715D+09
      C(581)= -0.114863978995043182D+10
      C(582)= -0.206099535300812411D+10
      C(583)= -0.977790379838855863D+09
      C(584)=  0.948247786730784893D+09
      C(585)=  0.958624283868839383D+09
      C(586)= -0.189960158565461755D+10
      C(587)= -0.238306959102375269D+10
      C(588)= -0.190662103040271625D+07
      C(589)=  0.715093884045501351D+09
      C(590)= -0.907758548083255053D+09
      C(591)=  0.603500347737260103D+09
      C(592)= -0.120979148418599695D+09
      C(593)=  0.205595159336935310D+07
      C(594)= -0.637165084313310623D+09
      C(595)=  0.720166631656499654D+08
      C(596)=  0.946854984874647856D+09
      C(597)=  0.612400837568606138D+09
      C(598)= -0.967842768961844563D+09
      C(599)= -0.143831561495054245D+10
      C(600)=  0.603450072919735074D+09
      C(601)=  0.448742859583371699D+09
      C(602)=  0.653956867393401504D+09
      C(603)=  0.814079456386528257D+07
      C(604)=  0.983219302382661581D+09
      C(605)= -0.314196980145991035D+08
      C(606)=  0.144975621336531162D+10
      C(607)= -0.259648384901462615D+09
      C(608)=  0.674575721288320899D+09
      C(609)= -0.149574113360720277D+09
      C(610)= -0.688548109557763696D+09
      C(611)=  0.712572337045879245D+09
      C(612)=  0.334948610566188872D+09
      C(613)= -0.263535441608430535D+09
      C(614)= -0.205479517167684937D+10
      C(615)=  0.756075260572801977D+08
      C(616)= -0.212300244718800068D+09
      C(617)=  0.565522674724939913D+08
      C(618)=  0.657884970148748040D+09
      C(619)= -0.566556417862892151D+09
      C(620)=  0.532472936937214613D+09
      C(621)= -0.465405455029226720D+09
      C(622)= -0.979201460348189950D+09
      C(623)= -0.212391873842859358D+09
      C(624)= -0.643704229661016107D+09
      C(625)= -0.750291616825325727D+09
      C(626)=  0.434728043380463541D+09
      C(627)= -0.472989600266833156D+08
      C(628)=  0.333102312986664951D+09
      C(629)=  0.154501073146552825D+10
      C(630)=  0.125999477502544284D+10
      C(631)= -0.122269303449967170D+10
      C(632)= -0.250683096502883993D+08
      C(633)= -0.293808769735802412D+10
      C(634)=  0.169304292751091421D+09
      C(635)=  0.250179480480356026D+10
      C(636)=  0.240567264115661287D+10
      C(637)=  0.657276814463840961D+09
      C(638)=  0.348358281614016831D+09
      C(639)= -0.407903818482457777D+06
      C(640)= -0.150171981021067882D+10
      C(641)= -0.412451248421922827D+10
      C(642)=  0.134385809146550751D+10
      C(643)= -0.420813153842235744D+09
      C(644)=  0.196638962675836372D+10
      C(645)= -0.891762518555347919D+09
      C(646)= -0.151488087464771318D+10
      C(647)= -0.308849595098289251D+10
      C(648)= -0.206664099975099039D+10
      C(649)= -0.101950240622414780D+10
      C(650)=  0.273095601567285769D+08
      C(651)=  0.420532635495630836D+10
      C(652)=  0.851852657863470793D+09
      C(653)= -0.259241072845063972D+10
      C(654)= -0.619600581373873726D+07
      C(655)=  0.256860261243888903D+10
      C(656)=  0.287224630788667488D+10
      C(657)=  0.175009458145557594D+10
      C(658)=  0.142172180584420156D+10
      C(659)=  0.222832370677686024D+10
      C(660)= -0.903907280132698417D+08
      C(661)=  0.100448752910501206D+10
      C(662)= -0.424320992468694508D+09
      C(663)=  0.195746270238362312D+10
      C(664)= -0.973970959846331596D+09
      C(665)= -0.680808028501452565D+09
      C(666)= -0.763340181374518633D+09
      C(667)= -0.289915319310639381D+09
      C(668)= -0.143898512211139417D+10
      C(669)= -0.110689332587893128D+10
      C(670)= -0.399939552496202171D+09
      C(671)= -0.751663164418962121D+09
      C(672)= -0.106776659807155013D+10
      C(673)= -0.980970701208772540D+09
      C(674)= -0.839660710905282378D+09
      C(675)=  0.285053565201601088D+09
      C(676)=  0.342607832743341625D+09
      C(677)= -0.468571637780271471D+09
      C(678)= -0.849974629520948052D+09
      C(679)=  0.124540587484237885D+10
      C(680)=  0.632927925207293510D+09
      C(681)=  0.381702479254240468D+08
      C(682)=  0.771554924440482855D+09
      C(683)=  0.573601472353042006D+09
      C(684)=  0.126274255377684264D+01
      C(685)= -0.164488228186675123D+01
      C(686)= -0.299903954629339518D-03
      C(687)=  0.232123854751777881D+01
      C(688)=  0.164711667732665629D+00
      C(689)= -0.142831935721023329D+01
      C(690)=  0.468722062016370677D-01
      C(691)=  0.624098340914674005D-01
      C(692)=  0.155882250676297601D+00
      C(693)= -0.997325635887207864D+04
      C(694)=  0.719928154703061707D+00
      C(695)=  0.915674932908943973D+00
      C(696)=  0.319158044948354203D+01
      C(697)=  0.713873621644589340D+00
      C(698)=  0.470961194131425509D+00
      C(699)=  0.699656896417086505D+00
      C(700)=  0.193122229129305234D+01
      C(701)=  0.155841904041989898D+01
      C(702)=  0.104013960809273964D+01
      C(703)=  0.414000505937959229D+00
      C(704)=  0.115902093262765749D+01
      C(705)=  0.129776430835251499D+01

      BSORDER(1)=10

      POLORDER=26

      DO I=1,DEG
        NCBAS(I)=0.0D+00
        NCBAS(I)=2*BSORDER(I)+2
      END DO

      CALL POLNC(POLORDER,MOLTYP,NCPOL)

      IF (NC.NE.(NCPOL+SUM(NCBAS))) THEN 
        WRITE(*,*) "PROBLEMS IN DEFINING PROPER NUMBER OF COEFFICIENTS"
        WRITE(*,*) "TOTAL NUMBER OF COEFFS ARE NOT SUMMING UP CORRECTLY"
        STOP
      END IF

      ALLOCATE(POL(NCPOL))

      DO I=1,NCPOL
        POL(I)=0.00D+00
      END DO

      DO I=1,NX
        Y(I)=0.0D+00
      END DO

      SUMC=NCPOL
      DO I=1,DEG
        ALLOCATE(BS(NCBAS(I)))
        DO J=1,NCBAS(I)
          BS(J)=0.0D+00
        END DO
        K=SUMC
        DO J=1,NCBAS(I)
          K=K+1
          BS(J)=C(K)
        END DO

        IF (DEG.EQ.1) THEN
          DO O=1,3 
            CALL BASIS_CONTRACT(1,BSORDER(1),BS,R(O),Y(O))
          END DO
        END IF

        SUMC=SUMC+NCBAS(I)
        DEALLOCATE(BS)
      END DO

      TOTNUM=0
      S=0
      DO M=0,POLORDER
        MNUM=0
        DO I=0,M
          DO J=0,(M-I)
            K=M-I-J
            S=I+J+K
            IF (S.NE.I .AND. S.NE.J .AND. S.NE.K) THEN

              IF (MOLTYP.EQ."ABC") THEN  

                TOTNUM=TOTNUM+1
                MNUM=MNUM+1

                CALL PERMUTABC(I,J,K,P,ID)

                DO L=1,ID

                  POL(TOTNUM)=POL(TOTNUM)+
     *               (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

                END DO

                POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

              ELSE IF (MOLTYP.EQ."AB2") THEN 

                IF (J.LE.K) THEN 

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTAB2(I,J,K,P,ID)
                   
                  DO L=1,ID

                    POL(TOTNUM)=POL(TOTNUM)+
     *                 (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

                  END DO

                  POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

                END IF    
         
              ELSE IF (MOLTYP.EQ."A3") THEN

                IF (I.LE.J .AND. J.LE.K) THEN

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTA3(I,J,K,P,ID)
                   
                  DO L=1,ID

                    POL(TOTNUM)=POL(TOTNUM)+
     *                 (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

                  END DO

                  POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

                END IF

              END IF

            END IF
          END DO
        END DO

      END DO

      POT=SUM(POL)*REPDAMP(NX,R)

      DEALLOCATE(POL)

      RETURN
      END

!######################################################################################################
! C2 SIGMA G PLUS SINGLET EXPERIMENTALLY-DETERMINED POTENTIAL ENERGY CURVE
!######################################################################################################

      SUBROUTINE CHIPR_C2SIGPS_EXP(R,POT)
      IMPLICIT NONE
      INTEGER, PARAMETER :: NC=18
      INTEGER :: I,J
      INTEGER :: BSORDER
      INTEGER :: POLORDER
      INTEGER :: NCBAS,NCPOL
      DOUBLE PRECISION, DIMENSION(2) :: Z
      DOUBLE PRECISION, DIMENSION(NC) :: C
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BS
      DOUBLE PRECISION :: Y   
      DOUBLE PRECISION :: R,POT

      Z( 1)= 0.600000000000000000D+01
      Z( 2)= 0.600000000000000000D+01

      C( 1)=-0.267086505870802468D-02
      C( 2)=-0.285280637765205680D+00
      C( 3)= 0.592287240803180115D+01
      C( 4)=-0.558544049425728630D+02
      C( 5)= 0.235039555064813953D+03
      C( 6)=-0.483801598735906566D+03
      C( 7)= 0.448105577174231883D+03
      C( 8)=-0.118126253532500471D+03
      C( 9)= 0.741626089532241739D+00
      C(10)=-0.177582090569614409D+00
      C(11)=-0.603744429598234399D-02
      C(12)=-0.225044902042792103D+03
      C(13)= 0.548451216510895345D+00
      C(14)= 0.935246919782962038D+00
      C(15)= 0.615983736625717704D+00
      C(16)= 0.412112345538274116D+01
      C(17)= 0.146596049324467104D+01
      C(18)= 0.946493553516430497D+00

      BSORDER=4

      POLORDER=8

      NCBAS=2*BSORDER+2

      NCPOL=POLORDER

      ALLOCATE(BS(NCBAS))

      DO I=1,NCBAS
        BS(I)=0.00D+00
      END DO

      DO I=1,NCBAS
        J=NCPOL+I
        BS(I)=C(J)
      END DO
 
      Y=0.00D+00
      POT=0.00D+00

      CALL BASIS_CONTRACT(1,BSORDER,BS,R,Y)

      DO I=1,POLORDER
        POT=POT+(Z(1)*Z(2)/R)*C(I)*Y**(DBLE(I))
      END DO

      DEALLOCATE(BS)

      RETURN
      END  

!######################################################################################################
! C2 PI U TRIPLET EXPERIMENTALLY-DETERMINED POTENTIAL ENERGY CURVE
!######################################################################################################

      SUBROUTINE CHIPR_C2PIUT_EXP(R,POT)
      IMPLICIT NONE
      INTEGER, PARAMETER :: NC=16
      INTEGER :: I,J
      INTEGER :: BSORDER
      INTEGER :: POLORDER
      INTEGER :: NCBAS,NCPOL
      DOUBLE PRECISION, DIMENSION(2) :: Z
      DOUBLE PRECISION, DIMENSION(NC) :: C
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BS
      DOUBLE PRECISION :: Y   
      DOUBLE PRECISION :: R,POT

      Z( 1)= 0.600000000000000000D+01
      Z( 2)= 0.600000000000000000D+01

      C( 1)= 0.673199875363553191D-01
      C( 2)= 0.354670818253029277D-01
      C( 3)= 0.454530262577290817D-01
      C( 4)= 0.896638948685439183D-02
      C( 5)=-0.557192541182933640D-01
      C( 6)= 0.396468034787315002D-01
      C( 7)= 0.137634242346898183D+01
      C( 8)=-0.644084104783897993D+00
      C( 9)= 0.415028815178197785D-01
      C(10)=-0.580664867995453619D+03
      C(11)= 0.132718504668728920D+01
      C(12)= 0.941815085429398735D+00
      C(13)= 0.108834374997582084D+01
      C(14)= 0.765033716358887284D+00
      C(15)= 0.240379211673489301D+01
      C(16)= 0.799892875462521968D+00

      BSORDER=4

      POLORDER=6

      NCBAS=2*BSORDER+2

      NCPOL=POLORDER

      ALLOCATE(BS(NCBAS))

      DO I=1,NCBAS
        BS(I)=0.00D+00
      END DO

      DO I=1,NCBAS
        J=NCPOL+I
        BS(I)=C(J)
      END DO
 
      Y=0.00D+00
      POT=0.00D+00

      CALL BASIS_CONTRACT(1,BSORDER,BS,R,Y)

      DO I=1,POLORDER
        POT=POT+(Z(1)*Z(2)/R)*C(I)*Y**(DBLE(I))
      END DO

      DEALLOCATE(BS)

      RETURN
      END 

!######################################################################################################
