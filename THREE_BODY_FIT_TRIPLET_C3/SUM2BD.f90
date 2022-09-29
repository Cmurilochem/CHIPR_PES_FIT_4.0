!###############################################################################
! THIS IS A USER SUPPLIED SUBROUTINE THAT CALCULATES THE SUM OF TWO-DOBY ENERGIES
! FOR A GIVEN TRIATOMIC MOLECULE
!###############################################################################

      SUBROUTINE SUM2BD(R1,R2,R3,SUMV2)
      IMPLICIT NONE
      INTEGER :: I
      DOUBLE PRECISION :: R1,R2,R3,SUMV2
      DOUBLE PRECISION, DIMENSION(3) :: V1,V2,VTOT
      DOUBLE PRECISION, DIMENSION(3) :: RE
      DOUBLE PRECISION, DIMENSION(3) :: RJAC
      DOUBLE PRECISION, DIMENSION(3) :: RB,RC
      DOUBLE PRECISION, DIMENSION(3) :: COSGAMMA
      DOUBLE PRECISION :: DIAT1,DIAT2,DIAT3
      DOUBLE PRECISION :: X0,GAM,RANG2BD
      V1(1)=DIAT1(R1)
      V1(2)=DIAT2(R2)
      V1(3)=DIAT3(R3)
      CALL CHIPR_C2SIGPS(R1,V2(1))
      CALL CHIPR_C2SIGPS(R2,V2(2))
      CALL CHIPR_C2SIGPS(R3,V2(3))
      CALL TRIANG2JACOBI(R1,R2,R3,RE,RJAC,COSGAMMA,RB,RC)
      X0=5.0D+00
      GAM=1.0D+00
      DO I=1,3
        VTOT(I)=(V2(I)+(V1(I)-V2(I))*RANG2BD(RJAC(I),X0,GAM))
      END DO
      SUMV2=SUM(VTOT)
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
      DOUBLE PRECISION :: R1,V21
!
!     DEFINE HERE YOUR SUBROUTINE FOR DIATOMIC R1 
!
      CALL TWOBD(R1,V21)

      DIAT1=V21

      RETURN 
      END

!###############################################################################

      DOUBLE PRECISION FUNCTION DIAT2(R2)
      IMPLICIT NONE
      DOUBLE PRECISION :: R2,V22
!
!     DEFINE HERE YOUR SUBROUTINE FOR DIATOMIC R2
!
      CALL TWOBD(R2,V22)

      DIAT2=V22

      RETURN 
      END

!###############################################################################

      DOUBLE PRECISION FUNCTION DIAT3(R3)
      IMPLICIT NONE
      DOUBLE PRECISION :: R3,V23
!
!     DEFINE HERE YOUR SUBROUTINE FOR DIATOMIC R3 
!
      CALL TWOBD(R3,V23)

      DIAT3=V23

      RETURN 
      END

!######################################################################################################
! RANG FUNCTION THAT CHANGES ONE DIATOMIC TO THE OTHER
!######################################################################################################

      DOUBLE PRECISION FUNCTION RANG2BD(R,X0,GAM)
      IMPLICIT NONE
      DOUBLE PRECISION :: R,X0,GAM
      RANG2BD=0.50D+00*(1.00D+00+TANH(GAM*(R-X0)))
      END FUNCTION RANG2BD

!######################################################################################################
! CONVERTING TRIANGULAR TO JACOBI COORDINATES
!######################################################################################################

      SUBROUTINE TRIANG2JACOBI(R1,R2,R3,RE,RJAC,COSGAMMA,RB,RC)
      IMPLICIT NONE
      DOUBLE PRECISION :: R1,R2,R3
      DOUBLE PRECISION, DIMENSION(3) :: RE
      DOUBLE PRECISION, DIMENSION(3) :: RJAC
      DOUBLE PRECISION, DIMENSION(3) :: RB,RC
      DOUBLE PRECISION, DIMENSION(3) :: SQUAREDRJAC
      DOUBLE PRECISION, DIMENSION(3) :: COSGAMMA
      DOUBLE PRECISION :: PI,ZERO
      INTEGER :: I

      PI=4.0D0*DATAN(1.0D+00)
      ZERO=1.0D-12

      DO I=1,3

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

            RJAC(I)=SQRT(SQUAREDRJAC(I)) 

          END IF
    
        COSGAMMA(I)=(0.5D+00*(RC(I)**2-RB(I)**2))/(RJAC(I)*RE(I))

          IF (COSGAMMA(I)>1.0D+00) THEN 

            COSGAMMA(I)=1.0D+00 

          ELSE IF (COSGAMMA(I)<-1.0D+00) THEN 
   
            COSGAMMA(I)=-1.0D+00

          END IF 

      END DO

      END SUBROUTINE TRIANG2JACOBI

!######################################################################################################
! TOTAL TWO-BODY POTENTIAL ENERGY CURVE
!######################################################################################################

      SUBROUTINE TWOBD(R,POT)
      IMPLICIT NONE
      INTEGER :: I,J 
      INTEGER, PARAMETER :: M=1
      DOUBLE PRECISION, DIMENSION(6*M+3) :: C
      DOUBLE PRECISION, DIMENSION(3,3) :: H,VEC
      DOUBLE PRECISION, DIMENSION(3) :: EIGVAL
      DOUBLE PRECISION :: R,POT,SECH
      DOUBLE PRECISION :: POT1,POT2,POT3
      DOUBLE PRECISION :: RC12,RC23,RC13,RHO

      RHO=1.0D-08

      C(1)=-0.491457630094914899D-07
      C(2)= 0.157491862773895264D+00
      C(3)= 0.543041520945708385D-04
      C(4)= 0.634928047657012939D+00
      C(5)=-0.140642678423768619D-02
      C(6)= 0.157993710041046143D+01

      C(6*M+1)=2.4474817237408617D+00
      C(6*M+2)=3.1694695847347925D+00
      C(6*M+3)=2.6955705977852991D+00

      DO I=1,3
        DO J=1,3
          H(I,J)=0.0D+00
        END DO
      END DO 

      CALL CHIPR_C2SIGPS(R,POT1)
      CALL CHIPR_C2PIUT(R,POT2)
      CALL CHIPR_C2SIGMT(R,POT3)

      H(1,1)=POT1
      H(2,2)=POT2
      H(3,3)=POT3

      RC12=C(6*M+1)
      RC23=C(6*M+2)
      RC13=C(6*M+3)

      DO I=1,M
        H(1,2)=H(1,2)+&
     &  0.5D+00*C(I)*C(I+M)*SECH(C(I+M)*(R-RC12))
      END DO

      H(2,1)=H(1,2)

      DO I=2*M+1,3*M
        H(2,3)=H(2,3)+&
     &  0.5D+00*C(I)*C(I+M)*SECH(C(I+M)*(R-RC23))
      END DO

      H(3,2)=H(2,3)

      DO I=4*M+1,5*M
        H(1,3)=H(1,3)+&
     &  0.5D+00*C(I)*C(I+M)*SECH(C(I+M)*(R-RC13))
      END DO

      H(3,1)=H(1,3)

      CALL DIAG(H,EIGVAL,VEC,3,RHO)

      POT=EIGVAL(1)

      RETURN
      END

!###############################################################################

      SUBROUTINE DIAG(A,EIG,VEC,NN,RHO)
      IMPLICIT NONE
      INTEGER, PARAMETER :: NDIM=3
      INTEGER :: I,I1,J,N,NR,NN,N1,N2,NPAS,NT 
      INTEGER :: ITEMP,LV,NRR,K,L,M,NP,NV,QJ  
      DOUBLE PRECISION :: EIG(NDIM)
      DOUBLE PRECISION :: W(NDIM),BETASQ(NDIM)
      DOUBLE PRECISION :: GAMMA(NDIM),BETA(NDIM)
      DOUBLE PRECISION :: A(NDIM,NDIM),VEC(NDIM,NDIM)
      DOUBLE PRECISION :: P(NDIM),Q(NDIM),IPOSV(NDIM)
      DOUBLE PRECISION :: IVPOS(NDIM),IORD(NDIM)
      DOUBLE PRECISION :: RHOSQ,RHO,SHIFT,S,SGN,D,TEMP,B
      DOUBLE PRECISION :: SQRTS,WTAW,SUM,COSA,G,PP,PPBS,PPBR
      DOUBLE PRECISION :: COSAP,SINA,SINA2,DIA,U,A2,R2,R1,R12
      DOUBLE PRECISION :: DIF,WJ
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

!######################################################################################################
! C2 PI U TRIPLET POTENTIAL ENERGY CURVE
!######################################################################################################

      SUBROUTINE CHIPR_C2PIUT(R,POT)
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

      C( 1)= 0.688368821899691119D-01
      C( 2)= 0.342351245803137610D-01
      C( 3)= 0.154614262783052108D-01
      C( 4)= 0.869956989359632789D-02
      C( 5)= 0.805774584332814674D-02
      C( 6)=-0.155565331355019400D-02
      C( 7)= 0.140576516199323032D+01
      C( 8)=-0.651695535962590333D+00
      C( 9)= 0.415028815178197785D-01
      C(10)=-0.580664867995453619D+03
      C(11)= 0.131877461311567812D+01
      C(12)= 0.950567375315507390D+00
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
! C2 SIGMA G SINGLET POTENTIAL ENERGY CURVE
!######################################################################################################

      SUBROUTINE CHIPR_C2SIGPS(R,POT)
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

      C( 1)=-0.156080848440132039D-02
      C( 2)=-0.146090727319725361D+00
      C( 3)= 0.289447560018073835D+01
      C( 4)=-0.329128960443330456D+02
      C( 5)= 0.149172064999103469D+03
      C( 6)=-0.310764847510739344D+03
      C( 7)= 0.267631937984524939D+03
      C( 8)=-0.414922464726189091D+02
      C( 9)= 0.741626089532241739D+00
      C(10)=-0.176756629350128786D+00
      C(11)=-0.724212331201012048D-02
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
! C2 SIGMA G SINGLET POTENTIAL ENERGY CURVE
!######################################################################################################

      SUBROUTINE CHIPR_C2SIGMT(R,POT)
      IMPLICIT NONE
      INTEGER, PARAMETER :: NC=20
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

      C( 1)=-0.139383051089885024D+00
      C( 2)= 0.379894280261107440D+02
      C( 3)=-0.868777870230472763D+04
      C( 4)= 0.617478594715823303D+06
      C( 5)=-0.224837692867042422D+08
      C( 6)= 0.442548539070751369D+09
      C( 7)=-0.431210636544939327D+10
      C( 8)= 0.148558838775511398D+11
      C( 9)=-0.317605248014272630D+09
      C(10)= 0.374148111458616577D+12
      C(11)= 0.266819262786576694D+00
      C(12)=-0.210741782847390147D+00
      C(13)= 0.101634586991167734D-01
      C(14)=-0.393424499706291385D+03
      C(15)= 0.298588226648417887D+00
      C(16)= 0.299430688776365861D+00
      C(17)= 0.799635468897632928D+00
      C(18)= 0.439191465042685469D+00
      C(19)= 0.158695418480965822D+01
      C(20)= 0.138708213744285547D+01

      BSORDER=4

      POLORDER=10

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
