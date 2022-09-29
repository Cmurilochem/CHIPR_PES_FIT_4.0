!###############################################################################
! MODULE THAT CALCULATES PROPERTIES FOR DIATOMICS
!###############################################################################
 
      SUBROUTINE PROP_DIAT
      USE COMMON_VAR
      IMPLICIT NONE
      INTEGER :: T,V,U,NPC,SUMC,O,E,DUM
      DOUBLE PRECISION :: BAS,BASNEW,BASOLD
      DOUBLE PRECISION :: R,POT,STEPRD
      DOUBLE PRECISION :: RDMIN,RDMAX
      DOUBLE PRECISION :: ORIG,R0,TWOBD,THREBD
      DOUBLE PRECISION :: R0NEW,R0OLD,REQ,DEQ,BAND
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: CINT,XINT,YINT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BS,BSNEW,BSOLD
      LOGICAL :: CHECK
      DOUBLE PRECISION, DIMENSION(NX) :: ARG
      DOUBLE PRECISION :: RE,DE,SD,FCT,VIB
      DOUBLE PRECISION, DIMENSION(NP) :: P
      DOUBLE PRECISION, DIMENSION(NC) :: DC
      INTEGER :: IIV(IVIBMX),IIJ(IVIBMX)
      DOUBLE PRECISION :: ESOLN(IVIBMX),IROTMAX(IVIBMX)
      DOUBLE PRECISION :: ROTCTE(7,IVIBMX),ZPE
      DOUBLE PRECISION :: ALLROVIB(0:IVIBMX,0:IVIBMX)
      INTEGER :: IAN1,IMN1,IAN2,IMN2,NLEV1,NJM,JDJR,IVMAX
 
      OPEN(UNIT=18,FILE='dev_basis_abinitio_points.res',IOSTAT=ERR,&
     &FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
      IF (ERR /= 0) THEN
      PRINT*,"PROBLEMS IN OPENING THE FILE DEV_BASIS_ABINITIO_POINTS.RES"
      END IF

      OPEN(UNIT=19,FILE='plot_basis.res',IOSTAT=ERR,&
     &FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
      IF (ERR /= 0) THEN
      PRINT*,"PROBLEMS IN OPENING THE FILE PLOT_BASIS.RES"
      END IF

      OPEN(UNIT=20,FILE='plot_basis_abinitio.res',IOSTAT=ERR,&
     &FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
      IF (ERR /= 0) THEN
      PRINT*,"PROBLEMS IN OPENING THE FILE PLOT_BASIS_ABINITIO.RES"
      END IF

      OPEN(UNIT=21,FILE='plot_nodes.res',IOSTAT=ERR,&
     &FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
      IF (ERR /= 0) THEN
      PRINT*,"PROBLEMS IN OPENING THE FILE PLOT_NODES.RES"
      END IF

      IF (JOBTYP.EQ."OPTBASIS") THEN
!
!     WRITTING FINAL DEVATIONS IN THE FILE DEV_BASIS_ABINITIO_POINTS.RES
!     WRITTING AB INITIO BASIS IN THE FILE PLOT_BASIS_ABINITIO.RES
!
        DO T=1,CONTDEG
          NC=NCBS(T)
          NP=NPBS(T)
          M=MBS(T)
          ALLOCATE(CINT(NC),XINT(NP),YINT(NP))
          DO V=1,NC
            CINT(V)=0.0D+00
          END DO
          DO V=1,NP
            XINT(V)=0.0D+00
            YINT(V)=0.0D+00
          END DO  
          DO V=1,NC
            CINT(V)=CBS(V,T)
          END DO 
          DO U=1,NP 
            XINT(U)=XBS(T,U,T)
            YINT(U)=YBS(U,T) 
            CALL BASIS_CONTRACT(T,M,CINT,XINT(U),BAS)
            WRITE(18,*) XINT(U),(YINT(U)-BAS)*EH2CM
            WRITE(20,*) XINT(U),YINT(U),(YYBS(U,T)-EZERO)
          END DO
          WRITE(18,*)
          WRITE(18,*)
          WRITE(20,*)
          WRITE(20,*)
          DEALLOCATE(CINT,XINT,YINT)
        END DO
!
!     WRITTING FITTED BASIS IN THE FILE PLOT_BASIS.RES 
!
        RDMIN=0.1D+00 
        RDMAX=100.00D+00
        NPC=500

        STEPRD=(RDMAX-RDMIN)/(DBLE(NPC)-1.000D+00)

        DO T=1,CONTDEG
          NC=NCBS(T)
          M=MBS(T)
          ALLOCATE(CINT(NC))
          DO V=1,NC
            CINT(V)=0.0D+00
          END DO 
          DO V=1,NC
            CINT(V)=CBS(V,T)
          END DO 
          DO U=1,NPC
            R=RDMIN+STEPRD*(DBLE(U)-1.000D+00)
            CALL BASIS_CONTRACT(T,M,CINT,R,BAS)
            WRITE(19,*) R,BAS
          END DO
          WRITE(19,*)
          WRITE(19,*)
          DEALLOCATE(CINT)
        END DO
!
!     CALCULATING THE ORIGINS OF EACH BASIS FUNCTION
!     PUT IN THE FILE PLOT_NODES.RES 
!
        DO T=1,CONTDEG
          NC=NCBS(T)
          M=MBS(T)
          ALLOCATE(CINT(NC))
          DO V=1,NC
            CINT(V)=0.0D+00
          END DO 
          DO V=1,NC
            CINT(V)=CBS(V,T)
          END DO 
          DO U=1,M
            R0=ORIG(U,ZETAFINAL(1,T),RREFFINAL(1,T))
            CALL BASIS_CONTRACT(T,M,CINT,R0,BAS)
            WRITE(21,*) R0,BAS
          END DO
          WRITE(21,*)
          WRITE(21,*)
          DEALLOCATE(CINT)
        END DO

        CALL SYSTEM ("gnuplot PLOT_BASIS_OPT_CHIPR_DIATOM.gnu")

      ELSE IF (JOBTYP.EQ."FITPOL" .OR. JOBTYP.EQ."DIRECTFIT") THEN
      
      IF (JOBTYP.EQ."DIRECTFIT") THEN 
        DUM=NABINITIO
      ELSE
        DUM=NP
      END IF
!
!     PRINTING BASIS SET FEATURES FIRST
!     WRITTING FINAL DEVATIONS IN THE FILE DEV_BASIS_ABINITIO_POINTS.RES
!     WRITTING AB INITIO BASIS IN THE FILE PLOT_BASIS_ABINITIO.RES
!
        SUMC=NCPOL
        DO I=1,DEG
          ALLOCATE(BSNEW(NCBS(I)),BSOLD(NCBS(I)))
          DO J=1,NCBS(I)
            BSNEW(J)=0.0D+00
            BSOLD(J)=0.0D+00
          END DO
          K=SUMC
          DO J=1,NCBS(I)
            K=K+1
            BSNEW(J)=C(K)
            BSOLD(J)=CBSINITPOL(J,I)
          END DO

          DO J=1,NPBS(I)
            CALL BASIS_CONTRACT(I,MBS(I),BSOLD,XBS(I,J,I),BASOLD)
            CALL BASIS_CONTRACT(I,MBS(I),BSNEW,XBS(I,J,I),BASNEW)
            WRITE(18,*) XBS(I,J,I),(YBS(J,I)-BASOLD)*EH2CM,(YBS(J,I)-BASNEW)*EH2CM
            WRITE(20,*) XBS(I,J,I),YBS(J,I),(YYBS(J,I)-EZERO)
          END DO

          SUMC=SUMC+NCBS(I)
          WRITE(18,*)
          WRITE(18,*)
          WRITE(20,*)
          WRITE(20,*)
          DEALLOCATE(BSNEW,BSOLD)
        END DO 
!
!     WRITTING NEW AND OLD BASIS SET IN THE FILE PLOT_POT.RES 
!
        RDMIN=0.5D+00 
        RDMAX=50.00D+00
        NPC=500

        STEPRD=(RDMAX-RDMIN)/(DBLE(NPC)-1.000D+00)

        SUMC=NCPOL
        DO I=1,DEG
          ALLOCATE(BSNEW(NCBS(I)),BSOLD(NCBS(I)))
          DO J=1,NCBS(I)
            BSNEW(J)=0.0D+00
            BSOLD(J)=0.0D+00
          END DO
          K=SUMC
          DO J=1,NCBS(I)
            K=K+1
            BSNEW(J)=C(K)
            BSOLD(J)=CBSINITPOL(J,I)
          END DO

          DO J=1,NPC
            R=RDMIN+STEPRD*(DBLE(J)-1.000D+00)
            CALL BASIS_CONTRACT(I,MBS(I),BSOLD,R,BASOLD)
            CALL BASIS_CONTRACT(I,MBS(I),BSNEW,R,BASNEW)
            WRITE(19,*) R,BASOLD,BASNEW
          END DO

          SUMC=SUMC+NCBS(I)
          WRITE(19,*)
          WRITE(19,*)
          DEALLOCATE(BSNEW,BSOLD)
        END DO
!
!     CALCULATING THE ORIGINS OF EACH BASIS FUNCTION
!     PUT IN THE FILE PLOT_NODES.RES 
!
        SUMC=NCPOL
        DO I=1,DEG
          ALLOCATE(BSNEW(NCBS(I)),BSOLD(NCBS(I)))
          DO J=1,NCBS(I)
            BSNEW(J)=0.0D+00
            BSOLD(J)=0.0D+00
          END DO
          K=SUMC
          DO J=1,NCBS(I)
            K=K+1
            BSNEW(J)=C(K)
            BSOLD(J)=CBSINITPOL(J,I)
          END DO

          DO J=1,MBS(I)
            R0OLD=ORIG(J,ZETAINITPOL(1,I),RREFINITPOL(1,I))
            R0NEW=ORIG(J,ZETAFINAL(1,I),RREFFINAL(1,I))
            CALL BASIS_CONTRACT(I,MBS(I),BSOLD,R0OLD,BASOLD)
            CALL BASIS_CONTRACT(I,MBS(I),BSNEW,R0NEW,BASNEW)
            CALL CHIPR_DIAT(MBS(I),L,C,R0NEW,POT,DC)
            WRITE(21,*) R0OLD,BASOLD,R0NEW,BASNEW,POT
          END DO

          SUMC=SUMC+NCBS(I)
          WRITE(21,*)
          WRITE(21,*)
          DEALLOCATE(BSNEW,BSOLD)
        END DO 

        OPEN(UNIT=22,FILE='plot_abinitio.res',IOSTAT=ERR,&
     &  FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
        IF (ERR /= 0) THEN
        PRINT*,"PROBLEMS IN OPENING THE FILE PLOT_ABINITIO.RES"
        END IF

        OPEN(UNIT=23,FILE='plot_pot.res',IOSTAT=ERR,&
     &  FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
        IF (ERR /= 0) THEN
        PRINT*,"PROBLEMS IN OPENING THE FILE PLOT_POT.RES"
        END IF

        OPEN(UNIT=24,FILE='dev_abinitio_points.res',IOSTAT=ERR,&
     &  FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
        IF (ERR /= 0) THEN
        PRINT*,"PROBLEMS IN OPENING THE FILE DEV_ABINITIO_POINTS.RES"
        END IF

        OPEN(UNIT=25,FILE='strat_rmsd.res',IOSTAT=ERR,&
     &  FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
        IF (ERR /= 0) THEN
        PRINT*,"PROBLEMS IN OPENING THE FILE STRAT_RMSD.RES"
        END IF

        OPEN(UNIT=26,FILE='inp_data.res',IOSTAT=ERR,&
     &  FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
        IF (ERR /= 0) THEN
        PRINT*,"PROBLEMS IN OPENING THE FILE INP_DATA.RES"
        END IF

        OPEN(UNIT=27,FILE='coeffs_tab_latex.res',IOSTAT=ERR,&
     &  FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
        IF (ERR /= 0) THEN
        PRINT*,"PROBLEMS IN OPENING THE FILE COEFFS_TAB_LATEX.RES"
        END IF
        
        IF (JOBTYP.EQ."DIRECTFIT") THEN
        
          OPEN(UNIT=28,FILE='plot_levels.res',IOSTAT=ERR,&
     &    FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
          IF (ERR /= 0) THEN
          PRINT*,"PROBLEMS IN OPENING THE FILE PLOT_LEVELS.RES"
          END IF        
        
        END IF
!
!     EVALUATING THE AB INITIO POINTS
!     PUT IN THE FILE PLOT_ABINITIO.RES 
!
        DO I=1,DUM
          ARG(1)=X(I,NX)
          CALL CHIPR_DIAT(MBS,L,C,ARG,POT,DC)
          WRITE(22,*) X(I,NX),Y(I),POT
        END DO
!
!     WRITTING FINAL DEVATIONS IN THE FILE DEV_ABINITIO_POINTS.RES
!
        DO I=1,DUM
          ARG(1)=X(I,NX)
          CALL CHIPR_DIAT(MBS,L,C,ARG,POT,DC)
          WRITE(24,*) X(I,NX),(Y(I)-POT)*EH2CM
        END DO
!
!     WRITTING THE FITTED POTENTIAL IN THE FILE PLOT_POT.RES 
!
        RDMIN=0.1D+00 
        RDMAX=50.00D+00
        NPC=500

        DO I=1,NPC
          R=RDMIN+STEPRD*(DBLE(I)-1.000D+00)
          CALL CHIPR_DIAT(MBS,L,C,R,POT,DC)
          WRITE(23,*) R,POT
        END DO

        M=MBS(1)
!
!     CALCULATING EQUILIBRIUM GEOMETRY
!     DISSOCIATION ENERGY AND  
!     HARMONIC VIBRATIONAL FREQUENCY
!
!        CALL MINDEDIAT(RE,DE,C)

!        CALL DIATSECONDDERIVATIVE(RE,SD,C)

!        CALL DIATHARMFREQ(SD,FCT,VIB)

!        WRITE(6,550) RE,DE,VIB
!
!     WRITTING STRATIFIED RMSDs IN THE FILE STRAT_RMSD.RES
!
        DO I=1,DUM
          ARG(1)=X(I,NX)
          P(I)=0.0D+00
          CALL CHIPR_DIAT(MBS,L,C,ARG,P(I),DC)
        END DO

        CALL STRAT_RMSD(P,Y,W)
!
!     WRITTING FINAL COEFFICIENTS IN THE FILE COEFFS_TAB_LATEX.RES
!
        J=0
        DO I=(L+1),(L+M)
          J=J+1
          WRITE(27,'("$c_{",I3,"}$",6X,"&",1X,D25.18,1X,"\\")') J,C(I)
        END DO

        J=0
        DO I=(L+M+1),(L+2*M)
          J=J+1
          WRITE(27,'("$\gamma_{",I3,"}$",1X,"&",1X,D25.18,1X,"\\")') J,C(I)
        END DO

        WRITE(27,'("$\zeta$",8X,"&",1X,D25.18,1X,"\\")') C(L+2*M+2)

        WRITE(27,'("$R^{\rm ref}$",2X,"&",1X,D25.18,1X,"\\")') C(L+2*M+1)

        DO I=1,L
          WRITE(27,'("$C_{",I3,"}$",6X,"&",1X,D25.18,1X,"\\")') I,C(I)
        END DO
!
!     WRITTING FINAL ANALYTIC FUNCTION TO THE FILE "CHIPR_DIAT_FUNC.f90"
!
        WRITE(26,'(6X,"BSORDER=",1X,I5)') MBS

        WRITE(26,'(6X,"POLORDER=",1X,I5)') L

        WRITE(26,'(6X,"NC=",1X,I5)') NC

        DO I=1,2
          WRITE(26,'(6X,"Z(",I3,")=",1X,D42.35)') I,Z(I)
        END DO

        DO I=1,NC
          WRITE(26,'(6X,"C(",I3,")=",1X,D42.35)') I,C(I)
        END DO

        CALL SYSTEM ("chmod +x ../CHIPR-4.0_SOURCE_CODE/SCRIPT_CHIPR_DIAT.csh")
        CALL SYSTEM ("cp ../CHIPR-4.0_SOURCE_CODE/SCRIPT_CHIPR_DIAT.csh .")
        CALL SYSTEM ("./SCRIPT_CHIPR_DIAT.csh")
        CALL SYSTEM ("rm SCRIPT_CHIPR_DIAT.csh")
!
!  PRINT RESULTS FROM DIRECT FIT 
!
        IF (JOBTYP.EQ."DIRECTFIT") THEN
        
          CALL DIRECTMINDEDIAT(REQ,DEQ,C)
          
          VLIM=DEQ*EH2CM

          IAN1=INT(Z(1))
          IMN1=INT(AMU(1)) 
          IAN2=INT(Z(2))
          IMN2=INT(AMU(2)) 
!
!  THIS CALCULATES ALL ROVIBRATIONAL LEVELS UP TO DISSOCIATION (VLIM)
!
          NLEV1=-9999!-(NROVIB-1)!-9999 
          NJM=9999!0!9999
          JDJR=1
!      
!  THESE SHOULD DEFINE THE LOWEST STATE OF THE MOLECULE
!
          IIV(1)=0
          IIJ(1)=IOMEGA
!
!  CALLING LEVEL PROGRAM 
!
          CALL LEVEL(IAN1,IMN1,IAN2,IMN2,ICHARGE,&
     &               IOMEGA,VLIM,NLEV1,NJM,JDJR,&
     &               IIV,IIJ,ESOLN,ROTCTE,IVMAX,&
     &               IROTMAX,ALLROVIB,C)
!
!  DEFINING ZPE
!
          ZPE=ALLROVIB(0,IOMEGA)
          
          DO I=1,NROVIB
          
            J=NABINITIO+I
            
            BAND=(ALLROVIB(IVEXP(I),IJEXP(I))-ZPE)
            
            WRITE(28,*) BAND, Y(J), (Y(J)-BAND)
          
          END DO
          
        END IF
!
!     PLOT THE RESULTS
!
        CALL SYSTEM ("gnuplot PLOT_POL_OPT_CHIPR_DIATOM.gnu")

      END IF

      CLOSE(UNIT=18)
      CLOSE(UNIT=19)
      CLOSE(UNIT=20)
      CLOSE(UNIT=21)

      INQUIRE(UNIT=22,OPENED=CHECK)
      IF (CHECK) THEN 
        CLOSE(UNIT=22)
      ELSE
      END IF 
      INQUIRE(UNIT=23,OPENED=CHECK)
      IF (CHECK) THEN 
        CLOSE(UNIT=23)
      ELSE
      END IF 
      INQUIRE(UNIT=24,OPENED=CHECK)
      IF (CHECK) THEN 
        CLOSE(UNIT=24)
      ELSE 
      END IF
      INQUIRE(UNIT=25,OPENED=CHECK)
      IF (CHECK) THEN 
        CLOSE(UNIT=25)
      ELSE 
      END IF
      INQUIRE(UNIT=26,OPENED=CHECK)
      IF (CHECK) THEN 
        CLOSE(UNIT=26)
      ELSE 
      END IF
      INQUIRE(UNIT=27,OPENED=CHECK)
      IF (CHECK) THEN 
        CLOSE(UNIT=27)
      ELSE 
      END IF
      INQUIRE(UNIT=28,OPENED=CHECK)
      IF (CHECK) THEN 
        CLOSE(UNIT=28)
      ELSE 
      END IF  
      INQUIRE(UNIT=2020,OPENED=CHECK)
      IF (CHECK) THEN 
        CLOSE(UNIT=2020)
      ELSE 
      END IF      
      
! 550  FORMAT (/,"EQUILIBRIUM GEOMETRY (IN BOHR):",&
!     &         F15.5,/,"DISSOCIATION ENERGY  (IN HARTREE):",F12.5,/,&
!     &         "HARMONIC FREQUENCY   (IN CM-1):",F15.5,/)
 
      END SUBROUTINE PROP_DIAT

!###############################################################################
! SUBROUTINE TO FIND EQUILIBRIUM DISTANCE AND DISSOCIATION ENERGY OF A DIATOMIC
!###############################################################################

      SUBROUTINE MINDEDIAT(ROPT,DEOPT,C)
      USE COMMON_VAR, ONLY : NC,M,L
      IMPLICIT NONE
      INTEGER :: I
      INTEGER, PARAMETER :: NPOINTSPOT=500
      DOUBLE PRECISION, PARAMETER :: RMIN=1.9D+0
      DOUBLE PRECISION, PARAMETER :: RMAX=3.3D+00
      DOUBLE PRECISION, DIMENSION(NC) :: C,DC
      DOUBLE PRECISION, DIMENSION(NPOINTSPOT) :: R,GRAD
      DOUBLE PRECISION, DIMENSION(NPOINTSPOT) :: MODGRAD
      INTEGER, DIMENSION(1) :: P
      DOUBLE PRECISION :: STEP,ROPT,DEOPT,POTOPT

      STEP=(RMAX-RMIN)/(DBLE(NPOINTSPOT)-1.0D+00)

      DO I=1,NPOINTSPOT
        R(I)=RMIN+STEP*(DBLE(I)-1.0D+00)
        CALL DIATFIRSTDERIVATIVE(R(I),GRAD(I),C)
        MODGRAD(I)=ABS(GRAD(I))
      END DO

      P=MINLOC(MODGRAD)

      ROPT=R(P(1))

      CALL CHIPR_DIAT(M,L,C,ROPT,POTOPT,DC)

      DEOPT=-POTOPT

      END SUBROUTINE MINDEDIAT

!###############################################################################
! CALCULATION OF THE ANALYTIC GRADIENT OF A DIATOMIC CURVE USING CENTRAL FINITE DIFERENCE 
!###############################################################################

      SUBROUTINE DIATFIRSTDERIVATIVE(RE,FD,C)
      USE COMMON_VAR, ONLY : NC,M,L
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(NC) :: C,DC
      DOUBLE PRECISION :: RE,H,FD
      DOUBLE PRECISION :: POT1,POT2
      H=0.0001D+00
      CALL CHIPR_DIAT(M,L,C,RE-H,POT1,DC)
      CALL CHIPR_DIAT(M,L,C,RE+H,POT2,DC)
      FD=(POT2-POT1)/2.00D+00*H
      END SUBROUTINE DIATFIRSTDERIVATIVE

!###############################################################################
! CALCULATION OF THE DIATOMIC FORCE CONSTANT USING CENTRAL FINITE DIFERENCE
!###############################################################################

      SUBROUTINE DIATSECONDDERIVATIVE(RE,SD,C)
      USE COMMON_VAR, ONLY : NC,M,L,Z
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(NC) :: C,DC
      DOUBLE PRECISION :: RE,H,SD
      DOUBLE PRECISION :: POT0,POT1,POT2
      H=0.0001D+00
      CALL CHIPR_DIAT(M,L,C,RE,POT0,DC)
      CALL CHIPR_DIAT(M,L,C,RE-H,POT1,DC)
      CALL CHIPR_DIAT(M,L,C,RE+H,POT2,DC)
      SD=(POT2-2.00D+00*POT0+POT1)/H**2
      END SUBROUTINE DIATSECONDDERIVATIVE

!###############################################################################
! CALCULATION OF THE HARMONIC VIBRATIONAL FRENQUENCY OF A DIATOMIC
!###############################################################################

      SUBROUTINE DIATHARMFREQ(SD,FCT,VIB)
      USE COMMON_VAR, ONLY : NC,M,L,AMU
      IMPLICIT NONE
      DOUBLE PRECISION :: REDMASS
      DOUBLE PRECISION :: SD, FCT, VIB
      DOUBLE PRECISION, PARAMETER :: C=2.99792458D+10             ! LIGHT VELOCIT IN cm.s-1
      DOUBLE PRECISION, PARAMETER :: B2A=0.5291772086D+00         ! BOHR TO ANGS
      DOUBLE PRECISION, PARAMETER :: JH=4.3597439D-18             ! CONVERTION FROM HARTREES TO JOULES
      DOUBLE PRECISION, PARAMETER :: AUKG=1.660538921E-27         ! CONVERTION FROM A.U TO Kg
      DOUBLE PRECISION, PARAMETER :: PI2=6.28318530718D+00        ! PI SQUARED
      REDMASS=((AMU(1)*AMU(2))/(AMU(1)+AMU(2)))*AUKG              ! REDUCED MASS IN KG
      FCT=SD*(JH/(B2A*1.00D-10)**2)                               ! FORCE CONSTANT IN J.m**-2
      VIB=(1.00D+00/PI2)*(1.00D+00/C)*SQRT(FCT/REDMASS)           ! VIBRATIONAL FREQUENCY IN cm-1
      END SUBROUTINE DIATHARMFREQ

!###############################################################################
