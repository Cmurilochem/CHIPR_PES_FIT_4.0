!###############################################################################
! MODULE THAT CALCULATES PROPERTIES FOR TRIATOMICS
!###############################################################################
 
      SUBROUTINE PROP_TRIAT
      USE COMMON_VAR
      IMPLICIT NONE
      INTEGER :: T,V,U,NPC,SUMC,O,E
      DOUBLE PRECISION :: BAS,BASNEW,BASOLD
      DOUBLE PRECISION :: R,POT,STEPRD
      DOUBLE PRECISION :: RDMIN,RDMAX
      DOUBLE PRECISION :: ORIG,R0,TWOBD,THREBD
      DOUBLE PRECISION :: R0NEW,R0OLD
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: CINT,XINT,YINT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BS,BSNEW,BSOLD
      LOGICAL :: CHECK
      INTEGER, DIMENSION(NCPOL,NX) :: ID
      DOUBLE PRECISION, DIMENSION(NX) :: ARG
      DOUBLE PRECISION, DIMENSION(NP) :: P
      DOUBLE PRECISION, DIMENSION(13) :: THETA
      DOUBLE PRECISION, DIMENSION(18) :: R1
      DOUBLE PRECISION :: R2MIN,R2MAX,STEPR2
      DOUBLE PRECISION :: R1MIN,R1MAX,STEPR1
      DOUBLE PRECISION :: PHI,D1,D2,D3
      DOUBLE PRECISION, DIMENSION(20) :: Z1,Z2,Z3
      DOUBLE PRECISION, DIMENSION(NC) :: DC
      INTEGER, DIMENSION(NCPOL) :: SUBORDER      

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
!     WRITTING FITTED CURVE IN THE FILE PLOT_POT.RES 
!
        RDMIN=0.0D+00 
        RDMAX=15.00D+00
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

!        CALL SYSTEM ("cp SOURCE/PLOT_BASIS_OPT_CHIPR_TRIATOM.gnu .")
        CALL SYSTEM ("gnuplot PLOT_BASIS_OPT_CHIPR_TRIATOM.gnu")
!        CALL SYSTEM ("okular TRIAT_CHIPR.ps &")
!        CALL SYSTEM ("rm PLOT_BASIS_OPT_CHIPR_TRIATOM.gnu")

      ELSE IF (JOBTYP.EQ."FITPOL") THEN
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
        RDMIN=0.0D+00 
        RDMAX=15.00D+00
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
            WRITE(21,*) R0OLD,BASOLD,R0NEW,BASNEW
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
!
!     EVALUATING THE AB INITIO POINTS
!     PUT IN THE FILE PLOT_ABINITIO.RES 
!     THE FITTED POTENTIAL AT THE ABINITIO POINTS IS PUT IN THE FILE PLOT_POT.RES
!     THE ABSOLUTE DEVIATIONS ARE PUT IN THE FILE DEV_ABINITIO_POINTS.RES
!
        DO I=1,NP
          ARG(1)=X(I,1)
          ARG(2)=X(I,2)
          ARG(3)=X(I,3)
          CALL SUM2BD(ARG(1),ARG(2),ARG(3),TWOBD)
          CALL CHIPR_TRIAT(MBS,L,C,ARG,THREBD,DC)
          POT=(TWOBD+THREBD)
          WRITE(22,*) X(I,1),X(I,2),X(I,3),YYY(I)
          WRITE(23,*) X(I,1),X(I,2),X(I,3),POT          
          WRITE(24,*) X(I,1),X(I,2),X(I,3),I,(YYY(I)-POT)*EH2CM 
        END DO
!
!     WRITTING STRATIFIED RMSDs IN THE FILE STRAT_RMSD.RES
!
        DO I=1,NP
          ARG(1)=X(I,1)
          ARG(2)=X(I,2)
          ARG(3)=X(I,3)
          CALL SUM2BD(ARG(1),ARG(2),ARG(3),TWOBD)
          CALL CHIPR_TRIAT(MBS,L,C,ARG,THREBD,DC)
          P(I)=0.0D+00
          P(I)=(TWOBD+THREBD)
        END DO

        CALL STRAT_RMSD(P,YYY,W)
!
!     WRITTING FINAL COEFFICIENTS IN THE FILE COEFFS_TAB_LATEX.RES
!
            DO I=1,DEG
              IF (I.EQ.1) THEN
             
                K=NCPOL
                DO J=1,MBS(I)
                  K=K+1
                  WRITE(27,'("$c_{"I1,",",I3,"}$",6X,"&",X,D25.18,X,"\\")') I,J,C(K)
                END DO

                K=NCPOL
                DO J=1,MBS(I)
                  K=K+1
                  WRITE(27,'("$\gamma_{"I1,",",I3,"}$",X,"&",X,D25.18,X,"\\")') I,J,C(MBS(I)+K)
                END DO

                K=NCPOL
                WRITE(27,'("$\zeta$",10X,"&",X,D25.18,X,"\\")') C(K+2*MBS(I)+2)

                WRITE(27,'("$R^{\rm ref}$",4X,"&",X,D25.18,X,"\\")') C(K+2*MBS(I)+1)

                WRITE(27,*) 

              ELSE IF (I.EQ.2) THEN

                K=NCPOL+2*MBS(1)+2
                DO J=1,MBS(I)
                  K=K+1
                  WRITE(27,'("$c_{"I1,",",I3,"}$",6X,"&",X,D25.18,X,"\\")') I,J,C(K)
                END DO

                K=NCPOL+2*MBS(1)+2
                DO J=1,MBS(I)
                  K=K+1
                  WRITE(27,'("$\gamma_{"I1,",",I3,"}$",X,"&",X,D25.18,X,"\\")') I,J,C(MBS(I)+K)
                END DO

                K=NCPOL+2*MBS(1)+2
                WRITE(27,'("$\zeta$",10X,"&",X,D25.18,X,"\\")') C(K+2*MBS(I)+2)

                WRITE(27,'("$R^{\rm ref}$",4X,"&",X,D25.18,X,"\\")') C(K+2*MBS(I)+1)

                WRITE(27,*) 

              ELSE IF (I.EQ.3) THEN

                K=NCPOL+2*MBS(1)+2+2*MBS(2)+2
                DO J=1,MBS(I)
                  K=K+1
                  WRITE(27,'("$c_{"I1,",",I3,"}$",6X,"&",X,D25.18,X,"\\")') I,J,C(K)
                END DO

                K=NCPOL+2*MBS(1)+2+2*MBS(2)+2
                DO J=1,MBS(I)
                  K=K+1
                  WRITE(27,'("$\gamma_{"I1,",",I3,"}$",X,"&",X,D25.18,X,"\\")') I,J,C(MBS(I)+K)
                END DO

                K=NCPOL+2*MBS(1)+2+2*MBS(2)+2
                WRITE(27,'("$\zeta$",10X,"&",X,D25.18,X,"\\")') C(K+2*MBS(I)+2)

                WRITE(27,'("$R^{\rm ref}$",4X,"&",X,D25.18,X,"\\")') C(K+2*MBS(I)+1)

                WRITE(27,*) 

              END IF

            END DO

            CALL POLINDEXES(L,MOLTYP,NCPOL,SUBORDER,ID)

            DO I=1,NCPOL
              WRITE(27,'(I3,1X,"&",I3,1X,"&",I3,1X,"&",1X,"$C_{",I3,",",I3,",",I3,"}$",&
       &               1X,"&",I3,1X,"&",1X,I4,1X,"&",1X,D25.18,1X,"\\")') ID(I,1),ID(I,2),ID(I,3),&
       &               ID(I,1),ID(I,2),ID(I,3),SUBORDER(I),I,C(I)
            END DO
!
!     WRITTING FINAL ANALYTIC FUNCTION TO THE FILE "CHIPR_TRIAT_FUNC.f90"
!
        WRITE(26,'(6X,"MOLTYP=",1X,A4)') MOLTYP

        WRITE(26,'(6X,"DEG=",1X,I1)') DEG

        WRITE(26,'(6X,"NC=",1X,I5)') NC

        WRITE(26,'(6X,"NX=",1X,I1)') NX

        DO I=1,DEG
          WRITE(26,'(6X,"BSORDER(",I1,")=",1X,I2)') I,MBS(I)
        END DO

        WRITE(26,'(6X,"POLORDER=",1X,I3)') L

        DO I=1,NC
          WRITE(26,'(6X,"C(",I3,")=",1X,D42.35)') I,C(I)
        END DO

        CALL SYSTEM ("chmod +x ../CHIPR-4.0_SOURCE_CODE/SCRIPT_CHIPR_TRIAT.csh")
        CALL SYSTEM ("cp ../CHIPR-4.0_SOURCE_CODE/SCRIPT_CHIPR_TRIAT.csh .")
        CALL SYSTEM ("./SCRIPT_CHIPR_TRIAT.csh")
        CALL SYSTEM ("rm SCRIPT_CHIPR_TRIAT.csh")
!
!     THE USER NOW CAN EDIT AT WILL
!     AND PLOT WHAT (S)HE WANTS 	
!
!     IN MY CASE OF A3, MAKE VALENCE COORDINATES PLOTS
!
        R2MIN=0.500D+00 
        R2MAX=20.00D+00
        NPC=500

        STEPR2=(R2MAX-R2MIN)/(DBLE(NPC)-1.000D+00)

        THETA=(/60.0D+00,70.0D+00,80.0D+00,90.0D+00,100.0D+00,&
     &          110.0D+00,120.0D+00,130.0D+00,140.0D+00,&
     &          150.0D+00,160.0D+00,170.0D+00,180.0D+00/)

           R1=(/1.8D+00,1.9D+00,2.0D+00,2.1D+00,2.2D+00,&
     &          2.3D+00,2.4D+00,2.5D+00,2.6D+00,2.7D+00,&
     &          2.8D+00,2.9D+00,3.0D+00,3.1D+00,3.2D+00,&
     &          3.3D+00,3.4D+00,3.5D+00/)

        DO I=1,13

          DO J=1,18

            DO K=1,NPC

              ARG(1)=R1(J)

              ARG(2)=R2MIN+STEPR2*(DBLE(K)-1.000D+00)

              ARG(3)=SQRT(ARG(1)**2+ARG(2)**2&
     &                -2.0D+00*ARG(1)*ARG(2)*COS(THETA(I)*PI/DEGS))


              CALL SUM2BD(ARG(1),ARG(2),ARG(3),TWOBD)
              CALL CHIPR_TRIAT(MBS,L,C,ARG,THREBD,DC)
              POT=(TWOBD+THREBD)

              WRITE(I+100,*) ARG(1), ARG(2), ARG(3), POT

            END DO

              WRITE(I+100,*)
              WRITE(I+100,*)

          END DO

        END DO
!
!     ENERGY AT EQUILIBRIUM GEOMETRY
!
        DO I=1,1
          ARG(1)=2.62856036D+00
          ARG(2)=2.62856036D+00
          ARG(3)=2.62856036D+00
          CALL SUM2BD(ARG(1),ARG(2),ARG(3),TWOBD)
          CALL CHIPR_TRIAT(MBS,L,C,ARG,THREBD,DC)
          POT=(TWOBD+THREBD)
          WRITE(6,*) ARG(1),ARG(2),ARG(3),TWOBD,THREBD,POT
        END DO
!
!     EXTRA PLOTS FOR FIXED VALENCE ANGLES
!
        PHI=180.00D+00*PI/DEGS
        R1MIN=2.35D+00
        R2MIN=2.35D+00  
        R1MAX=2.6D+00
        R2MAX=2.6D+00
        NPC=11

        STEPR1=(R1MAX-R1MIN)/(DBLE(NPC)-1.0D+00)

        STEPR2=(R2MAX-R2MIN)/(DBLE(NPC)-1.0D+00) 

        DO I=1,NPC

          DO J=1,NPC

            D1=R1MIN+STEPR1*(DBLE(I)-1.000D+00)
            D2=R2MIN+STEPR2*(DBLE(J)-1.000D+00)
            D3=SQRT((D1**2)+(D2**2)-2.0D+00*D1*D2*COS(PHI))
 
            ARG(1)=D1
            ARG(2)=D2
            ARG(3)=D3

            CALL SUM2BD(ARG(1),ARG(2),ARG(3),TWOBD)
            CALL CHIPR_TRIAT(MBS,L,C,ARG,THREBD,DC)
            POT=(TWOBD+THREBD)

            WRITE(201,*) ARG(1), ARG(2), ARG(3), POT

          END DO

            WRITE(201,*)
            WRITE(201,*)

        END DO
!
!
!
        PHI=107.50D+00*PI/DEGS
        R1MIN=2.425D+00
        R2MIN=2.425D+00  
        R1MAX=2.65D+00
        R2MAX=2.65D+00
        NPC=11

        STEPR1=(R1MAX-R1MIN)/(DBLE(NPC)-1.0D+00)

        STEPR2=(R2MAX-R2MIN)/(DBLE(NPC)-1.0D+00) 

        DO I=1,NPC

          DO J=1,NPC

            D1=R1MIN+STEPR1*(DBLE(I)-1.000D+00)
            D2=R2MIN+STEPR2*(DBLE(J)-1.000D+00)
            D3=SQRT((D1**2)+(D2**2)-2.0D+00*D1*D2*COS(PHI))
 
            ARG(1)=D1
            ARG(2)=D2
            ARG(3)=D3

            CALL SUM2BD(ARG(1),ARG(2),ARG(3),TWOBD)
            CALL CHIPR_TRIAT(MBS,L,C,ARG,THREBD,DC)
            POT=(TWOBD+THREBD)

            WRITE(202,*) ARG(1), ARG(2), ARG(3), POT

          END DO

            WRITE(202,*)
            WRITE(202,*)

        END DO
!
!
!
        PHI=60.0D+00*PI/DEGS
        R1MIN=2.465D+00
        R2MIN=2.465D+00  
        R1MAX=2.725D+00
        R2MAX=2.725D+00
        NPC=11

        STEPR1=(R1MAX-R1MIN)/(DBLE(NPC)-1.0D+00)

        STEPR2=(R2MAX-R2MIN)/(DBLE(NPC)-1.0D+00) 

        DO I=1,NPC

          DO J=1,NPC

            D1=R1MIN+STEPR1*(DBLE(I)-1.000D+00)
            D2=R2MIN+STEPR2*(DBLE(J)-1.000D+00)
            D3=SQRT((D1**2)+(D2**2)-2.0D+00*D1*D2*COS(PHI))
 
            ARG(1)=D1
            ARG(2)=D2
            ARG(3)=D3

            CALL SUM2BD(ARG(1),ARG(2),ARG(3),TWOBD)
            CALL CHIPR_TRIAT(MBS,L,C,ARG,THREBD,DC)
            POT=(TWOBD+THREBD)

            WRITE(203,*) ARG(1), ARG(2), ARG(3), POT

          END DO

            WRITE(203,*)
            WRITE(203,*)

        END DO
!
!
!
        PHI=65.0D+00*PI/DEGS
        R1MIN=2.465D+00
        R2MIN=2.465D+00  
        R1MAX=2.725D+00
        R2MAX=2.725D+00
        NPC=11

        STEPR1=(R1MAX-R1MIN)/(DBLE(NPC)-1.0D+00)

        STEPR2=(R2MAX-R2MIN)/(DBLE(NPC)-1.0D+00) 

        DO I=1,NPC

          DO J=1,NPC

            D1=R1MIN+STEPR1*(DBLE(I)-1.000D+00)
            D2=R2MIN+STEPR2*(DBLE(J)-1.000D+00)
            D3=SQRT((D1**2)+(D2**2)-2.0D+00*D1*D2*COS(PHI))
 
            ARG(1)=D1
            ARG(2)=D2
            ARG(3)=D3

            CALL SUM2BD(ARG(1),ARG(2),ARG(3),TWOBD)
            CALL CHIPR_TRIAT(MBS,L,C,ARG,THREBD,DC)
            POT=(TWOBD+THREBD)

            WRITE(208,*) ARG(1), ARG(2), ARG(3), POT

          END DO

            WRITE(208,*)
            WRITE(208,*)

        END DO
!
!
!
        PHI=70.0D+00*PI/DEGS
        R1MIN=2.465D+00
        R2MIN=2.465D+00  
        R1MAX=2.725D+00
        R2MAX=2.725D+00
        NPC=11

        STEPR1=(R1MAX-R1MIN)/(DBLE(NPC)-1.0D+00)

        STEPR2=(R2MAX-R2MIN)/(DBLE(NPC)-1.0D+00) 

        DO I=1,NPC

          DO J=1,NPC

            D1=R1MIN+STEPR1*(DBLE(I)-1.000D+00)
            D2=R2MIN+STEPR2*(DBLE(J)-1.000D+00)
            D3=SQRT((D1**2)+(D2**2)-2.0D+00*D1*D2*COS(PHI))
 
            ARG(1)=D1
            ARG(2)=D2
            ARG(3)=D3

            CALL SUM2BD(ARG(1),ARG(2),ARG(3),TWOBD)
            CALL CHIPR_TRIAT(MBS,L,C,ARG,THREBD,DC)
            POT=(TWOBD+THREBD)

            WRITE(205,*) ARG(1), ARG(2), ARG(3), POT

          END DO

            WRITE(205,*)
            WRITE(205,*)

        END DO
!
!
!
        PHI=130.0D+00*PI/DEGS
        R1MIN=2.2D+00
        R2MIN=2.2D+00  
        R1MAX=2.6D+00
        R2MAX=2.6D+00
        NPC=11

        STEPR1=(R1MAX-R1MIN)/(DBLE(NPC)-1.0D+00)

        STEPR2=(R2MAX-R2MIN)/(DBLE(NPC)-1.0D+00) 

        DO I=1,NPC

          DO J=1,NPC

            D1=R1MIN+STEPR1*(DBLE(I)-1.000D+00)
            D2=R2MIN+STEPR2*(DBLE(J)-1.000D+00)
            D3=SQRT((D1**2)+(D2**2)-2.0D+00*D1*D2*COS(PHI))
 
            ARG(1)=D1
            ARG(2)=D2
            ARG(3)=D3

            CALL SUM2BD(ARG(1),ARG(2),ARG(3),TWOBD)
            CALL CHIPR_TRIAT(MBS,L,C,ARG,THREBD,DC)
            POT=(TWOBD+THREBD)

            WRITE(206,*) ARG(1), ARG(2), ARG(3), POT

          END DO

            WRITE(206,*)
            WRITE(206,*)

        END DO
!
!
!
        PHI=140.0D+00*PI/DEGS
        R1MIN=2.2D+00
        R2MIN=2.2D+00  
        R1MAX=2.6D+00
        R2MAX=2.6D+00
        NPC=11

        STEPR1=(R1MAX-R1MIN)/(DBLE(NPC)-1.0D+00)

        STEPR2=(R2MAX-R2MIN)/(DBLE(NPC)-1.0D+00) 

        DO I=1,NPC

          DO J=1,NPC

            D1=R1MIN+STEPR1*(DBLE(I)-1.000D+00)
            D2=R2MIN+STEPR2*(DBLE(J)-1.000D+00)
            D3=SQRT((D1**2)+(D2**2)-2.0D+00*D1*D2*COS(PHI))
 
            ARG(1)=D1
            ARG(2)=D2
            ARG(3)=D3

            CALL SUM2BD(ARG(1),ARG(2),ARG(3),TWOBD)
            CALL CHIPR_TRIAT(MBS,L,C,ARG,THREBD,DC)
            POT=(TWOBD+THREBD)

            WRITE(207,*) ARG(1), ARG(2), ARG(3), POT

          END DO

            WRITE(207,*)
            WRITE(207,*)

        END DO
!
!
!
        Z1=(/2.35000000,2.35775632,2.36551263,&
       &       2.37326895,2.38102526,2.38878158,&
       &       2.39653789,2.40429421,2.41205053,& 
       &       2.41980684,2.42756316,2.43531947,& 
       &       2.44307579,2.45083210,2.45858842,& 
       &       2.46634474,2.47410105,2.48185737,&
       &       2.48961368,2.49737000/)

        Z2=(/2.52702540,2.52514503,2.52326465,2.52142429,&
       &     2.51958392,2.51782357,2.51630326,&
       &     2.51530306,2.51602321,2.51650330,& 
       &     2.51534307,2.51362272,2.51170234,& 
       &     2.50970194,2.50766153,2.50562112,& 
       &     2.50354071,2.50150030,2.49941988,2.49737000/)

        DO I=1,20

          Z3(I)=Z2(I)+Z1(I)

          ARG(1)=Z1(I)
          ARG(2)=Z2(I)
          ARG(3)=Z3(I)

          CALL SUM2BD(ARG(1),ARG(2),ARG(3),TWOBD)
          CALL CHIPR_TRIAT(MBS,L,C,ARG,THREBD,DC)
          POT=(TWOBD+THREBD)

          WRITE(204,*) ARG(1), ARG(2), ARG(3), POT

        END DO
!
!       FINISH HERE YOUR EXTRA PLOTS
!
!        CALL SYSTEM ("cp SOURCE/PLOT_POL_OPT_CHIPR_TRIATOM.gnu .")
        CALL SYSTEM ("gnuplot PLOT_POL_OPT_CHIPR_TRIATOM.gnu")
!        CALL SYSTEM ("okular TRIAT_CHIPR.ps &")
!        CALL SYSTEM ("rm PLOT_POL_OPT_CHIPR_TRIATOM.gnu")

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

      RETURN
      END

!################################################################################
! SUBROUTINE TO CALCULATE THE EXCITATION INDEXES FOR A GIVEN POLYNOMIAL
! OF ORDER "ORDER" AND MOLECULE TYPE "MOLTYP"
!################################################################################ 

      SUBROUTINE POLINDEXES(ORDER,MOLTYP,NCOEFFS,SUBORDER,ID)
      IMPLICIT NONE
      INTEGER :: I,J,K,L,M,S
      INTEGER :: NC,MNUM,ORDER,NCOEFFS
      INTEGER, DIMENSION(NCOEFFS,3) :: ID
      CHARACTER(LEN=3) :: MOLTYP
      INTEGER, DIMENSION(NCOEFFS) :: SUBORDER

      NC=0
      S=0

      DO I=1,NCOEFFS
        DO J=1,3
          ID(I,J)=0
        END DO
      END DO

      DO M=0,ORDER
        MNUM=0
        DO I=0,M
          DO J=0,(M-I)
            K=M-I-J
            S=I+J+K
            IF (S.NE.I .AND. S.NE.J .AND. S.NE.K) THEN
!
!     FOR THE CASE OF AN ABC MOLECULE 
! 
              IF (MOLTYP.EQ."ABC") THEN  

                NC=NC+1
                MNUM=MNUM+1

                ID(NC,1)=I
                ID(NC,2)=J
                ID(NC,3)=K
                
                SUBORDER(NC)=(I+J+K)
!
!     FOR THE CASE OF AN AB2 MOLECULE 
!
              ELSE IF (MOLTYP.EQ."AB2") THEN 

                IF (J.LE.K) THEN 

                  NC=NC+1
                  MNUM=MNUM+1

                  ID(NC,1)=I
                  ID(NC,2)=J
                  ID(NC,3)=K
                  
                  SUBORDER(NC)=(I+J+K)

                END IF    
!
!     FOR THE CASE OF AN A3 MOLECULE 
!         
              ELSE IF (MOLTYP.EQ."A3") THEN

                IF (I.LE.J .AND. J.LE.K) THEN

                  NC=NC+1
                  MNUM=MNUM+1

                  ID(NC,1)=I
                  ID(NC,2)=J
                  ID(NC,3)=K
                  
                  SUBORDER(NC)=(I+J+K)

                END IF
!
!     END OF CASES
!
              END IF

            END IF
          END DO
        END DO
      END DO
!
!     CHECKING - JUST IN CASE
!
      IF (NC.NE.NCOEFFS) THEN 
        WRITE(6,*) "PROBLEMS IN DEFINING PROPER NUMBER OF COEFFICIENTS"
        WRITE(6,*) "TOTAL NUMBER OF COEFFS ARE NOT SUMMING UP CORRECTLY"
        STOP
      END IF

      RETURN
      END
