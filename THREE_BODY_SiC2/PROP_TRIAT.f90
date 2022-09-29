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
! FOR THIS SYSTEM
      DOUBLE PRECISION, DIMENSION(7) :: rJAC1
      DOUBLE PRECISION, DIMENSION(9) :: JACANG1
      DOUBLE PRECISION, DIMENSION(6) :: rJAC2
      DOUBLE PRECISION, DIMENSION(15) :: JACANG2
      DOUBLE PRECISION :: BRMIN,BRMAX,STEP,BR
      INTEGER :: NPO,TMP
        
      

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

        CALL SYSTEM ("gnuplot PLOT_BASIS_OPT_CHIPR_TRIATOM.gnu")

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
          WRITE(24,*) X(I,1),X(I,2),X(I,3),I,(YYY(I)-POT)*EH2CM,YYY(I) 
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
                  WRITE(27,'("$c_{",I1,",",I3,"}$",6X,"&",1X,D25.18,1X,"\\")') I,J,C(K)
                END DO

                K=NCPOL
                DO J=1,MBS(I)
                  K=K+1
                  WRITE(27,'("$\gamma_{",I1,",",I3,"}$",1X,"&",1X,D25.18,1X,"\\")') I,J,C(MBS(I)+K)
                END DO

                K=NCPOL
                WRITE(27,'("$\zeta$",10X,"&",1X,D25.18,1X,"\\")') C(K+2*MBS(I)+2)

                WRITE(27,'("$R^{\rm ref}$",4X,"&",1X,D25.18,1X,"\\")') C(K+2*MBS(I)+1)

                WRITE(27,*) 

              ELSE IF (I.EQ.2) THEN

                K=NCPOL+2*MBS(1)+2
                DO J=1,MBS(I)
                  K=K+1
                  WRITE(27,'("$c_{",I1,",",I3,"}$",6X,"&",1X,D25.18,1X,"\\")') I,J,C(K)
                END DO

                K=NCPOL+2*MBS(1)+2
                DO J=1,MBS(I)
                  K=K+1
                  WRITE(27,'("$\gamma_{",I1,",",I3,"}$",1X,"&",1X,D25.18,1X,"\\")') I,J,C(MBS(I)+K)
                END DO

                K=NCPOL+2*MBS(1)+2
                WRITE(27,'("$\zeta$",10X,"&",1X,D25.18,1X,"\\")') C(K+2*MBS(I)+2)

                WRITE(27,'("$R^{\rm ref}$",4X,"&",1X,D25.18,1X,"\\")') C(K+2*MBS(I)+1)

                WRITE(27,*) 

              ELSE IF (I.EQ.3) THEN

                K=NCPOL+2*MBS(1)+2+2*MBS(2)+2
                DO J=1,MBS(I)
                  K=K+1
                  WRITE(27,'("$c_{",I1,",",I3,"}$",6X,"&",1X,D25.18,1X,"\\")') I,J,C(K)
                END DO

                K=NCPOL+2*MBS(1)+2+2*MBS(2)+2
                DO J=1,MBS(I)
                  K=K+1
                  WRITE(27,'("$\gamma_{",I1,",",I3,"}$",1X,"&",1X,D25.18,1X,"\\")') I,J,C(MBS(I)+K)
                END DO

                K=NCPOL+2*MBS(1)+2+2*MBS(2)+2
                WRITE(27,'("$\zeta$",10X,"&",1X,D25.18,1X,"\\")') C(K+2*MBS(I)+2)

                WRITE(27,'("$R^{\rm ref}$",4X,"&",1X,D25.18,1X,"\\")') C(K+2*MBS(I)+1)

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
        BRMIN=0.5D+00
        BRMAX=20.0D+00
        NPO=500
        STEP=(BRMAX-BRMIN)/(DBLE(NPO)-1.0D+00)
 
!########## C2-----Si CHANNEL

        JACANG1=(/90.0D+00,75.0D+00,60.0D+00,45.0D+00,&
     &          30.0D+00,15.0D+00,0.0D+00,67.5D+00,82.5D+00/)
     
        rJAC1=(/2.00000000D+00,2.10000000D+00,&
     &          2.30000000D+00,2.47913152D+00,2.72913152D+00,&
     &          2.97913152D+00,3.22913152D+00/)
     
        DO I=1,9   !FOR EACH ANGLE
          DO J=1,7 !FOR EACH r
            DO K=1,NPO
              TMP=2000+I
              BR=BRMIN+STEP*(DBLE(K)-1.000D+00)
              
              D1=rJAC1(J)
              D2=SQRT((0.5D+00*rJAC1(J))**2+BR**2&
     &           -2.0D+00*(0.5D+00*rJAC1(J))*BR*COS(JACANG1(I)*PI/DEGS))
              D3=SQRT((0.5D+00*rJAC1(J))**2+BR**2&
     &           +2.0D+00*(0.5D+00*rJAC1(J))*BR*COS(JACANG1(I)*PI/DEGS))
     
              ARG(1)=D1
              ARG(2)=D2
              ARG(3)=D3  
              
              CALL SUM2BD(ARG(1),ARG(2),ARG(3),TWOBD)
              CALL CHIPR_TRIAT(MBS,L,C,ARG,THREBD,DC)
              POT=(TWOBD+THREBD)

              WRITE(TMP,*) rJAC1(J),BR,JACANG1(I),POT,D1,D2,D3
              
            END DO
            WRITE(TMP,*)
            WRITE(TMP,*)
          END DO
        END DO

!########## SiC----C CHANNEL
 
        JACANG2=(/180.0D+00,165.0D+00,150.0D+00,135.0D+00,&
     &         120.0D+00,105.0D+00,90.0D+00,75.0D+00,&
     &         60.0D+00,45.0D+00,30.0D+00,15.0D+00,&
     &         0.0D+00,112.5D+00,127.5D+00/) 
     
        rJAC2=(/2.80000000D+00,2.96000000D+00,&
     &          3.26000000D+00,3.56000000D+00,&
     &          3.86000000D+00,4.15000000D+00/) 
     
        DO I=1,15   !FOR EACH ANGLE
          DO J=1,6  !FOR EACH r
            DO K=1,NPO
              TMP=3000+I
              BR=BRMIN+STEP*(DBLE(K)-1.000D+00)
              
              D2=rJAC2(J)
              D1=SQRT((0.5D+00*rJAC2(J))**2+BR**2&
     &           -2.0D+00*(0.5D+00*rJAC2(J))*BR*COS(JACANG2(I)*PI/DEGS))
              D3=SQRT((0.5D+00*rJAC2(J))**2+BR**2&
     &           +2.0D+00*(0.5D+00*rJAC2(J))*BR*COS(JACANG2(I)*PI/DEGS))
     
              ARG(1)=D1
              ARG(2)=D2
              ARG(3)=D3  
              
              CALL SUM2BD(ARG(1),ARG(2),ARG(3),TWOBD)
              CALL CHIPR_TRIAT(MBS,L,C,ARG,THREBD,DC)
              POT=(TWOBD+THREBD)

              WRITE(TMP,*) rJAC2(J),BR,JACANG2(I),POT,D1,D2,D3
              
            END DO
            WRITE(TMP,*)
            WRITE(TMP,*)
          END DO
        END DO     
 
!########## 
!
!       FINISH HERE YOUR EXTRA PLOTS
!
        CALL SYSTEM ("gnuplot PLOT_POL_OPT_CHIPR_TRIATOM.gnu")

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
