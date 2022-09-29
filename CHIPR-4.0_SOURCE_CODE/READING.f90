!###############################################################################
! THIS SUBROUTINE READS RELEVANT DATA TO THE MAIN FITTING PROGRAM
!###############################################################################

      SUBROUTINE READING
      USE COMMON_VAR
      IMPLICIT NONE
      INTEGER :: DUM
      LOGICAL :: CHECK
      DOUBLE PRECISION :: WREEXP,WDEEXP,WWEEXP
      DOUBLE PRECISION :: WEIGHT,SUMV2,SUMV23
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: ARG
!
!     FILE CONTAING INPUT PARAMETERS 
!
      OPEN(UNIT=1,FILE='coeffs_guess.txt',IOSTAT=ERR,&
     & FORM='FORMATTED',STATUS='OLD',ACTION='READ')
        IF (ERR /= 0) THEN
        PRINT*,"PROBLEMS IN OPENING THE FILE COEFFS_GUESS.TXT"
        END IF        
!
!     READ NUMER OF ATOMS
!
      READ(1,*) NATOM
!
!     CHECKING THE READING PROCESS
!      
      IF (NATOM.GT.4) THEN
        WRITE(6,16)
        WRITE(6,*) "AT ITS CURRENT VERSION, THE PROGRAM DOES NOT GO BEYOND TETRATOMICS"
        WRITE(6,*) "PLEASE, CHECK THE NUMBER OF INPUT ATOMS"
        WRITE(6,17)
        STOP
      END IF
!
!     NX IS THE NUMBER OF DEGRESS OF FREEDOM FOR THE N ATOM SYSTEM
!     IF NX=1, THE MOLECULE IS A DIATOMIC
!     IF NX=3, THE MOLECULE IS A TRIATOMIC
!     IF NX=6, THE MOLECULE IS A TETRATOMIC... ETC
!
      NX=NATOM*(NATOM-1)/2
!
!     READ ATOM TYPE, ATOMIC NUMBER AND AMUs
!
      ALLOCATE(ATOM(NATOM),Z(NATOM),AMU(NATOM))

      DO I=1,NATOM
       READ(1,*,IOSTAT=IOS) ATOM(I), Z(I), AMU(I)
!
!     CHECKING THE READING PROCESS
!
       IF (IOS.GT.0 .OR. IOS.LT.0) THEN 
         WRITE(6,16)
         WRITE(6,*) "PLEASE, CHECK THE NUMBER AND TYPE OF ATOMS"
         WRITE(6,FMT='(1X,A49,2X,I2)') "IN TOTAL, THE NUMBER OF ATOMS AND TYPE SHOULD BE:",NATOM
         WRITE(6,17)
         STOP
       END IF
!       
      END DO
!
!     READ ZERO OF ENERGY
!
      READ(1,*) EZERO

      IF (NX.EQ.1 .OR. NX.EQ.3 .OR. NX.EQ.6) THEN
!
!     READING TYPE OF MOLECULE
!
        IF (NX.EQ.3 .OR. NX.EQ.6) THEN 

          READ(1,*) MOLTYP

        END IF

        IF (NX.EQ.6) THEN
!
!     MOLTYP=A4   - THE MOLECULE IS OF A4-TYPE
!     MOLTYP=AB3  - THE MOLECULE IS OF AB3-TYPE
!     MOLTYP=A2B2 - THE MOLECULE IS OF A2B2-TYPE
!     MOLTYP=ABC2 - THE MOLECULE IS OF ABC2-TYPE
!     MOLTYP=ABCD - THE MOLECULE IS OF ABCD-TYPE
!

!
!     SETTING VARIABLES
!
          IF (MOLTYP.EQ."A4") THEN 
            DEG=1
          ELSE IF (MOLTYP.EQ."AB3") THEN 
            DEG=2
          ELSE IF (MOLTYP.EQ."A2B2") THEN 
            DEG=3
          ELSE IF (MOLTYP.EQ."ABC2") THEN 
            DEG=4
          ELSE IF (MOLTYP.EQ."ABCD") THEN 
            DEG=6
          ELSE
            WRITE(6,*) "PLEASE, GIVE A PROPER VALUE FOR THE MOLTYP VARIABLE"
            STOP
          END IF

        ELSE IF (NX.EQ.3) THEN
!
!     MOLTYP=A3  - THE MOLECULE IS OF A3-TYPE
!     MOLTYP=AB2 - THE MOLECULE IS OF AB2-TYPE
!     MOLTYP=ABC - THE MOLECULE IS OF ABC-TYPE
!

!
!     SETTING VARIABLES
!
          IF (MOLTYP.EQ."A3") THEN 
            DEG=1
          ELSE IF (MOLTYP.EQ."AB2") THEN 
            DEG=2
          ELSE IF (MOLTYP.EQ."ABC") THEN 
            DEG=3
          ELSE
            WRITE(6,*) "PLEASE, GIVE A PROPER VALUE FOR THE MOLTYP VARIABLE"
            STOP
          END IF

        ELSE IF (NX.EQ.1) THEN
          DEG=1
        END IF
!
!     JOBTYP=OPTBASIS       - OPTIMIZE BASIS SET FOR EACH DEGREE OF FREEDOM ONLY - TAKEN ALWAYS AS THE FIRST STEP IN THE FIT
!     JOBTYP=FITPOL         - OPTIMIZE BASIS SET AND THEN FIT THE POLYNOMIAL     - TAKEN ALWAYS AS THE SECOND STEP IN THE FIT
!     JOBTYP=DIRECTFIT      - DIRECT FIT OF DIATOMICS TO SPECTROSCOPIC AND ABINITIO DATA
!
        READ(1,*) JOBTYP
!
!     CHECKING INPUT
!
        IF (JOBTYP.EQ."OPTBASIS") THEN
          CONTINUE
        ELSE IF (JOBTYP.EQ."FITPOL") THEN
          CONTINUE
        ELSE IF (JOBTYP.EQ."DIRECTFIT") THEN
          IF (NATOM.GT.2) THEN
            WRITE(6,16)
            WRITE(6,*) "AT ITS CURRENT VERSION, THE PROGRAM DOES NOT DO DIRECT FITS FOR MOLECULES BEYOND DIATOMICS"
            WRITE(6,17)
            STOP
          ELSE
            CONTINUE
          END IF          
        ELSE 
          WRITE(6,*) "PLEASE, GIVE A PROPER VALUE FOR THE JOBTYP VARIABLE"
          STOP
        END IF

        IF (JOBTYP.EQ."OPTBASIS") THEN
!
!     IN CASE OF OPTBASIS
!
!     OPTTYP=FIXNONLIN   - OPTIMIZE ONLY THE LINEAR COEFFICIENTS IN THE CONTRACTION, FIX NON-LINEAR PARAMETERS 
!     OPTTYP=OPTNONLIN   - OPTIMIZE GAMMA, IN ADDITION TO THE LINEAR PARAMETERS IN THE CONTRACTION, FIX ZETA AND RREF
!     OPTTYP=OPTZETARREF - OPTIMIZE ZETA AND RREF, IN ADDITION TO THE LINEAR PARAMETERS IN THE CONTRACTION, FIX GAMMA
!     OPTTYP=OPTALL      - OPTIMIZE GAMMA, ZETA AND RREF, IN ADDITION TO THE LINEAR PARAMETERS IN THE CONTRACTION
!
          READ(1,*) OPTTYP
!
!     CHECKING INPUT
!
          IF (OPTTYP.EQ."FIXNONLIN") THEN
            CONTINUE
          ELSE IF (OPTTYP.EQ."OPTNONLIN") THEN
            CONTINUE
          ELSE IF (OPTTYP.EQ."OPTZETARREF") THEN
            CONTINUE
          ELSE IF (OPTTYP.EQ."OPTALL") THEN
            CONTINUE
          ELSE 
            WRITE(6,*) "PLEASE, GIVE A PROPER VALUE FOR THE OPTTYP VARIABLE"
            STOP
          END IF

        ELSE IF (JOBTYP.EQ."FITPOL" .OR. JOBTYP.EQ."DIRECTFIT") THEN
!     
!     IN CASE OF FITPOL
!
!     OPTTYP=FIXBAS         - OPTIMIZE ONLY THE LINEAR COEFFICIENTS IN THE POLYNOMIAL, FIX BASIS FROM PREVIOUS OPTIMIZATION
!     OPTTYP=FIXBASNONLIN   - OPTIMIZE THE LINEAR COEFFICIENTS OF BOTH POLYNOMIAL AND BASIS SET, FIX NON-LINEAR PARAMETERS
!     OPTTYP=OPTBASNONLIN   - OPTIMIZE GAMMA, IN ADDITION TO THE LINEAR COEFFICIENTS OF BOTH POLYNOMIAL AND BASIS SET, FIX ZETA AND RREF
!     OPTTYP=OPTBASZETARREF - OPTIMIZE ZETA AND RREF, IN ADDITION TO THE LINEAR COEFFICIENTS OF BOTH POLYNOMIAL AND BASIS SET, FIX GAMMA
!     OPTTYP=OPTALL         - OPTIMIZE GAMMA, ZETA AND RREF, IN ADDITION TO THE LINEAR COEFFICIENTS OF BOTH POLYNOMIAL AND BASIS SET
!
          READ(1,*) OPTTYP
!
!     CHECKING INPUT
! 
          IF (OPTTYP.EQ."FIXBAS") THEN
            CONTINUE
          ELSE IF (OPTTYP.EQ."FIXBASNONLIN") THEN
            CONTINUE
          ELSE IF (OPTTYP.EQ."OPTBASNONLIN") THEN
            CONTINUE
          ELSE IF (OPTTYP.EQ."OPTBASZETARREF") THEN
            CONTINUE
          ELSE IF (OPTTYP.EQ."OPTALL") THEN
            CONTINUE
          ELSE 
            WRITE(6,*) "PLEASE, GIVE A PROPER VALUE FOR THE OPTTYP VARIABLE"
            STOP
          END IF

        END IF
!
!     IN CASE OF FITPOL OR DIRECTFIT
!
        IF (JOBTYP.EQ."FITPOL" .OR. JOBTYP.EQ."DIRECTFIT") THEN
!
! FILE CONTAING AB INITIO DATA
!
        OPEN(UNIT=2,FILE='abinitio_data.txt',IOSTAT=ERR,&
     &   FORM='FORMATTED',STATUS='OLD',ACTION='READ')
          IF (ERR /= 0) THEN
          PRINT*,"PROBLEMS IN OPENING THE FILE ABINITIO_DATA.TXT"
          END IF
!          
! FILE CONTAING A GUESS GEOMETRY FOR THE GLOBAL MINIMUM
!
        OPEN(UNIT=3,FILE='mol_gm_guess.txt',IOSTAT=ERR,&
     &   FORM='FORMATTED',STATUS='OLD',ACTION='READ')
          IF (ERR /= 0) THEN
          PRINT*,"PROBLEMS IN OPENING THE FILE MOL_GM_GUESS.TXT"
          END IF          
!
!     IF CHIPR IS USED THEN
!     READ L (THE ORDER OF THE POLYNOMIAL)
!     AND THE TOTAL NUMBER OF AB INITIO POINTS
!
          READ(1,*) L, NP
          
            IF (NX.EQ.6 .AND. L.LT.4) THEN
          
              WRITE(6,16)
              WRITE(6,*) "FOR TETRATOMICS, THE LOWEST ORDER OF THE POLYNOMIAL IS 4"
              WRITE(6,17)
              STOP
            
            ELSE IF (NX.EQ.3 .AND. L.LT.2) THEN
          
              WRITE(6,16)
              WRITE(6,*) "FOR TRIATOMICS, THE LOWEST ORDER OF THE POLYNOMIAL IS 2"
              WRITE(6,17)
              STOP
          
            ELSE
          
            END IF
!
! FOR THE CASE OF DIRECTFIT ONLY
!
            IF (JOBTYP.EQ."DIRECTFIT") THEN 
!          
! FILE CONTAING SPECTROSCOPIC DATA
!
             OPEN(UNIT=4,FILE='spec_data.txt',IOSTAT=ERR,&
     &        FORM='FORMATTED',STATUS='OLD',ACTION='READ')
               IF (ERR /= 0) THEN
               PRINT*,"PROBLEMS IN OPENING THE FILE SPEC_DATA.TXT"
               END IF
               
             OPEN(UNIT=2020,FILE='level.res',IOSTAT=ERR,&
     &        FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
               IF (ERR /= 0) THEN
               PRINT*,"PROBLEMS IN OPENING THE FILE LEVEL.RES"
               END IF
             
               NABINITIO=NP
! READING CHARGE AND ELECTRONIC CONTRIBUTION TO THE ANGULAR MOMENTUM OF THE DIATOMIC               
               READ(4,*) ICHARGE, IOMEGA
! READING EXPERIMENTAL EQUILIBRIUM GEOMETRY IN A.U. AND WEIGHT             
               READ(4,*) REEXP, WREEXP
               NP=NP+1
! READING EXPERIMENTAL DISSOCIATION ENERGY IN A.U. AND WEIGHT               
               READ(4,*) DEEXP, WDEEXP
               NP=NP+1
! READING EXPERIMENTAL HARMONIC VIBRATIONAL FREQUENCY IN CM-1 AND WEIGHT               
               READ(4,*) WEEXP, WWEEXP
               NP=NP+1
! READING THE NUMBER OF (RO)VIBRATIONAL LEVELS
               READ(4,*) NROVIB

               NP=NP+NROVIB
!
! CONSIDERING THAT THE FIRST DERIVATIVE IS ZERO AT RE
!
               EXPFD=0.0000D+00
!
! CALCULATING SECOND DERIVATIVE FROM EXP. HAR. FREQ.
!
               M1=AMU(1)
               M2=AMU(2)

               CALL FORCE(M1,M2,WEEXP,EXPSD)
               
               ALLOCATE(IVEXP(NROVIB),IJEXP(NROVIB))
               
               DO I=1,NROVIB
                 IVEXP(I)=0
                 IJEXP(I)=0
               END DO
        
           END IF
           
        END IF
!
!     ALOCATING REQUIRED VECTORS
!
        ALLOCATE(MBS(DEG),NPBS(DEG),NCBS(DEG))
!
!     INICIALIZING VECTORS - JUST IN CASE
!
        DO I=1,DEG
          MBS(I)=0
          NPBS(I)=0
          NCBS(I)=0
        END DO

!
!     READ MBS (THE LENGHT OF CONTRACTION)
!     AND THE NUMBER OF AB INITIO POINTS TO CALIBRATE THE BASIS SET FOR EACH DEGREE OF FREEDOM
!     DEFINE ALSO THE NUMBER OF COEFFS NCBS FOR EACH BASIS SET CONTRACTION
!
        DO I=1,DEG
          READ(1,*) MBS(I),NPBS(I)
          NCBS(I)=2*MBS(I)+2
          IF (NCBS(I).GT.NPBS(I)) THEN
            WRITE(6,'("# OF POINTS (",I5,") <, # OF COEFF",I5)') NPBS(I), NCBS(I)
            STOP
          ENDIF
        END DO
!
!     ALOCATING REQUIRED VECTORS
!
        ALLOCATE(CBS(MAXVAL(NCBS),DEG),&
     &           XBS(NX,MAXVAL(NPBS),DEG),&
     &           YYBS(MAXVAL(NPBS),DEG),&
     &           YBS(MAXVAL(NPBS),DEG),&
     &           WWBS(MAXVAL(NPBS),DEG),&
     &           WBS(MAXVAL(NPBS),DEG),&
     &           GAMAFINAL(MAXVAL(MBS),DEG),&
     &           RREFFINAL(1,DEG),&
     &           ZETAFINAL(1,DEG))
!
!     INICIALIZING VECTORS - JUST IN CASE
!
        DO I=1,DEG
          DO J=1,NCBS(I)
            CBS(J,I)=0.0D+00
          END DO
          DO K=1,NPBS(I)
            DO Q=1,NX
             XBS(Q,K,I)=0.0D+00
            END DO
            YBS(K,I)=0.0D+00
            WBS(K,I)=0.0D+00
          END DO
        END DO

        DO I=1,DEG
          DO J=1,MBS(I)
            GAMAFINAL(J,I)=0.0D+00 
          END DO
          RREFFINAL(1,I)=0.0D+00 
          ZETAFINAL(1,I)=0.0D+00
        END DO        
!
!     READ GUESS COEFFS FOR BASIS SET CONTRACTION AND RELEVANT AB INITIO DATA
!
        DO I=1,DEG
          DO J=1,NCBS(I)
            READ(1,*,IOSTAT=IOS) CBS(J,I)
!
!     CHECKING THE READING PROCESS
!
            IF (IOS.GT.0 .OR. IOS.LT.0) THEN 
              WRITE(6,16)
              WRITE(6,*) "PLEASE, CHECK THE NUMBER OF BASIS SET COEFFICIENTS"
              WRITE(6,FMT='(1X,A27,2X,I4)') "IN TOTAL, IT SHOULD AMOUNT:",NCBS(I)
              WRITE(6,FMT='(I4,2X,A33)') MBS(I), "CONTRACTION (LINEAR) COEFFICIENTS"
              WRITE(6,FMT='(I4,2X,A21)') MBS(I), "NON-LINEAR PARAMETERS"
              WRITE(6,*) "PLUS ONE REFERENCE DISTANCE AND ONE ZETA PARAMETER" 
              WRITE(6,17)
              STOP
            END IF
!            
          END DO
!
!     READING BASIS SET AB INITIO DATA
!
          DO K=1,NPBS(I)

            READ(1,*,IOSTAT=IOS) (XBS(Q,K,I),Q=1,NX),YYBS(K,I),WWBS(K,I)
!
!     SINCE THERE ARE RELATIVELLY FEW POINTS AT THIS STAGE, THE WEIGHTS ARE GIVEN BY "HAND"
!
            WBS(K,I)=WWBS(K,I)
!
!     CHECKING THE READING PROCESS
!
            IF (IOS.GT.0 .OR. IOS.LT.0) THEN
              WRITE(6,16)
              WRITE(6,*) "PLEASE, PROVIDE THE PROPER NUMBER OF AB INITIO POINTS FOR THE BASIS SET"
              WRITE(6,FMT='(1X,A22,2X,I5)') "THIS NUMBER SHOULD BE:",NPBS(I)
              WRITE(6,17)
              STOP
            END IF
!
            IF (NX.EQ.6) THEN  
!
!     CALCULATING THE FOUR-BODY ENERGIES
!
!     THE USER SHOULD PROVIDE THE SUBROUTINES THAT
!     CALCULATES THE SUM OF TWO- AND THREE-BODY ENERGIES
!     SEE SUBROUTINE SUM23BD IN "SUM23BD.f90"
!
              EVALID=2

              CALL SUM23BD(XBS(1,K,I),XBS(2,K,I),XBS(3,K,I),&
     &             XBS(4,K,I),XBS(5,K,I),XBS(6,K,I),SUMV23)
     
!
!     SUMV23 DEFINES THE SUM OF TWO- AND THREE-BODY ENERGIES
!
              YBS(K,I)=((YYBS(K,I)-EZERO)-(SUMV23))
                   
            ELSE IF (NX.EQ.3) THEN  
!
!     CALCULATING THE THREE-BODY ENERGIES
!
!     THE USER SHOULD PROVIDE THE SUBROUTINES FOR THE
!     RELEVANT DIATOMICS IN EACH COORDINATE
!     THESE ARE DEFINED BY THE FUNCTIONS DIAT1, DIAT2, AND DIAT3
!     IN THE SUBROUTINE SUM2BD IN "SUM2BD.f90" 
!
              EVALID=2

              CALL SUM2BD(XBS(1,K,I),XBS(2,K,I),XBS(3,K,I),SUMV2)
!
!     SUMV2 DEFINES THE SUM OF TWO-BODY ENERGIES
!
              YBS(K,I)=((YYBS(K,I)-EZERO)-(SUMV2))

            ELSE IF (NX.EQ.1) THEN 

              YBS(K,I)=(YYBS(K,I)-EZERO)

            END IF
            
          END DO
                   
        END DO
!
!     END OF OPTBASIS INPUT READING
!     
!
!     IN CASE OF FITPOL
!
        IF (JOBTYP.EQ."FITPOL" .OR. JOBTYP.EQ."DIRECTFIT") THEN

          IF (NX.EQ.1) THEN

            NCPOL=L

            NC=NCBS(1)+NCPOL
           
          ELSE IF (NX.EQ.3) THEN
!
!     CALCULATES THE NUMER OF COEFFS NC FOR A GIVEN POLYNOMIAL
!     OF ORDER L AND MOLECULE TYPE "MOLTYP"
!
            CALL POLNC(L,MOLTYP,NCPOL)

            NC=SUM(NCBS)+NCPOL
            
          ELSE IF (NX.EQ.6) THEN 
!
!     CALCULATES THE NUMER OF COEFFS NC FOR A GIVEN POLYNOMIAL
!     OF ORDER L AND MOLECULE TYPE "MOLTYP"
!          
            CALL POLNCTETRA(L,MOLTYP,NCPOL)
          
            NC=SUM(NCBS)+NCPOL

          END IF
                    
          ALLOCATE(CPOL(NCPOL),C(NC))

          DO I=1,NCPOL
            CPOL(I)=0.0D+00
          END DO

          DO I=1,NC
            C(I)=0.0D+00
          END DO
!
!     READ GUESS COEFFS FOR THE POLYNOMIAL
!
          DO I=1,NCPOL
            READ(1,*,IOSTAT=IOS) CPOL(I)
!
!     CHECKING THE READING PROCESS
!
            IF (IOS.GT.0 .OR. IOS.LT.0) THEN
              WRITE(6,16)
              WRITE(6,*) "PLEASE, PROVIDE THE PROPER NUMBER OF GUESS COEFFICIENTS FOR THE POLYNOMIAL"
              IF (NX.EQ.1) THEN
                MOL="DIATOMIC"
              ELSE IF (NX.EQ.3) THEN
                MOL="TRIATOMIC"
              ELSE IF (NX.EQ.6) THEN
                MOL="TETRATOMIC"
              ELSE
                STOP
              END IF
              WRITE(6,FMT='(1X,A18,6X,A10)') "FOR THIS MOLECULE:", MOL
              IF (NX.EQ.3 .OR. NX.EQ.6) THEN
                WRITE(6,FMT='(1X,A8,16X,A4)') "OF TYPE:",MOLTYP
              END IF
              WRITE(6,FMT='(1X,A48,2X,I4)') "THE POLYNOMIAL ORDER CHOSEN REQUIRES A TOTAL OF:",NCPOL
              WRITE(6,17)
              STOP
            END IF            
          END DO       
!
!     DEFINING TOTAL GUESS COEFFS INCLUDING THOSE IN THE POLYNOMIAL AND BASIS SETS
!
          DO I=1,NCPOL
            C(I)=CPOL(I)
          END DO

          K=0
          DO I=1,DEG
            DO J=1,NCBS(I)
              K=K+1
              C(NCPOL+K)=CBS(J,I)
            END DO
          END DO
!
!     ALOCATING REQUIRED VECTORS
!
          ALLOCATE(X(NP,NX),Y(NP),YY(NP),YYY(NP),ARG(NX),&
     &             W(NP),WW(NP),&
     &             FVEC(NP),DFVEC(NP,NC),&
     &             GAMAINITPOL(MAXVAL(MBS),DEG),&
     &             RREFINITPOL(1,DEG),&
     &             ZETAINITPOL(1,DEG),&
     &             CBSINITPOL(MAXVAL(NCBS),DEG))
!
!     INICIALIZING VECTORS
!
          DO I=1,NP
            DO J=1,NX
              X(I,J)=0.000D+00
            END DO
          END DO

          DO I=1,NP
            Y(I)=0.000D+00
            YY(I)=0.000D+00
            YYY(I)=0.000D+00
            W(I)=0.0000D+00 
            WW(I)=0.000D+00
            FVEC(I)=0.00000D+00
          END DO

          DO I=1,DEG
            DO J=1,MBS(I)
              GAMAINITPOL(J,I)=0.0D+00
            END DO
            RREFINITPOL(1,I)=0.0D+00 
            ZETAINITPOL(1,I)=0.0D+00
          END DO

          DO I=1,DEG
            DO J=1,NCBS(I)
              CBSINITPOL(J,I)=0.0D+00
            END DO
          END DO
!
!     KEEPING INITIAL VALUES OF THE NON-LINEAR PARAMETERS
! 
          K=0
          DO I=1,DEG
            DO J=1,MBS(I)
              K=MBS(I)
              GAMAINITPOL(J,I)=CBS(K+J,I)
            END DO
            DO J=1,NCBS(I)
              CBSINITPOL(J,I)=CBS(J,I)
            END DO
            RREFINITPOL(1,I)=CBS(2*MBS(I)+1,I)
            ZETAINITPOL(1,I)=CBS(2*MBS(I)+2,I)
          END DO
!
!     READ AB INITIO DATA
!
          IF (JOBTYP.EQ."DIRECTFIT") THEN 
            DUM=NABINITIO
          ELSE
            DUM=NP
          END IF
!
          DO I=1,DUM
            READ(2,*,IOSTAT=IOS) (X(I,J),J=1,NX), YY(I), WW(I)
!
!     CHECKING THE READING PROCESS
!            
            IF (IOS.GT.0 .OR. IOS.LT.0) THEN
              WRITE(6,16)
              WRITE(6,*) "PLEASE, PROVIDE THE PROPER NUMBER OF AB INITIO POINTS TO CALIBRATE THE POLYNOMIAL"
              WRITE(6,FMT='(1X,A22,2X,I5)') "THIS NUMBER SHOULD BE:",NP
              WRITE(6,17)
              STOP
            END IF
          END DO          
!
!     INITIALIZING VECTOR TO READ GLOBAL MINIMUM GUESS
!
          DO I=1,3*NATOMMAX
            CARTX(I)=0.0D+00
          END DO
!
!     READING THE CARTESIAN COORDINATES (IN A.U) FOR THE GLOBAL MINIMUM GUESS
!
          READ(3,*,IOSTAT=IOS) (CARTX(I),I=1,3*NATOM)
!
!     CHECKING THE READING PROCESS
!            
          IF (IOS.GT.0 .OR. IOS.LT.0) THEN
            WRITE(6,16)
            WRITE(6,*) "PLEASE, PROVIDE THE PROPER NUMBER OF CARTESIAN COORDINATES FOR THE GUESS GLOBAL MINIMUM"
            WRITE(6,FMT='(1X,A22,2X,I5)') "THIS NUMBER SHOULD BE:",(3*NATOM)
            WRITE(6,17)
            STOP
          END IF
!
!     DEFINING PROPER VALUES OF Y AND WEIGHTS
!
          IF (JOBTYP.EQ."DIRECTFIT") THEN 
            DUM=NABINITIO
          ELSE
            DUM=NP
          END IF
          
          DO I=1,DUM

            IF (NX.EQ.1) THEN

              Y(I)=(YY(I)-EZERO)

              IF (WW(I).NE.(1.0D+00)) THEN

                W(I)=WW(I)

              ELSE
!
!     W(I) MUST BE DEFINED USING THE INPUT WW(I)
!     USE THE FUNCTION "WEIGHT" IN WEIGHT_FUNC.f90
!     TO DEFINE A CONVENIENT FORM
!     THE SIMPLEST ONE IS W(I)=WW(I)
!
                ARG(1)=X(I,NX)
                W(I)=WEIGHT(ARG,Y(I),WW(I))

              END IF

            ELSE IF (NX.EQ.3) THEN
!
!     CALCULATING THE THREE-BODY ENERGIES
!
!     THE USER SHOULD PROVIDE THE SUBROUTINES FOR THE
!     RELEVANT DIATOMICS IN EACH COORDINATE
!     THESE ARE DEFINED BY THE FUNCTIONS DIAT1, DIAT2, AND DIAT3
!     IN THE SUBROUTINE SUM2BD IN "SOURCE/SUM2BD.f90" 
!
              EVALID=2

              CALL SUM2BD(X(I,1),X(I,2),X(I,3),SUMV2)

              YYY(I)=(YY(I)-EZERO)

              Y(I)=(YYY(I)-(SUMV2))

              IF (WW(I).NE.(1.0D+00)) THEN

                W(I)=WW(I)

              ELSE
!
!     W(I) MUST BE DEFINED USING THE INPUT WW(I)
!     USE THE FUNCTION "WEIGHT" IN WEIGHT_FUNC.f90
!     TO DEFINE A CONVENIENT FORM
!     THE SIMPLEST ONE IS W(I)=WW(I)
!
                ARG(1)=X(I,1)
                ARG(2)=X(I,2)
                ARG(3)=X(I,3)
                W(I)=WEIGHT(ARG,YYY(I),WW(I))
            
              END IF
              
            ELSE IF (NX.EQ.6) THEN              
!
!     CALCULATING THE FOUR-BODY ENERGIES
!
!     THE USER SHOULD PROVIDE THE SUBROUTINES THAT
!     CALCULATES THE SUM OF TWO- AND THREE-BODY ENERGIES
!     SEE SUBROUTINE SUM23BD IN "SOURCE/SUM23BD.f90"
!
              EVALID=2

              CALL SUM23BD(X(I,1),X(I,2),X(I,3),&
     &             X(I,4),X(I,5),X(I,6),SUMV23)
     
!
!     SUMV23 DEFINES THE SUM OF TWO- AND THREE-BODY ENERGIES
!
              YYY(I)=(YY(I)-EZERO)

              Y(I)=(YYY(I)-(SUMV23))  
             
              IF (WW(I).NE.(1.0D+00)) THEN

                W(I)=WW(I)

              ELSE
!
!     W(I) MUST BE DEFINED USING THE INPUT WW(I)
!     USE THE FUNCTION "WEIGHT" IN WEIGHT_FUNC.f90
!     TO DEFINE A CONVENIENT FORM
!     THE SIMPLEST ONE IS W(I)=WW(I)
!
                ARG(1)=X(I,1)
                ARG(2)=X(I,2)
                ARG(3)=X(I,3)
                ARG(4)=X(I,4)
                ARG(5)=X(I,5)
                ARG(6)=X(I,6)
                W(I)=WEIGHT(ARG,YYY(I),WW(I))
            
              END IF              
              
            END IF

          END DO        
!
!  FOR THE CASE OF JOBTYPE=DIRECTFIT
!
          IF (JOBTYP.EQ."DIRECTFIT") THEN
          
            DO I=1,NROVIB
            
              J=NABINITIO+I
              
              READ(4,*,IOSTAT=IOS) IVEXP(I), IJEXP(I), Y(J), W(J)
!
! CHECKING THE READING PROCESS
!            
              IF (IOS.GT.0 .OR. IOS.LT.0) THEN
                WRITE(6,16)
                WRITE(6,*) "PLEASE, PROVIDE THE PROPER NUMBER OF EXPERIMENTAL BAND ORIGINS"
                WRITE(6,FMT='(1X,A22,2X,I5)') "THIS NUMBER SHOULD BE:",(NROVIB)
                WRITE(6,17)
                STOP
              END IF
            
            END DO
!
!  IMPOSING THE CONDITION THAT THE FIRST DERIVATIVE IS ALWAYS ZERO AT RE
!          
            Y(NABINITIO+NROVIB+1)=EXPFD
            W(NABINITIO+NROVIB+1)=WREEXP
!
!  IMPOSING THE CONDITION THAT THE SECOND DERIVATIVE IS ALWAYS EQUAL TO THE EXPERIMENTAL VALUE
!
            Y(NABINITIO+NROVIB+2)=EXPSD
            W(NABINITIO+NROVIB+2)=WWEEXP
!
!  IMPOSING THE CONDITION THAT THE DE IS ALWAYS EQUAL TO THE EXPERIMENTAL DE
!            
            Y(NABINITIO+NROVIB+3)=DEEXP 
            W(NABINITIO+NROVIB+3)=WDEEXP           
                      
          END IF
                    
        END IF

      ELSE

      STOP

      END IF

      IF (NP.GT.NPMAX) THEN
        WRITE(6,'("INCREASE NPMAX: ",I5,">",I5)') NP, NPMAX
        STOP
      ELSE IF (NC.GT.NCMAX) THEN
        WRITE(6,'("INCREASE NCMAX: ",I5,">",I5)') NC, NCMAX
        STOP
      ELSE IF (NC.GT.NP) THEN
        WRITE(6,'("# OF POINTS (",I5,") <, # OF COEFF",I5)') NP, NC
        STOP
      ENDIF

      CLOSE(UNIT=1)
      INQUIRE(UNIT=2,OPENED=CHECK)
      IF (CHECK) THEN 
        CLOSE(UNIT=2)
      ELSE 
      END IF
      INQUIRE(UNIT=3,OPENED=CHECK)
      IF (CHECK) THEN 
        CLOSE(UNIT=3)
      ELSE 
      END IF
      INQUIRE(UNIT=4,OPENED=CHECK)
      IF (CHECK) THEN 
        CLOSE(UNIT=4)
      ELSE 
      END IF      
      
 16   FORMAT (80("#"),/,80("#"),//,1X,"ERROR",/) 
 17   FORMAT (80("#"),/,80("#"),/) 

      END SUBROUTINE READING
