!###############################################################################
!
!    Copyright 2020 C. M. R. Rocha and A. J. C. Varandas
!
!    CHIPR-4.0 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.
!
!    CHIPR-4.0 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with CHIPR-4.0.  If not, see <https://www.gnu.org/licenses/>.
!
!###############################################################################
!
!    MAIN PROGRAM FOR THE CHIPR FITTING SUBROUTINE
!    
!    CHIPR-4.0 IS DESIGNED FOR FITTING GLOBAL POTENTIAL ENERGY  
!    SURFACES OF DIATOMIC, TRIATOMIC AND TETRATOMIC MOLECULES OF ANY TYPE USING 
!    THE CHIPR METHOD [J. CHEM. PHYS. 138, 054120 (2013); DOI: 10.1063/1.4788912; 
!    J. CHEM. PHYS. 138, 134117 (2013); DOI: 10.1063/1.4795826; 
!    J. PHYS. CHEM. A 123, 8154 (2019); DOI: 10.1021/acs.jpca.9b03194;
!    PHYS. CHEM. CHEM. PHYS.  21, 24406 (2019); DOI: 10.1039/c9cp04890a]. 
!    
!    WHEN USING THIS PROGRAM, THE USER MUST CITE:  
!
!    1). "A General Code for Fitting Global Potential Energy Surfaces via CHIPR Method: Triatomic Molecules" 
!    C. M. R. Rocha and A. J. C. Varandas
!    COMPUT. PHYS. COMMUN. 247, 106913 (2020) 
!
!    2). "A General Code for Fitting Global Potential Energy Surfaces via CHIPR Method: Direct-Fit Diatomics and Tetratomic Molecules" 
!    C. M. R. Rocha and A. J. C. Varandas
!    COMPUT. PHYS. COMMUN. 
!    to be submitted (2020)
!
!###############################################################################

      PROGRAM CHIPR_FIT
      USE COMMON_VAR
      IMPLICIT NONE
      INTEGER :: IFLAG,SUMC

      EXTERNAL TOFIT
!
!     REMOVING OLD FILES
!
      CALL SYSTEM ("rm *.res 2> /dev/null")
!
      CALL CPU_TIME(START)
      
      WRITE(6,5)
!
!     READING RELEVANT DATA
!
      CALL READING
!
!     MODULE TO FIT DIATOMICS AND TRIATOMICS
!
      IF (NX.EQ.1 .OR. NX.EQ.3 .OR. NX.EQ.6) THEN 

        IF (NX.EQ.1) THEN 

          WRITE(6,8)

          IF (JOBTYP.EQ."DIRECTFIT") THEN

            WRITE(6,15)

          END IF
        
        ELSE IF (NX.EQ.3) THEN 

          WRITE(6,7) MOLTYP

        ELSE IF (NX.EQ.6) THEN 

          WRITE(6,6) MOLTYP
 
        END IF

        IF (JOBTYP.EQ."OPTBASIS") THEN

          CONTDEG=0

          DO I=1,DEG

            NP=NPBS(I)
            NC=NCBS(I)
            M=MBS(I)

            WRITE(6,10) I
            WRITE(6,51)
            WRITE(6,52) NP,NC
         
            IF (NX.EQ.1) THEN

              NXRED=NX

            ELSE IF (NX.EQ.3) THEN 
!
!     REDUCING THE DIMENSIONALITY OF THE PROBLEM
!     IN ORDER TO FIT THE 1D BASIS SETS
!
              NXRED=NX-2

            ELSE IF (NX.EQ.6) THEN 
!
!     REDUCING THE DIMENSIONALITY OF THE PROBLEM
!     IN ORDER TO FIT THE 1D BASIS SETS
!
              NXRED=NX-5

            END IF
!
!     ALOCATING REQUIRED VECTORS
!
            ALLOCATE(C(NC),X(NP,NXRED),&
     &               Y(NP),W(NP),GAMAINIT(M),&
     &               FVEC(NP),DFVEC(NP,NC))
!
!     INICIALIZING VECTORS - JUST IN CASE
!
            C=0.000D+00
            X=0.000D+00
            Y=0.000D+00     
            W=0.0000D+00 
            FVEC=0.0000D+00  
            
            !DO J=1,NC
            !  C(J)=0.000D+00
            !END DO            

            !DO J=1,NP
            !  DO K=1,NXRED
            !    X(J,K)=0.000D+00
            !  END DO
            !END DO

            !DO J=1,NP
            ! Y(J)=0.000D+00          
            ! W(J)=0.0000D+00 
            ! FVEC(J)=0.00000D+00
            !END DO
!
!     PASSING ALL THE INFORMATION REGARDING TO THE 1D CUTS TO RELEVANT VECTORS
!
            DO J=1,NP
              DO K=1,NXRED
!
!     NOTE THAT THE COORDINATES TO BE FITTED HERE CORRESPONDS TO THE SAME VALUES OF DEG
!
!
!     FOR THE CASE OF A TETRATOMIC, PUT CONDITIONS IN WHAT COLUMNS CORRESPOND THE 
!     COORDINATES TO BE FITTED
!
                IF (NX.EQ.6) THEN
                  
                  IF (MOLTYP.EQ."A4") THEN
! R1 CORRESPONDS TO A-A DISTANCES AND WILL CALIBRATE Y1-Y6 BASIS
                    X(J,K)=XBS(I,J,I)
                  ELSE IF (MOLTYP.EQ."AB3") THEN
! R1 CORRESPONDS TO A-B DISTANCES AND WILL CALIBRATE Y1-Y3 BASIS
! R4 CORRESPONDS TO B-B DISTANCES AND WILL CALIBRATE Y4-Y6 BASIS
                    X(J,K)=XBS(((I-1)*I+I),J,I)
                  ELSE IF (MOLTYP.EQ."A2B2") THEN
! R1 CORRESPONDS TO A-A DISTANCES AND WILL CALIBRATE Y1 BASIS
! R3 CORRESPONDS TO A-B DISTANCES AND WILL CALIBRATE Y2-Y5 BASIS
! R6 CORRESPONDS TO B-B DISTANCES AND WILL CALIBRATE Y6 BASIS
                    X(J,K)=XBS(((I*(I+1))/2),J,I)
                  ELSE IF (MOLTYP.EQ."ABC2") THEN
! R1 CORRESPONDS TO C-C DISTANCES AND WILL CALIBRATE Y1 BASIS
! R2 CORRESPONDS TO C-A DISTANCES AND WILL CALIBRATE Y2 AND Y4 BASIS
! R3 CORRESPONDS TO C-B DISTANCES AND WILL CALIBRATE Y3 AND Y5 BASIS
! R6 CORRESPONDS TO A-B DISTANCES AND WILL CALIBRATE Y6 BASIS
                    IF (I.EQ.4) THEN
                      X(J,K)=XBS(6,J,I)
                    ELSE
                      X(J,K)=XBS(I,J,I)
                    END IF
!                      PRINT*,X(J,K)
                  ELSE IF (MOLTYP.EQ."ABCD") THEN
! EACH Ri WILL CALIBRATE Yi BASIS AND WILL BE READ IN SEQUENCE
                    X(J,K)=XBS(I,J,I)
                  ELSE 
                    STOP
                  END IF

                ELSE 
! FOR A3 MOLECULES
! R1 CORRESPONDS TO A-A DISTANCES AND WILL CALIBRATE Y1-Y3 BASIS
!
!                
! FOR AB2 MOLECULES
! R1 CORRESPONDS TO B-B DISTANCES AND WILL CALIBRATE Y1 BASIS
! R2 CORRESPONDS TO A-B DISTANCES AND WILL CALIBRATE Y2-Y3 BASIS
!
! 
! FOR ABC MOLECULES
! EACH Ri WILL CALIBRATE Yi BASIS AND WILL BE READ IN SEQUENCE
!
                  X(J,K)=XBS(I,J,I)

                END IF
                  
              END DO
            END DO

            DO J=1,NP
              Y(J)=YBS(J,I)         
              W(J)=WBS(J,I)
            END DO

            DO J=1,NC
              C(J)=CBS(J,I)      
            END DO
!
!     KEEPING INITIAL VALUES FOR THE NON-LINEAR PARAMETERS
! 
            DO J=1,M
              GAMAINIT(J)=C(M+J)
            END DO

            RREF0INIT=C(2*M+1)
            ZETAINIT=C(2*M+2)
!
!     EVALUATING THE RELIABLITY OF THE GUESS COEFFICIENTS
!            
            CONT=-1

            IFLAG=0

            EVALID=0

            CALL TOFIT(I,NP,NC,C,FVEC,DFVEC,IFLAG)            
!
!     CALLING MINPACK
!
            INFO=0

            CALL LSFIT(NP,NC,C,FVEC,DFVEC)

            IF (INFO.GE.1.AND.INFO.LE.3) THEN

              WRITE(6,100)
!
!      KEEPING THE OPTIMIZED VALUES AT THE SEPARATE VECTORS
!
              DO J=1,NC
                CBS(J,I)=C(J)
              END DO

              DO J=1,M
                GAMAFINAL(J,I)=C(M+J)
              END DO

              RREFFINAL(1,I)=C(2*M+1)
              ZETAFINAL(1,I)=C(2*M+2)
 
              CALL FINALANALYSIS
  
            ELSE
              WRITE(6,200)
            ENDIF

            DEALLOCATE(C,X,Y,W,GAMAINIT,FVEC,DFVEC)

          END DO

        ELSE IF (JOBTYP.EQ."FITPOL" .OR. JOBTYP.EQ."DIRECTFIT") THEN

          IF (OPTTYP.EQ."FIXBAS") THEN 

            WRITE(6,12)
            WRITE(6,51)
            WRITE(6,52) NP,(NC-SUM(NCBS))
   
          ELSE 

            WRITE(6,13)
            WRITE(6,51)
            IF (JOBTYP.EQ."FITPOL") THEN
              WRITE(6,52) NP,NC
            ELSE IF (JOBTYP.EQ."DIRECTFIT") THEN 
              WRITE(6,53) NABINITIO,(NROVIB+3),NC
            END IF

          END IF
!
!     EVALUATING THE RELIABLITY OF THE GUESS COEFFICIENTS
!
          CONT=-1

          IFLAG=0

          EVALID=0

          CALL TOFIT(DEG,NP,NC,C,FVEC,DFVEC,IFLAG)
          
          !STOP
!
!     CALLING MINPACK
!
          INFO=0
          
          CALL LSFIT(NP,NC,C,FVEC,DFVEC)
          
          IF (INFO.GE.1.AND.INFO.LE.3) THEN

            WRITE(6,100)
!
!      KEEPING THE OPTIMIZED VALUES AT THE SEPARATE VECTORS
!
            SUMC=NCPOL
            DO I=1,DEG
              K=SUMC
              DO J=1,NCBS(I)
                K=K+1
                CBS(J,I)=C(K)
              END DO

              DO J=1,MBS(I)
                GAMAFINAL(J,I)=CBS(MBS(I)+J,I)
              END DO

              RREFFINAL(1,I)=CBS(2*MBS(I)+1,I)
              ZETAFINAL(1,I)=CBS(2*MBS(I)+2,I)

              SUMC=SUMC+NCBS(I)
            END DO         

            CALL FINALANALYSIS

          ELSE
            WRITE(6,200)
          ENDIF

        END IF

      END IF

      CALL CPU_TIME(FINAL)
      
      WRITE(6,210) 
      WRITE(6,*) "CPU TIME (IN SECONDS):", (FINAL-START)
      WRITE(6,210)
      
      STOP

  5   FORMAT (80("#"),/,80("#"),//,6X,"CHIPR-4.0: A GENERAL CODE FOR FITTING GLOBAL PESs VIA CHIPR METHOD",//,&
 &            20X,"C. M. R. ROCHA AND A. J. C. VARANDAS",/,&
 &            19X,"COMPUT. PHYS. COMMUN. 258, 107556 (2021)"/)
  6   FORMAT (80("#"),/,80("#"),//,13X,"CHIPR PES FOR A TETRATOMIC MOLECULE OF TYPE -->",1X,A5,/)
  7   FORMAT (80("#"),/,80("#"),//,13X,"CHIPR PES FOR A TRIATOMIC MOLECULE OF TYPE -->",1X,A3,/)
  8   FORMAT (80("#"),/,80("#"),//,20X,"CHIPR CURVE FOR A DIATOMIC MOLECULE",/)
 15   FORMAT (80("#"),/,80("#"),//,15X,"DIRECT-FIT TO AB INITIO AND EXPERIMENTAL DATA",/)
 10   FORMAT (80("#"),/,80("#"),//,20X,"FITTING BASIS SET FOR COORDINATE",1X,"Y",I1,/)
 12   FORMAT (80("#"),/,80("#"),//,9X,"FITTING POLYNOMIAL WITH OPTIMUM BASIS SETS FROM PREVIOUS RUN",/)
 13   FORMAT (80("#"),/,80("#"),//,16X,"FITTING POLYNOMIAL AND OPTIMIZING BASIS SETS",/)
 51   FORMAT (80("#"),/,80("#"),//,"INPUT DATA")
 52   FORMAT ("NUMBER OF POINTS: ",I5,/,"NUMBER OF COEFFICIENTS: ",I5,/)
 53   FORMAT ("NUMBER OF AB INITIO POINTS: ",I5,/,"NUMBER OF EXPERIMENTAL DATA:",I5,/,"NUMBER OF COEFFICIENTS: ",I5,/)
 100  FORMAT (80("#"),/,80("#"),//,32X,"FIT CONVERGED",//,80("#"),/,80("#"),/)
 200  FORMAT (80("#"),/,80("#"),//,30X,"FIT NOT CONVERGED",//,80("#"),/,80("#"),/)
 210  FORMAT (/,80("#"),/,80("#"),/) 

      END PROGRAM CHIPR_FIT
