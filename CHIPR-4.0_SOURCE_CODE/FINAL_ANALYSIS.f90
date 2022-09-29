!###############################################################################
! PRINTING FINAL ATRIBUTES
!###############################################################################

      SUBROUTINE FINALANALYSIS
      USE COMMON_VAR
      IMPLICIT NONE
      DOUBLE PRECISION :: RMSDWW,RMSDNW
!
!     WRITE FINAL COEFFS INTO THE FILE OUT.RES
!
      OPEN(UNIT=11,FILE='out.res',IOSTAT=ERR,&
     & FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
        IF (ERR /= 0) THEN
        PRINT*,"PROBLEMS IN OPENING THE FILE OUT.RES"
        END IF

      WRITE(6,'("LMDIF INFO:")')
!
!     PRINTING LMDIF INFO
!
      CALL DIAGNOSTIC(INFO)
!
!     CALCULATING FINAL RMSDs
!
      CALL RMDS(NP,FVEC,W,RMSDWW,RMSDNW)

      WRITE(6,53) NP,NC,RMSDWW*EH2CM,&
     &            RMSDNW*EH2CM

      WRITE(11,53) NP,NC,RMSDWW*EH2CM,&
     &            RMSDNW*EH2CM

      WRITE(6,1001)
      WRITE(11,1001)
!
!     STARTING PLOTS FOR DIATOMICS, TRIATOMICS AND ETCS
!
      IF (NX.EQ.1 .OR. NX.EQ.3 .OR. NX.EQ.6) THEN

        IF (JOBTYP.EQ."OPTBASIS") THEN
!
!     THIS INDEX ALLOW FOR EVALUATION OF THE FUNCTION RATHER THAN FITTING
!
          EVALID=1

          CONTDEG=CONTDEG+1
          
          DO K=1,CONTDEG
            NC=NCBS(K)
            DO Q=1,NC
              WRITE(6,'(7X,"C(",I3,")=",1X,D42.35)') Q,CBS(Q,K)
              WRITE(11,'(7X,"C(",I3,")=",1X,D42.35)') Q,CBS(Q,K)
            END DO
           WRITE(6,*) 
           WRITE(11,*)
          END DO

        ELSE IF (JOBTYP.EQ."FITPOL" .OR. JOBTYP.EQ."DIRECTFIT") THEN
!
!     THIS INDEX ALLOW FOR EVALUATION OF THE FUNCTION RATHER THAN FITTING
!
          EVALID=1

          DO Q=1,NC
            WRITE(6,'(7X,"C(",I4,")=",1X,D42.35)') Q,C(Q)
            WRITE(11,'(7X,"C(",I4,")=",1X,D42.35)') Q,C(Q)
          END DO

        END IF
!
!  CALLING THE GLOBAL OPTIMIZATION PROGRAM AND FINAL ANALYSIS FOR EACH CASE
!
        IF (NX.EQ.1) THEN
        
          CALL GLOBALOPT
          
          CALL PROP_DIAT

        ELSE IF (NX.EQ.3) THEN
        
          CALL GLOBALOPT        

          CALL PROP_TRIAT

        ELSE IF (NX.EQ.6) THEN
        
          CALL GLOBALOPT      

          CALL PROP_TETRA

        END IF

      END IF 

      CLOSE(UNIT=11)

 53   FORMAT (/,"NUMBER OF POINTS:",&
     &         I10,/,"NUMBER OF PARAMETERS:",&
     &         I10,/,"FINAL RMSDW (IN CM-1):",&
     &         F15.5,/,"FINAL RMSD  (IN CM-1):",F15.5,/) 
 1001 FORMAT ("FINAL SET OF PARAMETERS:")

      RETURN 
      END

!###############################################################################
! THIS SUBROUTINE CALCULATES RMSDs
!###############################################################################

      SUBROUTINE RMDS(M,FVEC,W,RMSDWW,RMSDNW)
      IMPLICIT NONE
      INTEGER :: I,M
      DOUBLE PRECISION :: FVEC(M),W(M)
      DOUBLE PRECISION :: RMSDNW,RMSDWW      
      RMSDNW=0.0D0
      
      DO I=1,M
        RMSDNW=RMSDNW+(FVEC(I)/W(I))**2
      ENDDO
      
      RMSDWW=0.0D0
     
      DO I=1,M
        RMSDWW=RMSDWW+FVEC(I)**2
      ENDDO
      
      RMSDNW=SQRT(RMSDNW/DBLE(M-1.0D0))
      RMSDWW=SQRT(RMSDWW/DBLE(M-1.0D0))
      RETURN
      END

!###############################################################################
! PRINT INFO
!###############################################################################

      SUBROUTINE DIAGNOSTIC(INFO)
      INTEGER :: INFO
      CHARACTER*90 OUT(0:10)
      OUT(0)="IMPROPER INPUT PARAMETERS"
      OUT(1)="BOTH ACTUAL AND PREDICTED RELATIVE REDUCTIONS IN THE & 
     & SUM OF SQUARES ARE AT MOST FTOL."
      OUT(2)="RELATIVE ERROR BETWEEN TWO CONSECUTIVE ITERATES IS AT &
     & MOST XTOL."
      OUT(4)="THE COSINE OF THE ANGLE BETWEEN FVEC AND ANY COLUMN OF &
     & THE JACOBIAN IS AT MOST GTOL IN ABSOLUTE VALUE."
      OUT(5)="NUMBER OF CALLS TO FCN HAS REACHED OR EXCEEDED MAXFEC."
      OUT(6)="FTOL IS TOO SMALL. NO FURTHER REDUCTION IN THE SUM OF &
     & SQUARES IS POSSIBLE."
      OUT(7)="XTOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN THE &
     & APPROXIMATE SOLUTION X IS POSSIBLE."
      OUT(8)="GTOL IS TOO SMALL. FVEC IS ORTHOGONAL TO THE COLUMNS OF & 
     & THE JACOBIAN TO MACHINE PRECISION."
      IF(INFO.EQ.3)THEN
        WRITE(*,'(A)')OUT(1)
        WRITE(*,'("AND")')
        WRITE(*,'(A)')OUT(2)
        RETURN
      ELSE
        WRITE(*,'(A)')OUT(INFO)
        RETURN
      ENDIF
      END

!###############################################################################