!###############################################################################
! SUBROUTINE RESPONSIBLE FOR FITTING
!###############################################################################

      SUBROUTINE TOFIT(DEG,NP,NC,C,FVEC,DFVEC,IFLAG)
      USE COMMON_VAR, ONLY : NX,M,L,X,Y,VLIM
      USE COMMON_VAR, ONLY : JOBTYP,NXRED
      USE COMMON_VAR, ONLY : W,CONT,EH2CM
      USE COMMON_VAR, ONLY : MBS,MOLTYP,NCPOL
      USE COMMON_VAR, ONLY : NABINITIO,NROVIB
      USE COMMON_VAR, ONLY : IVIBMX,NPMAX,ITVIB
      USE COMMON_VAR, ONLY : AMU,Z,ICHARGE,IOMEGA
      USE COMMON_VAR, ONLY : IVEXP,IJEXP
      USE COMMON_VAR, ONLY : REEXP,EXPFD,EXPSD,DEEXP
      USE COMMON_VAR, ONLY : GAMAINIT,RREF0INIT,ZETAINIT
      USE COMMON_VAR, ONLY : CBSINITPOL,GAMAINITPOL
      USE COMMON_VAR, ONLY : RREFINITPOL,ZETAINITPOL
      USE COMMON_VAR, ONLY : OPTTYP,EVALID
      IMPLICIT NONE
      INTEGER :: I,J,NP,NC,IFLAG,DEG,DUM,O
      DOUBLE PRECISION, DIMENSION(NC) :: C,DIFFC 
      DOUBLE PRECISION, DIMENSION(NP) :: FVEC
      DOUBLE PRECISION, DIMENSION(NP,NC) :: DFVEC
      DOUBLE PRECISION, DIMENSION(NP) :: POT,YVAL
      DOUBLE PRECISION, DIMENSION(NX) :: ARG 
      DOUBLE PRECISION :: SUM,SUMW,FROOT,RE,DE,V,FE
      INTEGER :: IIV(IVIBMX),IIJ(IVIBMX)
      DOUBLE PRECISION :: ESOLN(IVIBMX),IROTMAX(IVIBMX)
      DOUBLE PRECISION :: ROTCTE(7,IVIBMX),ZPE
      DOUBLE PRECISION :: ALLROVIB(0:IVIBMX,0:IVIBMX)
      INTEGER :: IAN1,IMN1,IAN2,IMN2,NLEV1,NJM,JDJR,IVMAX
      DOUBLE PRECISION :: BAND,DER,DEOPT
      DOUBLE PRECISION :: SUMAB,SUMWAB
      DOUBLE PRECISION :: SUMRVIB,SUMWRVIB      
!
!     OBJECT FUNCTION FF:
!     NP IS THE NUMBER OF DATA POINTS
!     NC IS THE NUMBER OF VARIABLES TO WHICH THE FF OBJECT FUNCTION DEPENDS
!     C IS A VECTOR OF COEFFICIENTS
!     FVEC IS THE CHI SQUARED TO BE MINIMIZED
!
      IF (JOBTYP.EQ."DIRECTFIT") THEN 
        DUM=NABINITIO
      ELSE
        DUM=NP
      END IF
          
      SUM=0.0D+00
      SUMW=0.0D+00

      DO I=1,NC
        POT(I)=0.000D+00
        YVAL(I)=0.000D+00
      END DO
      
      DO I=1,NC
        DO J=1,NC
          DFVEC(I,J)=0.000D+00
        END DO
      END DO            
!
!     STARTING FVEC EVALUATION
!
      DO I=1,DUM

        IF (NX.EQ.1 .OR. NX.EQ.3 .OR. NX.EQ.6) THEN

          IF (JOBTYP.EQ."OPTBASIS") THEN
       
            CALL BASIS_CONTRACT(DEG,M,C,X(I,NXRED),YVAL(I))
            FVEC(I)=(Y(I)-YVAL(I))*W(I)

          ELSE IF (JOBTYP.EQ."FITPOL" .OR. JOBTYP.EQ."DIRECTFIT") THEN
          
            IF (NX.EQ.1) THEN

              ARG(1)=X(I,NX)

              CALL CHIPR_DIAT(MBS,L,C,ARG,POT(I),DIFFC)
              FVEC(I)=(Y(I)-POT(I))*W(I)
!
!  EVALUATING THE GRADIENT OF THE FVEC VECTOR ANALYTICALLY IN THE Cs COEFFS
!
              DO J=1,NCPOL
                DFVEC(I,J)=-W(I)*DIFFC(J)
              END DO

            ELSE IF (NX.EQ.3) THEN

              ARG(1)=X(I,1)
              ARG(2)=X(I,2)
              ARG(3)=X(I,3)

              CALL CHIPR_TRIAT(MBS,L,C,ARG,POT(I),DIFFC)
              FVEC(I)=(Y(I)-POT(I))*W(I)
!
!  EVALUATING THE GRADIENT OF THE FVEC VECTOR ANALYTICALLY IN THE Cs COEFFS
!
              DO J=1,NCPOL
                DFVEC(I,J)=-W(I)*DIFFC(J)
              END DO
              
            ELSE IF (NX.EQ.6) THEN 
            
              ARG(1)=X(I,1)
              ARG(2)=X(I,2)
              ARG(3)=X(I,3)
              ARG(4)=X(I,4)
              ARG(5)=X(I,5)
              ARG(6)=X(I,6)              
              
              CALL CHIPR_TETRA(MBS,L,C,ARG,POT(I),DIFFC)
              FVEC(I)=(Y(I)-POT(I))*W(I) 
!
!  EVALUATING THE GRADIENT OF THE FVEC VECTOR ANALYTICALLY IN THE Cs COEFFS
!
              DO J=1,NCPOL
                DFVEC(I,J)=-W(I)*DIFFC(J)
              END DO
                            
            END IF

          END IF

        END IF

      END DO
!
!  FOR THE CASE OF JOBTYPE=DIRECTFIT
!
      IF (JOBTYP.EQ."DIRECTFIT" .AND. NX.EQ.1) THEN
!
!  FINDING THE EQUILIBRIUM GEOMETRY AND DISSOCIATION ENERGY OF THE ACTUAL CURVE
!
        CALL DIRECTMINDEDIAT(RE,DE,C)
!
!  CALCULATING FIRST DERIVATIVE AT EXPERIMENTAL EQ. DISTANCE
!
        O=MBS(1)
 
        CALL DCHIPR_DIAT(O,L,C,REEXP,V,.TRUE.,DER)
        
        DEOPT=-V
!
!  CALCULATING FORCE CONSTANT AT EXPERIMENTAL EQ. DISTANCE
! 
        CALL DIRECTFC(REEXP,FE,C)
!
!  SETTING UP PARAMETERS FOR LEVEL
!        
        VLIM=DE*EH2CM

        IAN1=INT(Z(1))
        IMN1=INT(AMU(1)) 
        IAN2=INT(Z(2))
        IMN2=INT(AMU(2)) 
!
!  THIS CALCULATES ALL ROVIBRATIONAL LEVELS UP TO DISSOCIATION (VLIM)
!
        NLEV1=-9999 !-(NROVIB-1)!-9999 
        NJM=9999 !0!9999
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
     &             IOMEGA,VLIM,NLEV1,NJM,JDJR,&
     &             IIV,IIJ,ESOLN,ROTCTE,IVMAX,&
     &             IROTMAX,ALLROVIB,C)
!
!  DEFINING ZPE
!
        ZPE=ALLROVIB(0,IOMEGA)
!
!  FVEC EVALUATION FOR THE CALCULATED LEVELS
!
        DO I=1,NROVIB
        
          J=NABINITIO+I
          
          BAND=(ALLROVIB(IVEXP(I),IJEXP(I))-ZPE)
! FVEC IN A.U.          
          FVEC(J)=((Y(J)-BAND)/EH2CM)*W(J)
                  
        END DO
!
!  IMPOSING THE CONDITION THAT THE FIRST DERIVATIVE IS ALWAYS ZERO AT RE
!
        I=NABINITIO+NROVIB+1

        FVEC(I)=(Y(I)-DER)*W(I)
!
!  IMPOSING THE CONDITION THAT THE SECOND DERIVATIVE IS ALWAYS EQUAL TO THE EXPERIMENTAL VALUE
!        
        I=NABINITIO+NROVIB+2
        
        FVEC(I)=(Y(I)-FE)*W(I)
!
!  IMPOSING THE CONDITION THAT THE DE IS ALWAYS EQUAL TO THE EXPERIMENTAL DE
!        
        I=NABINITIO+NROVIB+3
        
        FVEC(I)=(Y(I)-DEOPT)*W(I)
                  
      END IF
!
!     CALCULATING TOTAL CHI**2
!
      DO I=1,NP
        SUMW=SUMW+FVEC(I)**2
        SUM=SUM+(FVEC(I)/W(I))**2
      ENDDO    
!
! SPLIT RMSDs BETWEEN DIFFERENT DATA SETS FOR JOBTYP.EQ."DIRECTFIT"
!      
      IF (JOBTYP.EQ."DIRECTFIT") THEN 
      
        SUMAB=0.0D+00
        SUMWAB=0.0D+00
        SUMRVIB=0.0D+00
        SUMWRVIB=0.0D+00
      
        DO I=1,NABINITIO
          SUMWAB=SUMWAB+FVEC(I)**2
          SUMAB=SUMAB+(FVEC(I)/W(I))**2
        ENDDO 
        
        I=0
        J=0
        
        DO J=1,NROVIB
          I=NABINITIO+J
          SUMWRVIB=SUMWRVIB+FVEC(I)**2
          SUMRVIB=SUMRVIB+(FVEC(I)/W(I))**2
        ENDDO                  
      
      END IF 
!
!     PRINTING RESULTS
!
      IF (IFLAG.EQ.0)THEN
      CONT=CONT+1

        IF(ABS(SUM).GT.1E6)THEN
          PRINT*,'STOP',SUM
        ENDIF

        IF (CONT.EQ.0) THEN 
        
          IF (JOBTYP.EQ."DIRECTFIT") THEN
          
            WRITE(6,40) CONT
            WRITE(6,52) SQRT(SUMW/DBLE(NP-1))*EH2CM,&
     &                  SQRT(SUM/DBLE(NP-1))*EH2CM
     
          ELSE
          
            WRITE(6,40) CONT
            WRITE(6,50) SQRT(SUMW/DBLE(NP-1))*EH2CM,&
     &                  SQRT(SUM/DBLE(NP-1))*EH2CM
      
          END IF
        
        ELSE
        
          IF (JOBTYP.EQ."DIRECTFIT") THEN        

            WRITE(6,40) CONT
            WRITE(6,53) SQRT(SUMW/DBLE(NP-1))*EH2CM,&
     &                  SQRT(SUM/DBLE(NP-1))*EH2CM
     
            WRITE(6,54) SQRT(SUMWAB/DBLE(NABINITIO-1))*EH2CM,&
     &                  SQRT(SUMAB/DBLE(NABINITIO-1))*EH2CM
     
            WRITE(6,55) SQRT(SUMWRVIB/DBLE(NROVIB-1))*EH2CM,&
     &                  SQRT(SUMRVIB/DBLE(NROVIB-1))*EH2CM           
     
          ELSE 
          
            WRITE(6,40) CONT
            WRITE(6,51) SQRT(SUMW/DBLE(NP-1))*EH2CM,&
     &                  SQRT(SUM/DBLE(NP-1))*EH2CM 
     
          END IF         
          
        END IF

      ENDIF 

 40   FORMAT (80("#"),//,"INTERATION NUMBER:",I4)

 50   FORMAT ("INITIAL RMSDW (IN CM-1):",&
     &         F15.5,/,"INITIAL RMSD  (IN CM-1):",F15.5,/)

 51   FORMAT ("RMSDW (IN CM-1):",&
     &         F15.5,/,"RMSD  (IN CM-1):",F15.5,/)
     
 52   FORMAT ("INITIAL RMSDW INCL. AB INITIO & SPECTROSC. DATA (IN CM-1):",&
     &         F15.5,/,"INITIAL RMSD  INCL. AB INITIO & SPECTROSC. DATA (IN CM-1):",F15.5,/)
     
 53   FORMAT ("RMSDW INCL. AB INITIO & SPECTROSC. DATA (IN CM-1):",&
     &         F15.5,/,"RMSD  INCL. AB INITIO & SPECTROSC. DATA (IN CM-1):",F15.5,/)   
     
 54   FORMAT ("RMSDW INCL. AB INITIO DATA ONLY (IN CM-1):",&
     &         F15.5,/,"RMSD  INCL. AB INITIO DATA ONLY (IN CM-1):",F15.5,/) 
     
 55   FORMAT ("RMSDW INCL. SPECTROSC. DATA ONLY (IN CM-1):",&
     &         F15.5,/,"RMSD  INCL. SPECTROSC. DATA ONLY (IN CM-1):",F15.5,/)            
     

      END SUBROUTINE TOFIT
