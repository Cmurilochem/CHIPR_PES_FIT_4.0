      PROGRAM TEST
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IVIBMX=400,NPMAX=100,ITVIB=10000)
      DIMENSION :: IIV(IVIBMX),IIJ(IVIBMX)
      DIMENSION :: ESOLN(IVIBMX),IROTMAX(IVIBMX)
      DIMENSION :: ROTCTE(7,IVIBMX)
      DIMENSION :: ALLROVIB(0:IVIBMX,0:IVIBMX)  
      DIMENSION :: X(NPMAX),Y(NPMAX),W(NPMAX)
      DIMENSION :: C(NPMAX),FVEC(NPMAX),V(NPMAX)
      DIMENSION :: SORTROVIB(ITVIB)
      CHARACTER :: ID

      EXTERNAL CLASSTURN
      COMMON/DATARC/X,Y,W,VLIM

      OPEN(UNIT=100,FILE='level-C2-order.res')    

      AU2CM=219474.6313702D+00
      
      IAN1=6 
      IMN1=12 
      IAN2=6
      IMN2=12 
      ICHARGE=0 

      IOMEGA=0
      VLIM=0.23234*AU2CM  

      NLEV1=-19!-9999
      NJM=0!9999	
      JDJR=1
C      
C     THESE SHOULD DEFINE THE LOWEST STATE OF THE MOLECULE
C
      IIV(1)=0
      IIJ(1)=IOMEGA
       
      CALL LEVEL(IAN1,IMN1,IAN2,IMN2,ICHARGE,
     *           IOMEGA,VLIM,NLEV1,NJM,JDJR,
     *           IIV,IIJ,ESOLN,ROTCTE,IVMAX,
     *           IROTMAX,ALLROVIB)
      
      WRITE(*,*) IVMAX,ALLROVIB(IIV(1),IIJ(1))
      WRITE(*,*) (IROTMAX(I),I=1,IVMAX)
        
      ICONT=0
      DO I=1,IVMAX
        IV=I-1
          IF (IROTMAX(I).NE.0) THEN
            DO K=1,IROTMAX(I)
              J=IIJ(1)+(K-1)*JDJR
              ICONT=ICONT+1
              SORTROVIB(ICONT)=ALLROVIB(IV,J)
              WRITE(*,*) IV,J,ALLROVIB(IV,J),ICONT,SORTROVIB(ICONT)
            END DO
          ELSE IF (IROTMAX(I).EQ.0) THEN
            WRITE(*,*) IV,ALLROVIB(IV,IIJ(1)),ROTCTE(1,I)
          END IF  
      END DO
 
      ID= 'I'

      CALL DLASRT(ID,ICONT,SORTROVIB,INFO)

      WRITE(100,FMT='(I12)') ICONT

      DO L=1,ICONT
        DO I=1,IVMAX
          IV=I-1
            IF (IROTMAX(I).NE.0) THEN
              DO K=1,IROTMAX(I)
                J=IIJ(1)+(K-1)*JDJR
                  IF (SORTROVIB(L).EQ.ALLROVIB(IV,J)) THEN
                    WRITE(100,10) IV,J,ALLROVIB(IV,J)
                  END IF
              END DO
            END IF
        END DO
      END DO
C     
C     TEMPERATURE
      T=5000.00D+00
C     SPIN MULTIPLICIY      
      GS=1.0D+00
C     ELECTRONIC DEGENERACY
      GE=1.0D+00
      
C      CALL TOTPARTFUNC(NLEV1,NJM,JDJR,IIV,IIJ,
C     *                 IVMAX,IROTMAX,ALLROVIB,
C     *                 GS,GE,T,Q)
     
C      WRITE(*,*) "THE TOTAL PARTITION 
C     *            FUNCTION AT TEMPERATURE",T,"IS",Q


C      CALL VIBPARTFUNC(NLEV1,NJM,JDJR,IIV,IIJ,
C     *                 IVMAX,IROTMAX,ALLROVIB,
C     *                 T,QV)

C      WRITE(*,*) "THE VIBRATIONAL PARTITION 
C     *            FUNCTION AT TEMPERATURE",T,"IS",QV

      RANDVIB=0.12783283209302385D+00
      RANDROT=0.46032184425332506D+00

C      CALL VIBCDF(NLEV1,NJM,JDJR,IIV,IIJ,
C     *            IVMAX,IROTMAX,ALLROVIB,
C     *            T,QV,RANDVIB,IVT)

C      CALL ROTPARTFUNC(NLEV1,NJM,JDJR,IIV,IIJ,
C     *                 IVMAX,IROTMAX,ALLROVIB,
C     *                 T,IVT,QR)

C      WRITE(*,*) "THE ROTATIONAL PARTITION 
C     *            FUNCTION AT TEMPERATURE",T,"IS",QR

C      CALL ROTCDF(NLEV1,NJM,JDJR,IIV,IIJ,
C     *            IVMAX,IROTMAX,ALLROVIB,
C     *            T,IVT,QR,RANDROT,IJT)

      NP=2
      NC=2
      M=NP
      N=NC 

      LWA=M*N+5*N+M
      TOL=SQRT(DPMPAR(1))

      C(1)=0.5D+00
      C(2)=1.5D+00

      W(1)=1.0D+00
      W(2)=1.0D+00

      X(1)=ALLROVIB(0,0)
      X(2)=ALLROVIB(0,0)

      Y(1)=0.0D+00
      Y(2)=0.0D+00

      CALL LMDIF1(CLASSTURN,M,N,C,FVEC,TOL,INFO,IWA,WA,LWA)

      CALL DIATPOT(C(1),V1,VLIM)

      CALL DIATPOT(C(2),V2,VLIM)          
      
!      WRITE(*,*) V1,C(1),V2,C(2)

   10 FORMAT(I3,I4,X,E14.8)
 
      END

C######################################################################################
C######################################################################################

      SUBROUTINE CLASSTURN(M,N,C,FVEC,IFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/DATARC/X,Y,W,VLIM
      PARAMETER(NPMAX=100)
      DIMENSION :: X(NPMAX),Y(NPMAX),W(NPMAX)
      DIMENSION :: C(N),FVEC(M),V(NPMAX)

      DO I=1,M
        CALL DIATPOT(C(I),V(I),VLIM)
        FF=X(I)-V(I)
        FVEC(I)=(Y(I)-FF)/W(I)
      END DO

      RETURN
      END 

C######################################################################################
C CALCULATING CUMULATIVE DISTRIBUTION FUNCTION FOR ROTATIONAL LEVELS USING BOLTZMANN DISTRIBUTION
C THIS SELECTS THE ROTATIONAL LEVEL OF A PREVIOSLY DETERMINED VIBRATIONAL LEVEL
C THAT MAKES THE CDF (SUM) TO BE APPROXIMATELY EQUAL TO A SELECTED RAND NUMBER
C###################################################################################### 

      SUBROUTINE ROTCDF(NLEV1,NJM,JDJR,IIV,IIJ,
     *                  IVMAX,IROTMAX,ALLROVIB,
     *                  T,IVT,QR,RANDROT,IJT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IVIBMX=400)
      DIMENSION :: IIV(IVIBMX),IIJ(IVIBMX)
      DIMENSION :: IROTMAX(IVIBMX)
      DIMENSION :: ALLROVIB(0:IVIBMX,0:IVIBMX)
      DIMENSION :: PROB(0:IVIBMX)
      DIMENSION :: CDFROT(0:IVIBMX)
      DIMENSION :: MINLOCROT(1)
C     BOLTZMANN CONSTANT IN J.K**-1 OR IN Kg.m**2.s**-2.K-1     
      BKB=1.38064852D-23
C     CONVERTING CM-1 TO J      
      CM2J=1.98644560548147D-23
C      
      SUM=0.0D+00

      CTE=1.0D+10
  
      DO I=0,IVIBMX
        PROB(I)=CTE
        CDFROT(I)=CTE
      END DO 
 
      DO K=1,IROTMAX(IVT+1)
        J=IIJ(1)+(K-1)*JDJR              
        GJ=(2.0D+00*DBLE(J)+1.0D+00)
        DELTA=ALLROVIB(IVT,J)-ALLROVIB(IVT,IIJ(1))
        SUM=SUM+GJ*EXP(-(DELTA*CM2J)/(BKB*T))
        PROB(J)=SUM/QR
        CDFROT(J)=ABS(PROB(J)-RANDROT)
        WRITE(*,*) IVT,J,ALLROVIB(IVT,J),
     *             GJ*EXP(-(DELTA*CM2J)/(BKB*T)),PROB(J),CDFROT(J)
      END DO
C
C     FINDS THE ROTATIONAL LEVEL IN THE VIBRATIONAL STATE PREVIOUSLY SELECTED 
C     THAT MAKES THE PROB CDF SUM APPROACHES THE RAND NUMBER  
C
      MINLOCROT=MINLOC(CDFROT,MASK=CDFROT.LT.CTE)-1
      IJT=MINLOCROT(1)

      WRITE(*,*) IVT,IJT,ALLROVIB(IVT,IJT)
        
      RETURN
      END 

C######################################################################################
C CALCULATING ROTATIONAL PARTITION FUNCTION FOR THE VIBRATIONAL LEVEL SELECT BY THE CDF
C######################################################################################      
      
      SUBROUTINE ROTPARTFUNC(NLEV1,NJM,JDJR,IIV,IIJ,
     *                       IVMAX,IROTMAX,ALLROVIB,
     *                       T,IVT,QR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IVIBMX=400)
      DIMENSION :: IIV(IVIBMX),IIJ(IVIBMX)
      DIMENSION :: IROTMAX(IVIBMX)
      DIMENSION :: ALLROVIB(0:IVIBMX,0:IVIBMX)
C     BOLTZMANN CONSTANT IN J.K**-1 OR IN Kg.m**2.s**-2.K-1     
      BKB=1.38064852D-23
C     CONVERTING CM-1 TO J      
      CM2J=1.98644560548147D-23
C      
      SUM=0.0D+00
     
      DO K=1,IROTMAX(IVT+1)
        J=IIJ(1)+(K-1)*JDJR              
        GJ=(2.0D+00*DBLE(J)+1.0D+00)
        DELTA=ALLROVIB(IVT,J)-ALLROVIB(IVT,IIJ(1))
        SUM=SUM+GJ*EXP(-(DELTA*CM2J)/(BKB*T))
      END DO

      QR=SUM
        
      RETURN
      END      

C######################################################################################
C CALCULATING CUMULATIVE DISTRIBUTION FUNCTION FOR VIBRATIONAL LEVELS USING BOLTZMANN DISTRIBUTION
C THIS SELECTS THE VIBRATIONAL LEVEL THAT MAKES THE CDF (SUM) TO BE APPROXIMATELY EQUAL TO 
C A SELECTED RAND NUMBER
C######################################################################################      
      
      SUBROUTINE VIBCDF(NLEV1,NJM,JDJR,IIV,IIJ,
     *                  IVMAX,IROTMAX,ALLROVIB,
     *                  T,QV,RANDVIB,IVT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IVIBMX=400)
      DIMENSION :: IIV(IVIBMX),IIJ(IVIBMX)
      DIMENSION :: IROTMAX(IVIBMX)
      DIMENSION :: ALLROVIB(0:IVIBMX,0:IVIBMX)
      DIMENSION :: PROB(0:IVMAX-1),CDFVIB(0:IVMAX-1)
      DIMENSION :: MINLOCVIB(1)
C     BOLTZMANN CONSTANT IN J.K**-1 OR IN Kg.m**2.s**-2.K-1     
      BKB=1.38064852D-23
C     CONVERTING CM-1 TO J      
      CM2J=1.98644560548147D-23
C      
      SUM=0.0D+00
C
C     THE CALCULATION ASSUMES THAT THE ENSEMBLE POPULATES ONLY THE GROUND ROTATIONAL LEVEL
C        
      DO I=1,IVMAX
        IV=I-1         
        DELTA=ALLROVIB(IV,IIJ(1))-ALLROVIB(IIV(1),IIJ(1))              
        SUM=SUM+EXP(-(DELTA*CM2J)/(BKB*T))
        PROB(IV)=SUM/QV
        CDFVIB(IV)=ABS(PROB(IV)-RANDVIB)
        WRITE(*,*) IV,ALLROVIB(IV,IIJ(1)),
     *             EXP(-(DELTA*CM2J)/(BKB*T)),PROB(IV),CDFVIB(IV) 
      END DO
C
C     FINDS THE VIBRATIONAL LEVEL THAT MAKES THE PROB CDF SUM APPROACHES THE RAND NUMBER  
C
      MINLOCVIB=MINLOC(CDFVIB)-1
      IVT=MINLOCVIB(1)

      WRITE(*,*) IVT,ALLROVIB(IVT,IIJ(1))

      RETURN
      END

C######################################################################################
C CALCULATING VIBRATIONAL PARTITION FUNCTION 
C######################################################################################      
      
      SUBROUTINE VIBPARTFUNC(NLEV1,NJM,JDJR,IIV,IIJ,
     *                       IVMAX,IROTMAX,ALLROVIB,
     *                       T,QV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IVIBMX=400)
      DIMENSION :: IIV(IVIBMX),IIJ(IVIBMX)
      DIMENSION :: IROTMAX(IVIBMX)
      DIMENSION :: ALLROVIB(0:IVIBMX,0:IVIBMX)
C     BOLTZMANN CONSTANT IN J.K**-1 OR IN Kg.m**2.s**-2.K-1     
      BKB=1.38064852D-23
C     CONVERTING CM-1 TO J      
      CM2J=1.98644560548147D-23
C      
      SUM=0.0D+00
C
C     THE CALCULATION ASSUMES THAT THE ENSEMBLE POPULATES ONLY THE GROUND ROTATIONAL LEVEL
C        
      DO I=1,IVMAX
        IV=I-1         
        DELTA=ALLROVIB(IV,IIJ(1))-ALLROVIB(IIV(1),IIJ(1))              
        SUM=SUM+EXP(-(DELTA*CM2J)/(BKB*T))
      END DO
      
      QV=SUM
        
      RETURN
      END

C######################################################################################
C CALCULATING TOTAL PARTITION FUNCTION (WITHOUT TRANSLATIONAL DEGREES OF FREEDOM)
C######################################################################################      
      
      SUBROUTINE TOTPARTFUNC(NLEV1,NJM,JDJR,IIV,IIJ,
     *                       IVMAX,IROTMAX,ALLROVIB,
     *                       GS,GE,T,Q)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IVIBMX=400)
      DIMENSION :: IIV(IVIBMX),IIJ(IVIBMX)
      DIMENSION :: IROTMAX(IVIBMX)
      DIMENSION :: ALLROVIB(0:IVIBMX,0:IVIBMX)
C     BOLTZMANN CONSTANT IN J.K**-1 OR IN Kg.m**2.s**-2.K-1     
      BKB=1.38064852D-23
C     CONVERTING CM-1 TO J      
      CM2J=1.98644560548147D-23
C      
      SUM=0.0D+00
      
      DO I=1,IVMAX
        IV=I-1
          IF (IROTMAX(I).NE.0) THEN
            DO K=1,IROTMAX(I)
            
              J=IIJ(1)+(K-1)*JDJR              
C     ROTATIONAL DEGENERACY FACTOR              
              GJ=(2.0D+00*DBLE(J)+1.0D+00)
C
              DELTA=ALLROVIB(IV,J)-ALLROVIB(IIV(1),IIJ(1))

              SUM=SUM+
     *        GE*GS*GJ*EXP(-(DELTA*CM2J)/(BKB*T))

            END DO
C      IN THE FOLLOWING CASE, WE HAVE ONLY THE VIBRATIONAL PART OF THE PARTITION FUNCTION             
          ELSE IF (IROTMAX(I).EQ.0) THEN
         
              DELTA=ALLROVIB(IV,IIJ(1))-ALLROVIB(IIV(1),IIJ(1))
              
              SUM=SUM+EXP(-(DELTA*CM2J)/(BKB*T))

          END IF  
      END DO
      
      Q=SUM
        
      RETURN
      END      

C######################################################################################
C DETERMINING THE DIATOMIC TO PERFORM THE ROVIBRATIONAL CALCULATIONS
C######################################################################################      
      
      SUBROUTINE DIATPOT(R,V,VLIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ANGS2BOHR=0.5291772D0    
      EH2CM=219474.6313702D+00 

      X=(R/ANGS2BOHR)
      
C
C     DEFINING THE POTENTIAL IN CM-1
C

      CALL CHIPR_DIAT(X,V,.FALSE.,DVDR)
C
      V=V*EH2CM+VLIM
     
      RETURN
      END 
