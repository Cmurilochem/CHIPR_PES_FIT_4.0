!###############################################################################
! CALCULATION OF THE FORCE CONSTANT FROM THE EXP. HARMONIC FRENQUENCY
!###############################################################################

      SUBROUTINE FORCE(M1,M2,VIB,SD)
      IMPLICIT NONE
      DOUBLE PRECISION :: M1, M2
      DOUBLE PRECISION :: REDMASS
      DOUBLE PRECISION :: SD, FCT, VIB
      DOUBLE PRECISION, PARAMETER :: C=2.99792458D+10             ! LIGHT VELOCIT IN cm.s-1
      DOUBLE PRECISION, PARAMETER :: B2A=0.5291772086D+00         ! BOHR TO ANGS
      DOUBLE PRECISION, PARAMETER :: JH=4.3597439D-18             ! CONVERTION FROM HARTREES TO JOULES
      DOUBLE PRECISION, PARAMETER :: AUKG=1.660538921E-27         ! CONVERTION FROM A.U TO Kg
      DOUBLE PRECISION, PARAMETER :: PI2=6.28318530718D+00        ! PI SQUARED
      REDMASS=((M1*M2)/(M1+M2))*AUKG                              ! REDUCED MASS IN KG
      FCT=((VIB*PI2*C)**2)*REDMASS                                ! FORCE CONSTANT IN J.m**-2
      SD=FCT/(JH/(B2A*1.00D-10)**2)                               ! SECOND DERIVATIVE IN Eh/BOHR**2
      END SUBROUTINE FORCE
    
!######################################################################################################     
! CALCULATES THE FORCE CONSTANT BY DIFFERENCIATING THE ANALYTIC FIRST DERIVATIVE
! IT USES THE RICHARDSON'S EXTRAPOLATION FORMULA
!###################################################################################################### 

      SUBROUTINE DIRECTFC(RE,FC,C)
      USE COMMON_VAR, ONLY : MBS,L,NC
      IMPLICIT NONE
      INTEGER :: M
      DOUBLE PRECISION, DIMENSION(NC) :: C
      DOUBLE PRECISION :: RE,FC,STEP
      DOUBLE PRECISION :: FMIN1, FMIN2
      DOUBLE PRECISION :: FMAX1, FMAX2
      DOUBLE PRECISION :: V1,V2,V3,V4
      LOGICAL :: DER
      DER=.TRUE.
      STEP=1.0D-4
      M=MBS(1)
      CALL DCHIPR_DIAT(M,L,C,(RE-STEP),V1,DER,FMIN1)
      CALL DCHIPR_DIAT(M,L,C,(RE+STEP),V2,DER,FMAX1)
      CALL DCHIPR_DIAT(M,L,C,(RE-2.0D+00*STEP),V3,DER,FMIN2)
      CALL DCHIPR_DIAT(M,L,C,(RE+2.0D+00*STEP),V4,DER,FMAX2) 
      FC=(FMIN2-8.0D+00*FMIN1+8.0D+00*FMAX1-FMAX2)/(12.0D+00*STEP)
      END SUBROUTINE DIRECTFC 
    
!######################################################################################################     
! CALCULATES THE EQUILIBRIUM DISTANCE AND DISSOCIATION ENERGY VIA BISECTION METHOD     
!######################################################################################################

      SUBROUTINE DIRECTMINDEDIAT(RE,DE,C)
      USE COMMON_VAR, ONLY : MBS,L,NC,REEXP
      IMPLICIT NONE 
      INTEGER :: IFLAG
      DOUBLE PRECISION :: FROOT,X1,X2,ROOT,EPS
      DOUBLE PRECISION :: RE,DE,POT
      DOUBLE PRECISION, DIMENSION(NC) :: C,DIFFC
      EXTERNAL FROOT
!      
! DEFINING THE INTERVAL TO FIND THE ROOT 
! SETTING THE GUESS AS THE EXPERIMENTAL DISTANCE
!
      X1=REEXP-0.5D+00
      X2=REEXP+0.5D+00
      EPS=1.0D-06
!
!  FINDING THE ROOT VIA BISECTION METHOD
!
      CALL BISECTION(FROOT,X1,X2,EPS,ROOT,IFLAG,C)
      IF (IFLAG.GT.0) THEN
        RE=ROOT
        CALL CHIPR_DIAT(MBS,L,C,RE,POT,DIFFC)
        DE=-POT
      ELSE
        WRITE(6,*) "PROBLEMS IN ROOT FINDING"
        WRITE(6,*) "SEE X1, X2 & EPS VARIABLES IN THE SUBROUTINE 'DIRECTMINDEDIAT' IN DIRECT_FIT.F90"
        STOP
      END IF
      END SUBROUTINE DIRECTMINDEDIAT
  
!######################################################################################################
! CALCULATES THE ANALYTIC FIRST DERIVATIVE AND FIND THE POINT AT WHICH IT IS ZERO
!######################################################################################################      
      
      DOUBLE PRECISION FUNCTION FROOT(R,C)
      USE COMMON_VAR, ONLY : MBS,L,NC
      IMPLICIT NONE
      INTEGER :: M
      LOGICAL :: DER
      DOUBLE PRECISION :: R,POT,DVDR
      DOUBLE PRECISION, DIMENSION(NC) :: C
      DER=.TRUE.
      M=MBS(1)
      CALL DCHIPR_DIAT(M,L,C,R,POT,DER,DVDR)
      FROOT=DVDR
      END FUNCTION      
      
!######################################################################################
! DETERMINING THE POTENTIAL TO PERFORM THE ROVIBRATIONAL CALCULATIONS IN LEVEL
! R IN ANGSTROM
! POTENTIAL IN CM-1
!######################################################################################      
      
      SUBROUTINE DIATPOT(R,V,VLIM,C)
      USE COMMON_VAR, ONLY : EH2CM,ANGS2AU
      USE COMMON_VAR, ONLY : MBS,L,NC
      IMPLICIT NONE
      DOUBLE PRECISION :: VLIM
      DOUBLE PRECISION :: R,X,V
      DOUBLE PRECISION, DIMENSION(NC) :: C,DIFFC
!
!     CONVERTING ANGSTROM TO BOHR
!
      X=R/ANGS2AU      
!
!     DEFINING THE POTENTIAL IN CM-1
!
      CALL CHIPR_DIAT(MBS,L,C,X,V,DIFFC)
!
      V=V*EH2CM+VLIM
     
      RETURN
      END       
           
!######################################################################################################
! CALCULATES ANALYTIC FIRST DERIVATIVE OF THE ACTUAL CURVE
!######################################################################################################

      SUBROUTINE DCHIPR_DIAT(BSORDER,POLORDER,C,R,POT,DER,DVDR)
      USE COMMON_VAR, ONLY : Z,DEG,NC,NX
      IMPLICIT NONE
      INTEGER :: I,J
      INTEGER :: BSORDER
      INTEGER :: POLORDER
      INTEGER :: NCBAS,NCPOL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BS
      DOUBLE PRECISION :: Y   
      DOUBLE PRECISION :: R,POT
      LOGICAL :: DER
      DOUBLE PRECISION :: DYDR
      DOUBLE PRECISION :: DVDRPART1
      DOUBLE PRECISION :: DVDRPART2
      DOUBLE PRECISION :: DVDR
      DOUBLE PRECISION, DIMENSION(NC) :: C

      NCBAS=2*BSORDER+2

      NCPOL=POLORDER
!
!     CHECKING - JUST IN CASE
!
      IF (NC.NE.(NCPOL+NCBAS)) THEN 
        WRITE(6,*) "PROBLEMS IN DEFINING PROPER NUMBER OF COEFFICIENTS"
        WRITE(6,*) "TOTAL NUMBER OF COEFFS ARE NOT SUMMING UP CORRECTLY"
        STOP
      END IF      
      
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

!###### CALCULATING ANALYTIC DERIVATIVES IF DER=.TRUE.

      IF (DER) THEN 
        CALL DBASIS_CONTRACT(1,BSORDER,BS,R,DYDR)
        DVDRPART1=-POT/R
        DVDRPART2=0.00D+00
        DO I=1,POLORDER
          DVDRPART2=DVDRPART2+(Z(1)*Z(2)/R)*C(I)*DBLE(I)*Y**(DBLE(I-1))
        END DO
        DVDRPART2=DVDRPART2*DYDR
        DVDR=DVDRPART1+DVDRPART2
      ELSE
        DVDR=1000
      ENDIF

!######

      DEALLOCATE(BS)

      RETURN
      END      
      
!######################################################################################################

      SUBROUTINE DBASIS_CONTRACT(DEG,M,C,R,DYDR)
      IMPLICIT NONE
      INTEGER :: I,J,DEG
      INTEGER :: M
      DOUBLE PRECISION :: RREF0,ZETA
      DOUBLE PRECISION, DIMENSION(2*M+2) :: C
      DOUBLE PRECISION, DIMENSION(M) :: GAMA,DPHIDR
      DOUBLE PRECISION :: R,DYDR

      DO I=1,M
        GAMA(I)=C(M+I)
      END DO

      RREF0=C(2*M+1)
      ZETA=C(2*M+2)

      DO I=1,M-1
        CALL DPHISECBASIS(DEG,I,RREF0,ZETA,GAMA(I),1,R,DPHIDR(I))
      END DO

      DO I=M,M
        CALL DPHICSECBASIS(DEG,I,RREF0,ZETA,GAMA(I),1,6,R,DPHIDR(I))
      END DO

      DYDR=0.00D+00

      DO J=1,M
        DYDR=DYDR+C(J)*DPHIDR(J)
      END DO

      RETURN
      END

!######################################################################################################

      SUBROUTINE DPHISECBASIS(DEG,IND,RREF0,ZETA,GAMA,ETA,R,DPHIDR)
      IMPLICIT NONE
      INTEGER :: DEG, IND, ETA
      DOUBLE PRECISION :: RREF0, ZETA, RREFIND 
      DOUBLE PRECISION :: GAMA, R, VAL, RHO
      DOUBLE PRECISION :: SECH, ORIG
      DOUBLE PRECISION :: PART1,PART2
      DOUBLE PRECISION :: DPHIDR
      RREFIND=ORIG(IND,ZETA,RREF0)
      RHO=(R-RREFIND)
      PART1=(SECH(GAMA*RHO))**(DBLE(ETA))
      PART2=(TANH(GAMA*RHO))
      DPHIDR=-(DBLE(ETA))*GAMA*PART1*PART2
      RETURN
      END

!######################################################################################################

      SUBROUTINE DPHICSECBASIS(DEG,IND,RREF0,ZETA,GAMA,ETA,LR,R,DPHIDR)
      IMPLICIT NONE
      INTEGER :: DEG, IND, LR, ETA
      DOUBLE PRECISION :: RREF0, ZETA, RREFIND 
      DOUBLE PRECISION :: GAMA, R, VAL, RHO
      DOUBLE PRECISION :: SECH, BETA, FAC, ORIG
      DOUBLE PRECISION :: FAC1PART1,FAC2PART1
      DOUBLE PRECISION :: FAC3PART1,FAC1PART2
      DOUBLE PRECISION :: FAC2PART2,FAC3PART2
      DOUBLE PRECISION :: PART1,PART2
      DOUBLE PRECISION :: DPHIDR
      BETA=1.00D+00/5.0D+00
      RREFIND=ORIG(IND,ZETA,RREF0)
      RHO=(R-RREFIND)
      FAC1PART1=(TANH(BETA*R)/R)**(DBLE(LR-1))
      FAC2PART1=(BETA*(SECH(BETA*R))**(2)/R)-(TANH(BETA*R)/R**2)
      FAC3PART1=(SECH(GAMA*RHO))**(DBLE(ETA))
      PART1=(DBLE(LR))*FAC1PART1*FAC2PART1*FAC3PART1
      FAC1PART2=(SECH(GAMA*RHO))**(DBLE(ETA))
      FAC2PART2=(TANH(GAMA*RHO))
      FAC3PART2=(TANH(BETA*R)/R)**(DBLE(LR))
      PART2=-(DBLE(ETA))*GAMA*FAC1PART2*FAC2PART2*FAC3PART2
      DPHIDR=PART1+PART2
      RETURN
      END      

!###############################################################################
! CALCULATION OF THE FORCE CONSTANT FROM THE HARMONIC VIBRATIONAL FRENQUENCY
!###############################################################################      
      SUBROUTINE BISECTION(F,X1,X2,EPS,ROOT,FLAG,COEFFS)
      USE COMMON_VAR, ONLY : NC
!============================================================
! Solutions of equation f(x)=0 on [x1,x2] interval
! Method: Bisectional (closed domain) (a single root)
! Alex G. January 2010
!------------------------------------------------------------
! input ...
! f   - function - evaluates f(x) for any x in [x1,x2]
! x1  - left endpoint of initial interval
! x2  - right endpoint of initial interval
! eps - desired uncertainity of the root as |b-a|<eps
! output ...
! Root  - root of the equation f(x)=0
! flag  - indicator of success
!         >0 - a single root found, flag=number of iterations
!          0 - no solutions for the bisectional method
!         <0 - not a root but singularity, flag=number of iterations
!
! Comments: Function f(x) has to change sign between x1 and x2
!           Max number of iterations - 200 (accuracy (b-a)/2**200)
!====================================================================
       implicit none
       double precision f, x1, x2, eps, Root
       double precision a, b, c
       integer i, flag
       integer, parameter:: iter=200
       DOUBLE PRECISION, DIMENSION(NC) :: COEFFS
       
       !* check the bisection condition
       if(f(x1,COEFFS)*f(x2,COEFFS)>0.0) then
         flag = 0
         return
       end if
       
       !* initialize calculations
       a=x1
       b=x2
       
       !* Iterative refining the solution 
       do i=1,iter
         c=(b+a)/2.0
         if(f(a,COEFFS)*f(c,COEFFS).le.0.0) then
             b = c
           else
             a=c
         end if
       ! condition(s) to stop iterations)
         if(abs(b-a)<= eps) exit  
       end do
       Root=(b+a)/2.0
       
       !* check if it is a root or singularity
       if (abs(f(Root,COEFFS)) < 1.0) then
         flag=i
         else
         flag = -i
       end if
       END SUBROUTINE BISECTION      
