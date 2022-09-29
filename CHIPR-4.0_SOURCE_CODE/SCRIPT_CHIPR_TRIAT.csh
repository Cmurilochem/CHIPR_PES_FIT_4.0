#!/bin/csh

set MOLTYP=`grep "MOLTYP" inp_data.res | awk '{print $2}'`

set DEG=`grep "DEG" inp_data.res | awk '{print $2}'`

set NC=`grep "NC" inp_data.res | awk '{print $2}'`

set NX=`grep "NX" inp_data.res | awk '{print $2}'`

cat <<EOF > CHIPR_TRIAT_FUNC.f90
!######################################################################################################

      SUBROUTINE CHIPR_TRIAT(R,POT,DER,DVDR)
      IMPLICIT NONE
      CHARACTER(LEN=3), PARAMETER :: MOLTYP="$MOLTYP"
      INTEGER, PARAMETER :: DEG=$DEG
      INTEGER, PARAMETER :: NC=$NC
      INTEGER, PARAMETER :: NX=$NX
      INTEGER :: I,J,K,L,M,S,O,ID
      INTEGER :: POLORDER
      INTEGER, DIMENSION(DEG) :: BSORDER,NCBAS
      DOUBLE PRECISION, DIMENSION(NC) :: C
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: POL,BS
      DOUBLE PRECISION, DIMENSION(3) :: R
      DOUBLE PRECISION :: POT,REPDAMP
      INTEGER :: NCPOL,NCTOTAL,SUMC
      DOUBLE PRECISION, DIMENSION(NX) :: Y
      INTEGER :: TOTNUM,MNUM
      INTEGER, DIMENSION(3,6) :: P
      LOGICAL :: DER
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DVDY1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DVDY2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DVDY3
      DOUBLE PRECISION, DIMENSION(NX) :: DVDR,DYDR

EOF

grep "C(" inp_data.res >> CHIPR_TRIAT_FUNC.f90

echo "" >> CHIPR_TRIAT_FUNC.f90

grep "BSORDER(" inp_data.res >> CHIPR_TRIAT_FUNC.f90

echo "" >> CHIPR_TRIAT_FUNC.f90

grep "POLORDER" inp_data.res >> CHIPR_TRIAT_FUNC.f90

cat <<EOF >> CHIPR_TRIAT_FUNC.f90

      DO I=1,DEG
        NCBAS(I)=0
        NCBAS(I)=2*BSORDER(I)+2
      END DO

      CALL POLNC(POLORDER,MOLTYP,NCPOL)

      IF (NC.NE.(NCPOL+SUM(NCBAS))) THEN 
        WRITE(*,*) "PROBLEMS IN DEFINING PROPER NUMBER OF COEFFICIENTS"
        WRITE(*,*) "TOTAL NUMBER OF COEFFS ARE NOT SUMMING UP CORRECTLY"
        STOP
      END IF

      ALLOCATE(POL(NCPOL))
      ALLOCATE(DVDY1(NCPOL),DVDY2(NCPOL),DVDY3(NCPOL))

      DO I=1,NCPOL
        POL(I)=0.00D+00
        DVDY1(I)=0.00D+00
        DVDY2(I)=0.00D+00
        DVDY3(I)=0.00D+00
      END DO

      DO I=1,NX
        Y(I)=0.0D+00
        DYDR(I)=0.0D+00
        DVDR(I)=0.0D+00
      END DO

      SUMC=NCPOL
      DO I=1,DEG
        ALLOCATE(BS(NCBAS(I)))
        DO J=1,NCBAS(I)
          BS(J)=0.0D+00
        END DO
        K=SUMC
        DO J=1,NCBAS(I)
          K=K+1
          BS(J)=C(K)
        END DO

EOF

if ($DEG == 1) then 

cat <<EOF >> CHIPR_TRIAT_FUNC.f90
!######
! A3-TYPE
!######
        IF (DEG.EQ.1) THEN
          DO O=1,3 
            CALL BASIS_CONTRACT(1,BSORDER(1),BS,R(O),Y(O))

!###### CALCULATING ANALYTIC DERIVATIVES OF Y WITH RESPECT TO R IF DER=.TRUE.
            
            IF (DER) THEN
              CALL DBASIS_CONTRACT(1,BSORDER(1),BS,R(O),DYDR(O))
            END IF

!######
          END DO
        END IF

EOF

else if ($DEG == 2) then 

cat <<EOF >> CHIPR_TRIAT_FUNC.f90
!######
! AB2-TYPE
!######
        IF (DEG.EQ.2) THEN
          IF (I.EQ.1) THEN
!
!           B-B BASIS
!
            CALL BASIS_CONTRACT(1,BSORDER(1),BS,R(1),Y(1))

!###### CALCULATING ANALYTIC DERIVATIVES OF Y WITH RESPECT TO R IF DER=.TRUE.
            
            IF (DER) THEN
              CALL DBASIS_CONTRACT(1,BSORDER(1),BS,R(1),DYDR(1))
            END IF

!###### 
          ELSE
!
!           A-B BASIS
!
            DO O=2,3 
              CALL BASIS_CONTRACT(2,BSORDER(2),BS,R(O),Y(O))

!###### CALCULATING ANALYTIC DERIVATIVES OF Y WITH RESPECT TO R IF DER=.TRUE.
            
              IF (DER) THEN
                CALL DBASIS_CONTRACT(2,BSORDER(2),BS,R(O),DYDR(O))
              END IF

!######
            END DO
          END IF            
        END IF

EOF

else if ($DEG == 3) then 

cat <<EOF >> CHIPR_TRIAT_FUNC.f90
!######
! ABC-TYPE
!######
        IF (DEG.EQ.3) THEN
          CALL BASIS_CONTRACT(I,BSORDER(I),BS,R(I),Y(I))

!###### CALCULATING ANALYTIC DERIVATIVES OF Y WITH RESPECT TO R IF DER=.TRUE.
            
          IF (DER) THEN
            CALL DBASIS_CONTRACT(I,BSORDER(I),BS,R(I),DYDR(I))
          END IF

!######
        END IF

EOF

endif

cat <<EOF >> CHIPR_TRIAT_FUNC.f90
        SUMC=SUMC+NCBAS(I)
        DEALLOCATE(BS)
      END DO

      TOTNUM=0
      S=0
      DO M=0,POLORDER
        MNUM=0
        DO I=0,M
          DO J=0,(M-I)
            K=M-I-J
            S=I+J+K
            IF (S.NE.I .AND. S.NE.J .AND. S.NE.K) THEN

              IF (MOLTYP.EQ."ABC") THEN  

                TOTNUM=TOTNUM+1
                MNUM=MNUM+1

                CALL PERMUTABC(I,J,K,P,ID)

                DO L=1,ID

                  POL(TOTNUM)=POL(TOTNUM)+&
     &               (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

!###### CALCULATING ANALYTIC DERIVATIVES OF V WITH RESPECT TO Y IF DER=.TRUE.

                      IF (DER) THEN
             
                        DVDY1(TOTNUM)=DVDY1(TOTNUM)+&
     &                 (DBLE(P(1,L))*Y(1)**(P(1,L)-1)*Y(2)**(P(2,L))*Y(3)**(P(3,L))) 

                        DVDY2(TOTNUM)=DVDY2(TOTNUM)+&
     &                 (Y(1)**(P(1,L))*DBLE(P(2,L))*Y(2)**(P(2,L)-1)*Y(3)**(P(3,L)))    

                        DVDY3(TOTNUM)=DVDY3(TOTNUM)+&
     &                 (Y(1)**(P(1,L))*Y(2)**(P(2,L))*DBLE(P(3,L))*Y(3)**(P(3,L)-1)) 
             
                      END IF

!######
                END DO

                POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

!###### MULTIPLYING THE CONTRIBUTIONS OF DVDY BY THE POLYNOMIAL COEFFS. IF DER=.TRUE.

                IF (DER) THEN
                  DVDY1(TOTNUM)=C(TOTNUM)*DVDY1(TOTNUM)/DBLE(ID)
                  DVDY2(TOTNUM)=C(TOTNUM)*DVDY2(TOTNUM)/DBLE(ID)
                  DVDY3(TOTNUM)=C(TOTNUM)*DVDY3(TOTNUM)/DBLE(ID)
                END IF

!######
              ELSE IF (MOLTYP.EQ."AB2") THEN 

                IF (J.LE.K) THEN 

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTAB2(I,J,K,P,ID)
                   
                  DO L=1,ID

                    POL(TOTNUM)=POL(TOTNUM)+&
     &                 (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

!###### CALCULATING ANALYTIC DERIVATIVES OF V WITH RESPECT TO Y IF DER=.TRUE.

                      IF (DER) THEN
             
                        DVDY1(TOTNUM)=DVDY1(TOTNUM)+&
     &                 (DBLE(P(1,L))*Y(1)**(P(1,L)-1)*Y(2)**(P(2,L))*Y(3)**(P(3,L))) 

                        DVDY2(TOTNUM)=DVDY2(TOTNUM)+&
     &                 (Y(1)**(P(1,L))*DBLE(P(2,L))*Y(2)**(P(2,L)-1)*Y(3)**(P(3,L)))    

                        DVDY3(TOTNUM)=DVDY3(TOTNUM)+&
     &                 (Y(1)**(P(1,L))*Y(2)**(P(2,L))*DBLE(P(3,L))*Y(3)**(P(3,L)-1)) 
             
                      END IF

!######
                  END DO

                  POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

!###### MULTIPLYING THE CONTRIBUTIONS OF DVDY BY THE POLYNOMIAL COEFFS. IF DER=.TRUE.

                  IF (DER) THEN
                    DVDY1(TOTNUM)=C(TOTNUM)*DVDY1(TOTNUM)/DBLE(ID)
                    DVDY2(TOTNUM)=C(TOTNUM)*DVDY2(TOTNUM)/DBLE(ID)
                    DVDY3(TOTNUM)=C(TOTNUM)*DVDY3(TOTNUM)/DBLE(ID)
                  END IF

!######
                END IF    
         
              ELSE IF (MOLTYP.EQ."A3") THEN

                IF (I.LE.J .AND. J.LE.K) THEN

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTA3(I,J,K,P,ID)
                   
                  DO L=1,ID

                    POL(TOTNUM)=POL(TOTNUM)+&
     &                 (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

!###### CALCULATING ANALYTIC DERIVATIVES OF V WITH RESPECT TO Y IF DER=.TRUE.

                      IF (DER) THEN
             
                        DVDY1(TOTNUM)=DVDY1(TOTNUM)+&
     &                 (DBLE(P(1,L))*Y(1)**(P(1,L)-1)*Y(2)**(P(2,L))*Y(3)**(P(3,L))) 

                        DVDY2(TOTNUM)=DVDY2(TOTNUM)+&
     &                 (Y(1)**(P(1,L))*DBLE(P(2,L))*Y(2)**(P(2,L)-1)*Y(3)**(P(3,L)))    

                        DVDY3(TOTNUM)=DVDY3(TOTNUM)+&
     &                 (Y(1)**(P(1,L))*Y(2)**(P(2,L))*DBLE(P(3,L))*Y(3)**(P(3,L)-1)) 
             
                      END IF

!######
                  END DO

                  POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

!###### MULTIPLYING THE CONTRIBUTIONS OF DVDY BY THE POLYNOMIAL COEFFS. IF DER=.TRUE.

                  IF (DER) THEN
                    DVDY1(TOTNUM)=C(TOTNUM)*DVDY1(TOTNUM)/DBLE(ID)
                    DVDY2(TOTNUM)=C(TOTNUM)*DVDY2(TOTNUM)/DBLE(ID)
                    DVDY3(TOTNUM)=C(TOTNUM)*DVDY3(TOTNUM)/DBLE(ID)
                  END IF

!######
                END IF

              END IF

            END IF
          END DO
        END DO

      END DO

      POT=SUM(POL)*REPDAMP(NX,R)

!###### CALCULATING ANALYTIC DERIVATIVES OF V WITH RESPECT TO R IF DER=.TRUE.

      IF (DER) THEN
        DVDR(1)=SUM(DVDY1)*DYDR(1)
        DVDR(2)=SUM(DVDY2)*DYDR(2)
        DVDR(3)=SUM(DVDY3)*DYDR(3)
      ELSE 
        DVDR(1)=1000
        DVDR(2)=1000
        DVDR(3)=1000
      END IF

!######

      DEALLOCATE(POL)
      DEALLOCATE(DVDY1,DVDY2,DVDY3)

      RETURN
      END

!####################################################################################

      SUBROUTINE CARTDERPES3BD(X,DVDX)
      IMPLICIT NONE 
      INTEGER :: I
      INTEGER, PARAMETER :: NATOM=3     
      DOUBLE PRECISION, DIMENSION(3*NATOM) :: X
      DOUBLE PRECISION, DIMENSION(3*NATOM) :: DVDX,O
      DOUBLE PRECISION, DIMENSION(NATOM) :: Y,R,DVDR
      DOUBLE PRECISION :: V
!
!     This subprogram uses the chain rule to calculate the derivatives of the
!     energy with respect to the cartesian coordinates from the derivatives
!     with respect to the internal coordinates for a three-body system.
!     The convention assumed in this program is as follows:
!
!     R(1) = R(A-B)
!     R(2) = R(A-C)
!     R(3) = R(B-C)
!     X(1) - X(3) : X, Y, Z for atom A
!     X(4) - X(6) : X, Y, Z for atom B
!     X(7) - X(9) : X, Y, Z for atom C
!
      CALL CART2INTER3BD(X,O,3*NATOM)

      R(1)=O(1)
      R(2)=O(2)
      R(3)=O(3)

      CALL CHIPR_TRIAT(R,V,.TRUE.,DVDR) 

      DO I=1,3
        Y(I)=DVDR(I)/R(I)
      END DO

      DVDX(1)=-Y(1)*(X(4)-X(1))-Y(2)*(X(7)-X(1)) 
      DVDX(2)=-Y(1)*(X(5)-X(2))-Y(2)*(X(8)-X(2)) 
      DVDX(3)=-Y(1)*(X(6)-X(3))-Y(2)*(X(9)-X(3)) 
      DVDX(4)= Y(1)*(X(4)-X(1))-Y(3)*(X(7)-X(4)) 
      DVDX(5)= Y(1)*(X(5)-X(2))-Y(3)*(X(8)-X(5))
      DVDX(6)= Y(1)*(X(6)-X(3))-Y(3)*(X(9)-X(6))
      DVDX(7)= Y(2)*(X(7)-X(1))+Y(3)*(X(7)-X(4))
      DVDX(8)= Y(2)*(X(8)-X(2))+Y(3)*(X(8)-X(5))
      DVDX(9)= Y(2)*(X(9)-X(3))+Y(3)*(X(9)-X(6))

      RETURN
      END

!####################################################################################

      SUBROUTINE CART2INTER3BD(X,R,NTA)
      IMPLICIT NONE
      INTEGER, PARAMETER :: NA=3,NT=3*NA
      INTEGER :: NTA,I,K,J
      DOUBLE PRECISION, DIMENSION (NT):: X,R
      K=0
      DO I=1,NTA,3
       DO J=I,NTA,3
        IF(I.NE.J)THEN
         K=K+1
         R(K)=SQRT((X(I)-X(J))**2+(X(I+1)-X(J+1))**2+&
     &   (X(I+2)-X(J+2))**2)
        END IF
       END DO
      END DO
      RETURN
      END

!################################################################################

      DOUBLE PRECISION FUNCTION REPDAMP(NX,R)
      IMPLICIT NONE
      INTEGER :: I,NX
      DOUBLE PRECISION, DIMENSION(NX) :: R,H
      DOUBLE PRECISION :: KAPPA, XI, R0
      R0=0.5D+00
      KAPPA=100.0D+00
      XI=10.0D+00
      DO I=1,NX
        H(I)=0.5D+00*(1.00D+00+TANH(KAPPA*(R(I)-R0)))
      END DO
      REPDAMP=(PRODUCT(H))**(XI)      
      END FUNCTION

!################################################################################      
      
      SUBROUTINE PERMUTABC(I,J,K,P,ID)
      IMPLICIT NONE
      INTEGER :: I,J,K
      INTEGER :: L,M,ID
      INTEGER, DIMENSION(3) :: INTER
      INTEGER, DIMENSION(3,6) :: P

      DO L=1,3 
        DO M=1,6
          P(L,M)=0
        END DO
      END DO

      INTER(1)=I
      INTER(2)=J
      INTER(3)=K

      ID=1

      P(1,1)=INTER(1)
      P(2,1)=INTER(2)
      P(3,1)=INTER(3) 

      RETURN
      END

!################################################################################      
      
      SUBROUTINE PERMUTAB2(I,J,K,P,ID)
      IMPLICIT NONE
      INTEGER :: I,J,K
      INTEGER :: L,M,ID
      INTEGER, DIMENSION(3) :: INTER
      INTEGER, DIMENSION(3,6) :: P

      DO L=1,3 
        DO M=1,6
          P(L,M)=0
        END DO
      END DO
      
      INTER(1)=I
      INTER(2)=J
      INTER(3)=K
      
      IF (J.NE.K) THEN  
 
        ID=2

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)      

        P(1,2)=INTER(1)
        P(2,2)=INTER(3)
        P(3,2)=INTER(2) 

       ELSE 

        ID=1

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)

       END IF

      RETURN 
      END 

!################################################################################      
      
      SUBROUTINE PERMUTA3(I,J,K,P,ID)
      IMPLICIT NONE
      INTEGER :: I,J,K
      INTEGER :: L,M,ID
      INTEGER, DIMENSION(3) :: INTER
      INTEGER, DIMENSION(3,6) :: P

      DO L=1,3 
        DO M=1,6
          P(L,M)=0
        END DO
      END DO
      
      INTER(1)=I
      INTER(2)=J
      INTER(3)=K
      
      IF (I.EQ.J.AND.J.EQ.K.AND.I.EQ.K) THEN
      
        ID=1

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)  
        
      ELSE IF (I.NE.J.AND.J.NE.K.AND.I.NE.K) THEN   

        ID=6

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)      

        P(1,2)=INTER(2)
        P(2,2)=INTER(1)
        P(3,2)=INTER(3)  

        P(1,3)=INTER(3)
        P(2,3)=INTER(2)
        P(3,3)=INTER(1) 

        P(1,4)=INTER(1)
        P(2,4)=INTER(3)
        P(3,4)=INTER(2) 

        P(1,5)=P(1,2)
        P(2,5)=P(3,2)
        P(3,5)=P(2,2)  

        P(1,6)=P(1,3)
        P(2,6)=P(3,3)
        P(3,6)=P(2,3)  
        
      ELSE IF (I.EQ.J) THEN       
     
        ID=3     

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3) 

        P(1,2)=INTER(3)
        P(2,2)=INTER(2)
        P(3,2)=INTER(1) 

        P(1,3)=INTER(1)
        P(2,3)=INTER(3)
        P(3,3)=INTER(2)
        
      ELSE IF (J.EQ.K) THEN       
                    
        ID=3     

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)      

        P(1,2)=INTER(2)
        P(2,2)=INTER(1)
        P(3,2)=INTER(3)  

        P(1,3)=INTER(3)
        P(2,3)=INTER(2)
        P(3,3)=INTER(1) 
        
      ELSE IF (I.EQ.K) THEN    
  
        ID=3  

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)      

        P(1,2)=INTER(2)
        P(2,2)=INTER(1)
        P(3,2)=INTER(3)  

        P(1,3)=INTER(1)
        P(2,3)=INTER(3)
        P(3,3)=INTER(2)        
    
      END IF

      RETURN 
      END 

!################################################################################ 

      SUBROUTINE POLNC(ORDER,MOLTYP,NC)
      IMPLICIT NONE
      INTEGER :: I,J,K,L,M,S
      INTEGER :: NC,MNUM,ORDER
      CHARACTER(LEN=3) :: MOLTYP

      NC=0
      S=0

      DO M=0,ORDER
        MNUM=0
        DO I=0,M
          DO J=0,(M-I)
            K=M-I-J
            S=I+J+K
            IF (S.NE.I .AND. S.NE.J .AND. S.NE.K) THEN
 
              IF (MOLTYP.EQ."ABC") THEN  

                NC=NC+1
                MNUM=MNUM+1

              ELSE IF (MOLTYP.EQ."AB2") THEN 

                IF (J.LE.K) THEN 

                  NC=NC+1
                  MNUM=MNUM+1

                END IF    
        
              ELSE IF (MOLTYP.EQ."A3") THEN

                IF (I.LE.J .AND. J.LE.K) THEN

                  NC=NC+1
                  MNUM=MNUM+1

                END IF

              END IF

            END IF
          END DO
        END DO
      END DO

      RETURN
      END

!######################################################################################################

      SUBROUTINE BASIS_CONTRACT(DEG,M,C,R,YVAL)
      IMPLICIT NONE
      INTEGER :: I,J,DEG
      INTEGER :: M
      DOUBLE PRECISION :: RREF0,ZETA
      DOUBLE PRECISION, DIMENSION(2*M+2) :: C
      DOUBLE PRECISION, DIMENSION(M) :: GAMA,VAL
      DOUBLE PRECISION :: R,YVAL

      DO I=1,M
        GAMA(I)=C(M+I)
      END DO

      RREF0=C(2*M+1)
      ZETA=C(2*M+2)

      DO I=1,M-1
        CALL PHISECBASIS(DEG,I,RREF0,ZETA,GAMA(I),1,R,VAL(I))
      END DO

      DO I=M,M
        CALL PHICSECBASIS(DEG,I,RREF0,ZETA,GAMA(I),1,6,R,VAL(I))
      END DO

      YVAL=0.00D+00

      DO J=1,M
        YVAL=YVAL+C(J)*VAL(J)
      END DO

      RETURN
      END

!######################################################################################################

      SUBROUTINE PHISECBASIS(DEG,IND,RREF0,ZETA,GAMA,ETA,R,VAL)
      IMPLICIT NONE
      INTEGER :: DEG, IND, ETA
      DOUBLE PRECISION :: RREF0, ZETA, RREFIND 
      DOUBLE PRECISION :: GAMA, R, VAL, RHO
      DOUBLE PRECISION :: SECH, ORIG
      RREFIND=ORIG(IND,ZETA,RREF0)
      RHO=(R-RREFIND)
      VAL=(SECH(GAMA*RHO))**(DBLE(ETA))
      RETURN
      END

!######################################################################################################

      SUBROUTINE PHICSECBASIS(DEG,IND,RREF0,ZETA,GAMA,ETA,LR,R,VAL)
      IMPLICIT NONE
      INTEGER :: DEG, IND, LR, ETA
      DOUBLE PRECISION :: RREF0, ZETA, RREFIND 
      DOUBLE PRECISION :: GAMA, R, VAL, RHO
      DOUBLE PRECISION :: SECH, BETA, FAC, ORIG
      BETA=1.00D+00/5.0D+00
      RREFIND=ORIG(IND,ZETA,RREF0)
      RHO=(R-RREFIND)
      FAC=(TANH(BETA*R)/R)**(DBLE(LR))
      VAL=FAC*(SECH(GAMA*RHO))**(DBLE(ETA))
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

!######################################################################################################

      DOUBLE PRECISION FUNCTION ORIG(IND,ZETA,RREF0)
      IMPLICIT NONE
      INTEGER :: IND
      DOUBLE PRECISION :: ZETA,RREF0
      ORIG=ZETA*(RREF0)**(DBLE(IND)-1.0D+00)
      RETURN
      END

!######################################################################################################

      DOUBLE PRECISION FUNCTION SECH(X)
      IMPLICIT NONE
      DOUBLE PRECISION :: X
      SECH=1.00D+00/(COSH(X))      
      RETURN
      END
EOF

exit
