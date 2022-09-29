#!/bin/csh

set MOLTYP=`grep "MOLTYP" inp_data.res | awk '{print $2}'`

set DEG=`grep "DEG" inp_data.res | awk '{print $2}'`

set NC=`grep "NC" inp_data.res | awk '{print $2}'`

set NX=`grep "NX" inp_data.res | awk '{print $2}'`

cat <<EOF > CHIPR_TETRA_FUNC.f90
!######################################################################################################

      SUBROUTINE CHIPR_TETRA(R,POT,DER,DVDR)
      IMPLICIT NONE
      CHARACTER(LEN=4), PARAMETER :: MOLTYP="$MOLTYP"
      INTEGER, PARAMETER :: DEG=$DEG
      INTEGER, PARAMETER :: NC=$NC
      INTEGER, PARAMETER :: NX=$NX
      INTEGER :: I,J,K,L,M,N,O,Q,S,LL,ID
      INTEGER :: POLORDER
      INTEGER, DIMENSION(DEG) :: BSORDER,NCBAS
      DOUBLE PRECISION, DIMENSION(NC) :: C
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: POL,BS
      DOUBLE PRECISION, DIMENSION(6) :: R
      DOUBLE PRECISION :: POT,REPDAMPTETRA
      INTEGER :: NCPOL,NCTOTAL,SUMC
      DOUBLE PRECISION, DIMENSION(NX) :: Y
      INTEGER :: TOTNUM,MNUM
      INTEGER, DIMENSION(6,24) :: P
      LOGICAL :: DER
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DVDY1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DVDY2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DVDY3
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DVDY4
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DVDY5
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DVDY6      
      DOUBLE PRECISION, DIMENSION(NX) :: DVDR,DYDR 

EOF

grep "C(" inp_data.res >> CHIPR_TETRA_FUNC.f90

echo "" >> CHIPR_TETRA_FUNC.f90

grep "BSORDER(" inp_data.res >> CHIPR_TETRA_FUNC.f90

echo "" >> CHIPR_TETRA_FUNC.f90

grep "POLORDER" inp_data.res >> CHIPR_TETRA_FUNC.f90

cat <<EOF >> CHIPR_TETRA_FUNC.f90

      DO I=1,DEG
        NCBAS(I)=0
        NCBAS(I)=2*BSORDER(I)+2
      END DO

      CALL POLNCTETRA(POLORDER,MOLTYP,NCPOL)

      IF (NC.NE.(NCPOL+SUM(NCBAS))) THEN 
        WRITE(*,*) "PROBLEMS IN DEFINING PROPER NUMBER OF COEFFICIENTS"
        WRITE(*,*) "TOTAL NUMBER OF COEFFS ARE NOT SUMMING UP CORRECTLY"
        STOP
      END IF

      ALLOCATE(POL(NCPOL))
      ALLOCATE(DVDY1(NCPOL),DVDY2(NCPOL),DVDY3(NCPOL))
      ALLOCATE(DVDY4(NCPOL),DVDY5(NCPOL),DVDY6(NCPOL))

      DO I=1,NCPOL
        POL(I)=0.00D+00
        DVDY1(I)=0.00D+00
        DVDY2(I)=0.00D+00
        DVDY3(I)=0.00D+00
        DVDY4(I)=0.00D+00
        DVDY5(I)=0.00D+00
        DVDY6(I)=0.00D+00        
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

cat <<EOF >> CHIPR_TETRA_FUNC.f90
!######
! A4-TYPE
!######
        IF (DEG.EQ.1) THEN
          DO O=1,6 
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

cat <<EOF >> CHIPR_TETRA_FUNC.f90
!######
! AB3-TYPE
!######
        IF (DEG.EQ.2) THEN
          IF (I.EQ.1) THEN
!
!           A-B BASIS
!
            DO O=1,3
              CALL BASIS_CONTRACT(1,BSORDER(1),BS,R(O),Y(O))

!###### CALCULATING ANALYTIC DERIVATIVES OF Y WITH RESPECT TO R IF DER=.TRUE.
            
              IF (DER) THEN
                CALL DBASIS_CONTRACT(1,BSORDER(1),BS,R(O),DYDR(O))
              END IF
!######
            END DO
!###### 
          ELSE
!
!           B-B BASIS
!
            DO O=4,6 
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

cat <<EOF >> CHIPR_TETRA_FUNC.f90
!######
! A2B2-TYPE
!######
        IF (DEG.EQ.3) THEN
          IF (I.EQ.1) THEN
!
!           A-A BASIS
!
            O=1             
            CALL BASIS_CONTRACT(1,BSORDER(1),BS,R(O),Y(O)) 

!###### CALCULATING ANALYTIC DERIVATIVES OF Y WITH RESPECT TO R IF DER=.TRUE.
            
            IF (DER) THEN
              CALL DBASIS_CONTRACT(1,BSORDER(1),BS,R(O),DYDR(O))
            END IF
!######            
          ELSE IF (I.EQ.2) THEN
!
!           A-B BASIS
! 
            DO O=2,5 
              CALL BASIS_CONTRACT(2,BSORDER(2),BS,R(O),Y(O))

!###### CALCULATING ANALYTIC DERIVATIVES OF Y WITH RESPECT TO R IF DER=.TRUE.
            
              IF (DER) THEN
                CALL DBASIS_CONTRACT(2,BSORDER(2),BS,R(O),DYDR(O))
              END IF
!######
            END DO
!######          
          ELSE
!
!           B-B BASIS
!          
            O=6             
            CALL BASIS_CONTRACT(3,BSORDER(3),BS,R(O),Y(O)) 
           
!###### CALCULATING ANALYTIC DERIVATIVES OF Y WITH RESPECT TO R IF DER=.TRUE.
            
            IF (DER) THEN
              CALL DBASIS_CONTRACT(3,BSORDER(3),BS,R(O),DYDR(O))
            END IF
!######            
          END IF 
        END IF

EOF

else if ($DEG == 4) then 

cat <<EOF >> CHIPR_TETRA_FUNC.f90
!######
! ABC2-TYPE
!######
        IF (DEG.EQ.4) THEN
          IF (I.EQ.1) THEN
!
!           C-C BASIS
!
            O=1             
            CALL BASIS_CONTRACT(1,BSORDER(1),BS,R(O),Y(O))

!###### CALCULATING ANALYTIC DERIVATIVES OF Y WITH RESPECT TO R IF DER=.TRUE.
            
            IF (DER) THEN
              CALL DBASIS_CONTRACT(1,BSORDER(1),BS,R(O),DYDR(O))
            END IF
!######             
          ELSE IF (I.EQ.2) THEN 
!
!           C-A BASIS
!
            DO O=2,4,2 
              CALL BASIS_CONTRACT(2,BSORDER(2),BS,R(O),Y(O))

!###### CALCULATING ANALYTIC DERIVATIVES OF Y WITH RESPECT TO R IF DER=.TRUE.
            
              IF (DER) THEN
                CALL DBASIS_CONTRACT(2,BSORDER(2),BS,R(O),DYDR(O))
              END IF
!######
            END DO
            
         ELSE IF (I.EQ.3) THEN   
!
!           C-B BASIS
!
            DO O=3,5,2 
              CALL BASIS_CONTRACT(3,BSORDER(3),BS,R(O),Y(O))

!###### CALCULATING ANALYTIC DERIVATIVES OF Y WITH RESPECT TO R IF DER=.TRUE.
            
              IF (DER) THEN
                CALL DBASIS_CONTRACT(3,BSORDER(3),BS,R(O),DYDR(O))
              END IF
!######
            END DO
            
         ELSE 
!
!           A-B BASIS
!          
            O=6             
            CALL BASIS_CONTRACT(4,BSORDER(4),BS,R(O),Y(O))

!###### CALCULATING ANALYTIC DERIVATIVES OF Y WITH RESPECT TO R IF DER=.TRUE.
            
            IF (DER) THEN
              CALL DBASIS_CONTRACT(4,BSORDER(4),BS,R(O),DYDR(O))
            END IF
!######                        
          END IF 
        END IF

EOF

else if ($DEG == 6) then 

cat <<EOF >> CHIPR_TETRA_FUNC.f90
!######
! ABCD-TYPE
!######
        IF (DEG.EQ.6) THEN
          CALL BASIS_CONTRACT(I,BSORDER(I),BS,R(I),Y(I))

!###### CALCULATING ANALYTIC DERIVATIVES OF Y WITH RESPECT TO R IF DER=.TRUE.
            
          IF (DER) THEN
            CALL DBASIS_CONTRACT(I,BSORDER(I),BS,R(I),DYDR(I))
          END IF
!######
        END IF

EOF

endif

cat <<EOF >> CHIPR_TETRA_FUNC.f90
        SUMC=SUMC+NCBAS(I)
        DEALLOCATE(BS)
      END DO

      TOTNUM=0
      S=0      
      DO M=0,POLORDER
        MNUM=0
        DO I=0,M
          DO J=0,(M-I)
            DO K=0,(M-I-J)
              DO L=0,(M-I-J-K)
                DO N=0,(M-I-J-K-L)
                  O=(M-I-J-K-L-N)
                  S=(I+J+K+L+N+O)
!
! FIRST CONDITION: DO NOT INCLUDE TWO-BODY TERMS
!
                  IF (S.NE.I .AND. S.NE.J .AND. S.NE.K & 
     &               .AND. S.NE.L .AND. S.NE.N .AND. S.NE.O) THEN
!      
! SECOND CONDITION: DO NOT INCLUDE THREE-BODY TERMS
!
                  IF (S.NE.(I+J+K) .AND. S.NE.(I+J+L) .AND. S.NE.(I+J+N) .AND. & 
     &                S.NE.(I+J+O) .AND. S.NE.(J+K+L) .AND. S.NE.(J+K+N) .AND. &
     &                S.NE.(J+K+O) .AND. S.NE.(K+L+N) .AND. S.NE.(K+L+O) .AND. & 
     &                S.NE.(L+N+O) .AND. S.NE.(I+K+L) .AND. S.NE.(I+K+N) .AND. & 
     &                S.NE.(I+K+O) .AND. S.NE.(I+L+N) .AND. S.NE.(I+L+O) .AND. &
     &                S.NE.(I+N+O) .AND. S.NE.(J+L+N) .AND. S.NE.(J+L+O) .AND. &
     &                S.NE.(K+N+O) .AND. S.NE.(J+N+O)) THEN
            
!
!     FOR THE CASE OF AN ABCD MOLECULE 
!      
                   IF (MOLTYP.EQ."ABCD") THEN  

                     TOTNUM=TOTNUM+1
                     MNUM=MNUM+1  
                     
                     CALL PERMUTABCD(I,J,K,L,N,O,P,ID)

                     DO LL=1,ID
                     
                       POL(TOTNUM)=POL(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*Y(2)**(P(2,LL))*Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*Y(5)**(P(5,LL))*Y(6)**(P(6,LL)))  

!###### CALCULATING ANALYTIC DERIVATIVES OF V WITH RESPECT TO Y IF DER=.TRUE.     
     
                       IF (DER) THEN

                         DVDY1(TOTNUM)=DVDY1(TOTNUM)+&
     &                   (DBLE(P(1,LL))*Y(1)**(P(1,LL)-1)*&
     &                    Y(2)**(P(2,LL))*&
     &                    Y(3)**(P(3,LL))*&
     &                    Y(4)**(P(4,LL))*&
     &                    Y(5)**(P(5,LL))*&
     &                    Y(6)**(P(6,LL)))
                         
                          DVDY2(TOTNUM)=DVDY2(TOTNUM)+&
     &                   (Y(1)**(P(1,LL))*&
     &                    DBLE(P(2,LL))*Y(2)**(P(2,LL)-1)*&
     &                    Y(3)**(P(3,LL))*&
     &                    Y(4)**(P(4,LL))*&
     &                    Y(5)**(P(5,LL))*&
     &                    Y(6)**(P(6,LL)))
                         
                          DVDY3(TOTNUM)=DVDY3(TOTNUM)+&
     &                   (Y(1)**(P(1,LL))*&
     &                    Y(2)**(P(2,LL))*&
     &                    DBLE(P(3,LL))*Y(3)**(P(3,LL)-1)*&
     &                    Y(4)**(P(4,LL))*&
     &                    Y(5)**(P(5,LL))*&
     &                    Y(6)**(P(6,LL)))

                          DVDY4(TOTNUM)=DVDY4(TOTNUM)+&
     &                   (Y(1)**(P(1,LL))*&
     &                    Y(2)**(P(2,LL))*&
     &                    Y(3)**(P(3,LL))*&
     &                    DBLE(P(4,LL))*Y(4)**(P(4,LL)-1)*&
     &                    Y(5)**(P(5,LL))*&
     &                    Y(6)**(P(6,LL))) 
     
                          DVDY5(TOTNUM)=DVDY5(TOTNUM)+&
     &                   (Y(1)**(P(1,LL))*&
     &                    Y(2)**(P(2,LL))*&
     &                    Y(3)**(P(3,LL))*&
     &                    Y(4)**(P(4,LL))*&
     &                    DBLE(P(5,LL))*Y(5)**(P(5,LL)-1)*&
     &                    Y(6)**(P(6,LL)))
     
                          DVDY6(TOTNUM)=DVDY6(TOTNUM)+&
     &                   (Y(1)**(P(1,LL))*&
     &                    Y(2)**(P(2,LL))*&
     &                    Y(3)**(P(3,LL))*&
     &                    Y(4)**(P(4,LL))*&
     &                    Y(5)**(P(5,LL))*&
     &                    DBLE(P(6,LL))*Y(6)**(P(6,LL)-1))     
     
                       END IF                      

!######                                              
                     END DO

                     POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)
                     
!###### MULTIPLYING THE CONTRIBUTIONS OF DVDY BY THE POLYNOMIAL COEFFS. IF DER=.TRUE. 

                     IF (DER) THEN
                       DVDY1(TOTNUM)=C(TOTNUM)*DVDY1(TOTNUM)/DBLE(ID)
                       DVDY2(TOTNUM)=C(TOTNUM)*DVDY2(TOTNUM)/DBLE(ID)
                       DVDY3(TOTNUM)=C(TOTNUM)*DVDY3(TOTNUM)/DBLE(ID)
                       DVDY4(TOTNUM)=C(TOTNUM)*DVDY4(TOTNUM)/DBLE(ID)
                       DVDY5(TOTNUM)=C(TOTNUM)*DVDY5(TOTNUM)/DBLE(ID)
                       DVDY6(TOTNUM)=C(TOTNUM)*DVDY6(TOTNUM)/DBLE(ID)
                     END IF
!######                     
!
!     FOR THE CASE OF AN AB3 MOLECULE 
!
                   ELSE IF (MOLTYP.EQ."AB3") THEN

                     IF (I.LE.J .AND. L.LE.N .AND. N.LE.O .AND. J.LE.K) THEN

                       TOTNUM=TOTNUM+1
                       MNUM=MNUM+1

                       CALL PERMUTAB3(I,J,K,L,N,O,P,ID)
                       
                       DO LL=1,ID
                       
                         POL(TOTNUM)=POL(TOTNUM)+&                         
     &                      (Y(1)**(P(1,LL))*Y(2)**(P(2,LL))*Y(3)**(P(3,LL))*&
     &                       Y(4)**(P(4,LL))*Y(5)**(P(5,LL))*Y(6)**(P(6,LL))) 
     
!###### CALCULATING ANALYTIC DERIVATIVES OF V WITH RESPECT TO Y IF DER=.TRUE.     
     
                         IF (DER) THEN
                         
                           DVDY1(TOTNUM)=DVDY1(TOTNUM)+&
     &                     (DBLE(P(1,LL))*Y(1)**(P(1,LL)-1)*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL)))
                           
                            DVDY2(TOTNUM)=DVDY2(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      DBLE(P(2,LL))*Y(2)**(P(2,LL)-1)*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL)))
                           
                            DVDY3(TOTNUM)=DVDY3(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      DBLE(P(3,LL))*Y(3)**(P(3,LL)-1)*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL)))
                         
                            DVDY4(TOTNUM)=DVDY4(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      DBLE(P(4,LL))*Y(4)**(P(4,LL)-1)*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL))) 
                         
                            DVDY5(TOTNUM)=DVDY5(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      DBLE(P(5,LL))*Y(5)**(P(5,LL)-1)*&
     &                      Y(6)**(P(6,LL)))
                         
                            DVDY6(TOTNUM)=DVDY6(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      DBLE(P(6,LL))*Y(6)**(P(6,LL)-1))     
                         
                         END IF                      

!######                                                                     
                       END DO

                       POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

!###### MULTIPLYING THE CONTRIBUTIONS OF DVDY BY THE POLYNOMIAL COEFFS. IF DER=.TRUE. 

                       IF (DER) THEN
                         DVDY1(TOTNUM)=C(TOTNUM)*DVDY1(TOTNUM)/DBLE(ID)
                         DVDY2(TOTNUM)=C(TOTNUM)*DVDY2(TOTNUM)/DBLE(ID)
                         DVDY3(TOTNUM)=C(TOTNUM)*DVDY3(TOTNUM)/DBLE(ID)
                         DVDY4(TOTNUM)=C(TOTNUM)*DVDY4(TOTNUM)/DBLE(ID)
                         DVDY5(TOTNUM)=C(TOTNUM)*DVDY5(TOTNUM)/DBLE(ID)
                         DVDY6(TOTNUM)=C(TOTNUM)*DVDY6(TOTNUM)/DBLE(ID)
                       END IF
!######                        
                     END IF
!
!     FOR THE CASE OF AN A2B2 MOLECULE 
!
                   ELSE IF (MOLTYP.EQ."A2B2") THEN

                     IF (J.LE.K .AND. K.LE.L .AND. L.LE.N) THEN
                    
                       TOTNUM=TOTNUM+1
                       MNUM=MNUM+1     
                       
                       CALL PERMUTA2B2(I,J,K,L,N,O,P,ID)
                       
                       DO LL=1,ID
                         
                         POL(TOTNUM)=POL(TOTNUM)+&
     &                      (Y(1)**(P(1,LL))*Y(2)**(P(2,LL))*Y(3)**(P(3,LL))*&
     &                       Y(4)**(P(4,LL))*Y(5)**(P(5,LL))*Y(6)**(P(6,LL)))  
 
!###### CALCULATING ANALYTIC DERIVATIVES OF V WITH RESPECT TO Y IF DER=.TRUE.     
     
                         IF (DER) THEN
                         
                           DVDY1(TOTNUM)=DVDY1(TOTNUM)+&
     &                     (DBLE(P(1,LL))*Y(1)**(P(1,LL)-1)*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL)))
                           
                            DVDY2(TOTNUM)=DVDY2(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      DBLE(P(2,LL))*Y(2)**(P(2,LL)-1)*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL)))
                           
                            DVDY3(TOTNUM)=DVDY3(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      DBLE(P(3,LL))*Y(3)**(P(3,LL)-1)*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL)))
                         
                            DVDY4(TOTNUM)=DVDY4(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      DBLE(P(4,LL))*Y(4)**(P(4,LL)-1)*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL))) 
                         
                            DVDY5(TOTNUM)=DVDY5(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      DBLE(P(5,LL))*Y(5)**(P(5,LL)-1)*&
     &                      Y(6)**(P(6,LL)))
                         
                            DVDY6(TOTNUM)=DVDY6(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      DBLE(P(6,LL))*Y(6)**(P(6,LL)-1))     
                         
                         END IF                      

!######     
                       END DO                       
                      
                       POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID) 

!###### MULTIPLYING THE CONTRIBUTIONS OF DVDY BY THE POLYNOMIAL COEFFS. IF DER=.TRUE. 

                       IF (DER) THEN
                         DVDY1(TOTNUM)=C(TOTNUM)*DVDY1(TOTNUM)/DBLE(ID)
                         DVDY2(TOTNUM)=C(TOTNUM)*DVDY2(TOTNUM)/DBLE(ID)
                         DVDY3(TOTNUM)=C(TOTNUM)*DVDY3(TOTNUM)/DBLE(ID)
                         DVDY4(TOTNUM)=C(TOTNUM)*DVDY4(TOTNUM)/DBLE(ID)
                         DVDY5(TOTNUM)=C(TOTNUM)*DVDY5(TOTNUM)/DBLE(ID)
                         DVDY6(TOTNUM)=C(TOTNUM)*DVDY6(TOTNUM)/DBLE(ID)
                       END IF
!######                            
                     END IF                       
!
!     FOR THE CASE OF AN ABC2 MOLECULE 
!
                   ELSE IF (MOLTYP.EQ."ABC2") THEN

                     IF (J.LE.L .AND. K.LE.N) THEN

                       TOTNUM=TOTNUM+1
                       MNUM=MNUM+1     
                       
                       CALL PERMUTABC2(I,J,K,L,N,O,P,ID)
                       
                       DO LL=1,ID
                         
                         POL(TOTNUM)=POL(TOTNUM)+&
     &                      (Y(1)**(P(1,LL))*Y(2)**(P(2,LL))*Y(3)**(P(3,LL))*&
     &                       Y(4)**(P(4,LL))*Y(5)**(P(5,LL))*Y(6)**(P(6,LL)))  
 
!###### CALCULATING ANALYTIC DERIVATIVES OF V WITH RESPECT TO Y IF DER=.TRUE.     
     
                         IF (DER) THEN
                         
                           DVDY1(TOTNUM)=DVDY1(TOTNUM)+&
     &                     (DBLE(P(1,LL))*Y(1)**(P(1,LL)-1)*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL)))
                           
                            DVDY2(TOTNUM)=DVDY2(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      DBLE(P(2,LL))*Y(2)**(P(2,LL)-1)*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL)))
                           
                            DVDY3(TOTNUM)=DVDY3(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      DBLE(P(3,LL))*Y(3)**(P(3,LL)-1)*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL)))
                         
                            DVDY4(TOTNUM)=DVDY4(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      DBLE(P(4,LL))*Y(4)**(P(4,LL)-1)*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL))) 
                         
                            DVDY5(TOTNUM)=DVDY5(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      DBLE(P(5,LL))*Y(5)**(P(5,LL)-1)*&
     &                      Y(6)**(P(6,LL)))
                         
                            DVDY6(TOTNUM)=DVDY6(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      DBLE(P(6,LL))*Y(6)**(P(6,LL)-1))     
                         
                         END IF
                         
!######
                       END DO                       

                       POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)
                       
!###### MULTIPLYING THE CONTRIBUTIONS OF DVDY BY THE POLYNOMIAL COEFFS. IF DER=.TRUE. 

                       IF (DER) THEN
                         DVDY1(TOTNUM)=C(TOTNUM)*DVDY1(TOTNUM)/DBLE(ID)
                         DVDY2(TOTNUM)=C(TOTNUM)*DVDY2(TOTNUM)/DBLE(ID)
                         DVDY3(TOTNUM)=C(TOTNUM)*DVDY3(TOTNUM)/DBLE(ID)
                         DVDY4(TOTNUM)=C(TOTNUM)*DVDY4(TOTNUM)/DBLE(ID)
                         DVDY5(TOTNUM)=C(TOTNUM)*DVDY5(TOTNUM)/DBLE(ID)
                         DVDY6(TOTNUM)=C(TOTNUM)*DVDY6(TOTNUM)/DBLE(ID)
                       END IF
!######
                     END IF
!
!     FOR THE CASE OF AN A4 MOLECULE 
!
                   ELSE IF (MOLTYP.EQ."A4") THEN

                     IF (I.LE.J .AND. J.LE.K .AND. K.LE.L .AND. L.LE.N .AND. N.LE.O) THEN
                     
                       TOTNUM=TOTNUM+1
                       MNUM=MNUM+1     
                       
                       CALL PERMUTAA4(I,J,K,L,N,O,P,ID)
                       
                       DO LL=1,ID
                         
                         POL(TOTNUM)=POL(TOTNUM)+&
     &                      (Y(1)**(P(1,LL))*Y(2)**(P(2,LL))*Y(3)**(P(3,LL))*&
     &                       Y(4)**(P(4,LL))*Y(5)**(P(5,LL))*Y(6)**(P(6,LL)))  
 
!###### CALCULATING ANALYTIC DERIVATIVES OF V WITH RESPECT TO Y IF DER=.TRUE.     
     
                         IF (DER) THEN
                         
                           DVDY1(TOTNUM)=DVDY1(TOTNUM)+&
     &                     (DBLE(P(1,LL))*Y(1)**(P(1,LL)-1)*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL)))
                           
                            DVDY2(TOTNUM)=DVDY2(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      DBLE(P(2,LL))*Y(2)**(P(2,LL)-1)*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL)))
                           
                            DVDY3(TOTNUM)=DVDY3(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      DBLE(P(3,LL))*Y(3)**(P(3,LL)-1)*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL)))
                         
                            DVDY4(TOTNUM)=DVDY4(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      DBLE(P(4,LL))*Y(4)**(P(4,LL)-1)*&
     &                      Y(5)**(P(5,LL))*&
     &                      Y(6)**(P(6,LL))) 
                         
                            DVDY5(TOTNUM)=DVDY5(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      DBLE(P(5,LL))*Y(5)**(P(5,LL)-1)*&
     &                      Y(6)**(P(6,LL)))
                         
                            DVDY6(TOTNUM)=DVDY6(TOTNUM)+&
     &                     (Y(1)**(P(1,LL))*&
     &                      Y(2)**(P(2,LL))*&
     &                      Y(3)**(P(3,LL))*&
     &                      Y(4)**(P(4,LL))*&
     &                      Y(5)**(P(5,LL))*&
     &                      DBLE(P(6,LL))*Y(6)**(P(6,LL)-1))     
                         
                         END IF
                         
!######             
                       END DO                       

!###### MULTIPLYING THE CONTRIBUTIONS OF DVDY BY THE POLYNOMIAL COEFFS. IF DER=.TRUE. 

                       IF (DER) THEN
                         DVDY1(TOTNUM)=C(TOTNUM)*DVDY1(TOTNUM)/DBLE(ID)
                         DVDY2(TOTNUM)=C(TOTNUM)*DVDY2(TOTNUM)/DBLE(ID)
                         DVDY3(TOTNUM)=C(TOTNUM)*DVDY3(TOTNUM)/DBLE(ID)
                         DVDY4(TOTNUM)=C(TOTNUM)*DVDY4(TOTNUM)/DBLE(ID)
                         DVDY5(TOTNUM)=C(TOTNUM)*DVDY5(TOTNUM)/DBLE(ID)
                         DVDY6(TOTNUM)=C(TOTNUM)*DVDY6(TOTNUM)/DBLE(ID)
                       END IF
!######                     
                     END IF
                
                   END IF
!
!     END OF MOLECULE CASES
!                  
                  END IF
!                  
!     END OF SECOND CONDITION
!
                  END IF
!                  
!     END OF FIRST CONDITION
!
                END DO
              END DO
            END DO
          END DO
        END DO

      END DO
      
      POT=SUM(POL)*REPDAMPTETRA(NX,R)
      
!###### CALCULATING ANALYTIC DERIVATIVES OF V WITH RESPECT TO R IF DER=.TRUE.      
      
      IF (DER) THEN
        DVDR(1)=SUM(DVDY1)*DYDR(1)
        DVDR(2)=SUM(DVDY2)*DYDR(2)
        DVDR(3)=SUM(DVDY3)*DYDR(3)        
        DVDR(4)=SUM(DVDY4)*DYDR(4)
        DVDR(5)=SUM(DVDY5)*DYDR(5)
        DVDR(6)=SUM(DVDY6)*DYDR(6)
      ELSE 
        DVDR(1)=1000
        DVDR(2)=1000
        DVDR(3)=1000
        DVDR(4)=1000
        DVDR(5)=1000
        DVDR(6)=1000        
      END IF  
      
!######      
      
      DEALLOCATE(POL)
      DEALLOCATE(DVDY1,DVDY2,DVDY3)
      DEALLOCATE(DVDY4,DVDY5,DVDY6)

      RETURN
      END      

!####################################################################################

      SUBROUTINE CARTDERPES4BD(X,DVDX)
      INTEGER :: I
      INTEGER, PARAMETER :: NATOM=4     
      DOUBLE PRECISION, DIMENSION(3*NATOM) :: X
      DOUBLE PRECISION, DIMENSION(3*NATOM) :: DVDX,O
      DOUBLE PRECISION, DIMENSION(3*NATOM-6) :: Y,R,DVDR
      DOUBLE PRECISION :: V
!
!     This subprogram uses the chain rule to calculate the derivatives of the
!     energy with respect to the cartesian coordinates from the derivatives
!     with respect to the internal coordinates for a three-body system.
!     The convention assumed in this subprogram is as follows:
!
!     R(1) : R(A-B)
!     R(2) : R(A-C)
!     R(3) : R(A-D)
!     R(4) : R(B-C)
!     R(5) : R(B-D)
!     R(6) : R(C-D)
!     X(1)  - X(3)  : X, Y, Z for atom A
!     X(4)  - X(6)  : X, Y, Z for atom B
!     X(7)  - X(9)  : X, Y, Z for atom C
!     X(10) - X(12) : X, Y, Z for atom D
!
!
      CALL CART2INTER4BD(X,O,3*NATOM)

      R(1)=O(1)
      R(2)=O(2)
      R(3)=O(3)
      R(4)=O(4)
      R(5)=O(5)
      R(6)=O(6)

      CALL CHIPR_TETRA(R,V,.TRUE.,DVDR) 

      DO I=1,6
        Y(I)=DVDR(I)/R(I)
      END DO

      DVDX( 1)=(X(1)-X(4))*Y(1)+(X(1)-X(7))*Y(2)+(X(1)-X(10))*Y(3) 
      DVDX( 2)=(X(2)-X(5))*Y(1)+(X(2)-X(8))*Y(2)+(X(2)-X(11))*Y(3) 
      DVDX( 3)=(X(3)-X(6))*Y(1)+(X(3)-X(9))*Y(2)+(X(3)-X(12))*Y(3) 
      DVDX( 4)=(X(4)-X(1))*Y(1)+(X(4)-X(7))*Y(4)+(X(4)-X(10))*Y(5) 
      DVDX( 5)=(X(5)-X(2))*Y(1)+(X(5)-X(8))*Y(4)+(X(5)-X(11))*Y(5)
      DVDX( 6)=(X(6)-X(3))*Y(1)+(X(6)-X(9))*Y(4)+(X(6)-X(12))*Y(5)
      DVDX( 7)=(X(7)-X(1))*Y(2)+(X(7)-X(4))*Y(4)+(X(7)-X(10))*Y(6)
      DVDX( 8)=(X(8)-X(2))*Y(2)+(X(8)-X(5))*Y(4)+(X(8)-X(11))*Y(6)
      DVDX( 9)=(X(9)-X(3))*Y(2)+(X(9)-X(6))*Y(4)+(X(9)-X(12))*Y(6)
      DVDX(10)=(X(10)-X(1))*Y(3)+(X(10)-X(4))*Y(5)+(X(10)-X(7))*Y(6)
      DVDX(11)=(X(11)-X(2))*Y(3)+(X(11)-X(5))*Y(5)+(X(11)-X(8))*Y(6)
      DVDX(12)=(X(12)-X(3))*Y(3)+(X(12)-X(6))*Y(5)+(X(12)-X(9))*Y(6) 

      RETURN
      END

!####################################################################################

      SUBROUTINE CART2INTER4BD(X,R,NTA)
      IMPLICIT NONE
      INTEGER, PARAMETER :: NA=4,NT=3*NA
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

      DOUBLE PRECISION FUNCTION REPDAMPTETRA(NX,R)
      IMPLICIT NONE
      INTEGER :: I,NX
      DOUBLE PRECISION, DIMENSION(NX) :: R,H
      DOUBLE PRECISION :: KAPPA, XI, R0
      R0=0.5D+00
      KAPPA=100.0D+00
      XI=1.0D+00
      DO I=1,NX
        H(I)=(0.5D+00*(1.00D+00+TANH(KAPPA*(R(I)-R0))))
      END DO
      REPDAMPTETRA=(PRODUCT(H))**(XI)    
      END FUNCTION 

!################################################################################      
      
      SUBROUTINE PERMUTABCD(I,J,K,R,S,T,P,ID)
      IMPLICIT NONE
      INTEGER :: I,J,K,R,S,T
      INTEGER :: L,M,ID
      INTEGER, DIMENSION(6) :: INTER
      INTEGER, DIMENSION(6,24) :: P

      DO L=1,6 
        DO M=1,24
          P(L,M)=0
        END DO
      END DO

      INTER(1)=I
      INTER(2)=J
      INTER(3)=K
      INTER(4)=R
      INTER(5)=S
      INTER(6)=T      

      ID=1
        
!     ACTION OF THE INDENTITY OPERATOR ONLY

      P(1,1)=INTER(1)
      P(2,1)=INTER(2)
      P(3,1)=INTER(3) 
      P(4,1)=INTER(4) 
      P(5,1)=INTER(5)
      P(6,1)=INTER(6) 

      RETURN
      END

!################################################################################      
      
      SUBROUTINE PERMUTAA4(I,J,K,R,S,T,P,ID)
      IMPLICIT NONE
      INTEGER :: I,J,K,R,S,T
      INTEGER :: L,M,ID,O,U,V,X,SS
      INTEGER, DIMENSION(6) :: INTER
      INTEGER, DIMENSION(6,24) :: P,PINT,PRED   
      
      DO L=1,6 
        DO M=1,24
          P(L,M)=0
          PRED(L,M)=0
        END DO
      END DO
!
!     GENERATING ALL POSSIBLE PERMUTATIONS
!     THE REDUNDANT ONES ARE LATTER EXCLUDED
!
      INTER(1)=I
      INTER(2)=J
      INTER(3)=K
      INTER(4)=R
      INTER(5)=S
      INTER(6)=T
        
!     ACTION OF THE INDENTITY OPERATOR

      P(1,1)=INTER(1)
      P(2,1)=INTER(2)
      P(3,1)=INTER(3)
      P(4,1)=INTER(4)
      P(5,1)=INTER(5)
      P(6,1)=INTER(6)

!     ACTION OF THE (3,4) OPERATOR - R2<->R3 AND R4<->R5

      P(1,2)=INTER(1)
      P(2,2)=INTER(3)
      P(3,2)=INTER(2)
      P(4,2)=INTER(5)
      P(5,2)=INTER(4)
      P(6,2)=INTER(6) 
      
!     ACTION OF THE (2,3) OPERATOR - R1<->R2 AND R5<->R6      
      
      P(1,3)=INTER(2)
      P(2,3)=INTER(1)
      P(3,3)=INTER(3)
      P(4,3)=INTER(4)
      P(5,3)=INTER(6)
      P(6,3)=INTER(5)       
      
!     ACTION OF THE (2,3,4) OPERATOR - 123456 -> 231645      
      
      P(1,4)=INTER(2)
      P(2,4)=INTER(3)
      P(3,4)=INTER(1)
      P(4,4)=INTER(6)
      P(5,4)=INTER(4)
      P(6,4)=INTER(5)       
      
!     ACTION OF THE (2,4,3) OPERATOR - 123456 -> 312564 
      
      P(1,5)=INTER(3)      
      P(2,5)=INTER(1)      
      P(3,5)=INTER(2)      
      P(4,5)=INTER(5)      
      P(5,5)=INTER(6)      
      P(6,5)=INTER(4)       
      
!     ACTION OF THE (2,4) OPERATOR - R1<->R3 AND R4<->R6      
      
      P(1,6)=INTER(3)      
      P(2,6)=INTER(2)      
      P(3,6)=INTER(1)      
      P(4,6)=INTER(6)      
      P(5,6)=INTER(5)      
      P(6,6)=INTER(4)    
   
!     ACTION OF THE (1,2) OPERATOR - R3<->R5 AND R4<->R2  

      P(1,7)=INTER(1)       
      P(2,7)=INTER(4)       
      P(3,7)=INTER(5)    
      P(4,7)=INTER(2)    
      P(5,7)=INTER(3)       
      P(6,7)=INTER(6) 
      
!     ACTION OF THE (1,2)(3,4) OPERATOR - R3<->R4 AND R2<->R5  

      P(1,8)=INTER(1)       
      P(2,8)=INTER(5)       
      P(3,8)=INTER(4)    
      P(4,8)=INTER(3)    
      P(5,8)=INTER(2)       
      P(6,8)=INTER(6)    
      
!     ACTION OF THE (1,2,3) OPERATOR - 123456->415263  

      P(1,9)=INTER(4)       
      P(2,9)=INTER(1)       
      P(3,9)=INTER(5)    
      P(4,9)=INTER(2)    
      P(5,9)=INTER(6)       
      P(6,9)=INTER(3)
      
!     ACTION OF THE (1,2,3,4) OPERATOR - 123456->451623  

      P(1,10)=INTER(4)       
      P(2,10)=INTER(5)       
      P(3,10)=INTER(1)    
      P(4,10)=INTER(6)    
      P(5,10)=INTER(2)       
      P(6,10)=INTER(3)  
      
!     ACTION OF THE (1,2,4,3) OPERATOR - 123456->514362  

      P(1,11)=INTER(5)       
      P(2,11)=INTER(1)       
      P(3,11)=INTER(4)    
      P(4,11)=INTER(3)    
      P(5,11)=INTER(6)       
      P(6,11)=INTER(2)  
      
!     ACTION OF THE (1,2,4) OPERATOR - 123456->541632  

      P(1,12)=INTER(5)       
      P(2,12)=INTER(4)       
      P(3,12)=INTER(1)    
      P(4,12)=INTER(6)    
      P(5,12)=INTER(3)       
      P(6,12)=INTER(2)       
      
!     ACTION OF THE (1,3,2) OPERATOR - 123456->246135

      P(1,13)=INTER(2)       
      P(2,13)=INTER(4)       
      P(3,13)=INTER(6)    
      P(4,13)=INTER(1)    
      P(5,13)=INTER(3)       
      P(6,13)=INTER(5)  
      
!     ACTION OF THE (1,3,4,2) OPERATOR - 123456->264315

      P(1,14)=INTER(2)       
      P(2,14)=INTER(6)       
      P(3,14)=INTER(4)    
      P(4,14)=INTER(3)    
      P(5,14)=INTER(1)       
      P(6,14)=INTER(5)        
      
!     ACTION OF THE (1,3) OPERATOR - R1<->R4 AND R3<->R6

      P(1,15)=INTER(4)       
      P(2,15)=INTER(2)       
      P(3,15)=INTER(6)    
      P(4,15)=INTER(1)    
      P(5,15)=INTER(5)       
      P(6,15)=INTER(3)       
      
!     ACTION OF THE (1,3,4) OPERATOR - 123456->462513

      P(1,16)=INTER(4)       
      P(2,16)=INTER(6)       
      P(3,16)=INTER(2)    
      P(4,16)=INTER(5)    
      P(5,16)=INTER(1)       
      P(6,16)=INTER(3)
      
!     ACTION OF THE (1,3)(2,4) OPERATOR - R1<->R6 AND R3<->R4

      P(1,17)=INTER(6)       
      P(2,17)=INTER(2)       
      P(3,17)=INTER(4)    
      P(4,17)=INTER(3)    
      P(5,17)=INTER(5)       
      P(6,17)=INTER(1)      
      
!     ACTION OF THE (1,3,2,4) OPERATOR - 123456->642531

      P(1,18)=INTER(6)       
      P(2,18)=INTER(4)       
      P(3,18)=INTER(2)    
      P(4,18)=INTER(5)    
      P(5,18)=INTER(3)       
      P(6,18)=INTER(1)       
      
!     ACTION OF THE (1,4,3,2) OPERATOR - 123456->356124

      P(1,19)=INTER(3)       
      P(2,19)=INTER(5)       
      P(3,19)=INTER(6)    
      P(4,19)=INTER(1)    
      P(5,19)=INTER(2)       
      P(6,19)=INTER(4)     
      
!     ACTION OF THE (1,4,2) OPERATOR - 123456->365214

      P(1,20)=INTER(3)       
      P(2,20)=INTER(6)       
      P(3,20)=INTER(5)    
      P(4,20)=INTER(2)    
      P(5,20)=INTER(1)       
      P(6,20)=INTER(4)  
      
!     ACTION OF THE (1,4,3) OPERATOR - 123456->536142

      P(1,21)=INTER(5)       
      P(2,21)=INTER(3)       
      P(3,21)=INTER(6)    
      P(4,21)=INTER(1)    
      P(5,21)=INTER(4)       
      P(6,21)=INTER(2)  
      
!     ACTION OF THE (1,4) OPERATOR - R2<->R6 AND R1<->R5

      P(1,22)=INTER(5)       
      P(2,22)=INTER(6)       
      P(3,22)=INTER(3)    
      P(4,22)=INTER(4)    
      P(5,22)=INTER(1)       
      P(6,22)=INTER(2)         
      
!     ACTION OF THE (1,4,2,3) OPERATOR - 123456->635241

      P(1,23)=INTER(6)       
      P(2,23)=INTER(3)       
      P(3,23)=INTER(5)    
      P(4,23)=INTER(2)    
      P(5,23)=INTER(4)       
      P(6,23)=INTER(1)  
      
!     ACTION OF THE (1,4)(2,3) OPERATOR - R1<->R6 AND R2<->R5

      P(1,24)=INTER(6)       
      P(2,24)=INTER(5)       
      P(3,24)=INTER(3)    
      P(4,24)=INTER(4)    
      P(5,24)=INTER(2)       
      P(6,24)=INTER(1)                   
                            
      DO L=1,6 
        DO M=1,24
          PINT(L,M)=P(L,M)
        END DO
      END DO
!
!     EXCLUDING REDUNDANT EXCITATIONS
!
!
!     FIRST, MARK THEN WITH THE NUMBER 1000
! 
      DO L=1,24
        DO M=(L+1),24
          SS=0
          DO U=1,6  
            IF (P(U,L).EQ.P(U,M)) THEN
            SS=SS+U
              IF (SS.EQ.21) THEN                  
                DO O=1,6
                  PINT(O,M)=1000
                END DO  
              END IF
            END IF 
          END DO  
        END DO
      END DO
!
!     SECOND, EXCLUDE THOSE MARKED FROM THE VECTOR
!      
      DO L=1,6 
        V=0
        DO M=1,24
          IF (PINT(L,M).NE.1000) THEN
            V=V+1
            PRED(L,V)=PINT(L,M)
          END IF          
        END DO
        ID=V
      END DO      
!
!     DEFINE A THE NEW PERMUTATION OPERATORS
!
      DO L=1,6 
        DO M=1,24
          P(L,M)=0
        END DO
      END DO
      
      DO L=1,6 
        DO M=1,24
          P(L,M)=PRED(L,M)
        END DO
      END DO

      RETURN 
      END

!################################################################################      
      
      SUBROUTINE PERMUTABC2(I,J,K,R,S,T,P,ID)
      IMPLICIT NONE
      INTEGER :: I,J,K,R,S,T
      INTEGER :: L,M,ID,O,U,V,X,SS
      INTEGER, DIMENSION(6) :: INTER
      INTEGER, DIMENSION(6,24) :: P,PINT,PRED   
      
      DO L=1,6 
        DO M=1,24
          P(L,M)=0
          PRED(L,M)=0
        END DO
      END DO
!
!     GENERATING ALL POSSIBLE PERMUTATIONS
!     THE REDUNDANT ONES ARE LATTER EXCLUDED
!
      INTER(1)=I
      INTER(2)=J
      INTER(3)=K
      INTER(4)=R
      INTER(5)=S
      INTER(6)=T
        
!     ACTION OF THE INDENTITY OPERATOR

      P(1,1)=INTER(1)
      P(2,1)=INTER(2)
      P(3,1)=INTER(3)
      P(4,1)=INTER(4)
      P(5,1)=INTER(5)
      P(6,1)=INTER(6)

!     ACTION OF THE (1,2) OPERATOR

      P(1,2)=INTER(1)
      P(2,2)=INTER(4)
      P(3,2)=INTER(5)
      P(4,2)=INTER(2)
      P(5,2)=INTER(3)
      P(6,2)=INTER(6) 
                            
      DO L=1,6 
        DO M=1,24
          PINT(L,M)=P(L,M)
        END DO
      END DO
!
!     EXCLUDING REDUNDANT EXCITATIONS
!
!
!     FIRST, MARK THEN WITH THE NUMBER 1000
!
      DO L=1,2
        DO M=(L+1),2
          SS=0
          DO U=2,5
            IF (P(U,L).EQ.P(U,M)) THEN
              SS=SS+U              
              IF (SS.EQ.14) THEN                  
                DO O=1,6
                  PINT(O,M)=1000
                END DO  
              END IF
            END IF 
          END DO  
        END DO
      END DO  
!
!     SECOND, EXCLUDE THOSE MARKED FROM THE VECTOR
!      
      DO L=1,6 
        V=0
        DO M=1,2
          IF (PINT(L,M).NE.1000) THEN
            V=V+1
            PRED(L,V)=PINT(L,M)
          END IF          
        END DO
        ID=V
      END DO      
!
!     DEFINE A THE NEW PERMUTATION OPERATORS
!
      DO L=1,6 
        DO M=1,24
          P(L,M)=0
        END DO
      END DO
      
      DO L=1,6 
        DO M=1,24
          P(L,M)=PRED(L,M)
        END DO
      END DO        

      RETURN 
      END

!################################################################################ 

      SUBROUTINE PERMUTA2B2(I,J,K,R,S,T,P,ID)
      IMPLICIT NONE
      INTEGER :: I,J,K,R,S,T
      INTEGER :: L,M,ID,O,U,V,X,SS
      INTEGER, DIMENSION(6) :: INTER
      INTEGER, DIMENSION(6,24) :: P,PINT,PRED   
      
      DO L=1,6 
        DO M=1,24
          P(L,M)=0
          PRED(L,M)=0
        END DO
      END DO
!
!     GENERATING ALL POSSIBLE PERMUTATIONS
!     THE REDUNDANT ONES ARE LATTER EXCLUDED
!
      INTER(1)=I
      INTER(2)=J
      INTER(3)=K
      INTER(4)=R
      INTER(5)=S
      INTER(6)=T
        
!     ACTION OF THE INDENTITY OPERATOR

      P(1,1)=INTER(1)
      P(2,1)=INTER(2)
      P(3,1)=INTER(3)
      P(4,1)=INTER(4)
      P(5,1)=INTER(5)
      P(6,1)=INTER(6)

!     ACTION OF THE (3,4) OPERATOR

      P(1,2)=INTER(1)
      P(2,2)=INTER(3)
      P(3,2)=INTER(2)
      P(4,2)=INTER(5)
      P(5,2)=INTER(4)
      P(6,2)=INTER(6) 
        
!     ACTION OF THE (1,2) OPERATOR

      P(1,3)=INTER(1)
      P(2,3)=INTER(4)
      P(3,3)=INTER(5)
      P(4,3)=INTER(2)
      P(5,3)=INTER(3)
      P(6,3)=INTER(6)        
        
!     ACTION OF THE (1,2)(3,4) OPERATOR

      P(1,4)=INTER(1)
      P(2,4)=INTER(5)
      P(3,4)=INTER(4)
      P(4,4)=INTER(3)
      P(5,4)=INTER(2)
      P(6,4)=INTER(6)  
                    
      DO L=1,6 
        DO M=1,24
          PINT(L,M)=P(L,M)
        END DO
      END DO
!
!     EXCLUDING REDUNDANT EXCITATIONS
!
!
!     FIRST, MARK THEN WITH THE NUMBER 1000
!
      DO L=1,4
        DO M=(L+1),4
          SS=0
          DO U=2,5
            IF (P(U,L).EQ.P(U,M)) THEN
              SS=SS+U              
              IF (SS.EQ.14) THEN                  
                DO O=1,6
                  PINT(O,M)=1000
                END DO  
              END IF
            END IF 
          END DO  
        END DO
      END DO  
!
!     SECOND, EXCLUDE THOSE MARKED FROM THE VECTOR
!      
      DO L=1,6 
        V=0
        DO M=1,4
          IF (PINT(L,M).NE.1000) THEN
            V=V+1
            PRED(L,V)=PINT(L,M)
          END IF          
        END DO
        ID=V
      END DO      
!
!     DEFINE A THE NEW PERMUTATION OPERATORS
!
      DO L=1,6 
        DO M=1,24
          P(L,M)=0
        END DO
      END DO
      
      DO L=1,6 
        DO M=1,24
          P(L,M)=PRED(L,M)
        END DO
      END DO        

      RETURN 
      END

!################################################################################ 

      SUBROUTINE PERMUTAB3(I,J,K,R,S,T,P,ID)
      IMPLICIT NONE
      INTEGER :: I,J,K,R,S,T
      INTEGER :: L,M,ID,O,U,V,X,SS
      INTEGER, DIMENSION(6) :: INTER
      INTEGER, DIMENSION(6,24) :: P,PINT,PRED      
      
      DO L=1,6 
        DO M=1,24
          P(L,M)=0
          PRED(L,M)=0
        END DO
      END DO
!
!     GENERATING ALL POSSIBLE PERMUTATIONS
!     THE REDUNDANT ONES ARE LATTER EXCLUDED
!
      INTER(1)=I
      INTER(2)=J
      INTER(3)=K
      INTER(4)=R
      INTER(5)=S
      INTER(6)=T
        
!     ACTION OF THE INDENTITY OPERATOR

      P(1,1)=INTER(1)
      P(2,1)=INTER(2)
      P(3,1)=INTER(3)
      P(4,1)=INTER(4)
      P(5,1)=INTER(5)
      P(6,1)=INTER(6)

!     ACTION OF THE (3,4) OPERATOR

      P(1,2)=INTER(1)
      P(2,2)=INTER(3)
      P(3,2)=INTER(2)
      P(4,2)=INTER(5)
      P(5,2)=INTER(4)
      P(6,2)=INTER(6) 
        
!     ACTION OF THE (2,3) OPERATOR

      P(1,3)=INTER(2)
      P(2,3)=INTER(1)
      P(3,3)=INTER(3)
      P(4,3)=INTER(4)
      P(5,3)=INTER(6)
      P(6,3)=INTER(5)        
        
!     ACTION OF THE (2,4) OPERATOR

      P(1,4)=INTER(3)
      P(2,4)=INTER(2)
      P(3,4)=INTER(1)
      P(4,4)=INTER(6)
      P(5,4)=INTER(5)
      P(6,4)=INTER(4)  
        
!     ACTION OF THE (2,3,4) OPERATOR

      P(1,5)=INTER(2)
      P(2,5)=INTER(3)
      P(3,5)=INTER(1)
      P(4,5)=INTER(6)
      P(5,5)=INTER(4)
      P(6,5)=INTER(5)
        
!     ACTION OF THE (2,4,3) OPERATOR

      P(1,6)=INTER(3)
      P(2,6)=INTER(1)
      P(3,6)=INTER(2)
      P(4,6)=INTER(5)
      P(5,6)=INTER(6)
      P(6,6)=INTER(4)                

      DO L=1,6 
        DO M=1,24
          PINT(L,M)=P(L,M)
        END DO
      END DO
!
!     EXCLUDING REDUNDANT EXCITATIONS
!
!
!     FIRST, MARK THEN WITH THE NUMBER 1000
!
      DO L=1,6
        DO M=(L+1),6
          SS=0
          DO U=1,6  
            IF (P(U,L).EQ.P(U,M)) THEN
            SS=SS+U
              IF (SS.EQ.21) THEN                  
                DO O=1,6
                  PINT(O,M)=1000
                END DO  
              END IF
            END IF 
          END DO  
        END DO
      END DO  
!
!     SECOND, EXCLUDE THOSE MARKED FROM THE VECTOR
!      
      DO L=1,6 
        V=0
        DO M=1,6
          IF (PINT(L,M).NE.1000) THEN
            V=V+1
            PRED(L,V)=PINT(L,M)
          END IF          
        END DO
        ID=V
      END DO      
!
!     DEFINE A THE NEW PERMUTATION OPERATORS
!
      DO L=1,6 
        DO M=1,24
          P(L,M)=0
        END DO
      END DO
      
      DO L=1,6 
        DO M=1,24
          P(L,M)=PRED(L,M)
        END DO
      END DO        

      RETURN 
      END

!################################################################################ 

      SUBROUTINE POLNCTETRA(ORDER,MOLTYP,NC)
      IMPLICIT NONE
      INTEGER :: I,J,K,L,M,N,O,Q,R,S,LL   
      INTEGER :: NC,MNUM,ORDER
      CHARACTER(LEN=4) :: MOLTYP

      NC=0
      S=0

      DO M=0,ORDER
        MNUM=0
        DO I=0,M
          DO J=0,(M-I)
            DO K=0,(M-I-J)
              DO L=0,(M-I-J-K)
                DO N=0,(M-I-J-K-L)
                  O=(M-I-J-K-L-N)
                  S=(I+J+K+L+N+O)
!
! FIRST CONDITION: DO NOT INCLUDE TWO-BODY TERMS
!
                  IF (S.NE.I .AND. S.NE.J .AND. S.NE.K & 
      &              .AND. S.NE.L .AND. S.NE.N .AND. S.NE.O) THEN
!      
! SECOND CONDITION: DO NOT INCLUDE THREE-BODY TERMS
!
                  IF (S.NE.(I+J+K) .AND. S.NE.(I+J+L) .AND. S.NE.(I+J+N) .AND. & 
      &               S.NE.(I+J+O) .AND. S.NE.(J+K+L) .AND. S.NE.(J+K+N) .AND. &
      &               S.NE.(J+K+O) .AND. S.NE.(K+L+N) .AND. S.NE.(K+L+O) .AND. & 
      &               S.NE.(L+N+O) .AND. S.NE.(I+K+L) .AND. S.NE.(I+K+N) .AND. & 
      &               S.NE.(I+K+O) .AND. S.NE.(I+L+N) .AND. S.NE.(I+L+O) .AND. &
      &               S.NE.(I+N+O) .AND. S.NE.(J+L+N) .AND. S.NE.(J+L+O) .AND. &
      &               S.NE.(K+N+O) .AND. S.NE.(J+N+O)) THEN
!
!     FOR THE CASE OF AN ABCD MOLECULE 
!      
                   IF (MOLTYP.EQ."ABCD") THEN  

                     NC=NC+1
                     MNUM=MNUM+1                      
!
!     FOR THE CASE OF AN AB3 MOLECULE 
!
                   ELSE IF (MOLTYP.EQ."AB3") THEN

                     IF (I.LE.J .AND. L.LE.N .AND. N.LE.O .AND. J.LE.K) THEN

                       NC=NC+1
                       MNUM=MNUM+1

                     END IF
!
!     FOR THE CASE OF AN A2B2 MOLECULE 
!
                   ELSE IF (MOLTYP.EQ."A2B2") THEN

                     IF (J.LE.K .AND. K.LE.L .AND. L.LE.N) THEN
                     
                       NC=NC+1
                       MNUM=MNUM+1                     
                                        
                     END IF
!
!     FOR THE CASE OF AN ABC2 MOLECULE 
!
                   ELSE IF (MOLTYP.EQ."ABC2") THEN

                     IF (J.LE.L .AND. K.LE.N) THEN
                     
                       NC=NC+1
                       MNUM=MNUM+1                     
                     
                     END IF
!
!     FOR THE CASE OF AN A4 MOLECULE 
!
                   ELSE IF (MOLTYP.EQ."A4") THEN

                     IF (I.LE.J .AND. J.LE.K .AND. K.LE.L .AND. L.LE.N .AND. N.LE.O) THEN
                     
                       NC=NC+1
                       MNUM=MNUM+1
                     
                     END IF
                
                   END IF
                  
                  END IF
!                  
! END OF SECOND CONDITION
!
                  END IF
!                  
! END OF FIRST CONDITION
!
                END DO
              END DO
            END DO
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
