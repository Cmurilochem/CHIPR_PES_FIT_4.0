!######################################################################################################
! CHIPR MODEL TRIATOMIC POTENTIAL
!######################################################################################################
      SUBROUTINE CHIPR_TRIAT(BSORDER,POLORDER,C,R,POT,DIFFC)
!######################################################################################################
! DEG: NUMBER OF DIFFERENT SYMMETRY UNRELATED DEGREES OF FREEDON: 
! DEG=1 FOR AN A3, DEG=2 FOR AN AB2, AND DEG=3 FOR AN ABC-TYPE MOLECULE
! BSORDER: IS A VECTOR OF DIMENSION "DEG" EXPECIFYING THE ORDER OF EACH BASIS SET
! POLORDER: DEFINES THE ORDER OF THE POLYNOMIAL
! C IS THE VECTOR OF COEFFICIENTS
! R: IS THE 3D VECTOR OF DISTANCES IN A.U. 
! POT: THE POTENTIAL IN A.U.
! DIFFC: DEVIVATIVE OF THE POLINOMIAL WITH RESPECT TO THE COEFFICIENTS TO BE USED IN LMDIF 
!######################################################################################################
      USE COMMON_VAR, ONLY : MOLTYP,DEG,NC,NX
      IMPLICIT NONE
      INTEGER :: I,J,K,L,M,S,O,ID
! POLORDER: DEFINES THE ORDER OF THE POLYNOMIAL
      INTEGER :: POLORDER
! BSORDER: DEFINES THE ORDER OF EACH BASIS SET
      INTEGER, DIMENSION(DEG) :: BSORDER,NCBAS
      DOUBLE PRECISION, DIMENSION(NC) :: C,DIFFC
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: POL,BS
      DOUBLE PRECISION, DIMENSION(3) :: R
      DOUBLE PRECISION :: POT,REPDAMP,F1,F2,F3,F4
      INTEGER :: NCPOL,NCTOTAL,SUMC
      DOUBLE PRECISION, DIMENSION(NX) :: Y
      INTEGER :: TOTNUM,MNUM
      INTEGER, DIMENSION(3,6) :: P
!
!     CALCULATING THE NUMBER OF COEFFICIENTS FOR EACH FOR EACH BASIS SET
!
      DO I=1,DEG
        NCBAS(I)=0
        NCBAS(I)=2*BSORDER(I)+2
      END DO
      
      F1=1.5D+00!0.65!0.56D+00!0.63836341849131550407747681674663909D+00 ! 12
      F2=1.5D+00!2.85!1.975D+00!0.28983388620287229286986985243856907D+01 ! 12
      F3=0.9D+00!0.525!0.565D+00!0.47992291126355957064930635169730522D+00 ! 12
      F4=1.0D+00!1.5!1.4D+00!0.27397796376384699890138563205255195D+01 ! 12
!
!     CALCULATE THE NUMBER OF COEFFICIENTS FOR EACH MOLECULE TYPE 
!     AND POLYNOMIALS ORDER
!
      CALL POLNC(POLORDER,MOLTYP,NCPOL)
!
!     CHECKING - JUST IN CASE
!
      IF (NC.NE.(NCPOL+SUM(NCBAS))) THEN 
        WRITE(6,*) "PROBLEMS IN DEFINING PROPER NUMBER OF COEFFICIENTS"
        WRITE(6,*) "TOTAL NUMBER OF COEFFS ARE NOT SUMMING UP CORRECTLY"
        STOP
      END IF
!
!     ALLOCATING REQUIRED VECTORS
!
      ALLOCATE(POL(NCPOL))
!
!     INICIALIZING VECTORS
!
      DO I=1,NCPOL
        POL(I)=0.00D+00
      END DO

      DO I=1,NX
        Y(I)=0.0D+00
      END DO

      DO I=1,NC
        DIFFC(I)=0.00D+00
      END DO
!
!     HERE, CALCULATE THE VALUES OF Y'S IN THE BASIS SET
!
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
!
!     IN THE CASE OF AN A3-TYPE MOLECULE
!
        IF (DEG.EQ.1) THEN
          DO O=1,3 
            CALL BASIS_CONTRACT(1,BSORDER(1),BS,R(O),Y(O))
          END DO
        END IF
!
!     IN THE CASE OF AN AB2-TYPE MOLECULE
!
        IF (DEG.EQ.2) THEN
          IF (I.EQ.1) THEN
!
!           B-B BASIS
!
             CALL BASIS_CONTRACT(1,BSORDER(1),BS,R(1),Y(1)) 
             !Y(1)=(0.5D+00*(1.0D+00-TANH(0.55D+00*(R(1)-2.5D+00))))
             !Y(1)=R(1)/(1.0D+00+EXP(F1*(R(1)-F2)))!1.0D+00/(1.0D+00+EXP(F1*(R(1)-F2))) 
          ELSE
            DO O=2,3
!
!           A-B BASIS
! 
              CALL BASIS_CONTRACT(2,BSORDER(2),BS,R(O),Y(O))
              !Y(O)=(0.5D+00*(1.0D+00-TANH(0.25D+00*(R(O)-1.5D+00))))
              !Y(O)=R(O)/(1.0D+00+EXP(F3*(R(O)-F4)))!1.0D+00/(1.0D+00+EXP(F3*(R(O)-F4)))
            END DO
          END IF            
        END IF
!
!     IN THE CASE OF AN ABC-TYPE MOLECULE
!
        IF (DEG.EQ.3) THEN
          CALL BASIS_CONTRACT(I,BSORDER(I),BS,R(I),Y(I))
        END IF
!
!     END OF BASIS SETS EVALUATION
!
        SUMC=SUMC+NCBAS(I)
        DEALLOCATE(BS)
      END DO
!
!     EVALUATING THE POLYNOMIALS
!
      TOTNUM=0
      S=0
      DO M=0,POLORDER
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

                TOTNUM=TOTNUM+1
                MNUM=MNUM+1

                CALL PERMUTABC(I,J,K,P,ID)

                DO L=1,ID

                  POL(TOTNUM)=POL(TOTNUM)+&
     &               (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

!                  WRITE(*,*) P(1,L),P(2,L),P(3,L),TOTNUM,MNUM,POL(TOTNUM)

                END DO
!
!     NORMALIZING COEFFICIENT
!
                POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)
!
!     DEVIVATIVE OF THE POLINOMIAL WITH RESPECT TO C(TOTNUM) 
!
                     DIFFC(TOTNUM)=POL(TOTNUM) 
!
!     FOR THE CASE OF AN AB2 MOLECULE 
!
              ELSE IF (MOLTYP.EQ."AB2") THEN 

                IF (J.LE.K) THEN 

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTAB2(I,J,K,P,ID)
                   
                  DO L=1,ID

                    POL(TOTNUM)=POL(TOTNUM)+&
     &                 (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

!                    WRITE(*,*) P(1,L),P(2,L),P(3,L),TOTNUM,MNUM,POL(TOTNUM)

                  END DO
!
!     NORMALIZING COEFFICIENT
!
                  POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)
!
!     DEVIVATIVE OF THE POLINOMIAL WITH RESPECT TO C(TOTNUM) 
!
                  DIFFC(TOTNUM)=POL(TOTNUM) 

                END IF    
!
!     FOR THE CASE OF AN A3 MOLECULE 
!         
              ELSE IF (MOLTYP.EQ."A3") THEN

                IF (I.LE.J .AND. J.LE.K) THEN

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTA3(I,J,K,P,ID)
                   
                  DO L=1,ID

                    POL(TOTNUM)=POL(TOTNUM)+&
     &                 (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

!                    WRITE(*,*) P(1,L),P(2,L),P(3,L),TOTNUM,MNUM,POL(TOTNUM)

                  END DO
!
!     NORMALIZING COEFFICIENT
!
                  POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)
!
!     DEVIVATIVE OF THE POLINOMIAL WITH RESPECT TO C(TOTNUM) 
!
                  DIFFC(TOTNUM)=POL(TOTNUM) 

                END IF
!
!     END OF CASES
!
              END IF

            END IF
          END DO
        END DO

      END DO

      POT=SUM(POL)*REPDAMP(NX,R)

      DEALLOCATE(POL)

      RETURN
      END

!################################################################################
! DAMPING THE 3-BODY ENERGY TERM AT SHORT INTERNUCLEAR DISTANCES
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
        H(I)=(0.5D+00*(1.00D+00+TANH(KAPPA*(R(I)-R0))))
      END DO
      REPDAMP=(PRODUCT(H))**(XI)      
      END FUNCTION

!################################################################################
! FINDING THE POSSIBLE PERMUTATIONS OF I J K FOR AN ABC-TYPE MOLECULE
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
        
!     ACTION OF THE INDENTITY OPERATOR ONLY

      P(1,1)=INTER(1)
      P(2,1)=INTER(2)
      P(3,1)=INTER(3) 

      RETURN
      END

!################################################################################
! FINDING THE POSSIBLE PERMUTATIONS OF I J K FOR AN AB2-TYPE MOLECULE
! Y1 IS ALWAYS THE B-B BOND
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

!     CALCULATE THE INDENTITY AND (2,3) OPERATORS ONLY 
      
        ID=2
      
!     ACTION OF THE INDENTITY OPERATOR

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)      

!     ACTION OF THE (2,3) OPERATOR

        P(1,2)=INTER(1)
        P(2,2)=INTER(3)
        P(3,2)=INTER(2) 

       ELSE 

        ID=1

!     ACTION OF THE INDENTITY OPERATOR ONLY

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)

       END IF

      RETURN 
      END 

!################################################################################
! FINDING THE POSSIBLE PERMUTATIONS OF I J K FOR AN A3-TYPE MOLECULE
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
        
!     ACTION OF THE INDENTITY OPERATOR ONLY

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)  
        
      ELSE IF (I.NE.J.AND.J.NE.K.AND.I.NE.K) THEN   

!     CALCULATE THE FULL S3 SYMMETRIC GROUP 
      
        ID=6
      
!     ACTION OF THE INDENTITY OPERATOR

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)      

!     ACTION OF THE (1,2) OPERATOR

        P(1,2)=INTER(2)
        P(2,2)=INTER(1)
        P(3,2)=INTER(3)  

!     ACTION OF THE (1,3) OPERATOR

        P(1,3)=INTER(3)
        P(2,3)=INTER(2)
        P(3,3)=INTER(1) 

!     ACTION OF THE (2,3) OPERATOR

        P(1,4)=INTER(1)
        P(2,4)=INTER(3)
        P(3,4)=INTER(2) 

!     ACTION OF THE (1,2,3) OPERATOR   

        P(1,5)=P(1,2)
        P(2,5)=P(3,2)
        P(3,5)=P(2,2)  
      
!     ACTION OF THE (1,3,2) OPERATOR 

        P(1,6)=P(1,3)
        P(2,6)=P(3,3)
        P(3,6)=P(2,3)  
        
      ELSE IF (I.EQ.J) THEN       
        
!     DO NOT ALLOW (1,2) OPERATIONS
!     ALLOW ONLY SINGLE ACTIONS     
      
        ID=3     
        
!     ACTION OF THE INDENTITY OPERATOR

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3) 

!     ACTION OF THE (1,3) OPERATOR

        P(1,2)=INTER(3)
        P(2,2)=INTER(2)
        P(3,2)=INTER(1) 

!     ACTION OF THE (2,3) OPERATOR

        P(1,3)=INTER(1)
        P(2,3)=INTER(3)
        P(3,3)=INTER(2)
        
      ELSE IF (J.EQ.K) THEN       
        
!     DO NOT ALLOW (2,3) OPERATIONS
!     ALLOW ONLY SINGLE ACTIONS       
      
        ID=3     
        
!     ACTION OF THE INDENTITY OPERATOR

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)      

!     ACTION OF THE (1,2) OPERATOR

        P(1,2)=INTER(2)
        P(2,2)=INTER(1)
        P(3,2)=INTER(3)  

!     ACTION OF THE (1,3) OPERATOR 

        P(1,3)=INTER(3)
        P(2,3)=INTER(2)
        P(3,3)=INTER(1) 
        
      ELSE IF (I.EQ.K) THEN    
      
!     DO NOT ALLOW (1,3) OPERATIONS 
!     ALLOW ONLY SINGLE ACTIONS       
      
        ID=3  
        
!     ACTION OF THE INDENTITY OPERATOR

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)      

!     ACTION OF THE (1,2) OPERATOR

        P(1,2)=INTER(2)
        P(2,2)=INTER(1)
        P(3,2)=INTER(3)  

!     ACTION OF THE (2,3) OPERATOR

        P(1,3)=INTER(1)
        P(2,3)=INTER(3)
        P(3,3)=INTER(2)        
    
      END IF

      RETURN 
      END 

!################################################################################
! SUBROUTINE TO CALCULATE THE NUMBER OF COEFFICIENTS NC FOR A GIVEN POLYNOMIAL
! OF ORDER "ORDER" AND MOLECULE TYPE "MOLTYP"
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
!
!     FOR THE CASE OF AN ABC MOLECULE 
! 
              IF (MOLTYP.EQ."ABC") THEN  

                NC=NC+1
                MNUM=MNUM+1
!
!     FOR THE CASE OF AN AB2 MOLECULE 
!
              ELSE IF (MOLTYP.EQ."AB2") THEN 

                IF (J.LE.K) THEN 

                  NC=NC+1
                  MNUM=MNUM+1

                END IF    
!
!     FOR THE CASE OF AN A3 MOLECULE 
!         
              ELSE IF (MOLTYP.EQ."A3") THEN

                IF (I.LE.J .AND. J.LE.K) THEN

                  NC=NC+1
                  MNUM=MNUM+1

                END IF
!
!     END OF CASES
!
              END IF

            END IF
          END DO
        END DO
      END DO

      RETURN
      END
