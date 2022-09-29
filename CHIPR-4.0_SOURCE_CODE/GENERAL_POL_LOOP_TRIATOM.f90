      PROGRAM TEST_PROGRAM
      IMPLICIT NONE
      INTEGER :: I,J,K,L,M,S,CHECK
      INTEGER :: TOTNUM,MNUM,ORDER
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: C,POL
      DOUBLE PRECISION :: V3
      INTEGER, DIMENSION(3,6) :: P
      INTEGER :: ID,NC
      CHARACTER(LEN=3) :: MOLTYP
      DOUBLE PRECISION :: Y1,Y2,Y3
!
!     INITIAL USER DEFINED OPTIONS
!
      ORDER=2
      MOLTYP="A3"
!
!     CALCULATE THE NUMBER OF COEFFICIENTS FOR EACH MOLECULE TYPE 
!     AND POLYNOMIALS ORDER
!
      CALL POLNC(ORDER,MOLTYP,NC)
!
!     ALLOCATING THE REQUIRED VECTORS
!
      ALLOCATE(C(NC),POL(NC))
!
!     FOR TESTING
!
      DO I=1,NC 
        C(I)=1.0D+00
        POL(I)=0.00D+00
      END DO
!
!     HERE, CALCULATE THE VALUES OF Y'S IN THE BASIS SET
!
      Y1=1.0D+00
      Y2=1.0D+00
      Y3=1.0D+00
!
!     
! 
      TOTNUM=0
      S=0
      CHECK=0

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

                TOTNUM=TOTNUM+1
                MNUM=MNUM+1

                CALL PERMUTABC(I,J,K,P,ID)

                DO L=1,ID

                  CHECK=CHECK+1

                  POL(TOTNUM)=POL(TOTNUM)+&
     &               (Y1**(P(1,L))*Y2**(P(2,L))*Y3**(P(3,L)))

                  WRITE(*,*) P(1,L),P(2,L),P(3,L),TOTNUM,MNUM,CHECK,POL(TOTNUM)

                END DO
!
!     NORMALIZING COEFFICIENT
!
                POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

!                WRITE(*,*) POL(TOTNUM)
!
!     FOR THE CASE OF AN AB2 MOLECULE 
!
              ELSE IF (MOLTYP.EQ."AB2") THEN 

                IF (J.LE.K) THEN 

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTAB2(I,J,K,P,ID)
                   
                  DO L=1,ID

                    CHECK=CHECK+1

                    POL(TOTNUM)=POL(TOTNUM)+&
     &                 (Y1**(P(1,L))*Y2**(P(2,L))*Y3**(P(3,L)))

                    WRITE(*,*) P(1,L),P(2,L),P(3,L),TOTNUM,MNUM,CHECK,POL(TOTNUM)

                  END DO
!
!     NORMALIZING COEFFICIENT
!
                  POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

!                  WRITE(*,*) POL(TOTNUM) 

                END IF    
!
!     FOR THE CASE OF AN A3 MOLECULE 
!         
              ELSE IF (MOLTYP.EQ."A3") THEN

                IF (I.LE.J .AND. J.LE.K) THEN

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTA3(I,J,K,P,ID)

                  WRITE(*,*) I,J,K
                   
                  DO L=1,ID

                    CHECK=CHECK+1

                    POL(TOTNUM)=POL(TOTNUM)+&
     &                 (Y1**(P(1,L))*Y2**(P(2,L))*Y3**(P(3,L)))

                    WRITE(*,*) P(1,L),P(2,L),P(3,L),TOTNUM,MNUM,CHECK,POL(TOTNUM)

                  END DO
!
!     NORMALIZING COEFFICIENT
!
                  POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

!                  WRITE(*,*) POL(TOTNUM)

                END IF
!
!     END OF CASES
!
              END IF

            END IF
          END DO
        END DO

        WRITE(*,*)
 
      END DO

      V3=SUM(POL)

!      WRITE(*,*) V3

      DEALLOCATE(C,POL)

      CONTAINS

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

!     CALCULATE THE FULL S3 PERMUTATION GROUP 
      
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

!     ACTION OF THE (1,2)(2,3) OPERATOR   

        P(1,5)=P(1,2)
        P(2,5)=P(3,2)
        P(3,5)=P(2,2)  
      
!     ACTION OF THE (1,3)(3,2) OPERATOR 

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

      END PROGRAM
