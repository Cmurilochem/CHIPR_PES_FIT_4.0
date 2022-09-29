      PROGRAM TEST_PROGRAM
      IMPLICIT NONE
      INTEGER :: I,J,K,L,M,N,O,Q,R,S,LL,CHECK
      INTEGER :: TOTNUM,MNUM,ORDER
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: C,POL
      DOUBLE PRECISION :: V4
      INTEGER, DIMENSION(6,24) :: P
      INTEGER :: ID,NC
      CHARACTER(LEN=4) :: MOLTYP
      DOUBLE PRECISION :: Y1,Y2,Y3,Y4,Y5,Y6
      REAL :: START,FINAL       
!
!     INITIAL USER DEFINED OPTIONS
!
      ORDER=4
      MOLTYP="A4"
!
!     CALCULATE THE NUMBER OF COEFFICIENTS FOR EACH MOLECULE TYPE 
!     AND POLYNOMIALS ORDER
!
      CALL POLNCTETRA(ORDER,MOLTYP,NC)
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
      Y4=1.0D+00
      Y5=1.0D+00
      Y6=1.0D+00      
!
!     
! 
      CALL CPU_TIME(START)

      TOTNUM=0
      S=0
      CHECK=0
!$OMP DO
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

                     TOTNUM=TOTNUM+1
                     MNUM=MNUM+1  
                     
                     CALL PERMUTABCD(I,J,K,L,N,O,P,ID)

                     DO LL=1,ID
                     
                       CHECK=CHECK+1
                       
                       POL(TOTNUM)=POL(TOTNUM)+&
     &                    (Y1**(P(1,LL))*Y2**(P(2,LL))*Y3**(P(3,LL))*&
     &                     Y4**(P(4,LL))*Y5**(P(5,LL))*Y6**(P(6,LL)))                       
 
                       PRINT*, P(1,LL),P(2,LL),P(3,LL),P(4,LL),P(5,LL),P(6,LL),TOTNUM,MNUM,CHECK,POL(TOTNUM)

                     END DO
!
!     NORMALIZING COEFFICIENT
!
                     POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)                     
!
!     FOR THE CASE OF AN AB3 MOLECULE 
!
                   ELSE IF (MOLTYP.EQ."AB3") THEN

                     IF (I.LE.J .AND. L.LE.N .AND. N.LE.O .AND. J.LE.K) THEN

                       TOTNUM=TOTNUM+1
                       MNUM=MNUM+1

                       CALL PERMUTAB3(I,J,K,L,N,O,P,ID)
                       
                       DO LL=1,ID
                       
                         CHECK=CHECK+1
                         
                         POL(TOTNUM)=POL(TOTNUM)+&
     &                      (Y1**(P(1,LL))*Y2**(P(2,LL))*Y3**(P(3,LL))*&
     &                       Y4**(P(4,LL))*Y5**(P(5,LL))*Y6**(P(6,LL)))  
 
                         PRINT*, P(1,LL),P(2,LL),P(3,LL),P(4,LL),P(5,LL),P(6,LL),TOTNUM,MNUM,CHECK,POL(TOTNUM)
                         
                       END DO
                       
!
!     NORMALIZING COEFFICIENT
!
                       POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID) 
                       
                       WRITE(*,*) POL(TOTNUM) 

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
                       
                         CHECK=CHECK+1
                         
                         POL(TOTNUM)=POL(TOTNUM)+&
     &                      (Y1**(P(1,LL))*Y2**(P(2,LL))*Y3**(P(3,LL))*&
     &                       Y4**(P(4,LL))*Y5**(P(5,LL))*Y6**(P(6,LL)))  
 
                         PRINT*, P(1,LL),P(2,LL),P(3,LL),P(4,LL),P(5,LL),P(6,LL),TOTNUM,MNUM,CHECK,POL(TOTNUM)
                         
                       END DO                       
!
!     NORMALIZING COEFFICIENT
!
                       POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID) 
                       
                       WRITE(*,*) POL(TOTNUM) 
                     
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
                       
                         CHECK=CHECK+1
                         
                         POL(TOTNUM)=POL(TOTNUM)+&
     &                      (Y1**(P(1,LL))*Y2**(P(2,LL))*Y3**(P(3,LL))*&
     &                       Y4**(P(4,LL))*Y5**(P(5,LL))*Y6**(P(6,LL)))  
 
                         PRINT*, P(1,LL),P(2,LL),P(3,LL),P(4,LL),P(5,LL),P(6,LL),TOTNUM,MNUM,CHECK,POL(TOTNUM)
                         
                       END DO 
!
!     NORMALIZING COEFFICIENT
!
                       POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID) 
                       
                       WRITE(*,*) POL(TOTNUM)                    
                     
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
                       
                         CHECK=CHECK+1
                         
                         POL(TOTNUM)=POL(TOTNUM)+&
     &                      (Y1**(P(1,LL))*Y2**(P(2,LL))*Y3**(P(3,LL))*&
     &                       Y4**(P(4,LL))*Y5**(P(5,LL))*Y6**(P(6,LL)))  
 
                         PRINT*, P(1,LL),P(2,LL),P(3,LL),P(4,LL),P(5,LL),P(6,LL),TOTNUM,MNUM,CHECK,POL(TOTNUM)
                         
                       END DO 
!
!     NORMALIZING COEFFICIENT
!
                       POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID) 
                       
                       WRITE(*,*) POL(TOTNUM)                      
                     
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

        WRITE(*,*)
        WRITE(*,*)
 
      END DO
!$OMP END DO

      CALL CPU_TIME(FINAL)

      V4=SUM(POL)

      WRITE(*,*) V4,(FINAL-START)

      DEALLOCATE(C,POL)

      CONTAINS
      
!################################################################################
! FINDING THE POSSIBLE PERMUTATIONS OF I J K R S T FOR AN ABCD-TYPE MOLECULE
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
! FINDING THE POSSIBLE PERMUTATIONS OF I J K R S T FOR AN A4-TYPE MOLECULE
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
! FINDING THE POSSIBLE PERMUTATIONS OF I J K R S T FOR AN ABC2-TYPE MOLECULE
! Y1 IS ALWAYS C-C DISTANCES
! Y2 AND Y4 ARE A-C DISTANCES  
! Y3 AND Y5 ARE B-C DISTANCES
! Y6 CORRESPONDS TO THE A-B BONDS
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
! FINDING THE POSSIBLE PERMUTATIONS OF I J K R S T FOR AN A2B2-TYPE MOLECULE
! Y1 IS ALWAYS A-A BONDS
! Y2-Y5 REPRESENT A-B BONDS
! Y6 CORRESPONDS TO THE B-B BONDS
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
! FINDING THE POSSIBLE PERMUTATIONS OF I J K R S T FOR AN AB3-TYPE MOLECULE
! Y1, Y2, AND Y3 IS ALWAYS THE A-B BONDS
! Y4, Y5, AND Y6 CORRESPONDS TO THE B-B BONDS
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
! SUBROUTINE TO CALCULATE THE NUMBER OF COEFFICIENTS NC FOR A GIVEN POLYNOMIAL
! OF ORDER "ORDER" AND MOLECULE TYPE "MOLTYP"
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
      
      END PROGRAM
