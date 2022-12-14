!######################################################################################################

      SUBROUTINE CHIPR_SiC(R,POT,DER,DVDR)
      IMPLICIT NONE
      INTEGER, PARAMETER :: NC=24
      INTEGER :: I,J
      INTEGER :: BSORDER
      INTEGER :: POLORDER
      INTEGER :: NCBAS,NCPOL
      DOUBLE PRECISION, DIMENSION(2) :: Z
      DOUBLE PRECISION, DIMENSION(NC) :: C
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BS
      DOUBLE PRECISION :: Y   
      DOUBLE PRECISION :: R,POT
      LOGICAL :: DER
      DOUBLE PRECISION :: DYDR
      DOUBLE PRECISION :: DVDRPART1
      DOUBLE PRECISION :: DVDRPART2
      DOUBLE PRECISION :: DVDR

      Z(  1)=  0.60000000000000000000000000000000000D+01
      Z(  2)=  0.14000000000000000000000000000000000D+02

      C(  1)=  0.74748573762006467052776415016523970D-02
      C(  2)=  0.13896486594960813509835872991970973D+00
      C(  3)=  0.44031454894917558817724057007580996D+01
      C(  4)= -0.13888389959650467062601819634437561D+04
      C(  5)= -0.39173482714575526188127696514129639D+05
      C(  6)= -0.43432062720686977263540029525756836D+06
      C(  7)= -0.21963770319024808704853057861328125D+07
      C(  8)= -0.49699009718037340790033340454101562D+07
      C(  9)= -0.17303256821271706372499465942382812D+08
      C( 10)= -0.10317362706312467157840728759765625D+09
      C( 11)= -0.20444473089682644605636596679687500D+09
      C( 12)= -0.25330130046192488074302673339843750D+09
      C( 13)=  0.95823067790508639812469482421875000D+09
      C( 14)=  0.11930971123467397689819335937500000D+11
      C( 15)= -0.12401090338018053160773490617430070D+00
      C( 16)=  0.12915091364993452199838586125224538D-01
      C( 17)=  0.37487079955892822769047967312872061D-01
      C( 18)= -0.17897493475958460749097866937518120D+04
      C( 19)=  0.46295762003455565025689111280371435D+00
      C( 20)=  0.13113020004292406106571888813050464D+01
      C( 21)=  0.84564735156863768406054759907419793D+00
      C( 22)=  0.42514362518898690668223139255132992D+00
      C( 23)=  0.12998712639233773735725208098301664D+01
      C( 24)=  0.15795636862866988536779899732209742D+01

      BSORDER=4

      POLORDER=14

      NCBAS=2*BSORDER+2

      NCPOL=POLORDER

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

      SUBROUTINE CHIPR_C2(R,POT,DER,DVDR)
      IMPLICIT NONE
      INTEGER, PARAMETER :: NC=22
      INTEGER :: I,J
      INTEGER :: BSORDER
      INTEGER :: POLORDER
      INTEGER :: NCBAS,NCPOL
      DOUBLE PRECISION, DIMENSION(2) :: Z
      DOUBLE PRECISION, DIMENSION(NC) :: C
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BS
      DOUBLE PRECISION :: Y   
      DOUBLE PRECISION :: R,POT
      LOGICAL :: DER
      DOUBLE PRECISION :: DYDR
      DOUBLE PRECISION :: DVDRPART1
      DOUBLE PRECISION :: DVDRPART2
      DOUBLE PRECISION :: DVDR

      Z(  1)=  0.60000000000000000000000000000000000D+01
      Z(  2)=  0.60000000000000000000000000000000000D+01

      C(  1)=  0.25400653490294206743316252072872885D-01
      C(  2)=  0.29205851999229807958169402581916074D-01
      C(  3)=  0.28585799210269775827431004699974437D-01
      C(  4)=  0.67709595572945802804953885356553656D-02
      C(  5)= -0.30171664669936672925620868568330479D-01
      C(  6)= -0.31193025388905736700051463117233652D-01
      C(  7)=  0.70567454605977492088086933108570520D-01
      C(  8)=  0.19492105542621596114827298151794821D+00
      C(  9)=  0.17268454528472909625946840606047772D+00
      C( 10)=  0.46786282573429871511905275838216767D-01
      C( 11)= -0.15662574345293944766410731972428039D-01
      C( 12)= -0.82749902555507275775381614835168875D-02
      C( 13)= -0.16136303779771361474604240981989278D-02
      C( 14)= -0.43664393003147911054284691090288106D+00
      C( 15)= -0.32332760335351227176658994721947238D+01
      C( 16)=  0.61344284466357297787908464670181274D+05
      C( 17)=  0.11859219475394039422511127668258268D+00
      C( 18)=  0.34535683308796008006424926861654967D+01
      C( 19)=  0.10842372088218057424313656156300567D+01
      C( 20)=  0.57107351824085317293366870217141695D+00
      C( 21)=  0.11869360172429590516429698254796676D+01
      C( 22)=  0.19410370886066858897578413234441541D+01

      BSORDER=4

      POLORDER=12

      NCBAS=2*BSORDER+2

      NCPOL=POLORDER

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
