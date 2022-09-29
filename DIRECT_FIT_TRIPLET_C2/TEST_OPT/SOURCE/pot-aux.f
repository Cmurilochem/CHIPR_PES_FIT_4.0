C234567
      function POTEN(R,NRD)
      implicit real*8(A-H,O-Z)
      CALL CHIPR_DIAT(R,POT,.FALSE.,DVDR)
      POTEN=POT
      RETURN
      END

c  !!!****************************************************************!!!
