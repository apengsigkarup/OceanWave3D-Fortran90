 FUNCTION GASDEV(IDUM)
   implicit none
   INTEGER IDUM
   REAL gasdev
! Uses ran1, returns a normally distributed deviate with zero mean and
! unit variance.
   Integer iset
   REAL :: fac,gset,rsq,v1,v2,ran1,zero=0.,one=1.,two=2.
   SAVE iset, gset
   DATA iset/0/
   IF (ISET.EQ.0) THEN
1     V1=two*RAN1(IDUM)-one
      V2=two*RAN1(IDUM)-one
      rsq=V1**2+V2**2
      IF(rsq >= one .or. rsq == zero)GO TO 1
      FAC=SQRT(-two*LOG(Rsq)/Rsq)
      GSET=V1*FAC
      GASDEV=V2*FAC
      ISET=1
   ELSE
      GASDEV=GSET
      ISET=0
   ENDIF
   RETURN
 END FUNCTION GASDEV
