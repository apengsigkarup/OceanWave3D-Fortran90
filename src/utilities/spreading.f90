SUBROUTINE  SPREADING(BETA,PI,S0,BETA0,SEED2,NCTA)
  !
  ! This is a very simple implementation of assigning a random heading angle
  ! to each wave component based on the cosine spreading function.
  !
   IMPLICIT none
   integer, parameter :: long=selected_real_kind(12,99)

   INTEGER S,I,J,NT,NCTA,SEED2

   REAL  ran1_b
   REAL(kind=long) :: CTAMAX,CTAMIN,DCTA,BETA0,IPOW,DS,SUM1,AAA,DCTA1,S0,PI,SPM
   REAL(kind=long) :: CTA(ncta),COSINT(ncta),CC1(ncta),BETA(NCTA),PS(ncta)

   EXTERNAL ran1_b


   
   IPOW=2*S0

   CTAMAX=BETA0+0.5*PI
   CTAMIN=BETA0-0.5*PI

   !FOR THE MULTIDIRECTIONAL SPREADING  'G0'

   DCTA1=(CTAMAX-CTAMIN)/FLOAT(NCTA-1)

   DO I=1,NCTA
      CC1(I)=CTAMIN+(I-1)*DCTA1
   END DO

   SUM1=0.0
   DO I=1,NCTA                            
      AAA=CC1(I)

      SUM1=SUM1+((COS((AAA-BETA0)/2.))**IPOW)*DCTA1
   END DO
   SPM=1./SUM1     




   CTA(1)=CTAMIN
   CTA(NCTA)=CTAMAX

   DS=0.0

   DCTA=PI/REAL(NCTA)
   !
   ! Build the discrete probability function as a function of heading angle.
   !
   DO I=1,NCTA

      CTA(I)=CTAMIN+(I-1)*DCTA

      DS=DS+((COS((CTA(I)-BETA0)/2.0))**IPOW)*DCTA

      COSINT(I)=DS*SPM
   END DO
   
   !
   ! Loop over each component, give it a random probability between 0 and 1, and find
   ! the associated heading angle.
   !
   DO J=1,NCTA 

      PS(J)=RAN1_b(SEED2)

      !
      ! This should be replaced by a more efficient root-finding routine.
      !
      DO I=1,NCTA

         IF (COSINT(I).GE.PS(J)) then
            GO TO 100

         end if

      END DO

100   BETA(J)=CTA(I)

   END DO


   return
 END SUBROUTINE SPREADING


