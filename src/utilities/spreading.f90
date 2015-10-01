SUBROUTINE  SPREADING(BETA1,PI,S0,CTA0,SEED2,NCTA)

INTEGER S,I,J,NT,NCTA,NCTA1,SEED2
REAL*8 CTAMAX,CTAMIN,DCTA,CTA0,IPOW,DS,SUM1,AAA,DCTA1,S0,PI,SPM
REAL*8 CTA(16384),COSINT(16384),CC1(16384),BETA1(NCTA),PS(16384)
 


  IPOW=2*S0
 
  NCTA1=16384
  
!!!!!  out put
!      OPEN(1110,FILE='spreading.dat')
 
!   write(1110,*) s0,cta0,ncta
 
!  close(1110)

  CTAMAX=CTA0+0.5*PI
  CTAMIN=CTA0-0.5*PI
 
!FOR THE MULTIDIRECTIONAL SPREADING  'G0'

 
       DCTA1=(CTAMAX-CTAMIN)/FLOAT(NCTA1-1)
 
	DO I=1,NCTA1
        CC1(I)=CTAMIN+(I-1)*DCTA1
        END DO  
 
        SUM1=0.0
        DO I=1,NCTA1                             
        AAA=CC1(I)

        SUM1=SUM1+((COS((AAA-CTA0)/2.))**IPOW)*DCTA1
        END DO                                     
        SPM=1./SUM1     

 


  CTA(1)=CTAMIN
  CTA(NCTA)=CTAMAX

  DS=0.0

  DCTA=PI/REAL(NCTA)
!  OPEN(1,FILE='1 ',status='unknown')
!  OPEN(2,FILE='2.TXT')
  DO I=1,NCTA

  CTA(I)=CTAMIN+(I-1)*DCTA

  DS=DS+((COS((CTA(I)-CTA0)/2.0))**IPOW)*DCTA

  COSINT(I)=DS*SPM
!   WRITE(*,*) i
!   WRITE(1,*)   CTA(I) ,i 
  END DO

!  OPEN(4,FILE='4.TXT')
!  OPEN(3,FILE='3.TXT')
  DO J=1,NCTA 

 
   PS(J)=RAN1_b(SEED2)
!   WRITE(2,*)   PS(J),SEED2,J,spm
  DO I=1,NCTA

  IF (COSINT(I).GE.PS(J)) then
!   WRITE(4,*)   CTA(I),COSINT(I),PS(J),I
  GO TO 100
    
  end if


  
  END DO
    
 100  BETA1(J)=CTA(I)
 

 !  WRITE(3,*)   CTA(I),COSINT(I),PS(J),I
 
  END DO

!!!!  output
!   OPEN(1110,FILE='anglespreading')
!  do i=1,ncta
!   write(1110,*) beta1(i)
!   end do
!  close(1110)
 
  return
  END SUBROUTINE SPREADING


