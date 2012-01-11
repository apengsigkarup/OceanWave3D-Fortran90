SUBROUTINE LocalSmoothing2D(nx,ny,itime,fop)
!
! This subroutine attempts to apply automatic localized smoothing to avoid any 
! breakdowns due to discontinuities in the free-surface.  We check the 
! surface acceleration everywhere and when it exceeds a given factor times g 
! we smooth a 10-point region centered at the point with the classical 3-point 
! filter [.25 .5 .25].  The feature is turned off by switching the threshold 
! percentage to something very large.  
!
! **  This still needs to be extended to 3D. **  -hbb
!
  USE GlobalVariables
  IMPLICIT NONE
  INTEGER :: nx, ny, itime, fop

! Local variables and automatic workspace
  INTEGER, SAVE :: i_first=0
  INTEGER ::  i, j, i_flag, N_smoothed, i_smoothed(nx)
  REAL(KIND=long)  accel_tol, dt_inv, Max_accel
  REAL(KIND=long), POINTER :: et(:,:), ph(:,:), x(:,:), EtatHist(:,:,:)
  REAL(KIND=long) :: etatem(nx,ny), Phitem(nx,ny), etatt(nx,ny), alpha_s(nx,ny)
!
! Initial set up of the output file
!
     IF(i_first==0)THEN
        i_first=1
        OPEN(fop,file='local.smoothing',status='unknown')
        WRITE(fop,10)accel_tol_fact
10      FORMAT('# Points where eta_tt exceeded',e10.2, &
             ' times g and local smoothing was applied:  t, x, eta_tt')
     END IF
!
! Do nothing for a factor larger than 100
!
  IF(accel_tol_fact .le. 100. .and. itime .ge. 3)THEN
!
! Search for high accelerations and smooth an area around points where they occur.  
!
     !
     ! Initialize
     !
     x => FineGrid%x; et => WaveField%E; ph => WaveField%P; EtatHist => WaveField%EtatHist
     etatem=et
     Phitem=ph
     alpha_s=zero
     !
     ! Take a 3-point one-sided difference to get the acceleration at the current 
     ! time step.  
     !
     dt_inv=one/dt
     j=1;             ! 2D code so far.
     DO i=1,nx
        etatt(i,j)=half*dt_inv*(three*EtatHist(i,j,1)-four*EtatHist(i,j,2)+EtatHist(i,j,3))
     END DO

     accel_tol=accel_tol_fact*g
     !
     ! First build a coefficient vector which is .25 at those points which surround a point of high 
     ! accelerations.  
     !
     i_flag=0
     N_smoothed=0; Max_accel=zero;
     DO i=6,nx-5
        IF(abs(etatt(i,j)) >= accel_tol)THEN
           i_flag=1
           N_smoothed=N_smoothed+1
           i_smoothed(N_smoothed)=i
           alpha_s(i-5:i+5,j)=0.25_long
!        ELSEIF(abs(etatt(i,j)) >= Max_accel) THEN
!           Max_accel=abs(etatt(i,j))
        END IF
     END DO
     IF(i_flag==1)THEN
        !
        ! Smooth the elevation and surface velocity.  
        !
        DO i=2,nx-1
           et(i,j)=alpha_s(i,j)*etatem(i-1,j)+(one-two*alpha_s(i,j))*etatem(i,j) &
                +alpha_s(i,j)*etatem(i+1,j)
           ph(i,j)=alpha_s(i,j)*Phitem(i-1,j)+(one-two*alpha_s(i,j))*Phitem(i,j) &
                +alpha_s(i,j)*Phitem(i+1,j)
 !          etatem(i,j)=alpha_s(i,j)*et(i-1,j)+(one-two*alpha_s(i,j))*et(i,j)+alpha_s(i,j)*et(i+1,j)
 !          Phitem(i,j)=alpha_s(i,j)*ph(i-1,j)+(one-two*alpha_s(i,j))*ph(i,j)+alpha_s(i,j)*ph(i+1,j)
 !          et(i,j)=alpha_s(i,j)*etatem(i-1,j)+(one-two*alpha_s(i,j))*etatem(i,j) &
 !               +alpha_s(i,j)*etatem(i+1,j)
 !          ph(i,j)=alpha_s(i,j)*Phitem(i-1,j)+(one-two*alpha_s(i,j))*Phitem(i,j) &
 !               +alpha_s(i,j)*Phitem(i+1,j)
        END DO
        WRITE(fop,11)
        DO i=1,N_smoothed
           WRITE(fop,12)(itime-1)*dt,x(i_smoothed(i),1),etatt(i_smoothed(i),1)
        END DO
11      FORMAT()
12      FORMAT(3e16.4)
     END IF

  END IF
  RETURN
END SUBROUTINE LocalSmoothing2D
