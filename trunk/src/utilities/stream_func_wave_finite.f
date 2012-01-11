
      
      subroutine stream_func_wave_finite(nfx,xp,tp,n_four_modes,
     $     zz,yy,wavenum,grav,etai,phi_s)

! A subroutine to compute the dimensional elevation and surface 
! potential  at a set of points (xp) on the free-surface at time tp, 
! based on the input stream function theory solution which has been 
! computed by the subroutine stream_func_coeffs.  The solution is 
! uniform in y.  

      IMPLICIT NONE
      INTEGER nfx, n_four_modes
C      REAL*8 xp(*), tp, zz(*), yy(*), etai(*), phi_s(*), wavenum, 
C     $     grav
C Test GD
      REAL*8 xp(nfx), tp, zz(2*n_four_modes+10), yy(2*n_four_modes+10), 
     $     etai(nfx), phi_s(nfx), wavenum, grav
C Harry's definition
C      REAL*8 xp(2,*), tp, zz(*), yy(*), etai(*), phi_s(*), wavenum, 
C     $     grav
C Local variables
      INTEGER j, k, m
      real*8 SLOPE,U,W,E,S,B,C,ufactor,t,elx,uxx,d, omega_t, kinv, kx,
     $     ph, pfactor, km, sinhzh
      kinv=1.D00/wavenum
      ufactor=sqrt(grav*kinv)
      pfactor=ufactor*kinv
      omega_t=wavenum*ufactor*zz(4)*tp
      DO j=1,nfx
         kx=wavenum*xp(j)-omega_t
C         kx=wavenum*xp(1,j)-omega_t
         km=DBLE(n_four_modes)
         etai(j)=0.5d0*(yy(n_four_modes))*cos(km*kx)
         slope=-0.5d0*km*yy(n_four_modes)*sin(km*kx)
         DO  m=1,n_four_modes-1
            km=DBLE(m)
! Non-dimensional eta
            etai(j)=etai(j)+yy(m)*cos(km*kx)
            slope=slope-km*yy(m)*sin(km*kx)
         ENDDO
         ph=0.D00
         w=0.d00
!         ph=(zz(4)-zz(7))*xp(1,j)
C Here phi, is computed on the free-surface etai.  To get values 
C at another level of z, replace etai with z in the following two loops.  
         DO m=1,n_four_modes
            km=DBLE(m)
            e=EXP(km*(zz(1)+etai(j)))
            b=zz(n_four_modes+m+10)*0.5d0*(e+1.d0/e)/COSH(km*zz(1))
            sinhzh=km*zz(n_four_modes+m+10)*0.5d0
     &                              *(e-1.d0/e)/COSH(km*zz(1))
            ph=ph+b*SIN(km*kx)
            w=w+sinhzh*SIN(km*kx)
         ENDDO
! Dimensionalize
         etai(j)=kinv*etai(j)
         phi_s(j)=ph*pfactor
!	Next line ADDED: APEK 01-07-2008, ! phi_s=uE*x+ph, uE=c-ubar, 
         phi_s(j)=phi_s(j) + xp(j)*(zz(5)*sqrt(grav/wavenum))   
         w=w*ufactor
!hbb temporary output for checking the GOperator
!         write(26,*)phi_s(j),w
      END DO

      return
      end subroutine stream_func_wave_finite

      subroutine stream_func_wave_deep(nfx,xp,tp,n_four_modes,
     $     zz,yy,wavenum,grav,etai,phi_s)

C A subroutine to compute the dimensional elevation and velocity at 
C a set of points (xp) on the free-surface at time tp, based on the 
C input stream function theory solution which has been computed by 
C the subroutine stream_func_coeffs.  The solution is uniform in 
C y.  

      IMPLICIT NONE
      INTEGER nfx, n_four_modes
      REAL*8 xp(*), tp, zz(*), yy(*), etai(*), phi_s(*), wavenum, 
     $     grav
C Local variables
      INTEGER j, k, m
      real*8 SLOPE,U,W,E,S,B,C,ufactor,t,elx,uxx,d, omega_t, kinv, kx,
     $     ph, pfactor, km

      kinv=1.D00/wavenum
      ufactor=sqrt(grav*kinv)
      pfactor=ufactor*kinv
      omega_t=wavenum*ufactor*zz(4)*tp
      DO j=1,nfx
         kx=wavenum*xp(j)-omega_t
         km=DBLE(n_four_modes)
         etai(j)=0.5d0*yy(n_four_modes)*cos(km*kx)
         slope=-0.5d0*km*yy(n_four_modes)*sin(km*kx)
         DO  m=1,n_four_modes-1
            km=DBLE(m)
            etai(j)=etai(j)+yy(m)*cos(km*kx)
            slope=slope-km*yy(m)*sin(km*kx)
         ENDDO
         ph=0.D00
!         ph=(zz(4)-zz(7))*xp(1,j)
C Here phi is computed on the free-surface etai.  To get values 
C at another level of z, replace etai with z in the following two loops.  
         DO m=1,n_four_modes
            km=DBLE(m)
            ph=ph+zz(n_four_modes+m+10)*EXP(km*etai(j))*SIN(km*kx)
         ENDDO
         etai(j)=kinv*etai(j)
         phi_s(j)=ph*pfactor
!	Next line ADDED: APEK 01-07-2008, ! phi_s=uE*x+ph, uE=c-ubar, 
         phi_s(j)=phi_s(j) + xp(j)*(zz(5)*sqrt(grav/wavenum))   

      END DO

      return
      end subroutine stream_func_wave_deep

	  subroutine weights (z,x,n,nd,m,c)
C---
C--- Input Parameters
C--- z location where approximations are to be accurate,
C--- x(0:nd) grid point locations, found in x(0:n)
C--- n one less than total number of grid points; n must
C--- not exceed the parameter nd below,
C--- nd dimension of x- and c-arrays in calling program
C--- x(0:nd) and c(0:nd,0:m), respectively,
C--- m highest derivative for which weights are sought,
C--- Output Parameter
C--- c(0:nd,0:m) weights at grid locations x(0:n) for derivatives
C--- of order 0:m, found in c(0:n,0:m)
C---
	  implicit real*8 (a-h,o-z)
	  dimension x(0:nd),c(0:nd,0:m)
	  c1 = 1.0d0
	  c4 = x(0)-z
	  do 10 k=0,m
		do 10 j=0,n
   10			c(j,k) = 0.0d0
		c(0,0) = 1.0d0
		do 50 i=1,n
			mn = min(i,m)
			c2 = 1.0d0
			c5 = c4
			c4 = x(i)-z
			do 40 j=0,i-1
				c3 = x(i)-x(j)
				c2 = c2*c3
				if (j.eq.i-1) then
					do 20 k=mn,1,-1
   20						c(i,k) = c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
					c(i,0) = -c1*c5*c(i-1,0)/c2
				endif
				do 30 k=mn,1,-1
   30					c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
   40			c(j,0) = c4*c(j,0)/c3
   50		c1 = c2
	  return
	  end subroutine weights

