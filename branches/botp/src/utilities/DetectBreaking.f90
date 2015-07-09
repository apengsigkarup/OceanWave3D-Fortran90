SUBROUTINE detect_breaking(fop,nfx,Wave,istep)

  ! A subroutine to identify and update surface rollers in 1 horizontal 
  ! dimension.  This is based on the Zelt model but using roller identification 
  ! based on the surface slope.  The geometry of the rollers found are output 
  ! to unit 'fop' every full time-step which is indicated by 'istep=1'. 

  USE Precision
  USE Constants
  USE DataTypes
  USE GlobalVariables, ONLY: FineGrid, BreakMod, g, rho, dt, dx, tstep, nsteps
  IMPLICIT NONE

  INTEGER fop, istep, nfx
  Type (Wavefield_FS) Wave

  ! Local workspace
  CHARACTER*3 CHOUT3
  logical :: first_call=.true.
  save    first_call
  INTEGER i, j, jr, jm, nrmax, ir, ir0, nx_ramp, nr, i_model, io, i_Us_max
  REAL(kind=long):: tphib, v, c1, c2, sxm, a, b, c, d, sr, dxr, dxm,   &
       x, xm, xr, dx_ramp, ramp, Usx, r1, x0, x1, UxStar, h_roller, bgamma, &
       d_max, quart=one/four, dx_inv

  ! Pointers to make the code more readable.  
  REAL(kind=long), POINTER :: et(:,:), Vs(:,:), Ws(:,:), Qr_x(:,:), Mr_t(:,:)
  ! local workspace
  REAL(kind=long) :: Qr(nfx), Mr(nfx), eta_x(nfx), tan_phi_old(nfx)
  !
  ! Local pointers to simplify the coding.
  !
  et => Wave%E; Vs => Wave%Px; Ws => Wave%W; 
  Mr_t=>Wave%Mr_t; Qr_x=>Wave%Qr_x;

  nrmax=nfx/10
  !
  IF(first_call)THEN
     ! Set up on the first call
     first_call=.false.
     ALLOCATE ( BreakMod%roller_thickness(nfx), BreakMod%tan_phi(nfx),         &
          BreakMod%i_roller(2,nrmax), BreakMod%eta_roller(nrmax),       &
          BreakMod%x_roller(nrmax) )
     OPEN(unit=fop,file='roller.dat',status='unknown')
     WRITE(fop,12)
12   FORMAT('# Breaking model data:  itime,n_rollers,(i_roller(1,i),i_roller(2,i),i=1,n_rollers),',/, &
          '# (x_roller(i),eta_roller(i),roller_thickness(jr-1),tan_phi(i_roller(2,i)),  i=1,n_rollers)',/, &
          '# i_roller(1,i)=number of grid points in the roller, i_roller(2,i)=jr, roller toe grid point')
     PRINT *, ' '
     PRINT *, '  The Zelt-type breaking model is being applied.  T_1/2=',BreakMod%T_half
     PRINT *, '  phi_b=',BreakMod%tan_phi_b,' phi_0=',BreakMod%tan_phi_0 
     PRINT *, '  f_delta=',BreakMod%del_fac,' f_c=',BreakMod%cel_fac
     PRINT *, '  gamma=',BreakMod%gamma_break
     PRINT *,' '
     IF(BreakMod%i_breaking==2)THEN
        PRINT *, '  Roller geometry will be computed for the whole run.'
        IF(BreakMod%i_break_time<nsteps)THEN
           PRINT *, '  The Zelt breaker model will be turned on at t=',&
                (BreakMod%i_break_time-1)*dt
           print *, ' '
        ELSE
           print*, ' The roller terms will not be included in the dynamics.'
           print *, ' '
        END IF
     END IF
     BreakMod%tan_phi_b=Tan(deg2rad*BreakMod%tan_phi_b)
     BreakMod%tan_phi_0=Tan(deg2rad*BreakMod%tan_phi_0)
     BreakMod%tan_phi=BreakMod%tan_phi_b
     BreakMod%eta_roller=zero
  END IF

  ! Initialization
  c1=one/six
  c2=two/three
  do j=1,nfx
     BreakMod%roller_thickness(j)=zero
     tan_phi_old(j)=BreakMod%tan_phi(j)
     Qr(j)=zero
     Mr(j)=zero
     Qr_x(j,1)=zero
     Mr_t(j,1)=zero
  end do
  do j=1,nrmax
     BreakMod%x_roller(j)=zero
  end do
  !
  dx_inv=one/dx;
  ! eta_x at j+1/2.  This is just a second-order differencing, and really should be 
  ! made consistent with the rest of the differencing -hbb
  do j=2,nfx-2
     eta_x(j)=dx_inv*(et(j+1,1)-et(j,1))
  end do

  ! We work back from the right end of the domain looking for roller toes.  
  ! When a roller toe is found, all the roller points for that roller are 
  ! adjusted, and the loop counter incremented to the end of the roller+1.  

  BreakMod%n_rollers=0
  DO  j = nfx-2,3,-1

     IF(j<2)EXIT

     tphib = BreakMod%tan_phi(j)

     IF (eta_x(j-1) .lt. -tphib .and. eta_x(j) .ge. -tphib             &
          .and. BreakMod%roller_thickness(j) <= 0.0 ) then
        ! If this grid point and the one behind it have slopes which straddle 
        ! -tan(phi_b), then the roller toe should be between j and j-1 (unless 
        ! a roller thickness has already been assigned to this point.)  

        BreakMod%n_rollers=BreakMod%n_rollers+1  ! a new roller

        ! Find the sub-grid position of the roller toe.  
        a   = half*dx_inv*(eta_x(j)-eta_x(j-1)) ! eta_xx/2 at j
        b   = half*(eta_x(j)+eta_x(j-1))    ! eta_x    at j
        c   = et(j,1)                         ! eta      at j
        ! The toe point is where eta_x + eta_xx * dx = -tan(phi_b)
        dxr = -(tphib+b)/(two*a) 
        xr  = real(j-1,long)*dx + dxr ! The x-position of the roller toe.
        jr  = int(xr*dx_inv) + 2 ! The grid point just in front of the roller toe
        ! The elevation at the toe is eta = eta(j) + eta_x(j)*dx + eta_xx(j)*dx^2/2.  
        ! sr  = a*dxr*dxr + b*dxr + c   
        sr = a*dxr*dxr + b*dxr +c
        BreakMod%x_roller(BreakMod%n_rollers)=xr
        BreakMod%eta_roller(BreakMod%n_rollers)=sr

        ! The simplest form of the wave celerity.
        v = BreakMod%cel_fac*sqrt(g*FineGrid%h(j,1))

        ! Work back along the roller to find:  the number of grid points in this 
        ! roller, and for each point:  the roller thickness, and the rotational 
        ! depth integrated flux and momentum Qr_x and Mr_t.  These are defined 
        ! by Q=(c-Us)del; M=(c^2-Us^2)del; Qr_x=d/dx Q; 
        ! Mr_t=1/(h+eta) d/dx(M - c Q).  

        sxm = eta_x(jr-1)
        jm  = jr - 1
        ! The roller pointer tells us the number of grid points in the roller and 
        ! the grid point just in front of the toe.  
        BreakMod%i_roller(1,BreakMod%n_rollers)=0  
        BreakMod%i_roller(2,BreakMod%n_rollers)=jr 
        h_roller=FineGrid%h(jr,1)
        UxStar=0.3_long*sqrt(g/h_roller)
        d_max=zero
        DO  i = (jr-1),3,-1
           d = BreakMod%del_fac*(et(i,1) - (sr+tphib*(xr-real(i-1,long)*dx)))
           ! The back of the roller is the point where the tangent line crosses eta.  
           IF (d < zero) EXIT 
           IF (d > (et(i,1)-sr)) d = et(i,1) - sr
           IF(d>d_max)d_max=d
           BreakMod%roller_thickness(i) = d
           BreakMod%i_roller(1,BreakMod%n_rollers)=BreakMod%i_roller(1,BreakMod%n_rollers)+1
           Usx = half*dx_inv*((Vs(i+1,1) - Wave%Ex(i+1,1) * Ws(i+1,1))- &
                (Vs(i-1,1) - Wave%Ex(i-1,1) * Ws(i-1,1)))
           Qr(i) = g*BreakMod%roller_thickness(i)
           IF (eta_x(i) .lt. sxm) THEN
              sxm = eta_x(i)
              jm  = i
           ENDIF
        END DO
        if(d_max==0)d_max=1;
        Mr(jr-1:i)= Mr(jr-1:i)/d_max
        ! The maximum negative slope point.
        ! dx^2/2 * eta_xxx at j+1/2
        a   = half*(eta_x(jm+1)-2.0e0*eta_x(jm)+eta_x(jm-1))
        b   = half*(eta_x(jm+1)-eta_x(jm-1))
        ! dx * eta_xx at j+1/2
        dxm = -b/(two*a)
        xm  = (jm+half+dxm)*dx

     ENDIF
  END do
  ! Updata tanphi for the next time step.  All values of tan(phi) from 3 grid points 
  ! ahead of the toe to the back of the roller get the new value.  
  IF(istep==1) THEN
     BreakMod%tan_phi=BreakMod%tan_phi_b
     DO j=1,BreakMod%n_rollers
        tphib=tan_phi_old(BreakMod%i_roller(2,j))
        DO i=1,BreakMod%i_roller(1,j)+3
           ! Linear decay of the roller angle
           BreakMod%tan_phi(BreakMod%i_roller(2,j)+4-i) = tphib + dt/BreakMod%T_half*(BreakMod%tan_phi_0-tphib)
           ! Exponential decay of the roller angle
           !            BreakMod%tan_phi(BreakMod%i_roller(2,j)+4-i) = tphib + LOG(two)*dt/BreakMod%T_half*(BreakMod%tan_phi_0-tphib)
        END DO
     END DO
  END IF

  IF(BreakMod%n_rollers>0)THEN
     IF(BreakMod%i_breaking==1)THEN
        ! Smooth the pressure.  
        DO j=3,nfx-2
           !            Mr_t(j,1)=half*dx_inv*(Mr(j+1)-Mr(j-1)) ! Second order deriv.
           ! Second order deriv. plus 3 point smoothing, alpha=1/4.
!           Mr_t(j,1)=eighth*dx_inv*(Mr(j+2)+two*Mr(j+1)-two*Mr(j-1)-Mr(j-2)) 
           Qr_x(j,1)=(quart*Qr(j+1)+half*Qr(j)+quart*Qr(j-1)) 
        END DO
     ELSE
        ! For BreakMod%i_breaking==2, we compute the breaker geometry, but don't include 
        ! the dynamics.  
     END IF
  END IF
  !
  ! Write out the roller geometry data and update the roller angle just once per 
  ! time step.
  IF(istep>0 .and. BreakMod%n_rollers>0)THEN

     WRITE(fop,*)tstep,BreakMod%n_rollers,                                                    &
          (BreakMod%i_roller(1,i),BreakMod%i_roller(2,i),i=1,BreakMod%n_rollers),             &
          (BreakMod%x_roller(i),BreakMod%eta_roller(i), BreakMod%roller_thickness(BreakMod%i_roller(2,i)-1),         &
          BreakMod%tan_phi(BreakMod%i_roller(2,i)),  i=1,BreakMod%n_rollers)
  END IF

  !print*, 'max. Qr_x, Mr_t=',maxval(abs(Qr_x)),maxval(abs(Mr_t))

  RETURN
END SUBROUTINE detect_breaking
