SUBROUTINE BeachGen2(h0,h1,l,x,h,hx,hxx,Nx)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: i, Nx
REAL(KIND=long) :: l, x(Nx), h0, h1, h(Nx), hx(Nx), hxx(Nx)
DO i = 1, Nx
    IF (x(i) > zero .AND. x(i) < l) THEN
        ! GD: Always use integer power when possible... (more precise and faster)
        h(i)   = h0 - (h0-h1)/two*(one+tanh(sin(pi*(x(i)-l/two)/l)/(one-(two*(x(i)-l/two)/l)**2)))
        hx(i)  = -(one/two*h0-one/two*h1)*(one-tanh(sin(pi*(x(i)-one/two*l)/l)/(one-(two*x(i)-l)**2/l**2))**2)*(cos(pi*(x(i)&
		-one/two*l)/l)*pi/l/(one-(two*x(i)-l)**2/l**2)+four*sin(pi*(x(i)-one/two*l)/l)/(one-(two*x(i)-l)**2/l**2)&
		**2*(two*x(i)-l)/l**2)
        hxx(i) = two*(half*h0-half*h1)*tanh(sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2))*(one-tanh(sin(pi*(x(i)-half*l)/l)&
		/(one-(two*x(i)-l)**2/l**2))**2)*(cos(pi*(x(i)-half*l)/l)*pi/l/(one-(two*x(i)-l)**2/l**2)+four*sin(pi*(x(i)-half*l)/l)&
		/(one-(two*x(i)-l)**2/l**2)**2*(two*x(i)-l)/l**2)**2-(half*h0-half*h1)*(one-tanh(sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)&
		-l)**2/l**2))**2)*(-sin(pi*(x(i)-half*l)/l)*pi**2/l**2/(one-(two*x(i)-l)**2/l**2)+eight*cos(pi*(x(i)-half*l)/l)*pi/l**3&
		/(one-(two*x(i)-l)**2/l**2)**2*(two*x(i)-l)+32.0_long*sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2)**3*(two*x(i)-l)&
		**2/l**4+eight*sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2)**2/l**2)
!		hxx(i) = two*(half*h0-half*h1)*tanh(sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2))*(one-tanh(sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2))**2)*(cos(pi*(x(i)-half*l)/l)*pi/l/(one-(two*x(i)-l)**2/l**2)+four*sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2)**2*(two*x(i)-l)/l**2)**2-(half*h0-half*h1)*(one-tanh(sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2))**2)*(-sin(pi*(x(i)-half*l)/l)*pi**2/l**2/(one-(two*x(i)-l)**2/l**2)+eight*cos(pi*(x(i)-half*l)/l)*pi/l**3/(one-(two*x(i)-l)**2/l**2)**2*(two*x(i)-l)+32.0_long*sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2)**3*(two*x(i)-l)**2/l**4+eight*sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2)**2/l**2)
    else
        h(i)   = h0
        hx(i)  = zero
        hxx(i) = zero
    endif
end do
i = Nx
h(i) = h0 - (h0-h1)
END SUBROUTINE BeachGen2

!
! GD: for flexibility (particularly to use relaxation zones on both ends)
!
SUBROUTINE BeachGen2_bis(h0,h1,ind1,ind2,x,h,hx,hxx,Nx)
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: i, Nx, ind1, ind2
REAL(KIND=long) :: l, x(Nx), h0, h1, h(Nx), hx(Nx), hxx(Nx), tmp
! Left part
h(1:ind1)  = h0
hx(1:ind1) = zero
hxx(1:ind1)= zero
! Right part
h(ind2:Nx)  = h1
hx(ind2:Nx) = zero
hxx(ind2:Nx)= zero
! length of the beach
l=x(ind2)-x(ind1)
DO i = ind1+1, ind2-1
  	tmp = x(i)-x(ind1)
    ! GD: Always use integer power when possible... (more precise and faster)
    h(i)   = h0 - (h0-h1)/two*(one+tanh(sin(pi*(tmp-l/two)/l)/(one-(two*(tmp-l/two)/l)**2)))
    hx(i)  = -(one/two*h0-one/two*h1)*(one-tanh(sin(pi*(tmp-one/two*l)/l)/(one-(two*tmp-l)**2/l**2))**2)*(cos(pi*(tmp&
    -one/two*l)/l)*pi/l/(one-(two*tmp-l)**2/l**2)+four*sin(pi*(tmp-one/two*l)/l)/(one-(two*tmp-l)**2/l**2)&
    **2*(two*tmp-l)/l**2)
    hxx(i) = two*(half*h0-half*h1)*tanh(sin(pi*(tmp-half*l)/l)/(one-(two*tmp-l)**2/l**2))*(one-tanh(sin(pi*(tmp-half*l)/l)&
    /(one-(two*tmp-l)**2/l**2))**2)*(cos(pi*(tmp-half*l)/l)*pi/l/(one-(two*tmp-l)**2/l**2)+four*sin(pi*(tmp-half*l)/l)&
    /(one-(two*tmp-l)**2/l**2)**2*(two*tmp-l)/l**2)**2-(half*h0-half*h1)*(one-tanh(sin(pi*(tmp-half*l)/l)/(one-(two*tmp&
    -l)**2/l**2))**2)*(-sin(pi*(tmp-half*l)/l)*pi**2/l**2/(one-(two*tmp-l)**2/l**2)+eight*cos(pi*(tmp-half*l)/l)*pi/l**3&
    /(one-(two*tmp-l)**2/l**2)**2*(two*tmp-l)+32.0_long*sin(pi*(tmp-half*l)/l)/(one-(two*tmp-l)**2/l**2)**3*(two*tmp-l)&
    **2/l**4+eight*sin(pi*(tmp-half*l)/l)/(one-(two*tmp-l)**2/l**2)**2/l**2)
!		hxx(i) = two*(half*h0-half*h1)*tanh(sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2))*(one-tanh(sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2))**2)*(cos(pi*(x(i)-half*l)/l)*pi/l/(one-(two*x(i)-l)**2/l**2)+four*sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2)**2*(two*x(i)-l)/l**2)**2-(half*h0-half*h1)*(one-tanh(sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2))**2)*(-sin(pi*(x(i)-half*l)/l)*pi**2/l**2/(one-(two*x(i)-l)**2/l**2)+eight*cos(pi*(x(i)-half*l)/l)*pi/l**3/(one-(two*x(i)-l)**2/l**2)**2*(two*x(i)-l)+32.0_long*sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2)**3*(two*x(i)-l)**2/l**4+eight*sin(pi*(x(i)-half*l)/l)/(one-(two*x(i)-l)**2/l**2)**2/l**2)
end do

END SUBROUTINE BeachGen2_bis
