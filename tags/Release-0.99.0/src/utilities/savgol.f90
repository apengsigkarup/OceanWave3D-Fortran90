SUBROUTINE savgol(c,np,nl,nr,ld,m)

! Numerical recipes routine for computing the Savitsky-Golay filter
! coefficients for smoothing the ld^th derivative of a function using
! np = nl+nr+1 points and an m^th order polynomial.  Coefficients are
! returned in c in wraparound order.
USE Precision
  IMPLICIT NONE
  INTEGER ld, m, nl, np, nr
  REAL(kind=long) :: c(np)
  INTEGER :: imj, ipj, j, k, kk, mm, indx(m+1)
  REAL(kind=long) :: d, fac, sum, a(m+1,m+1), b(m+1), zero=0._long, &
       one=1._long

  if (np<nl+nr+1 .or. nl<0 .or. nr<0 .or. ld>m .or. nl+nr<m) then
       write(*,'(A)') 'bad args in savgol'
       read*
! pause
  endif

  do ipj=0,2*m
     sum=zero
     if(ipj==0)sum=one
     do k=1,nr
        sum=sum+real(k,long)**ipj
     end do
     do k=1,nl
        sum=sum+real(-k,long)**ipj
     end do
     mm=min(ipj,2*m-ipj)
     do imj=-mm,mm,2
        a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum
     end do
  end do
  call ludcmp(a,m+1,m+1,indx,d)
  do j=1,m+1
     b(j)=zero
  end do
  b(ld+1)=one
  call lubksb(a,m+1,m+1,indx,b)
  do kk=1,np
     c(kk)=zero
  end do
  do k=-nl,nr
     sum=b(1)
     fac=one
     do mm=1,m
        fac=fac*k
        sum=sum+b(mm+1)*fac
     end do
     kk=mod(np-k,np)+1
     c(kk)=sum
  end do
  return
end SUBROUTINE savgol
