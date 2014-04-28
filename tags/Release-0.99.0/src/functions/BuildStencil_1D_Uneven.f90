SUBROUTINE BuildStencil_1D_Uneven(Nx,rank,x,fx)
!
! A function to compute finite-difference coefficients for the first 
! 2*alpha derivatives of a function on a variably spaced grid x(1:N) using 
! 2*alpha+1 neighboring points.  N  sets of 
! coefficients are returned where set 1,2,...,alpha are one-sided schemes 
! at the left end, alpha+1,...,N-alpha are centered schemes, and 
! N-alpha+1,...,N  are one-sided schemes at the right end.  
! 
! By Allan P. Engsig-Karup.
  USE Precision
  IMPLICIT NONE

  INTEGER :: Nx, n, m, ip, alpha, rank,  lwork, info
  INTEGER :: ipiv(rank)
  REAL(kind=long) :: one=1._long, dxmin, fact
  REAL(kind=long) :: fx(rank,rank,Nx), x(Nx), work(5*rank*rank), dx(Nx),   &
       factorial(0:rank-1), mat(rank,rank)

  lwork=5*rank*rank
  alpha=(rank-1)/2

  factorial(0)=one
  Do n=1,rank-1
     factorial(n)=factorial(n-1)*real(n,long)
  END Do
  ! Scale each column of the matrix by the appropriate power of dx_min.
  dxMin=1.E10
  do m=1,Nx-1
     dx(m)=x(m+1)-x(m);
     dxMin=min(abs(dx(m)),dxMin)
  end do
  fact=one/dxMin
  ! One-sided schemes for the left end-points.  
  do ip=1,alpha
     do m=1,rank
        do n=1,rank
           mat(m,n)=(fact*(x(m)-x(ip)))**(n-1)/factorial(n-1)
        end do
     end do

     CALL dgetrf(rank,rank,mat,rank,ipiv,info) 
     IF(info/=0)THEN
        PRINT *, 'Problems with L-U of the finite-difference matrix',info,info
        STOP
     END IF
     CALL dgetri(rank,mat,rank,ipiv,work,lwork,info)
     IF(info/=0)THEN                    ! and invert it
        PRINT *, 'Problems with inversion of the finite-difference matrix',info,info
        STOP
     END IF
     ! Re-scale by the dx factor and re-order.  
     Do n=2,rank
        Do m=1,rank
           fx(m,n,ip)=fact**(n-1)*mat(n,m)
        END Do
     END Do
  end do
  ! The centered schemes
  do ip=alpha+1,Nx-alpha
     DO m=-alpha,alpha
        Do n=1,rank
           mat(m+alpha+1,n)=(fact*(x(ip+m)-x(ip)))**(n-1)/factorial(n-1);
        end Do
     end DO
     CALL dgetrf(rank,rank,mat,rank,ipiv,info) 
     IF(info/=0)THEN
        PRINT *, 'Problems with L-U of the finite-difference matrix',info,info
        STOP
     END IF
     CALL dgetri(rank,mat,rank,ipiv,work,lwork,info)
     IF(info/=0)THEN                    ! and invert it
        PRINT *, 'Problems with inversion of the finite-difference matrix',info,info
        STOP
     END IF
     ! Re-scale by the dx factor.  
     Do n=2,rank
        Do m=1,rank
           fx(m,n,ip)=fact**(n-1)*mat(n,m)
        END Do
     END Do
  end do

  ! One-sided schemes for the right end-points.  
  Do ip=Nx-alpha+1,Nx
     Do n=1,rank
        Do m=1,rank
           mat(m,n)=(fact*(x(Nx-rank+m)-x(ip)))**(n-1)/factorial(n-1);
        end Do
     end Do
     CALL dgetrf(rank,rank,mat,rank,ipiv,info) 
     IF(info/=0)THEN
        PRINT *, 'Problems with L-U of the finite-difference matrix',info,info
        STOP
     END IF
     CALL dgetri(rank,mat,rank,ipiv,work,lwork,info)
     IF(info/=0)THEN                    ! and invert it
        PRINT *, 'Problems with inversion of the finite-difference matrix',info,info
        STOP
     END IF
! Re-scale by the dx factor.  
     Do n=2,rank
        Do m=1,rank
           fx(m,n,ip)=fact**(n-1)*mat(n,m)
        END Do
     END Do
  end do
  

  RETURN
END SUBROUTINE BuildStencil_1D_Uneven
