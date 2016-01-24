SUBROUTINE FilterInit(filtercoefficients,filtercoefficients2)

USE Precision
USE GlobalVariables, ONLY: tmpfilter,filterNP,filterALPHA,filterORDER,sigma_filt,fileop
USE Constants
!
! Set up the filtering coefficients.  The variables tmpfilter, filtercoefficients 
! and filtercoefficients2 have all been allocated and initialized to zero in 
! ReadInputFileParameters.f90.  
!
! filtercoefficients(1:filterNP):  The centered S-G filter
! filtercoefficients(1:filterNP,filterAlpha):  The off-centered coefficients towards a boundary 
!
! The off-centered coefficients of Berland et al JCP (2007) are applied to the
! first 3 boundary points with the weights given by sigma_filt, then centered S-G filters 
! are used up to point 6.  We assume that filterAlpha<=6.  
!
!  If sigma_filt(1)=0, then the boundary filtering is turned off.  

IMPLICIT NONE

INTEGER         :: i
REAL(KIND=long) :: filtercoefficients(filterNP),   &
     filtercoefficients2(max(filterNP,13),max(filterAlpha,6))
!
! Determine Savitzky-Golay filter
! filtercoefficients contains the centered SG-filter applied inside domain
!
CALL savgol(tmpfilter,filterNP,filterALPHA,filterALPHA,0,filterORDER)
! REORDER COEFFICIENTS
filtercoefficients(filterALPHA+1) = tmpfilter(1)
filtercoefficients(1:filterALPHA) = tmpfilter((filterALPHA+1):2:-1)
filtercoefficients((filterALPHA+2):filterNP) = tmpfilter(filterNP:(filterALPHA+2):-1)
!
! Different ways to treat the boundaries (case 4- seems to be the best choice)
! One has to adapt the filtering close to boundaries
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1- General off-centered Savitzky-Golay filter: NOT WORKING, amplification at some wavenumbers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$$$$$$ DO i=1,filterALPHA
!$$$$$$    ! Determine off-centered Savitzky-Golay filter
!$$$$$$    CALL savgol(tmpfilter,filterNP,i-1,filterNP-i,0,filterORDER)
!$$$$$$    ! REORDER COEFFICIENTS
!$$$$$$    filtercoefficients2(i,i)=tmpfilter(1)
!$$$$$$    filtercoefficients2(1:i-1,i)=tmpfilter(i:2:-1)
!$$$$$$    filtercoefficients2(i+1:filterNP,i)=tmpfilter(filterNP:i+1:-1)
!$$$$$$ ENDDO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2- Determine the coefficients assuming symetry on the boundary...
! FIXME: shall we use filtering on global variable ? ... (eg. for free propagation, not optimal)
! NOT WORKING in a general case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$$$$$$ filtercoefficients2(1:filterNP,1:filterALPHA) = zero
!$$$$$$ DO i=1,filterALPHA
!$$$$$$   filtercoefficients2(1:filterALPHA+i,i) = filtercoefficients(filterALPHA-i+2:filterNP)
!$$$$$$   DO j=1,filterALPHA-i+1
!$$$$$$     filtercoefficients2(1+j,i)=filtercoefficients2(1+j,i)+filtercoefficients(filterALPHA+1-j-i+1)
!$$$$$$   ENDDO
!$$$$$$ ENDDO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3- ! GD: optimize the filtering close to the boundaries
! We use the ghost points now, i.e. first point uses 1 point on the left and some points on the right...
! FIXME find a way to generalize the process... (plus, not working perfectly...)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$$$$$$ filterALPHA = 6
!$$$$$$ filterNP = 2*filterALPHA+1
!$$$$$$ filterORDER = 10
!$$$$$$ !
!$$$$$$ DEALLOCATE(filtercoefficients,tmpfilter,filtercoefficients2)
!$$$$$$ ALLOCATE(filtercoefficients(filterNP),tmpfilter(filterNP),filtercoefficients2(filterNP,filterALPHA))
!$$$$$$ ALLOCATE(filtercoefficients3(filterNP,filterALPHA))
!$$$$$$ filtercoefficients = zero; tmpfilter = zero; filtercoefficients2 = zero
!$$$$$$ !
!$$$$$$ ! Determine Savitzky-Golay filter : centered scheme
!$$$$$$ CALL savgol(tmpfilter,filterNP,filterALPHA,filterALPHA,0,filterORDER)
!$$$$$$ ! REORDER COEFFICIENTS
!$$$$$$ filtercoefficients(filterALPHA+1) = tmpfilter(1)
!$$$$$$ filtercoefficients(1:filterALPHA) = tmpfilter((filterALPHA+1):2:-1)
!$$$$$$ filtercoefficients((filterALPHA+2):filterNP) = tmpfilter(filterNP:(filterALPHA+2):-1)
!$$$$$$
!$$$$$$ !Off-centered scheme
!$$$$$$      ! First point no filter : test
!$$$$$$      i=1
!$$$$$$      ! First point
!$$$$$$      CALL savgol(tmpfilter,5,1,3,0,2)
!$$$$$$ !$$$$$$      CALL savgol(tmpfilter,13,1,11,0,2)
!$$$$$$      i=1
!$$$$$$      filtercoefficients2(i+1,i)=tmpfilter(1)
!$$$$$$      filtercoefficients2(1:i,i)=tmpfilter(i+1:2:-1)
!$$$$$$      filtercoefficients2(i+2:5,i)=tmpfilter(5:i+2:-1)
!$$$$$$ !$$$$$$      filtercoefficients2(i+2:13,i)=tmpfilter(13:i+2:-1)
!$$$$$$
!$$$$$$      ! Second point no filter : test
!$$$$$$      i=2
!$$$$$$      ! Second point
!$$$$$$      CALL savgol(tmpfilter,7,2,4,0,3)
!$$$$$$ !$$$$$$      CALL savgol(tmpfilter,13,2,10,0,2)
!$$$$$$      i=2
!$$$$$$      filtercoefficients2(i+1,i)=tmpfilter(1)
!$$$$$$      filtercoefficients2(1:i,i)=tmpfilter(i+1:2:-1)
!$$$$$$      filtercoefficients2(i+2:7,i)=tmpfilter(7:i+2:-1)
!$$$$$$      ! Third point
!$$$$$$      CALL savgol(tmpfilter,7,3,3,0,4)
!$$$$$$ !$$$$$$      CALL savgol(tmpfilter,13,3,9,0,2)
!$$$$$$      i=3
!$$$$$$      filtercoefficients2(i+1,i)=tmpfilter(1)
!$$$$$$      filtercoefficients2(1:i,i)=tmpfilter(i+1:2:-1)
!$$$$$$      filtercoefficients2(i+2:7,i)=tmpfilter(7:i+2:-1)
!$$$$$$ !$$$$$$      filtercoefficients2(i+2:13,i)=tmpfilter(13:i+2:-1)
!$$$$$$      ! Fourth point
!$$$$$$      CALL savgol(tmpfilter,9,4,4,0,5)
!$$$$$$ !$$$$$$      CALL savgol(tmpfilter,13,4,8,0,2)
!$$$$$$      i=4
!$$$$$$      filtercoefficients2(i+1,i)=tmpfilter(1)
!$$$$$$      filtercoefficients2(1:i,i)=tmpfilter(i+1:2:-1)
!$$$$$$      filtercoefficients2(i+2:9,i)=tmpfilter(9:i+2:-1)
!$$$$$$ !$$$$$$      filtercoefficients2(i+2:13,i)=tmpfilter(13:i+2:-1)
!$$$$$$      ! Fifth point
!$$$$$$      CALL savgol(tmpfilter,11,5,5,0,7)
!$$$$$$ !$$$$$$      CALL savgol(tmpfilter,13,5,7,0,2)
!$$$$$$      i=5
!$$$$$$      filtercoefficients2(i+1,i)=tmpfilter(1)
!$$$$$$      filtercoefficients2(1:i,i)=tmpfilter(i+1:2:-1)
!$$$$$$      filtercoefficients2(i+2:11,i)=tmpfilter(11:i+2:-1)
!$$$$$$      ! Six point and other points with centered schemes...
!$$$$$$      CALL savgol(tmpfilter,filterNP,filterALPHA,filterALPHA,0,filterORDER)
!$$$$$$      i=filterALPHA
!$$$$$$      filtercoefficients2(i+1,i)=tmpfilter(1)
!$$$$$$      filtercoefficients2(1:i,i)=tmpfilter(i+1:2:-1)
!$$$$$$      filtercoefficients2(i+2:filterNP,i)=tmpfilter(filterNP:i+2:-1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 4- other choice of filtering not taking into account ghost point...
! We use the optimized off-centered filters of Berland et al. (2007)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! other choice of filtering not taking into account ghost point...
!
! First point no filter : test
filtercoefficients2=zero
i=1
filtercoefficients2(i,i)=1.d0
! First point, filter from Berland et al. (JCP 2007)
tmpfilter=zero
tmpfilter(1) = 0.320882352941
tmpfilter(2) = -0.465
tmpfilter(3) = 0.179117647059
tmpfilter(4) = -0.035
filtercoefficients2(:,i)=filtercoefficients2(:,i)-sigma_filt(1)*tmpfilter(:)
tmpfilter=zero
! Second point : test
i=2
filtercoefficients2(i,i)=1.d0
! Second point, filter from Berland et al. (JCP 2007)
tmpfilter=zero
tmpfilter(1) = -0.085777408970
tmpfilter(2) = 0.277628171524
tmpfilter(3) = -0.356848072173
tmpfilter(4) = 0.223119093072
tmpfilter(5) = -0.057347064865
tmpfilter(6) = -0.000747264596
tmpfilter(7) = -0.000027453993
filtercoefficients2(:,i)=filtercoefficients2(:,i)-sigma_filt(2)*tmpfilter(:)
tmpfilter=zero

! Third point : test
i=3
filtercoefficients2(i,i)=1.d0
! Third point, filter from Berland et al. (JCP 2007)
tmpfilter=zero
tmpfilter(1) = 0.052523901012
tmpfilter(2) = -0.206299133811
tmpfilter(3) = 0.353527998250
tmpfilter(4) = -0.348142394842
tmpfilter(5) = 0.181481803619
tmpfilter(6) = 0.009440804370
tmpfilter(7) = -0.077675100452
tmpfilter(8) = 0.044887364863
tmpfilter(9) = -0.009971961849
tmpfilter(10)= 0.000113359420
tmpfilter(11)= 0.000113359420
filtercoefficients2(:,i)=filtercoefficients2(:,i)-sigma_filt(3)*tmpfilter(:)
tmpfilter=zero
! Fourth point
CALL savgol(tmpfilter,7,3,3,0,5)
i=4
filtercoefficients2(i,i)=tmpfilter(1)
filtercoefficients2(1:i-1,i)=tmpfilter(i:2:-1)
filtercoefficients2(i+1:7,i)=tmpfilter(7:i+1:-1)

! Fifth point
CALL savgol(tmpfilter,9,4,4,0,6)
i=5
filtercoefficients2(i,i)=tmpfilter(1)
filtercoefficients2(1:i-1,i)=tmpfilter(i:2:-1)
filtercoefficients2(i+1:9,i)=tmpfilter(9:i+1:-1)

! Six point
CALL savgol(tmpfilter,11,5,5,0,8)
i=6
filtercoefficients2(i,i)=tmpfilter(1)
filtercoefficients2(1:i-1,i)=tmpfilter(i:2:-1)
filtercoefficients2(i+1:11,i)=tmpfilter(11:i+1:-1)
!
DEALLOCATE(tmpfilter)
!
! Turn off the boundary filtering if sigma_filt(1)=0
!
!
if(abs(sigma_filt(1))<1.e-12)then
   print *, 'FilterInit:  Boundary filtering is turned off.'
   print *, ' ' 
   write(fileop(1),*) 'FilterInit:  Boundary filtering is turned off.'
   write(fileop(1),*) ' ' 
   filtercoefficients2=zero
   do i=1,max(filterALPHA,6)
      filtercoefficients2(i,i)=one
   end do
else
   print *, 'FilterInit:  Boundary filtering is turned on with values:'
   print *, sigma_filt
   print *, ' ' 
   write(fileop(1),*) 'FilterInit:  Boundary filtering is turned on with values:'
   write(fileop(1),*) sigma_filt
   write(fileop(1),*) ' ' 

end if


END SUBROUTINE FilterInit
