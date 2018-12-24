! Module containing all the mathematical functions needed for      13/12/2018
! SerraLINE
 module functions_mod

 use parms

 implicit none

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

 contains

!-----------------------------------------------------------------------
 !This subroutine is used in SerraLINE
 subroutine get_tangent_vectors(coords,N,frames,circle_str,tangents)
 
 implicit none
 
 logical, intent(in) :: circle_str
 integer, intent(in) :: N,frames
 real(dp), intent(in) :: coords(:,:,:)
 real(dp), intent(out) :: tangents(:,:,:)
 integer :: i,k

 do k=1,frames

  do i=1,N-1
    tangents(:,i,k) = coords(:,i+1,k) - coords(:,i,k)
    tangents(:,i,k) = normalize_vector(tangents(:,i,k),3)
  end do

  !If its a closed structure then the last bp has tangent vector
  if (circle_str) then
    tangents(:,N,k) = coords(:,1,k) - coords(:,N,k)
    tangents(:,N,k) = normalize_vector(tangents(:,N,k),3)
  end if

 end do

 end subroutine get_tangent_vectors
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !This subroutine calculates the bending angles between tangent vectors
 ! The output is in degrees
 subroutine get_bendings(ndim,frames,tangents,bends)
 
 implicit none
 
 integer, intent(in) :: ndim,frames
 real(dp), intent(in) :: tangents(:,:,:)
 real(dp), intent(out) :: bends(:,:)
 real(dp) :: a, b
 integer :: i, j, k, l

  do k=1,frames
    l=0
    do j=1,ndim-1
      do i=1,ndim-j
        l=l+1
        a = dot_product(tangents(:,i,k),tangents(:,i+j,k))
        b = dacos(a)
        bends(k,l) = rad_to_deg*b !Convert to degrees
      end do
    end do
  end do

 end subroutine get_bendings
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !This function normalizes vector v with dimension n
 function normalize_vector(v,n)
 implicit none
 real(dp) :: normalize_vector(n)
 real(dp) :: v(n)
 integer :: n
 integer :: i
 real(dp) :: v_sum

 v_sum=0._dp

 do i=1,n
  v_sum = v_sum + v(i)*v(i)
 end do

 v_sum=sqrt(v_sum)

 do i=1,n
  normalize_vector(i)=v(i)/v_sum
 end do

 end function 

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !Average and standard deviations are stored in the same arrays in SerraLINE
 function average_std(X,n)
 implicit none
 real(dp) :: average_std(2)
 integer :: n
 real(dp) :: X(n)
 real(dp) :: a !auxiliar variable
 integer :: i

 !Get Average
 average_std(1) = sum(X)/real(n,dp)

 !Now Standard deviation
 average_std(2)=0.0_dp

 do i=1,n
   a = X(i) - average_std(1)
   average_std(2) = average_std(2) + a*a
 end do
 average_std(2) = average_std(2)/real(n,dp)
 average_std(2) = sqrt(average_std(2))

 end function 

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Calculates the absolute value of vector a and stores it in b
 function absv(a) result (b)
 implicit none
 real(dp), intent(in) :: a(:)
 real(dp) :: b
 integer :: i
  b=0.
  do i=1,size(a)
   b=b+a(i)**2
  end do
  b=sqrt(b)
 end function absv

!-----------------------------------------------------------------------


 end module functions_mod
