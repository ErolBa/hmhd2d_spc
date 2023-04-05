#include"inc.h"
      MODULE ssub
      USE globals
      USE mp
      USE comm
      IMPLICIT NONE


! Prevent cray fortran compiler from "cloning" these subroutines
! Cloning these routines may cause segmentation fault in cce/8.6.x
! This may be fixed in later versions.

!DIR$ CLONENEVER bcxper, bcxsym, bcxasym, bcxext, bcx      
!DIR$ CLONENEVER bczper, bczsym, bczasym, bczext, bcz      
     
      
      CONTAINS
!**********************************************
      subroutine bcxper(dd,tag)
      REAL(RTYPE) :: dd(nx,nz)
      INTEGER :: tag

      if (nxproc.eq.1) then
        dd(1:2,:)=dd(nx-3:nx-2,:)
        dd(nx-1:nx,:)=dd(3:4,:)
      else
        call x_exchange(dd,tag)
      end if 
      return

      end subroutine bcxper
!**********************************************

!**********************************************
      subroutine bcxsym(dd,tag)

      REAL(RTYPE) :: dd(nx,nz)
      INTEGER :: tag
      INTEGER :: k
      if (nxproc.eq.1) then
        do k=1,nz
          dd(2,k)=dd(4,k)
          dd(1,k)=dd(5,k)
          dd(nx-1,k)=dd(nx-3,k)
          dd(nx-0,k)=dd(nx-4,k)
        end do
      elseif (iproc_x.eq.0) then
        call x_exchange_right(dd,tag)
        do k=1,nz
          dd(2,k)=dd(4,k)
          dd(1,k)=dd(5,k)
        end do
      elseif (iproc_x.eq.nxproc-1) then
        call x_exchange_left(dd,tag)
        do k=1,nz
          dd(nx-1,k)=dd(nx-3,k)
          dd(nx-0,k)=dd(nx-4,k)
        end do
      else
        call x_exchange(dd,tag)
      end if
      return
      end subroutine bcxsym
!**********************************************
!**********************************************

      subroutine bcxasym(dd,tag)

      REAL(RTYPE) :: dd(nx,nz)
      INTEGER :: tag

      INTEGER :: k   

      if (nxproc.eq.1) then 
        do k=1,nz
          dd(3,k)=zero
          dd(2,k)=-dd(4,k)
          dd(1,k)=-dd(5,k)
          dd(nx-2,k)=zero
          dd(nx-1,k)=-dd(nx-3,k)
          dd(nx-0,k)=-dd(nx-4,k)
        end do
      elseif (iproc_x.eq.0) then
        call x_exchange_right(dd,tag)
        do k=1,nz
          dd(3,k)=zero
          dd(2,k)=-dd(3,k)
          dd(1,k)=-dd(4,k)
        end do
      elseif (iproc_x.eq.nxproc-1) then
        call x_exchange_left(dd,tag) 
        do k=1,nz
          dd(nx-2,k)=zero
          dd(nx-1,k)=-dd(nx-2,k)
          dd(nx-0,k)=-dd(nx-3,k)
        end do
      else
        call x_exchange(dd,tag)
      end if
      return
      end subroutine bcxasym
!**********************************************
	
!**********************************************
      subroutine bcxext(dd,tag)

      REAL(RTYPE) :: dd(nx,nz)
      INTEGER :: tag

      INTEGER :: k   

      if (nxproc.eq.1) then 
        do k=1,nz
          dd(2,k)=dd(3,k)+(x(2)-x(3))*(dd(3,k)-dd(4,k))/(x(3)-x(4))
          dd(1,k)=dd(2,k)+(x(1)-x(2))*(dd(2,k)-dd(3,k))/(x(2)-x(3))
          dd(nx-1,k)=dd(nx-2,k) + (x(nx-1)-x(nx-2))  &
                                 *(dd(nx-2,k)-dd(nx-3,k)) &
                                 /(x(nx-2)-x(nx-3))
          dd(nx-0,k)=dd(nx-1,k) + (x(nx-0)-x(nx-1)) &
                                 *(dd(nx-1,k)-dd(nx-2,k)) &
                                 /(x(nx-1)-x(nx-2))
        end do
      elseif (iproc_x.eq.0) then
        call x_exchange_right(dd,tag)
        do k=1,nz
          dd(2,k)=dd(3,k)+(x(2)-x(3))*(dd(3,k)-dd(4,k))/(x(3)-x(4))
          dd(1,k)=dd(2,k)+(x(1)-x(2))*(dd(2,k)-dd(3,k))/(x(2)-x(3))
        end do
      elseif (iproc_x.eq.nxproc-1) then
        call x_exchange_left(dd,tag)
        do k=1,nz
          dd(nx-1,k)=dd(nx-2,k) + (x(nx-1)-x(nx-2))  &
                                 *(dd(nx-2,k)-dd(nx-3,k)) &
                                 /(x(nx-2)-x(nx-3))
          dd(nx-0,k)=dd(nx-1,k) + (x(nx-0)-x(nx-1)) &
                                 *(dd(nx-1,k)-dd(nx-2,k)) &
                                 /(x(nx-1)-x(nx-2))
        end do
      else
        call x_exchange(dd,tag)
      end if 
      return
      end subroutine bcxext
!**********************************************

!**********************************************
      subroutine bcx(dd, left, right, tag)

      REAL(RTYPE) :: dd(nx,nz)
      INTEGER :: tag

      CHARACTER(*) :: left, right
      REAL(RTYPE) :: d1, d2, d3, d4, d5, d6, d7, d8
      REAL(RTYPE) :: dm1, dm2, dp1, dp2, dp3

      INTEGER :: k   

! allow different boundary conditions on both sides
! '': exchange interior ghost cells only
! 'symh'  : symmetric, with boundary on half grid point
! 'symf'  : symmetric, with boundary on grid point
! 'asymh' : anti-symmetric, boundary on half grid point
! 'asymf' : anti-symmetric, with boundary on grid point
! 'ext'   : extrapolation
! 'exth1' : alternative extrapolation, such that d2dx=0, with boundary 
!            on half grid point
! 'extf1' : alternative extrapolation, such that d2dx=0, with boundary
!           on grid point      
! 'ext2'  : yet another extrapolation, use ddx and d2dx and 
!           Taylor expansion
! 'ext2_sym'  : use dd(4) and dd(5) to determine dd(3), dd(2), and
!               dd(1), assuming dd'(3)== and symmetry. 
! 'ext3'  : yet another extrapolation, use last 3 points to fit
!           a second order polynomial
! 'ext3_zero' : same as 'ext3' but set last point to zero
! 'ext3_zero_der' : 3rd order polynomial extrapolation, assuming ddx=0 at
!                  the boundary
! 'ext3_zero_der2' : 3rd order polynomial extrapolation, assuming d2dx2=0 at
!                  the boundary
! 'ext3_odd' : 3rd order polynomial extrapolation, assuming d2dx2=0 at
!                     the boundary and setting the last point to zero
! 'ext3_even' : 4th order polynomial extrapolation, assuming dxx=d3dx3=0 
!               the boundary 
! 'ext4'  : use last 4 points to fit a third order polynomial
! 'ext4_zero' : same as 'ext4' but set last point to zero

        
 
      if (iproc_x.eq.0) then
        if (nxproc>1) then 
          call x_exchange_right(dd,tag)
        end if   
        select case (left)
          case ('')
            continue
          case ('symh')
            do k=1,nz
              dd(2,k)=dd(3,k)
              dd(1,k)=dd(4,k)
            end do 
          case ('symf')
            do k=1,nz
              dd(2,k)=dd(4,k)
              dd(1,k)=dd(5,k)
            end do
          case ('asymh')
            do k=1,nz
              dd(2,k)=-dd(3,k)
              dd(1,k)=-dd(4,k)
            end do 
          case ('asymf')
            do k=1,nz
              dd(3,k)=zero    
              dd(2,k)=-dd(4,k)
              dd(1,k)=-dd(5,k)
            end do  
          case ('ext')
            do k=1,nz
              dd(2,k)=dd(3,k)+(x(2)-x(3))*(dd(3,k)-dd(4,k))/(x(3)-x(4))
              dd(1,k)=dd(2,k)+(x(1)-x(2))*(dd(2,k)-dd(3,k))/(x(2)-x(3))
            end do
          case ('exth1')
            do k=1,nz
              dd(2,k)=dd(3,k)
              dd(1,k)=two*dd(3,k)-dd(4,k)
            end do
          case ('extf1')
            do k=1,nz
              dd(2,k)=two*dd(3,k)-dd(4,k)
              dd(1,k)=two*dd(3,k)-dd(5,k)
            end do
          case ('ext2')
            do k=1,nz
              d1=ddx(dd,5,k)
              d2=d2dx(dd,5,k)
              d3=x(2)-x(5)
              dd(2,k)=dd(5,k)+d3*(d1+half*d2*d3)
              d3=x(1)-x(5)
              dd(1,k)=dd(5,k)+d3*(d1+half*d2*d3)
            end do
          case ('ext2_sym')
            dp1 = x(4)-x(3)
            dp2 = x(5)-x(3)
            dm1 = x(2)-x(3)
            dm2 = x(1)-x(3)
            d1 =  dp2**2/(dp2**2-dp1**2)
            d2 = -dp1**2/(dp2**2-dp1**2)
            do k=1,nz
              dd(3,k) = d1*dd(4,k) + d2*dd(5,k)
              d3 = (dd(5,k) - dd(4,k))/(dp2**2-dp1**2)
              dd(2,k) = dd(3,k) + d3*dm1**2
              dd(1,k) = dd(3,k) + d3*dm2**2
            end do
          case ('ext3')
            dp1 = x(4)-x(3)
            dp2 = x(5)-x(3)
            dm1 = x(2)-x(3)
            dm2 = x(1)-x(3)
            d2 = dm1*(dm1-dp2)/(dp1*(dp1-dp2))
            d3 = dm1*(dp1-dm1)/(dp2*(dp1-dp2))
            d5 = dm2*(dm2-dp2)/(dp1*(dp1-dp2))
            d6 = dm2*(dp1-dm2)/(dp2*(dp1-dp2))
            d1= one - d2 - d3
            d4= one - d5 - d6
            do k=1,nz
              dd(2,k)=dd(3,k)*d1+dd(4,k)*d2+dd(5,k)*d3  
              dd(1,k)=dd(3,k)*d4+dd(4,k)*d5+dd(5,k)*d6
            end do
          case ('ext3_zero')
            dp1 = x(4)-x(3)
            dp2 = x(5)-x(3)
            dm1 = x(2)-x(3)
            dm2 = x(1)-x(3)
            d2 = dm1*(dm1-dp2)/(dp1*(dp1-dp2))
            d3 = dm1*(dp1-dm1)/(dp2*(dp1-dp2))
            d5 = dm2*(dm2-dp2)/(dp1*(dp1-dp2))
            d6 = dm2*(dp1-dm2)/(dp2*(dp1-dp2))
            do k=1,nz
              dd(3,k) = zero
              dd(2,k) = dd(4,k)*d2+dd(5,k)*d3  
              dd(1,k) = dd(4,k)*d5+dd(5,k)*d6
            end do
          case ('ext3_zero_der')
            dp1 = x(4)-x(3)
            dp2 = x(5)-x(3)
            dm1 = x(2)-x(3)
            dm2 = x(1)-x(3)
            d2 = (dm1/dp1)**2*(dm1-dp2)/(dp1-dp2)
            d3 = (dm1/dp2)**2*(dp1-dm1)/(dp1-dp2)
            d5 = (dm2/dp1)**2*(dm2-dp2)/(dp1-dp2)
            d6 = (dm2/dp2)**2*(dp1-dm2)/(dp1-dp2)
            d1= one - d2 - d3
            d4= one - d5 - d6
            
            do k=1,nz
              dd(2,k)=dd(3,k)*d1+dd(4,k)*d2+dd(5,k)*d3  
              dd(1,k)=dd(3,k)*d4+dd(4,k)*d5+dd(5,k)*d6
            end do
          case ('ext3_zero_der2')
            dp1 = x(4)-x(3)
            dp2 = x(5)-x(3)
            dm1 = x(2)-x(3)
            dm2 = x(1)-x(3)
            d2 = dm1*(dm1*dm1 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d3 = dm1*(dp1*dp1 - dm1*dm1)/(dp2*(dp1*dp1 - dp2*dp2))
            d5 = dm2*(dm2*dm2 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d6 = dm2*(dp1*dp1 - dm2*dm2)/(dp2*(dp1*dp1 - dp2*dp2))
            d1= one - d2 - d3
            d4= one - d5 - d6
            
            do k=1,nz
              dd(2,k)=dd(3,k)*d1+dd(4,k)*d2+dd(5,k)*d3  
              dd(1,k)=dd(3,k)*d4+dd(4,k)*d5+dd(5,k)*d6
            end do
          case ('ext3_odd')
            dp1 = x(4)-x(3)
            dp2 = x(5)-x(3)
            dm1 = x(2)-x(3)
            dm2 = x(1)-x(3)
            d2 = dm1*(dm1*dm1 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d3 = dm1*(dp1*dp1 - dm1*dm1)/(dp2*(dp1*dp1 - dp2*dp2))
            d5 = dm2*(dm2*dm2 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d6 = dm2*(dp1*dp1 - dm2*dm2)/(dp2*(dp1*dp1 - dp2*dp2))
            
            do k=1,nz
              dd(3,k)=zero  
              dd(2,k)=dd(4,k)*d2+dd(5,k)*d3  
              dd(1,k)=dd(4,k)*d5+dd(5,k)*d6
            end do
          case ('ext3_even')
            dp1 = x(4)-x(3)
            dp2 = x(5)-x(3)
            dm1 = x(2)-x(3)
            dm2 = x(1)-x(3)
            d2 = dm1*dm1*(dm1*dm1 - dp2*dp2)/(dp1*dp1*(dp1*dp1 - dp2*dp2))
            d3 = dm1*dm1*(dp1*dp1 - dm1*dm1)/(dp2*dp2*(dp1*dp1 - dp2*dp2))
            d5 = dm2*dm2*(dm2*dm2 - dp2*dp2)/(dp1*dp1*(dp1*dp1 - dp2*dp2))
            d6 = dm2*dm2*(dp1*dp1 - dm2*dm2)/(dp2*dp2*(dp1*dp1 - dp2*dp2))
            d1= one - d2 - d3
            d4= one - d5 - d6
            
            do k=1,nz
              dd(2,k)=dd(3,k)*d1+dd(4,k)*d2+dd(5,k)*d3  
              dd(1,k)=dd(3,k)*d4+dd(4,k)*d5+dd(5,k)*d6
            end do

         case ('ext4')
            dp1 = x(4)-x(3)
            dp2 = x(5)-x(3)
            dp3 = x(6)-x(3)
            dm1 = x(2)-x(3)
            dm2 = x(1)-x(3)
            d2 = dm1*(dp2*dp3+dm1*(dm1-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d3 = dm1*(dm1*(dp3+dp1-dm1)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d4 = dm1*(dp1*dp2+dm1*(dm1-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))
            d6 = dm2*(dp2*dp3+dm2*(dm2-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d7 = dm2*(dm2*(dp3+dp1-dm2)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d8 = dm2*(dp1*dp2+dm2*(dm2-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))

            d1= one - d2 - d3 - d4
            d5= one - d6 - d7 - d8 

            do k=1,nz
              dd(2,k)=dd(3,k)*d1+dd(4,k)*d2+dd(5,k)*d3+dd(6,k)*d4  
              dd(1,k)=dd(3,k)*d5+dd(4,k)*d6+dd(5,k)*d7+dd(6,k)*d8
            end do

         case ('ext4_zero')
            dp1 = x(4)-x(3)
            dp2 = x(5)-x(3)
            dp3 = x(6)-x(3)
            dm1 = x(2)-x(3)
            dm2 = x(1)-x(3)
            d2 = dm1*(dp2*dp3+dm1*(dm1-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d3 = dm1*(dm1*(dp3+dp1-dm1)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d4 = dm1*(dp1*dp2+dm1*(dm1-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))
            d6 = dm2*(dp2*dp3+dm2*(dm2-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d7 = dm2*(dm2*(dp3+dp1-dm2)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d8 = dm2*(dp1*dp2+dm2*(dm2-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))

            do k=1,nz
              dd(3,k)=zero 
              dd(2,k)=dd(4,k)*d2+dd(5,k)*d3+dd(6,k)*d4  
              dd(1,k)=dd(4,k)*d6+dd(5,k)*d7+dd(6,k)*d8
            end do

          case default
            print*,'Invalid left boundary condition in x: ', left
            call exit(1)
        end select
      end if

      if (iproc_x.eq.nxproc-1) then
        if (nxproc>1) then
          call x_exchange_left(dd,tag)
        end if  
        select case (right)
          case ('')
            continue
          case ('symh')
            do k=1,nz
              dd(nx-1,k)=dd(nx-2,k)
              dd(nx-0,k)=dd(nx-3,k)
            end do
          case ('symf')
            do k=1,nz
              dd(nx-1,k)=dd(nx-3,k)
              dd(nx-0,k)=dd(nx-4,k)
            end do
          case ('asymh')
            do k=1,nz
              dd(nx-1,k)=-dd(nx-2,k)
              dd(nx-0,k)=-dd(nx-3,k)
            end do
          case ('asymf')
            do k=1,nz
              dd(nx-2,k)=zero    
              dd(nx-1,k)=-dd(nx-3,k)
              dd(nx-0,k)=-dd(nx-4,k)
            end do
          case ('ext')
            do k=1,nz
               dd(nx-1,k)=dd(nx-2,k) + (x(nx-1)-x(nx-2))  &
                                      *(dd(nx-2,k)-dd(nx-3,k)) &
                                      /(x(nx-2)-x(nx-3))
               dd(nx-0,k)=dd(nx-1,k) + (x(nx-0)-x(nx-1)) &
                                      *(dd(nx-1,k)-dd(nx-2,k)) &
                                      /(x(nx-1)-x(nx-2))
            end do
          case ('exth1')
            do k=1,nz
              dd(nx-1,k)=dd(nx-2,k)
              dd(nx-0,k)=two*dd(nx-2,k)-dd(nx-3,k)
            end do
          case ('extf1')
            do k=1,nz
              dd(nx-1,k)=two*dd(nx-2,k)-dd(nx-3,k)
              dd(nx-0,k)=two*dd(nx-2,k)-dd(nx-4,k)
            end do
          case ('ext2')
            do k=1,nz
              d1=ddx(dd,nx-4,k)
              d2=d2dx(dd,nx-4,k)
              d3=x(nx-1)-x(nx-4)
              dd(nx-1,k)=dd(nx-4,k)+d3*(d1+half*d2*d3)
              d3=x(nx)-x(nx-4)
              dd(nx  ,k)=dd(nx-4,k)+d3*(d1+half*d2*d3)
            end do
          case ('ext2_sym')
            dp1 = x(nx-3)-x(nx-2)
            dp2 = x(nx-4)-x(nx-2)
            dm1 = x(nx-1)-x(nx-2)
            dm2 = x(nx-0)-x(nx-2)
            d1 =  dp2**2/(dp2**2-dp1**2)
            d2 = -dp1**2/(dp2**2-dp1**2)
            do k=1,nz
              dd(nx-2,k) = d1*dd(nx-3,k) + d2*dd(nx-4,k)
              d3 = (dd(nx-4,k) - dd(nx-3,k))/(dp2**2-dp1**2)
              dd(nx-1,k) = dd(nx-2,k) + d3*dm1**2
              dd(nx  ,k) = dd(nx-2,k) + d3*dm2**2
            end do
          case ('ext3')
            dp1 = x(nx-3)-x(nx-2)
            dp2 = x(nx-4)-x(nx-2)
            dm1 = x(nx-1)-x(nx-2)
            dm2 = x(nx-0)-x(nx-2)
            d2 = dm1*(dm1-dp2)/(dp1*(dp1-dp2))
            d3 = dm1*(dp1-dm1)/(dp2*(dp1-dp2))
            d5 = dm2*(dm2-dp2)/(dp1*(dp1-dp2))
            d6 = dm2*(dp1-dm2)/(dp2*(dp1-dp2))
            d1= one - d2 - d3
            d4= one - d5 - d6
            
            do k=1,nz
              dd(nx-1,k)=dd(nx-2,k)*d1+dd(nx-3,k)*d2+dd(nx-4,k)*d3  
              dd(nx-0,k)=dd(nx-2,k)*d4+dd(nx-3,k)*d5+dd(nx-4,k)*d6
            end do
          case ('ext3_zero')
            dp1 = x(nx-3)-x(nx-2)
            dp2 = x(nx-4)-x(nx-2)
            dm1 = x(nx-1)-x(nx-2)
            dm2 = x(nx-0)-x(nx-2)
            d2 = dm1*(dm1-dp2)/(dp1*(dp1-dp2))
            d3 = dm1*(dp1-dm1)/(dp2*(dp1-dp2))
            d5 = dm2*(dm2-dp2)/(dp1*(dp1-dp2))
            d6 = dm2*(dp1-dm2)/(dp2*(dp1-dp2))
            do k=1,nz
              dd(nx-2,k)=zero                         
              dd(nx-1,k)=dd(nx-3,k)*d2+dd(nx-4,k)*d3  
              dd(nx-0,k)=dd(nx-3,k)*d5+dd(nx-4,k)*d6
            end do
          case ('ext3_zero_der')
            dp1 = x(nx-3)-x(nx-2)
            dp2 = x(nx-4)-x(nx-2)
            dm1 = x(nx-1)-x(nx-2)
            dm2 = x(nx-0)-x(nx-2)
            d2 = (dm1/dp1)**2*(dm1-dp2)/(dp1-dp2)
            d3 = (dm1/dp2)**2*(dp1-dm1)/(dp1-dp2)
            d5 = (dm2/dp1)**2*(dm2-dp2)/(dp1-dp2)
            d6 = (dm2/dp2)**2*(dp1-dm2)/(dp1-dp2)
            d1= one - d2 - d3
            d4= one - d5 - d6
            do k=1,nz
              dd(nx-1,k)=dd(nx-2,k)*d1+dd(nx-3,k)*d2+dd(nx-4,k)*d3  
              dd(nx-0,k)=dd(nx-2,k)*d4+dd(nx-3,k)*d5+dd(nx-4,k)*d6
            end do
          case ('ext3_zero_der2')
            dp1 = x(nx-3)-x(nx-2)
            dp2 = x(nx-4)-x(nx-2)
            dm1 = x(nx-1)-x(nx-2)
            dm2 = x(nx-0)-x(nx-2)
            d2 = dm1*(dm1*dm1 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d3 = dm1*(dp1*dp1 - dm1*dm1)/(dp2*(dp1*dp1 - dp2*dp2))
            d5 = dm2*(dm2*dm2 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d6 = dm2*(dp1*dp1 - dm2*dm2)/(dp2*(dp1*dp1 - dp2*dp2))
            d1= one - d2 - d3
            d4= one - d5 - d6

            do k=1,nz
              dd(nx-1,k)=dd(nx-2,k)*d1+dd(nx-3,k)*d2+dd(nx-4,k)*d3  
              dd(nx-0,k)=dd(nx-2,k)*d4+dd(nx-3,k)*d5+dd(nx-4,k)*d6
            end do
          case ('ext3_odd')
            dp1 = x(nx-3)-x(nx-2)
            dp2 = x(nx-4)-x(nx-2)
            dm1 = x(nx-1)-x(nx-2)
            dm2 = x(nx-0)-x(nx-2)
            d2 = dm1*(dm1*dm1 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d3 = dm1*(dp1*dp1 - dm1*dm1)/(dp2*(dp1*dp1 - dp2*dp2))
            d5 = dm2*(dm2*dm2 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d6 = dm2*(dp1*dp1 - dm2*dm2)/(dp2*(dp1*dp1 - dp2*dp2))

            do k=1,nz
              dd(nx-2,k)=zero  
              dd(nx-1,k)=dd(nx-3,k)*d2+dd(nx-4,k)*d3  
              dd(nx-0,k)=dd(nx-3,k)*d5+dd(nx-4,k)*d6
            end do

          case ('ext3_even')
            dp1 = x(nx-3)-x(nx-2)
            dp2 = x(nx-4)-x(nx-2)
            dm1 = x(nx-1)-x(nx-2)
            dm2 = x(nx-0)-x(nx-2)
            d2 = dm1*dm1*(dm1*dm1 - dp2*dp2)/(dp1*dp1*(dp1*dp1 - dp2*dp2))
            d3 = dm1*dm1*(dp1*dp1 - dm1*dm1)/(dp2*dp2*(dp1*dp1 - dp2*dp2))
            d5 = dm2*dm2*(dm2*dm2 - dp2*dp2)/(dp1*dp1*(dp1*dp1 - dp2*dp2))
            d6 = dm2*dm2*(dp1*dp1 - dm2*dm2)/(dp2*dp2*(dp1*dp1 - dp2*dp2))
            d1= one - d2 - d3
            d4= one - d5 - d6

            do k=1,nz
              dd(nx-1,k)=dd(nx-2,k)*d1+dd(nx-3,k)*d2+dd(nx-4,k)*d3  
              dd(nx-0,k)=dd(nx-2,k)*d4+dd(nx-3,k)*d5+dd(nx-4,k)*d6
            end do

          case ('ext4')
            dp1 = x(nx-3)-x(nx-2)
            dp2 = x(nx-4)-x(nx-2)
            dp3 = x(nx-5)-x(nx-2)
            dm1 = x(nx-1)-x(nx-2)
            dm2 = x(nx-0)-x(nx-2)
            
            d2 = dm1*(dp2*dp3+dm1*(dm1-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d3 = dm1*(dm1*(dp3+dp1-dm1)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d4 = dm1*(dp1*dp2+dm1*(dm1-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))
            d6 = dm2*(dp2*dp3+dm2*(dm2-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d7 = dm2*(dm2*(dp3+dp1-dm2)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d8 = dm2*(dp1*dp2+dm2*(dm2-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))

            d1= one - d2 - d3 - d4
            d5= one - d6 - d7 - d8 
            
            do k=1,nz
              dd(nx-1,k)= dd(nx-2,k)*d1+dd(nx-3,k)*d2 &
                         +dd(nx-4,k)*d3+dd(nx-5,k)*d4  
              dd(nx-0,k)= dd(nx-2,k)*d5+dd(nx-3,k)*d6 &
                         +dd(nx-4,k)*d7+dd(nx-5,k)*d8
            end do
 
          case ('ext4_zero')
            dp1 = x(nx-3)-x(nx-2)
            dp2 = x(nx-4)-x(nx-2)
            dp3 = x(nx-5)-x(nx-2)
            dm1 = x(nx-1)-x(nx-2)
            dm2 = x(nx-0)-x(nx-2)
            
            d2 = dm1*(dp2*dp3+dm1*(dm1-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d3 = dm1*(dm1*(dp3+dp1-dm1)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d4 = dm1*(dp1*dp2+dm1*(dm1-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))
            d6 = dm2*(dp2*dp3+dm2*(dm2-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d7 = dm2*(dm2*(dp3+dp1-dm2)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d8 = dm2*(dp1*dp2+dm2*(dm2-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))

            do k=1,nz
              dd(nx-2,k)= zero 
              dd(nx-1,k)= dd(nx-3,k)*d2 &
                         +dd(nx-4,k)*d3+dd(nx-5,k)*d4  
              dd(nx-0,k)= dd(nx-3,k)*d6 &
                         +dd(nx-4,k)*d7+dd(nx-5,k)*d8
            end do
           
          case default
            print*,'Invalid right boundary condition in x: ', right
            call exit(1)
        end select
      end if  
      if (iproc_x>0.and.iproc_x<nxproc-1) then
        call x_exchange(dd,tag)
      end if 
      return
      end subroutine bcx
!**********************************************

!**********************************************
      subroutine bczper(dd,tag)

      REAL(RTYPE) :: dd(nx,nz)
      INTEGER :: tag


      if (nzproc.eq.1) then
        dd(:,1:2)=dd(:,nz-3:nz-2) 
        dd(:,nz-1:nz)=dd(:,3:4)
      else
        call z_exchange(dd,tag)
      end if 
      return	
      end subroutine bczper
!**********************************************

!**********************************************
      subroutine bczsym(dd,tag)

      REAL(RTYPE) :: dd(nx,nz)
      INTEGER :: tag

      INTEGER :: i   

      if (nzproc.eq.1) then
        do i=1,nx
          dd(i,2)=dd(i,4)
          dd(i,1)=dd(i,5)
          dd(i,nz-1)=dd(i,nz-3)
          dd(i,nz-0)=dd(i,nz-4)
        end do
      elseif (iproc_z.eq.0) then
        call z_exchange_right(dd,tag)
        do i=1,nx
          dd(i,2)=dd(i,4)
          dd(i,1)=dd(i,5)
        end do
      elseif (iproc_z.eq.nzproc-1) then
        call z_exchange_left(dd,tag)
        do i=1,nx
          dd(i,nz-1)=dd(i,nz-3)
          dd(i,nz-0)=dd(i,nz-4)
        end do
      else
        call z_exchange(dd,tag)
      end if
      return
      end subroutine bczsym
!**********************************************

!**********************************************
      subroutine bczasym(dd,tag)
      REAL(RTYPE) :: dd(nx,nz)
      INTEGER :: tag

      INTEGER :: i   

      if (nzproc.eq.1) then 
        do i=1,nx
          dd(i,3)=zero
          dd(i,2)=-dd(i,4)
          dd(i,1)=-dd(i,5)
          dd(i,nz-2)=zero
          dd(i,nz-1)=-dd(i,nz-3)
          dd(i,nz-0)=-dd(i,nz-4)
        end do
      elseif (iproc_z.eq.0) then
        call z_exchange_right(dd,tag)
        do i=1,nx
          dd(i,3)=zero
          dd(i,2)=-dd(i,4)
          dd(i,1)=-dd(i,5)
        end do
      elseif (iproc_z.eq.nzproc-1) then
        call z_exchange_left(dd,tag)
        do i=1,nx
          dd(i,nz-2)=zero
          dd(i,nz-1)=-dd(i,nz-3)
          dd(i,nz-0)=-dd(i,nz-4)
        end do
      else
        call z_exchange(dd,tag)
      end if
      return
      end subroutine bczasym
!**********************************************


!**********************************************
      subroutine bczext(dd,tag)
      REAL(RTYPE) :: dd(nx,nz)
      INTEGER :: tag

      INTEGER :: i  
      if (nzproc.eq.1) then
        do i=1,nx
          dd(i,2)=dd(i,3)+(z(2)-z(3))*(dd(i,3)-dd(i,4))/(z(3)-z(4))
          dd(i,1)=dd(i,2)+(z(1)-z(2))*(dd(i,2)-dd(i,3))/(z(2)-z(3))
          dd(i,nz-1)=dd(i,nz-2) + (z(nz-1)-z(nz-2))  &
                                 *(dd(i,nz-2)-dd(i,nz-3)) &
                                 /(z(nz-2)-z(nz-3))
          dd(i,nz-0)=dd(i,nz-1) + (z(nz-0)-z(nz-1)) &
                                 *(dd(i,nz-1)-dd(i,nz-2)) &
                                 /(z(nz-1)-z(nz-2))
        end do
      elseif (iproc_z.eq.0) then
        call z_exchange_right(dd,tag)
        do i=1,nx
          dd(i,2)=dd(i,3)+(z(2)-z(3))*(dd(i,3)-dd(i,4))/(z(3)-z(4))
          dd(i,1)=dd(i,2)+(z(1)-z(2))*(dd(i,2)-dd(i,3))/(z(2)-z(3))
        end do
      elseif (iproc_z.eq.nzproc-1) then
        call z_exchange_left(dd,tag)
        do i=1,nx
          dd(i,nz-1)=dd(i,nz-2) + (z(nz-1)-z(nz-2))  &
                                 *(dd(i,nz-2)-dd(i,nz-3)) &
                                 /(z(nz-2)-z(nz-3))
          dd(i,nz-0)=dd(i,nz-1) + (z(nz-0)-z(nz-1)) &
                                 *(dd(i,nz-1)-dd(i,nz-2)) &
                                 /(z(nz-1)-z(nz-2))
        end do
      else
        call z_exchange(dd,tag)
      end if
      return
      end subroutine bczext
!**********************************************

!**********************************************
      subroutine bcz(dd, left, right, tag)

      REAL(RTYPE) :: dd(nx,nz)
      CHARACTER(*) :: left, right
      INTEGER :: tag

      INTEGER :: i   
      REAL(RTYPE) :: d1, d2, d3, d4, d5, d6, d7, d8
      REAL(RTYPE) :: dm1, dm2, dp1, dp2, dp3

! allow different boundary conditions on both sides
! '': exchange interior ghost cells only
! 'symh'  : symmetric, with boundary on half grid point
! 'symf'  : symmetric, with boundary on grid point
! 'asymh' : anti-symmetric, boundary on half grid point
! 'asymf' : anti-symmetric, with boundary on grid point
! 'ext'   : extrapolation
! 'exth1' : alternative extrapolation, such that d2dx=0, with boundary 
!            on half grid point
! 'extf1' : alternative extrapolation, such that d2dx=0, with boundary
!           on grid point      
! 'ext2'  : yet another extrapolation, use ddx and d2dx and 
!           Taylor expansion
! 'ext2_sym'  : use dd(4) and dd(5) to determine dd(3), dd(2), and
!               dd(1), assuming dd'(3)== and symmetry. 
! 'ext3'  : yet another extrapolation, use last 3 points to fit
!           a second order polynomial
! 'ext3_zero' : same as 'ext3' but set last point to zero
! 'ext3_zero_der' : 3rd order polynomial extrapolation, assuming ddx=0 at
!                  the boundary
! 'ext3_zero_der2' : 3rd order polynomial extrapolation, assuming d2dx2=0 at
!                  the boundary
! 'ext3_odd' : 3rd order polynomial extrapolation, assuming d2dx2=0 at
!                     the boundary and setting the last point to zero
! 'ext3_even' : 4th order polynomial extrapolation, assuming dxx=d3dx3=0 
!               the boundary 
! 'ext4'  : use last 4 points to fit a third order polynomial
! 'ext4_zero' : same as 'ext4' but set last point to zero
 

      if (iproc_z.eq.0) then
        if (nzproc.gt.1) then
          call z_exchange_right(dd,tag)
        end if   
        select case (left)
          case ('')
            continue
          case ('symh')
            do i=1,nx
              dd(i,2)=dd(i,3)
              dd(i,1)=dd(i,4)
            end do 
          case ('symf')
            do i=1,nx
              dd(i,2)=dd(i,4)
              dd(i,1)=dd(i,5)
            end do
          case ('asymh')
            do i=1,nx
              dd(i,2)=-dd(i,3)
              dd(i,1)=-dd(i,4)
            end do 
          case ('asymf')
            do i=1,nx
              dd(i,3)=zero
              dd(i,2)=-dd(i,4)
              dd(i,1)=-dd(i,5)
            end do 
          case ('ext')
            do i=1,nx
              dd(i,2)=dd(i,3)+(z(2)-z(3))*(dd(i,3)-dd(i,4))/(z(3)-z(4))
              dd(i,1)=dd(i,2)+(z(1)-z(2))*(dd(i,2)-dd(i,3))/(z(2)-z(3))
            end do
          case ('exth1')
            do i=1,nx
              dd(i,2)=dd(i,3)
              dd(i,1)=two*dd(i,3)-dd(i,4)
            end do
          case ('extf1')
            do i=1,nx
              dd(i,2)=two*dd(i,3)-dd(i,4)
              dd(i,1)=two*dd(i,3)-dd(i,5)
            end do
          case ('ext2')
            do i=1,nx
              d1=ddz(dd,i,5)
              d2=d2dz(dd,i,5)
              d3=z(2)-z(5)
              dd(i,2)=dd(i,5)+d3*(d1+half*d2*d3)
              d3=z(1)-z(5)
              dd(i,1)=dd(i,5)+d3*(d1+half*d2*d3)
            end do
          case ('ext2_sym')
            dp1 = z(4)-z(3)
            dp2 = z(5)-z(3)
            dm1 = z(2)-z(3)
            dm2 = z(1)-z(3)
            d1 =  dp2**2/(dp2**2-dp1**2)
            d2 = -dp1**2/(dp2**2-dp1**2)
            do i=1,nx
              dd(i,3) = d1*dd(i,4) + d2*dd(i,5)
              d3 = (dd(i,5) - dd(i,4))/(dp2**2-dp1**2)
              dd(i,2) = dd(i,3) + d3*dm1**2
              dd(i,1) = dd(i,3) + d3*dm2**2
            end do
           
          case ('ext3')
            dp1 = z(4)-z(3)
            dp2 = z(5)-z(3)
            dm1 = z(2)-z(3)
            dm2 = z(1)-z(3)
            d2 = dm1*(dm1-dp2)/(dp1*(dp1-dp2))
            d3 = dm1*(dp1-dm1)/(dp2*(dp1-dp2))
            d5 = dm2*(dm2-dp2)/(dp1*(dp1-dp2))
            d6 = dm2*(dp1-dm2)/(dp2*(dp1-dp2))
            d1= one - d2 - d3
            d4= one - d5 - d6
            
            do i=1,nx
              dd(i,2)=dd(i,3)*d1+dd(i,4)*d2+dd(i,5)*d3  
              dd(i,1)=dd(i,3)*d4+dd(i,4)*d5+dd(i,5)*d6
            end do
          case ('ext3_zero')
            dp1 = z(4)-z(3)
            dp2 = z(5)-z(3)
            dm1 = z(2)-z(3)
            dm2 = z(1)-z(3)
            d2 = dm1*(dm1-dp2)/(dp1*(dp1-dp2))
            d3 = dm1*(dp1-dm1)/(dp2*(dp1-dp2))
            d5 = dm2*(dm2-dp2)/(dp1*(dp1-dp2))
            d6 = dm2*(dp1-dm2)/(dp2*(dp1-dp2))
            do i=1,nx
              dd(i,3) = zero
              dd(i,2) = dd(i,4)*d2+dd(i,5)*d3  
              dd(i,1) = dd(i,4)*d5+dd(i,5)*d6
            end do
          case ('ext3_zero_der')
            dp1 = z(4)-z(3)
            dp2 = z(5)-z(3)
            dm1 = z(2)-z(3)
            dm2 = z(1)-z(3)
            d2 = (dm1/dp1)**2*(dm1-dp2)/(dp1-dp2)
            d3 = (dm1/dp2)**2*(dp1-dm1)/(dp1-dp2)
            d5 = (dm2/dp1)**2*(dm2-dp2)/(dp1-dp2)
            d6 = (dm2/dp2)**2*(dp1-dm2)/(dp1-dp2)
            d1= one - d2 - d3
            d4= one - d5 - d6
            
            do i=1,nx
              dd(i,2)=dd(i,3)*d1+dd(i,4)*d2+dd(i,5)*d3  
              dd(i,1)=dd(i,3)*d4+dd(i,4)*d5+dd(i,5)*d6
            end do
          case ('ext3_zero_der2')
            dp1 = z(4)-z(3)
            dp2 = z(5)-z(3)
            dm1 = z(2)-z(3)
            dm2 = z(1)-z(3)
            d2 = dm1*(dm1*dm1 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d3 = dm1*(dp1*dp1 - dm1*dm1)/(dp2*(dp1*dp1 - dp2*dp2))
            d5 = dm2*(dm2*dm2 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d6 = dm2*(dp1*dp1 - dm2*dm2)/(dp2*(dp1*dp1 - dp2*dp2))
            d1= one - d2 - d3
            d4= one - d5 - d6
 
            do i=1,nx
              dd(i,2)=dd(i,3)*d1+dd(i,4)*d2+dd(i,5)*d3  
              dd(i,1)=dd(i,3)*d4+dd(i,4)*d5+dd(i,5)*d6
            end do
          case ('ext3_odd')
            dp1 = z(4)-z(3)
            dp2 = z(5)-z(3)
            dm1 = z(2)-z(3)
            dm2 = z(1)-z(3)
            d2 = dm1*(dm1*dm1 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d3 = dm1*(dp1*dp1 - dm1*dm1)/(dp2*(dp1*dp1 - dp2*dp2))
            d5 = dm2*(dm2*dm2 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d6 = dm2*(dp1*dp1 - dm2*dm2)/(dp2*(dp1*dp1 - dp2*dp2))
 
            do i=1,nx
              dd(i,3)=zero
              dd(i,2)=dd(i,4)*d2+dd(i,5)*d3  
              dd(i,1)=dd(i,4)*d5+dd(i,5)*d6
            end do

          case ('ext3_even')
            dp1 = z(4)-z(3)
            dp2 = z(5)-z(3)
            dm1 = z(2)-z(3)
            dm2 = z(1)-z(3)
            d2 = dm1*dm1*(dm1*dm1 - dp2*dp2)/(dp1*dp1*(dp1*dp1 - dp2*dp2))
            d3 = dm1*dm1*(dp1*dp1 - dm1*dm1)/(dp2*dp2*(dp1*dp1 - dp2*dp2))
            d5 = dm2*dm2*(dm2*dm2 - dp2*dp2)/(dp1*dp1*(dp1*dp1 - dp2*dp2))
            d6 = dm2*dm2*(dp1*dp1 - dm2*dm2)/(dp2*dp2*(dp1*dp1 - dp2*dp2))

            d1= one - d2 - d3
            d4= one - d5 - d6
 
            do i=1,nx
              dd(i,2)=dd(i,3)*d1+dd(i,4)*d2+dd(i,5)*d3  
              dd(i,1)=dd(i,3)*d4+dd(i,4)*d5+dd(i,5)*d6
            end do

          case ('ext4')
            dp1 = z(4)-z(3)
            dp2 = z(5)-z(3)
            dp3 = z(6)-z(3)
            dm1 = z(2)-z(3)
            dm2 = z(1)-z(3)
            d2 = dm1*(dp2*dp3+dm1*(dm1-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d3 = dm1*(dm1*(dp3+dp1-dm1)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d4 = dm1*(dp1*dp2+dm1*(dm1-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))
            d6 = dm2*(dp2*dp3+dm2*(dm2-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d7 = dm2*(dm2*(dp3+dp1-dm2)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d8 = dm2*(dp1*dp2+dm2*(dm2-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))

            d1= one - d2 - d3 - d4
            d5= one - d6 - d7 - d8 

            do i=1,nx
              dd(i,2)=dd(i,3)*d1+dd(i,4)*d2+dd(i,5)*d3+dd(i,6)*d4  
              dd(i,1)=dd(i,3)*d5+dd(i,4)*d6+dd(i,5)*d7+dd(i,6)*d8
            end do
 
          case ('ext4_zero')
            dp1 = z(4)-z(3)
            dp2 = z(5)-z(3)
            dp3 = z(6)-z(3)
            dm1 = z(2)-z(3)
            dm2 = z(1)-z(3)
            d2 = dm1*(dp2*dp3+dm1*(dm1-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d3 = dm1*(dm1*(dp3+dp1-dm1)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d4 = dm1*(dp1*dp2+dm1*(dm1-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))
            d6 = dm2*(dp2*dp3+dm2*(dm2-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d7 = dm2*(dm2*(dp3+dp1-dm2)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d8 = dm2*(dp1*dp2+dm2*(dm2-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))

            do i=1,nx
              dd(i,3)=zero 
              dd(i,2)=dd(i,4)*d2+dd(i,5)*d3+dd(i,6)*d4  
              dd(i,1)=dd(i,4)*d6+dd(i,5)*d7+dd(i,6)*d8
            end do
             
           
          case default
            print*,'Invalid left boundary condition in z', left
            call exit(1)
        end select
      end if

      if (iproc_z.eq.nzproc-1) then
        if (nzproc.gt.1) then 
          call z_exchange_left(dd,tag)
        end if   
        select case (right)
          case ('')
            continue
          case ('symh')
            do i=1,nx
              dd(i,nz-1)=dd(i,nz-2)
              dd(i,nz-0)=dd(i,nz-3)
            end do
          case ('symf')
            do i=1,nx
              dd(i,nz-1)=dd(i,nz-3)
              dd(i,nz-0)=dd(i,nz-4)
            end do            
          case ('asymh')
            do i=1,nx
              dd(i,nz-1)=-dd(i,nz-2)
              dd(i,nz-0)=-dd(i,nz-3)
            end do
          case ('asymf')
            do i=1,nx
              dd(i,nz-2)=zero     
              dd(i,nz-1)=-dd(i,nz-3)
              dd(i,nz-0)=-dd(i,nz-4)
            end do
          case ('ext')
            do i=1,nx
              dd(i,nz-1)=dd(i,nz-2) + (z(nz-1)-z(nz-2))  &
                                      *(dd(i,nz-2)-dd(i,nz-3)) &
                                      /(z(nz-2)-z(nz-3))
              dd(i,nz-0)=dd(i,nz-1) + (z(nz-0)-z(nz-1)) &
                                      *(dd(i,nz-1)-dd(i,nz-2))  &
                                      /(z(nz-1)-z(nz-2))
            end do
          case ('exth1')
            do i=1,nx
              dd(i,nz-1)=dd(i,nz-2) 
              dd(i,nz-0)=two*dd(i,nz-2)-dd(i,nz-3) 
            end do
          case ('extf1')
            do i=1,nx
              dd(i,nz-1)=two*dd(i,nz-2)-dd(i,nz-3)
              dd(i,nz-0)=two*dd(i,nz-2)-dd(i,nz-4)
            end do
          case ('ext2')
            do i=1,nx
              d1=ddz(dd,i,nz-4)
              d2=d2dz(dd,i,nz-4)
              d3=z(nz-1)-z(nz-4)
              dd(i,nz-1)=dd(i,nz-4)+d3*(d1+half*d2*d3)
              d3=z(nz)-z(nz-4)
              dd(i,nz-0)=dd(i,nz-4)+d3*(d1+half*d2*d3)
            end do
          case ('ext2_sym')
            dp1 = z(nz-3)-z(nz-2)
            dp2 = z(nz-4)-z(nz-2)
            dm1 = z(nz-1)-z(nz-2)
            dm2 = z(nz-0)-z(nz-2)
            d1 =  dp2**2/(dp2**2-dp1**2)
            d2 = -dp1**2/(dp2**2-dp1**2)
            do i=1,nx
              dd(i,nz-2) = d1*dd(i,nz-3) + d2*dd(i,nz-4)
              d3 = (dd(i,nz-4) - dd(i,nz-3))/(dp2**2-dp1**2)
              dd(i,nz-1) = dd(i,nz-2) + d3*dm1**2
              dd(i,nz  ) = dd(i,nz-2) + d3*dm2**2
            end do
            
          case ('ext3')
            dp1 = z(nz-3)-z(nz-2)
            dp2 = z(nz-4)-z(nz-2)
            dm1 = z(nz-1)-z(nz-2)
            dm2 = z(nz-0)-z(nz-2)
            d2 = dm1*(dm1-dp2)/(dp1*(dp1-dp2))
            d3 = dm1*(dp1-dm1)/(dp2*(dp1-dp2))
            d5 = dm2*(dm2-dp2)/(dp1*(dp1-dp2))
            d6 = dm2*(dp1-dm2)/(dp2*(dp1-dp2))

            d1= one - d2 - d3
            d4= one - d5 - d6
            
            do i=1,nx
              dd(i,nz-1)=dd(i,nz-2)*d1+dd(i,nz-3)*d2+dd(i,nz-4)*d3  
              dd(i,nz-0)=dd(i,nz-2)*d4+dd(i,nz-3)*d5+dd(i,nz-4)*d6  
            end do
          case ('ext3_zero')
            dp1 = z(nz-3)-z(nz-2)
            dp2 = z(nz-4)-z(nz-2)
            dm1 = z(nz-1)-z(nz-2)
            dm2 = z(nz-0)-z(nz-2)
            d2 = dm1*(dm1-dp2)/(dp1*(dp1-dp2))
            d3 = dm1*(dp1-dm1)/(dp2*(dp1-dp2))
            d5 = dm2*(dm2-dp2)/(dp1*(dp1-dp2))
            d6 = dm2*(dp1-dm2)/(dp2*(dp1-dp2))
            do i=1,nx
              dd(i,nz-2)=zero                         
              dd(i,nz-1)=dd(i,nz-3)*d2+dd(i,nz-4)*d3  
              dd(i,nz-0)=dd(i,nz-3)*d5+dd(i,nz-4)*d6  
            end do
          case ('ext3_zero_der')
            dp1 = z(nz-3)-z(nz-2)
            dp2 = z(nz-4)-z(nz-2)
            dm1 = z(nz-1)-z(nz-2)
            dm2 = z(nz-0)-z(nz-2)
            d2 = (dm1/dp1)**2*(dm1-dp2)/(dp1-dp2)
            d3 = (dm1/dp2)**2*(dp1-dm1)/(dp1-dp2)
            d5 = (dm2/dp1)**2*(dm2-dp2)/(dp1-dp2)
            d6 = (dm2/dp2)**2*(dp1-dm2)/(dp1-dp2)

            d1= one - d2 - d3
            d4= one - d5 - d6
            
            do i=1,nx
              dd(i,nz-1)=dd(i,nz-2)*d1+dd(i,nz-3)*d2+dd(i,nz-4)*d3  
              dd(i,nz-0)=dd(i,nz-2)*d4+dd(i,nz-3)*d5+dd(i,nz-4)*d6  
            end do
          case ('ext3_zero_der2')
            dp1 = z(nz-3)-z(nz-2)
            dp2 = z(nz-4)-z(nz-2)
            dm1 = z(nz-1)-z(nz-2)
            dm2 = z(nz-0)-z(nz-2)
                       
            d2 = dm1*(dm1*dm1 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d3 = dm1*(dp1*dp1 - dm1*dm1)/(dp2*(dp1*dp1 - dp2*dp2))
            d5 = dm2*(dm2*dm2 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d6 = dm2*(dp1*dp1 - dm2*dm2)/(dp2*(dp1*dp1 - dp2*dp2))

            d1= one - d2 - d3
            d4= one - d5 - d6
            
            do i=1,nx
              dd(i,nz-1)=dd(i,nz-2)*d1+dd(i,nz-3)*d2+dd(i,nz-4)*d3  
              dd(i,nz-0)=dd(i,nz-2)*d4+dd(i,nz-3)*d5+dd(i,nz-4)*d6  
            end do
          case ('ext3_odd')
            dp1 = z(nz-3)-z(nz-2)
            dp2 = z(nz-4)-z(nz-2)
            dm1 = z(nz-1)-z(nz-2)
            dm2 = z(nz-0)-z(nz-2)
                       
            d2 = dm1*(dm1*dm1 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d3 = dm1*(dp1*dp1 - dm1*dm1)/(dp2*(dp1*dp1 - dp2*dp2))
            d5 = dm2*(dm2*dm2 - dp2*dp2)/(dp1*(dp1*dp1 - dp2*dp2))
            d6 = dm2*(dp1*dp1 - dm2*dm2)/(dp2*(dp1*dp1 - dp2*dp2))

            do i=1,nx
              dd(i,nz-2)=zero
              dd(i,nz-1)=dd(i,nz-3)*d2+dd(i,nz-4)*d3  
              dd(i,nz-0)=dd(i,nz-3)*d5+dd(i,nz-4)*d6  
            end do

          case ('ext3_even')
            dp1 = z(nz-3)-z(nz-2)
            dp2 = z(nz-4)-z(nz-2)
            dm1 = z(nz-1)-z(nz-2)
            dm2 = z(nz-0)-z(nz-2)
            d2 = dm1*dm1*(dm1*dm1 - dp2*dp2)/(dp1*dp1*(dp1*dp1 - dp2*dp2))
            d3 = dm1*dm1*(dp1*dp1 - dm1*dm1)/(dp2*dp2*(dp1*dp1 - dp2*dp2))
            d5 = dm2*dm2*(dm2*dm2 - dp2*dp2)/(dp1*dp1*(dp1*dp1 - dp2*dp2))
            d6 = dm2*dm2*(dp1*dp1 - dm2*dm2)/(dp2*dp2*(dp1*dp1 - dp2*dp2))

            d1= one - d2 - d3
            d4= one - d5 - d6

            do i=1,nx
              dd(i,nz-1)=dd(i,nz-2)*d1+dd(i,nz-3)*d2+dd(i,nz-4)*d3  
              dd(i,nz-0)=dd(i,nz-2)*d4+dd(i,nz-3)*d5+dd(i,nz-4)*d6  
            end do

          case ('ext4')
            dp1 = z(nz-3)-z(nz-2)
            dp2 = z(nz-4)-z(nz-2)
            dp3 = z(nz-5)-z(nz-2)
            dm1 = z(nz-1)-z(nz-2)
            dm2 = z(nz-0)-z(nz-2)

            d2 = dm1*(dp2*dp3+dm1*(dm1-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d3 = dm1*(dm1*(dp3+dp1-dm1)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d4 = dm1*(dp1*dp2+dm1*(dm1-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))
            d6 = dm2*(dp2*dp3+dm2*(dm2-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d7 = dm2*(dm2*(dp3+dp1-dm2)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d8 = dm2*(dp1*dp2+dm2*(dm2-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))


            d1= one - d2 - d3 - d4
            d5= one - d6 - d7 - d8 


            do i=1,nx
              dd(i,nz-1)= dd(i,nz-2)*d1+dd(i,nz-3)*d2 &
                         +dd(i,nz-4)*d3+dd(i,nz-5)*d4  
              dd(i,nz-0)= dd(i,nz-2)*d5+dd(i,nz-3)*d6 &
                         +dd(i,nz-4)*d7+dd(i,nz-5)*d8 
            end do
          case ('ext4_zero')
            dp1 = z(nz-3)-z(nz-2)
            dp2 = z(nz-4)-z(nz-2)
            dp3 = z(nz-5)-z(nz-2)
            dm1 = z(nz-1)-z(nz-2)
            dm2 = z(nz-0)-z(nz-2)
            
            d2 = dm1*(dp2*dp3+dm1*(dm1-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d3 = dm1*(dm1*(dp3+dp1-dm1)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d4 = dm1*(dp1*dp2+dm1*(dm1-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))
            d6 = dm2*(dp2*dp3+dm2*(dm2-dp2-dp3))/(dp1*(dp1-dp2)*(dp1-dp3))
            d7 = dm2*(dm2*(dp3+dp1-dm2)-dp1*dp3)/(dp2*(dp1-dp2)*(dp2-dp3))
            d8 = dm2*(dp1*dp2+dm2*(dm2-dp1-dp2))/(dp3*(dp1-dp3)*(dp2-dp3))

            do i=1,nx
              dd(i,nz-2)= zero
              dd(i,nz-1)= dd(i,nz-3)*d2 &
                         +dd(i,nz-4)*d3+dd(i,nz-5)*d4  
              dd(i,nz-0)= dd(i,nz-3)*d6 &
                         +dd(i,nz-4)*d7+dd(i,nz-5)*d8 
            end do

          case default
            print*,'Invalid right boundary condition in z: ', right
            call exit(1)
        end select
      end if 

      if (iproc_z>0.and.iproc_z<(nzproc-1)) then
        call z_exchange(dd,tag)
      end if 
      return
      end subroutine bcz
!**********************************************


!**********************************************
      pure function ddx(fn,i,k)
      USE comm 
      REAL(RTYPE) :: ddx
      REAL(RTYPE), intent(in) :: fn(nx, nz)
      INTEGER, intent(in)     :: i,k
#if NONUNIFX==1
      ddx=fn(i-2,k)*px51(i,1)+fn(i-1,k)*px51(i,2)  &
          +fn(i,k)*px51(i,3)+fn(i+1,k)*px51(i,4)   &
          +fn(i+2,k)*px51(i,5)
#else
      ddx=fn(i-2,k)*px51u(1)+fn(i-1,k)*px51u(2) &
          +fn(i+1,k)*px51u(4) &
          +fn(i+2,k)*px51u(5)
#endif
      return
      end function ddx
!**********************************************

!*********************************************
      pure function ddz(fn,i,k)
      USE comm
      REAL(RTYPE) :: ddz
      REAL(RTYPE), intent(in) :: fn(nx,nz)
      INTEGER, intent(in)     :: i,k
#if NONUNIFZ==1
      ddz=fn(i,k-2)*pz51(k,1)+fn(i,k-1)*pz51(k,2)  &
          +fn(i,k)*pz51(k,3)+fn(i,k+1)*pz51(k,4)   &
          +fn(i,k+2)*pz51(k,5) 
#else
      ddz=fn(i,k-2)*pz51u(1)+fn(i,k-1)*pz51u(2) &
          +fn(i,k+1)*pz51u(4) &
          +fn(i,k+2)*pz51u(5)
#endif
      return
      end function ddz
!**********************************************

!*********************************************
      pure function poisson(psi,phi,i,k)

      REAL(RTYPE) :: poisson
      REAL(RTYPE), intent(in) :: psi(nx,nz),phi(nx,nz)
      INTEGER, intent(in) :: i, k

      poisson=ddx(psi,i,k)*ddz(phi,i,k)  &
             -ddx(phi,i,k)*ddz(psi,i,k)

      return
      end function poisson
!**********************************************

!**********************************************
      pure function bgrad(bx,bz,phi,i,k)

      REAL(RTYPE) :: bgrad
      REAL(RTYPE), intent(in) :: bx(nx,nz),bz(nx,nz),phi(nx,nz)
      INTEGER, intent(in) :: i, k

      bgrad=bx(i,k)*ddx(phi,i,k)  & 
           +bz(i,k)*ddz(phi,i,k)

      return
      end function bgrad
!**********************************************

!**********************************************
      pure function d2dx(fn,i,k)
      USE comm
      REAL(RTYPE) :: d2dx
      REAL(RTYPE), intent(in) :: fn(nx,nz)
      INTEGER, intent(in)     :: i,k
#if NONUNIFX==1
      d2dx=fn(i-2,k)*px52(i,1)+fn(i-1,k)*px52(i,2)  &
           +fn(i,k)*px52(i,3)+fn(i+1,k)*px52(i,4)   &
           +fn(i+2,k)*px52(i,5) 
#else
!d^2/dx^2
      d2dx=fn(i-2,k)*px52u(1)+fn(i-1,k)*px52u(2) &
           +fn(i,k)*px52u(3)+fn(i+1,k)*px52u(4) &
           +fn(i+2,k)*px52u(5)

#endif
      return
      end function d2dx
!**********************************************

!**********************************************
      pure function d2dz(fn,i,k)
      USE comm
      REAL(RTYPE) :: d2dz
      REAL(RTYPE), intent(in) :: fn(nx,nz)
      INTEGER, intent(in)     :: i,k
#if NONUNIFZ==1
      d2dz=fn(i,k-2)*pz52(k,1)+fn(i,k-1)*pz52(k,2)  &
           +fn(i,k)*pz52(k,3)+fn(i,k+1)*pz52(k,4)   &
           +fn(i,k+2)*pz52(k,5)
#else 
      d2dz=fn(i,k-2)*pz52u(1)+fn(i,k-1)*pz52u(2) &
           +fn(i,k)*pz52u(3)+fn(i,k+1)*pz52u(4) &
           +fn(i,k+2)*pz52u(5)

#endif
      return
      end function d2dz
!**********************************************

!**********************************************
      pure function del2(fn,i,k)
      USE comm
      REAL(RTYPE) :: del2
      REAL(RTYPE), intent(in) :: fn(nx,nz)
      INTEGER, intent(in)     :: i, k

      del2=d2dx(fn,i,k)+d2dz(fn,i,k)

      return
      end function del2
!**********************************************

      pure function gauss(i,k,a0,alx,alz,xc,zc)
      USE comm
      REAL(RTYPE) :: gauss
      INTEGER, intent(in) :: i, k
      REAL(RTYPE), intent(in) :: a0, alx, alz, xc, zc

      REAL(RTYPE) :: dx, dz 

      dx=x(i)-xc
      dz=z(k)-zc
      gauss=a0*exp(-alx*dx*dx-alz*dz*dz)

      return
      end function gauss
!**********************************************
!
!**********************************************
      SUBROUTINE polint(xa,ya,x,y,dy)
      IMPLICIT NONE
      REAL(RTYPE), DIMENSION(:), INTENT(IN) :: xa,ya
      REAL(RTYPE), INTENT(IN) :: x
      REAL(RTYPE), INTENT(OUT) :: y,dy
!Given arrays xa and ya of length N, and given a value x, this routine returns a value y,
!and an error estimate dy. If P(x) is the polynomial of degree N-1 such that P(xai) =yai, 
!i = 1,...,N, then the returned value y = P(x).
      INTEGER :: m,n,ns
      REAL(RTYPE), DIMENSION(size(xa)) :: c,d,den,ho
      n=assert_eq2(size(xa),size(ya),'polint')
      c=ya     ! Initialize the tableau of c's and d's.
      d=ya
      ho=xa-x
      ns=iminloc(abs(x-xa))  !  Find index ns of closest table entry.
      y=ya(ns)               !  This is the initial approximation to y.
      ns=ns-1
      do m=1,n-1                          !For each column of the tableau,
        den(1:n-m)=ho(1:n-m)-ho(1+m:n)    !we loop over the current c's and d's and up
        if (any(den(1:n-m) == 0.0)) &     ! -date them.
        call nrerror('polint: calculation failure') !This error can occur only if two input x                                                    !a's are (to within roundoff) identical.
        den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
        d(1:n-m)=ho(1+m:n)*den(1:n-m)                ! Here the c's and d's are updated.
        c(1:n-m)=ho(1:n-m)*den(1:n-m)
        if (2*ns < n-m) then         ! After each column in the tableau is completed, we decide
                                     ! which correction, c or d, we want to add to our accumulating
                                     ! value of y, i.e., which path to take through the tableau   forking
                                     !  up or down. We do this in such a way as to take the most "straight
                                     !  line" route through the tableau to its apex, updating ns
                                     !  accordingly to keep track of where we are. This route keeps the
                                     !  partial approximations centered (insofar as possible) on the
                                     !  target x. The last dy added is thus the error indication.
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        end if
        y=y+dy
      end do
      END SUBROUTINE polint
!=============================================================
      SUBROUTINE polin2(x1a,x2a,ya,x1,x2,y,dy)
      IMPLICIT NONE
      REAL(RTYPE), DIMENSION(:), INTENT(IN) :: x1a,x2a
      REAL(RTYPE), DIMENSION(:,:), INTENT(IN) :: ya
      REAL(RTYPE), INTENT(IN) :: x1,x2
      REAL(RTYPE), INTENT(OUT) :: y,dy
! Given arrays x1a of length M and x2a of length N of independent variables, and an M*N array
! of function values ya, tabulated at the grid points defined by x1a and x2a, and given
! values x1 and x2 of the independent variables, this routine returns an interpolated function
! value y, and an accuracy indication dy (based only on the interpolation in the x1 direction, however).

      INTEGER :: j,m,ndum
      REAL(RTYPE), DIMENSION(size(x1a)) :: ymtmp
      REAL(RTYPE), DIMENSION(size(x2a)) :: yntmp

      m=assert_eq2(size(x1a),size(ya,1),'polin2: m')
      ndum=assert_eq2(size(x2a),size(ya,2),'polin2: ndum')
      do j=1,m    ! Loop over rows.
        yntmp=ya(j,:)  ! Copy row into temporary storage.
        call polint(x2a,yntmp,x2,ymtmp(j),dy)  ! Interpolate answer into temporary storage 
      end do
      call polint(x1a,ymtmp,x1,y,dy)           ! Do the final interpolation.
      END SUBROUTINE polin2
!=============================================================
!**********************************************
      SUBROUTINE polint4p(xa,ya,x,y,dy)
      IMPLICIT NONE
      REAL(RTYPE), DIMENSION(4), INTENT(IN) :: xa,ya
      REAL(RTYPE), INTENT(IN) :: x
      REAL(RTYPE), INTENT(OUT) :: y,dy
      INTEGER :: m,n,ns   
      REAL(RTYPE), DIMENSION(4):: c,d,den,ho
! 4-pt polynomial interpolation.  Since the no. of pts is fixed, this is more efficient
      
      n=assert_eq2(size(xa),size(ya),'polint4p')
      c=ya     ! Initialize the tableau of c's and d's.
      d=ya
      ho=xa-x
      ns=iminloc(abs(x-xa))
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
        den(1:n-m)=ho(1:n-m)-ho(1+m:n)
        if (any(den(1:n-m) == 0.0)) &
           call nrerror('polint4p: calculation failure')       
        den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
        d(1:n-m)=ho(1+m:n)*den(1:n-m)
        c(1:n-m)=ho(1:n-m)*den(1:n-m)
        if (2*ns < n-m) then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        end if
        y=y+dy
      end do
      END SUBROUTINE polint4p
!============================================================================
      FUNCTION assert_eq2(n1,n2,string)  
! Report and die if integers not all equal (used for size checking).
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER, INTENT(IN) :: n1,n2
      INTEGER :: assert_eq2
      if (n1 == n2) then
        assert_eq2=n1
      else
        write (*,*) 'nrerror: an assert_eq failed with this tag:', string
        STOP 'program terminated by assert_eq2'
      end if
      END FUNCTION assert_eq2
!======================================================================
      FUNCTION iminloc(arr)
!Index of minloc on an array.
      REAL(RTYPE), DIMENSION(:), INTENT(IN) :: arr
      INTEGER, DIMENSION(1) :: imin
      INTEGER :: iminloc
      imin=minloc(arr(:))
      iminloc=imin(1)
      END FUNCTION iminloc
!======================================================================
      SUBROUTINE nrerror(string)
!Report a message, then die.
      CHARACTER(LEN=*), INTENT(IN) :: string
      write (*,*) 'nrerror:',string
      STOP 'program terminated by nrerror'
      END SUBROUTINE nrerror
!======================================================================
      SUBROUTINE z_exchange(ff,tag)
      include 'mpif.h'

      REAL(RTYPE) :: ff(nx,nz)
      INTEGER :: tag

      INTEGER :: id, req(2), sreq(2), ierr, ss, tag1, tag2
      integer :: status(MPI_STATUS_SIZE,2)
      REAL(RTYPE) :: bufz_in1(nx*2), bufz_in2(nx*2)
      REAL(RTYPE) :: bufz_out1(nx*2), bufz_out2(nx*2)

      tag1=tag*2
      tag2=tag1+1

      id=proc_id(iproc_x,mod(iproc_z-1+nzproc,nzproc))
      call ireceive(bufz_in1,id,req(1),tag1)   ! prepost receive from left
      id=proc_id(iproc_x,mod(iproc_z+1,nzproc))
      call ireceive(bufz_in2,id,req(2),tag2)   ! prepost receive from right
!
      bufz_out2=reshape(ff(:,nz-3:nz-2),(/nx*2/))
      bufz_out1=reshape(ff(:,3:4),(/nx*2/))
!
      id=proc_id(iproc_x,mod(iproc_z+1,nzproc))
      call isend(bufz_out2,id,sreq(1),tag1)     ! sending right
!
      id=proc_id(iproc_x,mod(iproc_z-1+nzproc,nzproc))
      call isend(bufz_out1,id,sreq(2),tag2)    ! sending left
!
      call mpi_waitall(2,req,status,ierr)
!
      ff(:,1:2)=reshape(bufz_in1,(/nx,2/))
      ff(:,nz-1:nz)=reshape(bufz_in2,(/nx,2/))
!
      call mpi_waitall(2,sreq,status,ierr)
      return

      END SUBROUTINE z_exchange
!=====================================================================
      SUBROUTINE z_exchange_left(ff,tag)
      include 'mpif.h'

      REAL(RTYPE) :: ff(nx,nz)
      INTEGER :: tag

      INTEGER :: id, req, sreq, ierr, ss, tag1, tag2
      integer :: status(MPI_STATUS_SIZE)
      REAL(RTYPE) :: bufz_in1(nx*2), bufz_in2(nx*2)
      REAL(RTYPE) :: bufz_out1(nx*2), bufz_out2(nx*2)


      tag1=tag*2
      tag2=tag1+1

      id=proc_id(iproc_x,mod(iproc_z-1+nzproc,nzproc))
      call ireceive(bufz_in1,id,req,tag1)   ! prepost receive from left
!
      bufz_out1=reshape(ff(:,3:4),(/nx*2/))
      id=proc_id(iproc_x,mod(iproc_z-1+nzproc,nzproc))
      call isend(bufz_out1,id,sreq,tag2)    ! sending left
!
      call mpi_wait(req,status,ierr)
!
      ff(:,1:2)=reshape(bufz_in1,(/nx,2/))
!
      call mpi_wait(sreq,status,ierr)

      return

      END SUBROUTINE z_exchange_left
!======================================================================
      SUBROUTINE z_exchange_right(ff,tag)
      include 'mpif.h'

      REAL(RTYPE) :: ff(nx,nz)
      INTEGER :: tag

      INTEGER :: id, req, sreq,ierr, ss, tag1, tag2
      integer :: status(MPI_STATUS_SIZE)
      REAL(RTYPE) :: bufz_in1(nx*2), bufz_in2(nx*2)
      REAL(RTYPE) :: bufz_out1(nx*2), bufz_out2(nx*2)

      tag1=tag*2
      tag2=tag1+1

      id=proc_id(iproc_x,mod(iproc_z+1,nzproc))
      call ireceive(bufz_in2,id,req,tag2)   ! prepost receive from right
!
      bufz_out2=reshape(ff(:,nz-3:nz-2),(/nx*2/))

      id=proc_id(iproc_x,mod(iproc_z+1,nzproc))
      call isend(bufz_out2,id,sreq,tag1)     ! sending right
!
      call mpi_wait(req,status,ierr)
!
      ff(:,nz-1:nz)=reshape(bufz_in2,(/nx,2/))
!
      call mpi_wait(sreq,status,ierr)
      return

      END SUBROUTINE z_exchange_right

!======================================================================
      SUBROUTINE x_exchange(ff,tag)
      include 'mpif.h'
      REAL(RTYPE) :: ff(nx,nz)
      INTEGER :: tag

      INTEGER :: id, req(2), sreq(2), ierr, ss, tag1, tag2
      integer :: status(MPI_STATUS_SIZE,2)
      REAL(RTYPE) :: bufx_in1(nz*2), bufx_in2(nz*2)
      REAL(RTYPE) :: bufx_out1(nz*2), bufx_out2(nz*2)


      tag1=tag*2
      tag2=tag1+1

      id=proc_id(mod(iproc_x-1+nxproc,nxproc),iproc_z)
      call ireceive(bufx_in1,id,req(1),tag1)   ! prepost receive from left

      id=proc_id(mod(iproc_x+1,nxproc),iproc_z)
      call ireceive(bufx_in2,id,req(2),tag2)   ! prepost receive from right
!
      bufx_out2=reshape(ff(nx-3:nx-2,:),(/nz*2/))
      bufx_out1=reshape(ff(3:4,:),(/nz*2/))
!
      id=proc_id(mod(iproc_x+1,nxproc),iproc_z)
      call isend(bufx_out2,id,sreq(1),tag1)     ! sending right
!
      id=proc_id(mod(iproc_x-1+nxproc,nxproc),iproc_z)
      call isend(bufx_out1,id,sreq(2),tag2)    ! sending left

!
      call mpi_waitall(2,req,status,ierr)
!

      ff(1:2,:)=reshape(bufx_in1,(/2,nz/))
      ff(nx-1:nx,:)=reshape(bufx_in2,(/2,nz/))
!
      call mpi_waitall(2,sreq,status,ierr)
      return

      END SUBROUTINE x_exchange
!=====================================================================

      SUBROUTINE x_exchange_left(ff,tag)
      include 'mpif.h'

      REAL(RTYPE) :: ff(nx,nz)
      INTEGER :: tag

      INTEGER :: id, req, sreq, ierr, ss, tag1, tag2
      integer :: status(MPI_STATUS_SIZE)
      REAL(RTYPE) :: bufx_in1(nz*2), bufx_in2(nz*2)
      REAL(RTYPE) :: bufx_out1(nz*2), bufx_out2(nz*2)


      tag1=tag*2
      tag2=tag1+1

      id=proc_id(mod(iproc_x-1+nxproc,nxproc),iproc_z)
      call ireceive(bufx_in1,id,req,tag1)   ! prepost receive from left
!     
      bufx_out1=reshape(ff(3:4,:),(/nz*2/))

      id=proc_id(mod(iproc_x-1+nxproc,nxproc),iproc_z)
      call isend(bufx_out1,id,sreq,tag2)    ! sending left
!
      call mpi_wait(req,status,ierr)
!
      ff(1:2,:)=reshape(bufx_in1,(/2,nz/))
!
      call mpi_wait(sreq,status,ierr)
      return

      END SUBROUTINE x_exchange_left

!=====================================================================

      SUBROUTINE x_exchange_right(ff,tag)
      include 'mpif.h'

      REAL(RTYPE) :: ff(nx,nz)
      INTEGER :: tag

      INTEGER :: id, req, sreq, ierr, ss, tag1, tag2
      integer :: status(MPI_STATUS_SIZE)
      REAL(RTYPE) :: bufx_in1(nz*2), bufx_in2(nz*2)
      REAL(RTYPE) :: bufx_out1(nz*2), bufx_out2(nz*2)

      tag1=tag*2
      tag2=tag1+1

      id=proc_id(mod(iproc_x+1,nxproc),iproc_z)
      call ireceive(bufx_in2,id,req,tag2)   ! prepost receive from right
!
      bufx_out2=reshape(ff(nx-3:nx-2,:),(/nz*2/))

      id=proc_id(mod(iproc_x+1,nxproc),iproc_z)
      call isend(bufx_out2,id,sreq,tag1)     ! sending right
!
      call mpi_wait(req,status,ierr)  
!
      ff(nx-1:nx,:)=reshape(bufx_in2,(/2,nz/))
!
      call mpi_wait(sreq,status,ierr)
      return
      
      END SUBROUTINE x_exchange_right

!=====================================================================
      function integrate(ff)
      USE comm
      USE mp     

      REAL(RTYPE) :: integrate
      REAL(RTYPE) :: ff(nx,nz)
      REAL(RTYPE) :: ss, s1
      REAL(RTYPE) :: weightx(nx), weightz(nz)
      INTEGER :: i, k
! weights along x
      do i=3,nx-2
        weightx(i)=(x0(i+1)-x0(i))*dxdx0(i)
      end do
! weights along z
      do k=3,nz-2
        weightz(k)=(z0(k+1)-z0(k))*dzdz0(k)
      end do
! Now integrate
      ss=zero   
      !$OMP PARALLEL DO REDUCTION(+:ss)
      do k=3,nz-2
        s1=zero
        do i=3,nx-2
          s1=s1+ff(i,k)*weightx(i)         
        end do
        ss=ss+s1*weightz(k)
      end do
      !$OMP END PARALLEL DO
!
      if (nproc.gt.1) then
        call sum_allreduce(ss)
      end if 

      integrate=ss
      return
      end function integrate

!=====================================================================


      function proc_id(ip,kp)
      INTEGER :: ip, kp 
      INTEGER :: ib, kb, ir, kr
      INTEGER :: proc_id

! return the id of PE for the proc (ip,kp)

      ib=ip/block_dim_x
      kb=kp/block_dim_z
      ir=ip-ib*block_dim_x
      kr=kp-kb*block_dim_z

      proc_id=(ib+nxblock*kb)*blocksize+(ir+kr*block_dim_x)

      end function proc_id 
!=====================================================================
      END MODULE ssub
