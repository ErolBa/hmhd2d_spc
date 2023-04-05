














































 






 










 


















 








 























 

  










 
 






 
























 

 

 









 


























      MODULE stepping 
      USE globals, ONLY : nx,nz,RTYPE
      USE comm, ONLY : dt, half, two
      IMPLICIT NONE
      CONTAINS


!=====================================================================
      subroutine step2d(fi,f,df,s)
! time stepping for 2D array
      IMPLICIT NONE
      REAL(RTYPE) :: fi(nx,nz), f(nx,nz), df(nx,nz)
      INTEGER :: s     ! stage
      REAL(RTYPE) :: wk
      INTEGER :: i,k

! Trapezoidal Leapfrog
      if (s.eq.0)then
        !$OMP DO 
        do k=3,nz-2
        do i=3,nx-2
          wk=half*(f(i,k)+fi(i,k))+df(i,k)*dt
          f(i,k)=fi(i,k)
          fi(i,k)=wk 
        end do
        end do
        !$OMP END DO
      else
        !$OMP DO 
        do k=3,nz-2
        do i=3,nx-2
          fi(i,k)=f(i,k)+df(i,k)*dt
        end do
        end do
        !$OMP END DO
      end if





      end subroutine step2d
!=====================================================================
      subroutine step1d_x(fi,f,df,s)
! time stepping for 1D array along x
      IMPLICIT NONE
      REAL(RTYPE) :: fi(nx), f(nx), df(nx)
      INTEGER :: s     ! stage
      REAL(RTYPE) :: wk
      INTEGER :: i

! Trapezoidal Leapfrog
      if (s.eq.0)then
        !$OMP DO
        do i=3,nx-2
          wk=half*(f(i)+fi(i))+df(i)*dt
          f(i)=fi(i)
          fi(i)=wk 
        end do
        !$OMP END DO
      else
        !$OMP DO
        do i=3,nx-2
          fi(i)=f(i)+df(i)*dt
        end do
        !$OMP END DO
      end if





      end subroutine step1d_x
!=====================================================================
      subroutine step1d_z(fi,f,df,s)
! time stepping for 1D array along z
      IMPLICIT NONE
      REAL(RTYPE) :: fi(nz), f(nz), df(nz)
      INTEGER :: s     ! stage
      REAL(RTYPE) :: wk
      INTEGER :: k

! Trapezoidal Leapfrog
      if (s.eq.0)then
        !$OMP DO
        do k=3,nz-2
          wk=half*(f(k)+fi(k))+df(k)*dt
          f(k)=fi(k)
          fi(k)=wk 
        end do
        !$OMP END DO
      else
        !$OMP DO
        do k=3,nz-2
          fi(k)=f(k)+df(k)*dt
        end do
        !$OMP END DO
      end if





      end subroutine step1d_z
!=========================================================


      end module stepping
