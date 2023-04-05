#include"inc.h"
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

#if TIMESTEPPING == 0
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
#endif

#if TIMESTEPPING == 1
! SSP RK2
      if (s.eq.0)then
        !$OMP DO 
        do k=3,nz-2
        do i=3,nx-2
          wk=fi(i,k)+df(i,k)*dt
          f(i,k)=fi(i,k)
          fi(i,k)=wk 
        end do
        end do
        !$OMP END DO
      else
        !$OMP DO 
        do k=3,nz-2
        do i=3,nx-2
          fi(i,k)=half*(f(i,k)+fi(i,k)+df(i,k)*dt)
        end do
        end do
        !$OMP END DO
      end if
#endif


#if TIMESTEPPING == 2
! SSP RK3
      if (s.eq.0) then
        !$OMP DO 
        do k=3,nz-2
        do i=3,nx-2
          wk=fi(i,k)+df(i,k)*dt
          f(i,k)=fi(i,k)
          fi(i,k)=wk 
        end do
        end do
        !$OMP END DO
      else if (s.eq.1) then
        !$OMP DO 
        do k=3,nz-2
        do i=3,nx-2
          fi(i,k)=(3._RTYPE*f(i,k)+fi(i,k)+df(i,k)*dt)/4._RTYPE
        end do
        end do
        !$OMP END DO
      else 
        !$OMP DO 
        do k=3,nz-2
        do i=3,nx-2
          fi(i,k)=(f(i,k)+two*(fi(i,k)+df(i,k)*dt))/3._RTYPE
        end do
        end do
        !$OMP END DO
      end if
#endif

#if TIMESTEPPING == 3
! 4-stage SSP RK3
      if (s.eq.0) then
        !$OMP DO 
        do k=3,nz-2
        do i=3,nx-2
          wk=fi(i,k)+half*df(i,k)*dt
          f(i,k)=fi(i,k)
          fi(i,k)=wk
        end do
        end do
        !$OMP END DO
      else if (s.eq.1) then
        !$OMP DO 
        do k=3,nz-2
        do i=3,nx-2
          fi(i,k)=fi(i,k)+half*df(i,k)*dt
        end do 
        end do
        !$OMP END DO
      else if (s.eq.2) then
        !$OMP DO 
        do k=3,nz-2
        do i=3,nx-2
          fi(i,k)=(two*f(i,k)+fi(i,k)+half*df(i,k)*dt)/3._RTYPE
        end do
        end do
        !$OMP END DO
      else 
        !$OMP DO 
        do k=3,nz-2
        do i=3,nx-2
          fi(i,k)=fi(i,k)+half*df(i,k)*dt
        end do
        end do
        !$OMP END DO
      end if
#endif

      end subroutine step2d
!=====================================================================
      subroutine step1d_x(fi,f,df,s)
! time stepping for 1D array along x
      IMPLICIT NONE
      REAL(RTYPE) :: fi(nx), f(nx), df(nx)
      INTEGER :: s     ! stage
      REAL(RTYPE) :: wk
      INTEGER :: i

#if TIMESTEPPING == 0
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
#endif

#if TIMESTEPPING == 1
! SSP RK2
      if (s.eq.0)then
        !$OMP DO
        do i=3,nx-2
          wk=fi(i)+df(i)*dt
          f(i)=fi(i)
          fi(i)=wk 
        end do
        !$OMP END DO
      else
        !$OMP DO
        do i=3,nx-2
          fi(i)=half*(f(i)+fi(i)+df(i)*dt)
        end do
        !$OMP END DO
      end if
#endif


#if TIMESTEPPING == 2
! SSP RK3
      if (s.eq.0) then
        !$OMP DO
        do i=3,nx-2
          wk=fi(i)+df(i)*dt
          f(i)=fi(i)
          fi(i)=wk 
        end do
        !$OMP END DO
      else if (s.eq.1) then
        !$OMP DO
        do i=3,nx-2
          fi(i)=(3._RTYPE*f(i)+fi(i)+df(i)*dt)/4._RTYPE
        end do
        !$OMP END DO
      else 
        !$OMP DO
        do i=3,nx-2
          fi(i)=(f(i)+two*(fi(i)+df(i)*dt))/3._RTYPE
        end do
        !$OMP END DO
      end if
#endif

#if TIMESTEPPING == 3
! 4-stage SSP RK3
      if (s.eq.0) then
        !$OMP DO
        do i=3,nx-2
          wk=fi(i)+half*df(i)*dt
          f(i)=fi(i)
          fi(i)=wk
        end do
        !$OMP END DO
      else if (s.eq.1) then
        !$OMP DO
        do i=3,nx-2
          fi(i)=fi(i)+half*df(i)*dt
        end do
        !$OMP END DO 
      else if (s.eq.2) then
        !$OMP DO
        do i=3,nx-2
          fi(i)=(two*f(i)+fi(i)+half*df(i)*dt)/3._RTYPE
        end do
        !$OMP END DO
      else 
        !$OMP DO
        do i=3,nx-2
          fi(i)=fi(i)+half*df(i)*dt
        end do
        !$OMP END DO
      end if
#endif

      end subroutine step1d_x
!=====================================================================
      subroutine step1d_z(fi,f,df,s)
! time stepping for 1D array along z
      IMPLICIT NONE
      REAL(RTYPE) :: fi(nz), f(nz), df(nz)
      INTEGER :: s     ! stage
      REAL(RTYPE) :: wk
      INTEGER :: k

#if TIMESTEPPING == 0
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
#endif

#if TIMESTEPPING == 1
! SSP RK2
      if (s.eq.0)then
        !$OMP DO
        do k=3,nz-2
          wk=fi(k)+df(k)*dt
          f(k)=fi(k)
          fi(k)=wk 
        end do
        !$OMP END DO
      else
        !$OMP DO
        do k=3,nz-2
          fi(k)=half*(f(k)+fi(k)+df(k)*dt)
        end do
        !$OMP END DO
      end if
#endif


#if TIMESTEPPING == 2
! SSP RK3
      if (s.eq.0) then
        !$OMP DO
        do k=3,nz-2
          wk=fi(k)+df(k)*dt
          f(k)=fi(k)
          fi(k)=wk 
        end do
        !$OMP END DO
      else if (s.eq.1) then
        !$OMP DO
        do k=3,nz-2
          fi(k)=(3._RTYPE*f(k)+fi(k)+df(k)*dt)/4._RTYPE
        end do
        !$OMP END DO
      else 
        !$OMP DO
        do k=3,nz-2
          fi(k)=(f(k)+two*(fi(k)+df(k)*dt))/3._RTYPE
        end do
        !$OMP END DO
      end if
#endif

#if TIMESTEPPING == 3
! 4-stage SSP RK3
      if (s.eq.0) then
        !$OMP DO
        do k=3,nz-2
          wk=fi(k)+half*df(k)*dt
          f(k)=fi(k)
          fi(k)=wk
        end do
        !$OMP END DO
      else if (s.eq.1) then   
        !$OMP DO
        do k=3,nz-2
          fi(k)=fi(k)+half*df(k)*dt
        end do
        !$OMP END DO 
      else if (s.eq.2) then
        !$OMP DO
        do k=3,nz-2
          fi(k)=(two*f(k)+fi(k)+half*df(k)*dt)/3._RTYPE
        end do
        !$OMP END DO
      else 
        !$OMP DO
        do k=3,nz-2
          fi(k)=fi(k)+half*df(k)*dt
        end do
        !$OMP END DO
      end if
#endif

      end subroutine step1d_z
!=========================================================


      end module stepping
