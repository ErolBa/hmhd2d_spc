












! ================================================
! This module defines the problem to solve.
! Initial conditions and boundary conditions go here.
! Also the mesh is defined here.
! ================================================

! Setup of island saturation problem of Loureiro et al. PRL 2005
! qpa is the guide field
! qpb is the amplitude of the initial perturbation in Bz
! qpc is the amplitude of a perturbed flow

! Applied external E field to prevent the equilibrium from decaying.

! The standard is to set zl=2*pi. To be unstable needs xl>sqrt(5)*zl.
! Delta' = 2(5-k^2)(3+k^2)/k^2/sqrt(4+k^2)
! Here, k=2*pi/xl
! Note that this analystic formula assumes an infinite domain,
! where as the simulation domain is from -xl/2 to xl/2

! It is recommended to perturb flow instead of the magnetic flux,
! so that there is no seed island.

! xl       Delta'
! 2.9427   0.5
! 3.0      0.7135
! 3.0773   1.0
! 3.3491   2.0
! 3.89     4
! 4.0      4.413
! 4.9527   8.2
! 5.0      8.399
! 6.0      12.918
! 6.8577   17.3
! 7        18.075

! Important!! This version uses conducting Hard wall BC, which is
! different from Loureiro's periodic BC.




































































































































































































      MODULE prob
      USE globals
      USE comm
      USE ssub
      IMPLICIT NONE

      INTEGER, PARAMETER :: nhistmax=2002,  &
                            nhq=1            ! number of stored history quanitie
      REAL(RTYPE) ::  thist(nhistmax,nhq)

      REAL(RTYPE) :: Ey0(nz), Ex0(nz)  ! external E field

      CONTAINS
!#######################################################################
      subroutine init_prob
! problem dependent resources go here
      return
      end subroutine init_prob
!#######################################################################
      subroutine finish_prob
! free problem dependent resources
      return
      end subroutine finish_prob
!#######################################################################
      subroutine meshx

      INTEGER :: i
      REAL(RTYPE) :: dx, c1, c2, xmin, xmax

      xmin=-half
      xmax=half

! uniform grid
      dx=(xmax-xmin)/(nxtot-4)
      do i=1,nxtot
        x0tot(i)=dx*(i-3)+xmin
      end do
      x0=x0tot(ishift_x+1:ishift_x+nx)

      x=x0
      xtot=x0tot

      if (xl.ne.one) then
        x=xl*x     ! rescale
        xtot=xl*xtot
        x0=xl*x0
        x0tot=xl*x0tot
      end if

      end subroutine meshx
!#######################################################################
      subroutine meshz
      INTEGER :: k
      REAL(RTYPE) :: dz, c1, c2, zmin, zmax

      zmin=-half
      zmax=half
! uniform grid
      dz=(zmax-zmin)/(nztot-5)
      do k=1,nztot
        z0tot(k)=dz*(k-3)+zmin
      end do
      z0=z0tot(ishift_z+1:ishift_z+nz)

      z=z0
      ztot=z0tot

      if (zl.ne.one) then
        z=zl*z     ! rescale
        ztot=zl*ztot
        z0=zl*z0
        z0tot=zl*z0tot
      end if

      end subroutine meshz

!#######################################################################

      subroutine fin_primary
      USE comm
      IMPLICIT NONE
      INTEGER :: i, k

! Boundary conditions
      !$OMP SECTIONS
      !$OMP SECTION
      call bcxper(deni,1)
      call bcz(deni,'ext3','ext3',1)
      !$OMP SECTION
      call bcxper(puxi,2)
      call bcz(puxi,'ext3_zero_der','ext3_zero_der',2)
      !$OMP SECTION
      call bcxper(puzi,3)
      call bcz(puzi,'ext3_zero','ext3_zero',3)
      !$OMP SECTION
      call bcxper(psii,4)
      call bcz(psii,'ext3_zero','ext3_zero',4)
      !call bcz(psii,'ext3_odd','ext3_odd',4)

      !$OMP SECTION
      call bcxper(byi,5)
      call bcz(byi,'ext3','ext3',5)
      !call bcz(byi,'ext3_zero_der','ext3_zero_der',5)
      !$OMP SECTION
      call bcxper(puyi,6)
      call bcz(puyi,'ext3_zero_der','ext3_zero_der',6)
      !$OMP SECTION
      call bcxper(prei,7)
      call bcz(prei,'ext3','ext3',7)

!
      !$OMP END SECTIONS
      return
      end subroutine fin_primary

!########################################################################
      subroutine fin_aux
      USE comm
      IMPLICIT NONE
      INTEGER :: i, k

      !$OMP SECTIONS
      !$OMP END SECTIONS
      return
      end subroutine fin_aux

!########################################################################
      subroutine fin_b
      USE comm
      IMPLICIT NONE
! Have to fill in the ghost cells of bxi and bzi when step spuy
! or when maxwell stress tensor is used.
      !$OMP SECTIONS
      !$OMP SECTION
      call bcxper(bxi,5)
      call bcz(bxi,'ext3','ext3',5)
      !call bcz(bxi,'ext3_zero_der','ext3_zero_der',5)
      !$OMP SECTION
      call bcxper(bzi,6)
      !call bcz(bzi,'ext3_odd','ext3_odd',6)
      call bcz(bzi,'ext3','ext3',6)
      !$OMP END SECTIONS

!
      return
      end subroutine fin_b
!########################################################################
      subroutine fin_e
      USE comm
      IMPLICIT NONE
      INTEGER :: i, k
      !$OMP DO
      ! Add Ey0
      do k=3,nz-2
        do i=3, nx-2
          ey(i,k)=ey(i,k)+Ey0(k)
        end do
      end do
      !$OMP END DO
      !$OMP DO
      ! Add Ex0
      do k=3,nz-2
        do i=3, nx-2
          ex(i,k)=ex(i,k)+Ex0(k)
        end do
      end do
      !$OMP END DO
      !$OMP SECTIONS
      !$OMP SECTION
      call bcz(ex,'ext3_zero','ext3_zero',1)
      !$OMP SECTION
      call bcxper(ez,2)
      !$OMP END SECTIONS
      return
      end subroutine fin_e
!########################################################################

      subroutine initialize

      INTEGER :: i, k
      REAL(RTYPE) :: a

      a = 3.0_RTYPE*sqrt(3.0_RTYPE)/4.0_RTYPE

      do k=1,nz
        do i=1,nx
		psii(i,k)=1.00*a*(1.2*z(k)**3 + 1)/cosh(z(k))**2 + 0.269455*a
          ! add perturbation
          psii(i,k)=psii(i,k)-(xl/pi2)*qpb*cos(pi2*x(i)/xl)*cos(pi*z(k)/zl)
          psi0(i,k)= zero
          deni(i,k)= one
          prei(i,k)= T0*deni(i,k)
          pei(i,k)= prei(i,k)
          ! perturbed flow
          puxi(i,k)= qpc*(xl/zl)*sin(two*pi*x(i)/xl)*cos(2*pi*z(k)/zl)
          puyi(i,k)= zero
          puzi(i,k)= -qpc*cos(two*pi*x(i)/xl)*sin(2*pi*z(k)/zl)
		byi(i,k)=1.00*3.6*sqrt(0.0771604938271605*qpa**2 - (a*z(k)**2/cosh(z(k))**2 - 0.555555555555556*a*(1.2*z(k)**3 + 1)*sinh(z(k))/cosh(z(k))**3)**2)
        end do
      end do

      call restart

      end subroutine initialize
!########################################################################
      subroutine restart
! special requirement when restart goes here
      INTEGER :: i, k
      REAL(RTYPE) :: psiii(nx, nz), byii(nx, nz)
      REAL(RTYPE) :: a

      ! calculate Ex0, Ey0 to balance the resistive decay
      
      

      a = 3.0_RTYPE*sqrt(3.0_RTYPE)/4.0_RTYPE
      do k=1,nz
        do i=1,nx
      !     psiii(i,k)=a*(one/cosh(z(k))**2 - one/cosh(zl/2.0_RTYPE)**2)
      !     byii(i,k)= sqrt(qpa**2 - (a*two*tanh(z(k))/cosh(z(k))**2)**2)
		psiii(i,k)=1.00*a*(1.2*z(k)**3 + 1)/cosh(z(k))**2 + 0.269455*a
		byii(i,k)=1.00*3.6*sqrt(0.0771604938271605*qpa**2 - (a*z(k)**2/cosh(z(k))**2 - 0.555555555555556*a*(1.2*z(k)**3 + 1)*sinh(z(k))/cosh(z(k))**3)**2)
        end do
      end do

      ! write(*,*) "i=1,k,psii(1,k)"
      ! do k=1,nz
      !       write(*,*) k,psiii(1,k)
      ! end do
      ! stop


      do k=3, nz-2
        Ey0(k) = d2dz(psiii,3,k)*eta
        Ex0(k) = ddz(byii,3,k)*eta
      end do
      return
      end subroutine restart
!########################################################################
      subroutine diagnostics
! problem specific diagnostics go here
      return
      end subroutine diagnostics
!########################################################################
      subroutine tdiagcol
!
!-----------------------------------------------------------------------
!
! ****** Collect time history diagnostics.
!
!-----------------------------------------------------------------------
!
      thist(ihist,1)=time
      return
      end subroutine tdiagcol
!#######################################################################

      END MODULE prob
