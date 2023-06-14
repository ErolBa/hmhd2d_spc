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

#include"inc.h"


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

! #  if FDTYPE==0
! define(`DDZ', `($1(i,k-2)*pz51(k,1)+$1(i,k-1)*pz51(k,2) &
!                   +$1(i,k  )*pz51(k,3)+$1(i,k+1)*pz51(k,4) &
!                   +$1(i,k+2)*pz51(k,5))')
! #  else
! define(`DDZ', `($1(i,k-2)*pz51(k,1)+$1(i,k-1)*pz51(k,2) &
!                   +$1(i,k+1)*pz51(k,4) &
!                   +$1(i,k+2)*pz51(k,5))')
! #  endif

! uniform grid
dx=(xmax-xmin)/(nxtot-4)
do i=1,nxtot
  x0tot(i)=dx*(i-3)+xmin
end do
x0=x0tot(ishift_x+1:ishift_x+nx)

#if NONUNIFX==1
! nonuniform grid. Remap with 5-th order polynormial determined by 
! mesh size ratio

c1=15._RTYPE/(7._RTYPE*xmr+8._RTYPE)
c2=(one-c1)*40._RTYPE/7._RTYPE
  
do i=1,nxtot
  xtot(i)=c1*x0tot(i)+c2*x0tot(i)**3*(1-1.2_RTYPE*x0tot(i)**2)
end do

xtot(1   )=two*xmin-xtot(5)
xtot(2   )=two*xmin-xtot(4)
xtot(nxtot  )=two*xmax-xtot(nxtot-4)
xtot(nxtot-1)=two*xmax-xtot(nxtot-3)
 
x=xtot(ishift_x+1:ishift_x+nx)
#else
x=x0
xtot=x0tot
#endif

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

#if NONUNIFZ==1
! nonuniform grid. Remap with 5-th order polynormial determined by
! mesh size ratio

c1=15._RTYPE/(7._RTYPE*zmr+8._RTYPE)
c2=(one-c1)*40._RTYPE/7._RTYPE

do k=1,nztot
  ztot(k)=c1*z0tot(k)+c2*z0tot(k)**3*(1-1.2_RTYPE*z0tot(k)**2)
end do

ztot(1   )=two*zmin-ztot(5)
ztot(2   )=two*zmin-ztot(4)
ztot(nztot  )=two*zmax-ztot(nztot-4)
ztot(nztot-1)=two*zmax-ztot(nztot-3)
 
z=ztot(ishift_z+1:ishift_z+nz)
#else
z=z0
ztot=z0tot
#endif

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

#if TWO_AND_HALF==1 
!$OMP SECTION
call bcxper(byi,5)
call bcz(byi,'ext3','ext3',5)
!call bcz(byi,'ext3_zero_der','ext3_zero_der',5)
!$OMP SECTION
call bcxper(puyi,6)
call bcz(puyi,'ext3_zero_der','ext3_zero_der',6)
#endif
#if ISOTHERMAL==0
!$OMP SECTION
call bcxper(prei,7)
call bcz(prei,'ext3','ext3',7)
#  if (SPE==1 && TWOFLUID==1)
!$OMP SECTION
call bcxper(pei,8)
call bcz(pei,'ext3','ext3',8)
#  endif 
#endif

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
#if (TWOFLUID==1 || HYPER_RESISTIVITY==1)
! All ghost values of j's are needed to determine ve
! or when hyper-resistivity is used
!$OMP SECTION
call bcxper(jyi,1)
call bcz(jyi,'ext3','ext3',1)
#  if (TWO_AND_HALF==1)
!$OMP SECTION
call bcxper(jxi,2)
call bcz(jxi,'ext3','ext3',2)
!$OMP SECTION
call bcxper(jzi,3)
call bcz(jzi,'ext3','ext3',3) 
#  endif      
#endif
!$OMP END SECTIONS
return
end subroutine fin_aux

!########################################################################
subroutine fin_b
USE comm
IMPLICIT NONE
#if (TWO_AND_HALF==1 || JXB==0 || JY_CAL==1 || VDEPT_HDIF>2)
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
#endif      

!
return
end subroutine fin_b
!########################################################################
subroutine fin_e
USE comm
IMPLICIT NONE
INTEGER :: i, k, ip
!$OMP DO
! Add Ey0

! calc misc field, which is essentially E_bs
do k=3,nz-2
      do i=3, nx-2
            misc(i,k) =  bs_curr_const * ddz(prei, i, k) ! no eta, so misc includes jbs
      end do
end do



! if(iproc.eq.0) then
! ! If there is need to measure different components of Ey
!       do k=3,nz-2
!             do i=3, nx-2
!                   misc(i,k) = eta * bs_curr_const * ddz(prei, i, k)
!             end do
!       end do

!       ! write(*,*) "Max init (", iproc, ") ", MAXVAL(ABS(ey)), " ",NEW_LINE('a'), "Ey0", MAXVAL(ABS(Ey0)),  NEW_LINE('a'), "Jbs", max
      
!       write(*,*) "Ey0 ", NEW_LINE('a'),Ey0
!       ! write(*,*) "eta*J_bs ", NEW_LINE('a'),misc

! endif

do k=3,nz-2
  do i=3, nx-2
    ey(i,k) = ey(i,k)+Ey0(k) 
    ey(i,k) = ey(i,k) + eta * bs_curr_const * ddz(prei, i, k)! add bootstrap current term here(Erol) ! using ddx(p,i,k)
!     bs_curr_const should be on the order ~ 1e-3, 1e-4
!     if too large, HMHD crashes
  end do
end do
!$OMP END DO
#if TWO_AND_HALF==1 
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
#endif
return
end subroutine fin_e
!########################################################################

subroutine initialize

INTEGER :: i, k
REAL(RTYPE) :: a

a = 3.0_RTYPE*sqrt(3.0_RTYPE)/4.0_RTYPE
do k=1,nz
  do i=1,nx
psii(i,k)=1.00*-0.00744195*a + a/cosh(z(k))**2
    ! add perturbation
    psii(i,k)=psii(i,k)-(xl/pi2)*qpb*cos(pi2*x(i)/xl)*cos(pi*z(k)/zl)
    psi0(i,k)= zero
deni(i,k)=1.00*2.0
prei(i,k)=1.00*3.0*(z(k)+pi)/(2*pi) + 0.1
    pei(i,k)= prei(i,k)
!     deni(i,k) = 1.0
    ! perturbed flow
    puxi(i,k)= qpc*(xl/zl)*sin(two*pi*x(i)/xl)*cos(2*pi*z(k)/zl)
    puyi(i,k)= zero
    puzi(i,k)= -qpc*cos(two*pi*x(i)/xl)*sin(2*pi*z(k)/zl)
! byi(i,k)=1.00*sqrt(qpa**2 - 4*a**2*sinh(z(k))**2/cosh(z(k))**6 )!- 4.0 * prei(i,k) )
  end do
end do

do k=1,nz
      do i=1,nx
            byi(i,k)=1.00*sqrt(qpa**2 - (ddz(psii, i, k))**2 - (ddx(psii, i, k))**2 )!- 4.0 * prei(i,k) )
      end do
end do

call restart

! Calculate beta
! if(iproc.eq.0) then
!       write(*,*) "BETAS (",iproc,') '
! write(*,*) (2*MAXVAL(prei)) / qpa**2,","

! betatotal SUM()

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
psiii(i,k)=1.00*-0.00744195*a + a/cosh(z(k))**2
byii(i,k)=1.00*sqrt(qpa**2 - 4*a**2*sinh(z(k))**2/cosh(z(k))**6) !- 4.0 * prei(i,k) )

  end do
end do


do k=3, nz-2
  Ey0(k) = d2dz(psiii,3,k)*eta ! Important: use psiii, not psii (which also has a perturbation)
!   Ex0(k) = ddz(byi,3,k)*eta 
  Ex0(k) = 0.0 ! set to zero to simplify - Erol
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