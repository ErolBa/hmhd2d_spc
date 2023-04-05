#include"inc.h"

!############################################################
! Definitions for the m4 preprocessor
!############################################################

! Undefine intrinsic macros that conflict with f90

undefine(`include')   ! include statements
undefine(`len')       ! len defines length of character variables
undefine(`index')     ! dummy argument in sub set_fname
undefine(`format')

#if NONUNIFX==1
! d/dx
#  if FDTYPE==0
define(`DDX', `($1(i-2,k)*px51(i,1)+$1(i-1,k)*px51(i,2) &
                +$1(i,k)*px51(i,3)+$1(i+1,k)*px51(i,4) &
                +$1(i+2,k)*px51(i,5))')
#  else
define(`DDX', `($1(i-2,k)*px51(i,1)+$1(i-1,k)*px51(i,2) &
                +$1(i+1,k)*px51(i,4) &
                +$1(i+2,k)*px51(i,5))')
#  endif 
!d^2/dx^2
define(`D2DX', `($1(i-2,k)*px52(i,1)+$1(i-1,k)*px52(i,2) &
                 +$1(i  ,k)*px52(i,3)+$1(i+1,k)*px52(i,4) &
                 +$1(i+2,k)*px52(i,5))')
#else
! d/dx
define(`DDX', `($1(i-2,k)*px51u(1)+$1(i-1,k)*px51u(2) &
                +$1(i+1,k)*px51u(4) &
                +$1(i+2,k)*px51u(5))')
!d^2/dx^2
define(`D2DX', `($1(i-2,k)*px52u(1)+$1(i-1,k)*px52u(2) &
                 +$1(i  ,k)*px52u(3)+$1(i+1,k)*px52u(4) &
                 +$1(i+2,k)*px52u(5))')

#endif

#if NONUNIFZ==1
! d/dz
#  if FDTYPE==0
define(`DDZ', `($1(i,k-2)*pz51(k,1)+$1(i,k-1)*pz51(k,2) &
                +$1(i,k  )*pz51(k,3)+$1(i,k+1)*pz51(k,4) &
                +$1(i,k+2)*pz51(k,5))')
#  else
define(`DDZ', `($1(i,k-2)*pz51(k,1)+$1(i,k-1)*pz51(k,2) &
                +$1(i,k+1)*pz51(k,4) &
                +$1(i,k+2)*pz51(k,5))')
#  endif
! d^2/dz^2
define(`D2DZ', `($1(i,k-2)*pz52(k,1)+$1(i,k-1)*pz52(k,2) &
                 +$1(i,k  )*pz52(k,3)+$1(i,k+1)*pz52(k,4) &
                 +$1(i,k+2)*pz52(k,5))')
#else
! d/dz
define(`DDZ', `($1(i,k-2)*pz51u(1)+$1(i,k-1)*pz51u(2) &
                +$1(i,k+1)*pz51u(4) &
                +$1(i,k+2)*pz51u(5))')
! d^2/dz^2
define(`D2DZ', `($1(i,k-2)*pz52u(1)+$1(i,k-1)*pz52u(2) &
                 +$1(i,k  )*pz52u(3)+$1(i,k+1)*pz52u(4) &
                 +$1(i,k+2)*pz52u(5))')
#endif

! 3-point stencil finite difference

#if NONUNIFX==1
! d/dx
define(`DDX3', `($1(i-1,k)*px31(i,1) &
                +$1(i  ,k)*px31(i,2) &
                +$1(i+1,k)*px31(i,3))')
! d^2/dx^2
define(`D2DX3', `($1(i-1,k)*px32(i,1) &
                 +$1(i  ,k)*px32(i,2) &
                 +$1(i+1,k)*px32(i,3))')
#else
! d/dx
define(`DDX3', `($1(i-1,k)*px31u(1) &
                +$1(i+1,k)*px31u(3))')
! d^2/dx^2
define(`D2DX3', `($1(i-1,k)*px32u(1) &
                 +$1(i  ,k)*px32u(2) &
                 +$1(i+1,k)*px32u(3))')
#endif



#if NONUNIFZ==1
! d/dz
define(`DDZ3', `($1(i,k-1)*pz31(k,1) &
                +$1(i,k  )*pz31(k,2) &
                +$1(i,k+1)*pz31(k,3))')
! d^2/dz^2
define(`D2DZ3', `($1(i,k-1)*pz32(k,1) &
                 +$1(i,k  )*pz32(k,2) &
                 +$1(i,k+1)*pz32(k,3))')
#else
! d/dz
define(`DDZ3', `($1(i,k-1)*pz31u(1) &
                +$1(i,k+1)*pz31u(3))')
! d^2/dz^2
define(`D2DZ3', `($1(i,k-1)*pz32u(1) &
                 +$1(i,k  )*pz32u(2) &
                 +$1(i,k+1)*pz32u(3))')
#endif


      
module equation
USE comm
USE stepping
USE prob
IMPLICIT NONE

! The governing equations are defined here.    
CONTAINS
!========================================================================
subroutine firststep
  IMPLICIT NONE
  INTEGER :: f

! USE 1st order Euler for the first time step for trapzoidal leapfrog
  !$OMP WORKSHARE
  den=deni
  pux=puxi
  puy=puyi
  puz=puzi
  psi=psii
  pre=prei
#if TWOFLUID==1
  pe=pei
#endif
  by=byi  
  !$OMP END WORKSHARE      

  f=1
  
  call smom2(f)
#if ISOTHERMAL==0
  call spre2(f)
#endif
#if (TWOFLUID==1 && SPE==1 && ISOTHERMAL==0)
  call spe2(f) 
#endif
  call sden2(f)
  call spsi2(f)
#if TWO_AND_HALF==1 
  call sby2(f)
#endif

  !$OMP SINGLE
  time=time+dt
  !$OMP END SINGLE

  call fin

end subroutine firststep
!#######################################################################
subroutine step
  USE omp_lib

  IMPLICIT NONE
  INTEGER :: ns, f
  REAL(RTYPE) :: time_i(0:5)      ! intermediate time   

#if TIMESTEPPING == 0 
! Trapezoidal Leapfrog 
  ns=1
  time_i(0)=time+half*dt
  time_i(1)=time+dt
#endif 

#if TIMESTEPPING == 1
! SSP RK2
  ns=1
  time_i(0)=time+dt
  time_i(1)=time+dt
#endif

#if TIMESTEPPING == 2
! SSP RK3
  ns=2
  time_i(0)=time+dt
  time_i(1)=time+half*dt
  time_i(2)=time+dt
#endif 

#if TIMESTEPPING == 3
! 4-Stage SSP RK3
  ns=3
  time_i(0)=time+half*dt
  time_i(1)=time+dt
  time_i(2)=time+half*dt
  time_i(3)=time+dt
#endif

  
  do f=0, ns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call smom2(f)
#if ISOTHERMAL==0
    call spre2(f)
#endif
#if (TWOFLUID==1 && SPE==1 && ISOTHERMAL==0)
    call spe2(f)
#endif
    call sden2(f)
    call spsi2(f)
#if TWO_AND_HALF==1 
    call sby2(f)
#endif
    time=time_i(f)
    call fin

  end do
  return
end subroutine step
!#######################################################################
subroutine sden2(f)
  IMPLICIT NONE
  INTEGER :: f
  INTEGER :: i, k
        

  !$OMP DO 
  do k=1,nz
  do i=1,nx
    fx(i,k)=-vx(i,k)*deni(i,k)
    fz(i,k)=-vz(i,k)*deni(i,k)
  end do
  end do
  !$OMP END DO

  !$OMP DO 
  do k=3, nz-2
  do i=3, nx-2
    divf(i,k)=DDX(fx)+DDZ(fz)
! diffusion. Beware, this is not physical. Turn it off if you can.
#if DIF_DEN==1
! use 3-pt stencil
    divf(i,k)=divf(i,k)+difden*(D2DX3(deni)+D2DZ3(deni))
#elif DIF_DEN==2
! use 5-pt stencil. Note that higher order stencil is not always 
! "diffusive" and may generate spurious oscillations need steep jumps
    divf(i,k)=divf(i,k)+difden*(D2DX(deni)+D2DZ(deni))
#endif
  end do
  end do
  !$OMP END DO      

#if HDIF_DEN==1
!hyper diffusion
#  if (SIMPLE_LIMITER==0 || LIMITER_DEN==0)    
  call hyperdif(deni,hdifden,dif4den)
#  else
  call hyperdif_limited(deni,hdifden,dif4den)
#  endif
#endif

  call step2d(deni,den,divf,f)

  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
#if DEN_MIN==0
    if (deni(i,k).le.zero.or.deni(i,k).gt.den_max) then
      ifabort=.true.
      print*, 'x=',x(i),'  z=',z(k),'  den=',deni(i,k)
    end if
#elif DEN_MIN==1 
    deni(i,k)=max(deni(i,k),den_min) 
    if (deni(i,k).gt.den_max) then
      ifabort=.true.
      print*, 'x=',x(i),'  z=',z(k),'  den=',deni(i,k)
    end if
#endif
#if TEMP_MIN==1
    prei(i,k)=max(prei(i,k),deni(i,k)*T_min) 
#  if (TWOFLUID==1 && SPE==1)
    pei(i,k)=max(pei(i,k),deni(i,k)*T_min) 
#  endif
#endif
  end do
  end do
  !$OMP END DO


!
  return
end subroutine sden2
!*************************************************************
subroutine spre2(f)
  IMPLICIT NONE
  INTEGER :: f
  INTEGER :: i, k
  REAL(RTYPE) :: bb, bb2, bbx, bbz, k_para1, k_perp1, k_diff, &
                 dTx, dTz, nkbdgT, nkperp



  !$OMP DO 
  do k=1,nz
  do i=1,nx
    fx(i,k)=-vx(i,k)*prei(i,k)
    fz(i,k)=-vz(i,k)*prei(i,k)
  end do
  end do
  !$OMP END DO      

#if THERMAL_CONDUCT==1
! Calculate heat flux, use divfx, divfy, divfz as scratch
  !$OMP DO 
  do k=2,nz-1
  do i=2,nx-1
    bb2=bxi(i,k)**2+byi(i,k)**2+bzi(i,k)**2
    bb=max(sqrt(bb2),eps)
    bb2=max(bb2,eps)
    bbx=bxi(i,k)/bb
    bbz=bzi(i,k)/bb
      ! perpendicular = k_perp/B^2. If this becomes larger than k_para
      ! then isotropic
    k_para1=k_para*Ti(i,k)**2.5    ! parallel conductivity ~ T^(5/2)  
    k_perp1=k_perp*deni(i,k)*deni(i,k)/bb2/Ti(i,k)**0.5
    k_perp1=min(k_para1,k_perp1) 

    ! make thermal conductivity a constant to simplify things
    k_para1=1.00* 3e1
    k_perp1=1.00* 0

    k_diff=k_para1-k_perp1
    dTx=DDX3(Ti)
    dTz=DDZ3(Ti)
    nkbdgT=deni(i,k)*k_diff*(bbx*dTx+bbz*dTz)
    nkperp=deni(i,k)*k_perp1
    divfx(i,k)=-nkbdgT*bbx-nkperp*dTx
    divfz(i,k)=-nkbdgT*bbz-nkperp*dTz
  end do
  end do
  !$OMP END DO
#endif


  !$OMP DO 
  do k=3, nz-2
  do i=3, nx-2
    divf(i,k)=DDX(fx)+DDZ(fz)

!===============pre====================================
    divf(i,k)=divf(i,k)    &
              +omgamma*prei(i,k)*(DDX(vx)+DDZ(vz))
!======================================================
! Ohmic heating, the "half" comes from the fact that heat goes to electron
! and ion equally.  If PE is stepped, then heat goes to electron.
#if OHMIC_HEATING==1 && (TWOFLUID==0 || SPE==0)
    divf(i,k)=divf(i,k)-omgamma*half*electron_heating(i,k) 
#endif

! Viscous heating 
#if VISCOSITY==1 && VISCOUS_HEATING==1
#  if TWOFLUID==0 || SPE==0
! Heating goes to electron and ion equally --> half 
    divf(i,k)=divf(i,k)-omgamma*half*ion_heating(i,k) 
#  else    
    divf(i,k)=divf(i,k)-omgamma*ion_heating(i,k) 
#  endif
#endif




! Thermal exchange between species
#if THERMAL_EXCHANGE==1 && TWOFLUID==1 && SPE==1
    divf(i,k)=divf(i,k)-omgamma*heat_exchange(i,k)
#endif

! Thermal Conductivity
#if THERMAL_CONDUCT==1
#  if TWOFLUID==0 || SPE==0
    ! Distribute evenly between electron and ion -> half
    divf(i,k)=divf(i,k)+omgamma*half*( DDX3(divfx) &
                                          +DDZ3(divfz))
#  else
    divf(i,k)=divf(i,k)+omgamma*( DDX3(divfx) &
                                     +DDZ3(divfz))
#  endif 
#endif


  end do
  end do
  !$OMP END DO      

#if HDIF_PRE==1
!hyper diffusion
#  if (SIMPLE_LIMITER==0 || LIMITER_PRE==0)
  call hyperdif(prei,hdifpre,dif4pre)
#  else
  call hyperdif_limited(prei,hdifpre,dif4pre)
#  endif
!
#endif

  call step2d(prei,pre,divf,f)
!
#if POSITIVE_PRE==1
  !$OMP DO 
  do k=1,nz
  do i=1,nx
    prei(i,k)=max(prei(i,k),zero)
  end do
  end do
  !$OMP END DO
#endif 
  return
end subroutine spre2
!
!*************************************************************
subroutine spe2(f)
  IMPLICIT NONE
  INTEGER :: f
  INTEGER :: i, k

  !$OMP DO 
  do k=1,nz
  do i=1,nx
    fx(i,k)=-vex(i,k)*pei(i,k)
    fz(i,k)=-vez(i,k)*pei(i,k)
  end do
  end do
  !$OMP END DO      

  !$OMP DO 
  do k=3, nz-2
  do i=3, nx-2
    divf(i,k)=DDX(fx)+DDZ(fz)

!===============pe====================================
    divf(i,k)=divf(i,k)    &
              +omgamma*pei(i,k)*(DDX(vex)+DDZ(vez))
!======================================================
#if OHMIC_HEATING==1
    divf(i,k)=divf(i,k)-omgamma*electron_heating(i,k)
#endif
#if THERMAL_EXCHANGE==1
    divf(i,k)=divf(i,k)+omgamma*heat_exchange(i,k)                      
#endif
  end do
  end do
  !$OMP END DO

#if HDIF_PE==1
!hyper diffusion
#  if (SIMPLE_LIMITER==0 || LIMITER_PE==0)
  call hyperdif_e(pei,hdifpe,dif4pe)
#  else
  call hyperdif_e_limited(pei,hdifpe,dif4pe)
#  endif
!
#endif

  call step2d(pei,pe,divf,f)
!

#if POSITIVE_PE==1 
  !$OMP DO 
  do k=1,nz
  do i=1,nx
    pei(i,k)=max(pei(i,k),zero)
  end do
  end do
  !$OMP END DO
#endif

  return
end subroutine spe2
!

!*************************************************************
subroutine smom2(f)
!
  IMPLICIT NONE
  INTEGER :: f
  INTEGER :: i, k
  REAL(RTYPE) :: sft, dxvden, dzvden
! step momenta -- pux, puy, puz

! notice that pt=prei+pei
  !$OMP DO 
  do k=1,nz
  do i=1,nx
    fxx(i,k)=-vx(i,k)*puxi(i,k)-pt(i,k)
    fxz(i,k)=-vz(i,k)*puxi(i,k)
    fzz(i,k)=-vz(i,k)*puzi(i,k)-pt(i,k)
#if TWO_AND_HALF==1
    fxy(i,k)=-vy(i,k)*puxi(i,k)
    fyz(i,k)=-vz(i,k)*puyi(i,k)
#endif
 
#if ((TWOFLUID==0 || HALL_CAL==0) && JXB==0 && AMBIPOLAR==0)    
    ! Add Magnetic Stress Tensor if magnetic force is not pre-calculated
    fxx(i,k)=fxx(i,k)+( bxi(i,k)*bxi(i,k) &
                       -bzi(i,k)*bzi(i,k))*half
    fxz(i,k)=fxz(i,k)+bxi(i,k)*bzi(i,k)
    fzz(i,k)=fzz(i,k)+(-bxi(i,k)*bxi(i,k) &
                       +bzi(i,k)*bzi(i,k))*half
#  if TWO_AND_HALF==1
    fxx(i,k)=fxx(i,k)+(-byi(i,k)*byi(i,k))*half
    fzz(i,k)=fzz(i,k)+(-byi(i,k)*byi(i,k))*half
    fxy(i,k)=fxy(i,k)+bxi(i,k)*byi(i,k)
    fyz(i,k)=fyz(i,k)+byi(i,k)*bzi(i,k)
#  endif
#endif
  end do
  end do
  !$OMP END DO NOWAIT

#if VISCOSITY==1
! Viscous Stress Tensor
  !$OMP DO 
  do k=2,nz-1
  do i=2,nx-1
    vxx(i,k)=DDX3(vx)
    vxz(i,k)=(DDX3(vz)+DDZ3(vx))*half
    vzz(i,k)=DDZ3(vz)
#if TWO_AND_HALF==1
    vxy(i,k)=DDX3(vy)*half
    vyz(i,k)=DDZ3(vy)*half
#endif
    
#  if VISCOUS_HEATING==1
    ion_heating(i,k)= vxx(i,k)*vxx(i,k) &
                     +vzz(i,k)*vzz(i,k) &
                     +two*vxz(i,k)*vxz(i,k) 
#    if TWO_AND_HALF==1
    ion_heating(i,k)= ion_heating(i,k)   &
                     +two*( vxy(i,k)*vxy(i,k) &
                           +vyz(i,k)*vyz(i,k))
#    endif

    ion_heating(i,k)=ion_heating(i,k)*vden(i,k)
#  endif
  end do
  end do
  !$OMP END DO NOWAIT
#endif

#if RANDOM_FORCE==1
  sft=two*rfamp/sqrt(dt)
#endif

  !$OMP BARRIER

  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
    divfx(i,k)=DDX(fxx)+DDZ(fxz)
    divfz(i,k)=DDX(fxz)+DDZ(fzz)
#  if TWO_AND_HALF==1
    divfy(i,k)=DDX(fxy)+DDZ(fyz)
#  endif             

!===============  Momentum ====================================
#if (JXB==1 || AMBIPOLAR==1 || (TWOFLUID==1 && HALL_CAL==1))    
    divfx(i,k)=divfx(i,k)+fbx(i,k)
    divfz(i,k)=divfz(i,k)+fbz(i,k)
#  if TWO_AND_HALF==1
    divfy(i,k)=divfy(i,k)+fby(i,k)
#  endif             
#endif
!======================================================
! Viscosity & friction
#if VISCOSITY==1
    dxvden=DDX3(vden)
    dzvden=DDZ3(vden)
    divfx(i,k)=divfx(i,k) &
               +vden(i,k)*( D2DX3(vx) &
                           +half*( D2DZ3(vx) &
                                  +DDX3(vzz))) &
               +vxx(i,k)*dxvden &
               +vxz(i,k)*dzvden
    divfz(i,k)=divfz(i,k) &
               +vden(i,k)*( D2DZ3(vz) &
                           +half*( D2DX3(vz) &
                                  +DDZ3(vxx))) &
               +vxz(i,k)*dxvden &
               +vzz(i,k)*dzvden
#  if TWO_AND_HALF==1
    divfy(i,k)=divfy(i,k) &
               +vden(i,k)*half*(D2DX3(vy)+D2DZ3(vy)) &
               +vxy(i,k)*dxvden &
               +vyz(i,k)*dzvden
#  endif             
#endif
#if FRICTION==1
    divfx(i,k)=divfx(i,k)-fric*puxi(i,k)
    divfz(i,k)=divfz(i,k)-fric*puzi(i,k)
#  if TWO_AND_HALF==1
    divfy(i,k)=divfy(i,k)-fric*puyi(i,k)
#  endif    
#elif FRICTION==2    
    divfx(i,k)=divfx(i,k)-fric*fric_prof(i,k)*puxi(i,k)
    divfz(i,k)=divfz(i,k)-fric*fric_prof(i,k)*puzi(i,k)
#  if TWO_AND_HALF==1
    divfy(i,k)=divfy(i,k)-fric*fric_prof(i,k)*puyi(i,k)
#  endif    
#endif
#if RANDOM_FORCE==1
    divfx(i,k)=divfx(i,k)+sft*deni(i,k)*(ran(idum)-half)
    divfz(i,k)=divfz(i,k)+sft*deni(i,k)*(ran(idum)-half)
#  if TWO_AND_HALF==1
    divfy(i,k)=divfy(i,k)+sft*deni(i,k)*(ran(idum)-half)
#  endif    
#endif

! Gravity
#if GRAVITY==1
    divfz(i,k)=divfz(i,k)-deni(i,k)*gravity
#endif 

  end do
  end do
  !$OMP END DO


#if HDIF_MOM==1
!hyper diffusion
#  if (SIMPLE_LIMITER==0 || LIMITER_MOM==0)
  call hyperdif_vec(puxi,puyi,puzi,hdifmom,dif4mom)
#  else
  call hyperdif_vec_limited(puxi,puyi,puzi,hdifmom,dif4mom)
#  endif
#endif

!

  call step2d(puxi,pux,divfx,f)
  call step2d(puzi,puz,divfz,f)
#if TWO_AND_HALF==1
  call step2d(puyi,puy,divfy,f)
#endif
!
  return
end subroutine smom2
!



!*************************************************************
subroutine spux2(f)
!
  IMPLICIT NONE
  INTEGER :: f
  INTEGER :: i, k
  REAL(RTYPE) :: sft

! notice that ptb=prei+pei+byi^2/2
  !$OMP DO 
  do k=1,nz
  do i=1,nx
    fx(i,k)=-vx(i,k)*puxi(i,k)-pt(i,k)
    fz(i,k)=-vz(i,k)*puxi(i,k)
#if ((TWOFLUID==0 || HALL_CAL==0) && JXB==0)    
    ! Add Magnetic Stress Tensor if magnetic force is not pre-calculated
#  if TWO_AND_HALF==1    
    fx(i,k)=fx(i,k)+( bxi(i,k)*bxi(i,k)-byi(i,k)*byi(i,k) &
                       -bzi(i,k)*bzi(i,k))*half
    fz(i,k)=fz(i,k)+bxi(i,k)*bzi(i,k)
#  elif TWO_AND_HALF==0
    fx(i,k)=fx(i,k)+( bxi(i,k)*bxi(i,k) &
                       -bzi(i,k)*bzi(i,k))*half
    fz(i,k)=fz(i,k)+bxi(i,k)*bzi(i,k)
#  endif
#endif
    
  end do
  end do
  !$OMP END DO      

#if RANDOM_FORCE==1
  sft=two*rfamp/sqrt(dt)
#endif

  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
    divf(i,k)=DDX(fx)+DDZ(fz)

!===============pux====================================
#if (JXB==1 || (TWOFLUID==1 && HALL_CAL==1))    
    divf(i,k)=divf(i,k)+fbx(i,k)
#endif

!======================================================
! Viscosity & friction
#if VISCOSITY==1
    divf(i,k)=divf(i,k)+visc*(D2DX(puxi)+D2DZ(puxi)) 
#endif
#if FRICTION==1
    divf(i,k)=divf(i,k)-fric*puxi(i,k)
#elif FRICTION==2    
    divf(i,k)=divf(i,k)-fric_prof(i,k)*puxi(i,k)
#endif
#if RANDOM_FORCE==1
    divf(i,k)=divf(i,k)+sft*deni(i,k)*(ran(idum)-half)
#endif
  end do
  end do
  !$OMP END DO

#if HDIF_MOM==1
!hyper diffusion
#  if (SIMPLE_LIMITER==0 || LIMITER_MOM==0)
  call hyperdif(puxi,hdifmom,dif4mom)
#  else
  call hyperdif_limited(puxi,hdifmom,dif4mom)
#  endif
#endif

!
  call step2d(puxi,pux,divf,f)

!
  return
end subroutine spux2
!
!*************************************************************
subroutine spuy2(f)
  IMPLICIT NONE
  INTEGER :: f
  INTEGER :: i, k
  REAL(RTYPE) :: sft

  !$OMP DO 
  do k=1, nz
  do i=1, nx
    fx(i,k)=-vx(i,k)*puyi(i,k)
    fz(i,k)=-vz(i,k)*puyi(i,k)
#if ((TWOFLUID==0 || HALL_CAL==0) && JXB==0)    
    ! Add Magnetic Stress Tensor if magnetic force is not pre-calculated
    fx(i,k)=fx(i,k)+bxi(i,k)*byi(i,k)
    fz(i,k)=fz(i,k)+byi(i,k)*bzi(i,k)
#endif
    
  end do
  end do
  !$OMP END DO
#if RANDOM_FORCE==1
  sft=two*rfamp/sqrt(dt)
#endif

  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
    divf(i,k)=DDX(fx)+DDZ(fz)
!==================puy=========================================
#if (JXB==1 || (TWOFLUID==1 && HALL_CAL==1))    
    divf(i,k)=divf(i,k)+fby(i,k)
#endif
!==============================================================
! Viscosity & friction
#if VISCOSITY==1
    divf(i,k)=divf(i,k)+visc*(D2DX(puyi)+D2DZ(puyi)) 
#endif
#if FRICTION==1
    divf(i,k)=divf(i,k)-fric*puyi(i,k)
#elif FRICTION==2    
    divf(i,k)=divf(i,k)-fric_prof(i,k)*puyi(i,k)
#endif
#if RANDOM_FORCE==1
    divf(i,k)=divf(i,k)+sft*deni(i,k)*(ran(idum)-half)
#endif
  end do
  end do
!$OMP END DO

#if HDIF_MOM==1
!hyper diffusion
#  if (SIMPLE_LIMITER==0 || LIMITER_MOM==0)
  call hyperdif(puyi,hdifmom,dif4mom)
#  else
  call hyperdif_limited(puyi,hdifmom,dif4mom)
#  endif
#endif


  call step2d(puyi,puy,divf,f)
!
  return
end subroutine spuy2
!


!*************************************************************
subroutine spuz2(f)
  IMPLICIT NONE
  INTEGER :: f
  INTEGER :: i, k
  REAL(RTYPE) :: sft

! notice that ptb=prei+pei+byi^2/2

!$OMP DO 
  do k=1, nz
  do i=1, nx
    fx(i,k)=-vx(i,k)*puzi(i,k)
    fz(i,k)=-vz(i,k)*puzi(i,k)-pt(i,k)
#if ((TWOFLUID==0 || HALL_CAL==0) && JXB==0)    
    ! Add Magnetic Stress Tensor if magnetic force is not pre-calculated
#  if TWO_AND_HALF==1    
    fx(i,k)=fx(i,k)+bxi(i,k)*bzi(i,k)
    fz(i,k)=fz(i,k)+(-bxi(i,k)*bxi(i,k)-byi(i,k)*byi(i,k) &
                     +bzi(i,k)*bzi(i,k))*half
#  elif TWO_AND_HALF==0
    fx(i,k)=fx(i,k)+bxi(i,k)*bzi(i,k)
    fz(i,k)=fz(i,k)+(-bxi(i,k)*bxi(i,k) &
                     +bzi(i,k)*bzi(i,k))*half
#  endif
#endif
    
  end do
  end do
!$OMP END DO

#if RANDOM_FORCE==1
  sft=two*rfamp/sqrt(dt)
#endif

  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
    divf(i,k)=DDX(fx)+DDZ(fz)

!===============puz====================================
#if (JXB==1 || (TWOFLUID==1 && HALL_CAL==1))    
    divf(i,k)=divf(i,k)+fbz(i,k)
#endif
!======================================================
! Viscosity & friction
#if VISCOSITY==1
    divf(i,k)=divf(i,k)+visc*(D2DX(puzi)+D2DZ(puzi)) 
#endif
#if FRICTION==1
    divf(i,k)=divf(i,k)-fric*puzi(i,k)
#elif FRICTION==2    
    divf(i,k)=divf(i,k)-fric_prof(i,k)*puzi(i,k)
#endif
#if GRAVITY==1
    divf(i,k)=divf(i,k)-deni(i,k)*gravity
#endif 
#if RANDOM_FORCE==1
    divf(i,k)=divf(i,k)+sft*deni(i,k)*(ran(idum)-half)
#endif
  end do
  end do
  !$OMP END DO

#if HDIF_MOM==1
!hyper diffusion
#  if (SIMPLE_LIMITER==0 || LIMITER_MOM==0)
  call hyperdif(puzi,hdifmom,dif4mom)
#  else
  call hyperdif_limited(puzi,hdifmom,dif4mom)  
#  endif
#endif

  call step2d(puzi,puz,divf,f)

!

  return
end subroutine spuz2
!
!*************************************************************
subroutine sby2(f)
!
  IMPLICIT NONE
  INTEGER :: f
  INTEGER :: i, k

  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
!=============== curl E =====================================
    divf(i,k)=DDX(ez)-DDZ(ex)
!============================================================
#if (RESISTIVITY==1 && B_DIFF==1)
    divf(i,k)=divf(i,k)+eta*(D2DX(byi)+D2DZ(byi))
#elif (RESISTIVITY==1 && B_DIFF==2)
    divf(i,k)=divf(i,k)+eta*(D2DX3(byi)+D2DZ3(byi))
#endif

  end do
  end do
  !$OMP END DO

!hyper diffusion
#if (HDIF_B==1 && TWOFLUID==1)
#  if (SIMPLE_LIMITER==0 || LIMITER_B==0)  
  call hyperdif_e(byi,hdifb,dif4b)
#  else
  call hyperdif_e_limited(byi,hdifb,dif4b)
#  endif 
#elif (HDIF_B==1 && TWOFLUID==0)
#  if (SIMPLE_LIMITER==0 || LIMITER_B==0)  
  call hyperdif(byi,hdifb,dif4b)
#  else
  call hyperdif_limited(byi,hdifb,dif4b)
#  endif
#endif

  call step2d(byi,by,divf,f)

!
  return
end subroutine sby2
! 
!*************************************************************
subroutine spsi2(f)
  IMPLICIT NONE
  INTEGER :: f
  INTEGER :: i,k

  !$OMP DO 
  do k=3,nz-2
  do i=3, nx-2
    divf(i,k)=-ey(i,k)
#if (RESISTIVITY==1 && B_DIFF==1)
    divf(i,k)=divf(i,k)+eta*(D2DX(psii)+D2DZ(psii))
#elif (RESISTIVITY==1 && B_DIFF==2)
    divf(i,k)=divf(i,k)+eta*(D2DX3(psii)+D2DZ3(psii))
#endif
  end do
  end do
  !$OMP END DO

!hyper diffusion
#if (HDIF_B==1 && TWOFLUID==1)
#  if (SIMPLE_LIMITER==0 || LIMITER_B==0)  
  call hyperdif_e(psii,hdifb,dif4b)
#  else
  call hyperdif_e_limited(psii,hdifb,dif4b)
#  endif
#endif
#if (HDIF_B==1 && TWOFLUID==0) 
#  if (SIMPLE_LIMITER==0 || LIMITER_B==0)  
  call hyperdif(psii,hdifb,dif4b)
#  else
  call hyperdif_limited(psii,hdifb,dif4b)
#  endif
#endif

!
  call step2d(psii,psi,divf,f)

!
  return
end subroutine spsi2
!
!*************************************************************
subroutine hyperdif(foo,hypdif,dif4)
! add hyperdiffusion to divf
  IMPLICIT NONE
  REAL(RTYPE) :: foo(nx,nz), dif4, hypdif
  INTEGER :: i, k
  REAL(RTYPE) :: dif4x, dif4z

! Hyper diffusion, see Guzdar et. at. 1993.
! For non-uniform grids, we don't try to calculate it precisely, since
! it is hyper. We assume the grid size does not change much.
! hypdif is the coefficient for the velocity dependent hyper diffusion
! dif4 is the coefficient for 4th order hyper diffusion.
! The hyperdiffusion coefficient is approximately (dx^3)(2*|vx|*hypdif+dif4)/24
! along x and  (dz^3)(2*|vz|*hypdif+dif4)/24 along z

  !$OMP DO 
  do k=2,nz-2
    do i=2,nx-2
#if (VDEPT_HDIF>0 || HDIF4>0)                
#  if (VDEPT_HDIF>0 && HDIF4==1)
      dif4x=dif4+hypdif*vxm(i,k)
      dif4z=dif4+hypdif*vzm(i,k)
#  elif (VDEPT_HDIF>0 && HDIF4==2) 
      dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half+hypdif*vxm(i,k)
      dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half+hypdif*vzm(i,k)
#  elif (VDEPT_HDIF>0)
      dif4x=hypdif*vxm(i,k)
      dif4z=hypdif*vzm(i,k)
#  elif (HDIF4==1)
      dif4x=dif4
      dif4z=dif4
#  elif (HDIF4==2)
      dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half
      dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half
#  endif        
      fx(i,k)=foo(i+2,k)-foo(i-1,k)-3._RTYPE*(foo(i+1,k)-foo(i+0,k))
      fz(i,k)=foo(i,k+2)-foo(i,k-1)-3._RTYPE*(foo(i,k+1)-foo(i,k+0))
      fx(i,k)=fx(i,k)*dif4x
      fz(i,k)=fz(i,k)*dif4z
#else
      fx(i,k)=zero
      fz(i,k)=zero
#endif
    end do
  end do
  !$OMP END DO

! Although this part may be fused into the previous loop to 
! be more cache friendly, it will cause problem when OpenMP
! is used.

  !$OMP DO 
  do k=3,nz-2
    do i=3,nx-2
      divf(i,k)=divf(i,k)-((fx(i,k)-fx(i-1,k))*dx24(i)  &
                           +(fz(i,k)-fz(i,k-1))*dz24(k))
    end do
  end do
  !$OMP END DO


end subroutine hyperdif
!########################################################################

!*************************************************************
subroutine hyperdif_e(foo,hypdif,dif4)
! add hyperdiffusion to divf
! same as hyperdif except using ve instead of v.
  IMPLICIT NONE
  REAL(RTYPE) :: foo(nx,nz), dif4, hypdif
  INTEGER :: i, k
  REAL(RTYPE) :: dif4x, dif4z

! Hyper diffusion, see Guzdar et. at. 1993.
! For non-uniform grids, we don't try to calculate it precisely, since
! it is hyper. We assume the grid size does not change much.
! hypdif is the coefficient for the velocity dependent hyper diffusion
! dif4 is the coefficient for 4th order hyper diffusion.
! The hyperdiffusion coefficient is approximately (dx^3)(2*|vex|*hypdif+dif4)/24
! along x and  (dz^3)(2*|vez|*hypdif+dif4)/24 along z

  !$OMP DO 
  do k=2,nz-2
    do i=2,nx-2
#if (VDEPT_HDIF>0 || HDIF4>0)                
#  if (VDEPT_HDIF>0 && HDIF4==1)
      dif4x=dif4+hypdif*vexm(i,k)
      dif4z=dif4+hypdif*vezm(i,k)
#  elif (VDEPT_HDIF>0 && HDIF4==2) 
      dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half &
            +hypdif*vexm(i,k)
      dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half &
            +hypdif*vezm(i,k)
#  elif (VDEPT_HDIF>0)
      dif4x=hypdif*vexm(i,k)
      dif4z=hypdif*vezm(i,k)
#  elif (HDIF4==1)
      dif4x=dif4
      dif4z=dif4
#  elif (HDIF4==2)
      dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half
      dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half
#  endif        
      fx(i,k)=foo(i+2,k)-foo(i-1,k)-3._RTYPE*(foo(i+1,k)-foo(i+0,k))
      fz(i,k)=foo(i,k+2)-foo(i,k-1)-3._RTYPE*(foo(i,k+1)-foo(i,k+0))
      fx(i,k)=fx(i,k)*dif4x
      fz(i,k)=fz(i,k)*dif4z
#else
      fx(i,k)=zero
      fz(i,k)=zero
#endif
    end do
  end do
  !$OMP END DO      

! Although this part may be fused into the previous loop to 
! be more cache friendly, it will cause problem when OpenMP
! is used.

  !$OMP DO  
  do k=3,nz-2
    do i=3,nx-2
      divf(i,k)=divf(i,k)-((fx(i,k)-fx(i-1,k))*dx24(i)  &
                           +(fz(i,k)-fz(i,k-1))*dz24(k))
    end do
  end do
  !$OMP END DO

end subroutine hyperdif_e
!*************************************************************
subroutine hyperdif_limited(foo,hypdif,dif4)
! add hyperdiffusion to divf 
  IMPLICIT NONE
  REAL(RTYPE) :: foo(nx,nz), dif4, hypdif
  INTEGER :: i, k
  REAL(RTYPE) :: dif4x, dif4z

! Hyper diffusion, see Guzdar et. at. 1993.
! For non-uniform grids, we don't try to calculate it precisely, since
! it is hyper. We assume the grid size does not change much.
! hypdif is the coefficient for the velocity dependent hyper diffusion
! dif4 is the coefficient for 4th order hyper diffusion.
! The hyperdiffusion coefficient is approximately (dx^3)(2*|vx|*hypdif+dif4)/24
! along x and  (dz^3)(2*|vz|*hypdif+dif4)/24 along z
! Simple limiter is applied when SIMPLE_LIMITER>0 

  !$OMP DO 
  do k=2,nz-2
    do i=2,nx-2
#if (VDEPT_HDIF>0 || HDIF4>0)                
#  if (VDEPT_HDIF>0 && HDIF4==1)
      dif4x=dif4+hypdif*vxm(i,k)
      dif4z=dif4+hypdif*vzm(i,k)
#  elif (VDEPT_HDIF>0 && HDIF4==2) 
      dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half+hypdif*vxm(i,k)
      dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half+hypdif*vzm(i,k)
#  elif (VDEPT_HDIF>0)
      dif4x=hypdif*vxm(i,k)
      dif4z=hypdif*vzm(i,k)
#  elif (HDIF4==1)
      dif4x=dif4
      dif4z=dif4
#  elif (HDIF4==2)
      dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half
      dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half
#  endif        
      fx(i,k)=foo(i+2,k)-foo(i-1,k)-3._RTYPE*(foo(i+1,k)-foo(i+0,k))
      fz(i,k)=foo(i,k+2)-foo(i,k-1)-3._RTYPE*(foo(i,k+1)-foo(i,k+0))
      fx(i,k)=fx(i,k)*dif4x
      fz(i,k)=fz(i,k)*dif4z
#  if SIMPLE_LIMITER==1
      if (fx(i,k)*(foo(i+1,k)-foo(i,k)).gt.zero) then
        if (fx(i,k)*(foo(i+2,k)-foo(i+1,k)).lt.zero.or. &
            fx(i,k)*(foo(i,k)-foo(i-1,k)).lt.zero) then
          fx(i,k)=zero
        end if
      end if
      if (fz(i,k)*(foo(i,k+1)-foo(i,k)).gt.zero) then
        if (fz(i,k)*(foo(i,k+2)-foo(i,k+1)).lt.zero.or. &
            fz(i,k)*(foo(i,k)-foo(i,k-1)).lt.zero) then
          fz(i,k)=zero
        end if
      end if
#  elif SIMPLE_LIMITER==2
      if (fx(i,k)*(foo(i+1,k)-foo(i,k)).gt.zero) then
        fx(i,k)=zero
      end if
      if (fz(i,k)*(foo(i,k+1)-foo(i,k)).gt.zero) then
        fz(i,k)=zero
      end if
#  endif

#else
      fx(i,k)=zero
      fz(i,k)=zero
#endif
    end do
  end do
  !$OMP END DO

! Although this part may be fused into the previous loop to 
! be more cache friendly, it will cause problem when OpenMP
! is used.

  !$OMP DO 
  do k=3,nz-2
    do i=3,nx-2
      divf(i,k)=divf(i,k)-((fx(i,k)-fx(i-1,k))*dx24(i)  &
                           +(fz(i,k)-fz(i,k-1))*dz24(k))
    end do
  end do
  !$OMP END DO

end subroutine hyperdif_limited
!########################################################################

!*************************************************************
subroutine hyperdif_e_limited(foo,hypdif,dif4)
! add hyperdiffusion to divf
! same as hyperdif except using ve instead of v.
  IMPLICIT NONE
  REAL(RTYPE) :: foo(nx,nz), dif4, hypdif
  INTEGER :: i, k
  REAL(RTYPE) :: dif4x, dif4z

! Hyper diffusion, see Guzdar et. at. 1993.
! For non-uniform grids, we don't try to calculate it precisely, since
! it is hyper. We assume the grid size does not change much.
! hypdif is the coefficient for the velocity dependent hyper diffusion
! dif4 is the coefficient for 4th order hyper diffusion.
! The hyperdiffusion coefficient is approximately (dx^3)(2*|vex|*hypdif+dif4)/24
! along x and  (dz^3)(2*|vez|*hypdif+dif4)/24 along z
! Simple limiter is applied when SIMPLE_LIMITER>0

  !$OMP DO 
  do k=2,nz-2
    do i=2,nx-2
#if (VDEPT_HDIF>0 || HDIF4>0)                
#  if (VDEPT_HDIF>0 && HDIF4==1)
      dif4x=dif4+hypdif*vexm(i,k)
      dif4z=dif4+hypdif*vezm(i,k)
#  elif (VDEPT_HDIF>0 && HDIF4==2) 
      dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half &
            +hypdif*vexm(i,k)
      dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half &
            +hypdif*vezm(i,k)
#  elif (VDEPT_HDIF>0)
      dif4x=hypdif*vexm(i,k)
      dif4z=hypdif*vezm(i,k)
#  elif (HDIF4==1)
      dif4x=dif4
      dif4z=dif4
#  elif (HDIF4==2)
      dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half
      dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half
#  endif        
      fx(i,k)=foo(i+2,k)-foo(i-1,k)-3._RTYPE*(foo(i+1,k)-foo(i+0,k))
      fz(i,k)=foo(i,k+2)-foo(i,k-1)-3._RTYPE*(foo(i,k+1)-foo(i,k+0))
      fx(i,k)=fx(i,k)*dif4x
      fz(i,k)=fz(i,k)*dif4z
#  if SIMPLE_LIMITER==1
      if (fx(i,k)*(foo(i+1,k)-foo(i,k)).gt.zero) then
        if (fx(i,k)*(foo(i+2,k)-foo(i+1,k)).lt.zero.or. &
            fx(i,k)*(foo(i,k)-foo(i-1,k)).lt.zero) then
          fx(i,k)=zero
        end if
      end if
      if (fz(i,k)*(foo(i,k+1)-foo(i,k)).gt.zero) then
        if (fz(i,k)*(foo(i,k+2)-foo(i,k+1)).lt.zero.or. &
            fz(i,k)*(foo(i,k)-foo(i,k-1)).lt.zero) then
          fz(i,k)=zero
        end if
      end if
#  elif SIMPLE_LIMITER==2
      if (fx(i,k)*(foo(i+1,k)-foo(i,k)).gt.zero) then
        fx(i,k)=zero
      end if
      if (fz(i,k)*(foo(i,k+1)-foo(i,k)).gt.zero) then
        fz(i,k)=zero
      end if
#  endif

#else
      fx(i,k)=zero
      fz(i,k)=zero
#endif
    end do
  end do
  !$OMP END DO

! Although this part may be fused into the previous loop to 
! be more cache friendly, it will cause problem when OpenMP
! is used.

  !$OMP DO 
  do k=3,nz-2
    do i=3,nx-2
      divf(i,k)=divf(i,k)-((fx(i,k)-fx(i-1,k))*dx24(i)  &
                           +(fz(i,k)-fz(i,k-1))*dz24(k))
    end do
  end do
  !$OMP END DO

end subroutine hyperdif_e_limited

!*************************************************************
subroutine hyperdif_vec(fox,foy,foz,hypdif,dif4)
! Vector version of hyperdif  
! add hyperdiffusion to divfx, divfy, divfz
  USE comm
  IMPLICIT NONE
  REAL(RTYPE) :: fox(nx,nz), foy(nx,nz), foz(nx,nz)
  REAL(RTYPE) :: dif4, hypdif
  INTEGER :: i, k
  REAL(RTYPE) :: dif4x, dif4z


! Hyper diffusion, see Guzdar et. at. 1993.
! For non-uniform grids, we don't try to calculate it precisely, since
! it is hyper. We assume the grid size does not change much.
! hypdif is the coefficient for the velocity dependent hyper diffusion
! dif4 is the coefficient for 4th order hyper diffusion.
! The hyperdiffusion coefficient is approximately (dx^3)(2*|vx|*hypdif+dif4)/24
! along x, and (dz^3)(2*|vz|*hypdif+dif4)/24 along z


#if (VDEPT_HDIF>0 || HDIF4>0)                
  !$OMP DO 
  do k=2,nz-2
  do i=2,nx-2
#  if (VDEPT_HDIF>0 && HDIF4==1)
    dif4x=dif4+hypdif*vxm(i,k)
    dif4z=dif4+hypdif*vzm(i,k)
#  elif (VDEPT_HDIF>0 && HDIF4==2) 
    dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half &
          +hypdif*vxm(i,k)
    dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half &
          +hypdif*vzm(i,k)
#  elif (VDEPT_HDIF>0)
    dif4x=hypdif*vxm(i,k)
    dif4z=hypdif*vzm(i,k)
#  elif (HDIF4==1)
    dif4x=dif4
    dif4z=dif4
#  elif (HDIF4==2)
    dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half
    dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half
#  endif               
    fxx(i,k)=fox(i+2,k)-fox(i-1,k)-3._RTYPE*(fox(i+1,k) &
              -fox(i+0,k))
    fzx(i,k)=fox(i,k+2)-fox(i,k-1)-3._RTYPE*(fox(i,k+1) &
              -fox(i,k+0))
    fxx(i,k)=fxx(i,k)*dif4x
    fzx(i,k)=fzx(i,k)*dif4z

    fxy(i,k)=foy(i+2,k)-foy(i-1,k)-3._RTYPE*(foy(i+1,k) &
              -foy(i+0,k))
    fzy(i,k)=foy(i,k+2)-foy(i,k-1)-3._RTYPE*(foy(i,k+1) &
              -foy(i,k+0))
    fxy(i,k)=fxy(i,k)*dif4x
    fzy(i,k)=fzy(i,k)*dif4z

    fxz(i,k)=foz(i+2,k)-foz(i-1,k)-3._RTYPE*(foz(i+1,k) &
              -foz(i+0,k))
    fzz(i,k)=foz(i,k+2)-foz(i,k-1)-3._RTYPE*(foz(i,k+1) &
              -foz(i,k+0))
    fxz(i,k)=fxz(i,k)*dif4x
    fzz(i,k)=fzz(i,k)*dif4z
  end do
  end do
  !$OMP END DO

! Although this part may be fused into the previous loop to 
! be more cache friendly, it will cause problem when OpenMP
! is used.
  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
    divfx(i,k)=divfx(i,k)-((fxx(i,k)-fxx(i-1,k))*dx24(i)  &
                          +(fzx(i,k)-fzx(i,k-1))*dz24(k))
    divfz(i,k)=divfz(i,k)-((fxz(i,k)-fxz(i-1,k))*dx24(i)  &
                          +(fzz(i,k)-fzz(i,k-1))*dz24(k))
#  if TWO_AND_HALF==1
    divfy(i,k)=divfy(i,k)-((fxy(i,k)-fxy(i-1,k))*dx24(i)  &
                          +(fzy(i,k)-fzy(i,k-1))*dz24(k))
#  endif
  end do
  end do
  !$OMP END DO
#endif
end subroutine hyperdif_vec
!########################################################################
subroutine hyperdif_vec_e(fox,foy,foz,hypdif,dif4)
! Vector version of hyperdif_e  
! add hyperdiffusion to divfx, divfy, divfz
  USE comm
  IMPLICIT NONE
  REAL(RTYPE) :: fox(nx,nz), foy(nx,nz), foz(nx,nz)
  REAL(RTYPE) :: dif4, hypdif
 
  INTEGER :: i, k
  REAL(RTYPE) :: dif4x, dif4z


! Hyper diffusion, see Guzdar et. at. 1993.
! For non-uniform grids, we don't try to calculate it precisely, since
! it is hyper. We assume the grid size does not change much.
! hypdif is the coefficient for the velocity dependent hyper diffusion
! dif4 is the coefficient for 4th order hyper diffusion.
! The hyperdiffusion coefficient is approximately (dx^3)(2*|vex|*hypdif+dif4)/24
! along x, and (dz^3)(2*|vez|*hypdif+dif4)/24 along z


#if (VDEPT_HDIF>0 || HDIF4>0)                
  !$OMP DO 
  do k=2,nz-2
  do i=2,nx-2
#  if (VDEPT_HDIF>0 && HDIF4==1)
    dif4x=dif4+hypdif*vexm(i,k)
    dif4z=dif4+hypdif*vezm(i,k)
#  elif (VDEPT_HDIF>0 && HDIF4==2) 
    dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half &
          +hypdif*vexm(i,k)
    dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half &
          +hypdif*vezm(i,k)
#  elif (VDEPT_HDIF>0)
    dif4x=hypdif*vexm(i,k)
    dif4z=hypdif*vezm(i,k)
#  elif (HDIF4==1)
    dif4x=dif4
    dif4z=dif4
#  elif (HDIF4==2)
    dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half
    dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half
#  endif               
    fxx(i,k)=fox(i+2,k)-fox(i-1,k)-3._RTYPE*(fox(i+1,k) &
            -fox(i+0,k))
    fzx(i,k)=fox(i,k+2)-fox(i,k-1)-3._RTYPE*(fox(i,k+1) &
            -fox(i,k+0))
    fxx(i,k)=fxx(i,k)*dif4x
    fzx(i,k)=fzx(i,k)*dif4z

    fxy(i,k)=foy(i+2,k)-foy(i-1,k)-3._RTYPE*(foy(i+1,k) &
            -foy(i+0,k))
    fzy(i,k)=foy(i,k+2)-foy(i,k-1)-3._RTYPE*(foy(i,k+1) &
            -foy(i,k+0))
    fxy(i,k)=fxy(i,k)*dif4x
    fzy(i,k)=fzy(i,k)*dif4z

    fxz(i,k)=foz(i+2,k)-foz(i-1,k)-3._RTYPE*(foz(i+1,k) &
            -foz(i+0,k))
    fzz(i,k)=foz(i,k+2)-foz(i,k-1)-3._RTYPE*(foz(i,k+1) &
            -foz(i,k+0))
    fxz(i,k)=fxz(i,k)*dif4x
    fzz(i,k)=fzz(i,k)*dif4z
  end do
  end do
  !$OMP END DO

! Although this part may be fused into the previous loop to 
! be more cache friendly, it will cause problem when OpenMP
! is used.
  !$OMP DO
  do k=3,nz-2
  do i=3,nx-2
    divfx(i,k)=divfx(i,k)-((fxx(i,k)-fxx(i-1,k))*dx24(i)  &
                          +(fzx(i,k)-fzx(i,k-1))*dz24(k))
    divfz(i,k)=divfz(i,k)-((fxz(i,k)-fxz(i-1,k))*dx24(i)  &
                          +(fzz(i,k)-fzz(i,k-1))*dz24(k))
#  if TWO_AND_HALF==1
    divfy(i,k)=divfy(i,k)-((fxy(i,k)-fxy(i-1,k))*dx24(i)  &
                          +(fzy(i,k)-fzy(i,k-1))*dz24(k))
#  endif
  end do
  end do
  !$OMP END DO
#endif

end subroutine hyperdif_vec_e
!########################################################################
subroutine hyperdif_vec_limited(fox,foy,foz,hypdif,dif4)
  
! Vector version of hyperdif_limited  
! add hyperdiffusion to divfx, divfy, divfz

  USE comm
  IMPLICIT NONE
  REAL(RTYPE) :: fox(nx,nz), foy(nx,nz), foz(nx,nz)
  REAL(RTYPE) :: dif4, hypdif

  INTEGER :: i, k
  REAL(RTYPE) :: dif4x, dif4z


! Hyper diffusion, see Guzdar et. at. 1993.
! For non-uniform grids, we don't try to calculate it precisely, since
! it is hyper. We assume the grid size does not change much.
! hypdif is the coefficient for the velocity dependent hyper diffusion
! dif4 is the coefficient for 4th order hyper diffusion.
! The hyperdiffusion coefficient is approximately (dx^3)(2*|vx|*hypdif+dif4)/24
! along x, and (dz^3)(2*|vz|*hypdif+dif4)/24 along z
! Simple limiter is applied when SIMPLE_LIMITER>0

#if (VDEPT_HDIF>0 || HDIF4>0)                
  !$OMP DO 
  do k=2,nz-2
  do i=2,nx-2
#  if (VDEPT_HDIF>0 && HDIF4==1)
    dif4x=dif4+hypdif*vxm(i,k)
    dif4z=dif4+hypdif*vzm(i,k)
#  elif (VDEPT_HDIF>0 && HDIF4==2) 
    dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half &
          +hypdif*vxm(i,k)
    dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half &
          +hypdif*vzm(i,k)
#  elif (VDEPT_HDIF>0)
    dif4x=hypdif*vxm(i,k)
    dif4z=hypdif*vzm(i,k)
#  elif (HDIF4==1)
    dif4x=dif4
    dif4z=dif4
#  elif (HDIF4==2)
    dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half
    dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half
#  endif               
    fxx(i,k)=fox(i+2,k)-fox(i-1,k)-3._RTYPE*(fox(i+1,k) &
              -fox(i+0,k))
    fzx(i,k)=fox(i,k+2)-fox(i,k-1)-3._RTYPE*(fox(i,k+1) &
              -fox(i,k+0))
    fxx(i,k)=fxx(i,k)*dif4x
    fzx(i,k)=fzx(i,k)*dif4z

    fxy(i,k)=foy(i+2,k)-foy(i-1,k)-3._RTYPE*(foy(i+1,k) &
              -foy(i+0,k))
    fzy(i,k)=foy(i,k+2)-foy(i,k-1)-3._RTYPE*(foy(i,k+1) &
              -foy(i,k+0))
    fxy(i,k)=fxy(i,k)*dif4x
    fzy(i,k)=fzy(i,k)*dif4z

    fxz(i,k)=foz(i+2,k)-foz(i-1,k)-3._RTYPE*(foz(i+1,k) &
              -foz(i+0,k))
    fzz(i,k)=foz(i,k+2)-foz(i,k-1)-3._RTYPE*(foz(i,k+1) &
              -foz(i,k+0))
    fxz(i,k)=fxz(i,k)*dif4x
    fzz(i,k)=fzz(i,k)*dif4z

#  if SIMPLE_LIMITER==1
    if (fxx(i,k)*(fox(i+1,k)-fox(i,k)).gt.zero) then
      if (fxx(i,k)*(fox(i+2,k)-fox(i+1,k)).lt.zero.or. &
        fxx(i,k)*(fox(i,k)-fox(i-1,k)).lt.zero) then
        fxx(i,k)=zero
      end if
    end if
    if (fzx(i,k)*(fox(i,k+1)-fox(i,k)).gt.zero) then
      if (fzx(i,k)*(fox(i,k+2)-fox(i,k+1)).lt.zero.or. &
          fzx(i,k)*(fox(i,k)-fox(i,k-1)).lt.zero) then
        fzx(i,k)=zero
      end if
    end if

    if (fxy(i,k)*(foy(i+1,k)-foy(i,k)).gt.zero) then
      if (fxy(i,k)*(foy(i+2,k)-foy(i+1,k)).lt.zero.or. &
        fxy(i,k)*(foy(i,k)-foy(i-1,k)).lt.zero) then
        fxy(i,k)=zero
      end if
    end if
    if (fzy(i,k)*(foy(i,k+1)-foy(i,k)).gt.zero) then
      if (fzy(i,k)*(foy(i,k+2)-foy(i,k+1)).lt.zero.or. &
          fzy(i,k)*(foy(i,k)-foy(i,k-1)).lt.zero) then
        fzy(i,k)=zero
      end if
    end if

    if (fxz(i,k)*(foz(i+1,k)-foz(i,k)).gt.zero) then
      if (fxz(i,k)*(foz(i+2,k)-foz(i+1,k)).lt.zero.or. &
        fxz(i,k)*(foz(i,k)-foz(i-1,k)).lt.zero) then
        fxz(i,k)=zero
      end if
    end if
    if (fzz(i,k)*(foz(i,k+1)-foz(i,k)).gt.zero) then
      if (fzz(i,k)*(foz(i,k+2)-foz(i,k+1)).lt.zero.or. &
          fzz(i,k)*(foz(i,k)-foz(i,k-1)).lt.zero) then
        fzz(i,k)=zero
      end if
    end if
    
#  elif SIMPLE_LIMITER==2
    if (fxx(i,k)*(fox(i+1,k)-fox(i,k)).gt.zero) then
      fxx(i,k)=zero
    end if
    if (fzx(i,k)*(fox(i,k+1)-fox(i,k)).gt.zero) then
      fzx(i,k)=zero
    end if

    if (fxy(i,k)*(foy(i+1,k)-foy(i,k)).gt.zero) then
      fxy(i,k)=zero
    end if
    if (fzy(i,k)*(foy(i,k+1)-foy(i,k)).gt.zero) then
      fzy(i,k)=zero
    end if

    if (fxz(i,k)*(foz(i+1,k)-foz(i,k)).gt.zero) then
      fxz(i,k)=zero
    end if
    if (fzz(i,k)*(foz(i,k+1)-foz(i,k)).gt.zero) then
      fzz(i,k)=zero
    end if
    
#  endif
   
  end do
  end do
  !$OMP END DO

! Although this part may be fused into the previous loop to 
! be more cache friendly, it will cause problem when OpenMP
! is used.
  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
    divfx(i,k)=divfx(i,k)-((fxx(i,k)-fxx(i-1,k))*dx24(i)  &
                          +(fzx(i,k)-fzx(i,k-1))*dz24(k))
    divfz(i,k)=divfz(i,k)-((fxz(i,k)-fxz(i-1,k))*dx24(i)  &
                          +(fzz(i,k)-fzz(i,k-1))*dz24(k))
#  if TWO_AND_HALF==1
    divfy(i,k)=divfy(i,k)-((fxy(i,k)-fxy(i-1,k))*dx24(i)  &
                          +(fzy(i,k)-fzy(i,k-1))*dz24(k))
#  endif
  end do
  end do
  !$OMP END DO
#endif

end subroutine hyperdif_vec_limited
!########################################################################
subroutine hyperdif_vec_e_limited(fox,foy,foz,hypdif,dif4)
! Vector version of hyperdif_e_limited  
! add hyperdiffusion to divfx, divfy, divfz
  USE comm
  IMPLICIT NONE
  REAL(RTYPE) :: fox(nx,nz), foy(nx,nz), foz(nx,nz)
  REAL(RTYPE) :: dif4, hypdif

  INTEGER :: i, k
  REAL(RTYPE) :: dif4x, dif4z


! Hyper diffusion, see Guzdar et. at. 1993.
! For non-uniform grids, we don't try to calculate it precisely, since
! it is hyper. We assume the grid size does not change much.
! hypdif is the coefficient for the velocity dependent hyper diffusion
! dif4 is the coefficient for 4th order hyper diffusion.
! The hyperdiffusion coefficient is approximately (dx^3)(2*|vex|*hypdif+dif4)/24
! along x, and (dz^3)(2*|vez|*hypdif+dif4)/24 along z
! Simple limiter is applied when SIMPLE_LIMITER>0

#if (VDEPT_HDIF>0 || HDIF4>0)                
  !$OMP DO 
  do k=2,nz-2
  do i=2,nx-2
#  if (VDEPT_HDIF>0 && HDIF4==1)
    dif4x=dif4+hypdif*vexm(i,k)
    dif4z=dif4+hypdif*vezm(i,k)
#  elif (VDEPT_HDIF>0 && HDIF4==2) 
    dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half &
          +hypdif*vexm(i,k)
    dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half &
          +hypdif*vezm(i,k)
#  elif (VDEPT_HDIF>0)
    dif4x=hypdif*vexm(i,k)
    dif4z=hypdif*vezm(i,k)
#  elif (HDIF4==1)
    dif4x=dif4
    dif4z=dif4
#  elif (HDIF4==2)
    dif4x=dif4*(dif4_prof(i+1,k)+dif4_prof(i,k))*half
    dif4z=dif4*(dif4_prof(i,k+1)+dif4_prof(i,k))*half
#  endif               
    fxx(i,k)=fox(i+2,k)-fox(i-1,k)-3._RTYPE*(fox(i+1,k) &
              -fox(i+0,k))
    fzx(i,k)=fox(i,k+2)-fox(i,k-1)-3._RTYPE*(fox(i,k+1) &
              -fox(i,k+0))
    fxx(i,k)=fxx(i,k)*dif4x
    fzx(i,k)=fzx(i,k)*dif4z

    fxy(i,k)=foy(i+2,k)-foy(i-1,k)-3._RTYPE*(foy(i+1,k) &
              -foy(i+0,k))
    fzy(i,k)=foy(i,k+2)-foy(i,k-1)-3._RTYPE*(foy(i,k+1) &
              -foy(i,k+0))
    fxy(i,k)=fxy(i,k)*dif4x
    fzy(i,k)=fzy(i,k)*dif4z

    fxz(i,k)=foz(i+2,k)-foz(i-1,k)-3._RTYPE*(foz(i+1,k) &
              -foz(i+0,k))
    fzz(i,k)=foz(i,k+2)-foz(i,k-1)-3._RTYPE*(foz(i,k+1) &
              -foz(i,k+0))
    fxz(i,k)=fxz(i,k)*dif4x
    fzz(i,k)=fzz(i,k)*dif4z

#  if SIMPLE_LIMITER==1
    if (fxx(i,k)*(fox(i+1,k)-fox(i,k)).gt.zero) then
      if (fxx(i,k)*(fox(i+2,k)-fox(i+1,k)).lt.zero.or. &
        fxx(i,k)*(fox(i,k)-fox(i-1,k)).lt.zero) then
        fxx(i,k)=zero
      end if
    end if
    if (fzx(i,k)*(fox(i,k+1)-fox(i,k)).gt.zero) then
      if (fzx(i,k)*(fox(i,k+2)-fox(i,k+1)).lt.zero.or. &
          fzx(i,k)*(fox(i,k)-fox(i,k-1)).lt.zero) then
        fzx(i,k)=zero
      end if
    end if

    if (fxy(i,k)*(foy(i+1,k)-foy(i,k)).gt.zero) then
      if (fxy(i,k)*(foy(i+2,k)-foy(i+1,k)).lt.zero.or. &
        fxy(i,k)*(foy(i,k)-foy(i-1,k)).lt.zero) then
        fxy(i,k)=zero
      end if
    end if
    if (fzy(i,k)*(foy(i,k+1)-foy(i,k)).gt.zero) then
      if (fzy(i,k)*(foy(i,k+2)-foy(i,k+1)).lt.zero.or. &
          fzy(i,k)*(foy(i,k)-foy(i,k-1)).lt.zero) then
        fzy(i,k)=zero
      end if
    end if

    if (fxz(i,k)*(foz(i+1,k)-foz(i,k)).gt.zero) then
      if (fxz(i,k)*(foz(i+2,k)-foz(i+1,k)).lt.zero.or. &
        fxz(i,k)*(foz(i,k)-foz(i-1,k)).lt.zero) then
        fxz(i,k)=zero
      end if
    end if
    if (fzz(i,k)*(foz(i,k+1)-foz(i,k)).gt.zero) then
      if (fzz(i,k)*(foz(i,k+2)-foz(i,k+1)).lt.zero.or. &
          fzz(i,k)*(foz(i,k)-foz(i,k-1)).lt.zero) then
        fzz(i,k)=zero
      end if
    end if
    
#  elif SIMPLE_LIMITER==2
    if (fxx(i,k)*(fox(i+1,k)-fox(i,k)).gt.zero) then
      fxx(i,k)=zero
    end if
    if (fzx(i,k)*(fox(i,k+1)-fox(i,k)).gt.zero) then
      fzx(i,k)=zero
    end if

    if (fxy(i,k)*(foy(i+1,k)-foy(i,k)).gt.zero) then
      fxy(i,k)=zero
    end if
    if (fzy(i,k)*(foy(i,k+1)-foy(i,k)).gt.zero) then
      fzy(i,k)=zero
    end if

    if (fxz(i,k)*(foz(i+1,k)-foz(i,k)).gt.zero) then
      fxz(i,k)=zero
    end if
    if (fzz(i,k)*(foz(i,k+1)-foz(i,k)).gt.zero) then
      fzz(i,k)=zero
    end if
    
#  endif
   
  end do
  end do
  !$OMP END DO

! Although this part may be fused into the previous loop to 
! be more cache friendly, it will cause problem when OpenMP
! is used.
  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
    divfx(i,k)=divfx(i,k)-((fxx(i,k)-fxx(i-1,k))*dx24(i)  &
                          +(fzx(i,k)-fzx(i,k-1))*dz24(k))
    divfz(i,k)=divfz(i,k)-((fxz(i,k)-fxz(i-1,k))*dx24(i)  &
                          +(fzz(i,k)-fzz(i,k-1))*dz24(k))
#  if TWO_AND_HALF==1
    divfy(i,k)=divfy(i,k)-((fxy(i,k)-fxy(i-1,k))*dx24(i)  &
                          +(fzy(i,k)-fzy(i,k-1))*dz24(k))
#  endif
  end do
  end do
  !$OMP END DO
#endif


end subroutine hyperdif_vec_e_limited


!########################################################################
subroutine fin
! fill in ghost cells and calculate other fields


  call fin_primary
  call cal_b
  call fin_b
  call cal_aux_field
  call fin_aux
  call cal_e
#if TWO_AND_HALF==1 
  call fin_e
#endif  
#if VDEPT_HDIF>0
  call vdept_hdif_coef 
#endif       
end subroutine fin

!########################################################################
subroutine cal_b
  INTEGER :: i, k

#if PSI0==1
  !$OMP DO   
  do k=1,nz
  do i=1,nx
    psitot(i,k)=psii(i,k)+psi0(i,k)
  end do
  end do
  !$OMP END DO
#endif

  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
#if PSI0==1
    bxi(i,k)=-DDZ(psitot)
    bzi(i,k)=DDX(psitot) 
#else
    bxi(i,k)=-DDZ(psii)                 ! bxi=-dpsii/dz
    bzi(i,k)=DDX(psii)                  ! bzi=dpsii/dx
#endif
  end do
  end do
  !$OMP END DO
end subroutine cal_b

!########################################################################
subroutine cal_aux_field
! Calculate Auxiliary fields
  IMPLICIT NONE
  INTEGER :: i, k
  REAL(RTYPE) :: ideni, idt
!=========================================
  idt=-0.25_RTYPE/dt/omgamma
  !$OMP DO 
  do k=1,nz
  do i=1,nx
    ideni=one/deni(i,k)
#if TWO_AND_HALF==1
    vy(i,k)=puyi(i,k)*ideni
#endif
    vx(i,k)=puxi(i,k)*ideni    
    vz(i,k)=puzi(i,k)*ideni     

#if VISCOSITY==1
    vden(i,k)=visc*deni(i,k)
#endif


#if TWOFLUID==1
#  if ISOTHERMAL==0 && SPE==0
    pei(i,k)=prei(i,k)
    Ti(i,k)=prei(i,k)/deni(i,k)  
    Te(i,k)=Ti(i,k) 
#  elif ISOTHERMAL==0 && SPE==1
    Ti(i,k)=prei(i,k)/deni(i,k)  
    Te(i,k)=pei(i,k)/deni(i,k) 
#  elif ISOTHERMAL==1 
    prei(i,k)=T0*deni(i,k)
    pei(i,k)=prei(i,k)
    Ti(i,k)=T0
    Te(i,k)=T0
#  elif ISOTHERMAL==2
    prei(i,k)=T0*deni(i,k)
    pei(i,k)=Te0*deni(i,k)
    Ti(i,k)=T0
    Te(i,k)=Te0
#endif
    pt(i,k)=prei(i,k)+pei(i,k) 
#else
#  if ISOTHERMAL>0
    prei(i,k)=T0*deni(i,k)
    Ti(i,k)=T0
#  else
    Ti(i,k)=prei(i,k)/deni(i,k)  
#  endif
    pt(i,k)=two*prei(i,k)
#endif

      
#if RESISTIVITY==2 || RESISTIVITY==3
#  if ISOTHERMAL==0
#    if (TWOFLUID!=0 && SPE!=0)
    eta_perp(i,k)=eta*((T0/Te(i,k))**1.5_RTYPE)
#    else  
    eta_perp(i,k)=eta*((T0/Ti(i,k))**1.5_RTYPE)
#    endif
#  else
    eta_perp(i,k)=eta
#  endif          
#endif

#if ISOTHERMAL==0 && THERMAL_EXCHANGE==1 && TWOFLUID==1 && SPE==1 
#  if RESISTIVITY==1
    heat_exchange(i,k)=  &
      min(idt,(3._RTYPE*eta/(di*di))*deni(i,k))  &
        *(prei(i,k)-pei(i,k))
#  elif RESISTIVITY==2
    heat_exchange(i,k)=  &
      min(idt,(3._RTYPE*eta_perp(i,k)/(di*di))*deni(i,k))  &
        *(prei(i,k)-pei(i,k))
#  endif                      
#endif


#if (JXB==0)
#  if ((TWOFLUID==1 && HALL_CAL==1) || AMBIPOLAR==1)    
    ! Magnetic Stress Tensor
    fxx(i,k)=( bxi(i,k)*bxi(i,k) &
              -bzi(i,k)*bzi(i,k))*half
    fxz(i,k)=bxi(i,k)*bzi(i,k)
    fzz(i,k)=(-bxi(i,k)*bxi(i,k) &
              +bzi(i,k)*bzi(i,k))*half
#    if TWO_AND_HALF==1
    fxx(i,k)=fxx(i,k)+(-byi(i,k)*byi(i,k))*half
    fzz(i,k)=fzz(i,k)+(-byi(i,k)*byi(i,k))*half
    fxy(i,k)=bxi(i,k)*byi(i,k)
    fyz(i,k)=byi(i,k)*bzi(i,k)
#    endif
#  endif
#endif
 
  end do
  end do
  !$OMP END DO


  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
#if J_CAL==0    
! use 5-pt stencil to calculate J
#  if JY_CAL==0    
    jyi(i,k)=-(D2DX(psii)+D2DZ(psii))  ! jyi=-del2 psii
#  elif JY_CAL==1
    jyi(i,k)=DDZ(bxi)-DDX(bzi)         ! jyi=d Bx/dz -d Bz/dx
#  endif
#  if TWO_AND_HALF==1
    jxi(i,k)=-DDZ(byi)
    jzi(i,k)=DDX(byi)
#  endif
#elif J_CAL==1
! use 3-pt stencil to calculate J
#  if JY_CAL==0    
    jyi(i,k)=-(D2DX3(psii)+D2DZ3(psii))  ! jyi=-del2 psii
#  elif JY_CAL==1
    jyi(i,k)=DDZ3(bxi)-DDX3(bzi)         ! jyi=d Bx/dz -d Bz/dx
#  endif
#  if TWO_AND_HALF==1
    jxi(i,k)=-DDZ3(byi)
    jzi(i,k)=DDX3(byi)
#  endif
#endif
! pre-calculate magnetic force
#if (JXB==0)
#  if ((TWOFLUID==1 && HALL_CAL==1) || AMBIPOLAR==1) 
    fbx(i,k)=DDX(fxx)+DDZ(fxz)
    fbz(i,k)=DDX(fxz)+DDZ(fzz)
#    if TWO_AND_HALF==1    
    fby(i,k)=DDX(fxy)+DDZ(fyz)
#    endif    
#  endif    
#elif (JXB==1)
#  if TWO_AND_HALF==1
    fbx(i,k)=jyi(i,k)*bzi(i,k)-jzi(i,k)*byi(i,k)
    fby(i,k)=jzi(i,k)*bxi(i,k)-jxi(i,k)*bzi(i,k)
    fbz(i,k)=jxi(i,k)*byi(i,k)-jyi(i,k)*bxi(i,k)
#  else
    fbx(i,k)=jyi(i,k)*bzi(i,k)
    fbz(i,k)=-jyi(i,k)*bxi(i,k)
#  endif
#endif

  end do
  end do
  !$OMP END DO


end subroutine cal_aux_field
!########################################################################
subroutine cal_e
! Calculate E fields
  IMPLICIT NONE
  INTEGER :: i, k
  REAL(RTYPE) :: bb, bbx, bby, bbz, xi_e, F1, F2, delta
  REAL(RTYPE) :: c1, c2, c3, c1jpar
  REAL(RTYPE) :: e1, e2, e3
  REAL(RTYPE) :: dideni
#if TWOFLUID==1
  !$OMP DO 
  do k=1,nz
  do i=1,nx
    dideni=di/deni(i,k)
    vex(i,k)=vx(i,k)-dideni*jxi(i,k)
    vey(i,k)=vy(i,k)-dideni*jyi(i,k)
    vez(i,k)=vz(i,k)-dideni*jzi(i,k)
  end do
  end do
  !$OMP END DO
#endif

  !$OMP DO
  do k=3,nz-2
  do i=3,nx-2
#if TWOFLUID==0
!============= Ideal MHD part ================================
    ey(i,k)=vx(i,k)*bzi(i,k)-vz(i,k)*bxi(i,k)
#  if TWO_AND_HALF==1
    ! Ex and Ez only used in 3D or 2.5 D
    ez(i,k)=-vx(i,k)*byi(i,k)+bxi(i,k)*vy(i,k)
    ex(i,k)=+vz(i,k)*byi(i,k)-bzi(i,k)*vy(i,k)
#  endif
#elif TWOFLUID==1
!=========== Ideal + two-fluid Hall terms. ===============
! Note that two-fluid runs must be 2.5D
#  if HALL_CAL==0
! Calculate Hall term with JxB, using  -Ve x B
#    if ISOTHERMAL==0
    dideni=di/deni(i,k)
    ex(i,k)=vez(i,k)*byi(i,k)-vey(i,k)*bzi(i,k) &
            +dideni*(-DDX(pei))
    ey(i,k)=vex(i,k)*bzi(i,k)-vez(i,k)*bxi(i,k) 
    ez(i,k)=vey(i,k)*bxi(i,k)-vex(i,k)*byi(i,k) &
            +dideni*(-DDZ(pei))
#    else
! Isothermal EOS, ignore grad. pe term
! Note that, in this case ez and ex are not exactly the electric field,
! but the ignored pressure gradient term contributes nothing after taking 
! the curl
    ex(i,k)=vez(i,k)*byi(i,k)-vey(i,k)*bzi(i,k) 
    ey(i,k)=vex(i,k)*bzi(i,k)-vez(i,k)*bxi(i,k) 
    ez(i,k)=vey(i,k)*bxi(i,k)-vex(i,k)*byi(i,k) 
#    endif    
#  elif HALL_CAL==1
! Calculate Hall term with magnetic force (precalculated in CAL_AUX)
    dideni=di/deni(i,k)
#    if ISOTHERMAL==0
    ex(i,k)=vz(i,k)*byi(i,k)-vy(i,k)*bzi(i,k) &
            +dideni*(fbx(i,k)-DDX(pei))
    ey(i,k)=vx(i,k)*bzi(i,k)-vz(i,k)*bxi(i,k) &
            +dideni*(fby(i,k))
    ez(i,k)=vy(i,k)*bxi(i,k)-vx(i,k)*byi(i,k) &
              +dideni*(fbz(i,k)-DDZ(pei))
#    else
! Isothermal EOS, ignore grad. pe term
! Note that, in this case ez and ex are not exactly the electric field,
! but the ignored pressure gradient term contributes nothing after taking 
! the curl
    ex(i,k)=vz(i,k)*byi(i,k)-vy(i,k)*bzi(i,k) &
              +dideni*fbx(i,k)
    ey(i,k)=vx(i,k)*bzi(i,k)-vz(i,k)*bxi(i,k) &
              +dideni*fby(i,k)
    ez(i,k)=vy(i,k)*bxi(i,k)-vx(i,k)*byi(i,k) &
              +dideni*fbz(i,k)
#    endif            
#  endif 
#endif

#if AMBIPOLAR==1
!=================== Add Ambipolar diffusion term ============
    ey(i,k)=ey(i,k)-eta_ad*( fbz(i,k)*bxi(i,k) &
                            -fbx(i,k)*bzi(i,k))
#  if TWO_AND_HALF==1
    ex(i,k)=ex(i,k)-eta_ad*( fby(i,k)*bzi(i,k) &
                            -fbz(i,k)*byi(i,k))
    ez(i,k)=ez(i,k)-eta_ad*( fbx(i,k)*byi(i,k) &
                            -fby(i,k)*bxi(i,k))
#  endif 
#endif


#if RESISTIVITY==1 
#  if B_DIFF==0
!============ Add resistivity ====================================
    ey(i,k)=ey(i,k)+eta*jyi(i,k)
#    if TWO_AND_HALF==1 
! ex and ez are only used in sby
    ez(i,k)=ez(i,k)+eta*jzi(i,k)
    ex(i,k)=ex(i,k)+eta*jxi(i,k)
#    endif
#  endif
#  if OHMIC_HEATING==1
    electron_heating(i,k)=eta*( jxi(i,k)**2  &
                               +jyi(i,k)**2  &
                               +jzi(i,k)**2)

#  endif
#elif RESISTIVITY==2
!====== Isotropic but temperature dependent resistivity =======
    ey(i,k)=ey(i,k)+eta_perp(i,k)*jyi(i,k)
#  if TWO_AND_HALF==1 
! ex and ez are only used in sby
    ez(i,k)=ez(i,k)+eta_perp(i,k)*jzi(i,k)
    ex(i,k)=ex(i,k)+eta_perp(i,k)*jxi(i,k)
#  endif
#  if OHMIC_HEATING==1
    electron_heating(i,k)=eta_perp(i,k)*( jxi(i,k)**2  &
                                         +jyi(i,k)**2  &
                                         +jzi(i,k)**2)
#  endif
#elif RESISTIVITY==3
!====== Braginskii parallel and perpendicular resistivity ======       
    bb=sqrt(bxi(i,k)**2+byi(i,k)**2+bzi(i,k)**2)
    xi_e=bb*di/(eta_perp(i,k)*deni(i,k))
    delta=xi_e**4+14.79*xi_e**2+3.77
    F1=one-(6.42*xi_e**2+1.84)/delta
    c2=F1*eta_perp(i,k)
    if (bb.gt.zero) then
      F2=xi_e*(1.7*xi_e**2+0.78)/delta
      c1=(1.93/3.77-F1)*eta_perp(i,k)
      c3=F2*eta_perp(i,k) 
      bbx=bxi(i,k)/bb
      bby=byi(i,k)/bb
      bbz=bzi(i,k)/bb
      c1jpar=c1*(jxi(i,k)*bbx+jyi(i,k)*bby+jzi(i,k)*bbz)
      e1=c1jpar*bbx+c2*jxi(i,k)+c3*(jyi(i,k)*bbz-jzi(i,k)*bby)
      e2=c1jpar*bby+c2*jyi(i,k)+c3*(jzi(i,k)*bbx-jxi(i,k)*bbz)
      e3=c1jpar*bbz+c2*jzi(i,k)+c3*(jxi(i,k)*bby-jyi(i,k)*bbx) 
    else
      e1=c2*jxi(i,k)
      e2=c2*jyi(i,k)
      e3=c2*jzi(i,k)
    end if               
    ey(i,k)=ey(i,k)+e2
#  if TWO_AND_HALF==1
! ex and ez are only used in sby
    ez(i,k)=ez(i,k)+e3
    ex(i,k)=ex(i,k)+e1
#  endif
#  if OHMIC_HEATING==1
    electron_heating(i,k)=(e1*jxi(i,k)+e2*jyi(i,k)+e3*jzi(i,k))
#  endif
#endif

#if HYPER_RESISTIVITY==1
! Add hyper-resistivity
    ey(i,k)=ey(i,k)-etah*(D2DX(jyi)+D2DZ(jyi))
#  if TWO_AND_HALF==1 
    ex(i,k)=ex(i,k)-etah*(D2DX(jxi)+D2DZ(jxi))
    ez(i,k)=ez(i,k)-etah*(D2DX(jzi)+D2DZ(jzi))
#  endif
#endif
  end do
  end do
  !$OMP END DO

end subroutine cal_e
!########################################################################
subroutine vdept_hdif_coef
! Calculate velocity dependent hyperdiffusion coefficients
  IMPLICIT NONE
  INTEGER :: i, k

! average the flow in mid points.  This is to be used in velocity dependent
! hyper diffusion

#if (VDEPT_HDIF==2)
! Sound speed
  !$OMP DO 
  do k=2,nz-1
    do i=2,nx-1
      temp(i,k)=sqrt(gamma*pt(i,k)/deni(i,k))
    end do
  end do
  !$OMP END DO
#elif (VDEPT_HDIF==3)
! Alfven speed
  !$OMP DO 
  do k=2,nz-1
    do i=2,nx-1
      temp(i,k)=sqrt((bxi(i,k)**2+byi(i,k)**2+bzi(i,k)**2)/deni(i,k))
    end do
  end do
  !$OMP END DO
#elif (VDEPT_HDIF==4)
! fast magnetosonic wave speed
  !$OMP DO 
  do k=2,nz-1
    do i=2,nx-1
      temp(i,k)=sqrt((bxi(i,k)**2+byi(i,k)**2+bzi(i,k)**2 &
               +gamma*pt(i,k))/deni(i,k))   
    end do
  end do
  !$OMP END DO
#endif

  !$OMP DO 
  do k=2,nz-2
    do i=2,nx-2
#if (VDEPT_HDIF==1)
      vxm(i,k)=abs(vx(i+1,k))+abs(vx(i,k))
      vzm(i,k)=abs(vz(i,k+1))+abs(vz(i,k))
#elif (VDEPT_HDIF>1)
      vxm(i,k)=abs(vx(i+1,k))+abs(vx(i,k))+temp(i+1,k)+temp(i,k)
      vzm(i,k)=abs(vz(i,k+1))+abs(vz(i,k))+temp(i,k+1)+temp(i,k)
#endif
#if (TWOFLUID==1 && VDEPT_HDIF==1)
      vexm(i,k)=abs(vex(i+1,k))+abs(vex(i,k))
      vezm(i,k)=abs(vez(i,k+1))+abs(vez(i,k))
#elif (TWOFLUID==1 && VDEPT_HDIF>1)
      vexm(i,k)=abs(vex(i+1,k))+abs(vex(i,k))+temp(i+1,k)+temp(i,k)
      vezm(i,k)=abs(vez(i,k+1))+abs(vez(i,k))+temp(i,k+1)+temp(i,k)
#endif
    end do
  end do
  !$OMP END DO 
end subroutine vdept_hdif_coef 
!=========================================================================
end module equation
