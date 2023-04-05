














































 






 










 


















 








 























 

  










 
 






 
























 

 

 









 



























!############################################################
! Definitions for the m4 preprocessor
!############################################################

! Undefine intrinsic macros that conflict with f90

   ! include statements
       ! len defines length of character variables
     ! dummy argument in sub set_fname


! d/dx

!d^2/dx^2



! d/dz

! d^2/dz^2


! 3-point stencil finite difference

! d/dx

! d^2/dx^2




! d/dz

! d^2/dz^2



      
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
  by=byi  
  !$OMP END WORKSHARE      

  f=1
  
  call smom2(f)
  call spre2(f)
  call sden2(f)
  call spsi2(f)
  call sby2(f)

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

! Trapezoidal Leapfrog 
  ns=1
  time_i(0)=time+half*dt
  time_i(1)=time+dt




  
  do f=0, ns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call smom2(f)
    call spre2(f)
    call sden2(f)
    call spsi2(f)
    call sby2(f)
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
    divf(i,k)=(fx(i-2,k)*px51u(1)+fx(i-1,k)*px51u(2) &
                +fx(i+1,k)*px51u(4) &
                +fx(i+2,k)*px51u(5))+(fz(i,k-2)*pz51u(1)+fz(i,k-1)*pz51u(2) &
                +fz(i,k+1)*pz51u(4) &
                +fz(i,k+2)*pz51u(5))
! diffusion. Beware, this is not physical. Turn it off if you can.
  end do
  end do
  !$OMP END DO      

!hyper diffusion
  call hyperdif(deni,hdifden,dif4den)

  call step2d(deni,den,divf,f)

  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
    if (deni(i,k).le.zero.or.deni(i,k).gt.den_max) then
      ifabort=.true.
      print*, 'x=',x(i),'  z=',z(k),'  den=',deni(i,k)
    end if
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



  !$OMP DO 
  do k=3, nz-2
  do i=3, nx-2
    divf(i,k)=(fx(i-2,k)*px51u(1)+fx(i-1,k)*px51u(2) &
                +fx(i+1,k)*px51u(4) &
                +fx(i+2,k)*px51u(5))+(fz(i,k-2)*pz51u(1)+fz(i,k-1)*pz51u(2) &
                +fz(i,k+1)*pz51u(4) &
                +fz(i,k+2)*pz51u(5))

!===============pre====================================
    divf(i,k)=divf(i,k)    &
              +omgamma*prei(i,k)*((vx(i-2,k)*px51u(1)+vx(i-1,k)*px51u(2) &
                +vx(i+1,k)*px51u(4) &
                +vx(i+2,k)*px51u(5))+(vz(i,k-2)*pz51u(1)+vz(i,k-1)*pz51u(2) &
                +vz(i,k+1)*pz51u(4) &
                +vz(i,k+2)*pz51u(5)))
!======================================================
! Ohmic heating, the "half" comes from the fact that heat goes to electron
! and ion equally.  If PE is stepped, then heat goes to electron.

! Viscous heating 




! Thermal exchange between species

! Thermal Conductivity


  end do
  end do
  !$OMP END DO      

!hyper diffusion
  call hyperdif(prei,hdifpre,dif4pre)
!

  call step2d(prei,pre,divf,f)
!
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
    divf(i,k)=(fx(i-2,k)*px51u(1)+fx(i-1,k)*px51u(2) &
                +fx(i+1,k)*px51u(4) &
                +fx(i+2,k)*px51u(5))+(fz(i,k-2)*pz51u(1)+fz(i,k-1)*pz51u(2) &
                +fz(i,k+1)*pz51u(4) &
                +fz(i,k+2)*pz51u(5))

!===============pe====================================
    divf(i,k)=divf(i,k)    &
              +omgamma*pei(i,k)*((vex(i-2,k)*px51u(1)+vex(i-1,k)*px51u(2) &
                +vex(i+1,k)*px51u(4) &
                +vex(i+2,k)*px51u(5))+(vez(i,k-2)*pz51u(1)+vez(i,k-1)*pz51u(2) &
                +vez(i,k+1)*pz51u(4) &
                +vez(i,k+2)*pz51u(5)))
!======================================================
  end do
  end do
  !$OMP END DO

!hyper diffusion
  call hyperdif_e(pei,hdifpe,dif4pe)
!

  call step2d(pei,pe,divf,f)
!


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
    fxy(i,k)=-vy(i,k)*puxi(i,k)
    fyz(i,k)=-vz(i,k)*puyi(i,k)
 
    ! Add Magnetic Stress Tensor if magnetic force is not pre-calculated
    fxx(i,k)=fxx(i,k)+( bxi(i,k)*bxi(i,k) &
                       -bzi(i,k)*bzi(i,k))*half
    fxz(i,k)=fxz(i,k)+bxi(i,k)*bzi(i,k)
    fzz(i,k)=fzz(i,k)+(-bxi(i,k)*bxi(i,k) &
                       +bzi(i,k)*bzi(i,k))*half
    fxx(i,k)=fxx(i,k)+(-byi(i,k)*byi(i,k))*half
    fzz(i,k)=fzz(i,k)+(-byi(i,k)*byi(i,k))*half
    fxy(i,k)=fxy(i,k)+bxi(i,k)*byi(i,k)
    fyz(i,k)=fyz(i,k)+byi(i,k)*bzi(i,k)
  end do
  end do
  !$OMP END DO NOWAIT

! Viscous Stress Tensor
  !$OMP DO 
  do k=2,nz-1
  do i=2,nx-1
    vxx(i,k)=(vx(i-1,k)*px31u(1) &
                +vx(i+1,k)*px31u(3))
    vxz(i,k)=((vz(i-1,k)*px31u(1) &
                +vz(i+1,k)*px31u(3))+(vx(i,k-1)*pz31u(1) &
                +vx(i,k+1)*pz31u(3)))*half
    vzz(i,k)=(vz(i,k-1)*pz31u(1) &
                +vz(i,k+1)*pz31u(3))
    vxy(i,k)=(vy(i-1,k)*px31u(1) &
                +vy(i+1,k)*px31u(3))*half
    vyz(i,k)=(vy(i,k-1)*pz31u(1) &
                +vy(i,k+1)*pz31u(3))*half
    
  end do
  end do
  !$OMP END DO NOWAIT


  !$OMP BARRIER

  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
    divfx(i,k)=(fxx(i-2,k)*px51u(1)+fxx(i-1,k)*px51u(2) &
                +fxx(i+1,k)*px51u(4) &
                +fxx(i+2,k)*px51u(5))+(fxz(i,k-2)*pz51u(1)+fxz(i,k-1)*pz51u(2) &
                +fxz(i,k+1)*pz51u(4) &
                +fxz(i,k+2)*pz51u(5))
    divfz(i,k)=(fxz(i-2,k)*px51u(1)+fxz(i-1,k)*px51u(2) &
                +fxz(i+1,k)*px51u(4) &
                +fxz(i+2,k)*px51u(5))+(fzz(i,k-2)*pz51u(1)+fzz(i,k-1)*pz51u(2) &
                +fzz(i,k+1)*pz51u(4) &
                +fzz(i,k+2)*pz51u(5))
    divfy(i,k)=(fxy(i-2,k)*px51u(1)+fxy(i-1,k)*px51u(2) &
                +fxy(i+1,k)*px51u(4) &
                +fxy(i+2,k)*px51u(5))+(fyz(i,k-2)*pz51u(1)+fyz(i,k-1)*pz51u(2) &
                +fyz(i,k+1)*pz51u(4) &
                +fyz(i,k+2)*pz51u(5))

!===============  Momentum ====================================
!======================================================
! Viscosity & friction
    dxvden=(vden(i-1,k)*px31u(1) &
                +vden(i+1,k)*px31u(3))
    dzvden=(vden(i,k-1)*pz31u(1) &
                +vden(i,k+1)*pz31u(3))
    divfx(i,k)=divfx(i,k) &
               +vden(i,k)*( (vx(i-1,k)*px32u(1) &
                 +vx(i  ,k)*px32u(2) &
                 +vx(i+1,k)*px32u(3)) &
                           +half*( (vx(i,k-1)*pz32u(1) &
                 +vx(i,k  )*pz32u(2) &
                 +vx(i,k+1)*pz32u(3)) &
                                  +(vzz(i-1,k)*px31u(1) &
                +vzz(i+1,k)*px31u(3)))) &
               +vxx(i,k)*dxvden &
               +vxz(i,k)*dzvden
    divfz(i,k)=divfz(i,k) &
               +vden(i,k)*( (vz(i,k-1)*pz32u(1) &
                 +vz(i,k  )*pz32u(2) &
                 +vz(i,k+1)*pz32u(3)) &
                           +half*( (vz(i-1,k)*px32u(1) &
                 +vz(i  ,k)*px32u(2) &
                 +vz(i+1,k)*px32u(3)) &
                                  +(vxx(i,k-1)*pz31u(1) &
                +vxx(i,k+1)*pz31u(3)))) &
               +vxz(i,k)*dxvden &
               +vzz(i,k)*dzvden
    divfy(i,k)=divfy(i,k) &
               +vden(i,k)*half*((vy(i-1,k)*px32u(1) &
                 +vy(i  ,k)*px32u(2) &
                 +vy(i+1,k)*px32u(3))+(vy(i,k-1)*pz32u(1) &
                 +vy(i,k  )*pz32u(2) &
                 +vy(i,k+1)*pz32u(3))) &
               +vxy(i,k)*dxvden &
               +vyz(i,k)*dzvden

! Gravity

  end do
  end do
  !$OMP END DO


!hyper diffusion
  call hyperdif_vec(puxi,puyi,puzi,hdifmom,dif4mom)

!

  call step2d(puxi,pux,divfx,f)
  call step2d(puzi,puz,divfz,f)
  call step2d(puyi,puy,divfy,f)
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
    ! Add Magnetic Stress Tensor if magnetic force is not pre-calculated
    fx(i,k)=fx(i,k)+( bxi(i,k)*bxi(i,k)-byi(i,k)*byi(i,k) &
                       -bzi(i,k)*bzi(i,k))*half
    fz(i,k)=fz(i,k)+bxi(i,k)*bzi(i,k)
    
  end do
  end do
  !$OMP END DO      


  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
    divf(i,k)=(fx(i-2,k)*px51u(1)+fx(i-1,k)*px51u(2) &
                +fx(i+1,k)*px51u(4) &
                +fx(i+2,k)*px51u(5))+(fz(i,k-2)*pz51u(1)+fz(i,k-1)*pz51u(2) &
                +fz(i,k+1)*pz51u(4) &
                +fz(i,k+2)*pz51u(5))

!===============pux====================================

!======================================================
! Viscosity & friction
    divf(i,k)=divf(i,k)+visc*((puxi(i-2,k)*px52u(1)+puxi(i-1,k)*px52u(2) &
                 +puxi(i  ,k)*px52u(3)+puxi(i+1,k)*px52u(4) &
                 +puxi(i+2,k)*px52u(5))+(puxi(i,k-2)*pz52u(1)+puxi(i,k-1)*pz52u(2) &
                 +puxi(i,k  )*pz52u(3)+puxi(i,k+1)*pz52u(4) &
                 +puxi(i,k+2)*pz52u(5))) 
  end do
  end do
  !$OMP END DO

!hyper diffusion
  call hyperdif(puxi,hdifmom,dif4mom)

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
    ! Add Magnetic Stress Tensor if magnetic force is not pre-calculated
    fx(i,k)=fx(i,k)+bxi(i,k)*byi(i,k)
    fz(i,k)=fz(i,k)+byi(i,k)*bzi(i,k)
    
  end do
  end do
  !$OMP END DO

  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
    divf(i,k)=(fx(i-2,k)*px51u(1)+fx(i-1,k)*px51u(2) &
                +fx(i+1,k)*px51u(4) &
                +fx(i+2,k)*px51u(5))+(fz(i,k-2)*pz51u(1)+fz(i,k-1)*pz51u(2) &
                +fz(i,k+1)*pz51u(4) &
                +fz(i,k+2)*pz51u(5))
!==================puy=========================================
!==============================================================
! Viscosity & friction
    divf(i,k)=divf(i,k)+visc*((puyi(i-2,k)*px52u(1)+puyi(i-1,k)*px52u(2) &
                 +puyi(i  ,k)*px52u(3)+puyi(i+1,k)*px52u(4) &
                 +puyi(i+2,k)*px52u(5))+(puyi(i,k-2)*pz52u(1)+puyi(i,k-1)*pz52u(2) &
                 +puyi(i,k  )*pz52u(3)+puyi(i,k+1)*pz52u(4) &
                 +puyi(i,k+2)*pz52u(5))) 
  end do
  end do
!$OMP END DO

!hyper diffusion
  call hyperdif(puyi,hdifmom,dif4mom)


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
    ! Add Magnetic Stress Tensor if magnetic force is not pre-calculated
    fx(i,k)=fx(i,k)+bxi(i,k)*bzi(i,k)
    fz(i,k)=fz(i,k)+(-bxi(i,k)*bxi(i,k)-byi(i,k)*byi(i,k) &
                     +bzi(i,k)*bzi(i,k))*half
    
  end do
  end do
!$OMP END DO


  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
    divf(i,k)=(fx(i-2,k)*px51u(1)+fx(i-1,k)*px51u(2) &
                +fx(i+1,k)*px51u(4) &
                +fx(i+2,k)*px51u(5))+(fz(i,k-2)*pz51u(1)+fz(i,k-1)*pz51u(2) &
                +fz(i,k+1)*pz51u(4) &
                +fz(i,k+2)*pz51u(5))

!===============puz====================================
!======================================================
! Viscosity & friction
    divf(i,k)=divf(i,k)+visc*((puzi(i-2,k)*px52u(1)+puzi(i-1,k)*px52u(2) &
                 +puzi(i  ,k)*px52u(3)+puzi(i+1,k)*px52u(4) &
                 +puzi(i+2,k)*px52u(5))+(puzi(i,k-2)*pz52u(1)+puzi(i,k-1)*pz52u(2) &
                 +puzi(i,k  )*pz52u(3)+puzi(i,k+1)*pz52u(4) &
                 +puzi(i,k+2)*pz52u(5))) 
  end do
  end do
  !$OMP END DO

!hyper diffusion
  call hyperdif(puzi,hdifmom,dif4mom)

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
    divf(i,k)=(ez(i-2,k)*px51u(1)+ez(i-1,k)*px51u(2) &
                +ez(i+1,k)*px51u(4) &
                +ez(i+2,k)*px51u(5))-(ex(i,k-2)*pz51u(1)+ex(i,k-1)*pz51u(2) &
                +ex(i,k+1)*pz51u(4) &
                +ex(i,k+2)*pz51u(5))
!============================================================

  end do
  end do
  !$OMP END DO

!hyper diffusion
  call hyperdif(byi,hdifb,dif4b)

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
  end do
  end do
  !$OMP END DO

!hyper diffusion
  call hyperdif(psii,hdifb,dif4b)

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
      dif4x=hypdif*vxm(i,k)
      dif4z=hypdif*vzm(i,k)
      fx(i,k)=foo(i+2,k)-foo(i-1,k)-3._RTYPE*(foo(i+1,k)-foo(i+0,k))
      fz(i,k)=foo(i,k+2)-foo(i,k-1)-3._RTYPE*(foo(i,k+1)-foo(i,k+0))
      fx(i,k)=fx(i,k)*dif4x
      fz(i,k)=fz(i,k)*dif4z
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
      dif4x=hypdif*vexm(i,k)
      dif4z=hypdif*vezm(i,k)
      fx(i,k)=foo(i+2,k)-foo(i-1,k)-3._RTYPE*(foo(i+1,k)-foo(i+0,k))
      fz(i,k)=foo(i,k+2)-foo(i,k-1)-3._RTYPE*(foo(i,k+1)-foo(i,k+0))
      fx(i,k)=fx(i,k)*dif4x
      fz(i,k)=fz(i,k)*dif4z
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
! Simple limiter is applied when 0>0 

  !$OMP DO 
  do k=2,nz-2
    do i=2,nx-2
      dif4x=hypdif*vxm(i,k)
      dif4z=hypdif*vzm(i,k)
      fx(i,k)=foo(i+2,k)-foo(i-1,k)-3._RTYPE*(foo(i+1,k)-foo(i+0,k))
      fz(i,k)=foo(i,k+2)-foo(i,k-1)-3._RTYPE*(foo(i,k+1)-foo(i,k+0))
      fx(i,k)=fx(i,k)*dif4x
      fz(i,k)=fz(i,k)*dif4z

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
! Simple limiter is applied when 0>0

  !$OMP DO 
  do k=2,nz-2
    do i=2,nx-2
      dif4x=hypdif*vexm(i,k)
      dif4z=hypdif*vezm(i,k)
      fx(i,k)=foo(i+2,k)-foo(i-1,k)-3._RTYPE*(foo(i+1,k)-foo(i+0,k))
      fz(i,k)=foo(i,k+2)-foo(i,k-1)-3._RTYPE*(foo(i,k+1)-foo(i,k+0))
      fx(i,k)=fx(i,k)*dif4x
      fz(i,k)=fz(i,k)*dif4z

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


  !$OMP DO 
  do k=2,nz-2
  do i=2,nx-2
    dif4x=hypdif*vxm(i,k)
    dif4z=hypdif*vzm(i,k)
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
    divfy(i,k)=divfy(i,k)-((fxy(i,k)-fxy(i-1,k))*dx24(i)  &
                          +(fzy(i,k)-fzy(i,k-1))*dz24(k))
  end do
  end do
  !$OMP END DO
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


  !$OMP DO 
  do k=2,nz-2
  do i=2,nx-2
    dif4x=hypdif*vexm(i,k)
    dif4z=hypdif*vezm(i,k)
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
    divfy(i,k)=divfy(i,k)-((fxy(i,k)-fxy(i-1,k))*dx24(i)  &
                          +(fzy(i,k)-fzy(i,k-1))*dz24(k))
  end do
  end do
  !$OMP END DO

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
! Simple limiter is applied when 0>0

  !$OMP DO 
  do k=2,nz-2
  do i=2,nx-2
    dif4x=hypdif*vxm(i,k)
    dif4z=hypdif*vzm(i,k)
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
    divfy(i,k)=divfy(i,k)-((fxy(i,k)-fxy(i-1,k))*dx24(i)  &
                          +(fzy(i,k)-fzy(i,k-1))*dz24(k))
  end do
  end do
  !$OMP END DO

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
! Simple limiter is applied when 0>0

  !$OMP DO 
  do k=2,nz-2
  do i=2,nx-2
    dif4x=hypdif*vexm(i,k)
    dif4z=hypdif*vezm(i,k)
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
    divfy(i,k)=divfy(i,k)-((fxy(i,k)-fxy(i-1,k))*dx24(i)  &
                          +(fzy(i,k)-fzy(i,k-1))*dz24(k))
  end do
  end do
  !$OMP END DO


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
  call fin_e
  call vdept_hdif_coef 
end subroutine fin

!########################################################################
subroutine cal_b
  INTEGER :: i, k


  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
    bxi(i,k)=-(psii(i,k-2)*pz51u(1)+psii(i,k-1)*pz51u(2) &
                +psii(i,k+1)*pz51u(4) &
                +psii(i,k+2)*pz51u(5))                 ! bxi=-dpsii/dz
    bzi(i,k)=(psii(i-2,k)*px51u(1)+psii(i-1,k)*px51u(2) &
                +psii(i+1,k)*px51u(4) &
                +psii(i+2,k)*px51u(5))                  ! bzi=dpsii/dx
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
    vy(i,k)=puyi(i,k)*ideni
    vx(i,k)=puxi(i,k)*ideni    
    vz(i,k)=puzi(i,k)*ideni     

    vden(i,k)=visc*deni(i,k)


    Ti(i,k)=prei(i,k)/deni(i,k)  
    pt(i,k)=two*prei(i,k)

      



 
  end do
  end do
  !$OMP END DO


  !$OMP DO 
  do k=3,nz-2
  do i=3,nx-2
! use 5-pt stencil to calculate J
    jyi(i,k)=-((psii(i-2,k)*px52u(1)+psii(i-1,k)*px52u(2) &
                 +psii(i  ,k)*px52u(3)+psii(i+1,k)*px52u(4) &
                 +psii(i+2,k)*px52u(5))+(psii(i,k-2)*pz52u(1)+psii(i,k-1)*pz52u(2) &
                 +psii(i,k  )*pz52u(3)+psii(i,k+1)*pz52u(4) &
                 +psii(i,k+2)*pz52u(5)))  ! jyi=-del2 psii
    jxi(i,k)=-(byi(i,k-2)*pz51u(1)+byi(i,k-1)*pz51u(2) &
                +byi(i,k+1)*pz51u(4) &
                +byi(i,k+2)*pz51u(5))
    jzi(i,k)=(byi(i-2,k)*px51u(1)+byi(i-1,k)*px51u(2) &
                +byi(i+1,k)*px51u(4) &
                +byi(i+2,k)*px51u(5))
! pre-calculate magnetic force

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

  !$OMP DO
  do k=3,nz-2
  do i=3,nx-2
!============= Ideal MHD part ================================
    ey(i,k)=vx(i,k)*bzi(i,k)-vz(i,k)*bxi(i,k)
    ! Ex and Ez only used in 3D or 2.5 D
    ez(i,k)=-vx(i,k)*byi(i,k)+bxi(i,k)*vy(i,k)
    ex(i,k)=+vz(i,k)*byi(i,k)-bzi(i,k)*vy(i,k)



!============ Add resistivity ====================================
    ey(i,k)=ey(i,k)+eta*jyi(i,k)
! ex and ez are only used in sby
    ez(i,k)=ez(i,k)+eta*jzi(i,k)
    ex(i,k)=ex(i,k)+eta*jxi(i,k)

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


  !$OMP DO 
  do k=2,nz-2
    do i=2,nx-2
      vxm(i,k)=abs(vx(i+1,k))+abs(vx(i,k))
      vzm(i,k)=abs(vz(i,k+1))+abs(vz(i,k))
    end do
  end do
  !$OMP END DO 
end subroutine vdept_hdif_coef 
!=========================================================================
end module equation
