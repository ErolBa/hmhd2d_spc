#include"inc.h"
MODULE comm
USE globals
!
!-----------------------------------------------------------------------
!
!
! ****** Constants.
!

  REAL(RTYPE), PARAMETER ::  pi=3.141592653589793_RTYPE, pi2=2._RTYPE*pi
  REAL(RTYPE), PARAMETER ::  zero=0._RTYPE, one=1._RTYPE, &
                             two=2._RTYPE, half=0.5_RTYPE

!
! ****** File names.
!
  CHARACTER(LEN=72) :: infile,outfile,rsofile, &
                       rsifile=' ', hstroot, rsdir=''

  INTEGER, PARAMETER :: mxplist=30
  ! List of variables for floating point output
  CHARACTER(LEN=8) :: plotlist(mxplist)=reshape((/' '/),(/mxplist/),(/' '/))
  ! List of variables for two-byte binning
  CHARACTER(LEN=8) :: bin2(mxplist)=reshape((/' '/),(/mxplist/),(/' '/))
  ! List of variables for one-byte binning
  CHARACTER(LEN=8) :: bin1(mxplist)=reshape((/' '/),(/mxplist/),(/' '/))

  INTEGER :: nplist, nbin2, nbin1

!
! ****** Time step related variables.
!
  REAL(RTYPE) :: time, tmax
  REAL(RTYPE) :: dt=1.E-3_RTYPE
  INTEGER :: ntime, ntmax

  INTEGER :: ifrsout=0     ! 1 => dump restart  in the end

  LOGICAL ::  ifabort=.false.,   &
              ifend=.false.

!
! ******** Mesh related 
!


  REAL(RTYPE) :: xmr=1._RTYPE,    &  ! ratio of mesh size at edge to 
                 zmr=1._RTYPE        ! that at center

  REAL(RTYPE) :: xl=1._RTYPE, zl=1._RTYPE                
!
! ****** Physical parameters.
!

  REAL(RTYPE) :: &
    visc=1.e-3_RTYPE, &
    eta=1.e-3_RTYPE,  &
    etah=zero,  &  ! hyper-resistivity 
    eta_ad=zero, &  ! ambipolar diffusion 
      ! Thermal conductivity 
      ! assume perpendicular conductivity =min(k_perp/B^2,k_para)
    k_para=1.e-1_RTYPE, & 
    k_perp=1.e-3_RTYPE, & 
    ! k_para=0.0, & 
    ! k_perp=0.0, & 
                  
    di=zero,    &  ! ion skin depth
    gravity=zero,  &  ! gravity 
    difden=zero,        &  !Density diffusion, this is not physical
    fric=zero,          &    !friction coef.
    bs_curr_const = zero, & !Constant in front of the bootstrap current term in Ohm's law
    T0=one,    &   !temperature for isothermal EOS
    Te0=one,   &   !electron temperature for isothermal EOS
    rfamp=zero,     &   ! random forcing amplitude
    gamma=1.66666666666666667_RTYPE,  & ! adiabatic index
    omgamma                             ! 1-gamma


! Min & Max density and temperature

  REAL(RTYPE) :: den_min=zero, &   ! minimum density
                 T_min=zero      ! minimum temperature

  REAL(RTYPE) :: den_max=1.E10_RTYPE  ! Max density. 

!                    
! ****** Hyper diffusion ****************************
!
  REAL(RTYPE) ::  hdifden=0._RTYPE,  &  !v dependent hyperdiffusion for den
                  hdifpre=0._RTYPE,  &  !v dependent hyperdiffusion for pre
                  hdifpe =0._RTYPE,  &  !v dependent hyperdiffusion for pe
                  hdifmom=0._RTYPE,  &  !v dependent hyperdiffusion for momentum
                  hdifb  =0._RTYPE,  &  !v dependent hyperdiffusion for magnetic field
                  dif4den=0._RTYPE, &  !4th order diffusion for den
                  dif4pre=0._RTYPE, &  !4th order diffusion for pre
                  dif4pe =0._RTYPE, &  !4th order diffusion for pe
                  dif4mom=0._RTYPE, &  !4th order diffusion for momentum
                  dif4b  =0._RTYPE    !4th order diffusion for magnetic field

! ============ Large 2D arrays ==============================
!     Total 72 Arrays
!

#if ALLOC==0
  REAL(RTYPE), DIMENSION(nx, nz) :: &
#else
  REAL(RTYPE), DIMENSION(:, :), ALLOCATABLE :: &
#endif
    
! ****** Field Variables
    den, deni, &
    pux, puxi, &
    puz, puzi, &
    by, byi,   &
    puy, puyi, &
    psi, psii, &
    psi0, psitot, &  ! background and total
    pre, prei, &   ! ion pressure
    pe, pei, &  ! electron pressure
! Deduced variables
    pt,      &    ! total pressure
    bxi, bzi,  &
    jxi, jyi, jzi, &
    ex, ey, ez, &   ! E field
    Te, Ti, &
! Velocities
    vx, vy, vz,  &
    vxm, vzm, &
    vex, vey, vez, &
    vexm, vezm, &
!      
    heat_exchange, &  ! heat flux from electron to ion
    eta_perp,      &  ! perpendicular resistivity
    electron_heating, &
    ion_heating, &
!
    temp, &   ! Dummy array

! magnetic force       

    fbx, fby, fbz, & 

! ****** Flux ******

    fx, fy, fz,  &
    divf, &

! ****** visc*den *****

    vden, &


! ****** Stree Tensor  ******************               
      
    fxx, fxy, fxz, &
              fyz, &
    fzx, fzy, fzz, &   

! ****** Grad V Tensor  *********     
 
    vxx, vxy, vxz, &
              vyz, &
              vzz, &   

! ****  Div stress tensor                                              

    divfx, divfy, divfz, &

! ****** Hyper-diffusion profile ********************
 
    dif4_prof, & 

! ****** Friction  profile ********************
     
    fric_prof, &
    
    misc
  
! ==================================================== 
  
! ****** scratch matrices and 

  REAL(RTYPE) :: sx(3:nxtot-2), sz(3:nztot-2)
  REAL(RTYPE) :: ssx(3:nx-2), ssz(3:nz-2)


! ****** Coordinates
!
  REAL(RTYPE) :: x(nx),z(nz)        ! This can be nonunif. 
  REAL(RTYPE) :: x0(nx), z0(nz)     ! uniform mesh here 
  REAL(RTYPE) :: xtot(nxtot), ztot(nztot)
  REAL(RTYPE) :: x0tot(nxtot), z0tot(nztot)
  REAL(RTYPE) :: dxdx0(nx),dzdz0(nz)   ! dx/dx0, dz/dz0

!
! ***** Finite difference coefficients
!
! 5-points approx. d/dx, d/dy,  and d/dz
  REAL(RTYPE) :: px51(nx,5)
  REAL(RTYPE) :: pz51(nz,5)
! 5-points approx. d^2/dx^2, d^2/dy^2, and d^2/dz^2
  REAL(RTYPE) :: px52(nx,5)
  REAL(RTYPE) :: pz52(nz,5)
! 5-points approx. coefs. when grid are uniform
  REAL(RTYPE) :: px51u(5), px52u(5)
  REAL(RTYPE) :: pz51u(5), pz52u(5)
! 3-points approx. d/dx, d/dy,  and d/dz
  REAL(RTYPE) :: px31(nx,3)
  REAL(RTYPE) :: pz31(nz,3)
! 3-points approx. d^2/dx^2, d^2/dy^2, and d^2/dz^2
  REAL(RTYPE) :: px32(nx,3)
  REAL(RTYPE) :: pz32(nz,3)
! 3-points approx. coefs. when grid are uniform
  REAL(RTYPE) :: px31u(3), px32u(3)
  REAL(RTYPE) :: pz31u(3), pz32u(3)

!
! To be used in hyperdif
!
  
  REAL(RTYPE) :: dx24(nx), dz24(nz)

! ****** Time Limit
!

  REAL(DP) :: wall_clock_limit=1.e30  ! terminate the code and dump
                                      ! result when the execution time
                                      ! goes beyond this limit. In unit 
                                      ! second.  Need openmp as a timer.
  REAL(DP) :: tend, tstart, tend_kernel, tstart_kernel

! ***** Random seed

  INTEGER(ILWORD) :: idum=312560


!************ Free Parameters **********************
  REAL(RTYPE) :: qpa, qpb, qpc, qpd, qpe, qpf, qpg, qph

! ******* Misc *************************************
 
  INTEGER :: idebug
  LOGICAL :: ifdisp=.FALSE. ! whether display info in each time step

!
! ****** Control Plots, Dumps, and Time history 
!

  REAL(RTYPE) :: tpltxint=100._RTYPE,   &
                 tpltx_res=0._RTYPE,    &
                 thistint=0._RTYPE,     &
                 thist_res=0._RTYPE,    &
                 trsdump=0._RTYPE,      &
                 trs_res=0._RTYPE

  INTEGER ::  ipltrtype=1     ! plot data type, 2 ->double, 1 ->single 

  INTEGER ::  ipltxint=0,   &
              ipltx_res=0,  &
              ifplxstp, &
              ihistint=1, &
              ihist_res=0, &
              ifhststp,   &
              ihist=0, &
              irsdump=0, &
              irs_res=0

  INTEGER :: pdump_seq=0   ! output file sequence number

  LOGICAL :: pdump_last=.FALSE.  ! True -> dump the last time slice before quit

!
!  **** MPI related variables
!
  INTEGER :: iblock_x, iblock_z

  INTEGER :: iproc_x, iproc_z

  INTEGER :: ishift_x, ishift_z

  LOGICAL :: paralleldump=.TRUE.  ! use parallel hdf for dump file or not


CONTAINS
  
#if ALLOC!=0
  subroutine allocate_arrays()

    allocate(den(nx, nz))
    allocate(deni(nx, nz))
    allocate(pux(nx, nz))
    allocate(puxi(nx, nz))
    allocate(puz(nx, nz))
    allocate(puzi(nx, nz))  
    allocate(by(nx, nz)) 
    allocate(byi(nx, nz)) 
    allocate(puy(nx, nz))
    allocate(puyi(nx, nz))  
    allocate(psi(nx, nz))
    allocate(psii(nx, nz))  
    allocate(psi0(nx, nz))
    allocate(psitot(nx, nz)) 
    allocate(pre(nx, nz))
    allocate(prei(nx, nz))
    allocate(pe(nx, nz)) 
    allocate(pei(nx, nz))
    allocate(pt(nx, nz))    
    allocate(bxi(nx, nz))
    allocate(bzi(nx, nz))  
    allocate(jxi(nx, nz))
    allocate(jyi(nx, nz))
    allocate(jzi(nx, nz)) 
    allocate(ex(nx, nz)) 
    allocate(ey(nx, nz))
    allocate(ez(nx, nz)) 
    allocate(Te(nx, nz))
    allocate(Ti(nx, nz))  
    allocate(vx(nx, nz))
    allocate(vy(nx, nz)) 
    allocate(vz(nx, nz)) 
    allocate(vxm(nx, nz))
    allocate(vzm(nx, nz)) 
    allocate(vex(nx, nz))
    allocate(vey(nx, nz))
    allocate(vez(nx, nz)) 
    allocate(vexm(nx, nz))
    allocate(vezm(nx, nz))
    allocate(heat_exchange(nx, nz))
    allocate(eta_perp(nx, nz)) 
    allocate(electron_heating(nx, nz))
    allocate(ion_heating(nx, nz)) 
    allocate(temp(nx, nz))
    allocate(fbx(nx, nz)) 
    allocate(fby(nx, nz))
    allocate(fbz(nx, nz))  
    allocate(fx(nx, nz))
    allocate(fy(nx, nz))
    allocate(fz(nx, nz)) 
    allocate(divf(nx, nz)) 
    allocate(vden(nx, nz)) 
    allocate(fxx(nx, nz))
    allocate(fxy(nx, nz))
    allocate(fxz(nx, nz)) 
    allocate(fyz(nx, nz)) 
    allocate(fzx(nx, nz))
    allocate(fzy(nx, nz))
    allocate(fzz(nx, nz))    
    allocate(vxx(nx, nz))
    allocate(vxy(nx, nz))
    allocate(vxz(nx, nz))
    allocate(vyz(nx, nz))  
    allocate(vzz(nx, nz))   
    allocate(divfx(nx, nz))
    allocate(divfy(nx, nz))
    allocate(divfz(nx, nz)) 
    allocate(dif4_prof(nx, nz))   
    allocate(fric_prof(nx, nz))
    allocate(misc(nx, nz)) 
    return
  end subroutine allocate_arrays  
! ==========================================================
  subroutine deallocate_arrays()
    deallocate(den, deni, &
               pux, puxi, &
               puz, puzi, &
               by, byi,   &
               puy, puyi, &
               psi, psii, &
               psi0, psitot, &  ! background and total
               pre, prei, &   ! ion pressure
               pe, pei, &  ! electron pressure
               ! Deduced variables
               pt,      &    ! total pressure
               bxi, bzi,  &
               jxi, jyi, jzi, &
               ex, ey, ez, &   ! E field
               Te, Ti, &
               ! Velocities
               vx, vy, vz,  &
               vxm, vzm, &
               vex, vey, vez, &
               vexm, vezm, &
               !      
               heat_exchange, &  ! heat flux from electron to ion
               eta_perp,      &  ! perpendicular resistivity
               electron_heating, &
               ion_heating, &
               !
               temp, &   ! Dummy array
               ! magnetic force       
               fbx, fby, fbz, & 
               ! ****** Flux ******
               fx, fy, fz,  &
               divf, &
               ! ****** visc*den *****
               vden, &
               ! ****** Stree Tensor  ******************               
               fxx, fxy, fxz, &
                         fyz, &
               fzx, fzy, fzz, &   
               ! ****** Grad V Tensor  *********     
               vxx, vxy, vxz, &
                         vyz, &
                         vzz, &   
               ! ****  Div stress tensor                                 
               divfx, divfy, divfz, &
               ! ****** Hyper-diffusion profile ********************
               dif4_prof, & 
               ! ****** Friction  profile ********************
               fric_prof,  &

               misc &
              )  
    return
  end subroutine deallocate_arrays  
#endif
END MODULE comm
