MODULE nmlist

USE comm
!
!-----------------------------------------------------------------------
!
! ****** Namelist input variables.
!
      namelist /invars/ rsifile, rsdir, ifrsout , dt, ntmax, tmax, maxden, &
       visc, eta, k_perp, k_para, etah, eta_ad, fric, di, gravity,difden, &
       hdifden, hdifpre, &
       hdifpe, hdifmom, hdifb, &
       dif4den, dif4pre, dif4pe, dif4mom, dif4b, idum, &
       ipltxint,tpltxint,ihistint,thistint, ipltrtype, pdump_last, &
       plotlist,bin2, bin1, irsdump,trsdump, &
       gamma, T0, bs_curr_const, rfamp, wall_clock_limit, ifdisp, &
! min & max values
       den_min, den_max, T_min, &       
! free parameters 
       qpa, qpb, qpc, qpd, qpe, qpf, qpg, qph, &
! mesh related
       xmr,zmr, xl, zl, & 
!
       paralleldump
!
!-----------------------------------------------------------------------
!

END MODULE nmlist
