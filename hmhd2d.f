! ===================================================================
!  A 2-D nonlinear Hall MHD code in Cartesian Geometry
!  with the capability of nonuniform grids.
!
!  The computational domain is on the x-z plane, 
!  with y the symmetry direction.
!
!  Based on an earlier UMD code by A. B. Hassam,
!  P. N. Guzdar, et al.
!
!  Some subroutines from DEBSX by Dalton Schnack, Zoran Mikic, et al. 
!  are used.
!
!  Use the MPI module from GS2
!      
!  Authors:
!  Yi-Min Huang (Princeton University)
!
! ==================================================================

#include"inc.h"
      
      program HMHD2D

      USE comm
      USE ssub
      USE io
      USE mp
      USE util
      USE equation 
      USE prob
      USE omp_lib

      IMPLICIT NONE


      include 'mpif.h'
! Initialize parallel execution

#if HYBRID==0      
      call init_mp
#else 
      call init_mp_thread
#endif

      call check_mp(nxproc,1,nzproc)

! Wall time, use MPI
      tstart=MPI_Wtime()

      if (proc0) then    ! only PE0 executes this 
! Read input parameters, and set up output file.    
        call setup
      end if

! Broadcast 
      call bctparam
      call init

#if KERNEL_TIME==1     
! Measurn the time of ntmax timestep kernel
      call barrier
      tstart_kernel=MPI_Wtime()
!$OMP PARALLEL
      do ntime=1, ntmax

        call step
      end do
!$OMP END PARALLEL
      call barrier
      tend_kernel=MPI_Wtime()
      if (iproc.eq.0) then
        print *,'Timestepping Kernel time: ',tend_kernel-tstart_kernel
      end if
      call finish
#elif KERNEL_TIME==2
! Measurn the time of writing restart files
      call barrier
      tstart_kernel=MPI_Wtime()
      write (rsofile(lenstr(rsofile)+1:),'(a,i5.5)') '.',iproc
      call wrrsfile (rsofile)         
      call barrier
      tend_kernel=MPI_Wtime()
      if (iproc.eq.0) then
        print *,'Restart Output Kernel time: ',tend_kernel-tstart_kernel
      end if
      call finish

#elif KERNEL_TIME==3
! Measurn the time of writing HDF dump file
      call barrier
      tstart_kernel=MPI_Wtime()
      call pdump
      call barrier
      tend_kernel=MPI_Wtime()
      if (iproc.eq.0) then
        print *,'HDF Output Kernel time: ',tend_kernel-tstart_kernel
      end if
      call finish


#endif
      

!$OMP PARALLEL
      do while (.not.(ifabort.or.ifend)) 

!        call newdt
        call step

!$OMP SINGLE

        ntime=ntime+1
        call histchek         

        call diags
        call rsdump

        if (proc0) then

          if (ifdisp) print*, time, den(nx/2,nz/2),iproc

!
! ****** Check if finished.
!
! ****** Set logical flag IFEND=.true. when the run is over.
!
          if (ntime.ge.ntmax.or.time.ge.tmax) then 
            ifend=.true.
          end if

!
!  **** check if blowup  
!
!          if (isnan(deni(4,4))) then 
!          if (deni(4,4).le.zero.or.deni(4,4).ge.1.e30) then 
!            ifabort=.true.
!          end if
    

         
! ***** check wall time limit, use MPI_Wtime as timer*****

          tend=MPI_Wtime()
          if (tend-tstart.ge.wall_clock_limit) then 
            ifend=.TRUE.
          end if

        end if

! check blowup.  
        call or_allreduce(ifabort)
!
        call broadcast(ifend)
!$OMP END SINGLE
      enddo
!$OMP END PARALLEL

!
! ****** Write out the restart file if requested.
!

      if ((.not.ifabort).and.(ifrsout.ne.0)) then
        write (rsofile(lenstr(rsofile)+1:),'(a,i5.5)') '.',iproc
        call wrrsfile (rsofile)   
      else
        ! if blowup -- dump plot file
        ifplxstp=1
        call pdump
      end if
!
!
! ****** Plot spatial diagnostics (unless just plotted),
! ****** plot time histories, terminate the post-processor,
! ****** and exit.
!
      if (ifplxstp.eq.0.and.pdump_last) then
        call pdump
      end if
      if (proc0) then
        call dumphist
      end if
     
      call finish





!---------------------------------------------------------------------
      CONTAINS
!#########################################################################

      subroutine setup
!
!-----------------------------------------------------------------------
!
! ****** Set up the run, input/output, etc.
!
!-----------------------------------------------------------------------
!
      USE globals
      USE ident
      USE comm
      USE nmlist
      USE util, ONLY : alpha,letter, lenstr
      USE io, ONLY : rdlin
      IMPLICIT NONE

      
!
!-----------------------------------------------------------------------
!
      INTEGER :: ierr, lrid, i, length
      CHARACTER(LEN=80) :: line(20)

!-----------------------------------------------------------------------
! ****** Output line buffer.
!
!
!-----------------------------------------------------------------------
!
! ****** Parse the command line
!
      call parse (ierr)
      if (ierr.ne.0) call abort_mp
!
! ****** Create file names based on RUNID.
!
      lrid=lenstr(runid)
      if (lrid.le.0.or.lrid.gt.9) go to 150
      if (.not.letter(runid(1:1))) go to 150
      do 100 i=1,lrid
        if (.not.alpha(runid(i:i))) go to 150
  100 continue
      go to 200
  150 continue
      write (*,*)
      write (*,*) '### ERROR in SETUP:'
      write (*,*) '### Run ID is invalid.'
      write (*,*) 'RUNID = ',runid
      call abort_mp
  200 continue
!
      if (infile.eq.' ') infile='i'//runid(1:lrid)
      outfile='o'//runid(1:lrid)
      rsofile='rs'//runid(1:lrid)
      hstroot='h'//runid(1:lrid)
!
! ****** Open the necessary files.
!
      call ffopen (8,infile,'r')
      call ffopen (9,outfile,'rw')
!
! ****** Read the namelist parameters and write them out.
!
      write (9,210)
  210 format (/,' ### Parameter values:',/)
      read (8,invars)
      write (9,invars)
#if PROBVARS!=0       
! read in problem specific parameters
      read(8, probvars)
      write(9,probvars)
#endif         
      write (9,*)
      write (9,*)

!
! ****** Rewind the input file.
!
      rewind (8)
!


!
! ****** Write out identification and parameters.
!
      call date_and_time(rdate,rtime)
      call MPI_GET_PROCESSOR_NAME(run_machine,length,ierr)
!
      write (line,300) 'Code: ',idcode, &
                       'Version: ',vers, &
                       'Updated on: ',update, &
                       'Compiled on: ', compile_machine, &
                       'Run ID: ',runid, &
                       'Run started on: ',rdate, &
                       'Run started at: ',rtime, & 
                       'Run executed on: ',run_machine 
  300 format (8(/,a,a))
      call wrlins (line,9)
!
!       write (line,310) '### Parameters:', &
!                        'NXTOT       = ',nxtot, &
!                        'NZTOT       = ',nztot, &
!                        'NXBLOCK     = ', nxblock, &
!                        'NZBLOCK     = ', nzblock, &
!                        'BLOCK_DIM_X = ', block_dim_x, &
!                        'BLOCK_DIM_Z = ', block_dim_z                       
!   310 format (/,a,/,6(/,a,i6))
      ! call wrlins (line,9)
!
! ****** Echo the input file to the output file.
!
      ! write (line,330) '### Input file contents:'
  330 format (/,a,/)
      ! call wrlins (line,3)
!
!   400 continue
!       if (rdlin(8,line(1))) then
!         call wrlins (line,1)
!         go to 400
!       end if
!
      write (*,*)
      write (9,*)

      close (8)
      return
      end subroutine setup
!#######################################################################
      subroutine parse (ierr)
!
!-----------------------------------------------------------------------
!
! ****** Parse the command line arguments.
!
!-----------------------------------------------------------------------
!
! ****** Syntax is:
!
!     <executable> [-debug] <runid> [<infile>]
!
! ****** Return IERR=0 if the command line parameters are
! ****** valid.
!
!-----------------------------------------------------------------------
!
      USE globals
      USE ident
      USE comm
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER :: ierr
!
!
!-----------------------------------------------------------------------
!
! ****** Command-line arguments.
!
       INTEGER :: iargc
!
!-----------------------------------------------------------------------
!
      CHARACTER(LEN=80) :: strng
      INTEGER :: ifrid, ififil, iarg, narg
!
!-----------------------------------------------------------------------
!
      ierr=0
!
! ****** Set the defaults.
!
      idebug=0
      ifrid=0
      ififil=0
      iarg=1
!
      narg=iargc()
!
      if (narg.eq.0) go to 900
!
! ****** Parse argumets.
!
! ****** Get the debug switch.
! ****** Syntax: -debug
!
      if (narg.ge.1) then
        call getarg (iarg,strng)
        if (strng.eq.'-debug') then
          idebug=1
          iarg=iarg+1
          narg=narg-1
        end if
      end if
!
! ****** Get the run ID.
! ****** Syntax: <runid>
!
      if (narg.ge.1) then
        call getarg (iarg,runid)
        ifrid=1
        iarg=iarg+1
        narg=narg-1
      end if
!
! ****** Get the input file name.
! ****** Syntax: <infile>
!
      if (narg.ge.1) then
        call getarg (iarg,infile)
        ififil=1
        iarg=iarg+1
        narg=narg-1
      end if
!
! ****** There should be no further arguments left.
!
      if (narg.ne.0) go to 900
!
! ****** Parsing completed.  Load arguments.
!
! ****** Check that the run ID has been specified.
!
      if (ifrid.eq.0) go to 900
!
      if (ififil.eq.0) infile=' '
!
      return
!
  900 continue
!
! ****** Error exit.
!
      ierr=1
!
      write (*,'(3(/,1x,a,a))') 'Code: ',idcode, &
                                'Version: ',vers, &
                                'Updated on: ',update
!
      write (*,'(/,1x,a,/,6(/,1x,a,i6))') &
          '### Parameters:', &
          'NXTOT       = ',nxtot, &
          'NZTOT       = ',nztot, &
          'NXBLOCK     = ', nxblock, &
          'NZBLOCK     = ', nzblock, &
          'BLOCK_DIM_X = ', block_dim_x, &
          'BLOCK_DIM_Z = ', block_dim_z          
!
      write (*,*)
      write (*,*) '### ERROR in PARSE:'
      write (*,*) '### Command line syntax error.'
      write (*,*) 'The command line syntax is:'
      write (*,*)
      write (*,*) '<executable> [-debug] <runid> [<infile>]'
!
      return
      end subroutine parse 
!#######################################################################

      subroutine finish
!
!-----------------------------------------------------------------------
!
! ****** Finish up the run and close output files.
!
!-----------------------------------------------------------------------
!
      USE globals
      USE ident
      USE comm
      USE omp_lib
      USE util, ONLY : second
      IMPLICIT NONE
!-----------------------------------------------------------------------
      REAL(RTYPE) :: tused
!-----------------------------------------------------------------------
      CHARACTER(LEN=80) :: line(20)
      INTEGER :: ierr


      call barrier 


      if (iproc.eq.0) then 
! ****** Output line buffer.
!
        call date_and_time(rdate,rtime)
!


        write (line,100) 'Code: ',idcode, &
                         'Version: ',vers, &
                         'Updated on: ',update, &
                         'Run ID: ',runid, &
                         'Run ended on: ',rdate, &
                         'Run ended at: ',rtime
  100   format (6(/,a,a))
        call wrlins (line,7)
!
        write (line,110) '### End of run ...', &
                       'NTIME = ',ntime, &
                       'TIME = ',time
  110   format (/,a,//,a,i8,/,a,1pe14.7)
        call wrlins (line,5)
!  CPU time
        tused=second()
!  Wall time, use MPI   

        tend=mpi_wtime()

!
        write (line,120) '### CPU time used = ',tused,' seconds.'
        call wrlins (line,3)
        write (line,120) '### Run time used = ',tend-tstart,' seconds.'
        call wrlins (line,3)

  120   format (/,a,f10.2,a,/)
!
! ****** Close the output file.
!
        close (9)

      end if 
      call barrier
!
! ****** Exit.
!

#if ALLOC!=0
      call deallocate_arrays()
#endif





      if (ifabort) then
        call finish_prob
        call finish_mp
        call exit (1)
      else
        call finish_prob
        call finish_mp
        call exit (0)
      end if
!
      end subroutine finish
!#########################################################################


      subroutine init
      
      USE comm

      IMPLICIT NONE
      INTEGER :: i, k, ierr
      INTEGER :: rr, bb

#if ALLOC!=0
      call allocate_arrays()
#endif




!
! ****** Determine the number of quantities requested for
! ****** plotting.
!
      i=1
      do while (plotlist(i).ne.' ')
        i=i+1
      end do 

      nplist=i-1

      i=1
      do while (bin2(i).ne.' ')
        i=i+1
      end do
 
      nbin2=i-1

      i=1
      do while (bin1(i).ne.' ')
        i=i+1
      end do
 
      nbin1=i-1



! set up MPI block info.

      bb=(iproc/blocksize)
      rr=mod(iproc,blocksize)
      iblock_x=mod(bb,nxblock)
      iblock_z=bb/nxblock

      iproc_x=iblock_x*block_dim_x+mod(rr, block_dim_x)
      iproc_z=iblock_z*block_dim_z+rr/block_dim_x
     
      ishift_x=iproc_x*nx0
      ishift_z=iproc_z*nz0    

      call init_prob

      if (rsifile.eq.' ') then
! ===========setup mesh =============================
        call meshx
        call meshz

! setup finite difference coefficients

        call fdcoef



! ============ initial conditions =====================
        time=zero
        ntime=0


        call initialize

! BCs and others
!$OMP PARALLEL
        call fin
!$OMP END PARALLEL


        tpltx_res=tpltxint
        thist_res=thistint
        trs_res=trsdump
        ipltx_res=ipltxint
        ihist_res=ihistint
        irs_res=irsdump
        
#if KERNEL_TIME==0        
        ifhststp=1
        ifplxstp=1
        ihist=1
        call diags
#endif              

      else
! read restart file
        call rdrstrt


! setup mesh

        call meshx
        call meshz

! setup finite difference coefficients

        call fdcoef


! problem specific restart
 
        call restart 

! BCs and others -- must be called after mesh and fdcoef
!$OMP PARALLEL
        call fin
!$OMP END PARALLEL

      end if 



#if TIMESTEPPING == 0
!$OMP PARALLEL       
      call firststep
!$OMP END PARALLEL
      ntime=1
#  if KERNEL_TIME==0    
      call histchek
      call diags
      call rsdump
#  endif        
#endif

      return
      end subroutine init
!####################################################################
      subroutine fdcoef
! setup finite difference coefficients

      INTEGER :: i, k
      REAL(RTYPE) :: p2, p1, m2, m1, deno, dx, dz, dx0, dz0

! finite difference coefficients based on interpolation
      do i=3,nx-2
        p2=x(i+2)-x(i)
        p1=x(i+1)-x(i)
        m1=x(i-1)-x(i)
        m2=x(i-2)-x(i)
! 5pt coef. for i-2
        deno=(m2-m1)*m2*(m2-p1)*(m2-p2)
        px51(i,1)=-m1*p1*p2/deno
        px52(i,1)=two*(p1*p2+m1*p1+m1*p2)/deno
! 5pt coef. for i-1
        deno=(m1-m2)*m1*(m1-p1)*(m1-p2)
        px51(i,2)=-m2*p1*p2/deno
        px52(i,2)=two*(p1*p2+m2*p1+m2*p2)/deno
! 5pt coef. for i
        px51(i,3)=-one/m2-one/m1-one/p1-one/p2
        px52(i,3)=two/m1/m2+two/p1/m2+two/p2/m2+two/m1/p1+two/p2/m1+two/p1/p2
! 5pt coef. for i+1
        deno=(p1-m2)*(p1-m1)*p1*(p1-p2)
        px51(i,4)=-m2*m1*p2/deno
        px52(i,4)=two*(m1*p2+m2*m1+m2*p2)/deno
! 5pt coef. for i+2
        deno=(p2-m2)*(p2-m1)*p2*(p2-p1)
        px51(i,5)=-m2*m1*p1/deno
        px52(i,5)=two*(m2*p1+m2*m1+m1*p1)/deno
! The coefficients sum to zero -- this can never be done exactly
! but this step helps reducing the error.
        px51(i,3)=-(px51(i,1)+px51(i,2)+px51(i,4)+px51(i,5))
        px52(i,3)=-(px52(i,1)+px52(i,2)+px52(i,4)+px52(i,5))
      end do
!
      do k=3,nz-2
        p2=z(k+2)-z(k)
        p1=z(k+1)-z(k)
        m1=z(k-1)-z(k)
        m2=z(k-2)-z(k)
! 5pt coef. for k-2
        deno=(m2-m1)*m2*(m2-p1)*(m2-p2)
        pz51(k,1)=-m1*p1*p2/deno
        pz52(k,1)=two*(p1*p2+m1*p1+m1*p2)/deno
! 5pt coef. for k-1
        deno=(m1-m2)*m1*(m1-p1)*(m1-p2)
        pz51(k,2)=-m2*p1*p2/deno
        pz52(k,2)=two*(p1*p2+m2*p1+m2*p2)/deno
! 5pt coef. for k
        pz51(k,3)=-one/m2-one/m1-one/p1-one/p2
        pz52(k,3)=two/m1/m2+two/p1/m2+two/p2/m2+two/m1/p1+two/p2/m1+two/p1/p2
! 5pt coef. for k+1
        deno=(p1-m2)*(p1-m1)*p1*(p1-p2)
        pz51(k,4)=-m2*m1*p2/deno
        pz52(k,4)=two*(m1*p2+m2*m1+m2*p2)/deno
! 5pt coef. for k+2
        deno=(p2-m2)*(p2-m1)*p2*(p2-p1)
        pz51(k,5)=-m2*m1*p1/deno
        pz52(k,5)=two*(m2*p1+m2*m1+m1*p1)/deno
! The coefficients sum to zero -- this can never be done exactly
! but this step helps reducing the error.
        pz51(k,3)=-(pz51(k,1)+pz51(k,2)+pz51(k,4)+pz51(k,5))
        pz52(k,3)=-(pz52(k,1)+pz52(k,2)+pz52(k,4)+pz52(k,5))
      end do

      do i=2,nx-1
        p1=x(i+1)-x(i)
        m1=x(i-1)-x(i)
! 3pt coef. for i-1
        deno=-m1*(p1-m1)
        px31(i,1)=-p1/deno
        px32(i,1)=two/deno
! 3pt coef. for i
        deno=-m1*p1
        px31(i,2)=(p1+m1)/deno
        px32(i,2)=-two/deno
! 3pt coef. for i+1
        deno=p1*(p1-m1)
        px31(i,3)=-m1/deno
        px32(i,3)=two/deno
! The coefficients sum to zero -- this can never be done exactly
! but this step helps reducing the error.
        px31(i,2)=-(px31(i,1)+px31(i,3))
        px32(i,2)=-(px32(i,1)+px32(i,3))
      end do
!
      do k=2,nz-1
        p1=z(k+1)-z(k)
        m1=z(k-1)-z(k)
! 3pt coef. for k-1
        deno=-m1*(p1-m1)
        pz31(k,1)=-p1/deno
        pz32(k,1)=two/deno
! 3pt coef. for k
        deno=-m1*p1
        pz31(k,2)=(p1+m1)/deno
        pz32(k,2)=-two/deno
! 3pt coef. for k+1
        deno=p1*(p1-m1)
        pz31(k,3)=-m1/deno
        pz32(k,3)=two/deno
! The coefficients sum to zero -- this can never be done exactly
! but this step helps reducing the error.
        pz31(k,2)=-(pz31(k,1)+pz31(k,3))
        pz32(k,2)=-(pz32(k,1)+pz32(k,3))
      end do

! coefs. for uniform grid

      dx=(xtot(nxtot)-xtot(1))/(nxtot-1)
      dz=(ztot(nztot)-ztot(1))/(nztot-1)


      px51u(1)=one/(12._RTYPE*dx)
      px51u(2)=-two/(3._RTYPE*dx)
      px51u(3)=zero
      px51u(4)=two/(3._RTYPE*dx)
      px51u(5)=-one/(12._RTYPE*dx)

      px52u(1)=-one/(12._RTYPE*dx*dx)
      px52u(2)=4._RTYPE/(3._RTYPE*dx*dx)
      px52u(3)=-5._RTYPE/(two*dx*dx)
      px52u(4)=4._RTYPE/(3._RTYPE*dx*dx)
      px52u(5)=-one/(12._RTYPE*dx*dx)


      pz51u(1)=one/(12._RTYPE*dz)
      pz51u(2)=-two/(3._RTYPE*dz)
      pz51u(3)=zero
      pz51u(4)=two/(3._RTYPE*dz)
      pz51u(5)=-one/(12._RTYPE*dz)

      pz52u(1)=-one/(12._RTYPE*dz*dz)
      pz52u(2)=4._RTYPE/(3._RTYPE*dz*dz)
      pz52u(3)=-5._RTYPE/(two*dz*dz)
      pz52u(4)=4._RTYPE/(3._RTYPE*dz*dz)
      pz52u(5)=-one/(12._RTYPE*dz*dz)

      px31u(1)=-one/(two*dx)
      px31u(2)=zero
      px31u(3)=one/(two*dx)

      px32u(1)=one/(dx*dx)
      px32u(2)=-two/(dx*dx)
      px32u(3)=one/(dx*dx)

      pz31u(1)=-one/(two*dz)
      pz31u(2)=zero
      pz31u(3)=one/(two*dz)

      pz32u(1)=one/(dz*dz)
      pz32u(2)=-two/(dz*dz)
      pz32u(3)=one/(dz*dz)

      
! 5pt dx/dx0
      dx0=(x0tot(nxtot)-x0tot(1))/(nxtot-1)
      do i=3,nx-2
        dxdx0(i)=-(x(i+2)-x(i-2))/(12._RTYPE*dx0)   &
                 +(x(i+1)-x(i-1))*two/(3._RTYPE*dx0)
      end do
! 5pt dz/dz0
      dz0=(z0tot(nztot)-z0tot(1))/(nztot-1)
      do k=3,nz-2
        dzdz0(k)=-(z(k+2)-z(k-2))/(12._RTYPE*dz0)   &
                 +(z(k+1)-z(k-1))*two/(3._RTYPE*dz0)
      end do

! set up coefs to be used in hyperdif
      do i=3,nx-2
        dx24(i)=one/(24._RTYPE*dx0)/dxdx0(i)
      end do

      do k=3,nz-2
        dz24(k)=one/(24._RTYPE*dz0)/dzdz0(k)
      end do 

#if FDTYPE==1
! now use coordinate transfornm instead of interpolation
      do i=3,nx-2
! coef. for i-2
        px51(i,1)=one/(12._RTYPE*dx0)/dxdx0(i)
! coef. for i-1
        px51(i,2)=-two/(3._RTYPE*dx0)/dxdx0(i)
! coef. for i
        px51(i,3)=zero
! coef. for i+1
        px51(i,4)=-px51(i,2)
! coef. for i+2
        px51(i,5)=-px51(i,1)
      end do
!
      do k=3,nz-2
! coef. for k-2
        pz51(k,1)=one/(12._RTYPE*dz0)/dzdz0(k)
! coef. for k-1
        pz51(k,2)=-two/(3._RTYPE*dz0)/dzdz0(k)
! coef. for k
        pz51(k,3)=zero                        
! coef. for k+1
        pz51(k,4)=-pz51(k,2)
! coef. for k+2
        pz51(k,5)=-pz51(k,1)
      end do
#endif
      end subroutine fdcoef
!==============================================================================
      subroutine bctparam
      USE ident
! broadcast input parameters to all PEs
!!!
!!! most common error: adding a variable to inputs and forgetting to broadcast it
!!!
      call broadcast(runid)
      call broadcast(rsofile)
      call broadcast(rsifile)
      call broadcast(rsdir)
      call broadcast(ifrsout)
      call broadcast(dt)
      call broadcast(ntmax)
      call broadcast(tmax)
      call broadcast(den_max)
      call broadcast(den_min)
      call broadcast(T_min)      
      call broadcast(visc)
      call broadcast(eta)
      call broadcast(etah)
      call broadcast(eta_ad)
      call broadcast(k_perp)
      call broadcast(k_para)
      call broadcast(di)
      call broadcast(gravity)
      call broadcast(fric)
      call broadcast(difden)
      call broadcast(dif4den)
      call broadcast(dif4pre)
      call broadcast(dif4pe)
      call broadcast(dif4mom)
      call broadcast(dif4b)
      call broadcast(hdifden)
      call broadcast(hdifpre)
      call broadcast(hdifpe)
      call broadcast(hdifmom)
      call broadcast(hdifb)
      call broadcast(idum)
      !idum=idum+iproc   ! so that each proc has a different random seed
      call broadcast(ipltxint)
      call broadcast(tpltxint)
      call broadcast(ihistint)
      call broadcast(thistint)
      call broadcast(plotlist)
      call broadcast(bin1)
      call broadcast(bin2)
      call broadcast(irsdump)
      call broadcast(trsdump)
#if ISOTHERMAL==1
      gamma=one
#endif
      call broadcast(gamma)
      omgamma=one-gamma
      call broadcast(omgamma)
      call broadcast(T0)
      call broadcast(bs_curr_const)
      call broadcast(rfamp)
      call broadcast(wall_clock_limit)
      call broadcast(ifdisp)
      call broadcast(qpa)
      call broadcast(qpb)
      call broadcast(qpc)
      call broadcast(qpd)
      call broadcast(qpe)
      call broadcast(qpf)
      call broadcast(qpg)
      call broadcast(qph)
      call broadcast(xmr)
      call broadcast(zmr)
      call broadcast(xl)
      call broadcast(zl)
      call broadcast(paralleldump)
      call broadcast(ipltrtype)

      end subroutine bctparam
!==============================================================================

      end program HMHD2D      
