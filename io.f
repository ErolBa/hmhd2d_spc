#include"inc.h"

      MODULE io
! I/O subroutines
      USE mp

      CONTAINS
!#######################################################################
      subroutine wrlins (line,nlin)
      USE util, ONLY : lenstr
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! ****** Write NLIN lines from array LINE to the TTY and unit 9.
!
!-----------------------------------------------------------------------
!
      CHARACTER(LEN=*) :: line(*)
      INTEGER :: nlin
!----------------------------------------------------------------------
      INTEGER :: i, ll
!-----------------------------------------------------------------------
!
      do i=1,nlin
        ll=max0(1,lenstr(line(i)))
        write (*,*) line(i)(1:ll)
        if (proc0) then
          write (9,*) line(i)(1:ll)
        end if
      enddo
!
      return
      end subroutine wrlins
!#######################################################################
      subroutine rdrstrt
!
!-----------------------------------------------------------------------
!
! ****** Read the restart file.
!
!-----------------------------------------------------------------------
!
      USE globals
      USE comm
      USE util, ONLY : lenstr
      IMPLICIT NONE

      INTEGER :: id, lf, ierr
      CHARACTER(LEN=72) :: fname, fname1
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
! ****** Open the restart file.
!

      lf=lenstr(rsifile)
      write (fname,'(a,a,i5.5)') rsifile(1:lf),'.',iproc
      fname1=rsdir(1:lenstr(rsdir))//fname
      
      id=200+iproc
      ! open(unit=id,status='old',file=fname1,form='unformatted',iostat=ierr)
      open(unit=id,status='old',file=fname1,form='unformatted',iostat=ierr)

      if (ierr.ne.0) then
        if (proc0) then
          write (*,*)
          write (*,*) 'RSDIR may not exist. Try currernt directory instead.'
        end if
        
        open(unit=id,status='old',file=fname,form='unformatted',iostat=ierr)
      end if  

!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in RDRSTRT:'
        write (*,*) '### Error while opening the restart file at PROC ', iproc
        if (ierr.eq.1) then
          write (*,*) 'The file does not exist.'
        else
          write (*,*) 'IERR = ',ierr
        end if
        write (*,*) 'File name: ',rsifile(1:lf)
        call abort_mp
        call exit(1)
      end if
!
      if (proc0) then
        write (9,*)
        write (9,*) '### COMMENT from RDRSTRT:'
        write (9,*) 'Reading from restart file: ',rsifile(1:lf)
      end if
!

      read (id) time
      read (id) tpltx_res
      read (id) thist_res
      read (id) trs_res
      read (id) ipltx_res
      read (id) ihist_res
      read (id) irs_res
      read (id) pdump_seq
      read (id) ishift_x
      read (id) ishift_z
      read (id) x
      read (id) z
      read (id) xtot
      read (id) ztot
      read (id) deni
      read (id) puxi
      read (id) puzi
#if TWO_AND_HALF==1
      read (id) byi
      read (id) puyi
#endif
      read (id) psii
#if ISOTHERMAL==0
      read (id) prei
#endif     
#if (TWOFLUID==1 && ISOTHERMAL==0 && SPE==1)
      read (id) pei
#endif

!
      close(unit=id,iostat=ierr)

!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in RDRSTRT:'
        write (*,*) '### Error while closing the restart file at PROC ', iproc
        write (*,*) 'File name: ',rsifile(1:lf)
        write (*,*) 'IERR (from BCLOSE) = ',ierr
        call abort_mp
        call exit(1)
      end if
!

      return
      end subroutine rdrstrt 
!#######################################################################
      subroutine wrrsfile (fname)
!
!-----------------------------------------------------------------------
!
! ****** Write a restart file in file FNAME.
!
!-----------------------------------------------------------------------
!
      USE globals
      USE comm
      USE util, ONLY : lenstr
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      CHARACTER(LEN=*) :: fname 
      CHARACTER(LEN=72) :: fname1

      INTEGER :: lf, ierr, id
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
! ****** Create the restart file.
!
      fname1=rsdir(1:lenstr(rsdir))//fname
      lf=lenstr(fname)
!
      id=200+iproc
      open(unit=id,status='unknown',file=fname1,form='unformatted',iostat=ierr)
!
      if (ierr.ne.0) then
        if (proc0) then
          write (*,*)
          write (*,*) 'RSDIR may not exist. Use currernt directory instead.'
        end if  
        open(unit=id,status='unknown',file=fname,form='unformatted',iostat=ierr)
      end if  

      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRRSTRT:'
        write (*,*) '### Error while creating the restart file at PROC', iproc
        write (*,*) 'File name: ',fname(1:lf)
        write (*,*) 'IERR = ',ierr
        call abort_mp
        call exit(1)
      end if
!
      if (proc0) then
        write (9,*)
        write (9,*) '### COMMENT from WRRSTRT:'
        write (9,*) 'Writing to restart file: ',fname(1:lf)
        write (9,*) 'TIME = ',time
      end if 
      write (id) time
      write (id) tpltx_res
      write (id) thist_res
      write (id) trs_res
      write (id) ipltx_res
      write (id) ihist_res
      write (id) irs_res
      write (id) pdump_seq
      write (id) ishift_x
      write (id) ishift_z
      write (id) x
      write (id) z
      write (id) xtot
      write (id) ztot
      write (id) deni
      write (id) puxi
      write (id) puzi
#if TWO_AND_HALF==1
      write (id) byi
      write (id) puyi
#endif
      write (id) psii
#if ISOTHERMAL==0
      write (id) prei
#endif     
#if (TWOFLUID==1 && ISOTHERMAL==0 && SPE==1)
      write (id) pei
#endif
 

!
      close(unit=id,iostat=ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRRSTRT:'
        write (*,*) '### Error while closing the restart file at PROC', iproc
        write (*,*) 'File name: ',fname(1:lf)
        write (*,*) 'IERR (from BCLOSE) = ',ierr
        call abort_mp
        call exit(1)
      end if
!
      return
      end subroutine wrrsfile
!#######################################################################
      subroutine ffopen (iun,fname,mode)
      USE util, ONLY : lenstr
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! ****** Open file FNAME and link it to unit IUN.
!
!-----------------------------------------------------------------------
!
! ****** When MODE='r', the file must exist.
! ****** When MODE='w', the file is created.
! ****** When MODE='rw', the file must exist, but can be overwritten.
!
!-----------------------------------------------------------------------
!
      CHARACTER(LEN=*) :: fname,mode
      INTEGER :: iun

      INTEGER :: lf
!
!-----------------------------------------------------------------------
!
      lf=lenstr(fname)
!
      if (mode.eq.'r') then
        open (iun,file=fname(1:lf),status='old',err=900)
      else if (mode.eq.'rw') then
        open (iun,file=fname(1:lf),err=900)
      else if (mode.eq.'w') then
        open (iun,file=fname(1:lf),status='new',err=900)
      else
        write (*,*)
        write (*,*) '### ERROR in FFOPEN:'
        write (*,*) '### Invalid MODE requested.'
        write (*,*) 'MODE = ',mode
        write (*,*) 'File name = ',fname(1:lf)
        call abort_mp
        call exit(1)

      end if
!
      return
!
  900 continue
!
      write (*,*)
      write (*,*) '### ERROR in FFOPEN:'
      write (*,*) '### Error while opening file.'
      write (*,*) 'File name = ',fname(1:lf)
      write (*,*) 'MODE =  ',mode
      call abort_mp
      call exit(1)

      end subroutine ffopen
!#######################################################################
      function rdlin (iun,line)
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! ****** Read a line of text from unit IUN.
!
! ****** If a line was read successfully, return with
! ****** RDLIN=.true.; otherwise, set RDLIN=.false.
!
!-----------------------------------------------------------------------
!
      LOGICAL :: rdlin
      INTEGER, INTENT(IN) :: iun
      CHARACTER(LEN=*) ::  line
!
!-----------------------------------------------------------------------
!
      read (iun,'(a)',end=100,err=100) line
!
      rdlin=.true.
      return
!
  100 continue
!
! ****** Error while reading the line.
!
      rdlin=.false.
!
      return
      end function rdlin
!#######################################################################
      subroutine rsdump
!
!-----------------------------------------------------------------------
!
! ****** Write a restart file if on a dump interval.
!
!-----------------------------------------------------------------------
!
      USE globals
      USE ident
      USE comm
      USE util, ONLY : lenstr
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! ****** Restart file name and sequence number.
!
      CHARACTER(LEN=72) ::  rstname
      INTEGER, SAVE :: iseq=0
      
      INTEGER :: lrid, ifrsstp
      REAL(RTYPE) :: trem

!
!-----------------------------------------------------------------------
!
! ****** Check if this is a restart file dump interval.
!
! ****** Disable restart file dump at time step intervals
! ****** if IRSDUMP.le.0.
!

      ifrsstp=0

      if (irsdump>0) then
        irs_res=irs_res-1
        if (irs_res.eq.0) then
          ifrsstp=1
          irs_res=irsdump
        end if
      end if

!
! ****** Disable restart file dump at time intervals
! ****** if TRSDUMP.le.0.
!


      if (trsdump>zero) then
        trs_res=trs_res-dt
        if (trs_res<=zero) then
          ifrsstp=1
          trs_res=trs_res+trsdump
        end if
      end if
      
!
      if (ifrsstp.eq.0) return
!
! ****** Create a restart file name of the form: rs<runid>###.
!
      iseq=iseq+1
      rstname='rs'
      lrid=lenstr(runid)
      write (rstname(3:),'(a,i3.3,a,i5.5)') runid(1:lrid),iseq,'.',iproc
!
! ****** Write the restart file.
!
      call wrrsfile (rstname)
!
      return
      end subroutine rsdump
!#######################################################################
      subroutine pdump
!
!-----------------------------------------------------------------------
!
! ****** Dump requested fields at the present time.
!
!-----------------------------------------------------------------------
!
      USE globals
      USE comm
      USE util
      USE ssub 
      USE HDF5
      IMPLICIT NONE
      include 'mpif.h'
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      CHARACTER(LEN=32) :: fname
      CHARACTER(Len=4) :: ch3
!
      INTEGER :: ierr
      INTEGER(HID_T) :: file_id, plist_id

      INTEGER(HID_T) :: dset_id, dspace_id, memspace_id, rtype_id, rtype_id1
      INTEGER(HSIZE_T) :: data_dims(2), buf_dims(2), offset(2)=(/0,0/)
      INTEGER(HSIZE_T) :: data_dim(1)

      REAL(SP), ALLOCATABLE :: f1(:,:)
      REAL(DP), ALLOCATABLE :: f2(:,:)
      REAL(RTYPE), ALLOCATABLE :: f3(:), f4(:,:)
      INTEGER(IBYTE), ALLOCATABLE :: b1(:,:), b2(:,:)   ! for binning data
      INTEGER(IBYTE), ALLOCATABLE :: b3(:)  ! for binning data
      INTEGER(IWORD), ALLOCATABLE :: w1(:,:), w2(:,:)   ! for binning data
      INTEGER(IWORD), ALLOCATABLE :: w3(:)  ! for binning data
      INTEGER(HSIZE_T) :: stride(2)=(/1,1/), count(2)=(/1,1/) 
      INTEGER :: xs,xe,zs,ze,nnx,nnz

      CHARACTER(LEN=16) :: attr_name(4*mxplist)  ! store data name for offset, scale.
      REAL(RTYPE) :: hdf_tstart, hdf_tend, hdf_tend1, hdf_tend2, hdf_tend3      
      REAL(RTYPE) :: attr_value(4*mxplist)    ! the corresponding values
      INTEGER :: attr_id=0
      INTEGER :: i,k
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      call barrier
      hdf_tstart=MPI_Wtime()
      
      pdump_seq=pdump_seq+1
      if (pdump_seq.gt.9999) pdump_seq=pdump_seq-9999
      write (ch3,'(i4.4)') pdump_seq      

!
! Initialize FORTRAN interface
!
      call h5open_f(ierr)
!



! create the data file

!
! Setup file access property list with parallel I/O access.
!
      fname='data'//ch3//'.hdf'
      

      if (paralleldump) then 
        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
        CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
        !
        ! Create the file collectively.
        !
        call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,ierr,access_prp=plist_id)

        call h5pclose_f(plist_id,ierr)
!

        nnx=nx0
        nnz=nz0
        xs=3
        xe=nx-2
        zs=3
        ze=nz-2
        if (iproc_x.eq.0) then
          nnx=nnx+2
          xs=1 
          offset(1)=0
        else
          offset(1)=nx0*iproc_x+2
        end if

        if (iproc_x.eq.(nxproc-1)) then
          nnx=nnx+2
          xe=nx
        end if
 
        if (iproc_z.eq.0) then
          nnz=nnz+2
          zs=1
          offset(2)=0
        else
          offset(2)=nz0*iproc_z+2
        end if

        if (iproc_z.eq.(nzproc-1)) then
          nnz=nnz+2
          ze=nz
        end if 


        data_dims(1)=nxtot            ! data dimension in file
        data_dims(2)=nztot
        buf_dims(1)=nnx
        buf_dims(2)=nnz
        if (nplist.gt.0) then
          if (ipltrtype.eq.1) then
            allocate(f1(nnx,nnz),STAT=ierr)
            rtype_id=H5T_NATIVE_REAL
          else  
            allocate(f2(nnx,nnz),STAT=ierr)
            rtype_id=H5T_NATIVE_DOUBLE
          end if
        end if
        if (nbin1.gt.0) then 
          allocate(b1(nnx,nnz),STAT=ierr)
          allocate(b2(nx,nz),STAT=ierr)
        end if
        if (nbin2.gt.0) then 
          allocate(w1(nnx,nnz),STAT=ierr)
          allocate(w2(nx,nz),STAT=ierr)
        end if

      else
        if (nplist.gt.0) then
          ! only f3 in PE0 will be used.  But zaphod complains
          ! if other PEs do not allocate.
          allocate(f3(nx*nz*nproc),STAT=ierr)   
        end if
        if (nbin1.gt.0) then 
          allocate(b2(nx,nz),STAT=ierr) 
          allocate(b3(nx*nz*nproc),STAT=ierr) ! same here, only PE0 uses it
        end if 
        if (nbin2.gt.0) then 
          allocate(w2(nx,nz),STAT=ierr) 
          allocate(w3(nx*nz*nproc),STAT=ierr) ! same here, only PE0 uses it
        end if 
        if (proc0) then  ! PE0 opens the dump file
          call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,ierr)
          ! buffer for gathering data
          if (nplist.gt.0) then
            allocate(f4(nx,nz),STAT=ierr)
            if (ipltrtype.eq.1) then
              allocate(f1(nxtot,nztot),STAT=ierr)
              rtype_id=H5T_NATIVE_REAL
            else  
              allocate(f2(nxtot,nztot),STAT=ierr)
              rtype_id=H5T_NATIVE_DOUBLE
            end if
          end if 
          if (nbin1.gt.0) then 
            allocate(b1(nxtot,nztot),STAT=ierr)
          end if 
          if (nbin2.gt.0) then
            allocate(w1(nxtot,nztot),STAT=ierr)
          end if
          data_dims(1)=nxtot
          data_dims(2)=nztot
          buf_dims(1)=nxtot
          buf_dims(2)=nztot


        end if 
      end if


      if (RTYPE.eq.DP) then
        rtype_id1=H5T_NATIVE_DOUBLE
      else
        rtype_id1=H5T_NATIVE_REAL
      end if 

      call wr2d (deni,'den')
!
      call wr2d (puxi,'pux')
!
      call wr2d (puyi,'puy')
!
      call wr2d (puzi,'puz')
!
      call wr2d (vx,'vx')
!                     
      call wr2d (vy,'vy')
!                   
      call wr2d (vz,'vz')
!                
      call wr2d (vex,'vex')
!                    
      call wr2d (vey,'vey')
!                     
      call wr2d (vez,'vez')
!                                                                       
      call wr2d (bxi,'bx')
!
      call wr2d (byi,'by')
!
      call wr2d (bzi,'bz')
!
      call wr2d (psii,'psi')
!
      call wr2d (psitot,'psitot')
!
      call wr2d (jyi,'jy')
!
      call wr2d (jxi,'jx')
!
      call wr2d (jzi,'jz')
!
      call wr2d (prei,'pre')
!
      call wr2d (pei,'pe')
!
      call wr2d (Ti,'Ti')
!
      call wr2d (Te,'Te')
!
      call wr2d (ex,'ex')
!
      call wr2d (ey,'ey')
!
      call wr2d (ez,'ez')

      call wr2d (misc, 'misc')
! div v
      if (has(nplist,plotlist,'divv').or.  &
          has(nbin2,bin2,'divv').or.   &
          has(nbin1,bin1,'divv')) then
        !$OMP PARALLEL DO  
        do k=3,nz-2
        do i=3,nx-2
          temp(i,k)=ddx(vx,i,k)+ddz(vz,i,k)
        end do
        end do  
        !$OMP END PARALLEL DO  
        call wr2d(temp,'divv')
      end if
! vorticity in y
      if (has(nplist,plotlist,'vorty').or.  &
          has(nbin2,bin2,'vorty').or.   &
          has(nbin1,bin1,'vorty')) then
        !$OMP PARALLEL DO  
        do k=3,nz-2
        do i=3,nx-2
          temp(i,k)=ddz(vx,i,k)-ddx(vz,i,k)
        end do
        end do  
        !$OMP END PARALLEL DO  
        call wr2d(temp,'vorty')
      end if


!
      if (paralleldump) then 
!
        hdf_tend1=MPI_Wtime()
!        
        call h5fclose_f(file_id,ierr)
!
        if (nplist.gt.0) then
          if (ipltrtype.eq.1) then
            deallocate(f1,STAT=ierr)
          else
            deallocate(f2,STAT=ierr)
          end if
        end if

        if (nbin1.gt.0) then 
          deallocate(b1,STAT=ierr)
          deallocate(b2,STAT=ierr)
        end if
        if (nbin2.gt.0) then 
          deallocate(w1,STAT=ierr)
          deallocate(w2,STAT=ierr)
        end if
!   
      else
        if (nplist.gt.0) then
          deallocate(f3,STAT=ierr)
        end if
        if (nbin1.gt.0) then 
          deallocate(b3,STAT=ierr)
          deallocate(b2,STAT=ierr)
        end if
        if (nbin2.gt.0) then
          deallocate(w3,STAT=ierr)
          deallocate(w2,STAT=ierr)
        end if

        if (proc0) then  ! PE0 closes the dump file

          hdf_tend1=MPI_Wtime()

          call h5fclose_f(file_id,ierr)
          if (nplist.gt.0) then
            deallocate(f4,STAT=ierr)
            if (ipltrtype.eq.1) then
              deallocate(f1,STAT=ierr)
            else
              deallocate(f2,STAT=ierr)
            end if
          end if
          if (nbin1.gt.0) then 
            deallocate(b1,STAT=ierr)
          end if
          if (nbin2.gt.0) then
            deallocate(w1,STAT=ierr)
          end if
        end if 
      end if

      hdf_tend2=MPI_Wtime()

! write scale, offset, time, and coordinates by PE0
      if (iproc.eq.0) then

! create file to store time and space
!
        call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,ierr)
! write time
        data_dim(1)=1
        call h5screate_simple_f(1,data_dim,dspace_id,ierr)
        call h5dcreate_f(file_id,'ttime',rtype_id1,dspace_id, &
                       dset_id,ierr)
        call h5dwrite_f(dset_id,rtype_id1,time,data_dim,ierr)
        call h5sclose_f(dspace_id,ierr)
        call h5dclose_f(dset_id,ierr)

! write x

        data_dim(1)=nxtot
        call h5screate_simple_f(1,data_dim,dspace_id,ierr)
        call h5dcreate_f(file_id,'x',rtype_id1,dspace_id, &
                       dset_id,ierr)
        call h5dwrite_f(dset_id,rtype_id1,xtot,data_dim,ierr)
        call h5sclose_f(dspace_id,ierr)
        call h5dclose_f(dset_id,ierr)
! write z

        data_dim(1)=nztot
        call h5screate_simple_f(1,data_dim,dspace_id,ierr)
        call h5dcreate_f(file_id,'z',rtype_id1,dspace_id, &
                       dset_id,ierr)
        call h5dwrite_f(dset_id,rtype_id1,ztot,data_dim,ierr)
        call h5sclose_f(dspace_id,ierr)
        call h5dclose_f(dset_id,ierr)

! write scale and offset
! Write offset and scale
        data_dim(1)=1
        do while (attr_id.gt.0)
          call h5screate_simple_f(1,data_dim,dspace_id,ierr)
          call h5dcreate_f(file_id,attr_name(attr_id),rtype_id1,dspace_id, &
                       dset_id,ierr)
          call h5dwrite_f(dset_id,rtype_id1,attr_value(attr_id),data_dim,ierr)
          call h5sclose_f(dspace_id,ierr)
          call h5dclose_f(dset_id,ierr)
          attr_id=attr_id-1
        end do
        call h5fclose_f(file_id,ierr)

        
      end if

      hdf_tend3=MPI_Wtime()

!
!
! Close FORTRAN interface.
!
      call h5close_f(ierr)

!

      call barrier

      hdf_tend=MPI_Wtime()


!
      if (proc0) then
        write (9,*)
        write (9,*) '### COMMENT from PDUMP:'
        write (9,*) '### Wrote out the requested fields.'
        write (9,*) 'File sequence number = ', pdump_seq
        write (9,*) 'NTIME = ',ntime
        write (9,*) 'TIME = ',time
        write (9,*) 'WRITING TIME of main arrays (SEC) = ',hdf_tend1-hdf_tstart
        write (9,*) 'TIME for closing file (SEC) = ',hdf_tend2-hdf_tend1 
        write (9,*) 'WRITING TIME of auxiliary info. (SEC) = ',hdf_tend3-hdf_tend2
        write (9,*) 'Total HDF WRITING TIME (SEC) = ',hdf_tend-hdf_tstart
      end if
!
      return
!
      CONTAINS
!#######################################################################
      subroutine wr2d (ff,dname)
!
!-----------------------------------------------------------------------
!
! ****** Write a 2D field to a file.
!
!-----------------------------------------------------------------------
!
      USE globals
      USE comm
      USE util, ONLY : lenstr
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      REAL(RTYPE) :: ff(nx,nz)
      CHARACTER(len=*) :: dname
  


!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      if (has(nplist,plotlist,dname)) then
        if (paralleldump) then
          call wr2dhdf (ff,dname)
        else
          call wr2dhdf_seq (ff,dname)
        end if
      else if (has(nbin2,bin2,dname)) then
        if (paralleldump) then
          call wr2dhdf_bin (ff,dname,2)
        else
          call wr2dhdf_bin_seq (ff,dname,2)
        end if
      else if (has(nbin1,bin1,dname)) then
        if (paralleldump) then
          call wr2dhdf_bin (ff,dname,1)
        else
          call wr2dhdf_bin_seq (ff,dname,1)
        end if
      end if

!
      return
      end subroutine wr2d
!***********************************************************************
      subroutine wr2dhdf (ff,dname)
      USE globals
      USE comm
      USE HDF5
      IMPLICIT NONE

      REAL(RTYPE) ::  ff(nx,nz)
      CHARACTER(LEN=*) :: dname
      INTEGER :: ierr


!-----------------------------------------------------------------------
!
! ****** Write the 2D data in array F into an HDF file.
!
!   


      if (ipltrtype.eq.1) then
        f1=real(ff(xs:xe,zs:ze))
      else  
        f2=dble(ff(xs:xe,zs:ze))
      end if

        
      call h5screate_simple_f(2,data_dims,dspace_id,ierr) 
      call h5screate_simple_f(2,buf_dims,memspace_id,ierr)

      !
      ! Create chunked dataset.
      !
      call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, ierr)

      call h5dcreate_f(file_id,dname,rtype_id,dspace_id, &
                       dset_id,ierr,plist_id) 

      call h5sclose_f(dspace_id, ierr)
      ! 
      ! Select hyperslab in the file.
      !

      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, ierr, &
                                 stride, buf_dims)
      !
      ! Create property list for collective dataset write
      !

      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,ierr)
      !
      ! Write the dataset collectively. 
      !
      if (ipltrtype.eq.1) then 
        call h5dwrite_f(dset_id, rtype_id, f1, data_dims, ierr, &
                     file_space_id = dspace_id, mem_space_id = memspace_id, &
                     xfer_prp = plist_id)
      else
        call h5dwrite_f(dset_id, rtype_id, f2, data_dims, ierr, &
                     file_space_id = dspace_id, mem_space_id = memspace_id, &
                     xfer_prp = plist_id)
      end if 
      call h5sclose_f(dspace_id,ierr)
      call h5sclose_f(memspace_id,ierr)
      call h5dclose_f(dset_id,ierr)
      call h5pclose_f(plist_id,ierr)

      return
      end subroutine wr2dhdf
!---------------------------------------------------------------------
      subroutine wr2dhdf_seq (ff, dname)
      USE globals
      USE comm
      USE HDF5
      USE ssub
      IMPLICIT NONE

      REAL(RTYPE) ::  ff(nx,nz)
      CHARACTER(LEN=*) :: dname


      INTEGER :: ibx, ibz, i, ox, oz
      INTEGER :: xs,xe,zs,ze,nnx,nnz
      INTEGER :: ierr
!-----------------------------------------------------------------------
!
! ****** Write the 2D data in array F into an HDF file.
!
!   
      


! gather data to PE0

      call gather(reshape(ff,(/nx*nz/)),f3,0,MPI_COMM_WORLD)

! write data by PE0

      if (proc0) then

! reconstruct the array from gathered data
 
        do ibz=0, nzproc-1
          do ibx=0, nxproc-1
            i=proc_id(ibx,ibz)*nx*nz
            f4=reshape(f3(i+1:i+nx*nz),(/nx,nz/))
            nnx=nx0
            nnz=nz0
            xs=3
            xe=nx-2
            zs=3
            ze=nz-2
            if (ibx.eq.0) then
              nnx=nnx+2
              xs=1 
              ox=0
            else
              ox=nx0*ibx+2
            end if

            if (ibx.eq.(nxproc-1)) then
              nnx=nnx+2
              xe=nx
            end if
  
            if (ibz.eq.0) then
              nnz=nnz+2
              zs=1
              oz=0
            else
              oz=nz0*ibz+2
            end if

            if (ibz.eq.(nzproc-1)) then
              nnz=nnz+2
              ze=nz
            end if


            if (ipltrtype.eq.1) then
              f1(ox+1:ox+nnx,oz+1:oz+nnz)=real(f4(xs:xe,zs:ze))
            else
              f2(ox+1:ox+nnx,oz+1:oz+nnz)=dble(f4(xs:xe,zs:ze))
            end if
          end do
        end do

! write data



        call h5screate_simple_f(2,data_dims(1),dspace_id,ierr)
        call h5dcreate_f(file_id,dname,rtype_id,dspace_id, &
                       dset_id,ierr)
        call h5screate_simple_f(2,buf_dims,memspace_id,ierr)
        call h5sselect_hyperslab_f(memspace_id,H5S_SELECT_SET_F,offset, &
                                 data_dims,ierr)

        if (ipltrtype.eq.1) then
          call h5dwrite_f(dset_id,rtype_id,f1,buf_dims,ierr, &
                      memspace_id,dspace_id)
        else
          call h5dwrite_f(dset_id,rtype_id,f2,buf_dims,ierr, &
                      memspace_id,dspace_id)
        end if

        call h5sclose_f(dspace_id,ierr)
        call h5sclose_f(memspace_id,ierr)
        call h5dclose_f(dset_id,ierr)

      end if 

  
      call barrier

      return
      end subroutine wr2dhdf_seq
!-----------------------------------------------------------------------
      subroutine wr2dhdf_bin (ff,dname,nb)
      USE globals
      USE comm
      USE HDF5
      !USE hdfwrite 
      IMPLICIT NONE

      REAL(RTYPE) ::  ff(nx,nz)
      CHARACTER(LEN=*) :: dname
      INTEGER :: nb
      INTEGER :: ierr
      REAL(RTYPE) :: fmin,fmax
      REAL(RTYPE) :: ofst, scle
      INTEGER(HID_T) :: itype, otype
      INTEGER :: order


!-----------------------------------------------------------------------
!
! ****** Write the 2D data in array F into an HDF file.
!
!   

! This is more elegant, but not working on Franklin & Hopper
! get the native integer type
      
!      call h5tcopy_f(H5T_NATIVE_INTEGER,itype,ierr)
!      call h5tcopy_f(H5T_NATIVE_INTEGER,otype,ierr)
!      if (nb.eq.1) then
!        call h5tset_size_f(otype,1,ierr)
!        if (huge(0_BYTE).eq.127) then
!          call h5tset_size_f(itype,1,ierr)
!        else if (huge(0_BYTE).eq.32767) then
!          call h5tset_size_f(itype,2,ierr)
!        else if (huge(0_BYTE).eq.2147483647) then 
!          call h5tset_size_f(itype,4,ierr)
!        end if 
!      else if (nb.eq.2) then
!        call h5tset_size_f(otype,2,ierr)
!        if (huge(0_WORD).eq.32767) then
!          call h5tset_size_f(itype,2,ierr)
!        else if (huge(0_WORD).eq.2147483647) then 
!          call h5tset_size_f(itype,4,ierr)
!        end if 
!      end if

! get the order of native integer type
      
      call h5tget_order_f(H5T_NATIVE_INTEGER,order,ierr)

      if (order.eq.0) then  ! Little Endian
        if (nb.eq.1) then
          otype=H5T_STD_I8LE   ! output type
          if (huge(0_IBYTE).eq.127) then
            itype=H5T_STD_I8LE       ! input type 8-bit
          else if (huge(0_IBYTE).eq.32767) then
            itype=H5T_STD_I16LE       ! input type 16-bit
          else if (huge(0_IBYTE).eq.2147483647) then
            itype=H5T_STD_I32LE       ! input type 32-bit
          end if
        else if (nb.eq.2) then
          otype=H5T_STD_I16LE   ! output type
          if (huge(0_IWORD).eq.32767) then
            itype=H5T_STD_I16LE       ! input type 16-bit
          else if (huge(0_IWORD).eq.2147483647) then
            itype=H5T_STD_I32LE       ! input type 32-bit
          end if
        end if
      else    ! Big Endian
        if (nb.eq.1) then
          otype=H5T_STD_I8BE   ! output type
          if (huge(0_IBYTE).eq.127) then
            itype=H5T_STD_I8BE       ! input type 8-bit
          else if (huge(0_IBYTE).eq.32767) then
            itype=H5T_STD_I16BE       ! input type 16-bit
          else if (huge(0_IBYTE).eq.2147483647) then
            itype=H5T_STD_I32BE       ! input type 32-bit
          end if
        else if (nb.eq.2) then
          otype=H5T_STD_I16BE   ! output type
          if (huge(0_IWORD).eq.32767) then
            itype=H5T_STD_I16BE       ! input type 16-bit
          else if (huge(0_IWORD).eq.2147483647) then
            itype=H5T_STD_I32BE       ! input type 32-bit
          end if
        end if
      end if
       
 
      fmin=minval(ff)
      fmax=maxval(ff)
      call min_allreduce(fmin)
      call max_allreduce(fmax)

      ofst=fmin
      if (nb.eq.2) then  ! two-byte binning
        scle=(fmax-fmin)/65535._RTYPE
      else
        scle=(fmax-fmin)/255._RTYPE
      end if
      if (scle.lt.1.e-10_RTYPE) then
        scle=1.e-10_RTYPE
      end if
      
      if (nb.eq.2) then 
        w2=nint((ff-ofst)/scle-32768._RTYPE)
        ofst=ofst+scle*32768._RTYPE
        w1=w2(xs:xe,zs:ze)
      else
        b2=nint((ff-ofst)/scle-128._RTYPE)
        ofst=ofst+scle*128._RTYPE
        b1=b2(xs:xe,zs:ze)
      end if



      call h5screate_simple_f(2,data_dims,dspace_id,ierr) 
      call h5screate_simple_f(2,buf_dims,memspace_id,ierr)
      !
      ! Create chunked dataset.
      !
      call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, ierr)

      call h5dcreate_f(file_id,dname,otype,dspace_id, &
                       dset_id,ierr,plist_id) 

      call h5sclose_f(dspace_id, ierr)
      ! 
      ! Select hyperslab in the file.
      !

      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, ierr, &
                                 stride, buf_dims)
      !
      ! Create property list for collective dataset write
      !

      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,ierr)
      !
      ! Write the dataset collectively. 
      !

      if (nb.eq.1) then
        call h5dwrite_f(dset_id,itype, b1, data_dims, ierr, &
                   file_space_id = dspace_id, mem_space_id = memspace_id, &
                   xfer_prp = plist_id)
      else if (nb.eq.2) then
        call h5dwrite_f(dset_id,itype, w1, data_dims, ierr, &
                   file_space_id = dspace_id, mem_space_id = memspace_id, &
                   xfer_prp = plist_id)
      end if


      call h5sclose_f(dspace_id,ierr)
      call h5sclose_f(memspace_id,ierr)
      call h5dclose_f(dset_id,ierr)
      call h5pclose_f(plist_id,ierr)

      if (proc0) then
! store offset and scale
        attr_id=attr_id+1
        attr_name(attr_id)=dname//'_offset'
        attr_value(attr_id)=ofst
        attr_id=attr_id+1
        attr_name(attr_id)=dname//'_scale'
        attr_value(attr_id)=scle
      end if
      return
      end subroutine wr2dhdf_bin
!-----------------------------------------------------------------------
      subroutine wr2dhdf_bin_seq (ff,dname,nb)
      USE globals
      USE comm
      USE HDF5
      !USE hdfwrite 
      USE ssub
      IMPLICIT NONE

      REAL(RTYPE) ::  ff(nx,nz)
      CHARACTER(LEN=*) :: dname
      INTEGER :: nb

      REAL(RTYPE) :: fmin,fmax
      REAL(RTYPE) :: ofst, scle
      INTEGER :: ibx, ibz, i, ox, oz
      INTEGER :: xs,xe,zs,ze,nnx,nnz
      INTEGER :: ierr
     
      INTEGER(HID_T) :: itype, otype 
      INTEGER :: order 


!-----------------------------------------------------------------------
!
! ****** Write the 2D data in array F into an HDF file.
!
!   


! This is more elegant, but not working on Franklin & Hopper
! get the native integer type
      
!      call h5tcopy_f(H5T_NATIVE_INTEGER,itype,ierr)
!      call h5tcopy_f(H5T_NATIVE_INTEGER,otype,ierr)
!      if (nb.eq.1) then
!        call h5tset_size_f(otype,1,ierr)
!        if (huge(0_BYTE).eq.127) then
!          call h5tset_size_f(itype,1,ierr)
!        else if (huge(0_BYTE).eq.32767) then
!          call h5tset_size_f(itype,2,ierr)
!        else if (huge(0_BYTE).eq.2147483647) then 
!          call h5tset_size_f(itype,4,ierr)
!        end if 
!      else if (nb.eq.2) then
!        call h5tset_size_f(otype,2,ierr)
!        if (huge(0_WORD).eq.32767) then
!          call h5tset_size_f(itype,2,ierr)
!        else if (huge(0_WORD).eq.2147483647) then 
!          call h5tset_size_f(itype,4,ierr)
!        end if 
!      end if

! get the order of native integer type
      
      call h5tget_order_f(H5T_NATIVE_INTEGER,order,ierr)

      if (order.eq.0) then  ! Little Endian
        if (nb.eq.1) then
          otype=H5T_STD_I8LE   ! output type
          if (huge(0_IBYTE).eq.127) then
            itype=H5T_STD_I8LE       ! input type 8-bit
          else if (huge(0_IBYTE).eq.32767) then
            itype=H5T_STD_I16LE       ! input type 16-bit
          else if (huge(0_IBYTE).eq.2147483647) then
            itype=H5T_STD_I32LE       ! input type 32-bit
          end if
        else if (nb.eq.2) then
          otype=H5T_STD_I16LE   ! output type
          if (huge(0_IWORD).eq.32767) then
            itype=H5T_STD_I16LE       ! input type 16-bit
          else if (huge(0_IWORD).eq.2147483647) then
            itype=H5T_STD_I32LE       ! input type 32-bit
          end if
        end if
      else    ! Big Endian
        if (nb.eq.1) then
          otype=H5T_STD_I8BE   ! output type
          if (huge(0_IBYTE).eq.127) then
            itype=H5T_STD_I8BE       ! input type 8-bit
          else if (huge(0_IBYTE).eq.32767) then
            itype=H5T_STD_I16BE       ! input type 16-bit
          else if (huge(0_IBYTE).eq.2147483647) then
            itype=H5T_STD_I32BE       ! input type 32-bit
          end if
        else if (nb.eq.2) then
          otype=H5T_STD_I16BE   ! output type
          if (huge(0_IWORD).eq.32767) then
            itype=H5T_STD_I16BE       ! input type 16-bit
          else if (huge(0_IWORD).eq.2147483647) then
            itype=H5T_STD_I32BE       ! input type 32-bit
          end if
        end if
      end if
       
      
      fmin=minval(ff)
      fmax=maxval(ff)
      call min_allreduce(fmin)
      call max_allreduce(fmax)

      ofst=fmin
      if (nb.eq.2) then  ! two-byte binning
        scle=(fmax-fmin)/65535._RTYPE
      else
        scle=(fmax-fmin)/255._RTYPE
      end if
      if (scle.lt.1.e-10_RTYPE) then
        scle=1.e-10_RTYPE
      end if

      if (nb.eq.2) then 
        w2=nint((ff-ofst)/scle-32768._RTYPE)
        ofst=ofst+scle*32768._RTYPE
      else
        b2=nint((ff-ofst)/scle-128._RTYPE)
        ofst=ofst+scle*128._RTYPE
      end if
       
! gather data to PE0

      if (nb.eq.2) then 
        call gather(reshape(w2,(/nx*nz/)),w3,0,MPI_COMM_WORLD)
      else
        call gather(reshape(b2,(/nx*nz/)),b3,0,MPI_COMM_WORLD)
      end if


! write data by PE0

      if (proc0) then

! reconstruct the array from gathered data
 
        do ibz=0, nzproc-1
          do ibx=0, nxproc-1
            i=proc_id(ibx,ibz)*nx*nz
            if (nb.eq.2) then 
              w2=reshape(w3(i+1:i+nx*nz),(/nx,nz/))
            else
              b2=reshape(b3(i+1:i+nx*nz),(/nx,nz/))
            end if
            nnx=nx0
            nnz=nz0
            xs=3
            xe=nx-2
            zs=3
            ze=nz-2
            if (ibx.eq.0) then
              nnx=nnx+2
              xs=1 
              ox=0
            else
              ox=nx0*ibx+2
            end if

            if (ibx.eq.(nxproc-1)) then
              nnx=nnx+2
              xe=nx
            end if
  
            if (ibz.eq.0) then
              nnz=nnz+2
              zs=1
              oz=0
            else
              oz=nz0*ibz+2
            end if

            if (ibz.eq.(nzproc-1)) then
              nnz=nnz+2
              ze=nz
            end if
            if (nb.eq.2) then 
              w1(ox+1:ox+nnx,oz+1:oz+nnz)=w2(xs:xe,zs:ze)
            else
              b1(ox+1:ox+nnx,oz+1:oz+nnz)=b2(xs:xe,zs:ze)
            end if
            
          end do
        end do


! write data


        call h5screate_simple_f(2,data_dims(1),dspace_id,ierr)
        call h5dcreate_f(file_id,dname,otype,dspace_id, &
                       dset_id,ierr)
        call h5screate_simple_f(2,buf_dims,memspace_id,ierr)
        call h5sselect_hyperslab_f(memspace_id,H5S_SELECT_SET_F,offset, &
                                 data_dims,ierr)

        if (nb.eq.1) then
          call h5dwrite_f(dset_id,itype,b1,buf_dims,ierr, &
                      memspace_id,dspace_id)
        else if (nb.eq.2) then
          call h5dwrite_f(dset_id,itype,w1,buf_dims,ierr, &
                      memspace_id,dspace_id)
        end if


        call h5sclose_f(dspace_id,ierr)
        call h5sclose_f(memspace_id,ierr)
        call h5dclose_f(dset_id,ierr)

! store offset and scale
        attr_id=attr_id+1
        attr_name(attr_id)=dname//'_offset'
        attr_value(attr_id)=ofst
        attr_id=attr_id+1
        attr_name(attr_id)=dname//'_scale'
        attr_value(attr_id)=scle
      end if
  
      call barrier

      return
      end subroutine wr2dhdf_bin_seq
!-----------------------------------------------------------------------
      end subroutine pdump

!#######################################################################
      subroutine dumphist
!
!-----------------------------------------------------------------------
!
! ****** Dump time history data to a file.
!
!-----------------------------------------------------------------------
!
      USE globals
      USE comm
      USE util, ONLY : lenstr
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! ****** History file name and sequence number.
!
      CHARACTER(LEN=72) :: hstfile
      INTEGER, SAVE :: iseq=0
      INTEGER :: lf, nhist
!
!-----------------------------------------------------------------------
!
      if (ihist.le.1) return
!
! ****** Number of history intervals.
!
      nhist=ihist
!
! ****** Append a sequence number to the history file name root.
!
      iseq=iseq+1
      if (iseq.gt.99) iseq=iseq-99
      hstfile=hstroot
      lf=lenstr(hstfile)
      write (hstfile(lf+1:lf+2),'(i2.2)') iseq
      lf=lenstr(hstfile)
!
      hstfile(lf+1:)='.hdf'
      lf=lenstr(hstfile)
      call dumphisthdf(hstfile(1:lf),nhist)

!
      if (proc0) then
        write (9,*)
        write (9,*) '### Comment from DUMPHIST:'
        write (9,*) 'Wrote time histories to file: ',hstfile(1:lf)
        write (9,*) 'NTIME = ',ntime
        write (9,*) 'TIME = ',time
        write (9,*) 'NHIST = ',nhist
      end if 
!
      return
      end subroutine dumphist

!#######################################################################
      subroutine dumphisthdf(hstfile, nhist)
!history dump in HDF format
      USE comm
      USE prob
      USE HDF5
      USE util, ONLY : lenstr
      IMPLICIT NONE
      CHARACTER(LEN=*) :: hstfile
      INTEGER :: nhist
      INTEGER :: ierr
      INTEGER(HID_T) :: file_id
      INTEGER(HID_T) :: dset_id, dspace_id, memspace_id, rtype_id
      INTEGER(HSIZE_T) :: data_dims(2), buf_dims(2), offset(2)=(/0,0/)

!
! Initialize FORTRAN interface
!
      call h5open_f(ierr)
!
      if (RTYPE.eq.DP) then
        rtype_id=H5T_NATIVE_DOUBLE
      else
        rtype_id=H5T_NATIVE_REAL
      end if  

! create file
!
      call h5fcreate_f(hstfile,H5F_ACC_TRUNC_F,file_id,ierr)
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in DUMPHIST:'
        write (*,*) '### Error while opening a history file.'
        write (*,*) 'File name: ',hstfile
        write (*,*) 'IERR  = ',ierr
        write (*,*) 'NTIME = ',ntime
        write (*,*) 'TIME = ',time
        call wrrsfile (rsdir(1:lenstr(rsdir))//rsofile)
        call abort_mp
        call exit(1)

      end if

!
! write history data
!
      data_dims(1)=nhist
      data_dims(2)=nhq
      buf_dims(1)=nhistmax  ! dimension of the buffer
      buf_dims(2)=nhq
      call h5screate_simple_f(2,data_dims,dspace_id,ierr)     
      call h5dcreate_f(file_id,'thist',rtype_id,dspace_id, &
                       dset_id,ierr)
      call h5screate_simple_f(2,buf_dims,memspace_id,ierr)
! Select the subarray; offset=(0,0) means no offset. Notice this is more like 
! the C convention. The hdf5 manual uses "start" instead of "offset", which is misleading.
! A fortran array can start from any index, while a C array always starts from index 0.
! "start" is appropriate for C, while "offset" is appropriate for fortran.  (Yi-Min Huang)  
      call h5sselect_hyperslab_f(memspace_id,H5S_SELECT_SET_F,offset, &
                                 data_dims,ierr) 
      call h5dwrite_f(dset_id,rtype_id,thist,buf_dims,ierr,&
                      memspace_id,dspace_id)
      call h5sclose_f(dspace_id,ierr)
      call h5sclose_f(memspace_id,ierr)
      call h5dclose_f(dset_id,ierr)
!
! close file
!
      call h5fclose_f(file_id,ierr)
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in DUMPHIST:'
        write (*,*) '### Error while closing a history file.'
        write (*,*) 'File name: ',hstfile
        write (*,*) 'IERR  = ',ierr
        write (*,*) 'NTIME = ',ntime
        write (*,*) 'TIME = ',time
        call wrrsfile (rsdir(1:lenstr(rsdir))//rsofile)
        call abort_mp
        call exit(1)
      end if
!
! Close FORTRAN interface.
!
      call h5close_f(ierr)
      end subroutine dumphisthdf
!#########################################################################
      subroutine diags
!
!-----------------------------------------------------------------------
!
! ****** Collect time diagnostics.
!
!-----------------------------------------------------------------------
!
      USE globals
      USE comm
      USE prob, ONLY : tdiagcol, diagnostics
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!

      call diagnostics

      if (ifhststp.ne.0) call tdiagcol
!
! ****** Plot spatial diagnostics if on diagnostic interval.
!
      if (ifplxstp.ne.0) call pdump
!
      return
      end subroutine diags

!#######################################################################
      subroutine histchek
!
!-----------------------------------------------------------------------
!
! ****** Set switch for collection of time histories.
! ****** Also, set switch for plotting of spatial diagnostics.
!
!-----------------------------------------------------------------------
!
      USE globals
      USE comm
      USE prob
      IMPLICIT NONE
      REAL(RTYPE) :: trem
!
!-----------------------------------------------------------------------
!
! ****** Collect time-histories at intervals of IHISTINT timesteps,
! ****** or every time interval THISTINT.
!
! ****** Disable collection at timestep intervals if IHISTINT.le.0.
!

      ifhststp=0

      if (ihistint>0) then
        ihist_res=ihist_res-1
        if (ihist_res<=0) then
          ifhststp=1
          ihist_res=ihistint
        end if
      end if
      
!
! ****** Disable collection at time intervals if THISTINT.le.0.
!

      if (thistint>zero) then
        thist_res=thist_res-dt
        if (thist_res<=0.49*dt) then
          ifhststp=1
          thist_res=thist_res+thistint
        end if
      end if
      
!
! ****** Time history diagnostics are collected whenever IFHSTSTP.ne.0.
!
! ****** Check for overflowing history buffer.
!
      if (ihist.eq.nhistmax) then
        if (iproc.eq.0) call dumphist     ! PE0 dump history
        ihist=0
      end if
!
! ****** Increment the time-history index.
!
      if (ifhststp.ne.0) ihist=ihist+1
!
! ****** Set the switch for spatial diagnostics.
!
! ****** Plot x-diagnostics at intervals of IPLTXINT time steps,
! ****** and every TPLTXINT time interval.
!
! ****** Disable plotting at time step intervals if IPLTXINT.le.0.
!

      ifplxstp=0

      if (ipltxint>0) then
        ipltx_res=ipltx_res-1
        if (ipltx_res<=0) then
          ifplxstp=1
          ipltx_res=ipltxint
        end if
      end if
      
!
! ****** Disable plotting at time intervals if TPLTXINT.le.0.
!

      if (tpltxint>zero) then
        tpltx_res=tpltx_res-dt
        if (tpltx_res<=0.49*dt) then
          ifplxstp=1
          tpltx_res=tpltx_res+tpltxint
        end if
      end if
      
!
      return
      end subroutine histchek
!#######################################################################

      END MODULE io
