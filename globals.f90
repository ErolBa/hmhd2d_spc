MODULE globals
!
!-----------------------------------------------------------------------
!
      INTEGER, PARAMETER :: SP=kind(1.0),DP=kind(1.0D0)     ! Single and double precision
      INTEGER, PARAMETER :: RTYPE=DP   ! Accuracy of variables
      REAL(RTYPE), PARAMETER :: eps=tiny(1._RTYPE) ! smallest positive float
! Integer type
      INTEGER, PARAMETER :: IBYTE=selected_int_kind(2)    ! at least 1 byte
      INTEGER, PARAMETER :: IWORD=selected_int_kind(4)    ! at least 2 bytes
      INTEGER, PARAMETER :: ILWORD=selected_int_kind(8)   ! at least 4 bytes
! Resolutions

      INTEGER, PARAMETER :: nx0=30, &  !d=26! no. of real grids in each node 
                            nz0=38    !d=34

      INTEGER, PARAMETER :: nx=nx0+4,  &  ! resolution in x in each node
                            nz=nz0+4      ! resolution in z in each node  


      INTEGER, PARAMETER :: nxblock=1,   &    ! Number of blocks in x direction
                            nzblock=1         ! Number of blocks in z direction


      INTEGER, PARAMETER :: block_dim_x=3,  &  ! The size of each block 
                            block_dim_z=3      ! This may be utilize for better
                                               ! data proximity

      INTEGER, PARAMETER :: nxproc=nxblock*block_dim_x, &  ! no. procs in x
                            nzproc=nzblock*block_dim_z     ! no. procs in z

      INTEGER, PARAMETER :: blocksize=block_dim_x*block_dim_z 

      INTEGER, PARAMETER :: nxtot=nx0*nxproc+4, &   ! no. of total grid points in x
                            nztot=nz0*nzproc+4      ! no. of total grid points in z

!
!-----------------------------------------------------------------------
!

END MODULE globals
