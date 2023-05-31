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
! Resolution should be as small as possible, but if it's too small results will be strange and will crash the python analysis tools
	INTEGER, PARAMETER :: nx0=100, nz0=80

      INTEGER, PARAMETER :: nx=nx0+4,  &  ! resolution in x in each node
                            nz=nz0+4      ! resolution in z in each node  


      INTEGER, PARAMETER :: nxblock=1,   &    ! Number of blocks in x direction
                            nzblock=1         ! Number of blocks in z direction

      ! The size of each block may be utilize for better data proximity

	INTEGER, PARAMETER :: block_dim_x=2, block_dim_z=4

      INTEGER, PARAMETER :: nxproc=nxblock*block_dim_x, &  ! no. procs in x
                            nzproc=nzblock*block_dim_z     ! no. procs in z

      INTEGER, PARAMETER :: blocksize=block_dim_x*block_dim_z 

      INTEGER, PARAMETER :: nxtot=nx0*nxproc+4, &   ! no. of total grid points in x
                            nztot=nz0*nzproc+4      ! no. of total grid points in z

!
!-----------------------------------------------------------------------
!

END MODULE globals
