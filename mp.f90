
module mp
  implicit none
  include 'mpif.h'
  private

  public :: init_mp, init_mp_thread, check_mp, finish_mp, abort_mp
  public :: broadcast, sum_reduce, sum_allreduce
  public :: max_reduce, max_allreduce
  public :: min_reduce, min_allreduce
  public :: nproc, iproc, proc0, npe
  public :: send, receive
  public :: isend, ireceive
  public :: rsend
  public :: irsend
  public :: barrier
  public :: allgather
  public :: gather, scatter
  public :: gettime
  public :: or_allreduce    

  integer :: nproc, iproc
  logical :: proc0=.TRUE.

  integer :: npe 


  interface broadcast
     module procedure broadcast_integer, broadcast_integer_array
     module procedure broadcast_real,    broadcast_real_array
     module procedure broadcast_double,  broadcast_double_array
     module procedure broadcast_complex, broadcast_complex_array
     module procedure broadcast_double_complex, broadcast_double_complex_array
     module procedure broadcast_logical, broadcast_logical_array
     module procedure broadcast_character, broadcast_character_array
     module procedure bcastfrom_integer, bcastfrom_integer_array
     module procedure bcastfrom_real,    bcastfrom_real_array
     module procedure bcastfrom_double,  bcastfrom_double_array
     module procedure bcastfrom_complex, bcastfrom_complex_array
     module procedure bcastfrom_double_complex, bcastfrom_double_complex_array
     module procedure bcastfrom_character, bcastfrom_character_array
     module procedure bcastfrom_logical, bcastfrom_logical_array
     module procedure bcastfromto_integer, bcastfromto_integer_array
     module procedure bcastfromto_real,    bcastfromto_real_array
     module procedure bcastfromto_double,  bcastfromto_double_array
     module procedure bcastfromto_complex, bcastfromto_complex_array
     module procedure bcastfromto_double_complex, bcastfromto_double_complex_array
     module procedure bcastfromto_character, bcastfromto_character_array
     module procedure bcastfromto_logical, bcastfromto_logical_array
  end interface

  interface sum_reduce
     module procedure sum_reduce_integer, sum_reduce_integer_array
     module procedure sum_reduce_real,    sum_reduce_real_array
     module procedure sum_reduce_double,  sum_reduce_double_array
     module procedure sum_reduce_complex, sum_reduce_complex_array
     module procedure sum_reduce_double_complex, sum_reduce_double_complex_array
  end interface

  interface sum_allreduce
     module procedure sum_allreduce_integer, sum_allreduce_integer_array
     module procedure sum_allreduce_real,    sum_allreduce_real_array
     module procedure sum_allreduce_double,  sum_allreduce_double_array
     module procedure sum_allreduce_complex, sum_allreduce_complex_array
     module procedure sum_allreduce_double_complex, sum_allreduce_double_complex_array
  end interface

  interface max_reduce
     module procedure max_reduce_integer, max_reduce_integer_array
     module procedure max_reduce_real,    max_reduce_real_array
     module procedure max_reduce_double,  max_reduce_double_array
  end interface

  interface max_allreduce
     module procedure max_allreduce_integer, max_allreduce_integer_array
     module procedure max_allreduce_real,    max_allreduce_real_array
     module procedure max_allreduce_double,  max_allreduce_double_array
  end interface

  interface min_reduce
     module procedure min_reduce_integer, min_reduce_integer_array
     module procedure min_reduce_real,    min_reduce_real_array
     module procedure min_reduce_double,  min_reduce_double_array
  end interface

  interface min_allreduce
     module procedure min_allreduce_integer, min_allreduce_integer_array
     module procedure min_allreduce_real,    min_allreduce_real_array
     module procedure min_allreduce_double,    min_allreduce_double_array
  end interface

  interface send
     module procedure send_integer, send_integer_array
     module procedure send_real,    send_real_array
     module procedure send_double,    send_double_array
     module procedure send_complex, send_complex_array
     module procedure send_double_complex, send_double_complex_array
     module procedure send_logical, send_logical_array
  end interface

  interface receive
     module procedure receive_integer, receive_integer_array
     module procedure receive_real,    receive_real_array
     module procedure receive_double,    receive_double_array
     module procedure receive_complex, receive_complex_array
     module procedure receive_double_complex, receive_double_complex_array
     module procedure receive_logical, receive_logical_array
  end interface

  interface isend
     module procedure isend_integer, isend_integer_array
     module procedure isend_real,    isend_real_array
     module procedure isend_double,    isend_double_array
     module procedure isend_complex, isend_complex_array
     module procedure isend_double_complex, isend_double_complex_array
     module procedure isend_logical, isend_logical_array
  end interface

  interface ireceive
     module procedure ireceive_integer, ireceive_integer_array
     module procedure ireceive_real,    ireceive_real_array
     module procedure ireceive_double,    ireceive_double_array
     module procedure ireceive_complex, ireceive_complex_array
     module procedure ireceive_double_complex, ireceive_double_complex_array
     module procedure ireceive_logical, ireceive_logical_array
  end interface

  interface rsend
     module procedure rsend_integer, rsend_integer_array
     module procedure rsend_real,    rsend_real_array
     module procedure rsend_double,    rsend_double_array
     module procedure rsend_complex, rsend_complex_array
     module procedure rsend_double_complex, rsend_double_complex_array
     module procedure rsend_logical, rsend_logical_array
  end interface

  interface irsend
     module procedure irsend_integer, irsend_integer_array
     module procedure irsend_real,    irsend_real_array
     module procedure irsend_double,    irsend_double_array
     module procedure irsend_complex, irsend_complex_array
     module procedure irsend_double_complex, irsend_double_complex_array
     module procedure irsend_logical, irsend_logical_array
  end interface

  interface allgather
     module procedure allgather_real
     module procedure allgather_double
  end interface

  interface gather
     module procedure gather_real_array
     module procedure gather_double_array
     module procedure gather_integer_array
     module procedure gather_integer1_array
     module procedure gather_integer2_array
     module procedure gather_char_array
  end interface

  interface scatter
     module procedure scatter_real_array
     module procedure scatter_double_array
  end interface


contains

  function gettime()
     implicit none
     real :: gettime,oldtime,newtime
     save oldtime

     newtime = MPI_WTIME()
     gettime = newtime - oldtime
     oldtime = newtime

   end function gettime


  subroutine scatter_real_array(ii_send,ii_recv,root,comm)
    implicit none
    real(4), dimension(:), intent(in) :: ii_send
    real(4), dimension(:), intent(out) :: ii_recv
    integer :: root,comm
    integer :: ierror

    call mpi_scatter(ii_send,size(ii_recv),MPI_REAL,ii_recv,size(ii_recv), &
                   MPI_REAL,root,comm,ierror)

  end subroutine scatter_real_array

  subroutine scatter_double_array(ii_send,ii_recv,root,comm)
    implicit none
    real(8), dimension(:), intent(in) :: ii_send
    real(8), dimension(:), intent(out) :: ii_recv
    integer :: root,comm
    integer :: ierror

    call mpi_scatter(ii_send,size(ii_recv),MPI_DOUBLE_PRECISION,ii_recv,size(ii_recv), &
                   MPI_DOUBLE_PRECISION,root,comm,ierror)

  end subroutine scatter_double_array


  subroutine gather_real_array(ii,ii_gather,root,comm)
    implicit none
    real(4), dimension(:), intent(in) :: ii
    real(4), dimension(:), intent(out) :: ii_gather
    integer :: root, comm
    integer :: ierror
    
    call mpi_gather(ii,size(ii),MPI_REAL,ii_gather,size(ii), &
                   MPI_REAL,root,comm,ierror)

  end subroutine gather_real_array

  subroutine gather_double_array(ii,ii_gather,root,comm)
    implicit none
    real(8), dimension(:), intent(in) :: ii
    real(8), dimension(:), intent(out) :: ii_gather
    integer :: root, comm
    integer :: ierror
    call mpi_gather(ii,size(ii),MPI_DOUBLE_PRECISION,ii_gather,size(ii), &
                   MPI_DOUBLE_PRECISION,root,comm,ierror)
  end subroutine gather_double_array

  subroutine gather_integer_array(ii,ii_gather,root,comm)
    implicit none
    integer, dimension(:), intent(in) :: ii
    integer, dimension(:), intent(out) :: ii_gather
    integer :: root, comm
    integer :: ierror
     
    call mpi_gather(ii,size(ii),MPI_INTEGER,ii_gather,size(ii), &
                   MPI_INTEGER,root,comm,ierror)
  end subroutine gather_integer_array

  subroutine gather_integer1_array(ii,ii_gather,root,comm)
    implicit none
    integer*1, dimension(:), intent(in) :: ii
    integer*1, dimension(:), intent(out) :: ii_gather
    integer :: root, comm
    integer :: ierror
    
    call mpi_gather(ii,size(ii),MPI_INTEGER1,ii_gather,size(ii), &
                   MPI_INTEGER1,root,comm,ierror)
  end subroutine gather_integer1_array

  subroutine gather_integer2_array(ii,ii_gather,root,comm)
    implicit none
    integer*2, dimension(:), intent(in) :: ii
    integer*2, dimension(:), intent(out) :: ii_gather
    integer :: root, comm
    integer :: ierror
     
    call mpi_gather(ii,size(ii),MPI_INTEGER2,ii_gather,size(ii), &
                   MPI_INTEGER2,root,comm,ierror)
  end subroutine gather_integer2_array

  subroutine gather_char_array(ii,ii_gather,root,comm)
    implicit none
    character(len=1), dimension(:), intent(in) :: ii
    character(len=1), dimension(:), intent(out) :: ii_gather
    integer :: root, comm
    integer :: ierror

    call mpi_gather(ii,size(ii),MPI_CHARACTER,ii_gather,size(ii), &
                   MPI_CHARACTER,root,comm,ierror)

  end subroutine gather_char_array

  subroutine allgather_real(ii,ii_gather)
    implicit none
    real(4), intent(in) :: ii
    real(4), dimension(:), intent(out) :: ii_gather
    integer :: ierror

    call mpi_allgather(ii,1,MPI_REAL,ii_gather,size(ii_gather), &
                   MPI_REAL,MPI_COMM_WORLD,ierror)

  end subroutine allgather_real

  subroutine allgather_double(ii,ii_gather)
    implicit none
    real(8), intent(in) :: ii
    real(8), dimension(:), intent(out) :: ii_gather
    integer :: ierror

    call mpi_allgather(ii,1,MPI_DOUBLE_PRECISION,ii_gather,size(ii_gather), &
                   MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)

  end subroutine allgather_double


  subroutine init_mp
    implicit none
    integer :: ierror, rank, temp

    call mpi_init (ierror)
    call mpi_comm_size (mpi_comm_world, nproc, ierror)
    call mpi_comm_rank (mpi_comm_world, iproc, ierror)
    proc0 = iproc == 0
    npe = nproc
    temp = gettime()    !!! Initializes gettime

  end subroutine init_mp

  subroutine init_mp_thread
    implicit none
    integer :: ierror, rank, temp, request, provided

    request=MPI_THREAD_MULTIPLE
    call mpi_init_thread (request,provided,ierror)
    if (.not.request.eq.provided) then
      print*,'Insufficient MPI thread safety to run Hybrid Mode'
      print*,'Request:  ',request
      print*,'Provided: ',provided
      call finish_mp
      Stop 300
    end if
    call mpi_comm_size (mpi_comm_world, nproc, ierror)
    call mpi_comm_rank (mpi_comm_world, iproc, ierror)
    proc0 = iproc == 0
    npe = nproc
    temp = gettime()    !!! Initializes gettime

  end subroutine init_mp_thread



  subroutine check_mp(nxproc,nyproc,nzproc)
     implicit none
     integer, intent(in) :: nxproc,nyproc,nzproc
 
     if (nproc /= nxproc*nyproc*nzproc) then
       write(*,*) 'The number of procs ', nproc
       write(*,*) 'is different from NXPROC*NYPROC*NZPROC ', nxproc*nyproc*nzproc
       call finish_mp
       stop 342
     end if     
  end subroutine check_mp

  subroutine finish_mp
    implicit none
    integer :: ierror

    call mpi_finalize (ierror)
  end subroutine finish_mp

  subroutine abort_mp
    implicit none
    integer :: ierror

    call mpi_abort (MPI_COMM_WORLD,ierror)
  end subroutine abort_mp 

  subroutine broadcast_character (char)
    implicit none
    character(*), intent (in out) :: char
    integer :: ierror
    call mpi_bcast (char, len(char), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_character

  subroutine broadcast_character_array (char)
    implicit none
    character(*), dimension(:), intent (in out) :: char
    integer :: ierror
    call mpi_bcast (char, len(char)*size(char), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_character_array


  subroutine broadcast_integer (i)
    implicit none
    integer, intent (in out) :: i
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_integer

  subroutine broadcast_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_integer_array

  subroutine broadcast_real (x)
    implicit none
    real(4), intent (in out) :: x
    integer :: ierror
    call mpi_bcast (x, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_real

  subroutine broadcast_real_array (x)
    implicit none
    real(4), dimension (:), intent (in out) :: x
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_REAL, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_real_array

  subroutine broadcast_double (x)
    implicit none
    real(8), intent (in out) :: x
    integer :: ierror
    call mpi_bcast (x, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_double

  subroutine broadcast_double_array (x)
    implicit none
    real(8), dimension (:), intent (in out) :: x
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_double_array


  subroutine broadcast_complex (z)
    implicit none
    complex(4), intent (in out) :: z
    integer :: ierror
    call mpi_bcast (z, 1, MPI_COMPLEX, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_complex

  subroutine broadcast_complex_array (z)
    implicit none
    complex(4), dimension (:), intent (in out) :: z
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_COMPLEX, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_complex_array

  subroutine broadcast_double_complex (z)
    implicit none
    complex(8), intent (in out) :: z
    integer :: ierror
    call mpi_bcast (z, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_double_complex

  subroutine broadcast_double_complex_array (z)
    implicit none
    complex(8), dimension (:), intent (in out) :: z
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_double_complex_array

  subroutine broadcast_logical (f)
    implicit none
    logical, intent (in out) :: f
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_logical

  subroutine broadcast_logical_array (f)
    implicit none
    logical, dimension (:), intent (in out) :: f
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_logical_array

!
  subroutine bcastfrom_character (c, src)
    implicit none
    character(*), intent (in out) :: c
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (c, len(c), MPI_CHARACTER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_character

  subroutine bcastfrom_character_array (c, src)
    implicit none
    character(*), dimension(:) , intent (in out) :: c
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (c, len(c)*size(c), MPI_CHARACTER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_character_array

  subroutine bcastfrom_integer (i, src)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_integer

  subroutine bcastfrom_integer_array (i, src)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_integer_array

  subroutine bcastfrom_real (x, src)
    implicit none
    real(4), intent (in out) :: x
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (x, 1, MPI_REAL, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_real

  subroutine bcastfrom_real_array (x, src)
    implicit none
    real(4), dimension (:), intent (in out) :: x
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_REAL, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_real_array

  subroutine bcastfrom_double (x, src)
    implicit none
    real(8), intent (in out) :: x
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (x, 1, MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_double

  subroutine bcastfrom_double_array (x, src)
    implicit none
    real(8), dimension (:), intent (in out) :: x
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_double_array


  subroutine bcastfrom_complex (z, src)
    implicit none
    complex(4), intent (in out) :: z
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (z, 1, MPI_COMPLEX, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_complex

  subroutine bcastfrom_complex_array (z, src)
    implicit none
    complex(4), dimension (:), intent (in out) :: z
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_COMPLEX, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_complex_array

  subroutine bcastfrom_double_complex (z, src)
    implicit none
    complex(8), intent (in out) :: z
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (z, 1, MPI_DOUBLE_COMPLEX, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_double_complex

  subroutine bcastfrom_double_complex_array (z, src)
    implicit none
    complex(8), dimension (:), intent (in out) :: z
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_double_complex_array

  subroutine bcastfrom_logical (f, src)
    implicit none
    logical, intent (in out) :: f
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_logical

  subroutine bcastfrom_logical_array (f, src)
    implicit none
    logical, dimension (:), intent (in out) :: f
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_logical_array

!

  subroutine bcastfromto_character (c, src, comm)
    implicit none
    character(*), intent (in out) :: c
    integer, intent (in) :: src, comm
    integer :: ierror
    call mpi_bcast (c, len(c), MPI_CHARACTER, src, comm, ierror)
  end subroutine bcastfromto_character

  subroutine bcastfromto_character_array (c, src, comm)
    implicit none
    character(*), dimension(:) , intent (in out) :: c
    integer, intent (in) :: src, comm
    integer :: ierror
    call mpi_bcast (c, len(c)*size(c), MPI_CHARACTER, src, comm, ierror)
  end subroutine bcastfromto_character_array

  subroutine bcastfromto_integer (i, src, comm)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: src, comm
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, src, comm, ierror)
  end subroutine bcastfromto_integer

  subroutine bcastfromto_integer_array (i, src, comm)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: src, comm
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, src, comm, ierror)
  end subroutine bcastfromto_integer_array

  subroutine bcastfromto_real (x, src, comm)
    implicit none
    real(4), intent (in out) :: x
    integer, intent (in) :: src, comm
    integer :: ierror
    call mpi_bcast (x, 1, MPI_REAL, src, comm, ierror)
  end subroutine bcastfromto_real

  subroutine bcastfromto_real_array (x, src, comm)
    implicit none
    real(4), dimension (:), intent (in out) :: x
    integer, intent (in) :: src, comm
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_REAL, src, comm, ierror)
  end subroutine bcastfromto_real_array

  subroutine bcastfromto_double (x, src, comm)
    implicit none
    real(8), intent (in out) :: x
    integer, intent (in) :: src, comm
    integer :: ierror
    call mpi_bcast (x, 1, MPI_DOUBLE_PRECISION, src, comm, ierror)
  end subroutine bcastfromto_double

  subroutine bcastfromto_double_array (x, src, comm)
    implicit none
    real(8), dimension (:), intent (in out) :: x
    integer, intent (in) :: src, comm
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, src, comm, ierror)
  end subroutine bcastfromto_double_array


  subroutine bcastfromto_complex (z, src, comm)
    implicit none
    complex(4), intent (in out) :: z
    integer, intent (in) :: src, comm
    integer :: ierror
    call mpi_bcast (z, 1, MPI_COMPLEX, src, comm, ierror)
  end subroutine bcastfromto_complex

  subroutine bcastfromto_complex_array (z, src, comm)
    implicit none
    complex(4), dimension (:), intent (in out) :: z
    integer, intent (in) :: src, comm
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_COMPLEX, src, comm, ierror)
  end subroutine bcastfromto_complex_array

  subroutine bcastfromto_double_complex (z, src, comm)
    implicit none
    complex(8), intent (in out) :: z
    integer, intent (in) :: src, comm
    integer :: ierror
    call mpi_bcast (z, 1, MPI_DOUBLE_COMPLEX, src, comm, ierror)
  end subroutine bcastfromto_double_complex

  subroutine bcastfromto_double_complex_array (z, src, comm)
    implicit none
    complex(8), dimension (:), intent (in out) :: z
    integer, intent (in) :: src, comm
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, src, comm, ierror)
  end subroutine bcastfromto_double_complex_array

  subroutine bcastfromto_logical (f, src, comm)
    implicit none
    logical, intent (in out) :: f
    integer, intent (in) :: src, comm
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, src, comm, ierror)
  end subroutine bcastfromto_logical

  subroutine bcastfromto_logical_array (f, src, comm)
    implicit none
    logical, dimension (:), intent (in out) :: f
    integer, intent (in) :: src, comm
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, src, comm, ierror)
  end subroutine bcastfromto_logical_array



!
  subroutine sum_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_integer

  subroutine sum_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_integer_array

  subroutine sum_reduce_real (a, dest)
    implicit none
    real(4), intent (in out) :: a
    integer, intent (in) :: dest
    real(4) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_REAL, MPI_SUM, dest, MPI_COMM_WORLD, ierror)

  end subroutine sum_reduce_real

  subroutine sum_reduce_real_array (a, dest)
    implicit none
    real(4), dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    real(4), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_REAL, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_real_array

  subroutine sum_reduce_double (a, dest)
    implicit none
    real(8), intent (in out) :: a
    integer, intent (in) :: dest
    real(8) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)

  end subroutine sum_reduce_double

  subroutine sum_reduce_double_array (a, dest)
    implicit none
    real(8), dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    real(8), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_double_array

  subroutine sum_reduce_complex (z, dest)
    implicit none
    complex(4), intent (in out) :: z
    integer, intent (in) :: dest
    complex(4) :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, 1, MPI_COMPLEX, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_complex

  subroutine sum_reduce_complex_array (z, dest)
    implicit none
    complex(4), dimension (:), intent (in out) :: z
    integer, intent (in) :: dest
    complex(4), dimension (size(z)) :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, size(z), MPI_COMPLEX, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_complex_array

  subroutine sum_reduce_double_complex (z, dest)
    implicit none
    complex(8), intent (in out) :: z
    integer, intent (in) :: dest
    complex(8) :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_double_complex

  subroutine sum_reduce_double_complex_array (z, dest)
    implicit none
    complex(8), dimension (:), intent (in out) :: z
    integer, intent (in) :: dest
    complex(8), dimension (size(z)) :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, size(z), MPI_DOUBLE_COMPLEX, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_double_complex_array


  subroutine sum_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_integer

  subroutine sum_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_integer_array

  subroutine sum_allreduce_real (a)
    implicit none
    real(4), intent (in out) :: a
    real(4) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_real

  subroutine sum_allreduce_real_array (a)
    implicit none
    real(4), dimension (:), intent (in out) :: a
    real(4), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_real_array

  subroutine sum_allreduce_double (a)
    implicit none
    real(8), intent (in out) :: a
    real(8) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_double

  subroutine sum_allreduce_double_array (a)
    implicit none
    real(8), dimension (:), intent (in out) :: a
    real(8), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_double_array

  subroutine sum_allreduce_complex (z)
    implicit none
    complex(4), intent (in out) :: z
    complex(4) :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, 1, MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_complex

  subroutine sum_allreduce_complex_array (z)
    implicit none
    complex(4), dimension (:), intent (in out) :: z
    complex(4), dimension (size(z)) :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, size(z), MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_complex_array

  subroutine sum_allreduce_double_complex (z)
    implicit none
    complex(8), intent (in out) :: z
    complex(8) :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_double_complex

  subroutine sum_allreduce_double_complex_array (z)
    implicit none
    complex(8), dimension (:), intent (in out) :: z
    complex(8), dimension (size(z)) :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, size(z), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_double_complex_array


  subroutine barrier
    implicit none
    integer :: ierror
    call mpi_barrier (MPI_COMM_WORLD, ierror)
  end subroutine barrier

  subroutine send_integer (i, dest, tag)
    implicit none
    integer, intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, 1, MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_integer

  subroutine send_integer_array (i, dest, tag)
    implicit none
    integer, dimension (:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, size(i), MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_integer_array

  subroutine send_real (a, dest, tag)
    implicit none
    real(4), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, 1, MPI_REAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_real

  subroutine send_real_array (a, dest, tag)
    implicit none
    real(4), dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, size(a), MPI_REAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_real_array

  subroutine send_double (a, dest, tag)
    implicit none
    real(8), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, 1, MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_double

  subroutine send_double_array (a, dest, tag)
    implicit none
    real(8), dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, size(a), MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_double_array


  subroutine send_complex (z, dest, tag)
    implicit none
    complex(4), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, 1, MPI_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_complex

  subroutine send_complex_array (z, dest, tag)
    implicit none
    complex(4), dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, size(z), MPI_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_complex_array

  subroutine send_double_complex (z, dest, tag)
    implicit none
    complex(8), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, 1, MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_double_complex

  subroutine send_double_complex_array (z, dest, tag)
    implicit none
    complex(8), dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, size(z), MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_double_complex_array

  subroutine send_logical (f, dest, tag)
    implicit none
    logical, intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, 1, MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_logical

  subroutine send_logical_array (f, dest, tag)
    implicit none
    logical, dimension (:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, size(f), MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_logical_array

  subroutine receive_integer (i, src, tag)
    implicit none
    integer, intent (out) :: i
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, 1, MPI_INTEGER, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_integer

  subroutine receive_integer_array (i, src, tag)
    implicit none
    integer, dimension (:), intent (out) :: i
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, size(i), MPI_INTEGER, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_integer_array

  subroutine receive_real (a, src, tag)
    implicit none
    real(4), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, 1, MPI_REAL, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_real

  subroutine receive_real_array (a, src, tag)
    implicit none
    real(4), dimension (:), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, size(a), MPI_REAL, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_real_array

  subroutine receive_double (a, src, tag)
    implicit none
    real(8), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, 1, MPI_DOUBLE_PRECISION, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_double

  subroutine receive_double_array (a, src, tag)
    implicit none
    real(8), dimension (:), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, size(a), MPI_DOUBLE_PRECISION, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_double_array


  subroutine receive_complex (z, src, tag)
    implicit none
    complex(4), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, 1, MPI_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_complex

  subroutine receive_complex_array (z, src, tag)
    implicit none
    complex(4), dimension (:), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, size(z), MPI_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_complex_array

  subroutine receive_double_complex (z, src, tag)
    implicit none
    complex(8), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, 1, MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_double_complex

  subroutine receive_double_complex_array (z, src, tag)
    implicit none
    complex(8), dimension (:), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, size(z), MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_double_complex_array


  subroutine receive_logical (f, src, tag)
    implicit none
    logical, intent (out) :: f
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, 1, MPI_LOGICAL, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_logical

  subroutine receive_logical_array (f, src, tag)
    implicit none
    logical, dimension (:), intent (out) :: f
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, size(f), MPI_LOGICAL, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_logical_array


!=====================================================================

  subroutine isend_integer (i, dest, req, tag)
    implicit none
    integer, intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (i, 1, MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine isend_integer

  subroutine isend_integer_array (i, dest, req, tag)
    implicit none
    integer, dimension (:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (i, size(i), MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine isend_integer_array

  subroutine isend_real (a, dest, req, tag)
    implicit none
    real(4), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (a, 1, MPI_REAL, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine isend_real

  subroutine isend_real_array (a, dest, req, tag)
    implicit none
    real(4), dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (a, size(a), MPI_REAL, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine isend_real_array

  subroutine isend_double (a, dest, req, tag)
    implicit none
    real(8), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (a, 1, MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine isend_double

  subroutine isend_double_array (a, dest, req, tag)
    implicit none
    real(8), dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (a, size(a), MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine isend_double_array


  subroutine isend_complex (z, dest, req, tag)
    implicit none
    complex(4), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (z, 1, MPI_COMPLEX, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine isend_complex

  subroutine isend_complex_array (z, dest, req, tag)
    implicit none
    complex(4), dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (z, size(z), MPI_COMPLEX, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine isend_complex_array

  subroutine isend_double_complex (z, dest, req, tag)
    implicit none
    complex(8), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (z, 1, MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine isend_double_complex

  subroutine isend_double_complex_array (z, dest, req, tag)
    implicit none
    complex(8), dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (z, size(z), MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine isend_double_complex_array

  subroutine isend_logical (f, dest, req, tag)
    implicit none
    logical, intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (f, 1, MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine isend_logical

  subroutine isend_logical_array (f, dest, req, tag)
    implicit none
    logical, dimension (:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (f, size(f), MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine isend_logical_array

  subroutine ireceive_integer (i, src, req, tag)
    implicit none
    integer, intent (out) :: i
    integer, intent (in) :: src
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (i, 1, MPI_INTEGER, src, tagp, MPI_COMM_WORLD, &
         req, ierror)
  end subroutine ireceive_integer

  subroutine ireceive_integer_array (i, src, req, tag)
    implicit none
    integer, dimension (:), intent (out) :: i
    integer, intent (in) :: src
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (i, size(i), MPI_INTEGER, src, tagp, MPI_COMM_WORLD, &
         req, ierror)
  end subroutine ireceive_integer_array

  subroutine ireceive_real (a, src, req, tag)
    implicit none
    real(4), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (a, 1, MPI_REAL, src, tagp, MPI_COMM_WORLD, &
         req, ierror)
  end subroutine ireceive_real

  subroutine ireceive_real_array (a, src, req, tag)
    implicit none
    real(4), dimension (:), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (a, size(a), MPI_REAL, src, tagp, MPI_COMM_WORLD, &
        req, ierror)
  end subroutine ireceive_real_array

  subroutine ireceive_double (a, src, req, tag)
    implicit none
    real(8), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (a, 1, MPI_DOUBLE_PRECISION, src, tagp, MPI_COMM_WORLD, &
         req, ierror)
  end subroutine ireceive_double

  subroutine ireceive_double_array (a, src, req, tag)
    implicit none
    real(8), dimension (:), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (a, size(a), MPI_DOUBLE_PRECISION, src, tagp, MPI_COMM_WORLD, &
         req, ierror)
  end subroutine ireceive_double_array


  subroutine ireceive_complex (z, src, req, tag)
    implicit none
    complex(4), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (z, 1, MPI_COMPLEX, src, tagp, MPI_COMM_WORLD, &
         req, ierror)
  end subroutine ireceive_complex

  subroutine ireceive_complex_array (z, src, req, tag)
    implicit none
    complex(4), dimension (:), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (z, size(z), MPI_COMPLEX, src, tagp, MPI_COMM_WORLD, &
         req, ierror)
  end subroutine ireceive_complex_array

  subroutine ireceive_double_complex (z, src, req, tag)
    implicit none
    complex(8), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (z, 1, MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
         req, ierror)
  end subroutine ireceive_double_complex

  subroutine ireceive_double_complex_array (z, src, req, tag)
    implicit none
    complex(8), dimension (:), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (z, size(z), MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        req, ierror)
  end subroutine ireceive_double_complex_array


  subroutine ireceive_logical (f, src, req, tag)
    implicit none
    logical, intent (out) :: f
    integer, intent (in) :: src
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (f, 1, MPI_LOGICAL, src, tagp, MPI_COMM_WORLD, &
         req, ierror)
  end subroutine ireceive_logical

  subroutine ireceive_logical_array (f, src, req, tag)
    implicit none
    logical, dimension (:), intent (out) :: f
    integer, intent (in) :: src
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (f, size(f), MPI_LOGICAL, src, tagp, MPI_COMM_WORLD, &
        req, ierror)
  end subroutine ireceive_logical_array

!======================================================================

  subroutine rsend_integer (i, dest, tag)
    implicit none
    integer, intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_rsend (i, 1, MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine rsend_integer

  subroutine rsend_integer_array (i, dest, tag)
    implicit none
    integer, dimension (:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_rsend (i, size(i), MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine rsend_integer_array

  subroutine rsend_real (a, dest, tag)
    implicit none
    real(4), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_rsend (a, 1, MPI_REAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine rsend_real

  subroutine rsend_real_array (a, dest, tag)
    implicit none
    real(4), dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_rsend (a, size(a), MPI_REAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine rsend_real_array

  subroutine rsend_double (a, dest, tag)
    implicit none
    real(8), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_rsend (a, 1, MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine rsend_double

  subroutine rsend_double_array (a, dest, tag)
    implicit none
    real(8), dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_rsend (a, size(a), MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine rsend_double_array


  subroutine rsend_complex (z, dest, tag)
    implicit none
    complex(4), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_rsend (z, 1, MPI_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine rsend_complex

  subroutine rsend_complex_array (z, dest, tag)
    implicit none
    complex(4), dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_rsend (z, size(z), MPI_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine rsend_complex_array

  subroutine rsend_double_complex (z, dest, tag)
    implicit none
    complex(8), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_rsend (z, 1, MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine rsend_double_complex

  subroutine rsend_double_complex_array (z, dest, tag)
    implicit none
    complex(8), dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_rsend (z, size(z), MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine rsend_double_complex_array

  subroutine rsend_logical (f, dest, tag)
    implicit none
    logical, intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_rsend (f, 1, MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine rsend_logical

  subroutine rsend_logical_array (f, dest, tag)
    implicit none
    logical, dimension (:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_rsend (f, size(f), MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine rsend_logical_array
!======================================================================
  subroutine irsend_integer (i, dest, req, tag)
    implicit none
    integer, intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irsend (i, 1, MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine irsend_integer

  subroutine irsend_integer_array (i, dest, req, tag)
    implicit none
    integer, dimension (:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irsend (i, size(i), MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine irsend_integer_array

  subroutine irsend_real (a, dest, req, tag)
    implicit none
    real(4), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irsend (a, 1, MPI_REAL, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine irsend_real

  subroutine irsend_real_array (a, dest, req, tag)
    implicit none
    real(4), dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irsend (a, size(a), MPI_REAL, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine irsend_real_array

  subroutine irsend_double (a, dest, req, tag)
    implicit none
    real(8), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irsend (a, 1, MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine irsend_double

  subroutine irsend_double_array (a, dest, req, tag)
    implicit none
    real(8), dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irsend (a, size(a), MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine irsend_double_array


  subroutine irsend_complex (z, dest, req, tag)
    implicit none
    complex(4), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irsend (z, 1, MPI_COMPLEX, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine irsend_complex

  subroutine irsend_complex_array (z, dest, req, tag)
    implicit none
    complex(4), dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irsend (z, size(z), MPI_COMPLEX, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine irsend_complex_array

  subroutine irsend_double_complex (z, dest, req, tag)
    implicit none
    complex(8), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irsend (z, 1, MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine irsend_double_complex

  subroutine irsend_double_complex_array (z, dest, req, tag)
    implicit none
    complex(8), dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irsend (z, size(z), MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine irsend_double_complex_array

  subroutine irsend_logical (f, dest, req, tag)
    implicit none
    logical, intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irsend (f, 1, MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine irsend_logical

  subroutine irsend_logical_array (f, dest, req, tag)
    implicit none
    logical, dimension (:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (out) :: req
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irsend (f, size(f), MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, req, ierror)
  end subroutine irsend_logical_array

!======================================================================
  subroutine max_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_integer

  subroutine max_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_integer_array

  subroutine max_reduce_double (a, dest)
    implicit none
    real(8), intent (in out) :: a
    integer, intent (in) :: dest
    real(8) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_double

  subroutine max_reduce_double_array (a, dest)
    implicit none
    real(8), dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    real(8), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_double_array


  subroutine max_reduce_real (a, dest)
    implicit none
    real(4), intent (in out) :: a
    integer, intent (in) :: dest
    real(4) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_REAL, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_real

  subroutine max_reduce_real_array (a, dest)
    implicit none
    real(4), dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    real(4), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_REAL, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_real_array


  subroutine max_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_integer

  subroutine max_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_integer_array

  subroutine max_allreduce_real (a)
    implicit none
    real(4), intent (in out) :: a
    real(4) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_real

  subroutine max_allreduce_real_array (a)
    implicit none
    real(4), dimension (:), intent (in out) :: a
    real(4), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_real_array

  subroutine max_allreduce_double (a)
    implicit none
    real(8), intent (in out) :: a
    real(8) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_double

  subroutine max_allreduce_double_array (a)
    implicit none
    real(8), dimension (:), intent (in out) :: a
    real(8), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_double_array


  subroutine min_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_integer


  subroutine min_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_integer_array

  subroutine min_reduce_real (a, dest)
    implicit none
    real(4), intent (in out) :: a
    integer, intent (in) :: dest
    real(4) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_REAL, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_real

  subroutine min_reduce_real_array (a, dest)
    implicit none
    real(4), dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    real(4), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_REAL, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_real_array

  subroutine min_reduce_double (a, dest)
    implicit none
    real(8), intent (in out) :: a
    integer, intent (in) :: dest
    real(8) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_double

  subroutine min_reduce_double_array (a, dest)
    implicit none
    real(8), dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    real(8), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_double_array


  subroutine min_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_integer

  subroutine min_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_integer_array

  subroutine min_allreduce_real (a)
    implicit none
    real(4), intent (in out) :: a
    real(4) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_real

  subroutine min_allreduce_real_array (a)
    implicit none
    real(4), dimension (:), intent (in out) :: a
    real(4), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_REAL, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_real_array

  subroutine min_allreduce_double (a)
    implicit none
    real(8), intent (in out) :: a
    real(8) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_double

  subroutine min_allreduce_double_array (a)
    implicit none
    real(8), dimension (:), intent (in out) :: a
    real(8), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_double_array

  subroutine or_allreduce (a)
    implicit none
    logical :: a
    logical :: a1
    integer :: ierror
     
    a1=a 
    call MPI_Allreduce(a1,a,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierror)

  end subroutine or_allreduce
end module mp
