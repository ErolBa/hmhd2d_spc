      MODULE util
      USE globals
      CONTAINS
!#######################################################################
      function has (n,table,item)
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! ****** Check whether TABLE(N) contains ITEM.
!
!-----------------------------------------------------------------------
!
      LOGICAL :: has
      INTEGER :: n
      CHARACTER(LEN=*) :: table(n),item

      INTEGER :: i
      
!
!-----------------------------------------------------------------------
!
      has=.true.
!
      do 100 i=1,n
        if (table(i).eq.item) return
  100 continue
!
      has=.false.
!
      return
      end function has
!#######################################################################
      function letter (ch)
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! ****** Check if character CH is a letter.
!
!-----------------------------------------------------------------------
!
      LOGICAL :: letter
      CHARACTER :: ch
      INTEGER :: ich 
!
!-----------------------------------------------------------------------
!
      letter=.false.
!
      ich=ichar(ch)
!
      if ((ich.ge.65.and.ich.le. 90).or. &
          (ich.ge.97.and.ich.le.122)) letter=.true.
!
      return
      end function letter
!#######################################################################
      function digit (ch)
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! ****** Check if character CH is a digit.
!
!-----------------------------------------------------------------------
!
      LOGICAL :: digit
      CHARACTER :: ch
       
      INTEGER(ILWORD) :: ich
!
!-----------------------------------------------------------------------
!
      digit=.false.
!
      ich=ichar(ch)
!
      if (ich.ge.48.and.ich.le. 57) digit=.true.
!
      return
      end  function digit
!#######################################################################
      function alpha (ch)
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! ****** Check if CH is an alphanumeric character.
!
!-----------------------------------------------------------------------
!
      LOGICAL :: alpha
      CHARACTER, INTENT(IN) ::  ch
!
!-----------------------------------------------------------------------
!
      alpha=.false.
!
      if (letter(ch).or.digit(ch)) alpha=.true.
!
      return
      end  function alpha

!#######################################################################
      function lenstr (str)
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! ****** Get the length of string STR, ignoring right-filling blanks.
!
!-----------------------------------------------------------------------
!
      INTEGER :: lenstr
      CHARACTER(LEN=*) :: str

      INTEGER :: i 
!
!-----------------------------------------------------------------------
!
      do i=len(str),1,-1
        if (str(i:i).ne.' ') go to 100
      enddo
      i=0
  100 continue
!
      lenstr=i
!
      return
      end function lenstr

!#######################################################################
      subroutine adjstrl (str)
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! ****** Left-adjust the character string STR; i.e., shift the
! ****** characters so there are no blanks on the left end of
! ****** the string.
!
!-----------------------------------------------------------------------
!
      CHARACTER(LEN=*) :: str

      INTEGER :: lstr, i
!
!-----------------------------------------------------------------------
!
      lstr=len(str)
!
! ****** Find the first nonblank character.
!
      do i=1,lstr
        if (str(i:i).ne.' ') go to 200
      end do
!
! ****** All-blank string; leave untouched.
!
      return
!
! ****** Adjust to the left.
!
  200 continue
!
      if (i.gt.1) str(1:)=str(i:lstr)
!
      return
      end subroutine adjstrl

!######################################################################
      FUNCTION ran(idum)
      USE globals

      IMPLICIT NONE

      INTEGER(ILWORD), INTENT(INOUT) :: idum
      REAL(RTYPE) :: ran

!Minimal random number generator of Park and Miller combined 
!with a Marsaglia shift
!sequence. Returns a uniform random deviate between 0.0 and 1.0 
!(exclusive of the endpoint
!values). This fully portable, scalar generator has the
!traditional (not Fortran 90) calling
!sequence with a random deviate as the returned function 
!value: call with idum a negative
!integer to initialize; thereafter, do not alter idum except to reinitialize.
!The period of this generator is about 3.11E18.

      INTEGER(ILWORD), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
      REAL(RTYPE), SAVE :: am
      INTEGER(ILWORD), SAVE :: ix=-1,iy=-1,k

      if (idum <= 0 .or. iy < 0) then       !Initialize.
        am=nearest(1.0_RTYPE,-1.0_RTYPE)/IM
        iy=ior(ieor(888889999,abs(idum)),1)
        ix=ieor(777755555,abs(idum))
        idum=abs(idum)+1              !Set idum positive.
      end if
      ix=ieor(ix,ishft(ix,13))        !Marsaglia shift sequence with period 2^32-1.
      ix=ieor(ix,ishft(ix,-17))
      ix=ieor(ix,ishft(ix,5))
      k=iy/IQ                         !Park-Miller sequence by Schrages method,
      iy=IA*(iy-k*IQ)-IR*k            !   period 2^31-2.
      if (iy < 0) iy=iy+IM
      ran=am*ior(iand(IM,ieor(ix,iy)),1) !Combine the two generators with masking to
      END FUNCTION ran                   !   ensure nonzero value.


!#######################################################################
      subroutine getlun (ilun,ierr)
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! ****** Return an unused logical unit identifier.
!
!-----------------------------------------------------------------------
!
! ****** Upon successful completion (IERR=0), the first
! ****** unused logical unit number between MINLUN and
! ****** MAXLUN, inclusive, is returned in variable ILUN.
! ****** If all units between these limits are busy,
! ****** IERR=1 is returned.
!
!-----------------------------------------------------------------------
!
      INTEGER :: ilun,ierr
!
!-----------------------------------------------------------------------
!
! ****** Range of valid units.
!
      INTEGER, PARAMETER :: minlun=30,maxlun=49
      LOGICAL busy
      INTEGER :: i
!
!-----------------------------------------------------------------------
!
      ierr=0
!
! ****** Find an unused unit number.
!
      do 100 i=minlun,maxlun
        inquire (unit=i,opened=busy)
        if (.not.busy) go to 200
  100 continue
!
! ****** Fall through here if all units are busy.
!
      ierr=1
      return
!
! ****** An unused unit was found.
!
  200 continue
!
      ilun=i
!
      return
      end subroutine getlun
!#######################################################################
      function second()
      USE globals
      implicit none     
!******** RETURNS THE ELAPSED CPU TIME, IN SECONDS.
!******** THIS ROUTINE IS OPERATING SYSTEM-DEPENDENT.
     
!*        UNIX VERSION
      REAL(RTYPE) :: second
      REAL :: tarr
      call cpu_time(tarr)
      second = tarr
      return
      end function second
!#######################################################################


      END MODULE util
