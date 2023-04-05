










!******** Identification of the code **********************************
MODULE ident
!
! ****** Code name.
!
  CHARACTER(LEN=*), PARAMETER :: idcode='HMHD2D',   &
                                 vers='3.21',      &
                                 update='05/26/2022'  

!
! ****** Machine names
!
  ! compile machine is obtained from $HOSTNAME in makefile
  CHARACTER(LEN=*), PARAMETER :: compile_machine='spcpc603'
  ! run_machine is obtained in run time
  CHARACTER(LEN=30) :: run_machine                               

!
! ****** Run time and date.
!
  CHARACTER(LEN=8) :: rdate
  CHARACTER(LEN=10) :: rtime

!
! ****** Run ID.
!
  CHARACTER(LEN=16) :: runid


END MODULE ident


