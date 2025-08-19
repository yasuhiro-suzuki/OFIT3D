!=mgcpu.f90
!
!==Version
!
! $Revision: $
! $Id: $
!
!==Overview
!
!==Reference
!
!==Error Handlings
!
!==Known Bugs
!
!==Note
!
!==TODO
!
SUBROUTINE mgcpu (itime, mess & ! (in)
  &              )

  USE kind_spec
#ifdef MPI
  USE mpi_param_mod, ONLY : myrank
  USE mpi
#endif

  IMPLICIT NONE

  INTEGER :: i,     &
    &        itime, &
    &        min1,  &
    &        min2
  !INTEGER, SAVE :: it =  0
  INTEGER, SAVE :: it = -1
  REAL(DP) :: ccp,   &
    &         sec1,  &
    &         sec2,  &
    &         etime, &
    &         ta(2)
  REAL(DP), SAVE :: cpu1(100), &
    &               cpu2(100)
  CHARACTER(LEN=20) :: mess,      &
    &                  mess1(100)
  CHARACTER(LEN=100) :: fmt


#ifdef MPI
  IF(myrank == 0)THEN
#endif

  it  = it + 1

#ifdef MPI
  ccp =  MPI_WTIME()
  IF(itime == -1)THEN
    ta(1) =  ccp
    ta(2) =  ccp
    RETURN
  END IF
#else
!  ccp = etime(ta)
  CALL CPU_TIME(ccp)
#endif

  IF(it == 1)THEN
    cpu1(it) =  ccp
    cpu2(it) =  ccp
    cpu1(it) =  ccp
    cpu2(it) =  ccp
  ELSE
    cpu1(it) =  ccp
    cpu2(it) =  ccp - cpu1(it-1)
    !cpu2(it) =  ccp - cpu2(it-1)
    !cpu1(it) =  cpu1(it-1) + cpu2(it)
  END IF

  mess1(it) =  mess

  IF(itime == 0) RETURN

  PRINT *

  fmt = '(I4, 2X, A20, A15, I5, A4, F9.3, A4, A13, I5, A4, F9.3, A4)'

  loop010 : DO  i=1,it
    min1 =  cpu1(i) / 60
    min2 =  cpu2(i) / 60
    sec1 =  cpu1(i) - min1 * 60
    sec2 =  cpu2(i) - min2 * 60
    PRINT fmt, i, mess1(i), '; Exec Time ...', min2, ' min', sec2, ' sec', &
      &                     '; Total ... ',    min1, ' min', sec1, ' sec'
  END DO loop010

#ifdef MPI
  END IF
#endif


END SUBROUTINE mgcpu
