!=main.f90
!
!==Version
!
! $Revision$
! $Id$
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
PROGRAM GCR

  USE param2,                ONLY : lfile
  USE cylindrical_coord_mod, ONLY : free_mem_field
  USE info_mod
  USE mpi_param_mod,         ONLY : nprocs,        &
    &                               myrank
  USE mpi

  IMPLICIT NONE

  CHARACTER(LEN=100) :: fmt
!for MPI
  INTEGER :: ierr


  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

  ver_info =  '1.0.0'

  IF(myrank == 0)THEN
    PRINT *
    PRINT *, ' GCR (MPI) --- GUIDING-CENTER ORBIT FOLLOWING CODE, VERSION ', ver_info, ' ---'
    PRINT *
  END IF

  CALL readin
  CALL precal
  CALL driver

  CALL free_mem_field

  IF(lfile) CALL file_open

  tot_time =  read_time + pre_time + drive_time

  IF(myrank == 0)THEN

    fmt = '(A30, F15.4, A10)'

    PRINT *
    PRINT *
    PRINT fmt, 'TIME IN SUBROUTINE READIN  ', read_time,  'SECONDS'
    PRINT fmt, 'TIME IN SUBROUTINE PRECALC ', pre_time,   'SECONDS'
    PRINT fmt, 'TIME IN SUBROUTINE DRIVER  ', drive_time, 'SECONDS'
    PRINT fmt, 'TOTAL COMPUTATIONAL TIME   ', tot_time,   'SECONDS'

    fmt = '(2X, A6, A2, A1, A2, A1, A4, A6, A2, A1, A2, A1, A2)'

    CALL DATE_AND_TIME(date0, time0, zone0)

    PRINT *
    PRINT fmt, ' DATE= ', date0(5:6), '/', date0(7:8), '/', date0(1:4), &
               ' TIME= ', time0(1:2), ':', time0(3:4), ':', time0(5:6)
    PRINT *
    PRINT *

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL MPI_FINALIZE(ierr)


END PROGRAM GCR
