!=driver.f90
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
SUBROUTINE driver

  USE kind_spec
  USE param2,        ONLY : ldivertor,        &
    &                       run_mode
  USE vessel_mod,    ONLY : read_vessel,      &
    &                       free_mem_vessel
  USE divertor_mod,  ONLY : read_divertor,    &
    &                       free_mem_divertor
  USE info_mod,      ONLY : drive_time
  USE mpi_param_mod, ONLY : myrank
  USE mpi

  IMPLICIT NONE
!
  REAL(RP) :: t1, &
    &         t2
!for MPI
  INTEGER :: ierr


  t1         =  MPI_WTIME()

  IF(myrank == 0)THEN
    PRINT *
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *, '          SUBROUTINE DRIVER                                      '
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *
  END IF

  CALL read_vessel

  IF(ldivertor) CALL read_divertor

  SELECT CASE(run_mode)

    CASE('single', 'SINGLE')

      CALL single_particles

    CASE('multi', 'MULTI')

      CALL multi_particles

  END SELECT 

  CALL free_mem_vessel

  IF(ldivertor) CALL free_mem_divertor

  t2         =  MPI_WTIME()

  drive_time =  t2 - t1

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)


END SUBROUTINE driver
