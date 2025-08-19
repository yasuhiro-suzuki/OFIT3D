!=mpi_param_mod.f90
!
!==Version
!
! $Revision: $
! $Id: $
!
!==Overview
!
! SUBROUTINEs for MPI communication
!
!==Reference
!
! None
!
!==Error Handlings
!
! None
!
!==Known Bugs
!
! None
!
!==Note
!
! None
!
!==TODO
!
! Change to MPI2
!

MODULE mpi_param_mod

  USE kind_spec
  USE mpi

  IMPLICIT NONE

  PRIVATE

!  include 'mpif.h'
  INTEGER :: myrank, &
    &        nprocs, &
    &        nsta,   &
    &        nend

  PUBLIC :: myrank,    &
    &       nprocs,    &
    &       nsta,      &
    &       nend,      &
    &       para_range

CONTAINS

  SUBROUTINE para_range (m1, m2, mprocs, mrank, & ! (in)
    &                    msta, mend             & ! (out)
    &                   )

    IMPLICIT NONE
!Arguments
    INTEGER, INTENT(IN) :: m1,     &
      &                    m2,     &
      &                    mprocs, &
      &                    mrank
    INTEGER, INTENT(OUT) :: msta, &
      &                     mend
!Local variables
    INTEGER :: mwork1, &
      &        mwork2


    mwork1 = (m2 - m1 + 1) / mprocs
    mwork2 =  MOD(m2 - m1 + 1, mprocs)

    msta   =  mrank * mwork1 + m1 + MIN(mrank, mwork2)
    mend   =  msta + mwork1 - 1

    IF(mwork2 > mrank) mend =  mend + 1


    RETURN
  END SUBROUTINE para_range

 
END MODULE mpi_param_mod
