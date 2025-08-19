!=READin.f90
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
SUBROUTINE readin

  USE kind_spec
  USE param2,                  ONLY : ladaptive,      &
    &                                 lcoll,          &
    &                                 lmonitor,       &
    &                                 lpunc,          &
    &                                 lvessel,        &
    &                                 ldivertor,      &
    &                                 lfile,          &
    &                                 nparticles,     &
    &                                 interval,       &
    &                                 wk_in,          &
    &                                 wp_in,          &
    &                                 dt,             &
    &                                 quit_t,         &
    &                                 stalam,         &
    &                                 rstart,         &
    &                                 zstart,         &
    &                                 pstart,         &
    &                                 pstart1,        &
    &                                 mode,           &
    &                                 collision_type, &
    &                                 anomalous_type, &
    &                                 run_mode,       &
    &                                 multi_mode
  USE collision_parameter_mod, ONLY : collision_step, &
    &                                 n_e,            &
    &                                 t_bg,           &
    &                                 rho_col,        &
    &                                 atomic_mass,    &
    &                                 nu_d
  USE vessel_mod,              ONLY : vessel_model,   &
    &                                 lvessel_vtk
  USE divertor_mod,            ONLY : ldivertor_vtk
  USE cylindrical_coord_mod,   ONLY : mag_form,       &
    &                                 flx_form,       &
    &                                 version,        &
    &                                 lflux,          &
    &                                 cj,             &
    &                                 cturn,          &
    &                                 cfact,          &
    &                                 badjust,        &
    &                                 bnorm,          &
    &                                 igrid,          &
    &                                 mbound,         &
    &                                 kstep,          &
    &                                 nlinp_coil_dat
  USE info_mod, ONLY : read_time
  USE mpi

  IMPLICIT NONE

  REAL(RP) :: t1, &
    &         t2
!for MPI
  INTEGER :: ierr

  NAMELIST /nlinp1/ mode,           &
    &               run_mode,       &
    &               multi_mode,     &
    &               ladaptive,      &
    &               lvessel,        &
    &               ldivertor,      &
    &               lcoll,          &
    &               lmonitor,       &
    &               lpunc,          &
    &               lfile,          &
    &               nparticles,     &
    &               interval,       &
    &               wk_in,          &
    &               wp_in,          &
    &               dt,             &
    &               quit_t,         &
    &               stalam,         &
    &               rstart,         &
    &               zstart,         &
    &               pstart,         &
    &               pstart1,        &
    &               vessel_model,   &
    &               lvessel_vtk,    &
    &               ldivertor_vtk
  NAMELIST /nlinp2/ collision_type, &
    &               anomalous_type, &
    &               collision_step, &
    &               n_e,            &
    &               t_bg,           &
    &               rho_col,        &
    &               atomic_mass,    &
    &               nu_d



  t1        =  MPI_WTIME()

  READ(10, nlinp1)
  !WRITE(6, nlinp1)

  IF(lcoll)THEN
    READ(10, nlinp2)
    !WRITE(6, nlinp2)
  END IF

  READ(10, nlinp_coil_dat)

  IF(lfile) CALL file_open

!for fail safe
  lvessel   =  .true.
!  lflux     =  .true.

!default setting is the HINT for the vacuum
!  mag_form  =  'vac'
!  flx_form  =  'flx'

  t2        =  MPI_WTIME()

  read_time =  t2 - t1

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)


END SUBROUTINE readin
