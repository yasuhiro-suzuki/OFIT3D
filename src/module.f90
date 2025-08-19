!=module.f90
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
MODULE kind_spec
  INTEGER, PARAMETER :: RP =  SELECTED_REAL_KIND(15)
END MODULE kind_spec

MODULE param1
  USE kind_spec
  IMPLICIT NONE
  REAL(RP), PARAMETER ::  ee    =  1.602176565E-19_RP,                           & ! elementary electric charge
    &                     mi    =  1.672621777E-27_RP,                           & ! mass of proton
    &                     mn    =  1.674927351E-27_RP,                           & ! mass of nuetron
    &                     me    =  9.10938291E-31_RP,                            & ! mass of electron
    &                     epsi0 =  8.85418782E-12_RP,                            & ! permittivity
    &                     pi    =  3.141592653589793238462643383279502884197_RP, & ! pi
    &                     pi2   =  pi  + pi,                                     & ! 2 pi
    &                     pi4   =  pi2 + pi2                                       ! 4 pi
  REAL(RP) :: qa,    & ! charge of test particle
    &         qb,    & ! charge of background
    &         ma,    & ! mass of test particle
    &         mb,    & ! mass of background
    &         qom_a, & ! ratio of charge and mass for test particle
    &         qom_b    ! ratio of charge and mass for background
END MODULE param1

MODULE param2
  USE kind_spec
  IMPLICIT NONE
  LOGICAL :: ladaptive  =  .false.,    &
    &        lcoll      =  .false.,    &
    &        lano       =  .false.,    &
    &        lmonitor   =  .false.,    &
    &        lpunc      =  .false.,    &
    &        lvessel    =  .true.,     &
    &        ldivertor  =  .false.,    &
    &        lfile      =  .false.
  INTEGER :: nparticles =  100,        &
    &        interval   =  1
  REAL(RP) :: wk_in     =  1.0E+03_RP, & ! enegy [eV]
    &         wp_in     =  0.0E+00_RP, & ! potential [V]
    &         dt        =  1.0E-08_RP, & ! time step for DOE [s]
    &         quit_t    =  1.0E-03_RP, &
    &         stalam    =  0.0_RP,     & ! normalized pitch angle
    &         rstart    =  0.0_RP,     & ! strating point along R [m]
    &         pstart    =  0.0_RP,     & ! strating point along phi [deg]
    &         pstart1   =  0.0_RP,     & ! strating point along phi
    &         zstart    =  0.0_RP,     & ! strating point along Z [m]
    &         slimit    =  0.98_RP
  CHARACTER(LEN=20) :: mode           =  'proton',    &
    &                  collision_type =  'slow',      &
    &                  anomalous_type =  'classical', &
    &                  run_mode       =  'single',    &
    &                  multi_mode     =  ''
END MODULE param2

MODULE inv_val_mod
  USE kind_spec
  IMPLICIT NONE
  REAL(RP) :: mu
END MODULE inv_val_mod

MODULE collision_parameter_mod
  USE kind_spec
  IMPLICIT NONE
  INTEGER :: collision_step(2) =  1
  REAL(RP) :: n_e          =  1.0E+19_RP, &
    &         t_bg         =  1.0E+03_RP, &
    &         rho_col      =  1.2E-03_RP, &
    &         d_ano,                      &
    &         atomic_mass,                &
    &         ln_lamda,                   &
    &         v_th_a,                     &
    &         v_th_b,                     &
    &         tau,                        &
    &         nu_ab,                      &
    &         nu_d,                       &
    &         nu_dxtau,                   &
    &         nu_e,                       &
    &         nu_extau
END MODULE collision_parameter_mod

MODULE wall_mod
  USE kind_spec
  IMPLICIT NONE
  REAL(RP) :: rwall, &
    &         zwall, &
    &         awall
END MODULE wall_mod

MODULE file_name_mod
  USE kind_spec
  IMPLICIT NONE
  CHARACTER(LEN=200) :: mag_file,   &
    &                   flx_file,   &
    &                   gcr_file,   &
    &                   punc1_file, &
    &                   punc2_file
  NAMELIST /fopen/ mag_file,   &
    &              flx_file,   &
    &              gcr_file,   &
    &              punc1_file, &
    &              punc2_file
END MODULE file_name_mod

MODULE info_mod
  USE kind_spec
  IMPLICIT NONE
  REAL(RP) :: read_time  =  0.0_RP, &
    &         pre_time   =  0.0_RP, &
    &         drive_time =  0.0_RP, &
    &         tot_time   =  0.0_RP
  CHARACTER(LEN=10) :: date0   =  '', &
    &                  time0   =  '', &
    &                  zone0   =  ''
  CHARACTER(LEN=5) :: ver_info =  ''
END MODULE info_mod
