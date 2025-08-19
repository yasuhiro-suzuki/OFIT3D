!=precal.f90
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
SUBROUTINE precal

  USE kind_spec
  USE param1,                  ONLY : ee,             &
    &                                 mi,             &
    &                                 mn,             &
    &                                 me,             &
    &                                 epsi0,          &
    &                                 pi4,            &
    &                                 qa,             &
    &                                 qb,             &
    &                                 ma,             &
    &                                 mb,             &
    &                                 qom_a,          &
    &                                 qom_b
  USE param2,                  ONLY : lcoll,          &
    &                                 wk_in,          &
    &                                 dt,             &
    &                                 mode
  USE cylindrical_coord_mod,   ONLY : magset
  USE collision_parameter_mod, ONLY : collision_step, &
    &                                 n_e,            &
    &                                 t_bg,           &
    &                                 d_ano,          &
    &                                 rho_col,        &
    &                                 ln_lamda,       &
    &                                 v_th_a,         &
    &                                 v_th_b,         &
    &                                 nu_ab,          &
    &                                 tau,            &
    &                                 nu_d,           &
    &                                 nu_dxtau,       &
    &                                 nu_e,           &
    &                                 nu_extau
  USE info_mod,                ONLY : pre_time
  USE mpi_param_mod,           ONLY : myrank
  USE mpi
 
  IMPLICIT NONE

  REAL(RP) ::    &
    &         v,        &
    &         xa,       &
    &         xb,       &
    &         upsilon,  &
    &         xsi,      &
    &         dupsilon
  REAL(RP), EXTERNAL :: upsilon_cal, &
    &                   xsi_cal
  CHARACTER(LEN=100) :: fmt
!
  REAL(RP) :: t1, &
    &         t2
!for MPI
  INTEGER :: ierr


  t1       =  MPI_WTIME()

  IF(myrank == 0)THEN
    PRINT *
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *, '          SUBROUTINE PRECALC                                     '
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *
  END IF

!
! set test particle
!

  SELECT CASE(mode)
    CASE('proton')
      ma    =  mi
      mb    =  me
      qa    =  ee
      qb    = -ee
      qom_a =  ABS(qa) / ma
      qom_b =  ABS(qb) / mb
    CASE('deuteron')
      ma    =  mi + mn
      mb    =  me
      qa    =  ee
      qb    = -ee
      qom_a =  ABS(qa) / ma
      qom_b =  ABS(qb) / mb
    CASE('triton')
      ma    =  mi + mn + mn
      mb    =  me
      qa    =  ee
      qb    = -ee
      qom_a =  ABS(qa) / ma
      qom_b =  ABS(qb) / mb
    CASE('electron')
      ma    =  me
      mb    =  me
      qa    = -ee
      qb    = -ee
      qom_a =  ABS(qa) / ma
      qom_b =  ABS(qb) / mb
    CASE('alpha')
      ma    =  mi + mi + mn + mn
      mb    =  me
      qa    =  ee + ee
      qb    = -ee
      qom_a =  ABS(qa) / ma
      qom_b =  ABS(qb) / mb
    CASE('c6+', 'C6+')
      ma    =  6 * mi + 6 * mn
      mb    =  me
      qa    =  6 * ee
      qb    = -ee
      qom_a =  ABS(qa) / ma
      qom_b =  ABS(qb) / mb
    CASE DEFAULT
      ma    =  mi
      mb    =  me
      qa    =  ee
      qb    = -ee
      qom_a =  ABS(qa) / ma
      qom_b =  ABS(qb) / mb
  END SELECT

!
! read magnetic field
!

  CALL magset

!
! if collisionless case, return
!

  IF(lcoll)THEN

    ln_lamda =  15.2_RP - 0.5_RP * LOG(n_e / 1.0e+20_RP) + LOG(t_bg / 1.0e+03_RP)

    v_th_a   =  SQRT(2 * t_bg  * qom_a)
    v_th_b   =  SQRT(2 * t_bg  * qom_b)
    v        =  SQRT(2 * wk_in * qom_a)

    IF(myrank == 0)THEN

      fmt = '(A)'

      PRINT *
      PRINT fmt, ' PARAMETERS OF BACKGROUD PLASMA'
      PRINT fmt, '------------------------------------------------'
      PRINT fmt, '    Te[eV]         n_e[m^-3]     ln_lamda'

      fmt = '(3ES15.7)'

      PRINT fmt, t_bg, n_e, ln_lamda

    END IF

    nu_ab    =  n_e * qa**2 * qb**2 * ln_lamda / (pi4 * epsi0**2 * ma**2 * v_th_a**3)
    xa       =  v / v_th_a
    xb       =  v / v_th_b
    upsilon  =  upsilon_cal(xb)
    xsi      =  xsi_cal(xb)
    nu_d     =  nu_ab * (upsilon-xsi) / xa**3
    nu_e     =  2 * nu_ab * xsi / xa**3

    IF(myrank == 0)THEN

      fmt = '(A)'

      PRINT *
      PRINT fmt, ' COLLISION FREQUENCY'
      PRINT fmt, '--------------------------------------------------------------------------------------------'
      PRINT fmt, '    v[m/s]        v_th_a[m/s]    v_th_b[m/s]     nu_ab[1/s]     nu_d[1/s]      nu_e[1/s]'

      fmt = '(6ES15.7)'

      PRINT fmt, v, v_th_a, v_th_b, nu_ab, nu_d, nu_e

    END IF

    tau      = dt
    nu_dxtau = nu_d * tau
    nu_extau = nu_e * tau

    IF(myrank == 0)THEN

      fmt = '(A)'

      PRINT *
      PRINT fmt, ' MONTECARLO PARAMETERS'
      PRINT fmt, '------------------------------------------------'
      PRINT fmt, '    tau[s]       nu_d x tau     nu_e x tau'
  
      fmt = '(3ES15.7)'

      PRINT fmt, tau, nu_dxtau, nu_extau

    END IF

    collision_step(2) =  1.0_RP / nu_d
    !rho_col           =  1.2e-03_RP
    d_ano             =  rho_col**2 / (dt * collision_step(2))

    IF(myrank == 0)THEN

      fmt = '(A)'

      PRINT *
      PRINT fmt, ' ANOUMALOUS TRANSPORT'
      PRINT fmt, '------------------------------------------------'
      PRINT *, '  rho_ano[m]     tau_ano[s]     D_ano[m^2/s]'

      fmt = '(3ES15.7)'

      PRINT fmt, rho_col, dt * collision_step(2), d_ano

    END IF

  END IF

  IF(myrank == 0)THEN
    PRINT *
    PRINT *
  END IF

  t2       =  MPI_WTIME()

  pre_time =  t2 - t1

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)


END SUBROUTINE precal
