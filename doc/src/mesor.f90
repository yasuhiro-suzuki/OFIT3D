SUBROUTINE mesor ( rstart, pstart, zstart, w_in, stalam, dt_in, nstep, & !(in)
  &                loss, loss_time, trap                               & !(out)
  &              )

  USE kind_spec
  USE param1,                  ONLY : pi,             &
    &                                 qom_a
  USE param2,                  ONLY : lvessel,        &
    &                                 lmonitor,       &
    &                                 interval
  USE cylindrical_coord_mod,   ONLY : pi2m,           &
    &                                 mgval1
  USE vessel_mod,              ONLY : check_vessel,   &
    &                                 vessel_loss

  IMPLICIT NONE

!Arguments
  INTEGER, INTENT(IN) :: nstep
  INTEGER, INTENT(OUT) :: loss, &
    &                     trap
  REAL(RP), INTENT(IN) :: rstart, &
    &                     pstart, &
    &                     zstart, &
    &                     w_in,   &
    &                     stalam, &
    &                     dt_in
  REAL(RP), INTENT(OUT) :: loss_time
!Local variables
  INTEGER :: i
  REAL(RP) :: v,        &
    &         vpara,    &
    &         vperp,    &
    &         bdotv,    &
    &         vperp2,   &
    &         vperpr,   &
    &         vperpp,   &
    &         vperpz,   &
    &         br,       &
    &         bp,       &
    &         bz,       &
    &         bb,       &
    &         mu,       &
    &         wtotl,    &
    &         werror,   &
    &         wpara,    &
    &         wperp, rho_perp, rho(3)
  CHARACTER(LEN=100) :: fmt
!For Runge-Kutta
  INTEGER :: iout
  REAL(RP) :: t,   &
    &         y(6)
!External function
  EXTERNAL :: fun_mesor
!For loss detection
  INTEGER :: ivessel,  &
    &        idivertor
  REAL(RP) :: dt_loss,       &
    &         y0(3),         &
    &         loss_points(3)


  CALL mgval1(rstart, pstart, zstart, br, bp, bz, bb )

  IF(bb == 0.0_RP)THEN
    PRINT*,'Initial value(follow) is invalid.'
    RETURN
  END IF

  v     =  SQRT(2 * w_in * qom_a)
  vpara =  v * COS(stalam * pi)
  vperp =  v * SIN(stalam * pi)
 
  fmt  = '(A20, ES12.4)'

  rho_perp =  vperp / bb / qom_a

  PRINT *
  PRINT *,   '------------------------------------'
  PRINT fmt, '  gyro radius [m] : ', rho_perp

  y(1)  =  0.0_RP + vpara * br / bb
  y(2)  =  0.0_RP + vpara * bp / bb / rstart
  y(3)  = -vperp  + vpara * bz / bb
  y(4)  =  rstart
  y(5)  =  pstart
  y(6)  =  zstart
  t     =  0.0_RP
  loss  =  0
  trap  =  0

  rho(1) = -(rstart * y(3) * bz -          y(2) * bp) / bb**2 / qom_a
  rho(2) = -(         y(2) * br -          y(1) * bz) / bb**2 / qom_a
  rho(3) = -(         y(1) * bp - rstart * y(3) * br) / bb**2 / qom_a

  PRINT fmt, '  rho(R) [m]      : ', rho(1)
  PRINT fmt, '  rho(phi) [rad]  : ', rho(2) / rstart
  PRINT fmt, '  rho(Z) [m]      : ', rho(3)
  PRINT fmt, '  r_c(R) [m]      : ', rstart - rho(1)
  PRINT fmt, '  r_c(phi) [rad]  : ', pstart - rho(2) / rstart
  PRINT fmt, '  r_c(phi) [a.u.] : ', (pstart - rho(2) / rstart) / pi2m
  PRINT fmt, '  r_c(Z) [m]      : ', zstart - rho(3)
  PRINT *

  IF(lmonitor)then
    bdotv  = (br * y(1) + bp * y(2) * y(4) + bz * y(3)) / bb
    vperpr =  y(1) - bdotv * br / bb
    vperpp =  y(2) - bdotv * bp / (bb * y(4))
    vperpz =  y(3) - bdotv * bz / bb
    vperp2 =  vperpr**2 + vperpp**2 * y(4)**2 + vperpz**2
    mu     =  0.5_RP * vperp2 / (bb * qom_a)
    wpara  =  0.5_RP * bdotv * ABS(bdotv) / qom_a
    wperp  =  mu * bb
    WRITE(50,'(20(ES22.12, A1))') t, ',',                & ! 1: t=0 [s]
      &                           y(4) * COS(y(5)), ',', & ! 2: X@t=0 [m]
      &                           y(4) * SIN(y(5)), ',', & ! 3: Y@t=0 [m]
      &                           y(6), ',',             & ! 4: Z@t=0 [m]
      &                           y(4), ',',             & ! 5: R@t=0 [m]
      &                           y(5), ',',             & ! 6: phi@t=0 [rad]
      &                           mu, ',',               & ! 7: mu@t=0
      &                           bb, ',',               & ! 8: B@t=0 [T]
      &                           wpara, ',',            & ! 9: W_||@t=0 [eV]
      &                           wperp, ',',            & ! 10: W_perp@t=0 [eV]
      &                           w_in, ',',             & ! 11: initial W@t=0 [eV]
      &                           0.0_RP, ','              ! 12: relative error of W@t=0
  END IF

  loop10 : DO i=1,nstep

!  -------------------------------------------------
    CALL rkg6(t, dt_in, 6, fun_mesor, y(:), iout)
!  -------------------------------------------------

    !CALL check_vessel(y(4), y(5), y(6), ivessel)

    IF(ivessel == 1)THEN
      IF(lmonitor)THEN
        CALL mgval1(y(4), y(5), y(6), br, bp, bz, bb)
        bdotv  = (br * y(1) + bp * y(2) * y(4) + bz * y(3)) / bb
        vperpr =  y(1) - bdotv * br / bb 
        vperpp =  y(2) - bdotv * bp / (bb * y(4))
        vperpz =  y(3) - bdotv * bz / bb 
        vperp2 =  vperpr**2 + vperpp**2 * y(4)**2 + vperpz**2
        mu     =  0.5_RP * vperp2 / (bb * qom_a)
        wpara  =  0.5_RP * bdotv * ABS(bdotv) / qom_a
        wperp  =  mu * bb
        wtotl  =  ABS(wpara) + wperp
        werror = (wtotl - w_in) / w_in
        WRITE(50,'(20(ES22.12, A1))') t + dt_in, ',',        & ! 1: t [s]
          &                           y(4) * COS(y(5)), ',', & ! 2: X [m]
          &                           y(4) * SIN(y(5)), ',', & ! 3: Y [m]
          &                           y(6), ',',             & ! 4: Z [m]
          &                           y(4), ',',             & ! 5: R [m]
          &                           y(5), ',',             & ! 6: phi [rad]
          &                           mu, ',',               & ! 7: mu
          &                           bb, ',',               & ! 8: B [T]
          &                           wpara, ',',            & ! 9: W_|| [eV]
          &                           wperp, ',',            & ! 10: W_perp [eV]
          &                           wtotl, ',',            & ! 11: total W [eV]
          &                           werror, ','              ! 12: relative error of W
      END IF
      loss      =  1
      loss_time =  t + dt_in
      RETURN
    END IF

    IF(lmonitor)THEN
      IF(mod(i,interval) == 0)THEN
        CALL mgval1(y(4), y(5), y(6), br, bp, bz, bb)
        bdotv  = (br * y(1) + bp * y(2) * y(4) + bz * y(3)) / bb
        vperpr =  y(1) - bdotv * br / bb 
        vperpp =  y(2) - bdotv * bp / (bb * y(4))
        vperpz =  y(3) - bdotv * bz / bb 
        vperp2 =  vperpr**2 + vperpp**2 * y(4)**2 + vperpz**2
        mu     =  0.5_RP * vperp2 / (bb * qom_a)
        wpara  =  0.5_RP * bdotv * ABS(bdotv) / qom_a
        wperp  =  mu * bb
        wtotl  =  ABS(wpara) + wperp
        werror = (wtotl - w_in) / w_in
        WRITE(50,'(20(ES22.12, A1))') t + dt_in, ',',        & ! 1: t [s]
          &                           y(4) * COS(y(5)), ',', & ! 2: X [m]
          &                           y(4) * SIN(y(5)), ',', & ! 3: Y [m]
          &                           y(6), ',',             & ! 4: Z [m]
          &                           y(4), ',',             & ! 5: R [m]
          &                           y(5), ',',             & ! 6: phi [rad]
          &                           mu, ',',               & ! 7: mu
          &                           bb, ',',               & ! 8: B [T]
          &                           wpara, ',',            & ! 9: W_|| [eV]
          &                           wperp, ',',            & ! 10: W_perp [eV]
          &                           wtotl, ',',            & ! 11: total W [eV]
          &                           werror, ','              ! 12: relative error of W
      END IF
    END IF

    bdotv  = (br * y(1) + bp * y(2) * y(4) + bz * y(3)) / bb

    IF((i /= 1) .and. (bdotv * vpara < 0.0_RP))THEN
      trap =  1
    END IF

    t = t + dt_in
  END DO loop10

  loss_time = t


  RETURN
END SUBROUTINE mesor
