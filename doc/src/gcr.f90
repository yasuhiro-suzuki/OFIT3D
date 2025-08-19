!=gcr.f90
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
SUBROUTINE gcr(rstart, pstart, zstart, w_in, stalam, dt_in, nstep, & !(in)
  &            loss, loss_time, trap                               & !(out)
  &           )

  USE kind_spec
  USE param1,                  ONLY : pi,             &
    &                                 pi2,            &
    &                                 qom_a,          &
    &                                 ma,             &
    &                                 mb
  USE param2,                  ONLY : ladaptive,      &
    &                                 lcoll,          &
    &                                 lano,           &
    &                                 lmonitor,       &
    &                                 lpunc,          &
    &                                 lvessel,        &
    &                                 ldivertor,      &
    &                                 interval,       &
    &                                 collision_type, &
    &                                 anomalous_type, &
    &                                 slimit
  USE inv_val_mod,             ONLY : mu
  USE cylindrical_coord_mod,   ONLY : sedge,          &
    &                                 pi2m,           &
    &                                 mgval1,         &
    &                                 mgval3, mgval2
  USE collision_parameter_mod, ONLY : nu_dxtau,       &
    &                                 nu_extau,       &
    &                                 collision_step, &
    &                                 rho_col,        &
    &                                 v_th_a,         &
    &                                 v_th_b,         &
    &                                 nu_ab
  USE vessel_mod,              ONLY : check_vessel,   &
    &                                 vessel_loss
  USE divertor_mod,            ONLY : divertor_loss
  USE mt95

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
  INTEGER :: i, &
    &        k
  REAL(RP) :: dt,       &
    &         v,        &
    &         vpara,    &
    &         vperp,    &
    &         vpara1,   &
    &         vperp1,   &
    &         br,       &
    &         bz,       &
    &         bp,       &
    &         bb,       &
    &         s,        &
    &         wtotl,    &
    &         werror,   &
    &         wpara,    &
    &         wperp, bvec(4), dbdr(4), dbdp(4), dbdz(4)
!Collision parameters
  INTEGER :: it_array(3), &
    &        iseed
  REAL(RP) :: random1,   &
    &         random2,   &
    &         pitch,     &
    &         phase,     &
    &         p,         &
    &         gamma_ab,  &
    &         u,         &
    &         erfu,      &
    &         erfpu,     &
    &         dvv2_ave,  &
    &         dvv_ave,   &
    &         dvce2_ave, &
    &         muc,       &
    &         sigma,     &
    &         dvv,       &
    &         dvc,       &
    &         dve
!For Runge-Kutta
  INTEGER :: iout
  REAL(RP) :: t, y(4)
!External function
  EXTERNAL :: fun_gcr
!For loss detection
  INTEGER :: ivessel,  &
    &        idivertor
  REAL(RP) :: dt_loss,       &
    &         y0(3),         &
    &         loss_points(3)


  IF(lcoll)THEN
    CALL itime(it_array(:))
    iseed =  it_array(2) * it_array(3) + it_array(1)
    CALL genrand_init(iseed)
  END IF

  CALL mgval1(rstart, pstart, zstart, br, bp, bz, bb)
  CALL mgval3(rstart, pstart, zstart, s)

  IF(bb == 0.0_RP)THEN
    PRINT *, ' Initial value is invalid. '
    RETURN
  END IF

  dt    =  dt_in

  v     =  SQRT(2 * w_in * qom_a)
  vpara =  v * COS(stalam * pi)
  vperp =  v * SIN(stalam * pi)
  mu    =  vperp**2 / (bb + bb)

  y(1)  =  rstart
  y(2)  =  pstart
  y(3)  =  zstart
  y(4)  =  vpara
  t     =  0.0_RP
  loss  =  0
  trap  =  0

  IF(lmonitor)THEN
    wpara =  0.5_RP * vpara * ABS(vpara) / qom_a
    wperp =  mu * bb / qom_a
#ifdef BIN
    WRITE(50) t, y(1) * COS(y(2)), y(1) * SIN(y(2)), y(3), y(1), y(2), bb, s, wpara, wperp, w_in, 0.0_RP
#else
    WRITE(50,'(20(ES22.12, A1))') t, ',', y(1) * COS(y(2)), ',', y(1) * SIN(y(2)), ',', y(3), ',', y(1), ',', y(2), ',', bb, ',', s, ',', wpara, ',', wperp, ',', w_in, ',', 0.0_RP, ','
#endif
  END IF

  !IF(interval == 0) interval =  10
  !IF(nstep > 1000)  interval =  nstep / 1000
  
  DO i=1,nstep
    y0(1) = y(1)
    y0(2) = y(2)
    y0(3) = y(3)

!  --------------------------------------------
    CALL rkg6(t, dt, 4, fun_gcr, y(:), iout)
!  --------------------------------------------

    CALL mgval3(y(1), y(2), y(3), s)

    IF(s >= sedge)THEN

      IF(ldivertor)THEN
        CALL divertor_loss(y0(1), y0(2), y0(3), y(1), y(2), y(3), idivertor, dt_loss, loss_points)

        IF(idivertor == 1)THEN
#ifdef BIN
          WRITE(52) loss_points(1), loss_points(2), loss_points(3), SIGN(1.0_RP, y(4))
#else
          WRITE(52,'(4ES22.12)') loss_points(1), loss_points(2), loss_points(3), SIGN(1.0_RP, y(4))
#endif
          IF(lmonitor)THEN
            CALL mgval1(y(1), y(2), y(3), br, bp, bz, bb)
            wpara  =  0.5_RP * y(4) * ABS(y(4)) / qom_a
            wperp  =  mu * bb / qom_a
            wtotl  =  ABS(wpara) + wperp
            werror = (wtotl - w_in) / w_in
#ifdef BIN
            WRITE(50) t + dt_loss, y(1) * COS(y(2)), y(1) * SIN(y(2)), y(3), y(1), y(2), bb, s, wpara, wperp, wtotl, werror
#else
            WRITE(50,'(20(ES22.12, A1))') t + dt_loss, ',', y(1) * COS(y(2)), ',', y(1) * SIN(y(2)), ',', y(3), ',', y(1), ',', y(2), ',', bb, ',', s, ',', wpara, ',', wperp, ',', wtotl, ',', werror, ','
#endif
          END IF
          loss      =  1
          loss_time =  t + dt_loss
          RETURN
        END IF
      END IF

      CALL check_vessel(y(1), y(2), y(3), ivessel)

      IF(ivessel == 1)THEN
        CALL vessel_loss(y0(1), y0(2), y0(3), y(1), y(2), y(3), ivessel, dt_loss, loss_points)
#ifdef BIN
        WRITE(51) loss_points(1), loss_points(2), loss_points(3), SIGN(1.0_RP, y(4))
#else
        WRITE(51,'(4ES22.12)') loss_points(1), loss_points(2), loss_points(3), SIGN(1.0_RP, y(4))
#endif
        IF(lmonitor)THEN
          CALL mgval1(y(1), y(2), y(3), br, bp, bz, bb)
          wpara  =  0.5_RP * y(4) * ABS(y(4)) / qom_a
          wperp  =  mu * bb / qom_a
          wtotl  =  ABS(wpara) + wperp
          werror = (wtotl - w_in) / w_in
#ifdef BIN
          WRITE(50) t + dt, y(1) * COS(y(2)), y(1) * SIN(y(2)), y(3), y(1), y(2), bb, s, wpara, wperp, wtotl, werror
#else
          WRITE(50,'(20(ES22.12, A1))') t + dt, ',', y(1) * COS(y(2)), ',', y(1) * SIN(y(2)), ',', y(3), ',', y(1), ',', y(2), ',', bb, ',', s, ',', wpara, ',', wperp, ',', wtotl, ',', werror, ','
#endif
        END IF
        loss      =  1
        loss_time =  t + dt
        RETURN
      END IF

    ELSE

      IF(lcoll)THEN

        SELECT CASE(collision_type)

          CASE('pitch')

            CALL mgval1(y(1), y(2), y(3), br, bp, bz, bb)

            wtotl = (0.5_RP * y(4)**2 + mu * bb) / qom_a
            v     =  SQRT(2 * wtotl * qom_a)

            !pitch =  y(4) * bb / v
            pitch =  y(4)/ v

            CALL genrand_real1(random1)

            random1 =  random1 - 0.5_RP
            p       = (1.0_RP - pitch**2) * nu_dxtau
            pitch   =  pitch * ( 1.0_RP - nu_dxtau ) + SIGN(1.0_RP, random1) * SQRT(p)
            pitch   =  MAX(pitch, -1.0_RP) ! fail safe
            pitch   =  MIN(pitch,  1.0_RP) ! fail safe

            !CALL genrand_real1(random1)

            !wtotl   =  wtotl * (1 - random1 * dt)
            !v       =  SQRT(2.0_RP * wtotl * qom_a)

            y(4)   =  v * pitch
            vperp  =  v * SQRT(1.0_RP - pitch**2)
            mu     =  vperp**2 / (bb + bb)

          CASE('slow')

            vpara     =  y(4)
            vperp     =  SQRT(2 * mu * bb)

            CALL mgval1(y(1), y(2), y(3), br, bp, bz, bb)

            wtotl     = (0.5_RP * y(4)**2 + mu * bb) / qom_a
            v         =  SQRT(2 * wtotl * qom_a)

            gamma_ab  =  nu_ab * v_th_a**3
            u         =  v / (SQRT(2.0_RP) * v_th_b)
            erfu      =  ERF(u)
            erfpu     = (2 / SQRT(pi)) * EXP(-u**2)

            dvv2_ave  =  gamma_ab / (2 * v) * (erfu / u**2 - erfpu / u) * dt * 2
            dvv_ave   = -(1 + ma / mb) * gamma_ab / v**2 * (erfu - u * erfpu) * dt
            dvce2_ave =  gamma_ab / (4 * v) * ((2 - 1 / u**2) * erfu + erfpu /u ) * dt * 2

            muc       =  dvv_ave
            sigma     =  SQRT(dvv2_ave - dvv_ave**2)

            CALL genrand_real3(random1)
            CALL genrand_real1(random2)

            dvv       =  muc + sigma * SQRT(-2.0_RP * LOG(random1)) * SIN(pi2 * random2)

            sigma     =  SQRT(dvce2_ave)

            CALL genrand_real3(random1)
            CALL genrand_real1(random2)

            dvc       =  sigma * SQRT(-2.0_RP * LOG(random1)) * SIN(pi2 * random2)

            sigma     =  SQRT(dvce2_ave)

            CALL genrand_real3(random1)
            CALL genrand_real1(random2)

            dve       =  sigma * SQRT(-2.0_RP * LOG(random1)) * SIN(pi2 * random2)

            vpara1    =  vpara + dvv * vpara / v - dvc * vperp / v
            vperp1    =  SQRT((vperp + dvv * vperp / v + dvc * vpara / v)**2 + dve**2)
 
            vpara     =  vpara1
            vperp     =  vperp1
            mu        =  vperp**2 / (bb + bb)

            y(4)      =  vpara

        END SELECT
      END IF

    END IF

    IF(lano)THEN

      !SELECT CASE('anomalous_type')

      !  CASE('classical')

          IF(MOD(i, collision_step(2)) == 0)THEN

            CALL mgval1(y(1), y(2), y(3), br, bp, bz, bb)

            wtotl = (0.5_RP * y(4)**2 + mu * bb) / qom_a
            v     =  SQRT(2 * wtotl * qom_a)
            pitch =  y(4)/ v
            vperp =  v * SQRT(1.0_RP - pitch**2)

            CALL genrand_real1(random2)

            phase =  pi2 * random2
            y(1)  =  y(1) + rho_col * COS(phase)
            y(2)  =  y(2) + rho_col * SIN(phase)

            CALL mgval1(y(1), y(2), y(3), br, bp, bz, bb)

            mu =  vperp**2 / (bb + bb)

          END IF

      !END SELECT

    END IF

    IF(lmonitor)THEN
      IF(MOD(i, interval) == 0)THEN
        CALL mgval1(y(1), y(2), y(3), br, bp, bz, bb)
        CALL mgval2(y(1), y(2), y(3), bvec, dbdr, dbdp, dbdz)
        wpara  =  0.5_RP * y(4) * ABS(y(4)) / qom_a
        wperp  =  mu * bb / qom_a
        wtotl  =  ABS(wpara) + wperp
        werror = (wtotl - w_in) / w_in
#ifdef BIN
        WRITE(50) t, y(1) * COS(y(2)), y(1) * SIN(y(2)), y(3), y(1), y(2), bb, s, wpara, wperp, wtotl, werror
#else
        WRITE(50,'(20(ES22.12, A1))') t, ',', y(1) * COS(y(2)), ',', y(1) * SIN(y(2)), ',', y(3), ',', y(1), ',', y(2), ',', bb, ',', s, ',', wpara, ',', wperp, ',', wtotl, ',', werror, ','
#endif
      END IF
    END IF

    IF((i /= 1) .and. (y(4) * vpara < 0.0_RP))THEN
      trap =  1
    END IF

    t     =  t + dt

    IF(ladaptive)THEN
      wtotl = (0.5_RP * y(4)**2 + mu * bb) / qom_a
      v     =  SQRT(2 * wtotl * qom_a)
      !dt    =  0.002_RP / v
      dt    =  0.0095720401164064_RP / v
    END IF

  END DO

  loss_time = t


  RETURN
END SUBROUTINE gcr
