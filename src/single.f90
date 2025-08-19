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
SUBROUTINE single_particles

  USE kind_spec
  USE param1,        ONLY : pi
  USE param2,        ONLY : ladaptive, &
    &                       lcoll,     &
    &                       lano,      &
    &                       lmonitor,  &
    &                       lpunc,     &
    &                       lvessel,   &
    &                       ldivertor, &
    &                       lfile,     &
    &                       mode,      &
    &                       rstart,    &
    &                       zstart,    &
    &                       pstart,    &
    &                       pstart1,   &
    &                       wk_in,     &
    &                       stalam,    &
    &                       dt,        &
    &                       quit_t
  USE mpi_param_mod, ONLY : myrank

  IMPLICIT NONE

  INTEGER :: nend, &
    &        loss, &
    &        trap, ii, i
  REAL(RP) :: loss_time, &
    &         x,         &
    &         y, temp
  CHARACTER(LEN=100) :: fmt


  IF(myrank == 0)THEN

    PRINT *
    PRINT *, '  RUN CONTROL PARAMETERS '
    PRINT *
    PRINT *, '  ladaptive   lcoll     lano'
    PRINT *, '-------------------------------'

    fmt = '(3X, L5, 5X, L5, 5X, L5)'

    PRINT fmt, ladaptive, lcoll, lano

    PRINT *, '   lvessel  ldivertor'
    PRINT *, '------------------------'

    PRINT fmt, lvessel, ldivertor

    PRINT *, '   lmonitor   lpunc     lfile'
    PRINT *, '--------------------------------'

    PRINT fmt, lmonitor, lpunc, lfile

    PRINT *
    PRINT *
 
    PRINT *, '  RUN_MODE: SINGLE_PARTICLE'


    PRINT *
    PRINT *, '   mode'
    PRINT *, '-----------'

    fmt  = '(A10)'

    PRINT fmt, TRIM(mode)

  END IF

  DO i=1,231
  READ(60,*) ii, rstart, zstart, pstart, wk_in, stalam, temp

  pstart =  pstart1 * pi / 180.0_RP
  nend   =  quit_t / dt
  nend   =  quit_t / ABS(dt)

  IF(myrank == 0)THEN


    fmt  = '(A18, I12, A1, I12)'

    PRINT *
    PRINT *, '    R[m]      phi[deg]      Z[m]       lamda      E_K[keV]'
    PRINT *, '-------------------------------------------------------------'

    fmt  = '(5ES12.4)'

    PRINT fmt, rstart, pstart * 180.0_RP / pi, zstart, COS(stalam * pi), 1.0e-03_RP * wk_in

    PRINT *
    PRINT *, '    dt[s]     t_end[s]       nstep'
    PRINT *, '-------------------------------------'

    fmt  = '(2ES12.4, I12)'

    PRINT fmt, dt, quit_t, nend

  END IF

#ifdef GC
!------------------------------------------------------------------------------------
  CALL gcr(rstart, pstart, zstart, wk_in, stalam, dt, nend, loss, loss_time, trap)
!------------------------------------------------------------------------------------
#elif FULL
!------------------------------------------------------------------------------------
  CALL mesor(rstart, pstart, zstart, wk_in, stalam, dt, nend, loss, loss_time, trap)
!------------------------------------------------------------------------------------
#endif 

  IF(myrank == 0)THEN

    PRINT *
    PRINT *, '  t_loss[s]'
    PRINT *, '-----------------------------------'

    x =  COS(stalam * pi)
    y =  SIN(stalam * pi)

    fmt  = '(ES12.4, A22)'

#ifdef DIGOUT1
    IF((loss == 0) .AND. (trap == 0))THEN
      PRINT fmt, loss_time, "    'untrapped unloss'"
      WRITE(88,'(13(ES12.4, A), I3, A)') rstart, ',', pstart, ',', zstart, ',', stalam, ',', x, ',', y, ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', loss_time, ',', 1, ','
      WRITE(89,'(8(ES12.4, A))') rstart, ',', stalam, ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ','
    ELSE IF((loss == 0) .AND. (trap == 1))THEN
      PRINT fmt, loss_time, "    'trapped unloss'"
      WRITE(88,'(13(ES12.4, A), I3, A)') rstart, ',', pstart, ',', zstart, ',', stalam, ',', 0.0, ',',  0.0, ',', x, ',', y, ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', loss_time, ',', 2, ','
      WRITE(89,'(8(ES12.4, A))') 0.0, ',',  0.0, ',', rstart, ',', stalam, ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ','
    ELSE IF((loss == 1) .and. (trap == 0))THEN
      PRINT fmt, loss_time, "    'untrapped loss'"
      WRITE(88,'(13(ES12.4, A), I3, A)') rstart, ',', pstart, ',', zstart, ',', stalam, ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', x, ',', y, ',', 0.0, ',',  0.0, ',', loss_time, ',', 3, ','
      WRITE(89,'(8(ES12.4, A))') 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', rstart, ',', stalam, ',', 0.0, ',',  0.0, ','
    ELSE IF((loss == 1) .and. (trap == 1))THEN
      PRINT fmt, loss_time, "    'trapped loss'"
      WRITE(88,'(13(ES12.4, A), I3, A)') rstart, ',', pstart, ',', zstart, ',', stalam, ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', x, ',', y, ',', loss_time, ',', 4, ','
      WRITE(89,'(8(ES12.4, A))') 0.0, ',', 0.0, ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', rstart, ',', stalam, ','
    END IF
#endif

    PRINT *
    PRINT *

  END IF
  END DO


END SUBROUTINE single_particles
