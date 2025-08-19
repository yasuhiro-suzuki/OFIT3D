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
SUBROUTINE multi_particles

  USE kind_spec
  USE param1,                ONLY : pi,         &
    &                               qom_a
  USE param2,                ONLY : ladaptive,  &
    &                               lcoll,      &
    &                               lano,       &
    &                               lmonitor,   &
    &                               lpunc,      &
    &                               lvessel,    &
    &                               ldivertor,  &
    &                               lfile,      &
    &                               nparticles, &
    &                               interval,   &
    &                               dt,         &
    &                               quit_t,     &
    &                               mode,       &
    &                               multi_mode
  USE cylindrical_coord_mod, ONLY : pi2m
  USE mpi_param_mod,         ONLY : nprocs,     &
    &                               myrank,     &
    &                               nsta,       &
    &                               nend,       &
    &                               para_range
  USE mpi

  IMPLICIT NONE

  INTEGER :: nstep,    &
    &        loss,     &
    &        trap,     &
    &        loss_n_p, &
    &        loss_n_t, &
    &        loss_t_p, &
    &        loss_t_t, &
    &        loss_tol, &
    &        ii,       &
    &        n
  REAL(RP) :: loss_time, &
    &         v,         &
    &         x,         &
    &         y,         &
    &         temp
  REAL(RP), ALLOCATABLE :: rstart(:), &
    &                      zstart(:), &
    &                      pstart(:), &
    &                      wk_in(:),  &
    &                      stalam(:)
  CHARACTER(LEN=100) :: fmt
  CHARACTER(LEN=200) :: label
!for MPI
  INTEGER :: ierr


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

    PRINT *, '  RUN_MODE: MULTI_PARTICLES'


    PRINT *
    PRINT *, '   mode'
    PRINT *, '-----------'

    fmt  = '(A10)'

    PRINT fmt, TRIM(mode)

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  SELECT CASE(multi_mode)

    CASE DEFAULT

      nparticles =  231

      ALLOCATE(rstart(nparticles), zstart(nparticles), pstart(nparticles), wk_in(nparticles), stalam(nparticles))

      DO n=1,nparticles
        READ(60,*) ii, rstart(n), zstart(n), pstart(n), wk_in(n), stalam(n), temp
      END DO

      CALL para_range(1, nparticles, nprocs, myrank, nsta, nend)

      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      IF(myrank == 0)THEN
        PRINT *
        PRINT *, '   PARTICLE DISTRIBUTION ON EACH NODES'
        PRINT *, '----------------------------------------------------------------------------'
      END IF

      fmt = '(A, I10, A, I10, A, I10, I10, A, I10)'
  
      PRINT fmt, '   #rank =', myrank, ' nsta =', nsta, ' nend =', nend, nend - nsta + 1, '/', nparticles

      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      loss_n_p =  0
      loss_n_t =  0
      loss_t_p =  0
      loss_t_t =  0

      DO n=nsta,nend

        pstart(n) =  pi2m * pstart(n)

        v         =  SQRT(2 * wk_in(n) * qom_a)
        !dt        =  0.002_RP / v
        dt        =  0.0095720401164064_RP / v

        nstep     =  quit_t / dt
        interval  =  nstep / 100000


        WRITE(label, fmt='(A,I0,A)') 'multi.', n, '.dat'
        OPEN(11, FILE=label, FORM='formatted', STATUS='unknown')

        IF(lmonitor)THEN
          WRITE(label, fmt='(A,I0,A)') 'multi.', n, '.gcr'
#ifdef BIN
          OPEN(50, FILE=label, FORM='unformatted', STATUS='unknown')
#else
          OPEN(50, FILE=label, FORM='formatted', STATUS='unknown')
#endif
        END IF

        WRITE(label, fmt='(A,I0,A)') 'multi.', n, '.loss'
        OPEN(88, FILE=label, FORM='formatted', STATUS='unknown')

        fmt  = '(A18, I12, A1, I12, A7, I12)'

        WRITE(11,fmt) ' Particle Number =', n, '/', nparticles, ' Rank# ', myrank

        WRITE(11,*)
        WRITE(11,*) '    R[m]      phi[deg]      Z[m]       lamda      E_K[keV]'
        WRITE(11,*) '-------------------------------------------------------------'

        fmt  = '(5ES12.4)'

        WRITE(11,fmt) rstart(n), pstart(n) * 180.0_RP / pi, zstart(n), COS(stalam(n) * pi), 1.0E-03_RP * wk_in(n)

        WRITE(11,*)
        WRITE(11,*) '    dt[s]     t_end[s]       nstep       interval'
        WRITE(11,*) '---------------------------------------------------'

        fmt  = '(2ES12.4, 2I12)'

        WRITE(11,fmt) dt, quit_t, nstep, interval

#ifdef GC
!----------------------------------------------------------------------------------------------------------
        CALL gcr(rstart(n), pstart(n), zstart(n), wk_in(n), stalam(n), dt, nstep, loss, loss_time, trap)
!----------------------------------------------------------------------------------------------------------
#elif FULL
!----------------------------------------------------------------------------------------------------------
        CALL mesor(rstart(n), pstart(n), zstart(n), wk_in(n), stalam(n), dt, nstep, loss, loss_time, trap)
!----------------------------------------------------------------------------------------------------------
#endif

        WRITE(11,*)
        WRITE(11,*) '  t_loss[s]'
        WRITE(11,*) '-----------------------------------'

        x =  rstart(n)
        y =  COS(stalam(n) * pi)

        fmt  = '(ES12.4, A22)'

        IF((loss == 0) .AND. (trap == 0))THEN
          WRITE(11,fmt) loss_time, "    'untrapped unloss'"
          WRITE(88,'(13(ES12.4, A), I3, A)') rstart(n), ',', pstart(n), ',', zstart(n), ',', stalam(n), ',', x, ',', y, ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', loss_time, ',', 1, ','
        ELSE IF((loss == 0) .AND. (trap == 1))THEN
          WRITE(11,fmt) loss_time, "    'trapped unloss'"
          WRITE(88,'(13(ES12.4, A), I3, A)') rstart(n), ',', pstart(n), ',', zstart(n), ',', stalam(n), ',', 0.0, ',',  0.0, ',', x, ',', y, ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', loss_time, ',', 2, ','
        ELSE IF((loss == 1) .and. (trap == 0))THEN
          loss_n_p =  loss_n_p + 1
          WRITE(11,fmt) loss_time, "    'untrapped loss'"
          WRITE(88,'(13(ES12.4, A), I3, A)') rstart(n), ',', pstart(n), ',', zstart(n), ',', stalam(n), ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', x, ',', y, ',', 0.0, ',',  0.0, ',', loss_time, ',', 3, ','
        ELSE IF((loss == 1) .and. (trap == 1))THEN
          loss_n_t =  loss_n_t + 1 
          WRITE(11,fmt) loss_time, "    'trapped loss'"
          WRITE(88,'(13(ES12.4, A), I3, A)') rstart(n), ',', pstart(n), ',', zstart(n), ',', stalam(n), ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', 0.0, ',',  0.0, ',', x, ',', y, ',', loss_time, ',', 4, ','
        END IF

        CLOSE(88)
        CLOSE(50)
        CLOSE(11)

      END DO

      CALL MPI_REDUCE(loss_n_p, loss_t_p, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(loss_n_t, loss_t_t, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      loss_tol =  loss_t_p + loss_t_t

      fmt  = '(A40, I12, A, I12)'

      IF(myrank == 0)THEN
        PRINT *
        PRINT fmt, ' Number of Prompt Loss Particle  = ', loss_t_p, '/', nparticles
        PRINT fmt, ' Number of Trapped Loss Particle = ', loss_t_t, '/', nparticles
        PRINT fmt, ' Total Number of Lost Particle   = ', loss_tol, '/', nparticles
        PRINT *
      END IF

      DEALLOCATE(rstart, zstart, pstart, wk_in, stalam)


  END SELECT


END SUBROUTINE multi_particles
