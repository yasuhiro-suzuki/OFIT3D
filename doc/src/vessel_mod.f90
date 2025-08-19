!> @file vessel_mod.f90
!------------------------------------------------------------------------------
!
! MODULE: vessel_mod
!
!> @author
!> Yasuhiro Suzuki, National Institute for Fusion Science
!
! DESCRIPTION:
!> @brief
!>
!
! REVISION HISTORY:
!> @date 19 Apr 2020
!
!> @version Initial Version
!
!------------------------------------------------------------------------------
MODULE vessel_mod

  USE kind_spec
  USE param1,   ONLY : pi, &
    &                  pi2

  IMPLICIT NONE

  PRIVATE

  LOGICAL :: lvessel,       & !<
    &        lcheck_vessel, &
    &        lintersect,    &
    &        lvessel_vtk,   &
    &        lvessel_txt
  INTEGER :: mtor, &
    &        nu,   &
    &        nv
  CHARACTER(LEN=10) :: vessel_model = ''
  CHARACTER(LEN=300) :: vessel_file = 'vessel.dat'

  REAL(RP) :: pi2m, &
    &         dtor
  REAL(RP), ALLOCATABLE :: xx(:,:), &
    &                      yy(:,:), &
    &                      zz(:,:), &
    &                      rr(:,:)

  PUBLIC :: lvessel,         &
    &       lcheck_vessel,   &
    &       lintersect,      &
    &       lvessel_vtk,     &
    &       vessel_model,    &
    &       vessel_file,     &
    &       nu,              &
    &       nv,              &
    &       read_vessel,     &
    &       free_mem_vessel, &
    &       vessel,          &
    &       check_vessel,    &
    &       vessel_loss

CONTAINS

  SUBROUTINE make_mem_vessel

    IMPLICIT NONE


    ALLOCATE(xx(nu,nv), yy(nu,nv), zz(nu,nv), rr(nu,nv))


  END SUBROUTINE make_mem_vessel

  SUBROUTINE free_mem_vessel

    IMPLICIT NONE


    DEALLOCATE(xx, yy, zz, rr)


  END SUBROUTINE free_mem_vessel

  SUBROUTINE read_vessel

    IMPLICIT NONE


    SELECT CASE(TRIM(vessel_model))
      CASE("2d", "2D")
        CALL read_vessel_2d
      CASE("old", "OLD")
        CALL read_vessel_old
      CASE("ana", "ANA")
        CALL read_vessel_ana
      CASE DEFAULT
        CALL read_vessel_3d
    END SELECT

    IF(lvessel_vtk) CALL write_vessel_vtk


  END SUBROUTINE read_vessel

  SUBROUTINE read_vessel_3d

    IMPLICIT NONE
!Local variables
    INTEGER :: i, &
      &        j
    REAL(RP) :: phi


    READ(30,*) mtor, nu, nv

    pi2m =  pi2  / mtor
    dtor =  pi2m / nv

    nv   =  nv + 1

    CALL make_mem_vessel

    DO j=1,nv-1
      DO i=1,nu
        READ(30,*) rr(i,j), zz(i,j)
      END DO
    END DO

    DO i=1,nu
      rr(i,nv) =  rr(i,1)
      zz(i,nv) =  zz(i,1)
    END DO

    DO j=1,nv
      phi =  dtor * (j - 1)
      DO i=1,nu
        xx(i,j) =  rr(i,j) * COS(phi)
        yy(i,j) =  rr(i,j) * SIN(phi)
      END DO
    END DO


  END SUBROUTINE read_vessel_3d

  SUBROUTINE read_vessel_3d_old

    IMPLICIT NONE
!Local variables
    INTEGER :: i, &
      &        j
    REAL(RP) :: phi


    READ(30,*) mtor, nu, nv

    pi2m =  pi2  / mtor
    dtor =  pi2m / (nv - 1)

    CALL make_mem_vessel

    DO j=1,nv
      DO i=1,nu
        READ(30,*) rr(i,j), zz(i,j)
      END DO
    END DO

    DO j=1,nv
      phi =  dtor * (j - 1)
      DO i=1,nu
        xx(i,j) =  rr(i,j) * COS(phi)
        yy(i,j) =  rr(i,j) * SIN(phi)
      END DO
    END DO


  END SUBROUTINE read_vessel_3d_old

  SUBROUTINE read_vessel_old

    IMPLICIT NONE
!Local variables
    INTEGER :: i, &
      &        j


    READ(30,*) mtor, nu, nv

    pi2m =  pi2  / mtor
    dtor =  pi2m / (nv - 1)

    CALL make_mem_vessel

    DO j=1,nv
      DO i=1,nu
        READ(30,*) xx(i,j), yy(i,j), zz(i,j), rr(i,j)
      END DO
    END DO


  END SUBROUTINE read_vessel_old

  SUBROUTINE read_vessel_2d

    IMPLICIT NONE
!Local variables
    INTEGER :: i, &
      &        j
    REAL(RP) :: phi
    REAL(RP), ALLOCATABLE :: r1(:), &
      &                      z1(:)

    READ(30,*) nu

    ALLOCATE(r1(nu), z1(nu))

    mtor =  1
    nv   =  120

    pi2m =  pi2  / mtor
    dtor =  pi2m / (nv - 1)

    CALL make_mem_vessel

    DO i=1,nu
      READ(30,*) r1(i), z1(i)
    END DO

    DO j=1,nv
      phi =  dtor * (j - 1)
      DO i=1,nu
        xx(i,j) =  r1(i) * COS(phi)
        yy(i,j) =  r1(i) * SIN(phi)
        zz(i,j) =  z1(i)
        rr(i,j) =  r1(i)
      END DO
    END DO

    DEALLOCATE(r1, z1)


  END SUBROUTINE read_vessel_2d

  SUBROUTINE read_vessel_ana

    IMPLICIT NONE
!Local variables
    INTEGER :: i, &
      &        j
    REAL(RP), PARAMETER :: r0 =  3.0_RP, &
      &                    a0 =  1.2_RP
    REAL(RP) :: theta,  &
      &         dtheta, &
      &         phi,    &
      &         r,      &
      &         z



    mtor   =  1
    nv     =  121
    nu     =  361

    pi2m   =  pi2  / mtor
    dtor   =  pi2m / (nv - 1)
    dtheta =  pi2 / (nu - 1)

    CALL make_mem_vessel

    DO j=1,nv
      phi =  dtor * (j - 1)
      DO i=1,nu
        theta   =  dtheta * (i - 1)
        r       =  r0 + a0 * COS(theta)
        z       =       a0 * SIN(theta)
        xx(i,j) =  r * COS(phi)
        yy(i,j) =  r * SIN(phi)
        zz(i,j) =  z
        rr(i,j) =  r
      END DO
    END DO


  END SUBROUTINE read_vessel_ana

  SUBROUTINE write_vessel_vtk

    IMPLICIT NONE
!Local variables
    INTEGER :: i1, &
      &        i2, &
      &        i3, &
      &        i,  &
      &        j


    OPEN(300, FILE='vessel.vtk', FORM='formatted', STATUS='unknown')

    WRITE(300,'(A)') '# vtk DataFile Version 3.0'
    WRITE(300,'(A)') 'Unstructured Grid'
    WRITE(300,'(A)') 'ASCII'
    WRITE(300,*)
    WRITE(300,'(A)') 'DATASET UNSTRUCTURED_GRID'
    WRITE(300,'(A,I12,A)') 'POINTS ', nu * nv, ' double'

    DO j=1,nv
      DO i=1,nu
        WRITE(300,'(3(ES15.7))') xx(i,j), yy(i,j), zz(i,j)
      END DO
    END DO

    WRITE(300,*)
    WRITE(300,'(A,2I12)') 'CELLS', 2 * (nv - 1) * (nu - 1), 6 * (nv - 1) * (nu - 1) + 2 * (nv - 1) * (nu - 1)

    DO j=1,nv-1
      DO i=1,nu-1
        i1 =  i  + (j - 1) * nu 
        i2 =  i1 + 1
        i3 =  i  + j * nu
        WRITE(300,'(I3,3I12)') 3, i1 - 1, i2 - 1, i3 - 1
        i1 =  i  + (j - 1) * nu +1
        i2 =  i  + j * nu
        i3 =  i2 + 1 
        WRITE(300,'(I3,3I12)') 3, i1 - 1, i2 - 1, i3 - 1
      END DO
    END DO

    WRITE(300,*)
    WRITE(300,'(A,I12)') 'CELL_TYPES', 2 * (nv - 1) * (nu - 1)

    DO j=1,nv-1
      DO i=1,nu-1
        WRITE(300,'(I2)') 5
        WRITE(300,'(I2)') 5
      END DO
    END DO

    WRITE(300,'(A,I12)') 'POINT_DATA', nv * nu
    WRITE(300,'(A)') 'SCALARS scalars float 1'
    WRITE(300,'(A)') 'LOOKUP_TABLE default'

    DO j=1,nv
      DO i=1,nu
        WRITE(300,'(ES15.7)') 1.0_RP
      END DO
    END DO

    CLOSE(300)


  END SUBROUTINE write_vessel_vtk

  SUBROUTINE write_vessel_txt

    IMPLICIT NONE
!Local variables
    INTEGER :: i, &
      &        j


    DO j=1,nv-1
      DO i=1,nu-1
        WRITE(200,'(3ES15.7)') xx(i,j), yy(i,j), zz(i,j)
        WRITE(200,'(3ES15.7)') xx(i+1,j), yy(i+1,j), zz(i+1,j)
        WRITE(200,'(3ES15.7)') xx(i,j+1), yy(i,j+1), zz(i,j+1)
        WRITE(200,'(3ES15.7)') xx(i,j), yy(i,j), zz(i,j)
        WRITE(200,*)
        WRITE(200,'(3ES15.7)') xx(i+1,j), yy(i+1,j), zz(i+1,j)
        WRITE(200,'(3ES15.7)') xx(i,j+1), yy(i,j+1), zz(i,j+1)
        WRITE(200,'(3ES15.7)') xx(i+1,j+1), yy(i+1,j+1), zz(i+1,j+1)
        WRITE(200,'(3ES15.7)') xx(i+1,j), yy(i+1,j), zz(i+1,j)
        WRITE(200,*)
      END DO
      WRITE(200,*)
      WRITE(200,*)
    END DO


  END SUBROUTINE write_vessel_txt
 
  SUBROUTINE vessel (phi,   & !(in)
    &                r1, z1 & !(out)
    &               )


    IMPLICIT NONE
!Arguments
    REAL(RP), INTENT(IN) :: phi
    REAL(RP), INTENT(OUT) :: r1(nu), &
      &                      z1(nu)
!Local variables
    INTEGER :: iphi
    REAL(RP) :: phi0,  &
      &         phi1,  &
      &         dtor1, &
      &         alpha


    iphi  =  phi / pi2m
    phi1  =  phi - pi2m * iphi
    IF(phi1 <  0.0_RP) phi1 =  phi1 + pi2m
    IF(phi1 >= pi2m)   phi1 =  phi1 - pi2m

    iphi  =  phi1 / dtor
    phi0  =  dtor * iphi
    dtor1 =  phi1 - phi0

    iphi  =  iphi + 1

    alpha =  dtor1 / dtor

    r1(:) =  rr(:,iphi) + alpha * (rr(:,iphi+1) - rr(:,iphi))
    z1(:) =  zz(:,iphi) + alpha * (zz(:,iphi+1) - zz(:,iphi))


  END SUBROUTINE vessel

  SUBROUTINE check_vessel (r, phi, z,     & !(in)
    &                      iout           & !(out)
    &                     )

    IMPLICIT NONE
!Arguments
    INTEGER, INTENT(OUT) :: iout
    REAL(RP), INTENT(IN) :: r,   &
      &                     phi, &
      &                     z
!Local variables
    INTEGER :: iphi,   &
      &        itheta, &
      &        i
    REAL(RP) :: theta,  &
      &         theta0, &
      &         theta1, &
      &         dtheta, &
      &         r1(nu), &
      &         z1(nu)


    CALL vessel(phi, r1, z1)

    theta  =  0.0_RP
    theta0 =  ATAN2(z1(1) - z, r1(1) - r)
    IF(theta0 < 0.0_RP) theta0 =  theta0 + pi2
    DO i=2,nu
      theta1 =  ATAN2(z1(i) - z, r1(i) - r)
      IF(theta1 < 0.0_RP) theta1 =  theta1 + pi2
      dtheta =  theta1 - theta0
      IF(dtheta >  pi) dtheta =  dtheta - pi2
      IF(dtheta < -pi) dtheta =  dtheta + pi2
      theta  =  theta + dtheta
      theta0 =  theta1
    END DO
    theta =  theta / pi2
    IF(theta >= 0.0_RP)THEN
      theta =  theta + 1.0e-06_RP
    ELSE
      theta =  theta - 1.0e-06_RP
    END IF
    itheta =  theta

    iout   =  0
    IF(itheta == 0) iout =  1

 
  END SUBROUTINE check_vessel

  SUBROUTINE vessel_loss (r0, phi0, z0, r1, phi1, z1, &!(in)
    &                     intersect, beta, point      &!(out)
    &                    )

    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: intersect
    REAL(RP), INTENT(IN) :: r0,       &
      &                     phi0,     &
      &                     z0,       &
      &                     r1,       &
      &                     phi1,     &
      &                     z1
    REAL(RP), INTENT(OUT) :: beta,    &
      &                      point(3)

    INTEGER  :: i,        &
     &          j,        &
     &          k
    INTEGER  :: iphi0,    &
     &          iphi1
    REAL(RP) :: x0,       &
     &          y0,       &
     &          x1,       &
     &          y1,       &
     &          p0,       &
     &          p1,       &
     &          dis,      &
     &          dis_ps,   &
     &          dis_pe,   &
     &          dis_12,   &
     &          dis_13,   &
     &          nn(3),    &
     &          a(3),     &
     &          b(3),     &
     &          c(3),     &
     &          ab(3),    &
     &          ac(3),    &
     &          bc(3),    &
     &          ca(3),    &
     &          ps(3),    &
     &          pe(3),    &
     &          aps(3),   &
     &          ape(3),   &
     &          ap(3),    &
     &          bp(3),    &
     &          cp(3),    &
     &          abxap(3), &
     &          bcxbp(3), &
     &          caxcp(3)


    intersect =  0

    iphi0     =  phi0 / pi2m
    iphi1     =  phi1 / pi2m
    p0        =  phi0 - pi2m * iphi0
    p1        =  phi1 - pi2m * iphi1

    IF(p0 <  0.0_RP) p0 =  p0 + pi2m
    IF(p0 >= pi2m)   p0 =  p0 - pi2m
    IF(p1 <  0.0_RP) p1 =  p1 + pi2m
    IF(p1 >= pi2m)   p1 =  p1 - pi2m

    iphi0 =  p0 / dtor + 1
    iphi1 =  p1 / dtor + 1

    x0    =  r0 * COS(p0)
    y0    =  r0 * SIN(p0)
    x1    =  r1 * COS(p1)
    y1    =  r1 * SIN(p1)

    ps(1) =  x0
    ps(2) =  y0
    ps(3) =  z0
    pe(1) =  x1
    pe(2) =  y1
    pe(3) =  z1

    DO k=1,3

      IF(k == 1)THEN
        j =  iphi0
      ELSE IF(k == 2)THEN
        j =  iphi0 + 1
        IF(j == nv)THEN
          j =  1
        END IF
      ELSE
        j =  iphi0 - 1
        IF(j == 0)THEN
             j =  nv - 1
        END IF
      END IF

      DO i = 1, nu - 1
        a(1)     =  xx(i,j)
        a(2)     =  yy(i,j)
        a(3)     =  zz(i,j)
        b(1)     =  xx(i+1,j)
        b(2)     =  yy(i+1,j)
        b(3)     =  zz(i+1,j)
        c(1)     =  xx(i,j+1)
        c(2)     =  yy(i,j+1)
        c(3)     =  zz(i,j+1)

        ab(1:3)  =  b(1:3) - a(1:3)
        ac(1:3)  =  c(1:3) - a(1:3)

        aps(1:3) =  ps(1:3) - a(1:3)
        ape(1:3) =  pe(1:3) - a(1:3)

        nn(1)    =  ab(2) * ac(3) - ab(3) * ac(2)
        nn(2)    =  ab(3) * ac(1) - ab(1) * ac(3)
        nn(3)    =  ab(1) * ac(2) - ab(2) * ac(1)

        dis      =  SQRT(nn(1)**2 + nn(2)**2 + nn(3)**2)
        nn(1:3)  =  nn(1:3) / dis

        dis_ps   =  DOT_PRODUCT(aps(:), nn(:))
        dis_pe   =  DOT_PRODUCT(ape(:), nn(:))
        IF(dis_ps * dis_pe <= 0.0_RP)THEN
          beta       =  ABS(dis_ps) / (ABS(dis_ps) + ABS(dis_pe))
          point(1:3) =  ps(1:3) + beta * (pe(1:3) - ps(1:3))
          ap(1:3)    =  point(1:3) - a(1:3)
          bp(1:3)    =  point(1:3) - b(1:3)
          cp(1:3)    =  point(1:3) - c(1:3)
          ab(1:3)    =  b(1:3) - a(1:3)
          bc(1:3)    =  c(1:3) - b(1:3)
          ca(1:3)    =  a(1:3) - c(1:3)

          abxap(1)   =  ab(2) * ap(3) - ab(3) * ap(2)
          abxap(2)   =  ab(3) * ap(1) - ab(1) * ap(3)
          abxap(3)   =  ab(1) * ap(2) - ab(2) * ap(1)

          bcxbp(1)   =  bc(2) * bp(3) - bc(3) * bp(2)
          bcxbp(2)   =  bc(3) * bp(1) - bc(1) * bp(3)
          bcxbp(3)   =  bc(1) * bp(2) - bc(2) * bp(1)

          caxcp(1)   =  ca(2) * cp(3) - ca(3) * cp(2)
          caxcp(2)   =  ca(3) * cp(1) - ca(1) * cp(3)
          caxcp(3)   =  ca(1) * cp(2) - ca(2) * cp(1)

          dis_12     =  abxap(1) * bcxbp(1) + abxap(2) * bcxbp(2) + abxap(3) * bcxbp(3)
          dis_13     =  abxap(1) * caxcp(1) + abxap(2) * caxcp(2) + abxap(3) * caxcp(3)

          IF(dis_12 >= 0.0_RP .AND. dis_13 >= 0.0_RP)THEN
            intersect = 1
            RETURN
          END IF
        ELSE
          intersect = 0
        END IF

        a(1)     =  xx(i+1,j)
        a(2)     =  yy(i+1,j)
        a(3)     =  zz(i+1,j)
        b(1)     =  xx(i,j+1)
        b(2)     =  yy(i,j+1)
        b(3)     =  zz(i,j+1)
        c(1)     =  xx(i+1,j+1)
        c(2)     =  yy(i+1,j+1)
        c(3)     =  zz(i+1,j+1)

        ab(1:3)  =  b(1:3) - a(1:3)
        ac(1:3)  =  c(1:3) - a(1:3)

        aps(1:3) =  ps(1:3) - a(1:3)
        ape(1:3) =  pe(1:3) - a(1:3)

        nn(1)    =  ab(2) * ac(3) - ab(3) * ac(2)
        nn(2)    =  ab(3) * ac(1) - ab(1) * ac(3)
        nn(3)    =  ab(1) * ac(2) - ab(2) * ac(1)

        dis      =  SQRT(nn(1)**2 + nn(2)**2 + nn(3)**2)
        nn(1:3)  =  nn(1:3) / dis

        dis_ps   =  DOT_PRODUCT(aps(:), nn(:))
        dis_pe   =  DOT_PRODUCT(ape(:), nn(:))
        IF(dis_ps * dis_pe <= 0.0_RP)THEN
          beta       =  ABS(dis_ps) / (ABS(dis_ps) + ABS(dis_pe))
          point(1:3) =  ps(1:3) + beta * (pe(1:3) - ps(1:3))
          ap(1:3)    =  point(1:3) - a(1:3)
          bp(1:3)    =  point(1:3) - b(1:3)
          cp(1:3)    =  point(1:3) - c(1:3)
          ab(1:3)    =  b(1:3) - a(1:3)
          bc(1:3)    =  c(1:3) - b(1:3)
          ca(1:3)    =  a(1:3) - c(1:3)

          abxap(1)   =  ab(2) * ap(3) - ab(3) * ap(2)
          abxap(2)   =  ab(3) * ap(1) - ab(1) * ap(3)
          abxap(3)   =  ab(1) * ap(2) - ab(2) * ap(1)

          bcxbp(1)   =  bc(2) * bp(3) - bc(3) * bp(2)
          bcxbp(2)   =  bc(3) * bp(1) - bc(1) * bp(3)
          bcxbp(3)   =  bc(1) * bp(2) - bc(2) * bp(1)

          caxcp(1)   =  ca(2) * cp(3) - ca(3) * cp(2)
          caxcp(2)   =  ca(3) * cp(1) - ca(1) * cp(3)
          caxcp(3)   =  ca(1) * cp(2) - ca(2) * cp(1)

          dis_12     =  abxap(1) * bcxbp(1) + abxap(2) * bcxbp(2) + abxap(3) * bcxbp(3)
          dis_13     =  abxap(1) * caxcp(1) + abxap(2) * caxcp(2) + abxap(3) * caxcp(3)

          IF(dis_12 >= 0.0_RP .AND. dis_13 >= 0.0_RP)THEN
            intersect = 1
            RETURN
          END IF
        ELSE
          intersect = 0
        END IF

      END DO
    END DO  


  END SUBROUTINE vessel_loss

END MODULE vessel_mod
