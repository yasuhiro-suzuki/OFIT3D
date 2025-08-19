!> @file divertor_mod.f90
!------------------------------------------------------------------------------
!
! MODULE: divertor_mod
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
MODULE divertor_mod

  USE kind_spec
  USE param1,   ONLY : pi2

  IMPLICIT NONE

  PRIVATE

  LOGICAL :: ldivertor      = .false., &
    &        ldivertor_vtk  = .false.
  INTEGER :: num_nodes      =  18, &
    &        num_elements   =  20
  CHARACTER(LEN=10) :: divertor_model = ''
  CHARACTER(LEN=300) :: divertor_file = '', &
    &                   nodes_file    = '', &
    &                   elements_file = ''

  INTEGER, ALLOCATABLE :: def_elements(:,:)
  REAL(RP), ALLOCATABLE :: x_nodes(:), &
    &                      y_nodes(:), &
    &                      z_nodes(:)

  PUBLIC :: ldivertor,         &
    &       ldivertor_vtk,     &
    &       divertor_model,    &
    &       num_nodes,         &
    &       num_elements,      &
    &       divertor_file,     &
    &       nodes_file,        &
    &       elements_file,     &
    &       read_divertor,     &
    &       free_mem_divertor, &
    &       divertor_loss

CONTAINS

  SUBROUTINE make_mem_divertor

    IMPLICIT NONE


    ALLOCATE(def_elements(num_elements,3), x_nodes(num_nodes), y_nodes(num_nodes), z_nodes(num_nodes))


  END SUBROUTINE make_mem_divertor

  SUBROUTINE free_mem_divertor

    IMPLICIT NONE


    DEALLOCATE(def_elements, x_nodes, y_nodes, z_nodes)


  END SUBROUTINE free_mem_divertor

  SUBROUTINE read_divertor

    IMPLICIT NONE


    SELECT CASE(TRIM(divertor_model))
      CASE('w7code')
        CALL read_divertor_w7code
      CASE DEFAULT
        CALL make_mem_divertor
        CALL read_nodes
        CALL read_elements
    END SELECT

    IF(ldivertor_vtk) CALL write_divertor_vtk


  END SUBROUTINE read_divertor

  SUBROUTINE read_nodes

    IMPLICIT NONE

    INTEGER :: itemp, &
      &        i


    DO i=1,num_nodes
      READ(32,*) itemp, x_nodes(i), y_nodes(i), z_nodes(i)
    END DO


  END SUBROUTINE read_nodes

  SUBROUTINE read_elements

    IMPLICIT NONE

    INTEGER :: itemp, &
      &        i


    DO i=1,num_elements
      READ(33,*) itemp, def_elements(i,1), def_elements(i,2), def_elements(i,3)
    END DO

    def_elements(1:num_elements,1:3) =  def_elements(1:num_elements,1:3) - 1


  END SUBROUTINE read_elements

  SUBROUTINE read_divertor_w7code

    IMPLICIT NONE

     INTEGER :: mtor, &
       &        nu,   &
       &        nv
!Local variables
    INTEGER :: i1, &
      &        i2, &
      &        i3, &
      &        i,  &
      &        j,  &
      &        k
    REAL(RP) :: r,   &
      &         phi, &
      &         z


    READ(31,*)
    READ(31,*) nv, nu, mtor

    num_nodes      =  nv * nu
    num_elements   =  2 * (nv - 1) * (nu - 1)

    CALL make_mem_divertor

    k =  0
    DO j=1,nv
      READ(31,*) phi
      phi =  phi * pi2 / 360.0_RP
      DO i=1,nu
        k =  k + 1
        READ(31,*) r, z
        x_nodes(k) =  1.0E-02_RP * r * COS(phi)
        y_nodes(k) =  1.0E-02_RP * r * SIN(phi)
        z_nodes(k) =  1.0E-02_RP * z
      END DO
    END DO

    k =  0
    DO j=1,nv-1
      DO i=1,nu-1
        i1 =  i  + (j - 1) * nu
        i2 =  i1 + 1
        i3 =  i  + j * nu
        k  =  k + 1
        def_elements(k,1) =  i1 - 1
        def_elements(k,2) =  i2 - 1
        def_elements(k,3) =  i3 - 1
        i1 =  i  + (j - 1) * nu + 1
        i2 =  i  + j * nu
        i3 =  i2 + 1
        k  =  k + 1
        def_elements(k,1) =  i1 - 1
        def_elements(k,2) =  i2 - 1
        def_elements(k,3) =  i3 - 1
      END DO
    END DO


  END SUBROUTINE read_divertor_w7code

  SUBROUTINE write_divertor_vtk

    IMPLICIT NONE
!Local variables
    INTEGER :: i1, &
      &        i2, &
      &        i3, &
      &        i


    OPEN(301, FILE='divertor.vtk', FORM='formatted', STATUS='unknown')

    WRITE(301,'(A)') '# vtk DataFile Version 3.0'
    WRITE(301,'(A)') 'Unstructured Grid'
    WRITE(301,'(A)') 'ASCII'
    WRITE(301,*)
    WRITE(301,'(A)') 'DATASET UNSTRUCTURED_GRID'
    WRITE(301,'(A,I12,A)') 'POINTS ', num_nodes, ' float'

    DO i=1,num_nodes
      WRITE(301,'(3(ES12.4))') x_nodes(i), y_nodes(i), z_nodes(i)
    END DO

    WRITE(301,*)
    WRITE(301,'(A,2I12)') 'CELLS', num_elements, 4 * num_elements

    DO i=1,num_elements
     WRITE(301,'(I3,3I12)') 3, def_elements(i,1), def_elements(i,2), def_elements(i,3)
    END DO

    WRITE(301,*)
    WRITE(301,'(A,I12)') 'CELL_TYPES', num_elements

    DO i=1,num_elements
      WRITE(301,'(I2)') 5
    END DO

    WRITE(301,'(A,I12)') 'POINT_DATA', num_nodes
    WRITE(301,'(A)') 'SCALARS scalars float 1'
    WRITE(301,'(A)') 'LOOKUP_TABLE default'

    DO i=1,num_nodes
      WRITE(301,'(ES12.4)') 5.0_RP
    END DO

    CLOSE(301)


  END SUBROUTINE write_divertor_vtk

  SUBROUTINE divertor_loss (r0, phi0, z0, r1, phi1, z1, &!(in)
    &                       intersect, beta, point      &!(out)
    &                      )

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

    INTEGER  :: i1,       &
      &         i2,       &
      &         i3,       &
      &         i,        &
      &         j,        &
      &         k
    REAL(RP) :: x0,       &
      &         y0,       &
      &         x1,       &
      &         y1,       &
      &         p0,       &
      &         p1,       &
      &         dis,      &
      &         dis_ps,   &
      &         dis_pe,   &
      &         dis_12,   &
      &         dis_13,   &
      &         nn(3),    &
      &         a(3),     &
      &         b(3),     &
      &         c(3),     &
      &         ab(3),    &
      &         ac(3),    &
      &         bc(3),    &
      &         ca(3),    &
      &         ps(3),    &
      &         pe(3),    &
      &         aps(3),   &
      &         ape(3),   &
      &         ap(3),    &
      &         bp(3),    &
      &         cp(3),    &
      &         abxap(3), &
      &         bcxbp(3), &
      &         caxcp(3)


    intersect =  0

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

    DO i=1,num_elements
      i1       =  def_elements(i,1) + 1
      i2       =  def_elements(i,2) + 1
      i3       =  def_elements(i,3) + 1
      a(1)     =  x_nodes(i1)
      a(2)     =  y_nodes(i1)
      a(3)     =  z_nodes(i1)
      b(1)     =  x_nodes(i2)
      b(2)     =  y_nodes(i2)
      b(3)     =  z_nodes(i2)
      c(1)     =  x_nodes(i3)
      c(2)     =  y_nodes(i3)
      c(3)     =  z_nodes(i3)

      ab(1:3)  =  b(1:3) - a(1:3)
      ac(1:3)  =  c(1:3) - a(1:3)

      aps(1:3) =  ps(1:3) - a(1:3)
      ape(1:3) =  pe(1:3) - a(1:3)

      nn(1)    =  ab(2) * ac(3) - ab(3) * ac(2)
      nn(2)    =  ab(3) * ac(1) - ab(1) * ac(3)
      nn(3)    =  ab(1) * ac(2) - ab(2) * ac(1)

      dis      =  SQRT(nn(1)**2 + nn(2)**2 + nn(3)**2)
      nn(1:3)  =  nn(1:3) / dis

      dis_ps   =  DOT_PRODUCT(aps(1:3), nn(1:3))
      dis_pe   =  DOT_PRODUCT(ape(1:3), nn(1:3))
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

        IF((dis_12 >= 0.0_RP) .AND. (dis_13 >= 0.0_RP))THEN
          intersect = 1
          RETURN
        END IF
      ELSE
        intersect = 0
      END IF

    END DO


  END SUBROUTINE divertor_loss

END MODULE divertor_mod
