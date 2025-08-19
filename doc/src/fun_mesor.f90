SUBROUTINE fun_mesor ( t, yy, & !(in)
  &                    yp     & !(out)
  &                  )

  USE kind_spec
  USE param1,                ONLY : qom_a
  USE cylindrical_coord_mod, ONLY : mgval1

  IMPLICIT NONE

!Arguments
  REAL(RP), INTENT(IN) :: t,    &
    &                     yy(6)
  REAL(RP), INTENT(OUT) :: yp(6)
!Local variable
  REAL(RP) :: vr,  &
    &         vp,  &
    &         vz,  &
    &         br,  &
    &         bp,  &
    &         bz,  &
    &         bb,  &
    &         r,   &
    &         phi, &
    &         z


  vr  =  yy(1)
  vp  =  yy(2)
  vz  =  yy(3)

  r   =  yy(4)
  phi =  yy(5)
  z   =  yy(6)

  CALL mgval1(r, phi, z, br, bp, bz, bb)

  IF(bb == 0.0_RP)THEN
    yp(1:6) =  0.0_RP 
    RETURN
  END IF

  yp(1) =  qom_a * (r * vp * bz -     vz * bp)     + r * vp**2
  yp(2) =  qom_a * (    vz * br -     vr * bz) / r - 2 * vr * vp / r
  yp(3) =  qom_a * (    vr * bp - r * vp * br)

  yp(4) =  vr
  yp(5) =  vp
  yp(6) =  vz


  RETURN
END SUBROUTINE fun_mesor
