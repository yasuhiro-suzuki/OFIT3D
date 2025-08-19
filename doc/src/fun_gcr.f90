!=fun_gcr.f90
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
SUBROUTINE fun_gcr (t, y,    & ! (in)
  &                 yp, iout & ! (out)
  &                )

  USE kind_spec
  USE param1,                ONLY : qa,    &
    &                               ma
  USE inv_val_mod,           ONLY : mu
  USE cylindrical_coord_mod, ONLY : mgval2

  IMPLICIT NONE
!Arguments
  INTEGER, INTENT(OUT) :: iout
  REAL(RP), INTENT(IN) :: t,  &
    &                     y(4)
  REAL(RP), INTENT(OUT) :: yp(4)
!Local variables
  REAL(RP) :: vvv,      &
    &         bmod,     &
    &         dbdr,     &
    &         dbdp,     &
    &         dbdz,     &
    &         rotbr,    &
    &         rotbp,    &
    &         rotbz,    &
    &         bgradb,   &
    &         rbgrdb,   &
    &         efac,     &
    &         delta,    &
    &         vpara,    &
    &         vpara2,   &
    &         ombinv,   &
    &         eps,      &
    &         aa1,      &
    &         aa2,      &
    &         aa3,      &
    &         bvec(4),  &
    &         dbvdr(4), &
    &         dbvdp(4), &
    &         dbvdz(4), &
    &         s,        &
    &         grads(3)


  iout =  0

  CALL mgval2(y(1), y(2), y(3), bvec(:), dbvdr(:), dbvdp(:), dbvdz(:))

  bmod =  bvec(4)

  IF(bmod == 0.0_RP)THEN
    yp(:) =  0.0_RP
    iout  =  1
    RETURN
  END IF

  dbdr   = dbvdr(4)
  dbdp   = dbvdp(4) / y(1)
  dbdz   = dbvdz(4)

  rotbr  =  dbvdp(3) / y(1) - dbvdz(2)
  rotbp  =  dbvdz(1) - dbvdr(3)
  rotbz  =  dbvdr(2) - dbvdp(1) / y(1) + bvec(2) / y(1)

  delta  =  bvec(1) * rotbr + bvec(2) * rotbp + bvec(3) * rotbz
  vpara  =  y(4) 
  vpara2 =  y(4) * y(4)

  ombinv =  ma / (qa * bmod**2)
  aa1    =  vpara / bmod
  aa2    =  ombinv * (mu + vpara2 / bmod)
  aa3    =  ombinv * vpara2

  bgradb =  bvec(1) * dbdr + bvec(2) * dbdp + bvec(3) * dbdz
  rbgrdb =  rotbr   * dbdr + rotbp   * dbdp + rotbz   * dbdz ! corrected by J. Morimoto

  yp(1)  =  aa1 * bvec(1) + aa2 * (bvec(2) * dbdz - bvec(3) * dbdp) + aa3 * rotbr
  yp(2)  = (aa1 * bvec(2) + aa2 * (bvec(3) * dbdr - bvec(1) * dbdz) + aa3 * rotbp) / y(1)
  yp(3)  =  aa1 * bvec(3) + aa2 * (bvec(1) * dbdp - bvec(2) * dbdr) + aa3 * rotbz
  yp(4)  = -mu * (bgradb / bmod + rbgrdb * ombinv * y(4))

  eps    = (ombinv / bmod) * vpara * delta
  efac   =  1.0_RP / (1.0_RP + eps)
  yp(1)  =  efac * yp(1)  
  yp(2)  =  efac * yp(2)  
  yp(3)  =  efac * yp(3)
  yp(4)  =  efac * yp(4)


  RETURN
END SUBROUTINE fun_gcr
