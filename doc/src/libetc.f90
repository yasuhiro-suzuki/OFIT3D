FUNCTION upsilon_cal (x & !(in)
  &                  )

  IMPLICIT NONE
!Arguments
  REAL(8), INTENT(IN)  :: x
!Local variables
  REAL(8), PARAMETER :: pi = 3.141592653589793238462643383279502884197_8
  REAL(8) :: upsilon_cal, &
    &        integral
  REAL(8), EXTERNAL :: funphi


!----------------------------------------
  CALL qsimp(funphi, 0.0_8, x, integral)
!----------------------------------------

  upsilon_cal =  2.0_8 * integral / SQRT( pi )

 
  RETURN
END FUNCTION upsilon_cal

FUNCTION xsi_cal (x & !(in)
  &              )

  IMPLICIT NONE
!Arguments
  REAL(8), INTENT(IN)  :: x
!Local variables
  REAL(8), PARAMETER :: pi = 3.141592653589793238462643383279502884197_8
  REAL(8) :: xsi_cal, &
    &        upsilon, &
    &        dupsilon
  REAL(8), EXTERNAL :: upsilon_cal


  upsilon  =  upsilon_cal(x)
  dupsilon =  2.0_8 / (SQRT(pi) * EXP(x**2))
  xsi_cal  = (upsilon - x * dupsilon) / (2.0_8 * x**2)
 

  RETURN
END FUNCTION xsi_cal

FUNCTION funphi (t & !(in)
  &             )

  IMPLICIT NONE
!Arguments
  REAL(8), INTENT(IN) :: t
!Local varibales
  REAL(8) :: funphi


  funphi = 1.0_8 / EXP(t**2)


  RETURN
END FUNCTION funphi

SUBROUTINE qsimp (func, a, b, & ! (in)
  &               s           & !(out)
  &              )

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: a, &
    &                    b
  REAL(8), INTENT(OUT) :: s

  INTEGER, PARAMETER :: jmax =  20
  INTEGER :: j
  REAL(8), PARAMETER :: eps =  1.0e-06_8
  REAL(8) ::  os,  &
    &         ost, &
    &         st
  REAL(8), EXTERNAL :: func


  ost = -1.0e+30_8
  os  = -1.0e+30_8

  loop010 : DO j=1,jmax
    CALL trapzd(func, a, b, st, j)
    s = (4.0_8 * st - ost) / 3.0_8
    IF(j > 5)then
      IF((ABS(s - os) < eps * ABS(os)) .or. ((s == 0.0_8) .and. (os == 0.0_8 ))) RETURN
    END IF
    os  =  s
    ost =  st
  END DO loop010

  STOP ' too many steps in qsimp '


END SUBROUTINE qsimp

SUBROUTINE trapzd(func, a, b, s, n)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(8), INTENT(IN) :: a, &
    &                    b
  REAL(8), INTENT(OUT) :: s

  INTEGER :: it, &
    &        j
  REAL(8) :: del, &
    &        sum, &
    &        tnm, &
    &        x
  REAL(8), EXTERNAL :: func


  IF(n == 1)THEN
    s =  0.5_8 * (b - a) * (func(a) + func(b))
  ELSE
    it  =  2**(n - 2)
    tnm =  it
    del = (b - a) / tnm
    x   =  a + 0.5_8 * del
    sum =  0.0_8
    loop010 : DO j=1,it
      sum =  sum + func(x)
      x   =  x + del
    END DO loop010
    s =  0.5_8 * (s + (b - a) * sum / tnm)
  END IF


  RETURN
END SUBROUTINE trapzd
