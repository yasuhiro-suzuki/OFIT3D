SUBROUTINE rkg6 (x0, h, n, fun, & !(in)
  &              y0,            & !(inout)
  &              iout           & !(out)
  &             )

!
!    This subroutine is written by Y. Suzuki
!        at Graduate School of Energy Science (Kyoto Univ)
!         2002/12/28
!
!    Based program is wrtten by K. Hamamatsu
!        at Faculty of Science (Hiroshima Univ.)
!          1980/12/18
!
!  Runge-Kutta-Huta Formulas  ( Sixth order 8-stage )
!
!  <<< Reference >>>
!  " Improved Sixth-order Runge-kutta formulas and Approximate
!    Continuous Solution of Ordinary Differential Equation "
!  by D. Sarafyan:  J. Math. Anal. Appl. 40, 436-455 (1972)
!

  USE kind_spec

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(OUT) :: iout
  REAL(RP), INTENT(IN) :: x0, &
    &                     h
  REAL(RP), INTENT(INOUT) :: y0(n)

  REAL(RP) :: x1,  &
    &         c01, &
    &         c02, &
    &         c03, &
    &         c04, &
    &         c05, &
    &         c06, &
    &         c07, &
    &         c08, &
    &         c09, &
    &         c10, &
    &         c11, &
    &         c12, &
    &         c13, &
    &         c14, &
    &         c15, &
    &         c16, &
    &         c17, &
    &         c18, &
    &         c19, &
    &         c20, &
    &         c21, &
    &         c22, &
    &         c23, &
    &         c24, &
    &         c25, &
    &         c26, &
    &         c27
  REAL(RP) :: f(n,9)
  EXTERNAL :: fun


  c01 =  h / 9.0_RP
  c02 =  h * 0.4166666666666667e-01_RP
  c03 =  h * 0.125_RP
  c04 =  h / 6.0_RP
  c05 =  h * 0.5_RP
  c06 =  h * 0.6666666666666667_RP
  c07 =  h / 3.0_RP
  c08 =  h * 0.375_RP
  c09 =  h * 0.1333333333333333e+01_RP
  c10 =  h * 0.3333333333333333e+01_RP
  c11 =  h * 0.7e+01_RP
  c12 =  h * 0.9666666666666667e+01_RP
  c13 =  h * 0.1533333333333333e+02_RP
  c14 =  h * 0.6111111111111111_RP
  c15 =  h * 0.1166666666666667e+01_RP
  c16 =  h * 0.1375e+01_RP
  c17 =  h * 0.8333333333333333_RP
  c18 =  h * 0.4390243902439024_RP
  c19 =  h * 0.8780487804878049_RP
  c20 =  h * 0.1304878048780488e+01_RP
  c21 =  h * 0.2097560975609756e+01_RP
  c22 =  h * 0.2963414634146341e+01_RP
  c23 =  h * 0.4317073170731707e+01_RP
  c24 =  h * 0.3214285714285714e-01_RP
  c25 =  h * 0.4880952380952381e-01_RP
  c26 =  h * 0.2571428571428571_RP
  c27 =  h * 0.3238095238095238_RP

  f(:,:) =  0.0_RP
  f(:,1) =  y0(:)

!
!... 1-stage
!
  x1 =  x0
  CALL fun(x1, f(:,1), f(:,3), iout)
  f(:,2) =  c01 * f(:,3) &
    &    +        f(:,1)
!
!... 2-stage
!
  x1 =  x0 + c01
  CALL fun(x1, f(:,2), f(:,4), iout)
  f(:,2) =  c02 * f(:,3) &
    &    +  c03 * f(:,4) &
    &    +        f(:,1)
!
!... 3-stage
!
  x1 =  x0 + c04
  CALL fun(x1, f(:,2), f(:,5), iout)
  f(:,2) =  c04 * f(:,3) &
    &    -  c05 * f(:,4) &
    &    +  c06 * f(:,5) &
    &    +        f(:,1)
!
!... 4-stage
!
  x1 =  x0 + c07
  CALL fun(x1, f(:,2), f(:,6), iout)
  f(:,2) =  c03 * f(:,3) &
    &    +  c08 * f(:,6) &
    &    +        f(:,1)
!
!... 5-stage
!
  x1 =  x0 + c05
  CALL fun(x1, f(:,2), f(:,7), iout)
  f(:,2) = -c09 * f(:,3) &
    &    +  c10 * f(:,7) &
    &    -  c11 * f(:,4) &
    &    -  c12 * f(:,6) &
    &    +  c13 * f(:,5) &
    &    +        f(:,1)
!
!... 6-stage
!
  x1 =  x0 + c06
  CALL fun(x1, f(:,2), f(:,8), iout)
  f(:,2) = -c01 * f(:,3) &
    &    +  c03 * f(:,8) &
    &    +  c14 * f(:,7) &
    &    -  c15 * f(:,5) &
    &    +  c16 * f(:,4) &
    &    +        f(:,1)
!
!... 7-stage
!
  x1 =  x0 + c17
  CALL fun(x1, f(:,2), f(:,9), iout)
  f(:,2) = -c18 * f(:,8) &
    &    +  c19 * f(:,9) &
    &    +  c20 * f(:,3) &
    &    -  c21 * f(:,7) &
    &    -  c22 * f(:,4) &
    &    +  c23 * f(:,6) &
    &    +        f(:,1)
!
!... 8-stage
!
  x1 =  x0 + h
  CALL fun(x1, f(:,2), f(:,4), iout)
  f(:,1) =  c24 * (f(:,6) + f(:,8)) &
    &    +  c25 * (f(:,3) + f(:,4)) &
    &    +  c26 * (f(:,5) + f(:,9)) &
    &    +  c27 *  f(:,7)           &
    &    +         f(:,1)
  y0(:)  =  f(:,1)


  RETURN
END SUBROUTINE rkg6
