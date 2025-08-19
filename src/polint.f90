subroutine polint ( xa, ya, n, x, & ! (in)
  &                 y, dy         & ! (out)
  &               )

  implicit none

  integer, parameter :: ikind =  4, &
    &                   rkind =  8
!Arguments
  integer(kind=ikind), intent(in) :: n
  real(kind=rkind), intent(in) :: x,     &
    &                             xa(n), &
    &                             ya(n)
  real(kind=rkind), intent(out) :: y,  &
    &                              dy
!Local variables
  integer(kind=ikind) :: i,  &
    &                    m,  &
    &                    ns
  real(kind=rkind) :: den,  &
    &                 dif,  &
    &                 dift, &
    &                 ho,   &
    &                 hp,   &
    &                 w,    &
    &                 c(n), &
    &                 d(n)
  

  ns  =  1
  dif =  abs( x - xa(1) )

  loop100 : do i=1,n
    dift =  abs( x - xa(i) )
    if( dift < dif )then
      ns  =  i
      dif =  dift
    endif
    c(i) =  ya(i)
    d(i) =  ya(i)
  end do loop100

  y  =  ya(ns)
  ns =  ns - 1

  loop200 : do m=1,n-1
    loop210 : do i=1,n-m
      ho  =  xa(i)   - x
      hp  =  xa(i+m) - x
      w   =  c(i+1)  - d(i)
      den =  ho      - hp
      if( den == 0.0_rkind ) stop 'failure in polint'
      den  =  w  / den
      d(i) =  hp * den
      c(i) =  ho * den
    end do loop210
    if( 2 * ns < n - m )then
      dy =  c(ns+1)
    else
      dy =  d(ns)
      ns =  ns - 1
    end if
    y =  y + dy
  end do loop200


  return
end subroutine polint
