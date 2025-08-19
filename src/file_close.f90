!=file_close.f90
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
SUBROUTINE file_close

  USE param2,                ONLY : lmonitor, &
    &                               lpunc
  USE cylindrical_coord_mod, ONLY : lflux

  IMPLICIT NONE


  CLOSE(25)

  IF(lflux)THEN
    CLOSE(26)
  END IF

  IF(lmonitor)THEN
    CLOSE(50)
  END IF

  IF(lpunc)THEN
    CLOSE(61)
    CLOSE(62)
  END IF


END SUBROUTINE file_close
