!=file_open
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
SUBROUTINE file_open

  USE param2,                ONLY : lmonitor,   &
    &                               lpunc
  USE cylindrical_coord_mod, ONLY : lflux
  USE file_name_mod,         ONLY : fopen,      &
    &                               mag_file,   &
    &                               flx_file,   &
    &                               gcr_file,   &
    &                               punc1_file, &
    &                               punc2_file

  IMPLICIT NONE


  READ(10,fopen)
  WRITE(6,fopen)

  OPEN(25, FILE=mag_file, FORM='unformatted')

  IF(lflux)THEN
    OPEN(26, FILE=flx_file, FORM='formatted')
  END IF

  IF(lmonitor)THEN
    OPEN(50, FILE=gcr_file, FORM='formatted')
  END IF

  IF(lpunc)THEN
    OPEN(61, FILE=punc1_file, FORM='formatted')
    OPEN(62, FILE=punc2_file, FORM='formatted')
  END IF


END SUBROUTINE file_open
