!==================================================================     
! Include file for using FFTW transforms in Fortran 90
!

      INTEGER, PARAMETER  ::   FFTW_FORWARD = -1, FFTW_BACKWARD = 1

      INTEGER, PARAMETER  ::   FFTW_REAL_TO_COMPLEX = -1,   &
                              FFTW_COMPLEX_TO_REAL = 1

      INTEGER, PARAMETER  ::   FFTW_ESTIMATE = 0, FFTW_MEASURE = 1

      INTEGER, PARAMETER  ::   FFTW_OUT_OF_PLACE = 0,  &
                              FFTW_IN_PLACE = 8,      &
                              FFTW_USE_WISDOM = 16

      INTEGER, PARAMETER  ::  FFTW_THREADSAFE = 128
!==================================================================     
