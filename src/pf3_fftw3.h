*******************************************
*  Header file for pf_fftw3.f
*
*  created:   1/28/2004 (Chen Shen)
*******************************************

      INCLUDE 'fftw3.f'
      include 'parameterc.inc'
      INTEGER  PF_FFTW3_INIT_TRUE, PF_FFTW3_INIT_FALSE
      PARAMETER(PF_FFTW3_INIT_TRUE =9000)
      PARAMETER(PF_FFTW3_INIT_FALSE=9001)

      DOUBLE COMPLEX  pfFFTW3_data(nx,ny,nz)
      INTEGER*8       pfFFTW3_plan_fwd, pfFFTW3_plan_bkwd
      INTEGER         pfFFTW3_initflag

      COMMON /PF_FFTW3/  pfFFTW3_data,
     &     pfFFTW3_plan_fwd, pfFFTW3_plan_bkwd,
     &     pfFFTW3_initflag
