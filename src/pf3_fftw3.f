************************************
*  3D FFT with FFTW  
*
*  use FFTW library v.3.0.1
*
*  last modified:   2/3/2004 (CS)
*  created:        1/28/2004 (Chen Shen)
*
*  tested on: mse-uc02 (Redhat 9.0, gcc)             
************************************

!---------------------------------------
!-   Initialize FFT threaded plans.
!---------------------------------------
      subroutine FFT_Init_Threads(nthread)
      implicit none
      integer  nthread
      integer  iret
       call dfftw_init_threads(iret)
       print*,'iret=',iret
       call dfftw_plan_with_nthreads(nthread)
      end subroutine FFT_Init_Threads

c---------------------------------------
c-   Initialize FFT routines
c-   construct internal plan for c->c and c<-c
c---------------------------------------
      
      SUBROUTINE FFT_Init(option)
      implicit none

      INCLUDE 'pf3_fftw3.h'
      INTEGER option
    
c number of threads
      IF(pfFFTW3_initflag.eq.PF_FFTW3_INIT_TRUE) THEN
         print *, 'Warning: fft already initialized!'
      END IF

      IF(option.eq.1) THEN
         call dfftw_plan_dft_3d(pfFFTW3_plan_fwd, nx, ny, nz, 
     &        pfFFTW3_data, pfFFTW3_data,
     &        FFTW_FORWARD, FFTW_ESTIMATE)

         call dfftw_plan_dft_3d(pfFFTW3_plan_bkwd, nx, ny, nz, 
     &        pfFFTW3_data, pfFFTW3_data,
     &        FFTW_BACKWARD, FFTW_ESTIMATE)
      
      ELSE IF(option.eq.2) THEN
         call dfftw_plan_dft_3d(pfFFTW3_plan_fwd, nx, ny, nz,
     &        pfFFTW3_data, pfFFTW3_data,
     &        FFTW_FORWARD, FFTW_MEASURE)

         call dfftw_plan_dft_3d(pfFFTW3_plan_bkwd, nx, ny, nz,
     &        pfFFTW3_data, pfFFTW3_data,
     &        FFTW_BACKWARD, FFTW_MEASURE)

      ELSE IF(option.eq.3) THEN



         call dfftw_plan_dft_3d(pfFFTW3_plan_fwd, nx, ny, nz,
     &        pfFFTW3_data, pfFFTW3_data,
     &        FFTW_FORWARD, FFTW_PATIENT)



         call dfftw_plan_dft_3d(pfFFTW3_plan_bkwd, nx, ny, nz,
     &        pfFFTW3_data, pfFFTW3_data,
     &        FFTW_BACKWARD, FFTW_PATIENT)

      ELSE IF(option.eq.4) THEN
         call dfftw_plan_dft_3d(pfFFTW3_plan_fwd, nx, ny, nz,
     &        pfFFTW3_data, pfFFTW3_data,
     &        FFTW_FORWARD, FFTW_EXHAUSTIVE)

         call dfftw_plan_dft_3d(pfFFTW3_plan_bkwd, nx, ny, nz,
     &        pfFFTW3_data, pfFFTW3_data,
     &        FFTW_BACKWARD, FFTW_EXHAUSTIVE)
      ELSE
         STOP 'valid option (1-est,2-measure,3-patient,4-exhaust)!'
      END IF

      pfFFTW3_initflag =PF_FFTW3_INIT_TRUE

      RETURN
      END

c---------------------------------------
c-   Finalize FFT routines
c-   destroy internal plan for c->c and c<-c
c---------------------------------------
      
      SUBROUTINE FFT_Finish
      implicit none

      INCLUDE 'pf3_fftw3.h'

      IF(pfFFTW3_initflag.ne.PF_FFTW3_INIT_TRUE) THEN
         stop 'fft not initialized!'
      END IF

      call dfftw_destroy_plan(pfFFTW3_plan_fwd)
      call dfftw_destroy_plan(pfFFTW3_plan_bkwd)

      pfFFTW3_initflag =PF_FFTW3_INIT_FALSE

      RETURN
      END

      
c---------------------------------------
c-  forward 3D FFT , real to complex
c---------------------------------------

      subroutine fftrc3(in_data, out_data)
      implicit none

      INCLUDE 'pf3_fftw3.h'

c     subroutine arguments
      real  in_data(nx,ny,nz)
      complex out_data(nx,ny,nz)

c     local data
      integer i,j,k


      IF(pfFFTW3_initflag.ne.PF_FFTW3_INIT_TRUE) THEN
         stop 'fft not initialized!'
      END IF

c      convert input array to complex(8) type

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k)
!$OMP DO
      do k=1,nz
      do j=1,ny
      do i=1,nx
         pfFFTW3_data(i,j,k) =in_data(i,j,k)
      end do
      end do 
      end do
!$OMP END DO   
!$OMP END PARALLEL  

      call dfftw_execute(pfFFTW3_plan_fwd)

c     copy to OUT_DATA

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k)
!$OMP DO
      do k=1,nz
      do j=1,ny
      do i=1,nx
         out_data(i,j,k) =cmplx(pfFFTW3_data(i,j,k))
      end do
      end do
      end do
!$OMP END DO   
!$OMP END PARALLEL  

      return
      end


c---------------------------------------
c-   forward 3D FFT , complex to complex
c---------------------------------------

      subroutine fftcc3(in_data,out_data)
      implicit none

      INCLUDE 'pf3_fftw3.h'

*     subroutine arguments
      complex in_data(nx,ny,nz)
      complex out_data(nx,ny,nz)
*     local data
      integer i,j,k

      IF(pfFFTW3_initflag.ne.PF_FFTW3_INIT_TRUE) THEN
         stop 'fft not initialized!'
      END IF

*     convert input array to complex(8) type

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k)
!$OMP DO
      do k=1,nz
      do j=1,ny
      do i=1,nx
         pfFFTW3_data(i,j,k) =in_data(i,j,k)
      end do
      end do 
      end do
!$OMP END DO   
!$OMP END PARALLEL  
      
      call dfftw_execute(pfFFTW3_plan_fwd)

*     copy to OUT_DATA
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k)
!$OMP DO
      do k=1,nz
      do j=1,ny
      do i=1,nx
         out_data(i,j,k) =cmplx(pfFFTW3_data(i,j,k))
      end do
      end do
      end do
!$OMP END DO   
!$OMP END PARALLEL  
      return
      end

c---------------------------------------
c-   backward 3D FFT , complex to complex
c---------------------------------------

      subroutine ifftcc3(in_data,out_data)
      implicit none

      include 'pf3_fftw3.h'

*     subroutine arguments
      complex in_data(nx,ny,nz)
      complex out_data(nx,ny,nz)
*     local data
      integer i,j,k

      IF(pfFFTW3_initflag.ne.PF_FFTW3_INIT_TRUE) THEN
         stop 'fft not initialized!'
      END IF

*     convert input array to complex(8) type
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k)
!$OMP DO
      do k=1,nz
      do j=1,ny
      do i=1,nx
         pfFFTW3_data(i,j,k) =in_data(i,j,k)
      end do
      end do 
      end do
!$OMP END DO   
!$OMP END PARALLEL        
      call dfftw_execute(pfFFTW3_plan_bkwd)

*     copy to OUT_DATA
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k)
!$OMP DO
      do k=1,nz
      do j=1,ny
      do i=1,nx
         out_data(i,j,k) =cmplx(pfFFTW3_data(i,j,k))/nxyz
      end do
      end do
      end do
!$OMP END DO   
!$OMP END PARALLEL  
      return
      end

c---------------------------------------
c-   backward 3D FFT , complex to real
c---------------------------------------

      subroutine ifftcr3(in_data, out_data)
      implicit none
      include 'pf3_fftw3.h'

*     subroutine arguments
      complex in_data(nx,ny,nz)
      real  out_data(nx,ny,nz)
*     local data
      integer i,j,k


      IF(pfFFTW3_initflag.ne.PF_FFTW3_INIT_TRUE) THEN
         stop 'fft not initialized!'
      END IF

*     convert input array to complex(8) type
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k)
!$OMP DO

      do k=1,nz
      do j=1,ny
      do i=1,nx
         pfFFTW3_data(i,j,k) =in_data(i,j,k)
      end do
      end do 
      end do

!$OMP END DO   
!$OMP END PARALLEL 

      call dfftw_execute(pfFFTW3_plan_bkwd)

*     copy to OUT_DATA

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k)
!$OMP DO

      do k=1,nz
      do j=1,ny
      do i=1,nx
         out_data(i,j,k) =pfFFTW3_data(i,j,k)/nxyz
      end do
      end do
      end do

!$OMP END DO   
!$OMP END PARALLEL 

      return
      end


************************************
* end of file
************************************
