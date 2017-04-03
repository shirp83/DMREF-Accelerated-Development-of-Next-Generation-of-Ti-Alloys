      program dendrite
      implicit none
      include 'parameterc.inc'


      integer p,q, i, j, k, l, m, n, nout, amode, irun
      integer sa,sap
      integer ot1, OMP_GET_NUM_THREADS,ot2
      real tt, alpha, t_total
      real pi
      parameter (pi=3.1415926)
      real chcp1, chcp2, cbcc1, cbcc2
      real ca01,ca02,cb01,cb02
c  NZHOU ADD
      real randtest

      real tprr(nx,ny,nz)
      complex tprk(nx,ny,nz),eta2k(nx,ny,nz,ga,nvar)
      complex tpkk(nx,ny,nz)
      complex etak(nx,ny,nz,ga,nvar), dfek(nx,ny,nz,ga,nvar)

      real tprr1(nx,ny,nz)
      complex tprk1(nx,ny,nz)

      real temp

      real aaa,bbb
      
      real count1

      integer ix,iy,iz

      



      real sint, cost
      real rot(3,3), trot(3,3), kaba(3,3),kaba1(3,3)
      real inter_ener

      real max, value, sum

      count1=250.0*1.0

      open(10,file='dendrite_input',status='unknown')
      open(90,file='final.dat',status='unknown')

c   
c   dt          the integration time step,
c   t_total     total time
c   nout        number of output
c   amp0        noise amplification
c   amode       option to choose different way to describe anisotropy
c          1  gradient energy coefficient tensor and its rotation, 2-fold
c          2  general one, 2 or 4-fold, kapa=kapa0(1+a*cos(m*(theta-theta0)))
c             can be easily generalized to other fold symmetry
c          3  a special form of dependence of kapa on theta, equivalent 
c             option 1.
c  irun         initial configuration option
c  kapae       means kapa0
c  kapac       gradient energy coefficient for concentration, we set it to zero
c  mc          chemical mobility
c  agration   intensity of anisotropy
c  symt       if amode=2, symmetry must be specified. symt = 2 or 4.
c  theta0     rotation angle, if set to 0, no rotation, if other than zero,
c             the actual rotation angle = pi/theta0
c  c1, c2     initial concentration of phase1 and phase 2.


      read(10,*) dt, t_total, nout
      read(10,*) irun
      read(10,*) gamma, L0, lamda
      read(10,*) agratio
      read(10,*) chcp1, chcp2, cbcc1, cbcc2
      read(10,*) fc

      print *, dt, t_total, nout
      print *, irun
      print *, gamma, L0, lamda
      print *, agratio
      print *, chcp1, chcp2, cbcc1, cbcc2
      print *, fc 
      
      ca01=0.104821
      ca02=0.0238412
      cb01=0.0847926
      cb02=0.109545

      aaa=0.5

c amplitude of the randomnoise

c      fc=1.0e-4


      call FFT_Init_Threads(12)
      CALL fft_init(3)
      write(*,*) 'FFT initialized.'



      call interfacecalc()
      call get_ML(ca01,ca02,cb01,cb02)
      call sym_op()
      call polymisc()
      call grad()
      call elastic()
      call tieline_factor(ca01,ca02,cb01,cb02)

     

      if(theta0.ne.0.0) then
         theta0=pi/theta0
      endif

      tt = 0.0
      ot1 =0
      ot2 =0

**    initial configuration

      call initconfig(irun,chcp1, chcp2, cbcc1, cbcc2)
      write(*,*) "initialization done!"

      
     
      call outputfield(tt,irun)

      write(*,*) "outputfield done!"

**    first check if we need to calculate missing angle due to high
**    anisotropy

      write(*,*)"zero off virtral stress and strain"

***********************************************
***** start of anisotropy can be neglected in your program
**********************************************

**    time evolution

      do while (tt .le. t_total)
  


      if(tt.gt.(count1*dt*tao))then

            call inhom()

            else ! coherent

            call inhomc()

    
         endif

***********************************************
***** end of anisotropy can be neglected in your program
**********************************************

         call getdfce(ca01,ca02,cb01,cb02)
         
         write(*,*) "getfce..."
         
*        add the coherency strain energy
         do sa=1,ga
         do p=1,nvar

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k)
!$OMP DO

         do k=1,nz
         do j=1,ny
         do i=1,nx

            tprr(i,j,k)=eta(i,j,k,sa,p)

            tprr1(i,j,k)=dfe(i,j,k,sa,p)


         end do
         end do
         end do

!$OMP END DO   
!$OMP END PARALLEL  

         call fftrc3(tprr,tprk)
         call fftrc3(tprr1,tprk1)

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k)
!$OMP DO

         do k=1,nz
         do j=1,ny
         do i=1,nx

            etak(i,j,k,sa,p)=tprk(i,j,k)

            dfek(i,j,k,sa,p)=tprk1(i,j,k)


         end do
         end do
         end do

!$OMP END DO   
!$OMP END PARALLEL  



         end do ! end loop of p
         end do ! end loop of sa

! start evolution of eta field


         do sa=1,ga
         do p=1,nvar

         do i=1,3
            do j=1,3
               kaba1(i,j)=0.0
            end do
         end do

         kaba1(1,1)=1.0/9.0
         kaba1(2,2)=1.0
         kaba1(3,3)=4.0/9.0

         sint=sin(48/180.0*pi)
         cost=cos(48/180.0*pi)

        
         if (p.le.6) then
      
         do i=1,3
            do j=1,3
               rot(i,j)=0.0
               trot(i,j)=0.0
            end do
         end do

         rot(1,1)= cost
         rot(1,2)= sint
         rot(2,1)=-sint
         rot(2,2)= cost
         rot(3,3)= 1.0

         do i=1,3
         do j=1,3
            do k=1,3
            do l=1,3
               trot(i,j)=trot(i,j)+rot(i,k)*rot(j,l)*kaba1(k,l)
            end do
            end do       
         end do
         end do

         do i=1,3
         do j=1,3
            kaba1(i,j)=0
            do k=1,3
            do l=1,3
               kaba1(i,j)=kaba1(i,j)+s(k,i,p)*s(l,j,p)*trot(k,l)
            end do
            end do       
         end do
         end do

         else

         do i=1,3
            do j=1,3
               rot(i,j)=0.0
               trot(i,j)=0.0
            end do
         end do

         rot(1,1)= cost
         rot(1,2)=-sint
         rot(2,1)= sint
         rot(2,2)= cost
         rot(3,3)= 1.0

         do i=1,3
         do j=1,3
            do k=1,3
            do l=1,3
               trot(i,j)=trot(i,j)+rot(i,k)*rot(j,l)*kaba1(k,l)
            end do
            end do       
         end do
         end do

         do i=1,3
         do j=1,3
            kaba1(i,j)=0
            do k=1,3
            do l=1,3
               kaba1(i,j)=kaba1(i,j)+s(k,i,p-6)*s(l,j,p-6)*trot(k,l)
            end do
            end do       
         end do
         end do

         end if

         do i=1,3
         do j=1,3
            kaba(i,j)=0.0
            do k=1,3
            do l=1,3
               kaba(i,j)=kaba(i,j)+Qrot(k,i,sa)*Qrot(l,j,sa)*kaba1(k,l)
            end do
            end do       
         end do
         end do
         
        if((tt.le.(count1*dt*tao)).and.(fc.gt.0.0))then
c         if(1.gt.2)then
        call gaussian_series() 
        call fftrc3(randnum3, randnum3k)

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k,m,n,inter_ener)
!$OMP DO

         do k=1,nz
         do j=1,ny
         do i=1,nx

            inter_ener=0
            do m=1,3
            do n=1,3
               inter_ener=inter_ener+kaba(m,n)
     &              *g(m,i,j,k)*g(n,i,j,k)
            end do
            end do

            etak(i,j,k,sa,p)=(etak(i,j,k,sa,p)-dfek(i,j,k,sa,p)
     &       *me*dt+randnum3k(i,j,k)*dt-me*100.0
     &           *(appk(i,j,k,sa,p))*dt)/
     &           (1.0+kapae*inter_ener*me*dt)

         end do
         end do
         end do

!$OMP END DO   
!$OMP END PARALLEL

         else 

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k,m,n,inter_ener)
!$OMP DO

         do k=1,nz
         do j=1,ny
         do i=1,nx

            inter_ener=0
            do m=1,3
            do n=1,3
               inter_ener=inter_ener+kaba(m,n)
     &              *g(m,i,j,k)*g(n,i,j,k)
            end do
            end do

            etak(i,j,k,sa,p)=(etak(i,j,k,sa,p)-dfek(i,j,k,sa,p)
     &       *me*dt-me*100.0
     &           *(appk1(i,j,k,sa,p))*dt)/
     &(1.0+kapae*inter_ener*me*dt)




         end do
         end do
         end do

!$OMP END DO   
!$OMP END PARALLEL

         end if

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k)
!$OMP DO

         do k=1,nz
         do j=1,ny
         do i=1,nx

           tprk(i,j,k)=etak(i,j,k,sa,p)


         end do
         end do
         end do

!$OMP END DO   
!$OMP END PARALLEL  

         call ifftcr3(tprk,tprr)

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k)
!$OMP DO

         do k=1,nz
         do j=1,ny
         do i=1,nx

           eta(i,j,k,sa,p)=tprr(i,j,k)*shapeg(i,j,k,sa)


         end do
         end do
         end do

!$OMP END DO   
!$OMP END PARALLEL  



         
c        end if

         end do ! end loop of p
         end do ! end loop of sa

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k,p,sa,sum)
!$OMP DO


         do k=1,nz
         do j=1,ny
         do i=1,nx
            sum=0
            do p=1,nvar
            do sa=1,ga
               if (eta(i,j,k,sa,p).lt. epsilon) then
                  eta(i,j,k,sa,p)=0
               else if(eta(i,j,k,sa,p).gt. 1.0-epsilon) then
                  eta(i,j,k,sa,p)=1.0
               end if
               sum=sum+eta(i,j,k,sa,p)
            end do
            end do
            if (sum.gt.1) then
               do p=1,nvar
               do sa=1,ga
                  eta(i,j,k,sa,p)=eta(i,j,k,sa,p)/sum
               end do
               end do
            end if
         end do
         end do
         end do


!$OMP END DO   
!$OMP END PARALLEL


         call Euler()

         tt = tt + dt*tao
         ot1 =ot1+1

*        output intermediate composition field, nout of outputs
 66      format('t = ', I4)
         if(ot1 .ge. t_total/(dt*Tao)/nout*ot2)  then

            write(*,*) tt

            call outputfield(real(ot2),irun)
            ot2 = ot2+1
         endif

      end do

 166  format('Job finished!')
      write(6,166)



      end



      subroutine outputfield(time,irun)
      include 'parameterc.inc'
      
      real    time,sum
      integer i,j,k,irun,p
      integer sa,sap
      character filnam*20, ctime*20, i2c*80

      itime=time
      ctime=i2c(itime)

      filnam='restart_'
      filnam(lens(filnam)+1:)=ctime(1:lens(ctime))
      filnam(lens(filnam)+1:)='.dat'

      open(unit=20,file=filnam(1:lens(filnam)),status='unknown')

      do sa=1,ga ! output the distribution of each variant in each grain
      
      do k=1,nz
      do j=1,ny
      do i=1,nx
         id1=0
         id2=0
         value1=0
         value2=0

         do p=1,nvar
            if (eta(i,j,k,sa,p) .gt. value1) then
               value2=value1
               id2=id1
               id1=p
               value1=eta(i,j,k,sa,p)
            else if(eta(i,j,k,sa,p) .gt. value2) then
               value2=eta(i,j,k,sa,p)
               id2=p
            end if
         end do
         write(20,*) id1, value1, id2, value2
      end do
      end do
      end do

      end do ! end loop of sa
      
      close(20)


      filnam='field_'
      filnam(lens(filnam)+1:)=ctime(1:lens(ctime))
      filnam(lens(filnam)+1:)='.dat'
      
      open(unit=20,file=filnam(1:lens(filnam)),status='unknown')
      
      do k=1,nz
      do j=1,ny
      do i=1,nx
         sum=0
         do p=1,nvar
         do sa=1,ga
            sum=sum+eta(i,j,k,sa,p)*eta(i,j,k,sa,p)
         end do
         end do
         write(20,21) sum
      end do
      end do
      end do

      close(20)     

      filnam='c1_'
      filnam(lens(filnam)+1:)=ctime(1:lens(ctime))
      filnam(lens(filnam)+1:)='.dat'
      
      open(unit=30,file=filnam(1:lens(filnam)),status='unknown')
           
      do i=1,nx
         do j=1,ny
            write(30,21) (c1(i,j,k),k=1,nz)
         end do
      end do

      close(30)

      filnam='c2_'
      filnam(lens(filnam)+1:)=ctime(1:lens(ctime))
      filnam(lens(filnam)+1:)='.dat'
      
      open(unit=40,file=filnam(1:lens(filnam)),status='unknown')
           
      do i=1,nx
         do j=1,ny
            write(40,21) (c2(i,j,k),k=1,nz)
         end do
      end do

      close(40)
     
 21   format(1x,1024f10.4)

      return
      end


      


c...  to setup the initial configuration :
c...  the initial value of order parameter and composition

      subroutine initconfig(mode,chcp1, chcp2, cbcc1, cbcc2)

      include 'parameterc.inc'

      integer i,j,k,ir2,mode,ix,p
      integer sa,sap
      real chcp1, chcp2, cbcc1, cbcc2
      real a, b, sint, cost, x1, y1,z1
      real pi
      parameter(pi=3.14159265758)

      real rot(3,3), trot(3,3), rot2(3,3)
      real center(3,nvar)    

      sint=sin(52.0/180.0*pi)
      cost=cos(52.0/180.0*pi)




      do i=1,3
      do j=1,3
         rot(i,j)=0.0
         rot2(i,j)=0.0
      end do
      end do

      rot(1,1)= cost
      rot(1,2)=-sint
      rot(2,1)= sint
      rot(2,2)= cost
      rot(3,3)= 1.0

      rot2(1,1)= cost
      rot2(1,2)= sint
      rot2(2,1)=-sint
      rot2(2,2)= cost
      rot2(3,3)= 1.0
         
      if(mode.eq.0) then

         do k=1,nz
         do j=1,ny        
         do i=1,nx
            do p=1,nvar
            do sa=1,ga
               eta(i,j,k,sa,p)=0
            end do
            end do
         end do
         end do
         end do 

         open(unit=30,file='restart_input',status='unknown')

         do sa=1,ga

         do k=1,nz
         do j=1,ny        
         do i=1,nx
            read (30,*) id1, value1, id2, value2
            if (id1.ne.0) eta(i,j,k,sa,id1)=value1
            if (id2.ne.0) eta(i,j,k,sa,id2)=value2
         end do
         end do
         end do
         
         end do!end loop of sa

         close(30)

         open(unit=30,file='c1_input',status='unknown')     
         do i=1,nx   
         do j=1,ny         
            read (30,*) (c1(i,j,k),k=1,nz)
         end do
         end do  
         close(30)

         open(unit=30,file='c2_input',status='unknown')           
         do i=1,nx 
         do j=1,ny     
            read (30,*) (c2(i,j,k),k=1,nz)
         end do
         end do
         close(30)

      else if(mode.eq.1) then
c  NZHOU MODIFIED         
cNZ         i=sdrandom(1000)

         call init_random_seed()
         

         


         do i=1, nx
         do j=1, ny
         do k=1, nz
            do p=1,nvar
            do sa=1,ga
                eta(i,j,k,sa,p)=0.0
            end do
            end do
            c1(i,j,k)  = cbcc1
            c2(i,j,k)  = cbcc2
         end do
         end do
         end do         


      else if(mode.eq.2) then
         
         a=3.0

         do i=1, nx
         do j=1, ny
         do k=1, nz
            do p=1,nvar
            do sa=1,ga

                eta(i,j,k,sa,p)=0.0

            end do
            end do
            c1(i,j,k)  = cbcc1
            c2(i,j,k)  = cbcc2
         end do
         end do
         end do 


        do p=1,1!1,nvar
            do k=1, nz
            do j=1, ny
            do i=1, nx

               x1= i-nx/2
               y1= j-ny/2
               z1= k-nz/2

               if(sqrt(y1**2+z1**2+x1**2).lt.a) then 
                  eta(i,j,k,1,p) = 1.0
                  c1(i,j,k)  = chcp1
                  c2(i,j,k)  = chcp2
               endif
            end do
            end do
            end do
         end do


    
      end if

      return
      end

      integer function lens(line)
      character line*(*)
      ll=len(line)
      do 10 i=ll,1,-1
         i1=ichar(line(i:i))
         if(i1.gt.32) goto 20
 10   continue
      i=0
 20   lens=i
 900  return
      end

      character*80 function i2c(i)
      implicit none
      integer ia(80)
      integer i,j,k,it,nc

      i2c=''
      it=i
      do k=1,80
         ia(k)=mod(it,10)
         if(it.ge.10) then
            it=it/10
         else
            nc=k
            goto 10
         endif
      enddo
 10   continue
      do k=1,nc
         i2c(nc-k+1:nc-k+1)=char(ia(k)+48)
      enddo
      return
      end

 
c************************************************
c
c  initiate the tieline  calculation
c
c************************************************

      subroutine tieline_factor(ca01,ca02,cb01,cb02)
      implicit real (a-h,o-z)

      include "parameterc.inc"
      real ca01, ca02, cb01, cb02

      GHSERAL=-11277.683+188.661987*T-31.748192*T*LOG(T)
c     &     -1.234E+28*T**(-9)
     &     -1.234E+28*T**(-9)+137529.83
      GHSERTI=+908.837+67.048538*T-14.9466*T*LOG(T)-.0081465*T**2
c     &     +2.02715E-07*T**3-1477660*T**(-1)
     &     +2.02715E-07*T**3-1477660*T**(-1)+58166.677
      GHSERV=-7967.84+143.291*T-25.9*T*LOG(T)
c     &+6.25E-05*T**2-6.8E-07*T**3
     &+6.25E-05*T**2-6.8E-07*T**3+77277.217

      GBCCAL=10083-4.812*T+GHSERAL
      GBCCTI=+5758.6+38.38987*T-7.4305*T*LOG(T)+.0093636*T**2
     &     -1.04805E-06*T**3-525093*T**(-1)+GHSERTI

      GHCPAL=+5481-1.799*T+GHSERAL
      GHCPV=+4000+2.4*T+GHSERV


C     3=TI 1=AL 2=V

      BL13_0 = -128500+39*T
      BL13_1 = 6000
      BL13_2 = 21200

      BL132_0 = -8000

      BL12_0 = -98000+32*T
      BL12_1 = +3500-5*T

      BL32_0 = +10500.0-1.5*T
      BL32_1 = 2000.0
      BL32_2 = 1000.0

      HL13_0 = -133750+39*T
      HL13_1 = 250
      HL13_2 = 17250

      HL132_0 = -60000
      HL132_1 = -30000
      HL132_2 = 60000

      HL12_0 = -98000+32*T
      HL12_1 = +3500-5*T

      HL32_0 = +30250-10.0*T

c****** equilibrium tieline composition for hcp
    
      X1=ca01        
      X2=ca02      
      X3=1.0-X1-X2
      
      dfcc_hcp1=R*T*(1.0/X1+1.0/X3)+2.0*X2*HL12_1
     &     -2.0*(HL13_0+HL13_1*(x1-x3)+HL13_2*(X1-X3)*(X1-X3))
     &     +2.0*(X3-X1)*(2.0*HL13_1+4.0*HL13_2*(X1-X3))
     &     +X1*X3*8.0*HL13_2
     &     -X2*2.0*(X1*HL132_0+X3*HL132_1+X2*HL132_2)
     &     +X2*(X3-X1)*(HL132_0-HL132_1)
     &     +(X3-X1)*X2*(HL132_0-HL132_1)


      dfcc_hcp2=R*T*(1.0/X2+1.0/X3)
     &     -2.0*X1*HL12_1
     &     -2.0*HL32_0
     &     -2.0*X1*(HL13_1+2.0*HL13_2*(x1-x3))
     &     +X1*X3*2.0*HL13_2
     &     -X1*2.0*(X1*HL132_0+X3*HL132_1+X2*HL132_2)
     &     +2.0*X1*(X3-X2)*(HL132_2-HL132_1) 

     
      
c****** equilibrium tieline composition for bcc     
 
      X1=cb01
      X2=cb02
      X3=1.0-X1-X2      

      dfcc_bcc1=R*T*(1.0/x1+1.0/x3)+2.0*X2*BL12_1
     &     +2.0*X2*(BL32_1+2.0*BL32_2*(X3-X2))
     &     +X3*X2*2.0*BL32_2
     &     -2.0*(BL13_0+BL13_1*(x1-x3)+BL13_2*(X1-X3)*(X1-X3))
     &     +2.0*(X3-X1)*(2.0*BL13_1+4.0*BL13_2*(X1-X3))
     &     +X1*X3*8.0*BL13_2
     &     -2.0*BL132_0*X2


      dfcc_bcc2=R*T*(1.0/X2+1.0/X3)-2.0*X1*BL12_1
     &     -2.0*(BL32_0+BL32_1*(x3-x2)+BL32_2*(X3-X2)*(X3-X2))
     &     +2.0*(X3-X2)*(-2.0*BL32_1-4.0*BL32_2*(X3-X2))
     &     +X3*X2*8.0*BL32_2
     &     -2.0*X1*(BL13_1+2.0*BL13_2*(x1-x3))
     &     +X1*X3*2.0*BL13_2
     &     -2.0*BL132_0*X1


      print *, dfcc_bcc1,dfcc_bcc2,dfcc_hcp1,dfcc_hcp2

      return
      end


c.... this subroutine calculate the chemical potential and structure affinity.
c.... the model system is the cu-ni alloy consisting of liquid and fcc phases.
c.... non-equilibrium free energy constructed by following Wang et al. 1993.
c....
      subroutine getdfce(ca01,ca02,cb01,cb02)
      implicit real (a-h,o-z)

      include 'parameterc.inc'

      integer p,q,sa,sap
      real sum, df_deta(ga,nvar),dfe_bcc
      real ca01,ca02,cb01,cb02
      real sign


      GHSERAL=-11277.683+188.661987*T-31.748192*T*LOG(T)
     &     -1.234E+28*T**(-9)+137529.83
      GHSERTI=+908.837+67.048538*T-14.9466*T*LOG(T)-.0081465*T**2
     &     +2.02715E-07*T**3-1477660*T**(-1)+58166.677
      GHSERV=-7967.84+143.291*T-25.9*T*LOG(T)+6.25E-05*T**2-6.8E-07*T**3
     &     +77277.217

      GBCCAL=10083-4.812*T+GHSERAL
      GBCCTI=+5758.6+38.38987*T-7.4305*T*LOG(T)+.0093636*T**2
     &     -1.04805E-06*T**3-525093*T**(-1)+GHSERTI

      GHCPAL=+5481-1.799*T+GHSERAL
      GHCPV=+4000+2.4*T+GHSERV


C     3=TI 1=AL 2=V

      BL13_0 = -128500+39*T
      BL13_1 = 6000
      BL13_2 = 21200

      BL132_0 = -8000

      BL12_0 = -98000+32*T
      BL12_1 = +3500-5*T

      BL32_0 = +10500.0-1.5*T
      BL32_1 = 2000.0
      BL32_2 = 1000.0

      HL13_0 = -133750+39*T
      HL13_1 = 250
      HL13_2 = 17250

      HL132_0 = -60000
      HL132_1 = -30000
      HL132_2 = 60000

      HL12_0 = -98000+32*T
      HL12_1 = +3500-5*T

      HL32_0 = +30250-10.0*T


!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k,cx1,cx2,sa,sap,p,q,ett2,ett1,pfeta,X1,X2,X3,
!$OMP&Xbcc1,Xbcc2,GBCC,dGBCC_1,dGBCC_2,Xhcp1,Xhcp2,
!$OMP&GHCP,dGHCP_1,dGHCP_2,pfeta1,pfeta2,dpfeta2,dqfeta2,
!$OMP&dGall_1,dGall_2,sum,df_deta,dfe_bcc,sign)
!$OMP DO      

      do k=1,nz
      do j=1,ny
      do i=1,nx

         cx1=c1(i,j,k)
         cx2=c2(i,j,k)
         ett2=0

         do p=1,nvar
         do sa=1,ga
            ett2=ett2+eta(i,j,k,sa,p)
         end do    
         end do

         ett1=1.0-ett2      

         pfeta=ett2*ett2*ett2*(10.0-15.0*ett2+6.0*ett2*ett2)

c********* composition of bcc phase  phase on right with 1-p

         X1=(dfcc_hcp1*cx1-pfeta*(ca01*dfcc_hcp1
     &        -cb01*dfcc_bcc1))
     &        /(pfeta*dfcc_bcc1+(1-pfeta)*dfcc_hcp1)
         
         
         X2=(dfcc_hcp2*cx2-pfeta*(ca02*dfcc_hcp2
     &        -cb02*dfcc_bcc2))
     &        /(pfeta*dfcc_bcc2+(1-pfeta)*dfcc_hcp2)
         
         X3=1.0-X1-X2
         
         Xbcc1=X1
         Xbcc2=X2
         
         GBCC=x3*GBCCTI+x1*GBCCAL+x2*GHSERV+R*T*(x3
     $        *log(x3)+x1*log(x1)+x2*log(x2))+X1*X2
     &        *(BL12_0+BL12_1*(X1-X2))+x3*x2
     $        *(BL32_0+BL32_1*(x3-x2)+BL32_2*(x3-x2)*(x3-x2))+x1*x3
     $        *(BL13_0+BL13_1*(x1-x3)+BL13_2*(X1-X3)*(X1-X3))
     &        +x1*x2*x3*BL132_0
         
         
         dGBCC_1=-GBCCTI+GBCCAL+R*T*LOG(X1/X3)+X2
     &        *(BL12_0+BL12_1*(X1-X2))+X1*X2*BL12_1
     &        -X2*(BL32_0+BL32_1*(x3-x2)+BL32_2*(x3-x2)*(x3-x2))
     &        -X3*X2*(BL32_1+2.0*BL32_2*(X3-X2))
     &        +(X3-X1)*(BL13_0+BL13_1*(x1-x3)+BL13_2*(X1-X3)*(X1-X3))
     &        +X1*X3*(2.0*BL13_1+4*BL13_2*(X1-X3))
     &        +BL132_0*X2*(X3-X1)
         
         dGBCC_2=-GBCCTI+GHSERV+R*T*LOG(X2/X3)+X1
     &        *(BL12_0+BL12_1*(X1-X2))-X1*X2*BL12_1
     &        +(X3-X2)*(BL32_0+BL32_1*(x3-x2)+BL32_2*(X3-X2)*(X3-X2))
     &        +X3*X2*(-2.0*BL32_1-4*BL32_2*(X3-X2))
     &        -X1*(BL13_0+BL13_1*(x1-x3)+BL13_2*(x1-x3)*(x1-x3))
     &        +X1*X3*(BL13_1+2.0*BL13_2*(X1-X3))
     &        +BL132_0*X1*(X3-X2)
         
c*********composition in Hcp phase on left with p
         
         
         X1=(dfcc_bcc1*cx1-(1.0-pfeta)*(cb01*dfcc_bcc1
     &        -ca01*dfcc_hcp1))/(pfeta*dfcc_bcc1+(1.0-pfeta)
     &        *dfcc_hcp1)
         
         X2=(dfcc_bcc2*cx2-(1.0-pfeta)*(cb02*dfcc_bcc2
     &         -ca02*dfcc_hcp2))/(pfeta*dfcc_bcc2+(1.0-pfeta)
     &        *dfcc_hcp2)
         
         X3=1-X1-X2
         
         Xhcp1=X1
         Xhcp2=X2
         
         GHCP=x3*GHSERTI+x1*GHCPAL+x2*GHCPV+R*T*(x3
     $        *log(x3)+x1*log(x1)+x2*log(x2))+X1*X2
     &        *(HL12_0+HL12_1*(X1-X2))+x3*x2*HL32_0
     $        +x1*x3*(HL13_0+HL13_1*(x1-x3)+HL13_2*(X1-X3)*(X1-X3))
     &        +x1*x2*x3*(X1*HL132_0+X3*HL132_1+X2*HL132_2)
         
         
         dGHCP_1=-GHSERTI+GHCPAL+R*T*LOG(X1/X3)+X2
     &        *(HL12_0+HL12_1*(X1-X2))+X1*X2*HL12_1
     &        -X2*HL32_0
     &        +(X3-X1)*(HL13_0+HL13_1*(x1-x3)+HL13_2*(X1-X3)*(X1-X3))
     &        +X1*X3*(2.0*HL13_1+4*HL13_2*(X1-X3))
     &        +X2*(X3-X1)*(X1*HL132_0+X3*HL132_1+X2*HL132_2)
     &        +X1*X2*X3*(HL132_0-HL132_1)
         
         dGHCP_2=-GHSERTI+GHCPV+R*T*LOG(X2/X3)+X1
     &        *(HL12_0+HL12_1*(X1-X2))-X1*X2*HL12_1
     &        +(X3-X2)*HL32_0
     &        -X1*(HL13_0+HL13_1*(x1-x3)+HL13_2*(x1-x3)*(x1-x3))
     &        +X1*X3*(HL13_1+2.0*HL13_2*(X1-X3))
     &        +X1*(X3-X2)*(X1*HL132_0+X3*HL132_1+X2*HL132_2)
     &        +X1*X2*X3*(HL132_2-HL132_1)


         pfeta1=ett1*ett1*ett1*(10.0-15.0*ett1+6.0*ett1*ett1)
         pfeta2=ett2*ett2*ett2*(10.0-15.0*ett2+6.0*ett2*ett2)

         dGall_1=pfeta1*dGBCC_1+pfeta2*dGHCP_1
         dGall_2=pfeta1*dGBCC_2+pfeta2*dGHCP_2

         dfc1(i,j,k)=dGall_1/gnormal            
         dfc2(i,j,k)=dGall_2/gnormal
        
         dpfeta2=30.0*ett2*(1.0-ett2)*ett2*(1.0-ett2)

! now we are considering the polycrystal case
                 
         dfe_bcc=(GBCC-Xbcc1*dGHCP_1-Xbcc2*dGHCP_2)*dpfeta2
     &         +wmega*ett2

         
         do p=1,nvar
         do sa=1,ga
            df_deta(sa,p)=(GHCP-Xhcp1*dGHCP_1-Xhcp2*dGHCP_2)*dpfeta2
     &         +wmega*(1.0-eta(i,j,k,sa,p))
         end do
         end do

         do p=1,nvar
         do sa=1,ga
            dfe(i,j,k,sa,p)=0.0
            if ( eta(i,j,k,sa,p) .gt. epsilon) then
               
               if(ett1.gt.epsilon) 
     &              dfe(i,j,k,sa,p)=df_deta(sa,p)-dfe_bcc
               do q=1,nvar
               do sap=1,ga
                  if (eta(i,j,k,sap,q) .gt.epsilon) then
                     dfe(i,j,k,sa,p)=dfe(i,j,k,sa,p)+df_deta(sa,p)
     $                               -df_deta(sap,q)
                  end if
               end do ! end loop of sap
               end do ! end loop of q
            end if
            dfe(i,j,k,sa,p)=dfe(i,j,k,sa,p)/gnormal!/ga
         end do
         end do

      end do
      end do
      end do

!$OMP END DO   
!$OMP END PARALLEL


     

      return
      end
   

      subroutine euler ()
      implicit real (a-h,o-z)
      include 'parameterc.inc'
      integer left, right, down, up, top, bottom
      integer p
      integer sa,sap

      real M11(nx,ny,nz),M12(nx,ny,nz)
      real M21(nx,ny,nz),M22(nx,ny,nz)

*        euler method for concentration 
      QALTI1=204000
      QALTI2=96000
      FALTI1=5.51E-6
      FALTI2=5.19E-10

      QALAL1=19000
      FALAL1=5E-4

      QALV1=290000
      FALV1=1.85E-5

      QTITI1=121000
      QTITI2=237000
      FTITI1=1.47E-8
      FTITI2=5.91E-5

      QVTI1=107000
      QVTI2=212000
      FVTI1=1.55E-9
      FVTI2=2.54E-5

      QTIV1=320000
      QTIV2=544000
      FTIV1=1.62E-4
      FTIV2=51

      QVV1=306000
      QVV2=493000
      FVV1=2.73E-5
      FVV2=1.58
         
      ALTITIV0=111000-20*T
      ALTITIV1=-43000+10*T
      ALVTIV0=78000+11*T
      ALVTIV1=47000-22*T

      ALTIALTI0=-247000-83*T
      ALTIALTI1=-146000
      ALALALTI0=-1059000+373*T
      ALALALTI1=-921000+342*T

      ALALALV0=-4143000
      ALVALV0=-5752000

      ALALTIV0=+180000-277*T

      
      VVAL12=-0.12027*QALTI1*T**(-1)
      VVAL14=-0.12027*QALTI2*T**(-1)
      DTALTI=FALTI1*EXP(VVAL12)+FALTI2*EXP(VVAL14)

      DTALAL=-QALAL1+R*T*LOG(FALAL1)

      DTALV =-QALV1+R*T*LOG(FALV1)

      VV12=-0.12027*QTITI1*T**(-1)
      VV14=-0.12027*QTITI2*T**(-1)
      DTTITI=FTITI1*EXP(VV12)+FTITI2*EXP(VV14)
      VV22=-0.12027*QTIV1*T**(-1)
      VV24=-0.12027*QTIV2*T**(-1)
      DTTIV=FTIV1*EXP(VV22)+FTIV2*EXP(VV24)
      VV32=-0.12027*QVTI1*T**(-1)
      VV34=-0.12027*QVTI2*T**(-1)
      DTVTI=FVTI1*EXP(VV32)+FVTI2*EXP(VV34)
      VV42=-0.12027*QVV1*T**(-1)
      VV44=-0.12027*QVV2*T**(-1)
      DTVV=FVV1*EXP(VV42)+FVV2*EXP(VV44)

      AMTITIHCP=-3.03E5+R*T*LOG(1.35E-3)
      AMALTIHCP=-3.29E5+R*T*LOG(6.6E-3)
      AMVTIHCP=-2.50E5+R*T*LOG(1E-3)

c 1=al 2=v 3=ti
c+++++++++++++++++++++++++ mobility of Al (B2)+++++++++++++++++++++++

      G1_1=DTALAL
      G1_2=DTALV
      G1_3=R*T*LOG(DTALTI)
      
      G1_13_0=ALALALTI0
      G1_13_1=ALALALTI1
      G1_32_0=ALALTIV0
      G1_12_0=ALALALV0

c+++++++++++++++++++++++++ mobility of V (B2)+++++++++++++++++++++++

      G2_1=DTALAL
      G2_2=R*T*LOG(DTVV)
      G2_3=R*T*LOG(DTVTI)
         
      G2_12_0=ALVALV0
      G2_32_0=ALVTIV0
      G2_32_1=ALVTIV1

c+++++++++++++++++++++++++ mobility of TI (B3)+++++++++++++++++++++++

      G3_1=DTALAL
      G3_2=R*T*LOG(DTTIV)
      G3_3=R*T*LOG(DTTITI)

      G3_13_0=ALTIALTI0
      G3_13_1=ALTIALTI1
      G3_32_0=ALTITIV0
      G3_32_1=ALTITIV1


c+++++++++++++++++++++++++ mobility of Al (B2)+++++++++++++++++++++++

      HG1_1=AMALTIHCP
      HG1_2=AMALTIHCP
      HG1_3=AMALTIHCP
      
c+++++++++++++++++++++++++ mobility of V (B2)+++++++++++++++++++++++

      HG2_1=AMVTIHCP
      HG2_2=AMVTIHCP
      HG2_3=AMVTIHCP
         

c+++++++++++++++++++++++++ mobility of TI (B3)+++++++++++++++++++++++

      HG3_1=AMTITIHCP
      HG3_2=AMTITIHCP
      HG3_3=AMTITIHCP

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k,p,sa,x1,x2,x3,ett,DG1_B,B1_B,DG2_B,B2_B,DG3_B,B3_B,
!$OMP&DG1_H,B1_H,DG2_H,B2_H,DG3_H,B3_H, B1, B2, B3, AA, EE1)
!$OMP DO

      do k=1,nz
      do j=1,ny
      do i=1,nx

         x1=c1(i,j,k)
         x2=c2(i,j,k)
         x3=1-x1-x2
         ett=0
         do p=1,nvar
         do sa=1,ga 
            ett=ett+eta(i,j,k,sa,p)
         end do
         end do

         if (ett .gt. 0.95) ett=1.0


c------------------------- mobility of AL (B1) ---------------------------
            
         DG1_B=x1*G1_1+x2*G1_2+x3*G1_3+x1*x2*G1_12_0+x3*x2*G1_32_0
     &        +x1*x3*(G1_13_0+G1_13_1*(X1-X3))

         B1_B=exp(DG1_B/R/T)/R/T/Bnormal
c------------------------- mobility of V (B2) ---------------------------

         DG2_B=x1*G2_1+x2*G2_2+x3*G2_3+x1*x2*G2_12_0+x3*x2*(G2_32_0
     &        +G2_32_1*(X3-X2))
         
         B2_B=exp(DG2_B/R/T)/R/T/Bnormal
c
c------------------------- mobility of TI (B3) ---------------------------

         DG3_B=x1*G3_1+x2*G3_2+x3*G3_3+x3*x2*(G3_32_0+
     &        G3_32_1*(X3-X2))+x1*x3*(G3_13_0+G3_13_1*(X1-X3))

         B3_B=exp(DG3_B/R/T)/R/T/Bnormal

C-----------------          HCP           ---------------------------

         DG1_H=HG1_1

         B1_H=exp(DG1_H/R/T)/R/T/Bnormal
            
c------------------------- mobility of V (B2) ---------------------------

         DG2_H=HG2_1

         B2_H=exp(DG2_H/R/T)/R/T/Bnormal

c------------------------- mobility of TI (B3) ---------------------------

         DG3_H=HG3_1

         B3_H=exp(DG3_H/R/T)/R/T/Bnormal

         B1=EXP((DG1_H+ETT*(DG1_B-DG1_H))/R/T)/R/T/Bnormal
         B2=EXP((DG2_H+ETT*(DG2_B-DG2_H))/R/T)/R/T/Bnormal
         B3=EXP((DG3_H+ETT*(DG3_B-DG3_H))/R/T)/R/T/Bnormal

         B1=-B1+B1_H+B1_B
         B2=-B2+B2_H+B2_B
         B3=-B3+B3_H+B3_B

         AA=-(1-x1)*x1*B1*(1-B3/B1)+x1*x2*B2*(1-B3/B2)
         EE1=-(1-x2)*x2*B2*(1-B3/B2)+x1*x2*B1*(1-B3/B1)
         M11(i,j,k)=((1-x1)*x1*B1+AA*x1)
         M12(i,j,k)=(AA*x2-x1*x2*B2)
         M21(i,j,k)=(EE1*x1-x1*x2*B1)
         M22(i,j,k)=((1-x2)*x2*B2+EE1*x2)


      end do      
      end do
      end do

!$OMP END DO   
!$OMP END PARALLEL


!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&PRIVATE(i,j,k,left,right,down,up,bottom,top)
!$OMP DO

      do k=1,nz
      do j=1,ny
      do i=1,nx

c================ periodical boundary conditions ====================

         left=i-1
         right=i+1
         down=j-1
         up=j+1
         bottom=k-1
         top=k+1
         if(left.eq.0) left=nx
         if(down.eq.0) down=ny
         if(bottom.eq.0) bottom=nz
         if(right.eq.nx+1) right=1
         if(up.eq.ny+1) up=1
         if(top.eq.nz+1) top=1


         c1(i,j,k)=c1(i,j,k)
     &        +dt*((m11(right,j,k)+m11(i,j,k))/2.0
     &        *(dfc1(right,j,k)-dfc1(i,j,k))/dx
     &        -(m11(i,j,k)+m11(left,j,k))/2.0
     &        *(dfc1(i,j,k)-dfc1(left,j,k))/dx)/dx
     &        +dt*((m11(i,down,k)+m11(i,j,k))/2.0
     &        *(dfc1(i,down,k)-dfc1(i,j,k))/dx
     &        -(m11(i,j,k)+m11(i,up,k))/2.0
     &        *(dfc1(i,j,k)-dfc1(i,up,k))/dx)/dx
     &        +dt*((m11(i,j,top)+m11(i,j,k))/2.0
     &        *(dfc1(i,j,top)-dfc1(i,j,k))/dx
     &        -(m11(i,j,k)+m11(i,j,bottom))/2.0
     &        *(dfc1(i,j,k)-dfc1(i,j,bottom))/dx)/dx
     &        +dt*((m12(right,j,k)+m12(i,j,k))/2.0
     &        *(dfc2(right,j,k)-dfc2(i,j,k))/dx
     &        -(m12(i,j,k)+m12(left,j,k))/2.0
     &        *(dfc2(i,j,k)-dfc2(left,j,k))/dx)/dx
     &        +dt*((m12(i,down,k)+m12(i,j,k))/2.0
     &        *(dfc2(i,down,k)-dfc2(i,j,k))/dx
     &        -(m12(i,j,k)+m12(i,up,k))/2.0
     &        *(dfc2(i,j,k)-dfc2(i,up,k))/dx)/dx
     &        +dt*((m12(i,j,top)+m12(i,j,k))/2.0
     &        *(dfc2(i,j,top)-dfc2(i,j,k))/dx
     &        -(m12(i,j,k)+m12(i,j,bottom))/2.0
     &        *(dfc2(i,j,k)-dfc2(i,j,bottom))/dx)/dx


         c2(i,j,k)=c2(i,j,k)
     &        +dt*((m21(right,j,k)+m21(i,j,k))/2.0
     &        *(dfc1(right,j,k)-dfc1(i,j,k))/dx
     &        -(m21(i,j,k)+m21(left,j,k))/2.0
     &        *(dfc1(i,j,k)-dfc1(left,j,k))/dx)/dx
     &        +dt*((m21(i,down,k)+m21(i,j,k))/2.0
     &        *(dfc1(i,down,k)-dfc1(i,j,k))/dx
     &        -(m21(i,j,k)+m21(i,up,k))/2.0
     &        *(dfc1(i,j,k)-dfc1(i,up,k))/dx)/dx
     &        +dt*((m21(i,j,top)+m21(i,j,k))/2.0
     &        *(dfc1(i,j,top)-dfc1(i,j,k))/dx
     &        -(m21(i,j,k)+m21(i,j,bottom))/2.0
     &        *(dfc1(i,j,k)-dfc1(i,j,bottom))/dx)/dx
     &        +dt*((m22(right,j,k)+m22(i,j,k))/2.0
     &        *(dfc2(right,j,k)-dfc2(i,j,k))/dx
     &        -(m22(i,j,k)+m22(left,j,k))/2.0
     &        *(dfc2(i,j,k)-dfc2(left,j,k))/dx)/dx
     &        +dt*((m22(i,down,k)+m22(i,j,k))/2.0
     &        *(dfc2(i,down,k)-dfc2(i,j,k))/dx
     &        -(m22(i,j,k)+m22(i,up,k))/2.0
     &        *(dfc2(i,j,k)-dfc2(i,up,k))/dx)/dx
     &        +dt*((m22(i,j,top)+m22(i,j,k))/2.0
     &        *(dfc2(i,j,top)-dfc2(i,j,k))/dx
     &        -(m22(i,j,k)+m22(i,j,bottom))/2.0
     &        *(dfc2(i,j,k)-dfc2(i,j,bottom))/dx)/dx

            
      end do
      end do
      end do

!$OMP END DO   
!$OMP END PARALLEL

*     End of on step system evolution 

      return
      end

      function gasdev()

      data iset /0/
      data gset /0./

      common /iseedd/iseed1,iseed2

      if(iset.eq.0) then

cc  NZHOU MODIFIED              
 111     call random_number(randtest)
         v1=2*randtest-1
         call random_number(randtest)
         v2=2*randtest-1
         r=v1*v1+v2*v2
         if(r.ge.1) goto 111
         fac=sqrt(-2*log(r)/r)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
      else
         gasdev=gset
         iset=0
      end if
      return
      end


      subroutine interfacecalc()
      implicit real (a-h,o-z)
      include 'parameterc.inc'


      gamma=gamma/L0*1e-2
      kapae=lamda*gamma*8/3.1415926/3.1415926/gnormal
      wmega=gamma/lamda*4

      Tao=L0*L0/Bnormal/Gnormal*1.0e-12
   
      print *, kapae, wmega, gamma, lamda
      print *, 'time scale', Tao, 's'

      return
      end


      subroutine get_ML(ca01,ca02,cb01,cb02)
      implicit real (a-h,o-z)
      include "parameterc.inc"

      real kesi,kesi_co1, kesi_co2,sum
      real lamda1, mob11, mob22
      integer i,j
      real ca01,ca02,cb01,cb02



c+++++++++++++++++++++++++ mobility of Cr (B1)++++++++++++++++++++++

     

       QALTI1=204000
       QALTI2=96000
       FALTI1=5.51E-6
       FALTI2=5.19E-10

       QALAL1=19000
       FALAL1=5E-4

       QALV1=290000
       FALV1=1.85E-5

       QTITI1=121000
       QTITI2=237000
       FTITI1=1.47E-8
       FTITI2=5.91E-5

       QVTI1=107000
       QVTI2=212000
       FVTI1=1.55E-9
       FVTI2=2.54E-5

       QTIV1=320000
       QTIV2=544000
       FTIV1=1.62E-4
       FTIV2=51

       QVV1=306000
       QVV2=493000
       FVV1=2.73E-5
       FVV2=1.58

       ALTITIV0=111000-20*T
       ALTITIV1=-43000+10*T
       ALVTIV0=78000+11*T
       ALVTIV1=47000-22*T

       ALTIALTI0=-247000-83*T
       ALTIALTI1=-146000
       ALALALTI0=-1059000+373*T
       ALALALTI1=-921000+342*T

       ALALALV0=-4143000
       ALVALV0=-5752000

       ALALTIV0=+180000-277*T


       VVAL12=-0.12027*QALTI1*T**(-1)
       VVAL14=-0.12027*QALTI2*T**(-1)
       DTALTI=FALTI1*EXP(VVAL12)+FALTI2*EXP(VVAL14)

       DTALAL=-QALAL1+R*T*LOG(FALAL1)

       DTALV =-QALV1+R*T*LOG(FALV1)

       VV12=-0.12027*QTITI1*T**(-1)
       VV14=-0.12027*QTITI2*T**(-1)
       DTTITI=FTITI1*EXP(VV12)+FTITI2*EXP(VV14)
       VV22=-0.12027*QTIV1*T**(-1)
       VV24=-0.12027*QTIV2*T**(-1)
       DTTIV=FTIV1*EXP(VV22)+FTIV2*EXP(VV24)
       VV32=-0.12027*QVTI1*T**(-1)
       VV34=-0.12027*QVTI2*T**(-1)
       DTVTI=FVTI1*EXP(VV32)+FVTI2*EXP(VV34)
       VV42=-0.12027*QVV1*T**(-1)
       VV44=-0.12027*QVV2*T**(-1)
       DTVV=FVV1*EXP(VV42)+FVV2*EXP(VV44)

       AMTITIHCP=-3.03E5+R*T*LOG(1.35E-3)
       AMALTIHCP=-3.29E5+R*T*LOG(6.6E-3)
       AMVTIHCP=-2.50E5+R*T*LOG(1E-3)

c 1=al 2=v 3=ti
c+++++++++++++++++++++++++ mobility of Al (B2)+++++++++++++++++++++++

      G1_1=DTALAL
      G1_2=DTALV
      G1_3=R*T*LOG(DTALTI)
      
      G1_13_0=ALALALTI0
      G1_13_1=ALALALTI1
      G1_32_0=ALALTIV0
      G1_12_0=ALALALV0

c+++++++++++++++++++++++++ mobility of V (B2)+++++++++++++++++++++++

      G2_1=DTALAL
      G2_2=R*T*LOG(DTVV)
      G2_3=R*T*LOG(DTVTI)
         
      G2_12_0=ALVALV0
      G2_32_0=ALVTIV0
      G2_32_1=ALVTIV1

c+++++++++++++++++++++++++ mobility of TI (B3)+++++++++++++++++++++++

      G3_1=DTALAL
      G3_2=R*T*LOG(DTTIV)
      G3_3=R*T*LOG(DTTITI)

      G3_13_0=ALTIALTI0
      G3_13_1=ALTIALTI1
      G3_32_0=ALTITIV0
      G3_32_1=ALTITIV1


c+++++++++++++++++++++++++ mobility of Al (B2)+++++++++++++++++++++++

      HG1_1=AMALTIHCP
      HG1_2=AMALTIHCP
      HG1_3=AMALTIHCP
      
c+++++++++++++++++++++++++ mobility of V (B2)+++++++++++++++++++++++

      HG2_1=AMVTIHCP
      HG2_2=AMVTIHCP
      HG2_3=AMVTIHCP
         

c+++++++++++++++++++++++++ mobility of TI (B3)+++++++++++++++++++++++

      HG3_1=AMTITIHCP
      HG3_2=AMTITIHCP
      HG3_3=AMTITIHCP

      lamda1=3.1415926*0.5*sqrt(kapae*gnormal*0.5/wmega)
      
      sum=0
      kesi=0
      j=1
      do i=-10,10
         if ( abs(i) .le. lamda1+0.1) then
         ett=0.5*(1.0-sin(sqrt(2*wmega/gnormal/kapae)*i))
         pfeta=ett*ett*ett*(10.0-15.0*ett+6.0*ett*ett)
         x1=pfeta*ca01+(1.0-pfeta)*cb01
         x2=pfeta*ca02+(1.0-pfeta)*cb02
         x3=1.0-x1-x2

         DG1_B=x1*G1_1+x2*G1_2+x3*G1_3+x1*x2*G1_12_0+x3*x2*G1_32_0
     &        +x1*x3*(G1_13_0+G1_13_1*(X1-X3))

         B1_B=exp(DG1_B/R/T)/R/T/Bnormal
c------------------------- mobility of V (B2) ---------------------------

         DG2_B=x1*G2_1+x2*G2_2+x3*G2_3+x1*x2*G2_12_0+x3*x2*(G2_32_0
     &        +G2_32_1*(X3-X2))

         B2_B=exp(DG2_B/R/T)/R/T/Bnormal
c
c------------------------- mobility of TI (B3) ---------------------------

         DG3_B=x1*G3_1+x2*G3_2+x3*G3_3+x3*x2*(G3_32_0+
     &           G3_32_1*(X3-X2))+x1*x3*(G3_13_0+G3_13_1*(X1-X3))

         B3_B=exp(DG3_B/R/T)/R/T/Bnormal

C-----------------          HCP           ---------------------------

         DG1_H=HG1_1

         B1_H=exp(DG1_H/R/T)/R/T/Bnormal

c------------------------- mobility of V (B2) ---------------------------

         DG2_H=HG2_1

         B2_H=exp(DG2_H/R/T)/R/T/Bnormal

c------------------------- mobility of TI (B3) ---------------------------

         DG3_H=HG3_1

         B3_H=exp(DG3_H/R/T)/R/T/Bnormal


         B1=EXP((DG1_H+ETT*(DG1_B-DG1_H))/R/T)/R/T/Bnormal
         B2=EXP((DG2_H+ETT*(DG2_B-DG2_H))/R/T)/R/T/Bnormal
         B3=EXP((DG3_H+ETT*(DG3_B-DG3_H))/R/T)/R/T/Bnormal

         B1=-B1+B1_H+B1_B
         B2=-B2+B2_H+B2_B
         B3=-B3+B3_H+B3_B

         AA=-(1-x1)*x1*B1*(1-B3/B1)+x1*x2*B2*(1-B3/B2)
         EE1=-(1-x2)*x2*B2*(1-B3/B2)+x1*x2*B1*(1-B3/B1)
         mob11=((1-x1)*x1*B1+AA*x1)
         mob22=((1-x2)*x2*B2+EE1*x2)
    

         kesi=kesi+(cb01-ca01)*(x1-ca01)*pfeta/mob11
         kesi=kesi+(cb02-ca02)*(x2-ca02)*pfeta/mob22

      end if 
      end do

      me=3.1415926*0.125*sqrt(2*wmega/gnormal/kapae)/kesi
      print *, me

      return
      end 


c  add by Rongpei Shi (02/17/2009)

       subroutine gaussian_series()
  
       implicit real (a-h,o-z)
       include 'parameterc.inc'
  
       real u1,u2
       real g1,g3
       real tp,randtest
       integer i,j,k,i1,i2
     
         i1=1

c111	 call random_seed()
111      call random_number(randtest)
         u1=randtest
c         call random_seed()
         call random_number(randtest)
         u2=randtest


         tp=sqrt(-2.0*log(u1))

         g1=tp*cos(TWOPI*u2)
         g3=tp*sin(TWOPI*u2)
 
      randnum(i1)=g1
      i2=i1+1
      randnum(i2)=g3
      i1=i1+2
      if (i1.lt.nxyz) goto 111
 
      do k=1,nz
      do j=1,ny
	  do i=1,nx
	 
       randnum3(i,j,k)=fc*randnum((k-1)*ny*nx+(j-1)*nx+i)
	  end do
	  end do
	  end do
 
	  return 
	  end
      subroutine init_random_seed()

      integer i,n,clock
      integer seed(8)

      call system_clock(COUNT=clock)
      
      do i=1,8
         seed(i)=clock+37*(i-1)

      enddo

      call random_seed(put=seed)

      end subroutine
      
