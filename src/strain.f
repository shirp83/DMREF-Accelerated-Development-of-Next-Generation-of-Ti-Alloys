

      subroutine sym_op ()
      implicit none
      include 'parameterc.inc'

      integer i,j,k,p

      do p=1,6
      do j=1,3
      do i=1,3
         s(i,j,p)=0.0
      end do
      end do
      end do

      s(1,1,1)= 1.0
      s(2,2,1)= 1.0
      s(3,3,1)= 1.0

      s(1,1,2)= 1.0
      s(2,3,2)= 1.0
      s(3,2,2)=-1.0

      s(1,2,3)= sqrt(2.0)/2.0
      s(1,3,3)=-sqrt(2.0)/2.0
      s(2,1,3)=-sqrt(2.0)/2.0
      s(2,2,3)= 0.5
      s(2,3,3)= 0.5
      s(3,1,3)= sqrt(2.0)/2.0
      s(3,2,3)= 0.5
      s(3,3,3)= 0.5

      s(1,2,4)= sqrt(2.0)/2.0
      s(1,3,4)=-sqrt(2.0)/2.0
      s(2,1,4)= sqrt(2.0)/2.0
      s(2,2,4)= 0.5
      s(2,3,4)= 0.5
      s(3,1,4)= sqrt(2.0)/2.0
      s(3,2,4)=-0.5
      s(3,3,4)=-0.5      

      s(1,2,5)=-sqrt(2.0)/2.0
      s(1,3,5)=-sqrt(2.0)/2.0
      s(2,1,5)= sqrt(2.0)/2.0
      s(2,2,5)=-0.5
      s(2,3,5)= 0.5
      s(3,1,5)=-sqrt(2.0)/2.0
      s(3,2,5)=-0.5
      s(3,3,5)= 0.5 

      s(1,2,6)=-sqrt(2.0)/2.0
      s(1,3,6)=-sqrt(2.0)/2.0
      s(2,1,6)=-sqrt(2.0)/2.0
      s(2,2,6)=-0.5
      s(2,3,6)= 0.5
      s(3,1,6)=-sqrt(2.0)/2.0
      s(3,2,6)= 0.5
      s(3,3,6)=-0.5

      return
      end 


      subroutine polymisc()
      implicit none
      include 'parameterc.inc'
      
      real pi
      parameter(pi=3.1415926)
      integer i,j,k,l
      integer np,ip,p
      real phi1,phi2,phi
      integer ix,iy,iz
      real euler1(ga,3)



c======read euler angle

       OPEN(unit=10, file='euler.txt', status='unknown')

       do i=1,ga
          
          read(10,*)euler1(i,1),euler1(i,2),euler1(i,3)

       enddo


       close(10)

       print*,'euler angle read in'

      

*     applied stress, yy component (s22) only
      DO i=1,3
         DO j=1,3
         s_appij(i,j)=0.0
      end DO
      END DO
c        s_appij(1,1) =+0.027
        s_appij(3,3) = 0.977/2000

       DO i=1,3
          DO j=1,3
          write(*,*) 'applied stress: sigma(',i,j,')=',s_appij(i,j)
       END DO
       END DO
c       write(*,*)

        OPEN(unit=10, file='10graintotal', status='unknown')

      do k=1,nz
      do j=1,ny
      do i=1,nx


         do np=1,ga

            shapeg(i,j,k,np)=0.0

          enddo

         read(10,*)ip
         
         graineta(i,j,k)=ip

         shapeg(i,j,k,ip)=1.0

      enddo
      enddo
      enddo

      close(10)

       


c bicrystal

c$$$                do ix=1,nx
c$$$                do iy=1,ny
c$$$                do iz=1,nz
c$$$
c$$$                   if(ix.le.nx/3)then
c$$$
c$$$                      graineta(ix,iy,iz)=1
c$$$
c$$$                   elseif((ix.gt.nx/3).and.(ix.le.2*nx/3))then
c$$$
c$$$                      graineta(ix,iy,iz)=2
c$$$
c$$$                   else
c$$$                      
c$$$                      graineta(ix,iy,iz)=3
c$$$
c$$$                  
c$$$                  endif
c$$$
c$$$                  END DO
c$$$                  END DO
c$$$                  END DO


c assign rotation matrix for each grain
      do np=1,ga

      do j=1,3
      do i=1,3
         Qrot(i,j,np)=0
      end do
      end do

      end do

! can not be unit tensor?

      do i=1,ga


         phi1=euler1(i,1)*TWOPI/360.0
         
         phi = euler1(i,2)*TWOPI/360.0

         phi2 = euler1(i,3)*TWOPI/360.0

         

            Qrot(1,1,i)=cos(phi1)*cos(phi2)-sin(phi1)*sin(phi2)
     $                           *cos(phi)

            Qrot(1,2,i)=sin(phi1)*cos(phi2)+cos(phi1)*sin(phi2)
     $                           *cos(phi)

            Qrot(1,3,i)=sin(phi2)*sin(phi)

            Qrot(2,1,i)=-cos(phi1)*sin(phi2)-sin(phi1)
     $                           *cos(phi2)*cos(phi)

            Qrot(2,2,i)=-sin(phi1)*sin(phi2)+cos(phi1)
     $                           *cos(phi2)*cos(phi)

            Qrot(2,3,i)=cos(phi2)*sin(phi)

            Qrot(3,1,i)=sin(phi1)*sin(phi)

            Qrot(3,2,i)=-cos(phi1)*sin(phi)

            Qrot(3,3,i)=cos(phi)

      enddo

      return
      end


      subroutine elastic()
      implicit none
      include 'parameterc.inc'

      real c11,c12,c13,c33,c44
      real cp11,cp12,cp44,cp13,cp33
      real e0(3,3),e1(3,3),e(3,3,nvar),s0(3,3,nvar,ga), rot(3,3)
      real alpha,temp(nx,ny,nz)
      real omega(3,3),iomega(3,3)
      real interval, length, value(360), angle
      real self

      real eee(3,3,nvar)

      real e2(3,3),e3(3,3),s1(3,3,nvar,ga) 

      real T0(6,6),S0ij(6,6)
      integer ind(6)

      integer sa,sap,np,nq

      real grainrotation(ga,3,3)

      real euler1(ga,3)
      
      
      integer i,j,k,l
      integer p,q,m,n,pp
      integer ix, iy, iz
      integer index
      real phi1,phi,phi2
      real double

      real PI

      PI=3.14159265758
      
c Stress free strain for single crystal

      do i=1,3
         do j=1,3
            e0(i,j)=0
            e1(i,j)=0
            e2(i,j)=0
            e3(i,j)=0
         end do
      end do

c semi-coherent
      e0(1,1)= -0.049
      e0(1,2)= -3.04e-3
      e0(2,1)= -3.04e-3
      e0(2,2)= 0.067
      e0(3,3)=-2.998e-4

      e1(1,1)= -0.049
      e1(1,2)= -e0(1,2)
      e1(2,1)= -e0(2,1)
      e1(2,2)= e0(2,2)
      e1(3,3)=-2.998e-4

c coherent
      e2(1,1)= -0.083
      e2(1,2)= 9.486e-3
      e2(2,1)= 9.486e-3
      e2(2,2)= 0.123
      e2(3,3)= 0.035

      e3(1,1)= e2(1,1)
      e3(1,2)= -e2(1,2)
      e3(2,1)= -e2(2,1)
      e3(2,2)= e2(2,2)
      e3(3,3)= e2(3,3)



      do p=1,6
      do i=1,3
      do j=1,3
         e(i,j,p)=0.0
         eee(i,j,p)=0.0
         do k=1,3
         do l=1,3
            e(i,j,p)=e(i,j,p)+s(k,i,p)*s(l,j,p)*e0(k,l)
            eee(i,j,p)=eee(i,j,p)+s(k,i,p)*s(l,j,p)*e2(k,l)
         end do
         end do
      end do
      end do
      end do

      do p=7,12
      do i=1,3
      do j=1,3
         e(i,j,p)=0.0
         eee(i,j,p)=0.0
         do k=1,3
         do l=1,3
            e(i,j,p)=e(i,j,p)+s(k,i,p-6)*s(l,j,p-6)*e1(k,l)
            eee(i,j,p)=eee(i,j,p)+s(k,i,p-6)*s(l,j,p-6)*e3(k,l)
         end do
         end do
      end do
      end do
      end do


c SFTS of each grain in polycrystal

      do sa=1,ga

      do p=1,nvar
      do i=1,3
      do j=1,3

         ee(i,j,p,sa)=0.0
         ec(i,j,p,sa)=0.0


         do k=1,3
         do l=1,3
c semi-coherent
         ee(i,j,p,sa)=ee(i,j,p,sa)+Qrot(k,i,sa)*Qrot(l,j,sa)*e(k,l,p)
c coherent

         ec(i,j,p,sa)=ec(i,j,p,sa)+Qrot(k,i,sa)*Qrot(l,j,sa)*eee(k,l,p)

         end do
         end do

      end do
      end do
      end do !! end loop of p
      end do !! end loop of sa


* Elastic modulus

      DO l=1, 3
      DO k=1, 3
      DO j=1, 3
      DO i=1, 3
         C0ijkl(i,j,k,l) =0.0
         Cijkl_soft(i,j,k,l)=0.0
      END DO
      END DO
      END DO
      END DO  

*BCC
      c11=0.977
      c12=0.827
      c13=0.827
      c33=0.977
      c44=0.375

**   assign Cijkl  

      C0ijkl(1,1,1,1) =c11
      C0ijkl(2,2,2,2) =c11
      C0ijkl(3,3,3,3) =c33
      C0ijkl(1,1,2,2) =c12
      C0ijkl(2,2,3,3) =c13
      C0ijkl(3,3,1,1) =c13
      C0ijkl(1,2,1,2) =c44
      C0ijkl(2,3,2,3) =c44
      C0ijkl(3,1,3,1) =c44
      
      C0ijkl(2,2,1,1) =c12
      C0ijkl(3,3,2,2) =c13
      C0ijkl(1,1,3,3) =c13
      C0ijkl(2,1,2,1) =c44
      C0ijkl(1,2,2,1) =c44
      C0ijkl(2,1,1,2) =c44
      C0ijkl(3,2,3,2) =c44
      C0ijkl(3,2,2,3) =c44
      C0ijkl(2,3,3,2) =c44
      C0ijkl(1,3,1,3) =c44
      C0ijkl(1,3,3,1) =c44
      C0ijkl(3,1,1,3) =c44

      do i=1,3
      do j=1,3
         rot(i,j)=0
      end do
      end do

      rot(1,2)=-sqrt(2.0)/2.0
      rot(1,3)= sqrt(2.0)/2.0
      rot(2,1)= 1.0
      rot(3,2)= sqrt(2.0)/2.0
      rot(3,3)= sqrt(2.0)/2.0

      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
         Cijkl(i,j,k,l) =0.0
         do p=1,3
         do q=1,3
         do m=1,3
         do n=1,3
            Cijkl(i,j,k,l)=Cijkl(i,j,k,l)+
     &           rot(p,i)*rot(q,j)*rot(m,k)*rot(n,l)*
     &           C0ijkl(p,q,m,n)
         end do
         end do
         end do
         end do
      end do
      end do
      end do
      end do

c------reference modulus 2*Cijkl

      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3

         C0ijkl(i,j,k,l) =1.00*Cijkl(i,j,k,l)

      end do
      end do
      end do
      end do


c--------------calculate S0ijkl

      DO i=1, 3
      DO j=1, 3
      DO k=1, 3
      DO l=1, 3

          if(i.eq.1.and.j.eq.1) p=1
          if(i.eq.2.and.j.eq.2) p=2
          if(i.eq.3.and.j.eq.3) p=3

          if(i.eq.2.and.j.eq.3) p=4
          if(i.eq.3.and.j.eq.2) p=4

          if(i.eq.1.and.j.eq.3) p=5
          if(i.eq.3.and.j.eq.1) p=5

          if(i.eq.1.and.j.eq.2) p=6
          if(i.eq.2.and.j.eq.1) p=6


          if(k.eq.1.and.l.eq.1) q=1
          if(k.eq.2.and.l.eq.2) q=2
          if(k.eq.3.and.l.eq.3) q=3

          if(k.eq.2.and.l.eq.3) q=4
          if(k.eq.3.and.l.eq.2) q=4

          if(k.eq.1.and.l.eq.3) q=5
          if(k.eq.3.and.l.eq.1) q=5

          if(k.eq.1.and.l.eq.2) q=6
          if(k.eq.2.and.l.eq.1) q=6


         C0ij(p,q) =C0ijkl(i,j,k,l)
         Cij(p,q)=Cijkl(i,j,k,l)

      END DO
      END DO
      END DO
      END DO


** define Cij_soft as hexagonal

      do i=1,6
         do j=1,6
            Cij_soft(i,j)=1.0*Cij(i,j)
         enddo
      enddo
  
      DO i=1, 3
      DO j=1, 3
      DO k=1, 3
      DO l=1, 3

          if(i.eq.1.and.j.eq.1) p=1
          if(i.eq.2.and.j.eq.2) p=2
          if(i.eq.3.and.j.eq.3) p=3

          if(i.eq.2.and.j.eq.3) p=4
          if(i.eq.3.and.j.eq.2) p=4

          if(i.eq.1.and.j.eq.3) p=5
          if(i.eq.3.and.j.eq.1) p=5

          if(i.eq.1.and.j.eq.2) p=6
          if(i.eq.2.and.j.eq.1) p=6


          if(k.eq.1.and.l.eq.1) q=1
          if(k.eq.2.and.l.eq.2) q=2
          if(k.eq.3.and.l.eq.3) q=3

          if(k.eq.2.and.l.eq.3) q=4
          if(k.eq.3.and.l.eq.2) q=4

          if(k.eq.1.and.l.eq.3) q=5
          if(k.eq.3.and.l.eq.1) q=5

          if(k.eq.1.and.l.eq.2) q=6
          if(k.eq.2.and.l.eq.1) q=6

       Cijkl_soft(i,j,k,l)=Cij_soft(p,q)



      END DO
      END DO
      END DO
      END DO


c----------------------

* modulus start:stiffness C0

      write(*,*) 'C0'
      do i=1,6
           write(*,*) ( C0ij(i,j),j=1,6)
      END DO

      do i=1,6
          dO j=1,6
           T0(i,j)=C0ij(i,j)
          END DO
      END DO

*  do inverse of 6*6 matrix
      CALL MIGS (T0,6,S0ij,ind)

*  compliance S0
      write(*,*) 'S0'
      do i=1,6
            write(*,*) ( S0ij(i,j),j=1,6)
       END DO


c---------S0ijkl

      DO i=1, 3
      DO j=1, 3
      DO k=1, 3
      DO l=1, 3

          if(i.eq.1.and.j.eq.1) p=1
          if(i.eq.2.and.j.eq.2) p=2
          if(i.eq.3.and.j.eq.3) p=3

          if(i.eq.2.and.j.eq.3) p=4
          if(i.eq.3.and.j.eq.2) p=4

          if(i.eq.1.and.j.eq.3) p=5
          if(i.eq.3.and.j.eq.1) p=5

          if(i.eq.1.and.j.eq.2) p=6
          if(i.eq.2.and.j.eq.1) p=6


          if(k.eq.1.and.l.eq.1) q=1
          if(k.eq.2.and.l.eq.2) q=2
          if(k.eq.3.and.l.eq.3) q=3

          if(k.eq.2.and.l.eq.3) q=4
          if(k.eq.3.and.l.eq.2) q=4

          if(k.eq.1.and.l.eq.3) q=5
          if(k.eq.3.and.l.eq.1) q=5

          if(k.eq.1.and.l.eq.2) q=6
          if(k.eq.2.and.l.eq.1) q=6
          
          double=1
          if(p.gt.3) double=double*2
          if(q.gt.3)double=double*2

       S0ijkl(i,j,k,l)=S0ij(p,q)/double

      END DO
      END DO
      END DO
      END DO

c======read euler angle

       OPEN(unit=10, file='euler.txt', status='unknown')

       do i=1,ga
          
          read(10,*)euler1(i,1),euler1(i,2),euler1(i,3)

       enddo


       close(10)

       print*,'euler angle read in'


         do i=1,ga


         phi1=euler1(i,1)*TWOPI/360.0
         
         phi = euler1(i,2)*TWOPI/360.0

         phi2 = euler1(i,3)*TWOPI/360.0
     
  
         
        write(*,*),"orientation:",i,"  Euler angle:",phi1*180/PI,
     $  phi*180/PI,phi2*180/PI



            grainrotation(i,1,1)=cos(phi1)*cos(phi2)-sin(phi1)*sin(phi2)
     $                           *cos(phi)

            grainrotation(i,1,2)=sin(phi1)*cos(phi2)+cos(phi1)*sin(phi2)
     $                           *cos(phi)

            grainrotation(i,1,3)=sin(phi2)*sin(phi)

            grainrotation(i,2,1)=-cos(phi1)*sin(phi2)-sin(phi1)
     $                           *cos(phi2)*cos(phi)

            grainrotation(i,2,2)=-sin(phi1)*sin(phi2)+cos(phi1)
     $                           *cos(phi2)*cos(phi)

            grainrotation(i,2,3)=cos(phi2)*sin(phi)

            grainrotation(i,3,1)=sin(phi1)*sin(phi)

            grainrotation(i,3,2)=-cos(phi1)*sin(phi)

            grainrotation(i,3,3)=cos(phi)
        


       write(*,*) "rotation matrix of:",i

           write(*,*) grainrotation(i,1,1),grainrotation(i,1,2),
     $grainrotation(i,1,3)
           write(*,*) grainrotation(i,2,1),grainrotation(i,2,2),
     $grainrotation(i,2,3)
           write(*,*) grainrotation(i,3,1),grainrotation(i,3,2),
     $grainrotation(i,3,3)





                   do pp=1,3
                    do j=1,3
                    do k=1,3
                     do l=1,3

         Cijkl_grain(i,pp,j,k,l)=0.0


                        do p=1,3
                           do q=1,3
                              do m=1,3
                                 do n=1,3

c the Euler angle set /rotation matrix bring the sample axis 
c to be coincidence with crystal axis
c the rotation matrix should be inverse or transpose of the g
								 
         Cijkl_grain(i,pp,j,k,l)= Cijkl_grain(i,pp,j,k,l)+
     !  grainrotation(i,p,pp)* grainrotation(i,q,j)*
     ! grainrotation(i,m,k)* grainrotation(i,n,l)*Cijkl(p,q,m,n)


                     END DO
                     END DO
                     END DO
                     END DO


                     END DO
                     END DO
                     END DO
                     END DO


         END DO


c$$$         stop



      return
      end      

      SUBROUTINE Grad()
      IMPLICIT none
      INCLUDE 'parameterc.inc'

      INTEGER ix,iy,iz

      DO iz =1, nz
      DO iy =1, ny
      DO ix =1, nx

         if(ix .LE. nx/2+1) THEN
            g(1,ix,iy,iz) =TWOPI/nx/dx *(ix-1)
         ELSE
            g(1,ix,iy,iz) =TWOPI/nx/dx *(ix-1-nx)
         END IF
         if(iy .LE. ny/2+1) THEN
            g(2,ix,iy,iz) =TWOPI/ny/dx *(iy-1)
         ELSE
            g(2,ix,iy,iz) =TWOPI/ny/dx *(iy-1-ny)
         END IF
         if(iz .LE. nz/2+1) THEN
            g(3,ix,iy,iz) =TWOPI/nz/dx *(iz-1)
         ELSE
            g(3,ix,iy,iz) =TWOPI/nz/dx *(iz-1-nz)
         END IF

         g2(ix,iy,iz) =g(1,ix,iy,iz)**2 + g(2,ix,iy,iz)**2 +
     !                 g(3,ix,iy,iz)**2
      END DO
      END DO
      END DO  

      open(50,file='gradient.dat',status='unknown')
      do ix=1,nx
         do iy=1,ny
            do iz=1,nz
               write (50,*) g2(ix,iy,iz)
            end do
         end do
      end do

      close(50)





      return
      end

*-------------------------------------------------------
* subroutine DoInverse
*    3x3 matrix inverse
*-------------------------------------------------------

      SUBROUTINE DoInverse(a,ra) 

      real a(3,3),ra(3,3),b(3,3)
      real det_a
      integer i,j

      det_a = a(1,1)*a(2,2)*a(3,3)
     &      + a(2,1)*a(3,2)*a(1,3)+a(1,2)*a(2,3)*a(3,1)
     &      - a(3,1)*a(2,2)*a(1,3) 
     &      - a(2,1)*a(1,2)*a(3,3)-a(3,2)*a(2,3)*a(1,1)

      b(1,1)=a(2,2)*a(3,3)-a(3,2)*a(2,3)
      b(1,2)=a(1,2)*a(3,3)-a(3,2)*a(1,3)
      b(1,3)=a(1,2)*a(2,3)-a(2,2)*a(1,3)
      b(2,1)=a(2,1)*a(3,3)-a(3,1)*a(2,3)
      b(2,2)=a(1,1)*a(3,3)-a(3,1)*a(1,3)
      b(2,3)=a(1,1)*a(2,3)-a(2,1)*a(1,3)
      b(3,1)=a(2,1)*a(3,2)-a(3,1)*a(2,2)
      b(3,2)=a(1,1)*a(3,2)-a(3,1)*a(1,2)
      b(3,3)=a(1,1)*a(2,2)-a(2,1)*a(1,2)

      do i=1,3
      do j=1,3
         if(mod(i+j,2).eq.0)then
            ra(i,j)=b(i,j)/det_a 
         else
           ra(i,j)=-b(i,j)/det_a
         endif
      enddo
      enddo

      return
      end

c-------------     
!
      SUBROUTINE MIGS (A,N,X,INDX)
!
! Subroutine to invert matrix A(N,N) with the inverse stored
! in X(N,N) in the output.  Copyright (c) Tao Pang 2001.
!
       IMPLICIT NONE
       INTEGER N
       INTEGER I,J,K
       INTEGER  INDX(6)
         real   A(6,6)
             real  X(6,6)
         real    B(6,6)
!
    

       DO I = 1, N
       DO J = 1, N
        B(I,J) = 0.0
        END DO
       END DO
         DO I = 1, N
          B(I,I) = 1.0
        END DO
!   

        CALL ELGS (A,N,INDX)
!




              DO I = 1, N-1
          DO J = I+1, N
            DO K = 1, N
              B(INDX(J),K) = B(INDX(J),K)-A(INDX(J),I)*B(INDX(I),K)
            END DO
          END DO
        END DO
!
        DO I = 1, N
          X(N,I) = B(INDX(N),I)/A(INDX(N),N)
          DO J = N-1, 1, -1
            X(J,I) = B(INDX(J),I)
            DO K = J+1, N
              X(J,I) = X(J,I)-A(INDX(J),K)*X(K,I)
            END DO
                  X(J,I) =  X(J,I)/A(INDX(J),J)
          END DO
        END DO
      END SUBROUTINE MIGS
!




      SUBROUTINE ELGS (A,N,INDX)
!
! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed matrix
! plus the pivoting element ratios below the diagonal in the output.
! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
!
        IMPLICIT NONE
        INTEGER N
        INTEGER I,J,K,ITMP
        INTEGER INDX(6)
        REAL :: C1,PI,PI1,PJ,tp1,tp2
        real  A(6,6)
        real   C(6)
!
! Initialize the index
!

 
        DO I = 1, N
          INDX(I) = I
        END DO
!  

! Find the rescaling factors, one from each row
!
  

      DO I = 1, N
          C1= 0.0
          DO J = 1, N

             tp2=ABS(A(I,J))

            C1 = AMAX1(C1,tp2)
          END DO
          C(I) = C1
        END DO
!


! Search the pivoting (largest) element from each column
!
        DO J = 1, N-1
          PI1 = 0.0
          DO I = J, N
            PI = ABS(A(INDX(I),J))/C(INDX(I))
            IF (PI.GT.PI1) THEN
              PI1 = PI
              K   = I
            ENDIF
          END DO
!
! Interchange the rows via INDX(N) to record pivoting order
!
          ITMP    = INDX(J)
          INDX(J) = INDX(K)
                INDX(K) = ITMP
                DO I = J+1, N
            PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
            A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
            DO K = J+1, N
              A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
            END DO
          END DO
        END DO
!
      END SUBROUTINE ELGS

