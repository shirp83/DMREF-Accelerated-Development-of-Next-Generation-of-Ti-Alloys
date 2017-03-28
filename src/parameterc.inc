c     system size
      integer nx,ny,nz,nxyz
      parameter (nx=128,ny=128, nz=128)
      parameter (nxyz=nx*ny*nz)
      integer nvar, node,ga
      parameter (nvar=12)
      parameter (node=nvar)
      parameter (ga=10)
      real T,R
      parameter(T=1070.0, R=8.31451)
      real Gnormal, Bnormal
      parameter(Gnormal=5000.0, Bnormal=1.0E-18)
      real TWOPI
      parameter(TWOPI=3.14159265758*2.0)
      real dx
      parameter (dx=1.0)
      real epsilon
      parameter (epsilon=1.0e-4)

c     field variable 
      real c1(nx,ny,nz),dfc1(nx,ny,nz)
      real c2(nx,ny,nz),dfc2(nx,ny,nz)
      real eta(nx,ny,nz,ga,node)
       

      real dfe(nx,ny,nz,ga,node)
      real g2(nx,ny,nz), g(3,nx,ny,nz)
      real s(3,3,6)
      complex appk(nx,ny,nz,ga,nvar)
      complex appk1(nx,ny,nz,ga,nvar)

c     randnoise term
      
      real fc
      real randnum(nxyz)
      real randnum3(nx,ny,nz)
      complex randnum3k(nx,ny,nz)

c     Rotation Matrix assigned to each grain

      real Qrot(3,3,ga)

c     shape function of grain

      real shapeg(nx,ny,nz,ga)
      integer graineta(nx,ny,nz)


c     integration time step and total time
      real dt

c     gradient energy coefficient and order parameter mobility
      real me, df, kapac, mc

c     anisotropy degree, symmetry, rotation angle, missing angle
      real agratio, symt, theta0, thetam, sthetam, gthetam

      real dfcc_hcp1, dfcc_hcp2, dfcc_bcc1, dfcc_bcc2

      real wmega, kapae,tao
      real gamma, L0
      integer lamda

      real Cijkl(3,3,3,3), C0ijkl(3,3,3,3)
      real S0ijkl(3,3,3,3)
      real Cijkl_soft(3,3,3,3)
      real s_appij(3,3)
      real Cijkl_grain(ga,3,3,3,3)
      real ee(3,3,nvar,ga)
      real ec(3,3,nvar,ga)


      real C0ij(6,6),Cij(6,6),Cij_soft(6,6)

      real eij(3,3,nx,ny,nz)



      integer nthreads
      	
      common /poly1/ Qrot
      common /poly2/ shapeg,graineta,Cijkl_grain
      common /numerical/ dt
      common /field/ eta
      common /f1/  dfe
      common /f2/  c1
      common /f3/  c2
      common /f4/ dfc1
      common /f5/ dfc2
      common /anisot/ agratio, symt, theta0, thetam, sthetam, gthetam
      common /elas/  s, g, g2,S0ijkl,eij
      common /elas3/ Cijkl,C0ijkl,s_appij,ee,ec
      common /elas4/ Cijkl_soft,C0ij,Cij,Cij_soft
      common /eq/ dfcc_hcp1, dfcc_hcp2, dfcc_bcc1, dfcc_bcc2
      common /popenmp/ nthreads
      common /interface/ wmega, kapae, gamma, L0, tao,me, lamda
      common /applied1/ appk
      common /applied3/ appk1
      common /noise1/ fc,randnum
      common /noise2/ randnum3
      common /noise3/ randnum3k
