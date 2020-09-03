      subroutine spec(n,iz,xyz,echo,e,g)
      implicit none             
      real*8 autoang,autokcal
c coversion factors
      parameter (autoang =0.52917726d0)
      parameter (autokcal=627.509541d0)
      logical echo
c energy
      real*8 e      
c number of atoms
      integer n
c coordinates in au
      real*8 xyz(3,n)
c gradient
      real*8 g  (3,n)
c cardinal numbers of elements
      integer   iz(n)
c coordination numbers of the atoms
      real*8 cn(n)
c covalent radii
      real*8 rcov(94)
      real*8 rcovin(94)
c covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
c values for metals decreased by 10 %
      data rcovin /
     .  0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67
     ., 1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54
     ., 1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09
     ., 1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39
     ., 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26
     ., 1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57
     ., 1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53
     ., 1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32
     ., 1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58
     ., 1.52, 1.53, 1.54, 1.55 /
      real*8 p(0:94)
      real*8 k1
      integer i

      p=0
      open(unit=1,file='~/.ppot_param')
c     read(1,*) p(0)
      read(1,*) p(1)
      read(1,*) p(6)
      read(1,*) p(7)
      read(1,*) p(8)
      close(1)
      p(1)=p(1)/100.
      p(6)=p(6)/100.
      p(7)=p(7)/100.
      p(8)=p(8)/100.
      p(0)=10.0

c scale and convert to au
      rcov=(4./3.)*rcovin/autoang

      call cnspec(n,p(0),rcov,iz,xyz,cn)

      e=0
      do i=1,n
         write(*,*)i,iz(i),cn(i)
         e=e+cn(i)*p(iz(i))
      enddo

      end


      subroutine cnspec(natoms,k1,rcov,iz,xyz,cn)
      implicit none  
      integer iz(*),natoms,i
      real*8 xyz(3,*),cn(*),rcov(94),k1

      integer iat    
      real*8 dx,dy,dz,r,damp,xn,rr,rco,r2

      do i=1,natoms
      xn=0.0d0
      do iat=1,natoms
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r2=dx*dx+dy*dy+dz*dz 
            r=sqrt(r2)                  
c covalent distance in Bohr
            rco=rcov(iz(i))+rcov(iz(iat))
            rr=rco/r
c counting function exponential has a better long-range behavior than MHGs inverse damping
            damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
            xn=xn+damp
         endif
      enddo
      cn(i)=xn  
      enddo

      end
