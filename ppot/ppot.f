      subroutine ppot(edisp,egcp,ebas)
      use defvar
      implicit none
c cart coords
      real*8 ,allocatable ::xyz (:,:)
c grad      
      real*8 ,allocatable ::g   (:,:)
      real*8 ,allocatable ::g1  (:,:)
      real*8 ,allocatable ::g2  (:,:)
      real*8 ,allocatable ::g3  (:,:)
      real*8 ,allocatable ::g4  (:,:)
c atomic ordinal number
      integer,allocatable ::at  (:)
c # atoms
      integer n
c # logical if tm files are written
      logical tmgrad      
c energy 
      real*8 edisp,egcp,ebas,eat,etot      
c D3 parameter
      real*8 rs6,s18,rs18      
c conversion factor
      real*8 autokcal 
      parameter (autokcal=627.509541d0)
c for D3 (H-Pu)       
      integer max_elem
      parameter (max_elem=94)
c pbex hybrid ?
      logical pbex, echo

c local stuff
      character*80 fname      
      character*128 atmp,btmp       
      character*2 asym              
      integer i
      logical ex
 
      tmgrad=.false.
      pbex  =.false.
      echo=.false.
      if (doGrimme==2) echo=.true.

      write(*,*)'+--------------------------------------------+)'
      write(*,*)'|             P P O T  (SG, Mar 12)          |)'
      write(*,*)'|   corrects for D3, BSSE by gcp and basis   |)'
      write(*,*)'|                      -                     |)'
      write(*,*)'|     slightly modified version by RP for    |)'
      write(*,*)'|          interaction with DFT:IoG          |)'
      write(*,*)'+--------------------------------------------+)'
      write(*,*)'please cite'
      write(*,*)'a) DFT-D3:'
      write(*,*)'S. Grimme, J. Antony, S. Ehrlich and H. Krieg,'
      write(*,*)'J. Chem. Phys. 132 (2010), 154104'
      write(*,*)'S. Grimme, S. Ehrlich and L. Goerigk,'
      write(*,*)'J. Comput. Chem. 32 (2011), 1456-1465'
      write(*,*)'b) gCP:'
      write(*,*)'H. Kruse and S. Grimme,'
      write(*,*)'J. Chem. Phys. 134 (2012), 154101.'
      write(*,*)'c) SRB correction (HF-3c):'
      write(*,*)'R. Sure and S. Grimme,'
      write(*,*)'J. Comput. Chem., 34 (2013), 1672.'

c how many atoms?
      n=ncenter
      allocate(xyz(3,n),at(n),g(3,n)) 
      allocate(g1(3,n),g2(3,n),g3(3,n),g4(3,n)) 

c read coord in TM format


      do i=1,n
          at(i)=a(i)%index
          xyz(1,i)=a(i)%x/b2a
          xyz(2,i)=a(i)%y/b2a
          xyz(3,i)=a(i)%z/b2a
      end do

 
      if (echo) then
      write(*,*)
      write(*,*)'coordinates (Bohr)'
      do i=1,n
         write(*,'(i4,2x,a2,3x,3f12.5)')i,asym(at(i)),xyz(1:3,i)
      enddo
      end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3D3         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      rs6 =0.3981
      s18 =2.5
      rs18=4.4211  
      if(pbex)then
c pbe0/svx              
c     rs6=0.144
c     s18=0.701
c     rs18=5.431     
c pbe38/svx              
c     rs6=0.000
c     s18=0.708
c     rs18=6.270    
c pbe12/svx              
      rs6=0.000
      s18=2.570
      rs18=7.340    
      endif

      call dftd3(n,at,xyz,echo,.false.,rs6,s18,rs18,edisp,g1)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c gCPgCPgCPgCPgCPgCPgCPgCPgCPgCPgCPgCPgCPgCPgCPgCPgCPgCPgCPgCP      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      atmp='hf/minix'
      if(pbex)atmp='dft/svx'

      call gcp(n,at,xyz,echo,atmp,egcp,g2)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c BASBASBASBASBASBASBASBASBASBASBASBASBASBASBASBASBASBASBASBAS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(.not.pbex)then
      call basegrad(n,echo,max_elem,at,xyz,0.7d0,0.03d0,ebas,g3)
      else
c     call basegrad(n,max_elem,at,xyz,1.0d0,0.100d0,ebas,g3)
      g3=0
      ebas=0
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     call spec(n,at,xyz,.true.,eat,g4)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      g=g1+g2+g3

      etot=edisp+egcp+ebas

      write(*,*)
      write(*,"(' ---=< ppot results >=--- ')")
      write(*,'('' Edisp [a.u.]:'',f20.10)')edisp
      write(*,'('' Egcp  [a.u.]:'',f20.10)')egcp  
      write(*,'('' Ebas  [a.u.]:'',f20.10)')ebas 
c     write(*,'('' Eat   [a.u.]:'',f12.8)')eat  
c      write(*,*)
c      write(*,*)'|G|=',sum(abs(g(1:3,1:n)))
      write(*,"('              --------------------')")
      
      write(*,'('' Eppot  [a.u.]:'',f20.10)')etot 
      write(*,*)


      end subroutine
