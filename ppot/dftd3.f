      subroutine dftd3(n,iz,xyz,echo,grad,rs6,s18,rs18,disp,g)
      implicit none             
      integer max_elem,maxc                      
c conversion factors
      real*8 autoang,autokcal,c6conv
      parameter (max_elem=94)
c maximum coordination number references per element
      parameter (maxc    =5)
c coversion factors
      parameter (autoang =0.52917726d0)
      parameter (autokcal=627.509541d0)
c DFT-D version
      integer version
c number of atoms
      integer n
c coordinates in au
      real*8 xyz(3,n)
c gradient
      real*8 g  (3,n)
c cardinal numbers of elements
      integer   iz(n)
c cut-off radii for all element pairs
      real*8 r0ab(max_elem,max_elem)
c C6 for all element pairs 
      real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
c how many different C6 for one element
      integer mxc(max_elem)
c C6810 
      real*8 c6,c8,c10
c coordination numbers of the atoms
      real*8 cn(n)
c covalent radii
      real*8 rcov(max_elem)
      real*8 rcovin(max_elem)
c atomic <r^2>/<r^4> values
      real*8 r2r4(max_elem)
      real*8 r2r4in(max_elem)
c energies
      real*8 e6, e8, e10, e12, disp, e6abc        
c THE PARAMETERS OF THE METHOD (not all a "free")
      real*8 rs6, rs8, rs10, s6, s18, alp6, alp8, alp10, s42, rs18, alp
c printout option
      logical echo
c grad ?
      logical grad
c analyse results ?
      logical anal
c third-order term?
      logical noabc
c gradient calctype
      logical numgrad
c special parameters
      logical tz
c O(N^2)  neglect threshold (save choice, for O(N^3) threshold see gdisp.f)
      real*8 rthr,rthr2
      parameter (rthr=40000.0d0)
      parameter (rthr2=3000.0d0)

c local and dummy variables
      character*80 atmp,btmp,ctmp,dtmp,etmp,ftmp,func
      character*2  esym 
      integer i,j,z,nn,iat,jat,i1,i2,k
      integer ida(max_elem),ipot,scomp,list(6)  
      real*8  x,y,dispr,displ,gdsp,dum,xx(10),dum6(86)
      real*8  dum1,dum2
      logical ex,pot

c PBE0/def2-QZVP atomic values 
      data r2r4in /
     .  8.0589,  3.4698, 29.0974, 14.8517, 11.8799,  7.8715,  5.5588,
     .  4.7566,  3.8025,  3.1036, 26.1552, 17.2304, 17.7210, 12.7442,
     .  9.5361,  8.1652,  6.7463,  5.6004, 29.2012, 22.3934, 19.0598,
     . 16.8590, 15.4023, 12.5589, 13.4788, 12.2309, 11.2809, 10.5569,
     . 10.1428,  9.4907, 13.4606, 10.8544,  8.9386,  8.1350,  7.1251,
     .  6.1971, 30.0162, 24.4103, 20.3537, 17.4780, 13.5528, 11.8451,
     . 11.0355, 10.1997,  9.5414,  9.0061,  8.6417,  8.9975, 14.0834,
     . 11.8333, 10.0179,  9.3844,  8.4110,  7.5152, 32.7622, 27.5708,
     . 23.1671, 21.6003, 20.9615, 20.4562, 20.1010, 19.7475, 19.4828,
     . 15.6013, 19.2362, 17.4717, 17.8321, 17.4237, 17.1954, 17.1631,
     . 14.5716, 15.8758, 13.8989, 12.4834, 11.4421, 10.2671,  8.3549,
     .  7.8496,  7.3278,  7.4820, 13.5124, 11.6554, 10.0959,  9.7340,
     .  8.8584,  8.0125, 29.8135, 26.3157, 19.1885, 15.8542, 16.1305,
     . 15.6161, 15.1226, 16.1576 /                                       
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

c k1-k3
      include 'param'
c scale and convert to au
      rcov=k2*rcovin/autoang
c init
      version=4
      s6  =1.0d0
      alp =14.0d0

      pot =.false.
      anal=.false.
      noabc=.true. 
      numgrad=.false.
      tz=.false.
      func=' HF-D3(BJ)-gCP/MINIX '
c J/mol nm^6 - > au
      c6conv=1.d-3/2625.4999d0/((0.052917726d0)**6)

c set radii
      call setr0ab(max_elem,autoang,r0ab)

c C6 hard-coded (c6ab.dat not used)
c this is alternative to loadc6
      call copyc6(btmp,maxc,max_elem,c6ab,mxc)         

c scale r4/r2 values of the atoms by sqrt(Z) 
c sqrt is also globally close to optimum
c together with the factor 1/2 this yield reasonable
c c8 for he, ne and ar. for larger Z, C8 becomes too large
c which effectively mimics higher R^n terms neglected due
c to stability reasons

      do i=1,max_elem
         dum    =0.5*r2r4in(i)*dfloat(i)**0.5   
c store it as sqrt because the geom. av. is taken
         r2r4(i)=sqrt(dum)                         
      enddo

c for global ad hoc parameters see
c k3 in subroutine getc6, k1 and k2 in subroutine ncoord*
c fixed or dependent ones:
      rs8  = rs18       
      rs10 = rs18
      alp6 = alp
      alp8 = alp+2.
      alp10= alp8+2. 
c note: if version=4 (Becke-Johnson), a1=rs6 and a2=rs18
c       and alp* have no meaning

cccccccccccccc
c energy call
cccccccccccccc
      call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
     .     rs6,rs8,rs10,alp6,alp8,alp10,version,.true.,rthr,
     .     e6,e8,e10,e12,e6abc)

      e6   = e6   *s6

      e8   = e8   *s18

      disp =-e6-e8

c output
      if (echo) then
      write(*,*)'---------------------------------------------------'        
      if(version.lt.4)then
      write(*,'(''                  DFT-D V'',i1)') version       
      else
      write(*,'(''                  DFT-D V3(BJ)'')') 
      endif
      write(*,*)'---------------------------------------------------'        
      write(*,'('' parameters'')') 
      if(version.eq.2)then
         write(*,'('' s6       :'',f10.4)') s6            
         write(*,'('' alpha6   :'',f10.4)') alp6        
      endif
      if(version.eq.3)then
         write(*,'('' s6       :'',f10.4)') s6            
         write(*,'('' s8       :'',f10.4)') s18           
         write(*,'('' rs6      :'',f10.4)') rs6           
         write(*,'('' rs18     :'',f10.4)') rs18          
         write(*,'('' alpha6   :'',f10.4)') alp6        
         write(*,'('' alpha8   :'',f10.4)') alp8           
         write(*,'('' k1-k3    :'',3f10.4)') k1,k2,k3     
      endif
      if(version.eq.4)then
         write(*,'('' s6       :'',f10.4)') s6            
         write(*,'('' s8       :'',f10.4)') s18           
         write(*,'('' a1       :'',f10.4)') rs6           
         write(*,'('' a2       :'',f10.4)') rs18          
         write(*,'('' k1-k3    :'',3f10.4)') k1,k2,k3     
      endif
      write(*,'(/'' E6    /kcal :'',f11.4)')-e6*autokcal
      if(version.gt.2)then
      write(*,'('' E8    /kcal :'',f11.4)')-e8*autokcal 
      endif
      write(*,*)
      write(*,'('' Edisp /kcal,au:'',f11.4,f12.8)') disp*autokcal,disp
      endif
 
c     if(echo)write(*,*) 'normal termination of dftd3'
c      write(*,*)'---------------------------------------------------'        

      end

