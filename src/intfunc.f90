! INTFUNC: Integration of functions on a Lebedev grid
! This is effectively the main routine of this program (that's why it's here, rather than in sub)
!
subroutine intfunc(intgv)
use function
use util
implicit real*8 (a-h,o-z)
real*8 intvalxc(5),intvalx(5),intvalc(5),intvalold(5),funcval(radpot*sphpot,5),funcvalx(radpot*sphpot,5),funcvalc(radpot*sphpot,5),beckeweigrid(radpot*sphpot)
real*8 funcnuc(radpot*sphpot,5),intvaln(5),funckin(radpot*sphpot,5),intvalk(5)
real*8 intHF(5),intJ(5),intK(5),valJ(radpot*sphpot,5),valK(radpot*sphpot,5)
real*8 nure, intgv
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)

write(*,"(' Radial points:',i5,'	Angular points:',i5,'	Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
call gen1cintgrid(gridatmorg,iradcut)


!initialize some values:
intvalxc=0d0		!xc energy
intvalold=0d0	!xc energy for printing integration on each center
intvalx=0d0		!x energy
intvalc=0d0		!c energy
intvaln=0d0		!nuclear attraction
intvalk=0d0		!kinetic energy
intgv=0d0
fx = 0d0
fxx =0d0

intJ=0			!Coulomb
intK=0			!Exchange 
intHF=0			!HF energy

!if asking for K, J is free:
if (doK==1) doJ=1

!calculate nuclear repulsion energy
call calcnure(nure)

!starting the integration work
do iatm=1,ncenter
	write(*,"(' Processing center',i6,'(',a2,')	  /',i6)") iatm,a(iatm)%name,ncenter
	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(iatm)%y
	gridatm%z=gridatmorg%z+a(iatm)%z
nthreads=getNThreads()
!$OMP parallel do shared(funcval) private(i,rnowx,rnowy,rnowz) num_threads(nthreads)
	do i=1+iradcut*sphpot,radpot*sphpot
		rnowx=gridatm(i)%x
		rnowy=gridatm(i)%y
		rnowz=gridatm(i)%z
! xc vlaues
! USE LIBXC:
		funcvalx(i,1)=DFTxfunc(rnowx,rnowy,rnowz)
		if (ksnamec .NE. " ") funcvalc(i,1)=DFTcfunc(rnowx,rnowy,rnowz)
! USE EXTERNAL FUNCTIONAL (self coded):
!		funcvalx(i,1)=extDFTxfunc(rnowx,rnowy,rnowz)
!		funcvalc(i,1)=extDFTcfunc(rnowx,rnowy,rnowz)
		

! printout of some DFT-development things. This require some tinkering, but it's left here, just in case...
!        call gendensgradab(rnowx,rnowy,rnowz,rhoa,rhob,tdens,grada,gradb,tgrad,abgrad,ttau,atau,btau)
!        call sogga11x(Fx,rhoa,rhob,grada,gradb,taua,taub) !,  
!
!        Fxx=extDFTxfunc(rnowx,rnowy,rnowz)
!
!
!         funcvalx(i,1)=Fxx

!        if (ABS(Fx-fxx).gt.1d-3) then
!         write(*,"(' ----RPV----	   :',f10.5,f10.5,f10.5,f20.10,f20.10,f20.10)") rnowx,rnowy,rnowz,fx,fxx,ABS(Fx-fxx)
!        end if
        
        !write(*,"(' ----RPV----	   :',f10.5,f10.5,f10.5,f20.10,f20.10)") rnowx,rnowy,rnowz,fx,funcvalx(i,1)




!
		funcval(i,1)=funcvalx(i,1)+funcvalc(i,1)
! KS-DFT kinetic energy:
!		funckin(i,1)=Hamkin(rnowx,rnowy,rnowz)
! OF-DFT kinetic energy:
		funckin(i,1)=extOFDFTkfunc(rnowx,rnowy,rnowz)
! core hamiltonian
		funcnuc(i,1)=calcnuc(rnowx,rnowy,rnowz)
! Accumulating Coulomb and Exchange for HF (this is very time consuming)
		if (doJ==1) then
		  valJ(i,1)=calcJ(rnowx,rnowy,rnowz)
		  valK(i,1)=calcK(rnowx,rnowy,rnowz)
		endif
	end do
!$OMP end parallel do
	call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid)
	do i=1+iradcut*sphpot,radpot*sphpot

! xc integration
		intvalxc=intvalxc+funcval(i,:)*gridatmorg(i)%value*beckeweigrid(i)
		intvalx=intvalx+funcvalx(i,:)*gridatmorg(i)%value*beckeweigrid(i)
		intvalc=intvalc+funcvalc(i,:)*gridatmorg(i)%value*beckeweigrid(i)
! core hamiltonian
		intvalk=intvalk+funckin(i,:)*gridatmorg(i)%value*beckeweigrid(i)
		intvaln=intvaln+funcnuc(i,:)*gridatmorg(i)%value*beckeweigrid(i)
! coulomb & HF if necessary		   
		if (doJ==1) then
		  intJ=intJ+valJ(i,:)*gridatmorg(i)%value*beckeweigrid(i)
		  intK=intK+valK(i,:)*gridatmorg(i)%value*beckeweigrid(i)
		endif
	end do
	write(*,"(' Accumulated value:',f20.10,'  Current center:',f20.10)") intvalxc(1),intvalxc(1)-intvalold(1)
	intvalold=intvalxc
end do !integration on the grid

write(*,*)
write(*,"('---=< integration done, printing xc results >=---')")
if (ksnamec .NE. " ") then
	write(*,"(' Final (',a,' + ',a,') exchange-correlation energy:',f20.10)") trim(ksnamex),trim(ksnamec),intvalxc(1)
	write(*,"(' (',a,') exchange energy:',f20.10)") trim(ksnamex),intvalx(1)
	write(*,"(' (',a,') correlation energy:',f20.10)") trim(ksnamec),intvalc(1)
else 
	write(*,"(' Final (',a,') exchange-correlation energy:',f20.10)") trim(ksnamex),intvalxc(1)
end if 
write(*,*)
write(*,"('---=< printing final results >=---')")
write(*,"(' xc energy (from above):',f20.10)") intvalxc(1)
write(*,"(' Kinetic energy        :',f20.10)") intvalK(1)
write(*,"(' Nuclear attraction    :',f20.10)") intvaln(1)
write(*,"(' Nuclear repulsion     :',f20.10)") nure
if (doJ==1) write(*,"(' Coulomb integral (J)  :',f20.10)") intJ(1)
if (doK==1) then 
	write(*,"(' Exchange integral (K) :',f20.10)") intK(1)
	write(*,"(' HF energy             :',f20.10)") intJ(1)+intK(1)
end if
	write(*,"('                        --------------------')")
if (doJ==1) then 
	intgv=nure+intvaln(1)+intvalK(1)+intJ(1)+intvalxc(1)
	write(*,"('     +--> Total ENERGY :',f20.10' <--+',/)") intgv
else
	write(*,"(' Calculation of Coulomb term is off, cannot compute total energy (this is a feature, not a bug!)')")
end if 
write(*,*)

end subroutine















! subroutine sensan(intgv)
! use function
! use util
! implicit real*8 (a-h,o-z)
! real*8 intvalxc(5),intvals00(5),intvals05(5),intvals10(5),intvalold(5),funcval(radpot*sphpot,5),funcvalx(radpot*sphpot,5),funcvalc(radpot*sphpot,5),funcvals00(radpot*sphpot,5),funcvals05(radpot*sphpot,5),funcvals10(radpot*sphpot,5),beckeweigrid(radpot*sphpot)
! real*8 funcnuc(radpot*sphpot,5),intvaln(5),funckin(radpot*sphpot,5),intvalk(5)
! real*8 intHF(5),intJ(5),intK(5),valJ(radpot*sphpot,5),valK(radpot*sphpot,5)
! real*8 nure, intgv
! type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
! 
! write(*,"(' Radial points:',i5,'	Angular points:',i5,'	Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
! call gen1cintgrid(gridatmorg,iradcut)
! 
! 
! !initialize some values:
! intvalxc=0d0		!xc energy
! intvalold=0d0	!xc energy for printing integration on each center
! intvals00=0d0	
! intvals05=0d0	
! intvals10=0d0
! intvaln=0d0		!nuclear attraction
! intvalk=0d0		!kinetic energy
! intgv=0d0
! fx = 0d0
! fxx =0d0
! 
! intJ=0			!Coulomb
! intK=0			!Exchange 
! intHF=0			!HF energy
! 
! 
! !calculate nuclear repulsion energy
! call calcnure(nure)
! 
! !starting the integration work
! do iatm=1,ncenter
! 	write(*,"(' Processing center',i6,'(',a2,')	  /',i6)") iatm,a(iatm)%name,ncenter
! 	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
! 	gridatm%y=gridatmorg%y+a(iatm)%y
! 	gridatm%z=gridatmorg%z+a(iatm)%z
! nthreads=getNThreads()
! !$OMP parallel do shared(funcval) private(i,rnowx,rnowy,rnowz) num_threads(nthreads)
! 	do i=1+iradcut*sphpot,radpot*sphpot
! 		rnowx=gridatm(i)%x
! 		rnowy=gridatm(i)%y
! 		rnowz=gridatm(i)%z
! ! xc vlaues
! ! USE LIBXC:
! 		funcvalx(i,1)=DFTxfunc(rnowx,rnowy,rnowz)
! 		if (ksnamec .NE. " ") funcvalc(i,1)=DFTcfunc(rnowx,rnowy,rnowz)
! 		funcval(i,1)=funcvalx(i,1)+funcvalc(i,1)
! 		
!         funcvals00(i,1)=Ghoshentro(rnowx,rnowy,wnowz,1)
!         funcvals05(i,1)=2D0/3D0*lagkin(rnowx,rnowy,rnowz,0)/fdens(rnowx,rnowy,rnowz)
!         funcvals10(i,1)=sensfunc(rnowx,rnowy,rnowz,1d0)     
! 		
! ! OF-DFT kinetic energy:
! 		funckin(i,1)=extOFDFTkfunc(rnowx,rnowy,rnowz)
! ! core hamiltonian
! 		funcnuc(i,1)=calcnuc(rnowx,rnowy,rnowz)
! ! Accumulating Coulomb and Exchange for HF (this is very time consuming)
! 		if (doJ==1) then
! 		  valJ(i,1)=calcJ(rnowx,rnowy,rnowz)
! !		  valK(i,1)=calcK(rnowx,rnowy,rnowz)
! 		endif
! 	end do
! !$OMP end parallel do
! 	call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid)
! 	do i=1+iradcut*sphpot,radpot*sphpot
! 
! ! xc integration
! 		intvalxc=intvalxc+funcval(i,:)*gridatmorg(i)%value*beckeweigrid(i)
! 
! 		intvals00=intvals00+funcvals00(i,:)*gridatmorg(i)%value*beckeweigrid(i)
!         intvals05=intvals05+funcvals05(i,:)*gridatmorg(i)%value*beckeweigrid(i)
!         intvals10=intvals10+funcvals10(i,:)*gridatmorg(i)%value*beckeweigrid(i)
!         
!                 
!         
!         
! ! core hamiltonian
! 		intvalk=intvalk+funckin(i,:)*gridatmorg(i)%value*beckeweigrid(i)
! 		intvaln=intvaln+funcnuc(i,:)*gridatmorg(i)%value*beckeweigrid(i)
! ! coulomb & HF if necessary		   
! 		if (doJ==1) then
! 		  intJ=intJ+valJ(i,:)*gridatmorg(i)%value*beckeweigrid(i)
! !		  intK=intK+valK(i,:)*gridatmorg(i)%value*beckeweigrid(i)
! 		endif
! 	end do
! 	write(*,"(' Accumulated value:',f20.10,'  Current center:',f20.10)") intvalxc(1),intvalxc(1)-intvalold(1)
! 	intvalold=intvalxc
! end do !integration on the grid
! 
! write(*,*)
! write(*,"('---=< integration done, printing xc results >=---')")
! if (ksnamec .NE. " ") then
! 	write(*,"(' Final (',a,' + ',a,') exchange-correlation energy:',f20.10)") trim(ksnamex),trim(ksnamec),intvalxc(1)
! 
! 
! 
! 
! 	write(*,"(' u00 energy:',f20.10)") intvals00(1)
! 	write(*,"(' u05 energy:',f20.10)") intvals05(1)
!     write(*,"(' u10 energy:',f20.10)") intvals10(1)
! 
! 
!     write(*,"(' u00SD:',f10.5)") (intvals00(1)-intvalxc(1))/intvalxc(1)*100
!     write(*,"(' u05SD:',f10.5)") (intvals05(1)-intvalxc(1))/intvalxc(1)*100
!     write(*,"(' u10SD:',f10.5)") (intvals10(1)-intvalxc(1))/intvalxc(1)*100    
! 
! !	write(*,"(' (',a,') exchange energy:',f20.10)") trim(ksnamex),intvalx(1)
! !	write(*,"(' (',a,') correlation energy:',f20.10)") trim(ksnamec),intvalc(1)
! else 
! 	write(*,"(' Final (',a,') exchange-correlation energy:',f20.10)") trim(ksnamex),intvalxc(1)
! end if 
! write(*,*)
! write(*,"('---=< printing final results >=---')")
! write(*,"(' xc energy (from above):',f20.10)") intvalxc(1)
! write(*,"(' Kinetic energy        :',f20.10)") intvalK(1)
! write(*,"(' Nuclear attraction    :',f20.10)") intvaln(1)
! write(*,"(' Nuclear repulsion     :',f20.10)") nure
! if (doJ==1) write(*,"(' Coulomb integral (J)  :',f20.10)") intJ(1)
! !if (doK==1) then 
! !	write(*,"(' Exchange integral (K) :',f20.10)") intK(1)
! !	write(*,"(' HF energy             :',f20.10)") intJ(1)+intK(1)
! !end if
! 	write(*,"('                        --------------------')")
! if (doJ==1) then 
! 	intgv=nure+intvaln(1)+intvalK(1)+intJ(1)+intvalxc(1)
! 	write(*,"('     +--> Total ENERGY :',f20.10' <--+',/)") intgv
! else
! 	write(*,"(' Calculation of Coulomb term is off, cannot compute total energy (this is a feature, not a bug!)')")
! end if 
! write(*,*)
! 
! end subroutine





!not used at present
!
subroutine dblint(intgv)
use function
use util
implicit real*8 (a-h,o-z)
real*8 intvalxc(5),intvalx(5),intvalc(5),intvalold(5),funcval(radpot*sphpot,5),funcvalx(radpot*sphpot,5),funcvalc(radpot*sphpot,5),beckeweigrid(radpot*sphpot)
real*8 funcnuc(radpot*sphpot,5),intvaln(5),funckin(radpot*sphpot,5),intvalk(5)
real*8 intHF(5),intJ(5),intK(5),valJ(radpot*sphpot,5),valK(radpot*sphpot,5)
real*8 nure, intgv
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)

write(*,"(' Radial points:',i5,'	Angular points:',i5,'	Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
call gen1cintgrid(gridatmorg,iradcut)


!initialize some values:
intvalxc=0d0		!xc energy
intvalold=0d0	!xc energy for printing integration on each center
intvalx=0d0		!x energy
intvalc=0d0		!c energy
intvaln=0d0		!nuclear attraction
intvalk=0d0		!kinetic energy
intgv=0d0

intJ=0			!Coulomb
intK=0			!Exchange 
intHF=0			!HF energy

!if asking for K, J is free:
if (doK==1) doJ=1

!calculate nuclear repulsion energy
call calcnure(nure)

!starting the integration work
do iatm=1,ncenter
	write(*,"(' Processing center',i6,'(',a2,')	  /',i6)") iatm,a(iatm)%name,ncenter
	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(iatm)%y
	gridatm%z=gridatmorg%z+a(iatm)%z
nthreads=getNThreads()
!$OMP parallel do shared(funcval) private(i,rnowx,rnowy,rnowz) num_threads(nthreads)
	do i=1+iradcut*sphpot,radpot*sphpot
		rnowx=gridatm(i)%x
		rnowy=gridatm(i)%y
		rnowz=gridatm(i)%z
! xc vlaues
		funcval(i,1)=funcvalx(i,1)+funcvalc(i,1)
	end do
!$OMP end parallel do
	call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid)
	do i=1+iradcut*sphpot,radpot*sphpot
! xc integration
		intvalxc=intvalxc+funcval(i,:)*gridatmorg(i)%value*beckeweigrid(i)
	end do
	write(*,"(' Accumulated value:',f20.10,'  Current center:',f20.10)") intvalxc(1),intvalxc(1)-intvalold(1)
	intvalold=intvalxc
end do !integration on the grid

write(*,*)
write(*,"('---=< integration done, printing xc results >=---')")
if (ksnamec .NE. " ") then
	write(*,"(' Final (',a,' + ',a,') exchange-correlation energy:',f20.10)") trim(ksnamex),trim(ksnamec),intvalxc(1)
else 
	write(*,"(' Final (',a,') exchange-correlation energy:',f20.10)") trim(ksnamex),intvalxc(1)
end if 
write(*,*)

end subroutine







! Cuba interface in 3D
integer function integrand(ndim, xx, ncomp, func)
use function
implicit none
integer          :: i
integer          :: ndim         ! number of dimensions of the integral
integer          :: ncomp        ! number of components of the integral
double precision :: xx(ndim)      ! point in ndim space
double precision :: func(ncomp)  ! value of integrand in point x

real*8  :: rnowx, rnowy, rnowz, volume, box, minx, miny, minz
real*8  :: xmax,xmin,ymax,ymin,zmax,zmin,xlength,ylength,zlength
real*8  :: funcX, funcC
real*8  coord(3,ncenter)

do i=1,ncenter
    coord(1,i)=a(i)%x*b2a
    coord(2,i)=a(i)%y*b2a
    coord(3,i)=a(i)%z*b2a
end do

xmax= maxval(coord(1,1:ncenter))
xmin= minval(coord(1,1:ncenter))
ymax= maxval(coord(2,1:ncenter))
ymin= minval(coord(2,1:ncenter))
zmax= maxval(coord(3,1:ncenter))
zmin= minval(coord(3,1:ncenter))

xlength=boxscale+(xmax-xmin)
ylength=boxscale+(ymax-ymin)
zlength=boxscale+(zmax-zmin)
volume=xlength*ylength*zlength    

minx=0.5d0*(xmax+xmin-xlength)
miny=0.5d0*(ymax+ymin-ylength)
minz=0.5d0*(zmax+zmin-zlength)

rnowx=minx+xx(1)*xlength
rnowy=miny+xx(2)*ylength
rnowz=minz+xx(3)*zlength

funcX=volume*DFTxfunc(rnowx,rnowy,rnowz)
funcC=volume*DFTcfunc(rnowx,rnowy,rnowz)

func=funcX+funcC

!    write(*,*) '   ---func          ', minx,miny,minz,xlength,ylength,zlength

integrand = 0
   
return
end function integrand




! Cuba interface in 6D for J
integer function integrandJ(ndim2, xx, ncomp, func)
use function
implicit none
integer          :: i
integer          :: ndim2         ! number of dimensions of the integral
integer          :: ncomp        ! number of components of the integral
double precision :: xx(ndim2)      ! point in ndim space
double precision :: func(ncomp)  ! value of integrand in point x

real*8  :: rnowx1,rnowy1,rnowz1,rnowx2,rnowy2,rnowz2 
real*8  :: volume, box, minx, miny, minz
real*8  :: xmax,xmin,ymax,ymin,zmax,zmin,xlength,ylength,zlength
real*8  coord(3,ncenter)

do i=1,ncenter
    coord(1,i)=a(i)%x*b2a
    coord(2,i)=a(i)%y*b2a
    coord(3,i)=a(i)%z*b2a
end do

xmax= maxval(coord(1,1:ncenter))
xmin= minval(coord(1,1:ncenter))
ymax= maxval(coord(2,1:ncenter))
ymin= minval(coord(2,1:ncenter))
zmax= maxval(coord(3,1:ncenter))
zmin= minval(coord(3,1:ncenter))

xlength=boxscale+(xmax-xmin)
ylength=boxscale+(ymax-ymin)
zlength=boxscale+(zmax-zmin)
volume=xlength*ylength*zlength    

minx=0.5d0*(xmax+xmin-xlength)
miny=0.5d0*(ymax+ymin-ylength)
minz=0.5d0*(zmax+zmin-zlength)

rnowx1=minx+xx(1)*xlength
rnowy1=miny+xx(2)*ylength
rnowz1=minz+xx(3)*zlength

rnowx2=minx+xx(4)*xlength
rnowy2=miny+xx(5)*ylength
rnowz2=minz+xx(6)*zlength


func=volume*volume*funcJ(rnowx1,rnowy1,rnowz1,rnowx2,rnowy2,rnowz2)


!    write(*,*) '   ---func          ', rnowx1,rnowy1,rnowz1, volume, func

integrandJ = 0
   
return
end function integrandJ


! Cuba interface in 6D for K
integer function integrandK(ndim2, xx, ncomp, func)
use function
implicit none
integer          :: i
integer          :: ndim2         ! number of dimensions of the integral
integer          :: ncomp        ! number of components of the integral
double precision :: xx(ndim2)      ! point in ndim space
double precision :: func(ncomp)  ! value of integrand in point x

real*8  :: rnowx1,rnowy1,rnowz1,rnowx2,rnowy2,rnowz2 
real*8  :: volume, box, minx, miny, minz
real*8  :: xmax,xmin,ymax,ymin,zmax,zmin,xlength,ylength,zlength
real*8  coord(3,ncenter)

do i=1,ncenter
    coord(1,i)=a(i)%x*b2a
    coord(2,i)=a(i)%y*b2a
    coord(3,i)=a(i)%z*b2a
end do

xmax= maxval(coord(1,1:ncenter))
xmin= minval(coord(1,1:ncenter))
ymax= maxval(coord(2,1:ncenter))
ymin= minval(coord(2,1:ncenter))
zmax= maxval(coord(3,1:ncenter))
zmin= minval(coord(3,1:ncenter))

xlength=boxscale+(xmax-xmin)
ylength=boxscale+(ymax-ymin)
zlength=boxscale+(zmax-zmin)
volume=xlength*ylength*zlength    

minx=0.5d0*(xmax+xmin-xlength)
miny=0.5d0*(ymax+ymin-ylength)
minz=0.5d0*(zmax+zmin-zlength)

rnowx1=minx+xx(1)*xlength
rnowy1=miny+xx(2)*ylength
rnowz1=minz+xx(3)*zlength

rnowx2=minx+xx(4)*xlength
rnowy2=miny+xx(5)*ylength
rnowz2=minz+xx(6)*zlength


func=volume*volume*funcK(rnowx1,rnowy1,rnowz1,rnowx2,rnowy2,rnowz2)


!    write(*,*) '   ---func          ', rnowx1,rnowy1,rnowz1, volume, func

integrandK = 0
   
return
end function integrandK