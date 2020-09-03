!!---------- Return normalization coefficient for specific type of cartesian type GTF, see Levine 5ed p487
!The meaning of itype is defined in GTFtype2name
real*8 function normgau(itype,exp)
use defvar
use util
implicit real*8 (a-h,o-z)
ix=type2ix(itype)
iy=type2iy(itype)
iz=type2iz(itype)
normgau=(2*exp/pi)**0.75D0*dsqrt( (8*exp)**(ix+iy+iz)*ft(ix)*ft(iy)*ft(iz)/(ft(2*ix)*ft(2*iy)*ft(2*iz)) )
end function

real*8 function rnormgau_ORCA(exp,Lval)
real*8 exp,pi
integer Lval,n1,n2,nf
pi=acos(-1D0)
call genn1n2nf(Lval,n1,n2,nf)
rnormgau_ORCA=dsqrt(dsqrt(2**n1*exp**n2/(pi**3*nf*nf)))
end function

real*8 function doSint(iGTF,jGTF)
integer iGTF,jGTF
real*8,external :: doSintactual
doSint=doSintactual(iGTF,jGTF,0,0,0,0,0,0)
end function

!!!------------------ Evaluate overlap integral for two unnormalized GTFs
!~p arguments are the shifts of GTF index, used by doKint, doveloint but not by doSint
real*8 function doSintactual(iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p)
use util
use defvar
implicit real*8(a-h,o-z)
integer iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p
x1=a(b(iGTF)%center)%x
y1=a(b(iGTF)%center)%y
z1=a(b(iGTF)%center)%z
x2=a(b(jGTF)%center)%x
y2=a(b(jGTF)%center)%y
z2=a(b(jGTF)%center)%z
ee1=b(iGTF)%exp
ee2=b(jGTF)%exp
ep=ee1+ee2
sqrtep=dsqrt(ep)
px=(ee1*x1+ee2*x2)/ep
py=(ee1*y1+ee2*y2)/ep
pz=(ee1*z1+ee2*z2)/ep		 
expterm=dexp( -ee1*ee2*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)/ep )
ix1=type2ix(b(iGTF)%functype)+ix1p
iy1=type2iy(b(iGTF)%functype)+iy1p
iz1=type2iz(b(iGTF)%functype)+iz1p
ix2=type2ix(b(jGTF)%functype)+ix2p
iy2=type2iy(b(jGTF)%functype)+iy2p
iz2=type2iz(b(jGTF)%functype)+iz2p
!chen book,P103
numx=ceiling( (ix1+ix2+1)/2D0 ) !Need to calculate n points
sx=0.0D0
do i=1,numx
	tmp=Rhm(numx,i)/sqrtep+px
	term1=(tmp-x1)**ix1
	term2=(tmp-x2)**ix2
	sx=sx+Whm(numx,i)*term1*term2
end do
sx=sx/sqrtep

numy=ceiling( (iy1+iy2+1)/2D0 )
sy=0.0D0
do i=1,numy
	tmp=Rhm(numy,i)/sqrtep+py
	term1=(tmp-y1)**iy1
	term2=(tmp-y2)**iy2
	sy=sy+Whm(numy,i)*term1*term2
end do
sy=sy/sqrtep

numz=ceiling( (iz1+iz2+1)/2D0 )
sz=0.0D0
do i=1,numz
	tmp=Rhm(numz,i)/sqrtep+pz
	term1=(tmp-z1)**iz1
	term2=(tmp-z2)**iz2
	sz=sz+Whm(numz,i)*term1*term2
end do
sz=sz/sqrtep

doSintactual=sx*sy*sz*expterm
end function

!!!------------------ Generate overlap matrix between all GTFs
!nsize should be nprims*(nprims+1)/2
!------ A warpper of dodipoleint, used to directly get a single component of dipole moment integral. icomp=1/2/3 corresponds to X/Y/Z component
real*8 function dipintcomp(icomp,iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p)
integer icomp,iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p
real*8 xcomp,ycomp,zcomp
call dodipoleint(iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p,xcomp,ycomp,zcomp)
if (icomp==1) dipintcomp=xcomp
if (icomp==2) dipintcomp=ycomp
if (icomp==3) dipintcomp=zcomp
end function

!!!------------------- Evaluate kinetic integral (i.e. -(1/2)der2 )for two unnormalized GTFs, see Chen's book, p104
real*8 function doTint(iGTF,jGTF)
use defvar
implicit real*8(a-h,o-z)
integer iGTF,jGTF
real*8 term(4)
ee1=b(iGTF)%exp
ee2=b(jGTF)%exp
ix1=type2ix(b(iGTF)%functype)
iy1=type2iy(b(iGTF)%functype)
iz1=type2iz(b(iGTF)%functype)
ix2=type2ix(b(jGTF)%functype)
iy2=type2iy(b(jGTF)%functype)
iz2=type2iz(b(jGTF)%functype)
term=0
if(ix1>0.and.ix2>0)	 term(1)=	ix1*ix2*doSintactual(iGTF,jGTF,-1,0,0,-1,0,0)
if(ix1>0)			 term(2)=-2*ee2*ix1*doSintactual(iGTF,jGTF,-1,0,0, 1,0,0)
if(ix2>0)			 term(3)=-2*ee1*ix2*doSintactual(iGTF,jGTF, 1,0,0,-1,0,0)
					 term(4)= 4*ee1*ee2*doSintactual(iGTF,jGTF, 1,0,0, 1,0,0)
Tx=sum(term)
term=0
if(iy1>0.and.iy2>0)	 term(1)=	iy1*iy2*doSintactual(iGTF,jGTF,0,-1,0,0,-1,0)
if(iy1>0)			 term(2)=-2*ee2*iy1*doSintactual(iGTF,jGTF,0,-1,0,0, 1,0)
if(iy2>0)			 term(3)=-2*ee1*iy2*doSintactual(iGTF,jGTF,0, 1,0,0,-1,0)
					 term(4)= 4*ee1*ee2*doSintactual(iGTF,jGTF,0, 1,0,0, 1,0)
Ty=sum(term)
term=0
if(iz1>0.and.iz2>0)	 term(1)=	iz1*iz2*doSintactual(iGTF,jGTF,0,0,-1,0,0,-1)
if(iz1>0)			 term(2)=-2*ee2*iz1*doSintactual(iGTF,jGTF,0,0,-1,0,0, 1)
if(iz2>0)			 term(3)=-2*ee1*iz2*doSintactual(iGTF,jGTF,0,0, 1,0,0,-1)
					 term(4)= 4*ee1*ee2*doSintactual(iGTF,jGTF,0,0, 1,0,0, 1)
Tz=sum(term)
doTint=(Tx+Ty+Tz)/2
end function

