Subroutine sogga11x(Fx,RhoA,RhoB,GradA,GradB,taua,taub) !,  &
!                   ElocX,fa1,fa2,fa3,fa4,fa5,fb1,fb2,fb3,fb4,fb5)
Implicit Real*8(A-H,O-Z)
real*8 CxA(0:5),CxB(0:5)
real*8,parameter :: pi=3.141592653589793D0

!SOGGA11 exchange parameters:
CxA(0) = 0.0d0 ! 0.50000d0 !
CxA(1) = 0.0d0 !-2.95535d0 !
CxA(2) = 0.0d0 ! 15.7974d0 !
CxA(3) = 1.0d0 !-91.1804d0 !
CxA(4) = 0.0d0 ! 96.2030d0 !
CxA(5) = 0.0d0 ! 0.18683d0 !
CxB(0) = 0.0d0 ! 0.50000d0 !
CxB(1) = 0.0d0 ! 3.50743d0 !
CxB(2) = 0.0d0 !-12.9523d0 !
CxB(3) = 0.0d0 ! 49.7870d0 !
CxB(4) = 0.0d0 !-33.2545d0 !
CxB(5) = 0.0d0 !-11.1396d0 !
Cmu = 0.2236536053d0 !factor that multiplies s**2 in gamma
!some pi-dependent prefactors:
As = (24d0*pi**2)**(-1d0/3d0)
Ax = -(3d0/4d0)*(3d0/pi)**(1d0/3d0)

Rho = RhoA+RhoB
GRho = GradA+GradB

If(Rho.gt.1D-6) then
    Rho43 = Rho**(4d0/3d0)
    ElocX = Ax*Rho43
    s = As*Grho/Rho43
    fy = s*s
    PON = Cmu*fy
    Ffrac = 1d0-1d0/(1d0+PON)
    Fexp  = 1d0-exp(-PON)
    fa0 = CxA(0)
    fa1 = CxA(1) *Ffrac
    fa2 = CxA(2) *Ffrac**2
    fa3 = CxA(3) *Ffrac**3
    fa4 = CxA(4) *Ffrac**4
    fa5 = CxA(5) *Ffrac**5
    fb0 = CxB(0)
    fb1 = CxB(1) *Fexp
    fb2 = CxB(2) *Fexp**2
    fb3 = CxB(3) *Fexp**3
    fb4 = CxB(4) *Fexp**4
    fb5 = CxB(5) *Fexp**5
    !enhancement factor:
    Fggax = fa0+fa1+fa2+fa3+fa4+fa5+fb0+fb1+fb2+fb3+fb4+fb5
    !final energy:
    Fx=ElocX*Fggax
!   write(*,"(' ----RPV4----	   :',f20.10,f20.10,f20.10)") Ax,ElocX
end if
end 

Subroutine sogga11c(Fc,RhoA,RhoB,GradA,GradB,taua,taub)
Implicit Real*8(A-H,O-Z)
real*8 CcA(0:5),CcB(0:5)
real*8,parameter :: pi=3.141592653589793D0

! SOGGA11 correlation parameters:       
CcA(0) = 0D+00 ! 0.50000D+00 !1D+00 !
CcA(1) = 0D+00 !-4.62334D+00 !0D+00 !
CcA(2) = 0D+00 ! 8.00410D+00 !0D+00 !
CcA(3) = 1D+00 !-130.226D+00 !0D+00 !
CcA(4) = 0D+00 ! 38.2685D+00 !0D+00 !
CcA(5) = 0D+00 ! 69.5599D+00 !0D+00 !
CcB(0) = 0D+00 ! 0.50000D+00 !0D+00 !
CcB(1) = 0D+00 ! 3.62334D+00 !0D+00 !
CcB(2) = 0D+00 ! 9.36393D+00 !0D+00 !
CcB(3) = 0D+00 ! 34.5114D+00 !0D+00 !
CcB(4) = 0D+00 !-18.5684D+00 !0D+00 !
CcB(5) = 0D+00 !-0.16519D+00 !0D+00 !
XNu = 15.75592D0
CC0 = 0.004235D0
beta= XNu*CC0   
! some pi-dependent prefactors:
AsC=(3d0/(4d0*pi))**(1d0/3d0)
SKFac=2d0*((((3d0*pi*pi)**(1d0/3d0))/pi))**(1d0/2d0)

Rho = RhoA+RhoB
GRho = GradA+GradB

If(Rho.gt.1D-6) then
    rS = (3d0/(4d0*pi)/Rho)**(1d0/3d0)
    Zcor = (RhoA-RhoB)/Rho
    Call PWLc(rS,Zcor,PotLC)
    Call EvGZet(Zcor,Gcor)
    ElocC = Rho*PotLC
    Rho76 = Rho**(7d0/6d0)
    Tcor  = GRho/(2d0*SKFac*Rho76*Gcor)
    T2 = Tcor*Tcor
    G3 = Gcor**3
    PONC = (G3*beta*T2)/PotLC
    Ffracc = 1d0-1d0/(1d0-PONC)
    Fexpc  = 1d0-exp(PONC)
    fa0c = CcA(0)
    fa1c = CcA(1) *Ffracc
    fa2c = CcA(2) *Ffracc**2
    fa3c = CcA(3) *Ffracc**3
    fa4c = CcA(4) *Ffracc**4
    fa5c = CcA(5) *Ffracc**5
    fb0c = CcB(0)
    fb1c = CcB(1) *Fexpc
    fb2c = CcB(2) *Fexpc**2
    fb3c = CcB(3) *Fexpc**3
    fb4c = CcB(4) *Fexpc**4
    fb5c = CcB(5) *Fexpc**5
    !enhancement factor:
    FggaC = fa0c+fa1c+fa2c+fa3c+fa4c+fa5c+fb0c+fb1c+fb2c+fb3c+fb4c+fb5c
    !final energy:
    Fc = ElocC*FggaC   
!write(*,"(' ----RPV4----	   :',f20.10,f20.10,f20.10)") Gcor,Grho
end if
end

Subroutine EvGZet(Zeta,GZeta)
Implicit Real*8(A-H,O-Z)

Real*8 Nine
!Save Zero, One, Two, Three, Four, Nine, F27
Data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/, Three/3.0d0/, &
     Four/4.0d0/, Nine/9.0d0/, F8/8.0D0/, F27/27.0D0/
Small = 1D-10
GZeta = Zero
OMZ = One - Zeta
OPZ = One + Zeta
OMZ2 = OMZ**2
OPZ2 = OPZ**2
F13 = One / Three
F19 = One / Nine
F427 = Four / F27
If(OMZ.gt.Small) then
  OMZ3 = OMZ ** (-F13)
  GZeta = GZeta + OMZ*OMZ3
  endIf
If(OPZ.gt.Small) then
  OPZ3 = OPZ ** (-F13)
  GZeta = GZeta + OPZ*OPZ3
  endIf
GZeta = GZeta / Two
Return
End

Subroutine EvFZet(S,Zeta,FZeta)
Implicit Real*8(A-H,O-Z)
Real*8 Nine
Data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/, Three/3.0d0/, &
     Four/4.0d0/, Nine/9.0d0/, F8/8.0D0/, F27/27.0D0/
Small = 1d-10
FZeta = -Two
OMZ = One - Zeta
OPZ = One + Zeta
OMZ2 = OMZ**2
OPZ2 = OPZ**2
F13 = One / Three
F43 = Four / Three
F49 = Four / Nine
F827 = F8 / F27
If(OMZ.gt.Small) then
  OMZ3 = OMZ ** F13
  fZeta = fZeta + OMZ*OMZ3
  endIf
If(OPZ.gt.Small) then
  OPZ3 = OPZ ** F13
  fZeta = fZeta + OPZ*OPZ3
  endIf
fZeta = fZeta * S
Return
End

Subroutine EvPWLC(A,A1,B1,B2,B3,B4,RS,V)
Implicit Real*8(A-H,O-Z)
!Save F1, F2, F3, F4
Data F1/1.0d0/, F2/2.0d0/, F3/3.0d0/, F4/4.0d0/

Q0 = -F2*A*(F1+A1*RS)
RS12 = Sqrt(RS)
RS32 = RS*RS12
Q1 = F2*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RS)
Q2 = Log(F1+F1/Q1)
V = Q0*Q2
Return
End

Subroutine PWLc(RS,Zeta,PotLC)
Implicit Real*8(A-H,O-Z)
real*8,parameter :: pi=3.141592653589793D0
Dimension ECLPar(6,3)
!Save F1, F2, F3, F4, F6, F8, F9, F12, F24, F36, GammaB, ECLPar,
!     $  BadFZZ, BadEC13
Data F1/1.0d0/, F2/2.0d0/, F3/3.0d0/, F4/4.0d0/, F6/6.0d0/, &
     F8/8.0d0/, F9/9.0d0/, F12/12.0d0/, F24/24.0d0/, F36/36.0d0/, &
     GammaB/0.5198421D0/, BadFZZ/1.709921D0/,BadEC13/0.01688690D0/, &
     ECLPar/0.03109070D0,0.21370D0, 7.5957D0,3.5876D0,1.6382D0,0.49294D0, &
            0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,0.62517D0, &
            0.01688690D0,0.11125D0,10.3570D0,3.6231D0,0.88026D0,0.49671D0/

Third = F1 / F3
EClP13 = F1 / (F6*Pi*Pi)
FZZI = F9*(F2**Third-F1) / F4
GammaI = F1 / (F2*F2**Third-F2)
Call EvFZet(GammaI,Zeta,FZeta)
Call EvPWLC(ECLPar(1,1),ECLPar(2,1),ECLPar(3,1),ECLPar(4,1),ECLPar(5,1),ECLPar(6,1),RS,EU)
Call EvPWLC(ECLPar(1,2),ECLPar(2,2),ECLPar(3,2),ECLPar(4,2),ECLPar(5,2),ECLPar(6,2),RS,EP)
Call EvPWLC(ECLP13,ECLPar(2,3),ECLPar(3,3),ECLPar(4,3),ECLPar(5,3),ECLPar(6,3),RS,AlphaM)
Z2 = Zeta*Zeta
Z3 = Zeta*Z2
Z4 = Zeta*Z3
GZ = FZeta*Z4
HZ = FZZI*(FZeta-GZ)
PotLC = EU*(F1-GZ) + EP*GZ - AlphaM*HZ
Return
End



Subroutine PBEk(Fk,RhoA,RhoB,GradA,GradB,taua,taub)
Implicit Real*8(A-H,O-Z)
real*8 CkA(0:5),CkB(0:5)
real*8,parameter :: pi=3.141592653589793D0

! Becke's
Fmu = 0.2236536053d0
! PBE2
!Fmu = 0.2942d0
! PBE3
!Fmu = 4.1355d0
! PBE4
!Fmu = 1.7107d0
!
!PBE2
!CkA(0) =  1d0
!CkA(1) =  2.0309d0
!CkA(2) =  0d0
!CkA(3) =  0d0
!CkA(4) =  0d0
!CkA(5) =  0d0
!
!PBE3
!CkA(0) =  1d0
!CkA(1) = -3.7425d0
!CkA(2) =  50.258d0
!CkA(3) =  0d0
!CkA(4) =  0d0
!CkA(5) =  0d0
!
!PBE4
!CkA(0) =  1d0
!CkA(1) = -7.2333d0
!CkA(2) =  61.645d0
!CkA(3) = -93.683d0
!CkA(4) =  0d0
!CkA(5) =  0d0

!fa0 should be kept at 1
CkA(0) =  1d0
CkA(1) =  0d0
CkA(2) =  0d0
CkA(3) =  0d0
CkA(4) =  0d0
CkA(5) =  0d0
CkB(0) =  0d0
CkB(1) =  0d0
CkB(2) =  0d0
CkB(3) =  0d0
CkB(4) =  0d0
CkB(5) =  0d0
!some pi-dependent prefactors:
As = (24d0*pi**2)**(-1d0/3d0)
Ctf = (3d0/10d0)*(3d0*pi*pi)**(2d0/3d0) 

Rho = RhoA+RhoB
GRho = GradA+GradB
If(Rho.gt.1D-6) then
    Rho43 = Rho**(4d0/3d0)
    Rho53 = Rho**(5d0/3d0)
    ElocTF = Ctf*Rho53
    s = As*Grho/Rho43
    fy = s*s
    PON = fmu*fy
    Ffrac = 1d0-1d0/(1d0+PON)
    Fexp  = 1d0-exp(-PON)
    fa0 = CkA(0)
    fa1 = CkA(1) *Ffrac
    fa2 = CkA(2) *Ffrac**2
    fa3 = CkA(3) *Ffrac**3
    fa4 = CkA(4) *Ffrac**4
    fa5 = CkA(5) *Ffrac**5
    fb0 = CkB(0)
    fb1 = CkB(1) *Fexp
    fb2 = CkB(2) *Fexp**2
    fb3 = CkB(3) *Fexp**3
    fb4 = CkB(4) *Fexp**4
    fb5 = CkB(5) *Fexp**5
    !enhancement factor(s):
    fts = fa0+fa1+fa2+fa3+fa4+fa5+fb0+fb1+fb2+fb3+fb4+fb5
    Fggak = fts !- F5o3*fy
    !final energy:
    Fk=ElocTF*Fggak
end if
end

