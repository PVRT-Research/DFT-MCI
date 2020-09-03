module deftype
implicit none
type atomtype
character*2 name !name of atom
integer index !The index in periodic table, if ECP was used, charge will smaller than this value
real*8 x,y,z,charge !Coordinate(Bohr) and charge of atoms
end type

type primtype
integer center,functype !The number of nuclei that the basis function centered on and its function type
real*8 exp !Exponent
end type

type content !Type for grid data points
real*8 x,y,z,value
end type

end module


module defvar
use deftype
implicit none
integer :: ido
!============ Store globally shared information
!
! Relevant for DFT:IoG (usually these are read from a file):
character*20 :: ksnamex="xc_gga_x_b88",ksnamec="xc_gga_c_lyp", density=""
integer :: doINT=1,doJ=1,doJMC=0,doK=1,doGrimme=0,doDFTc=0 !doINT =0 no integration, =1 do grid integration, =2 do MC integration, =3 do both
integer :: iautointgrid=1,radpot=72,sphpot=302 !sphpot=230/302/434/590/770, coarse is 72*302, ultrafine is 99*590
!
! Variables for MC integration
integer :: MClvl=4,seed=0,flags=0,key1=10,key2=1,key3=1,maxpass=5,ldxgiven=0,nnew=1000,key=0 !(MClvl selects the MC algorithm: 1=vegas, 2=divonne, 3=suave, 4=cuhre)
integer :: nstart=1000,nincrease=500,nbatch=1000,gridno=0
integer*8 ::  ngiven=0,nextra=0,mineval=96,maxeval=10000000 !(*this one should go in the input file)
double precision :: epsrel=1.d-16, epsabs=1.d-03 !(*these determine the precision of the integration in hartrees, must go in input)
double precision :: userdata=0.d0,border=0.d0,maxchisq=1.d0,mindeviation=.25d0, flatness=25.d0, boxscale=10.d0 !(scale of integration box must go in input)
!
!-------- Variables for wfn information(_org means using for backuping the first loaded molecule)
integer :: ibasmode=0 !0/1 = GTO/STO is used in current wavefunction
integer :: nmo=0,nprims=0,ncenter=0,ncenter_org=0,nmo_org=0,nprims_org=0 !Number of orbitals, primitive functions, nuclei
integer :: idxHOMO=0 !For fch and molden, record the index of original HOMO, this will be used to calculate linear response kernel for pi-electrons
integer :: ifiletype=0 !plain text=0, fch=1, wfn=2, wfx=3, chg=4, pdb/xyz=5, NBO .31=6, cub=7, grd=8, molden=9, gms=10, MDL mol=11
integer :: wfntype=0 !0/1/2/3/4 means R/U/ROHF /R/U-Post-HF wavefunction
real*8 :: totenergy=0,virialratio=2,nelec=0,naelec=0,nbelec=0
!-------- Variables for nuclei & basis function & Molecular orbital
type(atomtype),allocatable :: a(:),a_org(:),a_tmp(:)
type(primtype),allocatable :: b(:),b_org(:),b_tmp(:)
integer,allocatable :: connmat(:,:) !Connectivity matrix
real*8,allocatable :: MOocc(:),MOocc_org(:),MOene(:) !Occupation number & energy of orbital
integer,allocatable :: MOtype(:) !The type of orbitals, (alpha&beta)=0/alpha=1/beta=2, not read from .wfn directly
character*4,allocatable :: MOsym(:) !The symmetry of orbitals, meaningful when .molden/.gms is used since it sometimes records irreducible representation
real*8,allocatable :: CO(:,:),CO_org(:,:),CO_tmp(:,:) !Coefficient matrix of primitive basis functions, including both normalization and contraction coefficients
!
!
real*8,parameter :: pi=3.141592653589793D0,b2a=0.529177249D0 !1 Bohr = 0.529177249 Angstrom
real*8,parameter :: au2kcal=627.51D0,au2KJ=2625.5D0,au2eV=27.2113838D0
real*8,parameter :: masse=9.10938215D-31,chge=1.602176487D0,lightc=2.99792458D8,au2debye=2.5417462D0 !masse/chge: Mass/charge of an electron
real*8,parameter :: planckc=6.62606896D-34,h_bar=1.054571628D-34,amu2kg=1.66053878D-27
real*8,parameter :: boltzc=1.3806488D-23,boltzcau=3.1668114D-6,boltzceV=8.6173324D-5 !in J/K, in Hartree/K and in eV/K, respectively
real*8,parameter :: avogacst=6.02214179D23
integer,parameter :: nelesupp=150 !The number of elements supported, ghost(index=0) is not taken into account
real*8 ctrval(1000) !Value of contour lines

!Store important calculated data
real*8,allocatable :: cubmat(:,:,:) !cubmat, store density/laplacian...3D-matrix
real*8,allocatable :: cubmattmp(:,:,:) !For cube data exchanging, position of points must be identical to cubmat, so don't use type(content) for lowering memory consumption
real*8,allocatable :: cubmatvec(:,:,:,:) !Used to store vector field
real*8,allocatable :: gradd1(:,:),gradd2(:,:) !Gradient in direction1/2 for gradient line plot
real*8,allocatable :: distmat(:,:) !Distance matrix, in Bohr
character*200 filename,firstfilename
character*80 firstfileonlyname,extctrsetting
character,allocatable :: custommapname(:)*200,customop(:) !Custom operation for custom map/cube file
logical alive
integer,allocatable :: fragatm(:),fragatmbackup(:) !Store the index of atoms in fragment, has no relationship with frag1/frag2. fragatmbackup is used to backup fragatm during custom operation
integer,allocatable :: frag1(:),frag2(:) !These two fragments are only used for bond order analysis/composition analysis etc., store index of basis functions or atoms. Their size just fit their content
integer :: ncustommap=0,imodwfn=0 !if 1, means occupation number or orbital type or basis function information have been modified
integer :: iorbsel=1 !Which orbital is selected, and its value will be calculated by fmo and calchessmat_mo
integer :: iorbsel2=0 !Which orbital will be plotted together with iorbsel in plane map


!The name for superheavy atoms are consistent with Stuttgart PP website: http://www.tc.uni-koeln.de/PP/clickpse.en.html
character*2 :: ind2name(0:nelesupp)=(/ "Bq","H ","He", &   !X(number O) is ghost atom
"Li","Be","B ","C ","N ","O ","F ","Ne", & !3~10
"Na","Mg","Al","Si","P ","S ","Cl","Ar", & !11~18
"K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", & !19~36
"Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe", & !37~54
"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", & !55~71
"Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", & !72~86
"Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr", & !87~103
"Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Ut","Fl","Up","Lv","Us","Uo","Un","Ux",("??",ido=121,nelesupp) /) !104~all  Such as Uuo/Uup is replaced by Uo/Up
character*2 :: ind2name_up(0:nelesupp)=(/ "BQ","H ","HE", & !Same as ind2name, but all characters are upper case, to cater to .pdb file
"LI","BE","B ","C ","N ","O ","F ","NE", & !3~10
"NA","MG","AL","SI","P ","S ","CL","AR", & !11~18
"K ","CA","SC","TI","V ","CR","MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR","KR", & !19~36
"RB","SR","Y ","ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN","SN","SB","TE","I ","XE", & !37~54
"CS","BA","LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB","LU", & !55~71
"HF","TA","W ","RE","OS","IR","PT","AU","HG","TL","PB","BI","PO","AT","RN", & !72~86
"FR","RA","AC","TH","PA","U ","NP","PU","AM","CM","BK","CF","ES","FM","MD","NO","LR", & !87~103
"RF","DB","SG","BH","HS","MT","DS","RG","CN","UT","FL","UP","LV","US","UO","UN","UX",("??",ido=121,nelesupp) /) !104~all
!Bondi vdW radius, from J.Phys.Chem.,1964,68(3),441-451, unit is Angstrom, will be convert to Bohr when multiwfn start
real*8 :: vdwr(0:nelesupp)=(/ 0.4D0,1.2D0,1.4D0,& !Ghost,H,He
1.82D0,1.77D0,1.74D0,1.7D0,1.55D0,1.52D0,1.47D0,1.54D0,& !Li~Ne
2.27D0,1.73D0,1.73D0,2.1D0,1.8D0,1.8D0,1.75D0,1.88D0,& !Na~Ar
(2.0D0,ido=19,27),1.63D0,1.4D0,1.39D0,1.87D0,2.0D0,1.85D0,1.9D0,1.85D0,2.02D0,& ! Ni~Kr(28~36)
(2.0D0,ido=37,45),1.63D0,1.72D0,1.58D0,1.93D0,2.17D0,2.0D0,2.06D0,1.98D0,2.16D0,& !Pd~Xe(46~54)
(2.0D0,ido=55,77),1.72D0,1.66D0,1.55D0,1.96D0,2.02D0,(2.0D0,ido=83,nelesupp) /) !Pt~Pb(78~82)
!##No use currently!## Modified Bondi vdW radii, but for all main group (except for H,He), use IVA radius in corresponding row. Specifically used to molecular surface decomposition
real*8 :: vdwr_tianlu(0:nelesupp)=(/ 0.4D0,1.7D0,1.7D0,& !Ghost,H,Ne   H and Ne are augmented to carbon radius
(1.7D0,ido=3,10),& !Li~Ne
(2.1D0,ido=11,18),& !Na~Ar
1.87D0,1.87D0,	(2.0D0,ido=21,27),1.63D0,1.40D0,1.39D0,	 (1.87D0,ido=31,36),& !K ,Ca,  Ni~Zn(21~30),  Ga~Kr(31,37)
1.93D0,1.93D0,	(2.0D0,ido=39,45),1.63D0,1.72D0,1.58D0,	 (1.93D0,ido=49,54),& !Rb,Sr,  Y ~Cd(39~48),  In~Xe(49~54)
1.96D0,1.96D0,	(2.0D0,ido=57,77),1.72D0,1.66D0,1.55D0,	 (1.96D0,ido=81,86),& !Cs,Ba,  La~Hg(57~80),  Tl~Rn(81~86)
(2.0D0,ido=87,nelesupp) /) !Rn~Mt(87~109,~all)

!Covalent radius, from "Dalton Trans., 2008, 2832-2838", unit is Angstrom, will be convert to Bohr when multiwfn start
real*8 :: covr(0:nelesupp)=(/ 0.1D0,0.31D0,0.28D0,& !Ghost,H,Ne(1~2)
1.28D0,0.96D0,0.84D0,0.76D0,0.71D0,0.66D0,0.57D0,0.58D0,& !Li~Ne(3~10)	   here C is sp3
1.66D0,1.41D0,1.21D0,1.11D0,1.07D0,1.05D0,1.02D0,1.06D0,& !Na~Ar(11~18)
2.03D0,1.76D0,1.70D0,1.60D0,1.53D0,1.39D0,1.39D0,1.32D0,1.26D0,& !K~Co(19~27)  here MnD0,FeD0,Co is low-spinD0, high spin is 1.61D0,1.52D0,1.50
1.24D0,1.32D0,1.22D0,1.22D0,1.20D0,1.19D0,1.20D0,1.20D0,1.16D0,& !Ni~Kr(28~36)
2.20D0,1.95D0,1.90D0,1.75D0,1.64D0,1.54D0,1.47D0,1.46D0,1.42D0,& !Rb~Rh(37~45)
1.39D0,1.45D0,1.44D0,1.42D0,1.39D0,1.39D0,1.38D0,1.39D0,1.40D0,& !Pd~Xe(46~54)
2.44D0,2.15D0,2.07D0,2.04D0,2.03D0,2.01D0,1.99D0,1.98D0,1.98D0,1.96D0,1.94D0,& !Cs~Tb(55~65)
1.92D0,1.92D0,1.89D0,1.90D0,1.87D0,1.87D0,1.75D0,1.70D0,1.62D0,1.51D0,1.44D0,1.41D0,& !Dy~Ir(66~77)
1.36D0,1.36D0,1.32D0,1.45D0,1.46D0,1.48D0,1.40D0,1.50D0,1.50D0,2.60D0,2.21D0,& !Pt~Ra(78~88)
2.15D0,2.06D0,2.00D0,1.96D0,1.90D0,1.87D0,1.80D0,1.69D0,(1.5D0,ido=97,nelesupp) /) !Ac~Cm(89~96),~all
!(Covalent) radius proposed by Suresh, from J. Phys. Chem. A 2001, 105, 5940-5944. For missing values (including all noble gases and very heavy elements), the ones in covr array are used
!Unit is Angstrom, will be convert to Bohr when multiwfn start
real*8 :: covr_Suresh(0:nelesupp)=(/ 0.1D0,0.327D0,0.28D0,& !Ghost,H,Ne(1~2)
1.219D0,0.911D0,0.793D0,0.766D0,0.699D0,0.658D0,0.633D0,0.58D0,& !Li~Ne(3~10)
1.545D0,1.333D0,1.199D0,1.123D0,1.11D0,1.071D0,1.039D0,1.06D0,& !Na~Ar(11~18)
1.978D0,1.745D0,1.337D0,1.274D0,1.236D0,1.128D0,1.18D0,1.091D0,1.089D0,& !K~Co(19~27)
1.077D0,1.146D0,1.187D0,1.199D0,1.179D0,1.209D0,1.201D0,1.201D0,1.16D0,& !Ni~Kr(28~36)
2.217D0,1.928D0,1.482D0,1.377D0,1.353D0,1.24D0,1.287D0,1.212D0,1.229D0,& !Rb~Rh(37~45)
1.24D0,1.362D0,1.429D0,1.385D0,1.38D0,1.421D0,1.4D0,1.397D0,1.40D0,& !Pd~Xe(46~54)
2.442D0,2.149D0,1.653D0,2.04D0,2.03D0,2.01D0,1.99D0,1.98D0,1.98D0,1.96D0,1.94D0,& !Cs~Tb(55~65)
1.92D0,1.92D0,1.89D0,1.90D0,1.87D0,1.87D0,1.364D0,1.346D0,1.256D0,1.258D0,1.222D0,1.227D0,& !Dy~Ir(66~77)
1.227D0,1.273D0,1.465D0,1.531D0,1.434D0,1.496D0,1.40D0,1.50D0,1.50D0,2.60D0,2.21D0,& !Pt~Ra(78~88)
2.15D0,2.06D0,2.00D0,1.96D0,1.90D0,1.87D0,1.80D0,1.69D0,(1.5D0,ido=97,nelesupp) /) !Ac~Cm(89~96),~all
!Covalent radius, from Pyykko "Chem. Eur. J.,15,186 (2009)", unit is in Angstrom, will be convert to Bohr when multiwfn start
real*8 :: covr_pyy(0:nelesupp)=(/ 0.1D0,0.32D0,0.46D0,& !Ghost,H,Ne(1~2)
1.33D0,1.02D0,0.85D0,0.75D0,0.71D0,0.63D0,0.64D0,0.67D0,& !Li~Ne(3~10)
1.55D0,1.39D0,1.26D0,1.16D0,1.11D0,1.03D0,0.99D0,0.96D0,& !Na~Ar(11~18)
1.96D0,1.71D0,1.48D0,1.36D0,1.34D0,1.22D0,1.19D0,1.16D0,1.11D0,& !K~Co(19~27)
1.10D0,1.12D0,1.18D0,1.24D0,1.21D0,1.21D0,1.16D0,1.14D0,1.17D0,& !Ni~Kr(28~36)
2.10D0,1.85D0,1.63D0,1.54D0,1.47D0,1.38D0,1.28D0,1.25D0,1.25D0,& !Rb~Rh(37~45)
1.20D0,1.28D0,1.36D0,1.42D0,1.40D0,1.40D0,1.36D0,1.33D0,1.31D0,& !Pd~Xe(46~54)
2.32D0,1.96D0,1.80D0,1.63D0,1.76D0,1.74D0,1.73D0,1.72D0,1.68D0,1.69D0,1.68D0,& !Cs~Tb(55~65)
1.67D0,1.66D0,1.65D0,1.64D0,1.70D0,1.62D0,1.52D0,1.46D0,1.37D0,1.31D0,1.29D0,1.22D0,& !Dy~Ir(66~77)
1.23D0,1.24D0,1.34D0,1.44D0,1.44D0,1.51D0,1.45D0,1.47D0,1.42D0,2.23D0,2.01D0,& !Pt~Ra(78~88)
1.86D0,1.75D0,1.69D0,1.70D0,1.71D0,1.72D0,1.66D0,1.66D0,1.68D0,1.68D0,1.65D0,1.67D0,1.73D0,1.76D0,1.61D0,& !Ac~Lr(89~103)
1.57D0,1.49D0,1.43D0,1.41D0,1.34D0,1.29D0,1.28D0,1.21D0,1.22D0,1.36D0,1.43D0,1.62D0,1.75D0,1.65D0,1.57D0,(1.5D0,ido=119,nelesupp)  /) !Rf~118(104~118),~all
real*8 :: covr_tianlu(0:nelesupp)=(/ 0.1D0,0.31D0,0.28D0,& !H,Ne(1~2) !Based on CSD radii, but for all main group (except for H,He), use IVA radius in corresponding row
(0.76D0,ido=3,10),& !Li~Ne(3~10)
(1.11D0,ido=11,18),1.2D0,1.2D0,& !Na~Ar(11~18),K,Ca
1.70D0,1.60D0,1.53D0,1.39D0,1.39D0,1.32D0,1.26D0,1.24D0,1.32D0,1.22D0,& !Sc~Zn(21~30)  here MnD0,FeD0,Co is low-spinD0, high spin is 1.61D0,1.52D0,1.50
(1.2D0,ido=31,36),1.42D0,1.42D0,& !Ga~Kr(31~36),Rb,Sr
1.90D0,1.75D0,1.64D0,1.54D0,1.47D0,1.46D0,1.42D0,1.39D0,1.45D0,1.44D0,& !Y~Cd(39~48)
(1.39D0,ido=49,54),1.46D0,1.46D0,& !In~Xe(49~54),Cs,Ba
2.07D0,2.04D0,2.03D0,2.01D0,1.99D0,1.98D0,1.98D0,1.96D0,1.94D0,1.92D0,1.92D0,1.89D0,1.90D0,1.87D0,1.87D0,& !La~Lu(57~71)
1.75D0,1.70D0,1.62D0,1.51D0,1.44D0,1.41D0,1.36D0,1.36D0,1.32D0,& !Hf~Hg(72~80)
(1.46D0,ido=81,86),1.46D0,1.46D0,&!Tl~Rn(81~86),Fr(still 1.46),Ra(still 1.46)
2.15D0,2.06D0,2.00D0,1.96D0,1.90D0,1.87D0,1.80D0,1.69D0,(1.5D0,ido=97,nelesupp) /) !Ac~Cm(89~96),~all
!Radii proposed in Chem. Phys. Lett., 480 (2009) 127-131, the unit is Bohr!
real*8 :: radii_Hugo(0:nelesupp)=(/ 0.10D0,1.00D0,0.74D0,& !Ghost,H,Ne(1~2)
1.59D0,1.21D0,1.28D0,1.10D0,0.97D0,1.00D0,0.88D0,0.79D0,& !Li~Ne(3~10)
1.63D0,1.33D0,1.51D0,1.29D0,1.14D0,1.15D0,1.02D0,0.93D0,& !Na~Ar(11~18)
1.77D0,1.49D0,1.44D0,1.41D0,1.42D0,1.42D0,1.35D0,1.31D0,1.31D0,1.33D0,1.33D0,1.20D0,1.51D0,1.31D0,1.18D0,1.18D0,1.07D0,0.99D0,& !K~Kr
1.80D0,1.55D0,1.48D0,1.43D0,1.42D0,1.38D0,1.37D0,1.36D0,1.35D0,1.28D0,1.34D0,1.23D0,1.53D0,1.36D0,1.26D0,1.23D0,1.14D0,1.06D0,1.87D0,1.62D0,& !Rb~Xe,Cs,Ba
1.56D0,1.57D0,1.58D0,1.57D0,1.56D0,1.55D0,1.55D0,1.49D0,1.52D0,1.51D0,1.50D0,1.49D0,1.48D0,1.47D0,1.58D0,& !La~Lu
1.41D0,1.34D0,1.31D0,1.32D0,1.27D0,1.23D0,1.23D0,1.21D0,1.14D0,1.49D0,1.35D0,1.37D0,1.27D0,1.21D0,1.12D0,1.83D0,1.16D0,& !Hf~Rn,Fr,Ra
1.62D0,1.47D0,1.52D0,1.48D0,1.47D0,1.50D0,1.51D0,1.51D0,1.48D0,1.47D0,1.46D0,1.45D0,1.44D0,1.43D0,1.67D0,1.51D0,(1.5D0,ido=105,nelesupp) /) !Ac~Rf,~all

real*8 :: YWTatomcoeff(18,3)=reshape((/ & !Coef. of fitting B3LYP/6-31G* density by Weitao Yang group for the first three rows, see supporting info. of JACS,132,6498
0.2815D0,2.437D0,11.84D0,31.34D0,67.82D0,120.2D0,190.9D0,289.5D0,406.3D0,561.3D0,760.8D0,1016.0D0,1319.0D0,1658.0D0,2042.0D0,2501.0D0,3024.0D0,3625.0D0, &
0.0D0,0.0D0,0.06332D0,0.3694D0,0.8527D0,1.172D0,2.247D0,2.879D0,3.049D0,6.984D0,22.42D0,37.17D0,57.95D0,87.16D0,115.7D0,158.0D0,205.5D0,260.0D0, &
0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.06358D0,0.3331D0,0.8878D0,0.7888D0,1.465D0,2.17D0,3.369D0,5.211D0 /),(/18,3/))
real*8 :: YWTatomexp(18,3)=reshape((/ & !Corresponding exponent of YWTatom, the value setted to 1.0 don't have any meaning, only for avoiding divide zero
0.5288D0,0.3379D0,0.1912D0,0.139D0,0.1059D0,0.0884D0,0.0767D0,0.0669D0,0.0608D0,0.0549D0,0.0496D0,0.0449D0,0.0411D0,0.0382D0,0.0358D0,0.0335D0,0.0315D0,0.0296D0, &
1.0D0,1.0D0,0.9992D0,0.6945D0,0.53D0,0.548D0,0.4532D0,0.3974D0,0.3994D0,0.3447D0,0.2511D0,0.215D0,0.1874D0,0.1654D0,0.1509D0,0.1369D0,0.1259D0,0.1168D0, &
1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0236D0,0.7753D0,0.5962D0,0.6995D0,0.5851D0,0.5149D0,0.4974D0,0.4412D0 /),(/18,3/))
!Atomic weights, from http://www.chem.qmul.ac.uk/iupac/AtWt/, the data is mainly based on Pure Appl. Chem., 81, 2131-2156 (2009)
real*8 :: atmwei(0:nelesupp)=(/ 0D0,1.00794D0,4.0026D0,6.941D0,9.01218D0,10.811D0,12.0107D0,14.0067D0,15.9994D0,18.9984D0,20.1797D0,& !1~10
22.98977D0,24.305D0,26.98154D0,28.0855D0,30.97376D0,32.065D0,35.453D0,39.948D0,39.0983D0,40.078D0,& !11~20
44.95591D0,47.867D0,50.9415D0,51.9961D0,54.93805D0,55.845D0,58.93319D0,58.6934D0,63.546D0,65.38D0,& !21~30
69.723D0,72.64D0,74.9216D0,78.96D0,79.904D0,83.798D0,85.4678D0,87.62D0,88.90585D0,91.224D0,& !31~40
92.90638D0,95.96D0,98D0,101.07D0,102.9055D0,106.42D0,107.8682D0,112.411D0,114.818D0,118.71D0,& !41~50
121.76D0,127.6D0,126.90447D0,131.293D0,132.90545D0,137.327D0,138.90547D0,140.116D0,140.90765D0,144.242D0,& !51~60
145D0,150.36D0,151.964D0,157.25D0,158.92535D0,162.5D0,164.93032D0,167.259D0,168.93421D0,173.054D0,& !61~70
174.9668D0,178.49D0,180.94788D0,183.84D0,186.207D0,190.23D0,192.217D0,195.084D0,196.96657D0,200.59D0,& !71~80
204.3833D0,207.2D0,208.9804D0,209D0,210D0,222D0,223D0,226D0,227D0,232.03806D0,& !71~90
231.03588D0,238.02891D0,237D0,244D0,243D0,247D0,247D0,251D0,252D0,257D0,258D0,259D0,262D0,265D0,268D0,271D0,272D0,270D0,276D0,& !91~109
281D0,282D0,285D0,285D0,289D0,289D0,293D0,294D0,294D0,& !110~118
(0D0,ido=119,nelesupp) /) !119~all
 
!Series of Lebedev-Laikov routines
integer :: Lebelist(32)=(/ 6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810 /)
integer :: fact(0:10)=(/ 1,1,2,6,24,120,720,5040,40320,362880,3628800 /) ! Store factorials from 0~10 
integer :: isphergau=0 !By default, all basis functions are cartesian type, =1 means spherical (but some of them can be still cartesian type)
character*5 :: GTFtype2name(-32:56)=(/ & !The definition of such as G-4, H+5 can be found in http://sobereva.com/97
"H 0  ","H+1  ","H-1  ","H+2  ","H-2  ","H+3  ","H-3  ","H+4  ","H-4  ","H+5  ","H-5  ", & !-32:-22
"G 0  ","G+1  ","G-1  ","G+2  ","G-2  ","G+3  ","G-3  ","G+4  ","G-4  ", & !-21:-13
"F 0  ","F+1  ","F-1  ","F+2  ","F-2  ","F+3  ","F-3  ","D 0  ","D+1  ","D-1  ","D+2  ","D-2  ", & !-12:-6,-5:-1
"     ","S    ","X    ","Y    ","Z    ","XX   ","YY   ","ZZ   ","XY   ","XZ   ","YZ   ", & !0~10
"XXX  ","YYY  ","ZZZ  ","XXY  ","XXZ  ","YYZ  ","XYY  ","XZZ  ","YZZ  ","XYZ  ", & !f 11~20
"ZZZZ ","YZZZ ","YYZZ ","YYYZ ","YYYY ","XZZZ ","XYZZ ","XYYZ ","XYYY ","XXZZ ","XXYZ ","XXYY ","XXXZ ","XXXY ","XXXX ", & !g 21~35
"ZZZZZ","YZZZZ","YYZZZ","YYYZZ","YYYYZ","YYYYY","XZZZZ","XYZZZ","XYYZZ","XYYYZ","XYYYY","XXZZZ","XXYZZ","XXYYZ","XXYYY","XXXZZ","XXXYZ","XXXYY","XXXXZ","XXXXY","XXXXX" /) !h 36~56
character*5 :: type2ang(56)=(/ &
"S    ","P    ","P    ","P    ","D    ","D    ","D    ","D    ","D    ","D    ", & !0~10
"F    ","F    ","F    ","F    ","F    ","F    ","F    ","F    ","F    ","F    ", & !f 11~20
"G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ", & !g 21~35
"H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    " /) !h 36~56
!Here s,p,d sequences are identical to .wfn, .wfx, .fch, .molden  !Note: Sequence in .fch = sequence in Gaussian
!Here f sequence is identical to .wfn, .wfx, but not identical to .fch and .molden
!Here g sequence is identical to .fch, .wfn does not support higher than f function, not identical to .wfx and .molden
!here h sequence is identical to .wfx and .fch, .molden doesn't support h
!Notice: The .wfn produced by G09 B.01 and later supports g and h, the definition is identical to here, and thus can be normally loaded
!Overall, spd: Multiwfn=wfn=wfx=fch=molden	 f: Multiwfn=wfn=wfx!=fch=molden   g: Multiwfn=fch!=wfx=molden=Molden2AIM	h: Multiwfn=wfx=fch
integer :: type2ix(56)=(/ 0, 1,0,0, 2,0,0,1,1,0, 3,0,0,2,2,0,1,1,0,1, 0,0,0,0,0,1,1,1,1,2,2,2,3,3,4, 0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5 /)
integer :: type2iy(56)=(/ 0, 0,1,0, 0,2,0,1,0,1, 0,3,0,1,0,2,2,0,1,1, 0,1,2,3,4,0,1,2,3,0,1,2,0,1,0, 0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0 /)
integer :: type2iz(56)=(/ 0, 0,0,1, 0,0,2,0,1,1, 0,0,3,0,1,1,0,2,2,1, 4,3,2,1,0,3,2,1,0,2,1,0,1,0,0, 5,4,3,2,1,0,4,3,2,1,0,3,2,1,0,2,1,0,1,0,0 /)
!Negative value means the shell use spherical gauss function. -1=SP (also known as L in GAMESS), and impossible be used in Multiwfn (when detect it, split it as S and P)
character :: shtype2name(-5:5)=(/ "H","G","F","D","L","S","P","D","F","G","H" /)
!Convert shell type to the number of basis functions in the shell: 0=s,1=p,-1=sp,2=6d,-2=5d,3=10f,-3=7f,4=15g,-4=9g,5=21h,-5=11h
integer :: shtype2nbas(-5:5)=(/ 11,9,7,5,4,1,3,6,10,15,21 /) 

														!Note: Row/column of CO denote MO/GTF respectively, in contrary to convention
!-------- Describe inner electron density in EDF section
type(primtype),allocatable :: b_EDF(:)
real*8,allocatable :: CO_EDF(:)
integer :: nEDFprims=0,nEDFelec=0 !Electrons represented by EDF
!-------- Variables when basis functions are basis rather than primitive function as basis
integer :: nbasis=0,nshell=0,nprimshell=0 !The number of basis, basis shell and primitive shell. SP shell is counted as S and P shell separately
integer,allocatable :: shtype(:),shcon(:),shcen(:) !Type, contraction degree and attributed center of a basis shell
real*8,allocatable :: primshexp(:),primshcoeff(:) !Exponent and contraction coefficient of a primitive shell
integer,allocatable :: basshell(:) !The ith element is the shell index that the ith basis attributed to
integer,allocatable :: bascen(:),bastype(:) !Center/type of basis, definition is the same as GTF
integer,allocatable :: basstart(:),basend(:) !The ith element means the basis from where to where is attributed to the ith atom
integer,allocatable :: primstart(:),primend(:)	!The ith element means the GTF from where to where is attributed to the ith basis function
real*8,allocatable :: primconnorm(:) !element i means the contract. coeff. * normalization coeff. of GTF i, can be used for e.g. constructing basis integral from GTF integral
real*8,allocatable :: Sbas(:,:),Sbas_org(:,:) !Overlap matrix and its backup
real*8,allocatable :: Dbas(:,:,:) !Dipole moment integral matrix, the first index 1,2,3=X,Y,Z
real*8,allocatable :: Tbas(:,:) !Kinetic energy integral matrix
real*8,allocatable :: Vbas(:,:) !Nuclear attraction potential integral matrix
real*8,allocatable :: Velbas(:,:,:) !Velocity integral matrix, the first index 1,2,3=X,Y,Z
real*8,allocatable :: Magbas(:,:,:) !Magnetic integral matrix, the first index 1,2,3=X,Y,Z
!Coefficient matrix for alpha/beta orbital, CObasa(i,j) means the coefficient of ith basis in the jth orbital, differ to CO(:,:)
real*8,allocatable,target :: CObasa(:,:),CObasb(:,:) !wfntype==0,2,3 only allocate CObasa(nbasis,nmo), ==1,4 also allocate CObasb, dimension of both cobasa and cobasb would be (nbasis,nbasis)
real*8,allocatable,target :: CObasa_org(:,:),CObasb_org(:,:)
real*8,allocatable,target :: Ptot(:,:),Palpha(:,:),Pbeta(:,:) !Density matrix of total/alpha/beta, for wfntype==0.or.wfntype==3, only Ptot is filled, for others, all of Ptot,Palpha and Pbeta are filled
real*8,allocatable :: Palpha_org(:,:),Pbeta_org(:,:) !Backup P, e.g. for Wiberg bond order calculation
! real*8,allocatable :: twoPDM(:) !Store two-particle density matrix by one-dimension array, Not use currently

!-------- Trajectory
integer :: nframetraj=0 !The number of frames in the trajectory
real*8,allocatable :: traj(:,:,:) !traj(1/2/3,a,i) corresponds to x/y/z of the ath atom in frame i
!-------- Points loaded from external file
integer :: numextpt=0
real*8,allocatable :: extpt(:,:),extpttmp(:) !extpt(i,1:4) corresponds to X/Y/Z/value of point i, length unit is bohr. extpttmp only records function value
!-------- Atomic radial densities, originally loaded from .rad file
integer,allocatable :: atmradnpt(:) !How many radial points that each atom has
real*8 atmradpos(200) !Position of radial points, shared by all atoms, since it is generated by the same rule
real*8,allocatable :: atmradrho(:,:) !(iatm,j) corresponds to density value at j radial point of iatm atom
!---------- Used for passing data of spectrum plotting, as well as used by DOS
real*8,allocatable :: datax(:),str(:),FWHM(:) !Transition energy, strength and FWHM loaded from only one file


!!!!!!!!!!!!!!!!!!!!!! Parameter !!!!!!!!!!!!!!!!!!!!!!
!For passing Dislin main parent GUI and draw widget identifier
integer idissetlight1,idissetlight2,idissetlight3,idissetlight4,idissetlight5,idissetlightall0,idissetlightall1,idissetangle,idissetplaneXVU,idissetplaneYVU
integer idisgraph,idiszoomin,idiszoomout,idisisosurscl,idisscrval,idisshowbothsign,idisshowisosur,idisshowdatarange,idisshowmol,idisisosursec,iorbseltext,iorbtxt,iorblis
integer idisisosurquality,idisisosurnumpt,idisorbinfo2
integer idisshowatmlab,idisshowaxis,idisbondradius,idislabelsize,idisbondcrit,idisatmsize,idisshowpathlab !In draw mol GUI
integer idisshowattlab,idisdrawinternalbasin,idisattsize !Draw basin GUI
integer idisshow3n3,idisshow3n1,idisshow3p1,idisshow3p3,idisshowCPlab,idisshowpath,idisshowbassurf
integer idisshowlocminlab,idisshowlocmaxlab,idisshowlocminpos,idisshowlocmaxpos !For molecular surface analysis
!For setting isosurface style, colors
integer idisisosur1style,idisisosur1solid,idisisosur1mesh,idisisosur1point,idisisosur1solidmesh,idisisosur1tpr,idisisosur1opa
integer idisisosur2style,idisisosur2solid,idisisosur2mesh,idisisosur2point,idisisosur2solidmesh,idisisosur2tpr,idisisosur2opa
integer idisisosurallstyle,idisisosurallsolid,idisisosurallmesh,idisisosurallpoint,idisisosurallsolidmesh,idisisosuralltpr
integer GUI_mode !=1: Show mol and orbitals =2: Show plane =3: Show isosurface =4: Show mol and CPs =5: Show mol and surface extrema =6: Show basin or domain space

!Plotting external parameter, can be set in settings.ini
character :: graphformat*4="png " !ps/eps/pdf/wmf/gif/tiff/bmp
integer :: graph1Dwidth=1280,graph1Dheight=800,graph2Dwidth=1280,graph2Dheight=1200,graph3Dwidth=1400,graph3Dheight=1400
integer :: itickreverse=0,iticks=2,symbolsize=8,ilenunit1D=1,ilenunit2D=1,iatmlabtype=1,iatmlabtype3D=3,iplaneextdata=0
integer :: numdigx=2,numdigy=2,numdigz=3,numdiglinex=3,numdigliney=3,numdigctr=3
real*8 :: planestpx=1.5D0,planestpy=1.5D0,planestpz=0.1D0
integer :: fillcoloritpx=5,fillcoloritpy=3,pleatmlabsize=50
real*8 :: disshowlabel=0.5D0
real*8 :: bondclrR=0.1D0,bondclrG=1.0D0,bondclrB=0.1D0,atmlabclrR=0D0,atmlabclrG=0D0,atmlabclrB=0D0
real*8 :: CP3n3RGB(3)=(/0.72D0,0D0,0.72D0/),CP3n1RGB(3)=(/1D0,0.5D0,0D0/),CP3p1RGB(3)=(/1D0,1D0,0D0/),CP3p3RGB(3)=(/0D0,1D0,0D0/)
real*8 :: atm3Dclr(0:nelesupp,3) !Colors of the atom spheres shown in 3D plots, set in "loadsetting" routine

!Plotting Internal parameter
integer :: imodlayout=0
integer :: idrawbasinidx=0,idrawinternalbasin=0 !Draw which basin. If draw interal part of the basin
integer :: ifixorbsign=0 !if 1, during generating orbital isosurface by drawmolgui, most part will always be positive (namely if sum(cubmat)<0 or sum(cubmattmp)<0, the data sign will be inverted)
integer :: iatom_on_contour,iatom_on_contour_far=0,plesel,IGRAD_ARROW=0,ILABEL_ON_CONTOUR,LASTCTRVAL,ictrlabsize=20,ivdwctrlabsize=0,iwidthvdwctr=10,iwidthposctr=1,iwidthnegctr=1,iwidthgradline=1
integer :: iclrindctrpos=5,iclrindctrneg=5,ivdwclrindctr=3,iclrindgradline=6,vdwctrstyle(2)=(/ 1,0 /),ctrposstyle(2)=(/ 1,0 /),ctrnegstyle(2)=(/ 10,15 /)
integer :: isavepic=0,icurve_vertlinex=0,iclrindatmlab=1,imarkrefpos=0,ilog10y=0,iclrcurve=1,inowhiteblack=0
integer :: ifragcontri=0,nfragatmnum,nfragatmnumbackup,ipromol,inucespplot=0,idrawmol=1,idrawisosur=0,isosursec=0,idrawtype=1,idrawcontour=1
integer :: iinvgradvec=0,icolorvecfield=0,vecclrind=30,idrawplanevdwctr=0,iplaneoutall=0,icurvethick=5
real*8 :: surcolorzmin,surcolorzmax !fillctr is the contour value will be draw on fillcolor map
real*8 :: curve_vertlinex=0D0,curvexyratio=0.618D0 !Gold partition
real*8 :: gradplotstep=0.002D0,gradplotdis=0.01D0,gradplottest=0.2D0,cutgradvec=0.3D0
real*8 :: clrRcub1same=0.3D0,clrGcub1same=0.75D0,clrBcub1same=0.3D0,clrRcub1oppo=0.3D0,clrGcub1oppo=0.45D0,clrBcub1oppo=0.9D0 !Color for isosurface 1 with solid style
real*8 :: clrRcub2same=0.4D0,clrGcub2same=0.5D0,clrBcub2same=0.0D0,clrRcub2oppo=0.35D0,clrGcub2oppo=0.1D0,clrBcub2oppo=0.9D0 !Color for isosurface 2 with solid style
real*8 :: clrRcub1samemeshpt=0.3D0,clrGcub1samemeshpt=0.75D0,clrBcub1samemeshpt=0.3D0,clrRcub1oppomeshpt=0.3D0,clrGcub1oppomeshpt=0.45D0,clrBcub1oppomeshpt=0.9D0 !Color for isosurface 1 with solid style
real*8 :: clrRcub2samemeshpt=0.4D0,clrGcub2samemeshpt=0.5D0,clrBcub2samemeshpt=0.0D0,clrRcub2oppomeshpt=0.35D0,clrGcub2oppomeshpt=0.1D0,clrBcub2oppomeshpt=0.9D0 !Color for isosurface 2 with solid style
real*8 :: opacitycub1=0.7D0,opacitycub2=0.7D0 !Opacity for isosurface 1 and 2 with transparent style
!About topology information on plane
integer :: imark3n3=1,imark3n1=1,imark3p1=1,imark3p3=1,imarkpath=1,sizemarkcp=30,sizemarkpath=5,sizemark3n1path=5,idrawintbasple=0,isurfstyle=2
real*8 :: clrRpath=0.3D0,clrGpath=0.1D0,clrBpath=0D0,clrR3n1path=0.0D0,clrG3n1path=0D0,clrB3n1path=0.5D0
integer,allocatable :: boldlinelist(:)
character*3 :: drawsurmesh="ON "
!Parameter for drawing molecular structure or 3D map
integer :: ienablelight1=1,ienablelight2=1,ienablelight3=1,ienablelight4=0,ienablelight5=0 !If enable lighting 1~5
integer :: ishowatmlab=1,ishowCPlab=0,ishowpathlab=0,ishowaxis=1,isosurshowboth=1,ishowdatarange=0,idrawpath=1,idrawbassurf=1,ishowattlab=0,ishowatt=0
integer :: isosur1style=1,isosur2style=1 !isosurface style,1/2/3/4/5=solid,mesh,points,solid+mesh,transparent
integer :: ishowlocminlab=0,ishowlocmaxlab=0,ishowlocminpos=0,ishowlocmaxpos=0 !For molecular surface analysis
integer :: ishow3n3=0,ishow3n1=0,ishow3p1=0,ishow3p3=0
real*8 :: bondcrit=1.15D0,textheigh=38.0D0,ratioatmsphere=1.0D0,ratioCPsphere=0.8D0,bondradius=0.2D0,attsphsize=0.1D0
real*8 :: XVU=150.0D0,YVU=30.0D0,ZVU=6.0D0 !3D view angle
!Parameter for drawing domain defined by isosurfaces as grids
integer :: idrawdomainidx=0,idrawdomain=0,ndomain=0
real*8,allocatable :: gridxyz(:,:) !XYZ coordinate of grid that statisfied criterion
integer,allocatable :: domainsize(:) !The number of grids contained in each domain
integer,allocatable :: domaingrid(:,:) !The grid indices contained in each domain
!For passing ploting parameter from GUI routine to their call-back routine
!sur_value: The value of isosurface will be plot by drawmol routine when idrawisosur=1
real*8 :: dp_init1,dp_end1,dp_init2,dp_end2,dp_init3,dp_end3,sur_value=0.05D0

!!! Other external parameter !!!
integer :: ispecial=0 !=0: Normal, =1 specific for Chunying Rong, =2 for Shubin's 2nd project
integer :: isys=2 !isys=1/2/3: Windows/Linux/MacOS
integer :: igenDbas=0,igenMagbas=0,igenP=1,iwfntmptype=1,outmedinfo=0,intmolcust=0,isilent=0,idelvirorb=1,ifchprog=1,iloadascart=0
!integer :: iuserfunc=1000,iDFTxcsel=84,ispheratm=1,ishowchgtrans=0,SpherIVgroup=0,MCvolmethod=2,readEDF=1,isupplyEDF=2,ishowptESP=1,imolsurparmode=1
integer :: NICSnptlim=8000
real*8 :: bndordthres=0.05D0,compthres=0.5D0,compthresCDA=1D0,expcutoff=-40D0,espprecutoff=0D0
integer :: nthreads,ompstacksize=100000000
integer :: iniNThreads=0
character :: lastfile*200="",gaupath*80=""
!About function calculation, external or internal parameters
integer :: RDG_addminimal=1,ELF_addminimal=1,num1Dpoints=3000,atomdenscut=1,nprevorbgrid=120000
integer :: ipolarpara=0,iALIEdecomp=0,iskipnuc=0
real*8 :: RDG_maxrho=0.05D0,RDGprodens_maxrho=0.1D0,aug1D=1.5D0,aug2D=4.5D0,aug3D=6.0D0,radcut=10.0D0
real*8 :: refx=0D0,refy=0D0,refz=0D0
real*8 :: pleA=0D0,pleB=0D0,pleC=0D0,pleD=0D0 !!ABCD of the plane defined by main function 1000, used for special aims
real*8 :: globaltmp=0 !A variable can be used anywhere and can be set by option 5 of main function 1000, for debugging purpose avoiding re-compile code
!About line/plane/grid calculation, inner parameter
!For 3D grid data
real*8 :: orgx,orgy,orgz,endx,endy,endz,dx,dy,dz !Origin, end point and translation length
integer :: nx=80,ny=80,nz=80 !The number of points in three directions
!For 2D plane map
real*8 :: v1x,v1y,v2x,v2y,v1z,v2z,a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,d1,d2 !Translation vector 1 and 2, three point in self-defined plane for projecting label, d1,d2=Length of v1,v2
real*8 :: orgx2D,orgy2D,orgz2D !Origin
integer :: ngridnum1=100,ngridnum2=100 !The number of points in two directions
!Specific for Shubin's project
real*8 :: steric_addminimal=1D-4,steric_potcutrho=0D0,steric_potcons=0D0
!Other
integer :: ifirstMultiwfn=1 !If 1, means we re-load file via main function -11 and don't need to do some initializations

!Used for EDR(r;d) and D(r)
integer,parameter :: max_edr_exponents=50 !Maximum EDR exponents used to calculate EDR(r;d) and D(r). 
real*8 :: dedr,edrastart,edrainc	!Length scale to define EDR(r;d), start and increment in exponents to evaluate D(r) 
real*8 :: wrtexpo(max_edr_exponents) !For users write the number of EDR exponents that will be used in calculation. 
integer nedr	!No of EDR exponents used to evaluate D(r). 

contains
  integer function getNThreads()
	IMPLICIT NONE
	INTEGER currNThreads, TID, OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
!$OMP PARALLEL PRIVATE(TID) SHARED(currNThreads)
	  TID = OMP_GET_THREAD_NUM()
	  !Only master thread does this
!	   PRINT *, 'Hello from thread ', TID
	  IF (TID .EQ. 0) THEN
		currNThreads = OMP_GET_MAX_THREADS()
!		 PRINT *, 'OMP_NUM_THREADS = ', currNThreads
	  END IF
!$OMP END PARALLEL
	IF (iniNThreads .NE. 0) THEN
	 currNThreads = iniNThreads
	END IF
!	 PRINT *, 'Number of threads = ', getNThreads
	getNThreads = currNThreads
  END FUNCTION
end module


module util
implicit real*8(a-h,o-z)

interface sort
	module procedure sortr8
	module procedure sorti4
end interface
interface invarr
	module procedure invarrr8
	module procedure invarri4
end interface
!!------------------- Root and weight of hermite polynomial
real*8 Rhm(10,10),Whm(10,10)
data Rhm(1,1) /	 0.0D0						/
data Rhm(2,1) / -0.70710678118654752440D+00 /
data Rhm(2,2) /	 0.70710678118654752440D+00 /
data Rhm(3,1) / -1.22474487139158904910D+00 /
data Rhm(3,2) /	 0.0D0						/
data Rhm(3,3) /	 1.22474487139158904910D+00 /
data Rhm(4,1) / -1.65068012388578455588D+00 /
data Rhm(4,2) / -0.52464762327529031788D+00 /
data Rhm(4,3) /	 0.52464762327529031788D+00 /
data Rhm(4,4) /	 1.65068012388578455588D+00 /
data Rhm(5,1) / -2.02018287045608563293D+00 /
data Rhm(5,2) / -0.95857246461381850711D+00 /
data Rhm(5,3) /	 0.0D0						/
data Rhm(5,4) /	 0.95857246461381850711D+00 /
data Rhm(5,5) /	 2.02018287045608563293D+00 /
data Rhm(6,1) / -2.35060497367449222283D+00 /
data Rhm(6,2) / -1.33584907401369694971D+00 /
data Rhm(6,3) / -0.43607741192761650868D+00 /
data Rhm(6,4) /	 0.43607741192761650868D+00 /
data Rhm(6,5) /	 1.33584907401369694971D+00 /
data Rhm(6,6) /	 2.35060497367449222283D+00 /
data Rhm(7,1) / -2.65196135683523349245D+00 /
data Rhm(7,2) / -1.67355162876747144503D+00 /
data Rhm(7,3) / -0.81628788285896466304D+00 /
data Rhm(7,4) /	 0.0D0						/
data Rhm(7,5) /	 0.81628788285896466304D+00 /
data Rhm(7,6) /	 1.67355162876747144503D+00 /
data Rhm(7,7) /	 2.65196135683523349245D+00 /
data Rhm(8,1) / -2.93063742025724401922D+00 /
data Rhm(8,2) / -1.98165675669584292585D+00 /
data Rhm(8,3) / -1.15719371244678019472D+00 /
data Rhm(8,4) / -0.38118699020732211685D+00 /
data Rhm(8,5) /	 0.38118699020732211685D+00 /
data Rhm(8,6) /	 1.15719371244678019472D+00 /
data Rhm(8,7) /	 1.98165675669584292585D+00 /
data Rhm(8,8) /	 2.93063742025724401922D+00 /
data Rhm(9,1) / -3.19099320178152760723D+00 /
data Rhm(9,2) / -2.26658058453184311180D+00 /
data Rhm(9,3) / -1.46855328921666793167D+00 /
data Rhm(9,4) / -0.72355101875283757332D+00 /
data Rhm(9,5) /	 0.0D0						/
data Rhm(9,6) /	 0.72355101875283757332D+00 /
data Rhm(9,7) /	 1.46855328921666793167D+00 /
data Rhm(9,8) /	 2.26658058453184311180D+00 /
data Rhm(9,9) /	 3.19099320178152760723D+00 /
data Rhm(10,1) /  -3.43615911883773760333D+00 /
data Rhm(10,2) /  -2.53273167423278979641D+00 /
data Rhm(10,3) /  -1.75668364929988177345D+00 /
data Rhm(10,4) /  -1.03661082978951365418D+00 /
data Rhm(10,5) /  -0.34290132722370460879D+00 /
data Rhm(10,6) /   0.34290132722370460879D+00 /
data Rhm(10,7) /   1.03661082978951365418D+00 /
data Rhm(10,8) /   1.75668364929988177345D+00 /
data Rhm(10,9) /   2.53273167423278979641D+00 /
data Rhm(10,10) /  3.43615911883773760333D+00 /
data Whm(1,1) / 1.77245385090551602730D+00 / ! SQRT(PI)
data Whm(2,1) / 8.86226925452758013649D-01 /
data Whm(2,2) / 8.86226925452758013649D-01 /
data Whm(3,1) / 2.95408975150919337883D-01 /
data Whm(3,2) / 1.18163590060367735153D+00 /
data Whm(3,3) / 2.95408975150919337883D-01 /
data Whm(4,1) / 8.13128354472451771430D-02 /
data Whm(4,2) / 8.04914090005512836506D-01 /
data Whm(4,3) / 8.04914090005512836506D-01 /
data Whm(4,4) / 8.13128354472451771430D-02 /
data Whm(5,1) / 1.99532420590459132077D-02 /
data Whm(5,2) / 3.93619323152241159828D-01 /
data Whm(5,3) / 9.45308720482941881226D-01 /
data Whm(5,4) / 3.93619323152241159828D-01 /
data Whm(5,5) / 1.99532420590459132077D-02 /
data Whm(6,1) / 4.53000990550884564086D-03 /
data Whm(6,2) / 1.57067320322856643916D-01 /
data Whm(6,3) / 7.24629595224392524092D-01 /
data Whm(6,4) / 7.24629595224392524092D-01 /
data Whm(6,5) / 1.57067320322856643916D-01 /
data Whm(6,6) / 4.53000990550884564086D-03 /
data Whm(7,1) / 9.71781245099519154149D-04 /
data Whm(7,2) / 5.45155828191270305922D-02 /
data Whm(7,3) / 4.25607252610127800520D-01 /
data Whm(7,4) / 8.10264617556807326765D-01 /
data Whm(7,5) / 4.25607252610127800520D-01 /
data Whm(7,6) / 5.45155828191270305922D-02 /
data Whm(7,7) / 9.71781245099519154149D-04 /
data Whm(8,1) / 1.99604072211367619206D-04 /
data Whm(8,2) / 1.70779830074134754562D-02 /
data Whm(8,3) / 2.07802325814891879543D-01 /
data Whm(8,4) / 6.61147012558241291030D-01 /
data Whm(8,5) / 6.61147012558241291030D-01 /
data Whm(8,6) / 2.07802325814891879543D-01 /
data Whm(8,7) / 1.70779830074134754562D-02 /
data Whm(8,8) / 1.99604072211367619206D-04 /
data Whm(9,1) / 3.96069772632643819046D-05 /
data Whm(9,2) / 4.94362427553694721722D-03 /
data Whm(9,3) / 8.84745273943765732880D-02 /
data Whm(9,4) / 4.32651559002555750200D-01 /
data Whm(9,5) / 7.20235215606050957124D-01 /
data Whm(9,6) / 4.32651559002555750200D-01 /
data Whm(9,7) / 8.84745273943765732880D-02 /
data Whm(9,8) / 4.94362427553694721722D-03 /
data Whm(9,9) / 3.96069772632643819046D-05 /
data Whm(10,1) /  7.64043285523262062916D-06 /
data Whm(10,2) /  1.34364574678123269220D-03 /
data Whm(10,3) /  3.38743944554810631362D-02 /
data Whm(10,4) /  2.40138611082314686417D-01 /
data Whm(10,5) /  6.10862633735325798784D-01 /
data Whm(10,6) /  6.10862633735325798784D-01 /
data Whm(10,7) /  2.40138611082314686417D-01 /
data Whm(10,8) /  3.38743944554810631362D-02 /
data Whm(10,9) /  1.34364574678123269220D-03 /
data Whm(10,10) / 7.64043285523262062916D-06 /

contains
!Content sequences:
!!Geometry operation
!!Array, Vector
!!String process
!!Matrix calculation
!!Misc

!===============================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Geometry operation !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================!


!!---------- Get angle (degree) between two vectors
real*8 function vecang(vec1x,vec1y,vec1z,vec2x,vec2y,vec2z)
real*8 vec1x,vec1y,vec1z,vec2x,vec2y,vec2z
pi=3.141592653589793D0
rnorm1=dsqrt(vec1x**2+vec1y**2+vec1z**2)
rnorm2=dsqrt(vec2x**2+vec2y**2+vec2z**2)
costheta=(vec1x*vec2x+vec1y*vec2y+vec1z*vec2z)/rnorm1/rnorm2
if (costheta>1D0) costheta=1
vecang=acos(costheta)/pi*180
end function

!!---------- Get distance of point 0 to plane(defined by point 1,2,3) 
real*8 function potpledis(x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0)
real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,prjx,prjy,prjz
call pointprjple(x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,prjx,prjy,prjz)
potpledis=dsqrt((x0-prjx)**2+(y0-prjy)**2+(z0-prjz)**2)
end function

!!---------- Project a point(x0,y0,z0) to a plane defined by x/y/z-1/2/3, prjx/y/z are results
subroutine pointprjple(x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,prjx,prjy,prjz)
real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,prjx,prjy,prjz,A,B,C,D,t
call pointABCD(x1,y1,z1,x2,y2,z2,x3,y3,z3,A,B,C,D)
! (x0-x)/A=(y0-y)/B=(z0-z)/C ---> x=x0-t*A y=y0-t*B z=z0-t*C , substitute into Ax+By+Cz+D=0 solve t
t=(D+A*x0+B*y0+C*z0)/(A**2+B**2+C**2)
prjx=x0-t*A
prjy=y0-t*B
prjz=z0-t*C
end subroutine

!!-------- Input three points, get ABCD of Ax+By+Cz+D=0
subroutine pointABCD(x1,y1,z1,x2,y2,z2,x3,y3,z3,A,B,C,D)
real*8 v1x,v1y,v1z,v2x,v2y,v2z,x1,y1,z1,x2,y2,z2,x3,y3,z3,A,B,C,D
v1x=x2-x1
v1y=y2-y1
v1z=z2-z1
v2x=x3-x1
v2y=y3-y1
v2z=z3-z1
! Solve determinant(Vector multiply) to get the normal vector (A,B,C):
!  i   j   k   //unit vector
! v1x v1y v1z
! v2x v2y v2z
A=v1y*v2z-v1z*v2y
B=-(v1x*v2z-v1z*v2x)
C=v1x*v2y-v1y*v2x
D=A*(-x1)+B*(-y1)+C*(-z1)
end subroutine

!!-------- Input three points 0,1,2, get the vertical projection point of 0 to the line linking 1-2
subroutine pointprjline(x0,y0,z0,x1,y1,z1,x2,y2,z2,prjx,prjy,prjz)
real*8 x0,y0,z0,x1,y1,z1,x2,y2,z2,prjx,prjy,prjz,v12x,v12y,v12z,t
v12x=x2-x1
v12y=y2-y1
v12z=z2-z1
!Since prjx=x1+t*v12x prjy=y1+t*v12y prjz=z1+t*v12z
!So v12x*(x0-prjx)+v12y*(y0-prjy)+v12z*(z0-prjz)=0
!v12x*(x0-x1)-v12x*v12x*t + v12y*(y0-y1)-v12y*v12y*t + v12z*(z0-z1)-v12z*v12z*t =0
t=( v12x*(x0-x1)+v12y*(y0-y1)+v12z*(z0-z1) )/(v12x**2+v12y**2+v12z**2)
prjx=x1+t*v12x
prjy=y1+t*v12y
prjz=z1+t*v12z
end subroutine

!!---------- Get distance of point 0 to the line 1-2
real*8 function potlinedis(x0,y0,z0,x1,y1,z1,x2,y2,z2)
real*8 x0,y0,z0,x1,y1,z1,x2,y2,z2,prjx,prjy,prjz
call pointprjline(x0,y0,z0,x1,y1,z1,x2,y2,z2,prjx,prjy,prjz)
potlinedis=dsqrt((x0-prjx)**2+(y0-prjy)**2+(z0-prjz)**2)
end function

!!--------- Input three points, return angle between 1-2 and 2-3 (in degree)
real*8 function xyz2angle(x1,y1,z1,x2,y2,z2,x3,y3,z3)
real*8 :: x1,y1,z1,x2,y2,z2,x3,y3,z3,pi=3.141592653589793D0
vec1x=x1-x2
vec1y=y1-y2
vec1z=z1-z2
vec2x=x3-x2
vec2y=y3-y2
vec2z=z3-z2
dotprod=vec1x*vec2x+vec1y*vec2y+vec1z*vec2z
rnormv1=dsqrt( vec1x**2+vec1y**2+vec1z**2 )
rnormv2=dsqrt( vec2x**2+vec2y**2+vec2z**2 )
xyz2angle=acos(dotprod/(rnormv1*rnormv2))/pi*180
end function

!!--------- Input four points, return dihedral angle (in degree)
real*8 function xyz2dih(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
real*8 :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,pi=3.141592653589793D0,phi
v12x=x1-x2
v12y=y1-y2
v12z=z1-z2
v23x=x2-x3
v23y=y2-y3
v23z=z2-z3
v34x=x3-x4
v34y=y3-y4
v34z=z3-z4
call vecprod(v12x,v12y,v12z,v23x,v23y,v23z,p1x,p1y,p1z)
call vecprod(v23x,v23y,v23z,v34x,v34y,v34z,p2x,p2y,p2z)
phi=acos( (p1x*p2x+p1y*p2y+p1z*p2z)/(sqrt(p1x*p1x+p1y*p1y+p1z*p1z)*sqrt(p2x*p2x+p2y*p2y+p2z*p2z)) )
xyz2dih=phi/pi*180
end function

!!--------------- Get area of a triangle, need input coordinates of three points
function gettriangarea(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz)
implicit real*8 (a-h,o-z)
real*8 gettriangarea,pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz
! a---b				va=pb-pa vb=pc-pa
! |
! V
! c
va1=pbx-pax
va2=pby-pay
va3=pbz-paz
vb1=pcx-pax
vb2=pcy-pay
vb3=pcz-paz
call vecprod(va1,va2,va3,vb1,vb2,vb3,vc1,vc2,vc3)  !vc=va¡Ávb=|va||vb|sin¦È*i  where i is unit vector perpendicular to va and vb
absvc=dsqrt(vc1**2+vc2**2+vc3**2)
gettriangarea=0.5D0*absvc
end function

!!--------------- Get volume of a tetrahedron, need input coordinates of four points
function gettetravol(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,pdx,pdy,pdz)
implicit real*8 (a-h,o-z)
real*8 pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,pdx,pdy,pdz
! real*8 volmat(4,4)
! volmat(:,1)=(/ pax,pbx,pcx,pdx /)
! volmat(:,2)=(/ pay,pby,pcy,pdy /)
! volmat(:,3)=(/ paz,pbz,pcz,pdz /)
! volmat(:,4)=1D0
! gettetravol=abs(detmat(volmat))/6D0
! call showmatgau(volmat)
!vol=abs( (a-d)¡¤((b-d)¡Á(c-d)) )/6,  see http://en.wikipedia.org/wiki/Tetrahedron
vec1x=pax-pdx
vec1y=pay-pdy
vec1z=paz-pdz
vec2x=pbx-pdx
vec2y=pby-pdy
vec2z=pbz-pdz
vec3x=pcx-pdx
vec3y=pcy-pdy
vec3z=pcz-pdz
call vecprod(vec2x,vec2y,vec2z,vec3x,vec3y,vec3z,vec2x3x,vec2x3y,vec2x3z)
gettetravol=abs(vec1x*vec2x3x+vec1y*vec2x3y+vec1z*vec2x3z)/6D0
end function




!===============================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! Array, Vector !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================!

!!--------- Invert integer array, from istart to iend
!Real*8 version
subroutine invarrr8(array,istart,iend)
integer istart,iend
real*8 array(:)
len=iend-istart+1
do i=0,int(len/2D0)-1
	tmp=array(istart+i)
	array(istart+i)=array(iend-i)
	array(iend-i)=tmp
end do
end subroutine
!Integer 4 version
subroutine invarri4(array,istart,iend)
integer istart,iend,array(:)
len=iend-istart+1
do i=0,int(len/2D0)-1
	tmp=array(istart+i)
	array(istart+i)=array(iend-i)
	array(iend-i)=tmp
end do
end subroutine

!-------- Sort value from small to big by Bubble method
!mode=abs: sort by absolute value, =val: sort by value. Default is by value
!Real*8 version
subroutine sortr8(array,inmode)
integer N,i,j,mode
character,optional :: inmode*3
real*8 array(:),temp
N=size(array)
mode=1
if (present(inmode)) then
	if (inmode=="abs") mode=2
end if
if (mode==1) then
	do i=1,N
		do j=i+1,N
			if (array(i)>array(j)) then
				temp=array(i)
				array(i)=array(j)
				array(j)=temp
			end if
		end do
	end do
else if (mode==2) then
	do i=1,N
		do j=i+1,N
			if (abs(array(i))>abs(array(j))) then
				temp=array(i)
				array(i)=array(j)
				array(j)=temp
			end if
		end do
	end do
end if
end subroutine
!Integer 4 version
subroutine sorti4(array,inmode)
integer N,i,j,mode
character,optional :: inmode*3
integer array(:),temp
N=size(array)
mode=1
if (present(inmode)) then
	if (inmode=="abs") mode=2
end if
if (mode==1) then
	do i=1,N
		do j=i+1,N
			if (array(i)>array(j)) then
				temp=array(i)
				array(i)=array(j)
				array(j)=temp
			end if
		end do
	end do
else if (mode==2) then
	do i=1,N
		do j=i+1,N
			if (abs(array(i))>abs(array(j))) then
				temp=array(i)
				array(i)=array(j)
				array(j)=temp
			end if
		end do
	end do
end if
end subroutine

!!--------- Evaluate standard deviation of array elements
real*8 function stddevarray(array)
real*8 array(:),avg
avg=sum(array)/size(array)
stddevarray=dsqrt(sum((array-avg)**2)/size(array))
end function

!!--------- Evaluate covariant of two array elements
real*8 function covarray(array1,array2)
real*8 array1(:),array2(:),avg1,avg2
avg1=sum(array1)/size(array1)
avg2=sum(array2)/size(array2)
covarray=sum((array1-avg1)*(array2-avg2))/size(array1)
end function

!--- Vector/cross product, input two vectors, return a new vector (x,y,z)
subroutine vecprod(x1,y1,z1,x2,y2,z2,x,y,z)
real*8 x1,y1,z1,x2,y2,z2,x,y,z
! |i  j	 k |
! |x1 y1 z1|
! |x2 y2 z2|
x=	y1*z2-z1*y2
y=-(x1*z2-z1*x2)
z=	x1*y2-y1*x2
end subroutine

!--- Generate full arrangement arry
!ncol is the number of element, nrow=ncol!
!example: nrow=3!=5, ncol=3, the arr will be:
! 1			  3			  2
! 2			  1			  3
! 2			  3			  1
! 3			  1			  2
! 3			  2			  1
! 1			  2			  3
subroutine fullarrange(arr,nrow,ncol)
integer nrow,ncol,arr(nrow,ncol),seq(ncol)
seq=(/ (i,i=1,ncol) /)
arr(1,:)=seq !The first array will be 1,2,3,4...
do icyc=2,nrow
	do i=ncol-1,1,-1
		if (seq(i)<seq(i+1)) exit
	end do
	do j=ncol,1,-1
		if (seq(j)>seq(i)) exit
	end do
	itmp=seq(i)
	seq(i)=seq(j)
	seq(j)=itmp
	call invarr(seq,i+1,ncol)
	arr(icyc,:)=seq
end do
end subroutine



!===============================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! String process !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================!

!-------- Parse inputted integer string (e.g. 3,4,5,6,7-25,32-35,55,88) to array
!nelement will be return, which is the number of indices in the string
!if "array" is present, the terms will be written into it, else if "array" is not present, only count the number of terms as "nelement"
!array must be large enough to contain the elements
subroutine str2arr(inpstr,nelement,array)
character(len=*) inpstr
character c80tmp*80
integer nelement
integer,optional :: array(:)
if (present(array)) array=0
icommaold=0
nelement=0
do ipos=1,len(inpstr)
	if (inpstr(ipos:ipos)==','.or.ipos==len(inpstr)) then
		icommanew=ipos
		if (ipos==len(inpstr)) icommanew=ipos+1
		read(inpstr(icommaold+1:icommanew-1),*) c80tmp
		if (index(c80tmp,'-')==0) then
			nelement=nelement+1
			if (present(array)) read(c80tmp,*) array(nelement)
		else
			do iheng=1,len_trim(c80tmp)
				if (c80tmp(iheng:iheng)=='-') exit
			end do
			read(c80tmp(1:iheng-1),*) ilow
			read(c80tmp(iheng+1:),*) ihigh
			do itmp=ilow,ihigh
				nelement=nelement+1
				if (present(array)) array(nelement)=itmp
			end do
		end if
		icommaold=icommanew
	end if
end do
end subroutine

!---------Input path name, e.g. c:\ltwd\MIO.wfn , output file name, e.g. MIO
subroutine path2filename(pathnamein,filenameout)
character(len=*) pathnamein,filenameout
do i=len_trim(pathnamein),1,-1
	if (pathnamein(i:i)=='.') then
		iend=i-1
		exit
	end if
end do
istart=1
do i=iend,1,-1
	if (pathnamein(i:i)=='/'.or.pathnamein(i:i)=='\') then
		istart=i+1
		exit
	end if
end do
filenameout=' '
filenameout(1:iend-istart+1)=pathnamein(istart:iend)
end subroutine

!!--------- Convert a character to lower case
subroutine uc2lc(inc)
character*1 inc
itmp=ichar(inc)
if (itmp>=65.and.itmp<=90) itmp=itmp+32
inc=char(itmp)
end subroutine

!!--------- Convert a character to upper case
subroutine lc2uc(inc)
character*1 inc
itmp=ichar(inc)
if (itmp>=97.and.itmp<=122) itmp=itmp-32
inc=char(itmp)
end subroutine

!!--------- Convert a string to lower case
subroutine struc2lc(str)
character(len=*) str
do i=1,len_trim(str)
	call uc2lc(str(i:i))
end do
end subroutine

!!--------- Convert a string to upper case
subroutine strlc2uc(str)
character(len=*) str
do i=1,len_trim(str)
	call lc2uc(str(i:i))
end do
end subroutine



!===============================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Matrix calculation !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================!


!!----- Make the matrix to upper trigonal matrix
subroutine ratio_upper(mat)
real*8 :: mat(:,:),m,st
real*8,allocatable :: temp(:),s(:),divided(:)
integer :: a,i,j,t
a=size(mat,1)
allocate(temp(a))
allocate(s(a))
allocate(divided(a))
do i=1,a
	s(i)=maxval(abs(mat(i,1:a)))
end do
do i=1,a-1
	divided(i:a)=mat(i:a,i)/s(i:a)
	t=maxloc(abs(divided(i:a)),dim=1)
	temp(:)=mat(i,:)
	mat(i,:)=mat(i+t-1,:)
	mat(i+t-1,:)=temp(:)
	st=s(i)
	s(i)=s(i+t-1)
	s(i+t-1)=st
	do j=i+1,a
		m=mat(j,i)/mat(i,i)
		mat(j,i:a)=mat(j,i:a)-mat(i,i:a)*m
	end do
end do
deallocate(temp,s,divided)
end subroutine

!!----- Get value of determinant of a matrix
real*8 function detmat(mat)
real*8 mat(:,:)
real*8,allocatable :: mattmp(:,:)
isizemat=size(mat,1)
detmat=1D0
NOTlowertri=0
NOTuppertri=0
outter1: do i=1,isizemat !Check if already is lower-trigonal matrix
	do j=i+1,isizemat
		if (mat(i,j)>1D-12) then
			NOTlowertri=1 !There are at least one big value at upper trigonal part, hence not lower trigonal matrix
			exit outter1
		end if
	end do
end do outter1
outter2: do i=1,isizemat !Check if already is upper-trigonal matrix
	do j=1,i-1
		if (mat(i,j)>1D-12) then
			NOTuppertri=1 !There are at least one big value at lower trigonal part, hence not upper trigonal matrix
			exit outter2
		end if
	end do
end do outter2

if (NOTlowertri==0.or.NOTuppertri==0) then !Is lower or upper trigonal matrix, don't need to convert to trigonal matrix
	do i=1,isizemat
		detmat=detmat*mat(i,i)
	end do
else !Not upper or lower trigonal matrix
	allocate(mattmp(isizemat,isizemat))
	mattmp=mat
	call ratio_upper(mattmp)
	detmat=1D0
	do i=1,isizemat
		detmat=detmat*mattmp(i,i)
	end do
end if
end function

!!-------- Get trace of a matrix
real*8 function mattrace(mat)
real*8 mat(:,:)
mattrace=0
do i=1,size(mat,1)
	mattrace=mattrace+mat(i,i)
end do
end function

!!--- Use Jacobi method to diagonalize matrix, simple, but much slower than diagsymat and diaggemat if the matrix is large
subroutine diagmat(mat,S,eigval,inmaxcyc,inthres)
! mat: input and will be diagonalized matrix, S:eigenvector matrix(columns correspond to vectors), eigval:eigenvalue vector
! inmaxcyc: max cycle, inthres: expected threshold
implicit real*8 (a-h,o-z)
integer,optional :: inmaxcyc
real*8,optional :: inthres
real*8 thres,mat(:,:),S(:,:),eigval(:)
real*8,allocatable :: R(:,:)
n=size(mat,1)
allocate(R(n,n))
maxcyc=200
thres=1D-9
if (present(inmaxcyc)) maxcyc=inmaxcyc
if (present(inthres)) thres=inthres
S=0
do i=1,n
	S(i,i)=1.0D0
end do
do k=1,maxcyc+1
	R=0
	do i=1,n
		R(i,i)=1.0D0
	end do
	i=1
	j=2
	do ii=1,n
		do jj=ii+1,n
			if (abs(mat(ii,jj))>abs(mat(i,j))) then
				i=ii
				j=jj
			end if
		end do
	end do
	if (abs(mat(i,j))<thres) exit
	if (k==maxcyc+1) write(*,*) "Matrix diagonalization exceed max cycle before converge"
	phi=atan(2*mat(i,j)/(mat(i,i)-mat(j,j)))/2.0D0
	R(i,i)=cos(phi)
	R(j,j)=R(i,i)
	R(i,j)=-sin(phi)
	R(j,i)=-R(i,j)
	mat=matmul(matmul(transpose(R),mat),R)
	S=matmul(S,R)
end do
do i=1,n
	eigval(i)=mat(i,i)
end do
end subroutine

!!------------ Diagonalize a symmetry matrix 
!Repack the extremely complex "DSYEV" routine in lapack to terse form
!if istat/=0, means error occurs
subroutine diagsymat(mat,eigvecmat,eigvalarr,istat)
integer istat
real*8 mat(:,:),eigvecmat(:,:),eigvalarr(:)
real*8,allocatable :: lworkvec(:)
isize=size(mat,1)
allocate(lworkvec(3*isize-1))
call DSYEV('V','U',isize,mat,isize,eigvalarr,lworkvec,3*isize-1,istat)
eigvecmat=mat
mat=0D0
forall (i=1:isize) mat(i,i)=eigvalarr(i)
end subroutine

!!------------ Diagonalize a general matrix 
!Repack the extremal complex "DGEEV" routine in lapack to terse form
!eigvecmat is right eigenvector matrix
!eigvalarr is real part of eigenvalue, imaginary parts are discarded
!if istat/=0, means error appears
subroutine diaggemat(mat,eigvecmat,eigvalarr,istat)
integer istat,lwork
real*8 mat(:,:),eigvecmat(:,:),eigvalarr(:),tmpmat(1,1)
real*8,allocatable :: lworkvec(:),eigvalimgarr(:)
isize=size(mat,1)
lwork=8*isize !4*isize is enough, but for better performance we use larger size
allocate(lworkvec(lwork),eigvalimgarr(isize))
call DGEEV('N','V',isize,mat,isize,eigvalarr,eigvalimgarr,tmpmat,1,eigvecmat,isize,lworkvec,lwork,istat)
mat=0D0
forall (i=1:isize) mat(i,i)=eigvalarr(i)
end subroutine

!!--------------- A function to output inverted matrix, inputted matrix will not be affected. Essentially is a warpper of KROUT
function invmat(mat,N)
integer N,ierr
real*8 :: mat(N,N),invmat(N,N),tmpvec(N)
invmat=mat
call KROUT(0,N,0,invmat,N,tmpvec,N,ierr)
end function
!!--------------- A routine to invert matrix. The inputted matrix will be taken placed by inverted matrix. Essentially is a warpper of KROUT
subroutine invmatsub(mat,N)
integer N,ierr
real*8 :: mat(N,N),tmpvec(N)
call KROUT(0,N,0,mat,N,tmpvec,N,ierr)
end subroutine
!Taken and adapted from krout.f, which can be downloaded at http://jblevins.org/mirror/amiller/
!-----------------------------------------------------------------------
!  CROUT PROCEDURE FOR INVERTING MATRICES AND SOLVING EQUATIONS
!-----------------------------------------------------------------------
!  A IS A MATRIX OF ORDER N WHERE N IS GREATER THAN OR EQUAL TO 1.
!  IF MO = 0 THEN THE INVERSE OF A IS COMPUTED AND STORED IN A.
!  IF MO IS NOT 0 THEN THE INVERSE IS NOT COMPUTED.

!  IF M IS GREATER THAN 0 THEN B IS A MATRIX HAVING N ROWS AND M COLUMNS.
!  IN THIS CASE AX = B IS SOLVED AND THE SOLUTION X IS STORED IN B.
!  IF M=0 THEN THERE ARE NO EQUATIONS TO BE SOLVED.
!  N.B. B is passed as a VECTOR not as a matrix.

!  KA = THE LENGTH OF THE COLUMNS OF THE ARRAY A
!  KB = THE LENGTH OF THE COLUMNS OF THE ARRAY B (IF M > 0)

!  IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS. WHEN
!  THE ROUTINE TERMINATES IERR HAS ONE OF THE FOLLOWING VALUES ...
!	  IERR =  0	  THE REQUESTED TASK WAS PERFORMED.
!	  IERR = -1	  EITHER N, KA, OR KB IS INCORRECT.
!	  IERR =  K	  THE K-TH PIVOT ELEMENT IS 0.
!-----------------------------------------------------------------------
! Adapted from the routine KROUT in the NSWC Math. Library by Alan Miller
! Latest revision - 3 August 1998
subroutine KROUT(MO, N, M, A, KA, B, KB, IERR)
implicit none
integer, intent(in)							   :: MO
integer, intent(in)							   :: N
integer, intent(in)							   :: M
real*8, intent(in out), dimension(:,:) :: A		! a(ka,n)
integer, intent(in)							   :: KA
real*8, intent(in out), dimension(:)   :: B
integer, intent(in)							   :: KB
integer, intent(out)						   :: IERR
integer, allocatable, dimension(:)		   :: INDX
real*8, allocatable, dimension(:)  :: TEMP
integer		   :: I, J, JP1, K, KJ, KM1, KP1, L, LJ, MAXB, NJ, NMJ, NMK, NM1, ONEJ
real*8 :: D, DSUM, P, T
real*8, parameter :: ZERO = 0.0D0, ONE = 1.0D0
if (N < 1 .or. KA < N) then
  IERR = -1
  return
end if
if (M > 0 .and. KB < N) then
  IERR = -1
  return
end if
IERR = 0
if (N < 2) then
  D = A(1,1)
  if (D == ZERO) then
	IERR = N
	return
  end if
  if (MO == 0) A(1,1) = ONE / D
  if (M <= 0) return
  MAXB = KB*M
  do KJ = 1,MAXB,KB
	B(KJ) = B(KJ)/D
  end do
  return
end if
if (MO == 0) then
  allocate( INDX(N-1), TEMP(N) )
end if
NM1 = N - 1
do K = 1,NM1
  KP1 = K + 1
  P = abs(A(K,K))
  L = K
  do I = KP1,N
	T = abs(A(I,K))
	if (P >= T) cycle
	P = T
	L = I
  end do
  if (P == ZERO) then
	IERR = K
	return
  end if
  P = A(L,K)
  if (MO == 0) then
	INDX(K) = L
  end if
  if (K /= L) then
	do J = 1,N
	  T = A(K,J)
	  A(K,J) = A(L,J)
	  A(L,J) = T
	end do
	if (M > 0) then
	  KJ = K
	  LJ = L
	  do J = 1,M
		T = B(KJ)
		B(KJ) = B(LJ)
		B(LJ) = T
		KJ = KJ + KB
		LJ = LJ + KB
	  end do
	end if
  end if
  if (K <= 1) then
	do J = KP1,N
	  A(K,J) = A(K,J)/P
	end do
  else
	do J = KP1,N
	  DSUM = A(K,J) - dot_product( A(K,1:KM1), A(1:KM1,J) )
	  A(K,J) = DSUM / P
	end do
  end if
  do I = KP1,N
	DSUM = A(I,KP1) - dot_product( A(I,1:K), A(1:K,KP1) )
	A(I,KP1) = DSUM
  end do
  KM1 = K
end do
if (A(N,N) == ZERO) then
  IERR = N
  return
end if
if (M > 0) then
  MAXB = KB*M
  do ONEJ = 1,MAXB,KB
	KJ = ONEJ
	B(KJ) = B(KJ)/A(1,1)
	do K = 2,N
	  KJ = KJ + 1
	  DSUM = B(KJ)
	  KM1 = K - 1
	  LJ = ONEJ
	  do L = 1,KM1
		DSUM = DSUM - A(K,L)*B(LJ)
		LJ = LJ + 1
	  end do
	  B(KJ) = DSUM / A(K,K)
	end do
  end do
  do NJ = N,MAXB,KB
	KJ = NJ
	do NMK = 1,NM1
	  K = N - NMK
	  LJ = KJ
	  KJ = KJ - 1
	  DSUM = B(KJ)
	  KP1 = K + 1
	  do L = KP1,N
		DSUM = DSUM - A(K,L)*B(LJ)
		LJ = LJ + 1
	  end do
	  B(KJ) = DSUM
	end do
  end do
end if
if (MO /= 0) return
do J = 1,NM1
  A(J,J) = ONE / A(J,J)
  JP1 = J + 1
  do I = JP1,N
	DSUM = dot_product( A(I,J:I-1), A(J:I-1,J) )
	A(I,J) = -DSUM / A(I,I)
  end do
end do
A(N,N) = ONE / A(N,N)
do NMK = 1,NM1
  K = N - NMK
  KP1 = K + 1
  do J = KP1,N
	TEMP(J) = A(K,J)
	A(K,J) = ZERO
  end do
  do J = 1,N
	DSUM = A(K,J) - dot_product( TEMP(KP1:N), A(KP1:N,J) )
	A(K,J) = DSUM
  end do
end do
do NMJ = 1,NM1
  J = N - NMJ
  K = INDX(J)
  if (J == K) cycle
  do I = 1,N
	T = A(I,J)
	A(I,J) = A(I,K)
	A(I,K) = T
  end do
end do
if (MO == 0) deallocate( INDX, TEMP )
end subroutine


!------- Calculate how much is a matrix deviates from identity matrix
!error=¡Æ[i,j]abs( abs(mat(i,j))-¦Ä(i,j) )
real*8 function identmaterr(mat)
implicit real*8 (a-h,o-z)
real*8 mat(:,:)
nsize=size(mat,1)
identmaterr=0D0
do i=1,nsize
	do j=1,nsize
		if (i==j) then
			identmaterr=identmaterr+abs(abs(mat(i,j))-1D0)
		else
			identmaterr=identmaterr+abs(mat(i,j))
		end if
	end do
end do
end function


!----- Convert a square matrix to an array. imode=1/2/3: Full matrix; Lower half matrix; Upper half matrix
!For mode=1,2, "arr" should be nsize*(nsize+1)/2
subroutine mat2arr(mat,arr,imode)
implicit real*8 (a-h,o-z)
real*8 mat(:,:),arr(:)
nsize=size(mat,1)
itmp=0
if (imode==1) then !Full matrix
	do i=1,nsize
		do j=1,nsize
			itmp=itmp+1
			arr(itmp)=mat(i,j)
		end do
	end do
else if (imode==2) then !Lower half matrix
	do i=1,nsize
		do j=1,i
			itmp=itmp+1
			arr(itmp)=mat(i,j)
		end do
	end do
else !Upper half matrix
	do i=1,nsize
		do j=i,nsize
			itmp=itmp+1
			arr(itmp)=mat(i,j)
		end do
	end do
end if
end subroutine


!===============================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Misc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================!

!!---------- Get current time in second, the difference between two times of invoking this routine is consumed wall clock time
subroutine walltime(inow)
character nowdate*20,nowtime*20
integer inow
call date_and_time(nowdate,nowtime)
read(nowtime(1:2),*) inowhour
read(nowtime(3:4),*) inowminute
read(nowtime(5:6),*) inowsecond
inow=inowhour*3600+inowminute*60+inowsecond
end subroutine

!!----- Find the position of specific value in cube
subroutine findvalincub(cubfile,value,i,j,k)
real*8 cubfile(:,:,:),value
integer i,j,k
do ii=1,size(cubfile,1)
	do jj=1,size(cubfile,2)
		do kk=1,size(cubfile,3)
			if (cubfile(ii,jj,kk)==value) then
				i=ii
				j=jj
				k=kk
			end if
		end do
	end do
end do
end subroutine

!!------ Display matrix similar to gaussian program
subroutine showmatgau(mat,label,insemi,form,fileid,useri1,useri2,inncol,titlechar)
!The number of columns is always 5 (unadjustable)
!"Label" is the title, if the content is empty, title will not be printed
!If semi==1, only lower and diagonal element will be shown
!"form" is the format to show data, default is D14.6, can pass into such as "f14.8", total width should be 14 characters
!fildid is output destination, 6 corresponds to outputting to screen
!useri1 and useri2 define the dimension of the matrix, default or -1 means determining automatically
!inncol seems controls spacing between number labels of each frame
!Default titlechar is "i8,6x", if you manually set inncol, you also set this to broaden or narrow
implicit real*8(a-h,o-z)
real*8 :: mat(:,:)
character(*),optional :: label,form,titlechar
integer,optional :: insemi,fileid,useri1,useri2,inncol
integer :: semi,ides,ncol
semi=0
ides=6
ncol=5
i1=size(mat,1)
i2=size(mat,2)
if (present(useri1).and.useri1/=-1) i1=useri1
if (present(useri2).and.useri1/=-1) i2=useri2
if (present(insemi)) semi=insemi
if (present(fileid)) ides=fileid
if (present(inncol)) ncol=inncol
if (present(label).and.label/='') write(ides,*) "************ ",label," ************"
nt=ceiling(i2/float(ncol))
do i=1,nt !How many frame
	ns=(i-1)*5+1 !This frame start from where
	if (i/=nt) ne=(i-1)*ncol+ncol !This frame end to where
	if (i==nt) ne=i2
	!Write basis number in separate line
	write(ides,"(6x)",advance='no')
	do j=ns,ne
		if (present(titlechar)) then
			write(ides,'('//titlechar//')',advance='no') j
		else
			write(ides,"(i8,6x)",advance='no') j
		end if
	end do
	write(ides,*)
	!Write content in each regular line
	do k=1,i1
		if (k<ns.and.semi==1) cycle !The lines have been outputted are skipped
		write(ides,"(i6)",advance='no') k
		do j=ns,ne
			if (semi==1.and.k<j) cycle !Upper trigonal element were passed
			if (present(form)) then
				write(ides,'('//form//')',advance='no') mat(k,j)			
			else
				write(ides,"(D14.6)",advance='no') mat(k,j)
			end if
		end do
		write(ides,*) !Change to next line
	end do
end do
end subroutine

!!------- Read matrix outputted in the gaussian .out file, such as outputted by iop(3/33=1)
subroutine readmatgau(fileid,mat,semi,inform,inskipcol,inncol,innspace)
!e.g. trigonal: readmatgau(10,tempmat,1,"D14.6",7,5)   full: readmatgau(10,tempmat,0,"D13.6",7,5)
!inform is the format used to read data, default is D14.6, you can input "f8.4 " (Note: Must be 5 character!)
!semi=1 means the read matrix is lower trigonal matrix
!inskipcol is the skipped column (marker information) in front of each row, default is 7
!inncol is data in each row, default is 5
!innspace is space lines in between each frame, default is 1

!Before use, loclabel should be used to move reading position into title line of the matrix£¬namely move to "*** Overlap ***"
!This routine will skip title line, and skip one line in each frame reading
! *** Overlap ***
!				  1				2			  3				4			  5
!		1  0.100000D+01
!		2  0.236704D+00	 0.100000D+01
implicit real*8(a-h,o-z)
real*8 :: mat(:,:)
character,optional :: inform*5
character :: form*7,c80tmp*79
integer,optional :: inskipcol,inncol,semi,innspace
integer :: skipcol,ncol,fileid
form="(D14.6)" !Suitable for Sbas,Kinene,Potene,Hcore by 3/33=1 ,Fockmat,Densmat by 5/33=3
skipcol=7
ncol=5
nspace=1
if (present(inform)) form='('//inform//')'
if (present(inskipcol)) skipcol=inskipcol
if (present(inncol)) ncol=inncol
if (present(innspace)) nspace=innspace
read(fileid,*) !!!!!!!!!!!! Skip title line
i1=size(mat,1)
i2=size(mat,2)
nt=ceiling(i2/float(ncol))
mat=0D0
do i=1,nt !Number of frames
	ns=(i-1)*ncol+1
	if (i/=nt) ne=(i-1)*ncol+ncol
	if (i==nt) ne=i2
	do ii=1,nspace
		read(fileid,*) !!!!!!!!!!!! Skip number line when reading each new frame
	end do
	do k=1,i1 !Scan rows in each frame
!		  read(fileid,"(a)") c80tmp
!		  write(15,"(a)") c80tmp
!		  backspace(fileid)
		
		if (k<ns.and.present(semi).and.semi==1) cycle
		do iii=1,skipcol !Skip marker columns in each row
			read(fileid,"(1x)",advance='no')
		end do
		do j=ns,ne !Scan elements in each row
			if (present(semi).and.semi==1.and.k<j) cycle
			read(fileid,form,advance='no') mat(k,j)
!			  write(*,*) i,k,j,mat(k,j)
		end do
		read(fileid,*)
	end do
end do
if (present(semi).and.semi==1) then !When read is lower trigonal matrix, we assume it is a symmetric matrix
	mat=mat+transpose(mat)
	do i=1,i1
		mat(i,i)=mat(i,i)/2D0
	end do
end if
end subroutine

!!!------------------------- Determine how many lines in the fileid
!If imode==1, space line will be regarded as the sign of end of file. If imode==2, will count actual number of lines in the file
integer function totlinenum(fileid,imode)
integer fileid,ierror,imode
character*80 c80
totlinenum=0
do while(.true.)
	read(fileid,"(a)",iostat=ierror) c80
	if (imode==1) then
		if (ierror/=0.or.c80==" ") exit
	else if (imode==2) then
		if (ierror/=0) exit
	end if
	totlinenum=totlinenum+1
end do
rewind(fileid)
end function

!!!-------- Locate the line where the label first appears in fileid
!Return ifound=1 if found the label, else return 0
!Default is rewind, if irewind=0 then will not rewind
!If the current line just has the label, calling this subroutine will do nothing
!maxline define the maximum number of lines that will be searched, default is search the whole file
subroutine loclabel(fileid,label,ifound,irewind,maxline)
integer fileid,ierror
integer,optional :: ifound,irewind,maxline
character*200 c200
CHARACTER(LEN=*) label
if ((.not.present(irewind)).or.(present(irewind).and.irewind==1)) rewind(fileid)
if (.not.present(maxline)) then
	do while(.true.)
		read(fileid,"(a)",iostat=ierror) c200
		if (index(c200,label)/=0) then
			backspace(fileid)
			if (present(ifound)) ifound=1 !Found result
			return
		end if
		if (ierror/=0) exit
	end do
else
	do iline=1,maxline
		read(fileid,"(a)",iostat=ierror) c200
		if (index(c200,label)/=0) then
			backspace(fileid)
			if (present(ifound)) ifound=1 !Found result
			return
		end if
		if (ierror/=0) exit
	end do
end if
if (present(ifound)) ifound=0
end subroutine


!!!----------- Skip specific number of lines in specific fileid
subroutine skiplines(id,nskip)
integer id,nskip
do i=1,nskip
	read(id,*)
end do
end subroutine


!!!---------------- Calculate factorial
integer function ft(i)
integer i
ft=i
if (i==0) ft=1
do j=i-1,1,-1
	ft=ft*j
end do
end function

!---- Calculate gamma(Lval+1/2), see http://en.wikipedia.org/wiki/Gamma_function
real*8 function gamma_ps(n)
use defvar
integer n
gamma_ps=ft(2*n)*dsqrt(pi)/4**n/ft(n)
end function

!!-------- Get all combinations of any ncomb elements of array, which length is ntot
!outarray(A,B) is output array, A is the length and must be ntot!/ncomb!/(ntot-ncomb)!, B is generated array, should be equal to ncomb
subroutine combarray(array,ntot,ncomb,outarray)
integer array(ntot),ntot,ncomb,idxarr(ncomb),outarray(:,:)
ipos=ncomb !Current position in the array
forall (i=1:ncomb) idxarr(i)=i !Used to record index
ioutput=1
ncount=0
do while(ipos>0)
	if (ioutput==1) then
		ncount=ncount+1
		outarray(ncount,:)=array(idxarr(:))
	end if
	ioutput=0
	idxarr(ipos)=idxarr(ipos)+1
	if (idxarr(ipos)>ntot) then
		ipos=ipos-1 !Go back to last position
		cycle
	end if
	if (ipos<ncomb) then
		ipos=ipos+1
		idxarr(ipos)=idxarr(ipos-1)
		cycle
	end if
	if (ipos==ncomb) ioutput=1
end do
end subroutine


!!--------- Convert XY scatter data to density distribution
!xarr and yarr records the points. nlen is array length
!mat is the outputted matrix, matnx and matny are its number of element in X and Y
!x/ymin, x/ymax are lower and upper limit of "mat", the X and Y ranges contain nvalx and nvaly data
!e.g. n=5
!	|	1	|	2	|	3	|	4	|	5	 |
! xmin-------------------------------------xmax
subroutine xypt2densmat(xarr,yarr,nlen,mat,nvalx,nvaly,xmin,xmax,ymin,ymax)
integer nlen,nvalx,nvaly
real*8 xarr(nlen),yarr(nlen),mat(nvalx,nvaly),xmin,xmax,ymin,ymax
mat=0D0
spcx=(xmax-xmin)/nvalx
spcy=(ymax-ymin)/nvaly
!If enable parallel, program often prompts memory is not enough, I don't know how to solve this
! !$OMP PARALLEL DO SHARED(mat) PRIVATE(ix,iy,ipt,xlow,xhigh,ylow,yhigh) schedule(dynamic) NUM_THREADS(nthreads)
do ix=1,nvalx
	xlow=xmin+(ix-1)*spcx
	xhigh=xmin+ix*spcx
	do iy=1,nvaly
		ylow=ymin+(iy-1)*spcy
		yhigh=ymin+iy*spcy
		do ipt=1,nlen
			if (xarr(ipt)>xlow.and.xarr(ipt)<=xhigh.and.yarr(ipt)>ylow.and.yarr(ipt)<=yhigh) then
				mat(ix,iy)=mat(ix,iy)+1D0
			end if
		end do
	end do
end do
! !$OMP END PARALLEL DO
end subroutine

end module

module function
use defvar
implicit real*8 (a-h,o-z)

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate wavefunction value of a range of orbitals and their derivatives at a given point, up to third-order
!! istart and iend is the range of the orbitals will be calculated, to calculate all orbitals, use 1,nmo
!! runtype=1: value	 =2: value+dx/y/z  =3: value+dxx/yy/zz(diagonal of hess)  =4: value+dx/y/z+Hessian	
!!		  =5: value+dx/y/z+hess+3-order derivative tensor 
subroutine orbderv(runtype,istart,iend,x,y,z,wfnval,grad,hess,tens3)
real*8 x,y,z,wfnval(nmo)
real*8,optional :: grad(3,nmo),hess(3,3,nmo),tens3(3,3,3,nmo)
integer runtype,istart,iend

wfnval=0D0
if (present(grad)) grad=0D0
if (present(hess)) hess=0D0
if (present(tens3)) tens3=0D0
lastcen=-1 !arbitrary value

! if the center/exp of current GTF is the same as previous, then needn't recalculate them
do j=1,nprims
	ix=type2ix(b(j)%functype)
	iy=type2iy(b(j)%functype)
	iz=type2iz(b(j)%functype)
	ep=b(j)%exp
	
	if (b(j)%center/=lastcen) then
		sftx=x-a(b(j)%center)%x
		sfty=y-a(b(j)%center)%y
		sftz=z-a(b(j)%center)%z
		sftx2=sftx*sftx
		sfty2=sfty*sfty
		sftz2=sftz*sftz
		rr=sftx2+sfty2+sftz2
	end if
	if (expcutoff>0.or.-ep*rr>expcutoff) then
		expterm=exp(-ep*rr)
	else
		expterm=0D0
	end if
	lastcen=b(j)%center
!	  expterm=exp(-ep*dsqrt(rr))
	if (expterm==0D0) cycle
	
	!Calculate value for current GTF
	if (b(j)%functype==1) then !Some functype use manually optimized formula for cutting down computational time
	GTFval=expterm
	else if (b(j)%functype==2) then
	GTFval=sftx*expterm
	else if (b(j)%functype==3) then
	GTFval=sfty*expterm
	else if (b(j)%functype==4) then
	GTFval=sftz*expterm
	else if (b(j)%functype==5) then
	GTFval=sftx2*expterm
	else if (b(j)%functype==6) then
	GTFval=sfty2*expterm
	else if (b(j)%functype==7) then
	GTFval=sftz2*expterm
	else if (b(j)%functype==8) then
	GTFval=sftx*sfty*expterm
	else if (b(j)%functype==9) then
	GTFval=sftx*sftz*expterm
	else if (b(j)%functype==10) then
	GTFval=sfty*sftz*expterm
	else !If above condition is not satisfied(Angular moment higher than f), the function will calculated explicitly
	GTFval=sftx**ix *sfty**iy *sftz**iz *expterm
	end if
	!Calculate orbital wavefunction value
	do imo=istart,iend
		wfnval(imo)=wfnval(imo)+co(imo,j)*GTFval
	end do
	
	if (runtype>=2) then
		!Calculate 1-order derivative for current GTF
		tx=0.0D0
		ty=0.0D0
		tz=0.0D0
		if (ix/=0) tx=ix*sftx**(ix-1)
		if (iy/=0) ty=iy*sfty**(iy-1)
		if (iz/=0) tz=iz*sftz**(iz-1)
		GTFdx=sfty**iy *sftz**iz *expterm*(tx-2*ep*sftx**(ix+1))
		GTFdy=sftx**ix *sftz**iz *expterm*(ty-2*ep*sfty**(iy+1))
		GTFdz=sftx**ix *sfty**iy *expterm*(tz-2*ep*sftz**(iz+1))
		!Calculate 1-order derivative for orbitals
		do imo=istart,iend
			grad(1,imo)=grad(1,imo)+co(imo,j)*GTFdx
			grad(2,imo)=grad(2,imo)+co(imo,j)*GTFdy
			grad(3,imo)=grad(3,imo)+co(imo,j)*GTFdz
		end do

		if (runtype>=3) then
			!Calculate 2-order derivative for current GTF
			txx=0.0D0
			tyy=0.0D0
			tzz=0.0D0
			if (ix>=2) txx=ix*(ix-1)*sftx**(ix-2)
			if (iy>=2) tyy=iy*(iy-1)*sfty**(iy-2)
			if (iz>=2) tzz=iz*(iz-1)*sftz**(iz-2)
			GTFdxx=sfty**iy *sftz**iz *expterm*( txx + 2*ep*sftx**ix*(-2*ix+2*ep*sftx2-1) )
			GTFdyy=sftx**ix *sftz**iz *expterm*( tyy + 2*ep*sfty**iy*(-2*iy+2*ep*sfty2-1) )
			GTFdzz=sftx**ix *sfty**iy *expterm*( tzz + 2*ep*sftz**iz*(-2*iz+2*ep*sftz2-1) )
			ttx=tx-2*ep*sftx**(ix+1)
			tty=ty-2*ep*sfty**(iy+1)
			ttz=tz-2*ep*sftz**(iz+1)
			GTFdxy=sftz**iz *expterm*ttx*tty
			GTFdyz=sftx**ix *expterm*tty*ttz
			GTFdxz=sfty**iy *expterm*ttx*ttz
			!Calculate diagonal Hessian elements for orbitals
			do imo=istart,iend
				hess(1,1,imo)=hess(1,1,imo)+co(imo,j)*GTFdxx !dxx
				hess(2,2,imo)=hess(2,2,imo)+co(imo,j)*GTFdyy !dyy
				hess(3,3,imo)=hess(3,3,imo)+co(imo,j)*GTFdzz !dzz
			end do
			if (runtype>=4) then !Also process nondiagonal elements
				do imo=istart,iend
					hess(1,2,imo)=hess(1,2,imo)+co(imo,j)*GTFdxy !dxy
					hess(2,3,imo)=hess(2,3,imo)+co(imo,j)*GTFdyz !dyz
					hess(1,3,imo)=hess(1,3,imo)+co(imo,j)*GTFdxz !dxz
				end do
				hess(2,1,:)=hess(1,2,:)
				hess(3,2,:)=hess(2,3,:)
				hess(3,1,:)=hess(1,3,:)
			end if
			
			if (runtype>=5) then
				!Calculate 3-order derivative for current GTF
				ep2=ep*2D0
				ep4=ep*4D0
				epep4=ep2*ep2
				epep8=epep4*2D0
				!dxyz
				a1=0D0
				b1=0D0
				c1=0D0
				if (ix>=1) a1=ix*sftx**(ix-1)
				if (iy>=1) b1=iy*sfty**(iy-1)
				if (iz>=1) c1=iz*sftz**(iz-1)
				a2=-ep2*sftx**(ix+1)
				b2=-ep2*sfty**(iy+1)
				c2=-ep2*sftz**(iz+1)
				GTFdxyz=(a1+a2)*(b1+b2)*(c1+c2)*expterm
				!dxyy,dxxy,dxxz,dxzz,dyzz,dyyz
				atmp=0D0
				btmp=0D0
				ctmp=0D0
				if (ix>=2) atmp=ix*(ix-1)*sftx**(ix-2)
				if (iy>=2) btmp=iy*(iy-1)*sfty**(iy-2)
				if (iz>=2) ctmp=iz*(iz-1)*sftz**(iz-2)
				GTFdxyy=(a1+a2)*sftz**iz *expterm*(-ep4*iy*sfty**iy+epep4*sfty**(iy+2)+btmp-ep2*sfty**iy) !ok
				GTFdxxy=(b1+b2)*sftz**iz *expterm*(-ep4*ix*sftx**ix+epep4*sftx**(ix+2)+atmp-ep2*sftx**ix) !=dyxx
				GTFdxxz=(c1+c2)*sfty**iy *expterm*(-ep4*ix*sftx**ix+epep4*sftx**(ix+2)+atmp-ep2*sftx**ix) !=dzxx
				GTFdxzz=(a1+a2)*sfty**iy *expterm*(-ep4*iz*sftz**iz+epep4*sftz**(iz+2)+ctmp-ep2*sftz**iz)
				GTFdyzz=(b1+b2)*sftx**ix *expterm*(-ep4*iz*sftz**iz+epep4*sftz**(iz+2)+ctmp-ep2*sftz**iz)
				GTFdyyz=(c1+c2)*sftx**ix *expterm*(-ep4*iy*sfty**iy+epep4*sfty**(iy+2)+btmp-ep2*sfty**iy) !=dzyy,ok
				!dxxx,dyyy,dzzz
				aatmp1=0D0
				bbtmp1=0D0
				cctmp1=0D0
				if (ix>=1) aatmp1=ep2*ix*sftx**(ix-1)
				if (iy>=1) bbtmp1=ep2*iy*sfty**(iy-1)
				if (iz>=1) cctmp1=ep2*iz*sftz**(iz-1)
				aatmp2=0D0
				bbtmp2=0D0
				cctmp2=0D0
				if (ix>=2) aatmp2=ep2*ix*(ix-1)*sftx**(ix-1)
				if (iy>=2) bbtmp2=ep2*iy*(iy-1)*sfty**(iy-1)
				if (iz>=2) cctmp2=ep2*iz*(iz-1)*sftz**(iz-1)
				aatmp3=0D0
				bbtmp3=0D0
				cctmp3=0D0
				if (ix>=3) aatmp3=ix*(ix-1)*(ix-2)*sftx**(ix-3)
				if (iy>=3) bbtmp3=iy*(iy-1)*(iy-2)*sfty**(iy-3)
				if (iz>=3) cctmp3=iz*(iz-1)*(iz-2)*sftz**(iz-3)
				GTFdxxx=sfty**iy*sftz**iz*expterm*( (-2*ix+ep2*sftx2-1)*(-epep4*sftx**(ix+1) + aatmp1) - aatmp2 + epep8*sftx**(ix+1) + aatmp3 )
				GTFdyyy=sftx**ix*sftz**iz*expterm*( (-2*iy+ep2*sfty2-1)*(-epep4*sfty**(iy+1) + bbtmp1) - bbtmp2 + epep8*sfty**(iy+1) + bbtmp3 )
				GTFdzzz=sfty**iy*sftx**ix*expterm*( (-2*iz+ep2*sftz2-1)*(-epep4*sftz**(iz+1) + cctmp1) - cctmp2 + epep8*sftz**(iz+1) + cctmp3 )
				
				!Calculate 3-order derivative tensor for orbital wavefunction
				do imo=istart,iend
					tens3(1,1,1,imo)=tens3(1,1,1,imo)+co(imo,j)*GTFdxxx !dxxx
					tens3(2,2,2,imo)=tens3(2,2,2,imo)+co(imo,j)*GTFdyyy !dyyy
					tens3(3,3,3,imo)=tens3(3,3,3,imo)+co(imo,j)*GTFdzzz !dzzz
					tens3(1,2,2,imo)=tens3(1,2,2,imo)+co(imo,j)*GTFdxyy !dxyy*
					tens3(1,1,2,imo)=tens3(1,1,2,imo)+co(imo,j)*GTFdxxy !dxxy*
					tens3(1,1,3,imo)=tens3(1,1,3,imo)+co(imo,j)*GTFdxxz !dxxz*
					tens3(1,3,3,imo)=tens3(1,3,3,imo)+co(imo,j)*GTFdxzz !dxzz*
					tens3(2,3,3,imo)=tens3(2,3,3,imo)+co(imo,j)*GTFdyzz !dyzz*
					tens3(2,2,3,imo)=tens3(2,2,3,imo)+co(imo,j)*GTFdyyz !dyyz*
					tens3(1,2,3,imo)=tens3(1,2,3,imo)+co(imo,j)*GTFdxyz !dxyz
				end do
				tens3(1,2,1,:)=tens3(1,1,2,:) !dxyx=dxxy
				tens3(1,3,1,:)=tens3(1,1,3,:) !dxzx=dxxz
				tens3(1,3,2,:)=tens3(1,2,3,:) !dxzy=dxyz
				tens3(2,1,1,:)=tens3(1,1,2,:) !dyxx=dxxy
				tens3(2,1,2,:)=tens3(1,2,2,:) !dyxy=dxyy
				tens3(2,1,3,:)=tens3(1,2,3,:) !dyxz=dxyz
				tens3(2,2,1,:)=tens3(1,2,2,:) !dyyx=dxyy
				tens3(2,3,1,:)=tens3(1,2,3,:) !dyzx=dxyz
				tens3(2,3,2,:)=tens3(2,2,3,:) !dyzy=dyyz
				tens3(3,1,1,:)=tens3(1,1,3,:) !dzxx=dxxz
				tens3(3,1,2,:)=tens3(1,2,3,:) !dzxy=dxyz
				tens3(3,1,3,:)=tens3(1,3,3,:) !dzxz=dxzz
				tens3(3,2,1,:)=tens3(1,2,3,:) !dzyx=dxyz
				tens3(3,2,2,:)=tens3(2,2,3,:) !dzyy=dyyz
				tens3(3,2,3,:)=tens3(2,3,3,:) !dzyz=dyzz
				tens3(3,3,1,:)=tens3(1,3,3,:) !dzzx=dxzz
				tens3(3,3,2,:)=tens3(2,3,3,:) !dzzy=dyzz
			end if !end runtype>=5
			
		end if !end runtype>=3
	end if !end runtype>=2
end do
end subroutine


!!!----------- Calculate contribution from EDFs (recorded in wfx file) to density and corresponding derivatives (up to third-order)
!Only S-type GTFs are supported
! In wfx files, GTFs are used to expand core density
! runtype=1: Only calculate rho, =2: rho+dx/dy/dz =3: rho+dx/dy/dz+dxx/dyy/dzz
!		 =4: rho+dx/dy/dz+full Hessian =5: rho+dx/dy/dz+full Hessian+tens3
subroutine EDFrho(runtype,x,y,z,value,grad,hess,tens3)
integer runtype
real*8 x,y,z,value
real*8,optional :: grad(3),hess(3,3),tens3(3,3,3)
value=0D0
if (present(grad)) grad=0D0
if (present(hess)) hess=0D0
if (present(tens3)) tens3=0D0

do i=1,nEDFprims
	sftx=x-a(b_EDF(i)%center)%x
	sfty=y-a(b_EDF(i)%center)%y
	sftz=z-a(b_EDF(i)%center)%z
	sftx2=sftx*sftx
	sfty2=sfty*sfty
	sftz2=sftz*sftz
	rr=sftx2+sfty2+sftz2
	ep=b_EDF(i)%exp
	expterm=exp(-ep*rr)
	value=value+CO_EDF(i)*expterm
!	  write(11,"(i5,3D20.10)") i,value,CO_EDF(i),expterm
	if (runtype>=2) then
		tmp=2*CO_EDF(i)*expterm*ep
		grad(1)=grad(1)-tmp*sftx
		grad(2)=grad(2)-tmp*sfty
		grad(3)=grad(3)-tmp*sftz
		if (runtype>=3) then
			hess(1,1)=hess(1,1)+tmp*(2*ep*sftx2-1)
			hess(2,2)=hess(2,2)+tmp*(2*ep*sfty2-1)
			hess(3,3)=hess(3,3)+tmp*(2*ep*sftz2-1)
			if (runtype>=4) then
				epep4=ep*ep*4
				tmp2=CO_EDF(i)*epep4*expterm
				hess(1,2)=hess(1,2)+tmp2*sftx*sfty
				hess(1,3)=hess(1,3)+tmp2*sftx*sftz
				hess(2,3)=hess(2,3)+tmp2*sfty*sftz
				hess(2,1)=hess(1,2)
				hess(3,1)=hess(1,3)
				hess(3,2)=hess(2,3)
				if (runtype>=5) then
					tmp3=CO_EDF(i)*epep4*expterm
					tens3(1,1,1)=tens3(1,1,1)+tmp3*sftx*(3-2*ep*sftx2)
					tens3(2,2,2)=tens3(2,2,2)+tmp3*sfty*(3-2*ep*sfty2)
					tens3(3,3,3)=tens3(3,3,3)+tmp3*sftz*(3-2*ep*sftz2)
					tens3(1,2,2)=tens3(1,2,2)+tmp3*sftx*(1-2*ep*sfty2)
					tens3(1,1,2)=tens3(1,1,2)+tmp3*sfty*(1-2*ep*sftx2)
					tens3(1,1,3)=tens3(1,1,3)+tmp3*sftz*(1-2*ep*sftx2)
					tens3(1,3,3)=tens3(1,3,3)+tmp3*sftx*(1-2*ep*sftz2)
					tens3(2,3,3)=tens3(2,3,3)+tmp3*sfty*(1-2*ep*sftz2)
					tens3(2,2,3)=tens3(2,2,3)+tmp3*sftz*(1-2*ep*sfty2)
					tens3(1,2,3)=tens3(1,2,3)-CO_EDF(i)*8*ep**3*sftx*sfty*sftz*expterm
					tens3(1,2,1)=tens3(1,1,2) !dxyx=dxxy
					tens3(1,3,1)=tens3(1,1,3) !dxzx=dxxz
					tens3(1,3,2)=tens3(1,2,3) !dxzy=dxyz
					tens3(2,1,1)=tens3(1,1,2) !dyxx=dxxy
					tens3(2,1,2)=tens3(1,2,2) !dyxy=dxyy
					tens3(2,1,3)=tens3(1,2,3) !dyxz=dxyz
					tens3(2,2,1)=tens3(1,2,2) !dyyx=dxyy
					tens3(2,3,1)=tens3(1,2,3) !dyzx=dxyz
					tens3(2,3,2)=tens3(2,2,3) !dyzy=dyyz
					tens3(3,1,1)=tens3(1,1,3) !dzxx=dxxz
					tens3(3,1,2)=tens3(1,2,3) !dzxy=dxyz
					tens3(3,1,3)=tens3(1,3,3) !dzxz=dxzz
					tens3(3,2,1)=tens3(1,2,3) !dzyx=dxyz
					tens3(3,2,2)=tens3(2,2,3) !dzyy=dyyz
					tens3(3,2,3)=tens3(2,3,3) !dzyz=dyzz
					tens3(3,3,1)=tens3(1,3,3) !dzzx=dxzz
					tens3(3,3,2)=tens3(2,3,3) !dzzy=dyzz
				end if
			end if
		end if
	end if
end do
end subroutine

real*8 function nucesp(x,y,z)
nucesp=0D0
do i=1,nfragatmnum
	dist2mpx=(x-a(fragatm(i))%x)**2
	dist2mpy=(y-a(fragatm(i))%y)**2
	dist2mpz=(z-a(fragatm(i))%z)**2
	dist2=dist2mpx+dist2mpy+dist2mpz
	if (dist2==0D0) then
		 nucesp=1D3
		 return
	end if
	nucesp=nucesp+a(fragatm(i))%charge/dsqrt(dist2)
end do
end function

! Kinetic Energy
real*8 function Hamkin(x,y,z)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo)
call orbderv(3,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
hamx=sum( MOocc(1:nmo)*wfnhess(1,1,1:nmo)*wfnval(1:nmo) )
hamy=sum( MOocc(1:nmo)*wfnhess(2,2,1:nmo)*wfnval(1:nmo) )
hamz=sum( MOocc(1:nmo)*wfnhess(3,3,1:nmo)*wfnval(1:nmo) )
Hamkin=hamx+hamy+hamz
Hamkin=-Hamkin/2D0
end function

! Nuclear Attraction
real*8 function calcnuc(x,y,z)
real*8 x,y,z
calcnuc=-nucesp(x,y,z)*fdens(x,y,z)
end function

!!--- Simultaneously generate electron density, gradient norm for alpha and beta electrons, as well as dot product between grada and gradb
!---- Mainly used to evalute DFT functional. EDF is not taken into account
!adens/bdens/tdens means the density of alpha/beta/total density, similar for *grad
subroutine gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,abgrad,tgrad,atau,btau,ttau,alap,blap,tlap)
real*8 x,y,z,adens,bdens,agrad,bgrad,abgrad,wfnval(nmo),wfnderv(3,nmo),gradrhoa(3),gradrhob(3),gradrhot(3),tmparr(3),tmparr2(3),tmparr3(3,3)
real*8 wfnhess(3,3,nmo),EDFgrad(3),EDFhess(3,3),ttau,atau,btau,lagkina(3),lagkinb(3),alap,blap,tlap,hessrhoa(3,3),hessrhob(3,3)

call orbderv(4,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
adens=0D0
bdens=0D0
gradrhoa=0D0
gradrhob=0D0
atau=0D0
btau=0D0
lagkina=0D0
lagkinb=0D0
hessrhoa=0d0
hessrhob=0d0
alap=0d0
blap=0d0
do i=1,nmo
	if (MOtype(i)==1) then
		adens=adens+MOocc(i)*wfnval(i)**2
		gradrhoa(:)=gradrhoa(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
		lagkina(:)=lagkina(:)+MOocc(i)*wfnderv(:,i)**2
        hessrhoa(1,1)=hessrhoa(1,1)+MOocc(i)*( wfnderv(1,i)**2 + wfnval(i)*wfnhess(1,1,i) )
        hessrhoa(2,2)=hessrhoa(2,2)+MOocc(i)*( wfnderv(2,i)**2 + wfnval(i)*wfnhess(2,2,i) )
        hessrhoa(3,3)=hessrhoa(3,3)+MOocc(i)*( wfnderv(3,i)**2 + wfnval(i)*wfnhess(3,3,i) )
!        hessrhoa(1,2)=hessrhoa(1,2)+MOocc(i)*( wfnderv(1,i)*wfnderv(2,i)+wfnhess(1,2,i)*wfnval(i) )
!        hessrhoa(2,3)=hessrhoa(2,3)+MOocc(i)*( wfnderv(2,i)*wfnderv(3,i)+wfnhess(2,3,i)*wfnval(i) )
!        hessrhoa(1,3)=hessrhoa(1,3)+MOocc(i)*( wfnderv(1,i)*wfnderv(3,i)+wfnhess(1,3,i)*wfnval(i) )
	else if (MOtype(i)==2) then
		bdens=bdens+MOocc(i)*wfnval(i)**2
		gradrhob(:)=gradrhob(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
		lagkinb(:)=lagkinb(:)+MOocc(i)*wfnderv(:,i)**2
        hessrhob(1,1)=hessrhob(1,1)+MOocc(i)*( wfnderv(1,i)**2 + wfnval(i)*wfnhess(1,1,i) )
        hessrhob(2,2)=hessrhob(2,2)+MOocc(i)*( wfnderv(2,i)**2 + wfnval(i)*wfnhess(2,2,i) )
        hessrhob(3,3)=hessrhob(3,3)+MOocc(i)*( wfnderv(3,i)**2 + wfnval(i)*wfnhess(3,3,i) )
!        hessrhob(1,2)=hessrhob(1,2)+MOocc(i)*( wfnderv(1,i)*wfnderv(2,i)+wfnhess(1,2,i)*wfnval(i) )
!        hessrhob(2,3)=hessrhob(2,3)+MOocc(i)*( wfnderv(2,i)*wfnderv(3,i)+wfnhess(2,3,i)*wfnval(i) )
!        hessrhob(1,3)=hessrhob(1,3)+MOocc(i)*( wfnderv(1,i)*wfnderv(3,i)+wfnhess(1,3,i)*wfnval(i) )		
	else if (MOtype(i)==0) then
		tmpval=MOocc(i)/2D0*wfnval(i)**2
		adens=adens+tmpval
		bdens=bdens+tmpval
		tmparr(:)=MOocc(i)/2D0*wfnval(i)*wfnderv(:,i)
		gradrhoa(:)=gradrhoa(:)+tmparr(:)
		gradrhob(:)=gradrhob(:)+tmparr(:)
		tmparr2(:)=MOocc(i)/2D0*wfnderv(:,i)**2
		lagkina(:)=lagkina(:)+tmparr2(:)
		lagkinb(:)=lagkinb(:)+tmparr2(:)
        tmparr3(1,1)=MOocc(i)/2D0*( wfnderv(1,i)**2 + wfnval(i)*wfnhess(1,1,i) )
        tmparr3(2,2)=MOocc(i)/2D0*( wfnderv(2,i)**2 + wfnval(i)*wfnhess(2,2,i) )
        tmparr3(3,3)=MOocc(i)/2D0*( wfnderv(3,i)**2 + wfnval(i)*wfnhess(3,3,i) )
!        tmparr3(1,2)=MOocc(i)/2D0*( wfnderv(1,i)*wfnderv(2,i)+wfnhess(1,2,i)*wfnval(i) )
!        tmparr3(2,3)=MOocc(i)/2D0*( wfnderv(2,i)*wfnderv(3,i)+wfnhess(2,3,i)*wfnval(i) )
!        tmparr3(1,3)=MOocc(i)/2D0*( wfnderv(1,i)*wfnderv(3,i)+wfnhess(1,3,i)*wfnval(i) )
        hessrhoa(1,1)=hessrhoa(1,1)+tmparr3(1,1)
        hessrhoa(2,2)=hessrhoa(2,2)+tmparr3(2,2)
        hessrhoa(3,3)=hessrhoa(3,3)+tmparr3(3,3)
!        hessrhoa(1,2)=hessrhoa(1,2)+tmparr3(1,2)
!        hessrhoa(2,3)=hessrhoa(2,3)+tmparr3(2,3)
!        hessrhoa(1,3)=hessrhoa(1,3)+tmparr3(1,3)
        hessrhob(1,1)=hessrhob(1,1)+tmparr3(1,1)
        hessrhob(2,2)=hessrhob(2,2)+tmparr3(2,2)
        hessrhob(3,3)=hessrhob(3,3)+tmparr3(3,3)
!        hessrhob(1,2)=hessrhob(1,2)+tmparr3(1,2)
!        hessrhob(2,3)=hessrhob(2,3)+tmparr3(2,3)
!        hessrhob(1,3)=hessrhob(1,3)+tmparr3(1,3)		
	end if
end do
tdens=adens+bdens
gradrhoa=gradrhoa*2
gradrhob=gradrhob*2
gradrhot=gradrhoa+gradrhob
agrad=dsqrt(sum(gradrhoa**2))
bgrad=dsqrt(sum(gradrhob**2))
tgrad=dsqrt(sum(gradrhot**2))
abgrad=sum(gradrhoa*gradrhob)
atau=sum(lagkina)/2.0D0
btau=sum(lagkinb)/2.0D0
ttau=atau+btau
alap=hessrhoa(1,1)+hessrhoa(2,2)+hessrhoa(3,3)
blap=hessrhob(1,1)+hessrhob(2,2)+hessrhob(3,3)
tlap=alap+blap

!write(*,"(' ----RPV----	   :',f10.3,f10.3,f10.3,f40.10,f40.10,f40.10)") x,y,z,agrad**2,alap
!if (tdens.gt.1D-6) write(*,"(' ----RPV----	   :',f10.3,f10.3,f10.3,f40.10,f40.10,f40.10)") x,y,z,tdens,tgrad**2,ttau

end subroutine




!        !!!------------------------- Output Laplacian of electron density at a point
!        !label=x/y/z output 2-order derivative of electron density respect to xx/yy/zz, =t get their summing
!        real*8 function flapl(x,y,z,label)
!        real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),laplx,laply,laplz,EDFgrad(3),EDFhess(3,3)
!        character label
!        call orbderv(3,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
!        laplx=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
!        laply=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
!        laplz=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
!        !Add in contribution of electron density function, assume EDFs are S type
!        if (allocated(b_EDF)) then
!            call EDFrho(3,x,y,z,EDFdens,EDFgrad,EDFhess)
!            laplx=laplx+EDFhess(1,1)
!            laply=laply+EDFhess(2,2)
!            laplz=laplz+EDFhess(3,3)
!        end if
!        if (label=='x') flapl=laplx
!        if (label=='y') flapl=laply
!        if (label=='z') flapl=laplz
!        if (label=='t') flapl=laplx+laply+laplz
!        flapl=flapl*laplfac !laplfac is an external variable
!        end function






!!!---- Calculate the DFT exchange-correlation functionals (two ways: 1-via libxc; 2-via internal routine)
!Note that the inner core density represented by EDF field is not taken into account
!
!!!!! EXCHANGE
!
real*8 function DFTxfunc(x,y,z)
real*8 x,y,z
if (wfntype==0.or.wfntype==3) then !Close-shell
	DFTxfunc=DFTxfunc_close(x,y,z)
else !Open-shell
	DFTxfunc=DFTxfunc_open(x,y,z)
end if
end function

real*8 function DFTxfunc_close(x,y,z)
real*8 x,y,z
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,abgrad,tgrad,atau,btau,ttau,alap,blap,tlap)
call rks_libxc_x(tdens,tgrad**2,ttau,DFTxfunc_close)
end function
real*8 function DFTxfunc_open(x,y,z)
real*8 x,y,z
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,abgrad,tgrad,atau,btau,ttau,alap,blap,tlap)
call uks_libxc_x(adens,bdens,agrad**2,bgrad**2,abgrad,atau,btau,DFTxfunc_open)
end function

subroutine rks_libxc_x(rho,sigma,tau,zk)
   use xc_f90_types_m
   use xc_f90_lib_m
   implicit real*8 (a-h,o-z)
   real*8 rho,sigma,lapl,tau,zk
   TYPE(xc_f90_pointer_t) :: xc_func
   TYPE(xc_f90_pointer_t) :: xc_info
   nx=xc_f90_functional_get_number(ksnamex)
   call xc_f90_func_init(xc_func, xc_info, nx, XC_UNPOLARIZED)
   select case (xc_f90_info_family(xc_info))
   case(XC_FAMILY_LDA)
	 call xc_f90_lda_exc(xc_func, 1, rho, ex)
   case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
	 call xc_f90_gga_exc(xc_func, 1, rho, sigma, ex)
   case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
	 call xc_f90_mgga_exc(xc_func, 1, rho, sigma, lapl, tau, ex)
   end select
   zk=rho*ex
   call xc_f90_func_end(xc_func)
end

subroutine uks_libxc_x(rhoa,rhob,sigmaaa,sigmabb,sigmaab,atau,btau,zk)
   use xc_f90_types_m
   use xc_f90_lib_m
   implicit real*8 (a-h,o-z)
   real*8 rhoa,rhob,sigmaaa,sigmabb,sigmaab,atau,btau,zk
   real*8 rho(2),sigma(3),tau(2),lapl(2)
   TYPE(xc_f90_pointer_t) :: xc_func
   TYPE(xc_f90_pointer_t) :: xc_info
   nx=xc_f90_functional_get_number(ksnamex)
   rho(1)=rhoa
   rho(2)=rhob
   sigma(1)=sigmaaa
   sigma(2)=sigmabb
   sigma(3)=sigmaab
   tau(1)=atau
   tau(2)=btau
   call xc_f90_func_init(xc_func, xc_info, nx, XC_POLARIZED)
   select case (xc_f90_info_family(xc_info))
   case(XC_FAMILY_LDA)
	 call xc_f90_lda_exc(xc_func, 1, rho(1), ex)
   case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
	 call xc_f90_gga_exc(xc_func, 1, rho(1),sigma(1), ex)
   case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
	 call xc_f90_mgga_exc(xc_func, 1, rho(1),sigma(1), lapl(1), tau(1), ex)
   end select
   zk=(rhoa+rhob)*ex
   call xc_f90_func_end(xc_func)
end
!
!!!!! CORRELATION
!
real*8 function DFTcfunc(x,y,z)
real*8 x,y,z
if (wfntype==0.or.wfntype==3) then !Close-shell
	DFTcfunc=DFTcfunc_close(x,y,z)
else !Open-shell
	DFTcfunc=DFTcfunc_open(x,y,z)
end if
end function

real*8 function DFTcfunc_close(x,y,z)
real*8 x,y,z
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,abgrad,tgrad,atau,btau,ttau,alap,blap,tlap)
call rks_libxc_c(tdens,tgrad**2,ttau,DFTcfunc_close)
end function

real*8 function DFTcfunc_open(x,y,z)
real*8 x,y,z
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,abgrad,tgrad,atau,btau,ttau,alap,blap,tlap)
call uks_libxc_c(adens,bdens,agrad**2,bgrad**2,abgrad,atau,btau,DFTcfunc_open)
end function


subroutine rks_libxc_c(rho,sigma,tau,zk)
   use xc_f90_types_m
   use xc_f90_lib_m
   implicit real*8 (a-h,o-z)
   real*8 rho,sigma,lapl,tau,zk
   TYPE(xc_f90_pointer_t) :: xc_func
   TYPE(xc_f90_pointer_t) :: xc_info
   nc=xc_f90_functional_get_number(ksnamec)	  
   call xc_f90_func_init(xc_func, xc_info, nc, XC_UNPOLARIZED)
   select case (xc_f90_info_family(xc_info))
   case(XC_FAMILY_LDA)
	 call xc_f90_lda_exc(xc_func, 1, rho, ec)
   case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
	 call xc_f90_gga_exc(xc_func, 1, rho, sigma, ec)
   case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
	 call xc_f90_mgga_exc(xc_func, 1, rho, sigma, lapl, tau, ec)
   end select
   zk=rho*ec
   call xc_f90_func_end(xc_func)
end

subroutine uks_libxc_c(rhoa,rhob,sigmaaa,sigmabb,sigmaab,atau,btau,zk)
   use xc_f90_types_m
   use xc_f90_lib_m
   implicit real*8 (a-h,o-z)
   real*8 rhoa,rhob,sigmaaa,sigmabb,sigmaab,atau,btau,zk
   real*8 rho(2),sigma(3),tau(2),lapl(2)
   TYPE(xc_f90_pointer_t) :: xc_func
   TYPE(xc_f90_pointer_t) :: xc_info
   nc=xc_f90_functional_get_number(ksnamec)	  
   rho(1)=rhoa
   rho(2)=rhob
   sigma(1)=sigmaaa
   sigma(2)=sigmabb
   sigma(3)=sigmaab
   tau(1)=atau
   tau(2)=btau
   call xc_f90_func_init(xc_func, xc_info, nc, XC_POLARIZED)
   select case (xc_f90_info_family(xc_info))
   case(XC_FAMILY_LDA)
	 call xc_f90_lda_exc(xc_func, 1, rho(1), ec)
   case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
	 call xc_f90_gga_exc(xc_func, 1, rho(1),sigma(1), ec)
   case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
	 call xc_f90_mgga_exc(xc_func, 1, rho(1),sigma(1), lapl(1), tau(1), ec)
   end select
   zk=(rhoa+rhob)*ec
   call xc_f90_func_end(xc_func)
end







!!!------------------------- Output gradient of rho and RDG(reduced density gradient) at a point
!label=x/y/z output 1-order derivation of x/y/z, =t get norm, =r get RDG
real*8 function fgrad(x,y,z,label)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),gradrho(3),EDFgrad(3),sumgrad2
character label
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
rho=0D0
gradrho=0D0
do i=1,nmo
    rho=rho+MOocc(i)*wfnval(i)**2
    gradrho(:)=gradrho(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
end do
gradrho=2*gradrho
! Add in contribution of Electron density function
if (allocated(b_EDF)) then
    call EDFrho(2,x,y,z,EDFdens,EDFgrad)
    rho=rho+EDFdens
    gradrho=gradrho+EDFgrad
end if
if (label=='x') fgrad=gradrho(1)
if (label=='y') fgrad=gradrho(2)
if (label=='z') fgrad=gradrho(3)
if (label=='t') fgrad=dsqrt( sum(gradrho(:)**2) )
if (label=='r') then
    sumgrad2=sum(gradrho(:)**2)
    if (RDG_maxrho/=0.0D0.and.rho>=RDG_maxrho) then
        fgrad=100D0
!This occurs at distant region when exponent cutoff is used, the actual value should be very large. In order to avoid denominator become zero, we set it artifically to a big value
    else if (sumgrad2==0D0.or.rho==0D0) then
        RDG=999D0
    else
        fgrad=0.161620459673995D0*dsqrt(sumgrad2)/rho**(4.0D0/3.0D0) !0.161620459673995D0=1/(2*(3*pi**2)**(1/3))
    end if
end if
end function

!!!------------------------- Output Laplacian of electron density at a point
!label=x/y/z output 2-order derivative of electron density respect to xx/yy/zz, =t get their summing
real*8 function flapl(x,y,z,label)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),laplx,laply,laplz,EDFgrad(3),EDFhess(3,3)
character label
call orbderv(3,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
laplx=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
laply=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
laplz=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
!Add in contribution of electron density function, assume EDFs are S type
if (allocated(b_EDF)) then
    call EDFrho(3,x,y,z,EDFdens,EDFgrad,EDFhess)
    laplx=laplx+EDFhess(1,1)
    laply=laply+EDFhess(2,2)
    laplz=laplz+EDFhess(3,3)
end if
if (label=='x') flapl=laplx
if (label=='y') flapl=laply
if (label=='z') flapl=laplz
if (label=='t') flapl=laplx+laply+laplz
flapl=flapl*laplfac !laplfac is an external variable
end function

!!!------------------------- Output spin or Alpha or Beta electron density at a point
!itype='s' output spin density, ='a' output alpha density, ='b' output beta density
real*8 function fspindens(x,y,z,itype)
real*8 :: x,y,z,wfnval(nmo)
character itype
call orbderv(1,1,nmo,x,y,z,wfnval)
adens=0.0D0
bdens=0.0D0
do i=1,nmo
    if (MOtype(i)==1) then
        adens=adens+MOocc(i)*wfnval(i)**2
    else if (MOtype(i)==2) then
        bdens=bdens+MOocc(i)*wfnval(i)**2
    else if (MOtype(i)==0) then
        adens=adens+MOocc(i)/2D0*wfnval(i)**2
        bdens=bdens+MOocc(i)/2D0*wfnval(i)**2
    end if
end do
if (itype=='s') then
    fspindens=adens-bdens
    if (ipolarpara==1) fspindens=fspindens/(adens+bdens)
else if (itype=='a') then
    fspindens=adens
else if (itype=='b') then
    fspindens=bdens
end if
end function


!!!-----Output ELF(Electron Localization Function) or LOL(Localized orbital locator) and similar function value at a point
!label="ELF"/"LOL"
real*8 function ELF_LOL(x,y,z,label)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo)
real*8 D,Dh,gradrho(3),gradrhoa(3),gradrhob(3),rho,rhoa,rhob,rhospin,MOoccnow
real*8 :: Fc=2.871234000D0 ! Fermi constant = (3/10)*(3*Pi^2)**(2/3) = 2.871234, 1/2.871234=0.34828
real*8 :: Fc_pol=4.557799872D0 ! Fermi constant for spin polarized = (3/10)*(6*Pi^2)**(2/3) = 4.5578, 1/4.5578=0.2194
character*3 label

!The ELF defined by Tsirelson, CPL,351,142.
!The LOL defined by Tsirelson, Acta Cryst. (2002). B58, 780-785.
!These functions are not important, so I don't write a special code
!for it, since rho, nebla-rho, nebla^2-rho support EDF, this function also support EDF
if (ELFLOL_type==1) then
    rho=fdens(x,y,z)
    if (wfntype==0.or.wfntype==3) then !close shell cases
        Dh=Fc*rho**(5.0D0/3.0D0)
    else if (wfntype==1.or.wfntype==2.or.wfntype==4) then !Open shell cases
        rhospin=fspindens(x,y,z,'s') !rhospin=rhoa-rhob, rho=rhoa+rhob
        rhoa=(rhospin+rho)/2D0
        rhob=(rho-rhospin)/2D0
        Dh=Fc_pol*(rhoa**(5.0D0/3.0D0)+rhob**(5.0D0/3.0D0)) !kinetic energy of HEG
    end if
    if (label=="ELF") then
        !Restrictly speaking, the kinetic energy expansion should be replace by polarized form for open-shell
        D=Dh-(1/9D0)*fgrad(x,y,z,'t')**2/rho+(1/6D0)*flapl(x,y,z,'t')
        ELF_LOL=1/(1+(D/Dh)**2)
    else if (label=="LOL") then
        D=Dh+(1/72D0)*fgrad(x,y,z,'t')**2/rho+(1/6D0)*flapl(x,y,z,'t')
        t=Dh/D
        ELF_LOL=1D0/(1D0/t+1)
    end if
    if (ELF_LOL<ELFLOL_cut) ELF_LOL=0
    return
end if

!For Becke's and Lu Tian's definition below
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)

D=0.0D0
rho=0.0D0
rhoa=0D0
rhob=0D0
gradrho=0D0
gradrhoa=0D0
gradrhob=0D0
if (label=="ELF") then
    if (wfntype==0.or.wfntype==3) then !spin-unpolarized case
        do i=1,nmo
            rho=rho+MOocc(i)*wfnval(i)**2
            gradrho(:)=gradrho(:)+2.0D0*MOocc(i)*wfnval(i)*wfnderv(:,i)
            D=D+MOocc(i)*(sum(wfnderv(:,i)**2)) !Calculate actual kinetic term
        end do        
        Dh=Fc*rho**(5.0D0/3.0D0) !Thomas-Fermi uniform electron gas kinetic energy
        if (rho==0D0) then
            ELF_LOL=0.0D0
            return
        end if
        D=D/2.0D0-(sum(gradrho(:)**2))/rho/8D0
    else if (wfntype==1.or.wfntype==2.or.wfntype==4) then !spin-polarized case
        do i=1,nmo
            MOoccnow=MOocc(i)
            if (MOtype(i)==0) MOoccnow=MOocc(i)/2D0 !double occupied, present when wfntype==2 (ROHF), alpha and beta get half part
            if (MOtype(i)==1.or.MOtype(i)==0) then
                rhoa=rhoa+MOoccnow*wfnval(i)**2
                gradrhoa(:)=gradrhoa(:)+2.0D0*MOoccnow*wfnval(i)*wfnderv(:,i)
            end if
            if (MOtype(i)==2.or.MOtype(i)==0) then
                rhob=rhob+MOoccnow*wfnval(i)**2
                gradrhob(:)=gradrhob(:)+2.0D0*MOoccnow*wfnval(i)*wfnderv(:,i)
            end if
            D=D+MOocc(i)*(sum(wfnderv(:,i)**2)) !Calculate actual kinetic term
        end do
        Dh=Fc_pol*(rhoa**(5.0D0/3.0D0)+rhob**(5.0D0/3.0D0))
        D=D/2.0D0
        if (rhoa/=0D0) D=D-(sum(gradrhoa(:)**2))/rhoa/8
        if (rhob/=0D0) D=D-(sum(gradrhob(:)**2))/rhob/8
    end if
    if (ELF_addminimal==1) D=D+1D-5 !add 1D-5 to avoid D become zero, leads to unity in infinite
    if (ELFLOL_type==0) ELF_LOL=1/(1+(D/Dh)**2)
    if (ELFLOL_type==2) ELF_LOL=1/(1+(D/Dh)) !New formalism defined by Tian Lu
    if (ELFLOL_type==995) ELF_LOL=Dh !Thomas-Fermi kinetic energy density
    if (ELFLOL_type==996) ELF_LOL=D/Dh !For test
    if (ELFLOL_type==997) ELF_LOL=D !For test
    if (ELFLOL_type==998) ELF_LOL=1/(1+D) !For test
    if (ELFLOL_type==999) ELF_LOL=1/(1+D**2) !For test
!----------------------------
else if (label=="LOL") then
    t=0.0D0
    if (wfntype==0.or.wfntype==3) then !Spin-unpolarized case
        do i=1,nmo !Store actual kinetic energy to t first
            rho=rho+MOocc(i)*wfnval(i)**2
            t=t+MOocc(i)*(sum(wfnderv(:,i)**2))
        end do
        t=t/2D0
        Dh=Fc*rho**(5.0D0/3.0D0)
    else if (wfntype==1.or.wfntype==2.or.wfntype==4) then !spin-polarized case
        do i=1,nmo !Store actual kinetic energy to t first
            MOoccnow=MOocc(i)
            if (MOtype(i)==0) MOoccnow=MOocc(i)/2D0 !double occupied, present when wfntype==2 (ROHF), alpha and beta get half part
            if (MOtype(i)==1.or.MOtype(i)==0) rhoa=rhoa+MOoccnow*wfnval(i)**2
            if (MOtype(i)==2.or.MOtype(i)==0) rhob=rhob+MOoccnow*wfnval(i)**2
            t=t+MOocc(i)*(sum(wfnderv(:,i)**2))
        end do
        t=t/2D0
        Dh=Fc_pol*(rhoa**(5.0D0/3.0D0)+rhob**(5.0D0/3.0D0))
    end if
!--------- A new definition of LOL, however the value range is not as good as LOL
! ELF_LOL=t-Dh
! if (ELF_LOL>0) then
!     ELF_LOL=1D0/(1D0+1D0/ELF_LOL)
! else if (ELF_LOL<0) then
!     ELF_LOL=-1D0/(1D0+1D0/abs(ELF_LOL))
! end if
!-------------
!If there is very long distance between molecule and current point, t (above) is zero,
!and t=Dh/t is also zero (because rho converges faster), but can't be calculate directly, so simply skip
    if (t/=0.0D0) t=Dh/t
    if (ELFLOL_type==0) ELF_LOL=1D0/(1D0/t+1) !namely t/(1+t). This is default case
    if (ELFLOL_type==2) ELF_LOL=1D0/((1D0/t)**2+1) !New form defined by Tian Lu
end if
if (ELF_LOL<ELFLOL_cut) ELF_LOL=0
end function


!!!----- Calculate ELF/LOL, its gradient and Hessian matrix at x,y,z, store to hess 
!!!!!!!!!!!! currently can not calculate Hessian
!funsel="ELF" or "LOL"
!itype=1 Only calculate value and gradient, not Hessian
!itype=2 Calculate value, gradient and Hessian (Not available)
subroutine calchessmat_ELF_LOL(itype,x,y,z,value,grad,hess,funsel)
use util
integer itype
real*8 x,y,z,value,grad(3),hess(3,3),MOoccnow
real*8 wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),wfntens3(3,3,3,nmo)
real*8 rho,gradrho(3),hessrho(3,3),Dh,gradDh(3),hessDh(3,3),Ts,gradTs(3),hessTs(3,3),Wei,gradWei(3),hessWei(3,3)
real*8 rhoa,rhob,gradrhoa(3),gradrhob(3),hessrhoa(3,3),hessrhob(3,3),Dha,Dhb,gradDha(3),gradDhb(3),Weia,Weib,gradWeia(3),gradWeib(3)
! real*8 hessDha(3,3),hessDhb(3,3),hessWeia(3,3),hessWeib(3,3)
real*8 :: Fc=2.871234000D0 ! Fermi constant = (3/10)*(3*Pi^2)**(2/3) = 2.871234, 1/2.871234=0.34828
real*8 :: Fc_pol=4.557799872D0 ! Fermi constant for spin polarized = (3/10)*(6*Pi^2)**(2/3) = 4.5578, 1/4.5578=0.2194
real*8 :: corrELF=1D-5
character funsel*3

if (itype==1) call orbderv(4,1,nmo,x,y,z,wfnval,wfnderv,wfnhess) !Get Hessian of GTF, needn't 3-order tensor
if (itype==2) call orbderv(5,1,nmo,x,y,z,wfnval,wfnderv,wfnhess,wfntens3)

!spin-unpolarized case
if (wfntype==0.or.wfntype==3) then
    rho=0D0
    gradrho=0D0
    Ts=0D0
    gradTs=0D0
    do i=1,nmo
        rho=rho+MOocc(i)*wfnval(i)**2
        gradrho(:)=gradrho(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
        Ts=Ts+MOocc(i)*(sum(wfnderv(:,i)**2))
    end do
    gradrho=2*gradrho
    Ts=Ts/2D0
    Dh=Fc*rho**(5D0/3D0)
    gradDh(:)=5D0/3D0*Fc*rho**(2D0/3D0)*gradrho(:)
    do i=1,nmo
        gradTs(1)=gradTs(1)+MOocc(i)*(wfnderv(1,i)*wfnhess(1,1,i)+wfnderv(2,i)*wfnhess(1,2,i)+wfnderv(3,i)*wfnhess(1,3,i))
        gradTs(2)=gradTs(2)+MOocc(i)*(wfnderv(2,i)*wfnhess(2,2,i)+wfnderv(1,i)*wfnhess(2,1,i)+wfnderv(3,i)*wfnhess(2,3,i))
        gradTs(3)=gradTs(3)+MOocc(i)*(wfnderv(3,i)*wfnhess(3,3,i)+wfnderv(1,i)*wfnhess(3,1,i)+wfnderv(2,i)*wfnhess(3,2,i))
    end do
    if (funsel=="ELF") then
        !Calculate Hessian for rho
        hessrho(1,1)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
        hessrho(2,2)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
        hessrho(3,3)=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
        hessrho(1,2)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(2,1:nmo)+wfnhess(1,2,1:nmo)*wfnval(1:nmo) ) )
        hessrho(2,3)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)*wfnderv(3,1:nmo)+wfnhess(2,3,1:nmo)*wfnval(1:nmo) ) )
        hessrho(1,3)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(3,1:nmo)+wfnhess(1,3,1:nmo)*wfnval(1:nmo) ) )
        hessrho(2,1)=hessrho(1,2)
        hessrho(3,1)=hessrho(1,3)
        hessrho(3,2)=hessrho(2,3)
        !Calculate Weizsacker functional and its derivatives
        Wei=sum(gradrho(:)**2)/8D0/rho
        D=Ts-Wei+corrELF
        chi=D/Dh
        value=1D0/(1D0+chi**2)
        gradWei(1)= 0.25D0/rho*( gradrho(1)*hessrho(1,1)+gradrho(2)*hessrho(1,2)+gradrho(3)*hessrho(1,3) ) - wei/rho*gradrho(1)
        gradWei(2)= 0.25D0/rho*( gradrho(2)*hessrho(2,2)+gradrho(1)*hessrho(2,1)+gradrho(3)*hessrho(2,3) ) - wei/rho*gradrho(2)
        gradWei(3)= 0.25D0/rho*( gradrho(3)*hessrho(3,3)+gradrho(2)*hessrho(3,2)+gradrho(1)*hessrho(3,1) ) - wei/rho*gradrho(3)
        chidx=(gradTs(1)-gradWei(1))/Dh - gradDh(1)/Dh**2 *D
        chidy=(gradTs(2)-gradWei(2))/Dh - gradDh(2)/Dh**2 *D
        chidz=(gradTs(3)-gradWei(3))/Dh - gradDh(3)/Dh**2 *D
        grad(1)=-2D0*chi/(1+chi**2)**2 * chidx
        grad(2)=-2D0*chi/(1+chi**2)**2 * chidy
        grad(3)=-2D0*chi/(1+chi**2)**2 * chidz
    else if (funsel=="LOL") then
        value=1D0/(1D0+Ts/Dh)
        tmp=-1D0/Dh/(1D0+Ts/Dh)**2
        grad(1)=tmp*(gradTs(1)-Ts/Dh*gradDh(1))
        grad(2)=tmp*(gradTs(2)-Ts/Dh*gradDh(2))
        grad(3)=tmp*(gradTs(3)-Ts/Dh*gradDh(3))
    end if
!spin-polarized case
else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
    rhoa=0D0
    rhob=0D0
    gradrhoa=0D0
    gradrhob=0D0
    Ts=0D0
    gradTs=0D0
    do i=1,nmo
        MOoccnow=MOocc(i)
        if (MOtype(i)==0) MOoccnow=MOocc(i)/2D0
        if (MOtype(i)==1.or.MOtype(i)==0) then
            rhoa=rhoa+MOoccnow*wfnval(i)**2
            gradrhoa(:)=gradrhoa(:)+MOoccnow*wfnval(i)*wfnderv(:,i)
        end if
        if (MOtype(i)==2.or.MOtype(i)==0) then
            rhob=rhob+MOoccnow*wfnval(i)**2
            gradrhob(:)=gradrhob(:)+MOoccnow*wfnval(i)*wfnderv(:,i)            
        end if
        Ts=Ts+MOocc(i)*(sum(wfnderv(:,i)**2))
    end do
    gradrhoa=2*gradrhoa
    gradrhob=2*gradrhob
    Ts=Ts/2D0
    Dha=Fc_pol*rhoa**(5D0/3D0)
    Dhb=Fc_pol*rhob**(5D0/3D0)
    Dh=Dha+Dhb
    gradDha(:)=5D0/3D0*Fc_pol*rhoa**(2D0/3D0)*gradrhoa(:)
    gradDhb(:)=5D0/3D0*Fc_pol*rhob**(2D0/3D0)*gradrhob(:)
    gradDh=gradDha+gradDhb
    do i=1,nmo
        gradTs(1)=gradTs(1)+MOocc(i)*(wfnderv(1,i)*wfnhess(1,1,i)+wfnderv(2,i)*wfnhess(1,2,i)+wfnderv(3,i)*wfnhess(1,3,i))
        gradTs(2)=gradTs(2)+MOocc(i)*(wfnderv(2,i)*wfnhess(2,2,i)+wfnderv(1,i)*wfnhess(2,1,i)+wfnderv(3,i)*wfnhess(2,3,i))
        gradTs(3)=gradTs(3)+MOocc(i)*(wfnderv(3,i)*wfnhess(3,3,i)+wfnderv(1,i)*wfnhess(3,1,i)+wfnderv(2,i)*wfnhess(3,2,i))
    end do
    
    if (funsel=="ELF") then
!         !Calculate Hessian for rho
        hessrhoa=0D0
        hessrhob=0D0
        do i=1,nmo
            MOoccnow=MOocc(i)
            if (MOtype(i)==0) MOoccnow=MOocc(i)/2D0 !double occupied, present when wfntype==2 (ROHF), alpha and beta get half part
            if (MOtype(i)==1.or.MOtype(i)==0) then
                hessrhoa(1,1)=hessrhoa(1,1)+MOoccnow*( wfnderv(1,i)**2 + wfnval(i)*wfnhess(1,1,i) )
                hessrhoa(2,2)=hessrhoa(2,2)+MOoccnow*( wfnderv(2,i)**2 + wfnval(i)*wfnhess(2,2,i) )
                hessrhoa(3,3)=hessrhoa(3,3)+MOoccnow*( wfnderv(3,i)**2 + wfnval(i)*wfnhess(3,3,i) )
                hessrhoa(1,2)=hessrhoa(1,2)+MOoccnow*( wfnderv(1,i)*wfnderv(2,i)+wfnhess(1,2,i)*wfnval(i) )
                hessrhoa(2,3)=hessrhoa(2,3)+MOoccnow*( wfnderv(2,i)*wfnderv(3,i)+wfnhess(2,3,i)*wfnval(i) )
                hessrhoa(1,3)=hessrhoa(1,3)+MOoccnow*( wfnderv(1,i)*wfnderv(3,i)+wfnhess(1,3,i)*wfnval(i) )
            end if
            if (MOtype(i)==2.or.MOtype(i)==0) then
                hessrhob(1,1)=hessrhob(1,1)+MOoccnow*( wfnderv(1,i)**2 + wfnval(i)*wfnhess(1,1,i) )
                hessrhob(2,2)=hessrhob(2,2)+MOoccnow*( wfnderv(2,i)**2 + wfnval(i)*wfnhess(2,2,i) )
                hessrhob(3,3)=hessrhob(3,3)+MOoccnow*( wfnderv(3,i)**2 + wfnval(i)*wfnhess(3,3,i) )
                hessrhob(1,2)=hessrhob(1,2)+MOoccnow*( wfnderv(1,i)*wfnderv(2,i)+wfnhess(1,2,i)*wfnval(i) )
                hessrhob(2,3)=hessrhob(2,3)+MOoccnow*( wfnderv(2,i)*wfnderv(3,i)+wfnhess(2,3,i)*wfnval(i) )
                hessrhob(1,3)=hessrhob(1,3)+MOoccnow*( wfnderv(1,i)*wfnderv(3,i)+wfnhess(1,3,i)*wfnval(i) )
            end if
        end do
        hessrhoa=hessrhoa*2
        hessrhob=hessrhob*2
        hessrhoa(2,1)=hessrhoa(1,2)
        hessrhoa(3,1)=hessrhoa(1,3)
        hessrhoa(3,2)=hessrhoa(2,3)
        hessrhob(2,1)=hessrhob(1,2)
        hessrhob(3,1)=hessrhob(1,3)
        hessrhob(3,2)=hessrhob(2,3)
!         !Calculate Weizsacker functional and its derivatives
        Weia=sum(gradrhoa(:)**2)/8D0/rhoa
        Weib=sum(gradrhob(:)**2)/8D0/rhob
        Wei=Weia+Weib
        D=Ts-Wei+corrELF
        chi=D/Dh
        value=1D0/(1D0+chi**2)
        gradWeia(1)= 0.25D0/rhoa*( gradrhoa(1)*hessrhoa(1,1)+gradrhoa(2)*hessrhoa(1,2)+gradrhoa(3)*hessrhoa(1,3) ) - weia/rhoa*gradrhoa(1)
        gradWeia(2)= 0.25D0/rhoa*( gradrhoa(2)*hessrhoa(2,2)+gradrhoa(1)*hessrhoa(2,1)+gradrhoa(3)*hessrhoa(2,3) ) - weia/rhoa*gradrhoa(2)
        gradWeia(3)= 0.25D0/rhoa*( gradrhoa(3)*hessrhoa(3,3)+gradrhoa(2)*hessrhoa(3,2)+gradrhoa(1)*hessrhoa(3,1) ) - weia/rhoa*gradrhoa(3)
        gradWeib(1)= 0.25D0/rhob*( gradrhob(1)*hessrhob(1,1)+gradrhob(2)*hessrhob(1,2)+gradrhob(3)*hessrhob(1,3) ) - weib/rhob*gradrhob(1)
        gradWeib(2)= 0.25D0/rhob*( gradrhob(2)*hessrhob(2,2)+gradrhob(1)*hessrhob(2,1)+gradrhob(3)*hessrhob(2,3) ) - weib/rhob*gradrhob(2)
        gradWeib(3)= 0.25D0/rhob*( gradrhob(3)*hessrhob(3,3)+gradrhob(2)*hessrhob(3,2)+gradrhob(1)*hessrhob(3,1) ) - weib/rhob*gradrhob(3)
        gradWei=gradWeia+gradWeib
        chidx=(gradTs(1)-gradWei(1))/Dh - gradDh(1)/Dh**2 *D
        chidy=(gradTs(2)-gradWei(2))/Dh - gradDh(2)/Dh**2 *D
        chidz=(gradTs(3)-gradWei(3))/Dh - gradDh(3)/Dh**2 *D
        grad(1)=-2D0*chi/(1+chi**2)**2 * chidx
        grad(2)=-2D0*chi/(1+chi**2)**2 * chidy
        grad(3)=-2D0*chi/(1+chi**2)**2 * chidz
    else if (funsel=="LOL") then
        value=1D0/(1D0+Ts/Dh)
        tmp=-1D0/Dh/(1D0+Ts/Dh)**2
        grad(1)=tmp*(gradTs(1)-Ts/Dh*gradDh(1))
        grad(2)=tmp*(gradTs(2)-Ts/Dh*gradDh(2))
        grad(3)=tmp*(gradTs(3)-Ts/Dh*gradDh(3))
    end if
end if

! if (itype==1) return
! ! Calculate Hessian for LOL, also need Hessian for rho
! hessrho(1,1)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
! hessrho(2,2)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
! hessrho(3,3)=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
! hessrho(1,2)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(2,1:nmo)+wfnhess(1,2,1:nmo)*wfnval(1:nmo) ) )
! hessrho(2,3)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)*wfnderv(3,1:nmo)+wfnhess(2,3,1:nmo)*wfnval(1:nmo) ) )
! hessrho(1,3)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(3,1:nmo)+wfnhess(1,3,1:nmo)*wfnval(1:nmo) ) )
! hessrho(2,1)=hessrho(1,2)
! hessrho(3,1)=hessrho(1,3)
! hessrho(3,2)=hessrho(2,3)
! 
! hessTs=0D0
! do i=1,nmo
!     hessTs(1,1)=hessTs(1,1)+MOocc(i)*( wfnhess(1,1,i)**2 + wfnderv(1,i)*wfntens3(1,1,1,i) + wfnhess(1,2,i)**2 + wfnderv(2,i)*wfntens3(1,1,2,i) + wfnhess(1,3,i)**2 + wfnderv(3,i)*wfntens3(1,1,3,i) )
!     hessTs(2,2)=hessTs(2,2)+MOocc(i)*( wfnhess(2,2,i)**2 + wfnderv(2,i)*wfntens3(2,2,2,i) + wfnhess(2,1,i)**2 + wfnderv(1,i)*wfntens3(2,2,1,i) + wfnhess(2,3,i)**2 + wfnderv(3,i)*wfntens3(2,2,3,i) )
!     hessTs(3,3)=hessTs(3,3)+MOocc(i)*( wfnhess(3,3,i)**2 + wfnderv(3,i)*wfntens3(3,3,3,i) + wfnhess(3,2,i)**2 + wfnderv(2,i)*wfntens3(3,3,2,i) + wfnhess(3,1,i)**2 + wfnderv(1,i)*wfntens3(3,3,1,i) )
!     hessTs(1,2)=hessTs(1,2)+MOocc(i)*( wfnhess(1,2,i)*wfnhess(1,1,i) + wfnderv(1,i)*wfntens3(1,1,2,i) + wfnhess(2,2,i)*wfnhess(1,2,i) + wfnderv(2,i)*wfntens3(1,2,2,i) + wfnhess(2,3,i)*wfnhess(1,3,i) + wfnderv(3,i)*wfntens3(1,2,3,i) )
!     hessTs(2,3)=hessTs(2,3)+MOocc(i)*( wfnhess(2,3,i)*wfnhess(2,2,i) + wfnderv(2,i)*wfntens3(2,2,3,i) + wfnhess(3,3,i)*wfnhess(2,3,i) + wfnderv(3,i)*wfntens3(2,3,3,i) + wfnhess(3,1,i)*wfnhess(2,1,i) + wfnderv(1,i)*wfntens3(2,3,1,i) )
!     hessTs(1,3)=hessTs(1,3)+MOocc(i)*( wfnhess(1,3,i)*wfnhess(1,1,i) + wfnderv(1,i)*wfntens3(1,1,3,i) + wfnhess(3,3,i)*wfnhess(1,3,i) + wfnderv(3,i)*wfntens3(1,3,3,i) + wfnhess(3,2,i)*wfnhess(1,2,i) + wfnderv(2,i)*wfntens3(1,3,2,i) )
! end do
! hessTs(2,1)=hessTs(1,2)
! hessTs(3,1)=hessTs(1,3)
! hessTs(3,2)=hessTs(2,3)
! 
! tmp1=10D0/9D0*Fc/rho**(1D0/3D0)
! tmp2=5D0/3D0*Fc*rho**(2D0/3D0)
! hessDh(1,1)=tmp1*gradrho(1)**2 + tmp2*hessrho(1,1)
! hessDh(2,2)=tmp1*gradrho(2)**2 + tmp2*hessrho(2,2)
! hessDh(3,3)=tmp1*gradrho(3)**2 + tmp2*hessrho(3,3)
! hessDh(1,2)=tmp1*gradrho(1)*gradrho(2) + tmp2*hessrho(1,2)
! hessDh(2,3)=tmp1*gradrho(2)*gradrho(3) + tmp2*hessrho(2,3)
! hessDh(1,3)=tmp1*gradrho(1)*gradrho(3) + tmp2*hessrho(1,3)
! hessDh(2,1)=hessDh(1,2)
! hessDh(3,1)=hessDh(1,3)
! hessDh(3,2)=hessDh(2,3)
! 
! !Diagonal of Hessian of LOL
! apre=1/Dh**2/(1+Ts/Dh)**3
! bpre=-1/Dh/(1+Ts/Dh)**2
! apartxx=(gradDh(1)+2*gradTs(1)-Ts/Dh*gradDh(1)) * (gradTs(1)-Ts/Dh*gradDh(1))
! bpartxx=hessTs(1,1)-( gradDh(1)*gradTs(1)-Ts/Dh*gradDh(1)**2+Ts*hessDh(1,1) )/Dh
! hess(1,1)=apre*apartxx+bpre*bpartxx
! apartyy=(gradDh(2)+2*gradTs(2)-Ts/Dh*gradDh(2)) * (gradTs(2)-Ts/Dh*gradDh(2))
! bpartyy=hessTs(2,2)-( gradDh(2)*gradTs(2)-Ts/Dh*gradDh(2)**2+Ts*hessDh(2,2) )/Dh
! hess(2,2)=apre*apartyy+bpre*bpartyy
! apartzz=(gradDh(3)+2*gradTs(3)-Ts/Dh*gradDh(3)) * (gradTs(3)-Ts/Dh*gradDh(3))
! bpartzz=hessTs(3,3)-( gradDh(3)*gradTs(3)-Ts/Dh*gradDh(3)**2+Ts*hessDh(3,3) )/Dh
! hess(3,3)=apre*apartzz+bpre*bpartzz
! !Non-diagonal of Hessian of LOL
! bpartxy=hessTs(1,2)-( gradDh(1)*gradTs(2)-Ts/Dh*gradDh(1)*gradDh(2)+Ts*hessDh(1,2) )/Dh
! hess(1,2)=apre*apartxx+bpre*bpartxy
! bpartyz=hessTs(2,3)-( gradDh(2)*gradTs(3)-Ts/Dh*gradDh(2)*gradDh(3)+Ts*hessDh(2,3) )/Dh
! hess(2,3)=apre*apartyy+bpre*bpartyz
! bpartxz=hessTs(1,3)-( gradDh(1)*gradTs(3)-Ts/Dh*gradDh(1)*gradDh(3)+Ts*hessDh(1,3) )/Dh
! hess(1,3)=apre*apartxx+bpre*bpartxz
! hess(2,1)=hess(1,2)
! hess(3,1)=hess(1,3)
! hess(3,2)=hess(2,3)
end subroutine









!!!--------------- Output Shannon information entropy function at a point
!itype=1 rho/N*ln(rho/N), this is normal definition
!itype=2 rho*ln(rho), this is Shannon information density, see J. Chem. Phys., 126, 191107
real*8 function infoentro(itype,x,y,z)
real*8 x,y,z,rho
integer itype
if (nelec==0D0) then
    infoentro=0D0
else
    rho=fdens(x,y,z)
    if (itype==1) rho=fdens(x,y,z)/nelec
    if (rho<=1D-100) then
        infoentro=0.0D0
    else
        infoentro=-rho*log(rho)
    end if
end if
end function
!!!----- Ghosh entropy density, PNAS, 81, 8028
!If itype==1, G(r) will be used as kinetic energy density
!If itype==2, G(r)-der2rho/8 will be used instead, which is the kinetic energy density exactly corresponding to Eq. 22 of PNAS, 81, 8028.
real*8 function Ghoshentro(x,y,z,itype)
integer itype
real*8 kintot,x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo)
if (itype==1) call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
if (itype==2) call orbderv(3,1,nmo,x,y,z,wfnval,wfnderv,wfnhess) !If K(r) is used, use this!!!
rho=0D0
do i=1,nmo
    rho=rho+MOocc(i)*wfnval(i)**2
end do
ck=2.871234D0
TFkin=ck*rho**(5D0/3D0)
kintot=0D0
!   If we use Hamiltonian kinetic density
! hamx=sum( MOocc(1:nmo)*wfnhess(1,1,1:nmo)*wfnval(1:nmo) )
! hamy=sum( MOocc(1:nmo)*wfnhess(2,2,1:nmo)*wfnval(1:nmo) )
! hamz=sum( MOocc(1:nmo)*wfnhess(3,3,1:nmo)*wfnval(1:nmo) )
! kintot=-(hamx+hamy+hamz)/2
!   If we use Lagrangian kinetic density G(r)
do i=1,nmo
    kintot=kintot+MOocc(i)*sum(wfnderv(:,i)**2)
end do
kintot=kintot/2D0
if (itype==2) then
    xlapl=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
    ylapl=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
    zlapl=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
    kintot=kintot-(xlapl+ylapl+zlapl)/8
end if
! if (kintot<0) then
!     write(*,"(5f16.10)") kintot,TFkin,x,y,z
!     read (*,*)
! end if
rlambda=5D0/3D0+log(4D0*pi*ck/3D0)
if (kintot<0) then
    rlogterm=0
else
    rlogterm=log(kintot/TFkin)
end if
Ghoshentro=1.5D0*rho*(rlambda+rlogterm)
end function






!!!------------------------- Output Lagrangian kinetic G(r) at a point. idir=0/1/2/3 means total/x/y/z kinetic energy
real*8 function Lagkin(x,y,z,idir)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo)
integer idir
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
lagkin=0D0
if (idir==0) then
    do imo=1,nmo
        lagkin=lagkin+MOocc(imo)*sum(wfnderv(:,imo)**2)
    end do
else
    do imo=1,nmo
        lagkin=lagkin+MOocc(imo)*wfnderv(idir,imo)**2
    end do
end if
lagkin=lagkin/2D0
end function





real*8 function sensfunc(x,y,z,u)
real*8 x,y,z,u,elf
real*8 param,augment,gam,tueg,s,uu
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,abgrad,tgrad,atau,btau,ttau,alap,blap,tlap)
!call uspline(tdens,tgrad,ttau,u,sensfunc)
sensfunc=1.0d0


elf=2D0/3D0*lagkin(x,y,z,0)/fdens(x,y,z)
!Ghoshentro(x,y,z,1)
!ELF_LOL(x,y,z,"LOL")


!write(*,*),'RPV here',s,u,elf


param=0.05d0

gam = 1d0
s = tgrad/(2*(3*pi)**(1/3)*tdens**(4/3))
uu = gam*s**2/(1.0d0+gam*s**2)
tueg=3/10*(3*pi**2)**(2/3)*tdens**(5/3)
w = (ttau/tueg-1)/(ttau/tueg+1)

augment=1.1d0
!if((uu.gt.(u-param)).and.(uu.lt.(u+param))) then
if(elf.lt.0.2) then
    sensfunc=augment
!    write(*,*),'RPV here',s,u,elf
else
    sensfunc=1.0d0
end if
end function

















subroutine uspline(dens,grad,tau,u_val,splval)
real*8 dens,grad,u,u_val,splval,cutoff,s,augment,param,w,tueg,tau



param=0.05d0

gam = 1d0
s = grad/(2*(3*pi)**(1/3)*dens**(4/3))
u = gam*s**2/(1.0d0+gam*s**2)
tueg=3/10*(3*pi**2)**(2/3)*rho**(5/3)
w = (tau/tueg-1)/(tau/tueg+1)



augment=1.1d0


if((u.gt.(u_val-param)).and.(u.lt.(u_val+param))) then
    splval=augment
!    if(u_val.eq.1.0) write(*,*),'RPV here',s,u,u_val,splval
else
    splval=1.0d0
end if
end






























!-RPV- These are the general boilerplate routines for implementing and developing new XC approximations (outside libxc)
real*8 function extDFTxfunc(x,y,z)
real*8 x,y,z
extDFTxfunc=0d0
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,abgrad,tgrad,atau,btau,ttau,alap,blap,tlap)
call sogga11x(extDFTxfunc,adens,bdens,agrad,bgrad,atau,btau)
end function

real*8 function extDFTcfunc(x,y,z)
real*8 x,y,z
extDFTcfunc=0d0
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,abgrad,tgrad,atau,btau,ttau,alap,blap,tlap)
call sogga11c(extDFTcfunc,adens,bdens,agrad,bgrad,atau,btau)
end function

real*8 function extOFDFTkfunc(x,y,z)
real*8 x,y,z
extOFDFTkfunc=0d0
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,abgrad,tgrad,atau,btau,ttau,alap,blap,tlap)
call PBEk(extOFDFTkfunc,adens,bdens,agrad,bgrad,atau,btau)
end function









! The following is kept in here for future usage with B13 (eventually)
!
!!!---- Calculate atomic density based on STO fitted or radial density
!!if indSTO==0, then all atom densities will be evaluated based on interpolation. if indSTO=18, then use STO fitted atomic density for element <18
!real*8 function calcatmdens(iatm,x,y,z,indSTO)
!use util
!real*8 rho,x,y,z,posarr(200),rhoarr(200)
!integer iatm,indSTO
!calcatmdens=0
!r=dsqrt( (a(iatm)%x-x)**2 + (a(iatm)%y-y)**2 + (a(iatm)%z-z)**2 )
!ind=a(iatm)%index
!! if (r>6*vdwr(ind)) return !Doesn't improve speed evidently but deteriorate result in rare cases
!if (ind<=indSTO) then !H~Ar, use STO fitted density. This is faster than using Lagrange interpolation technique, but not normalized to expected electron number
!	 do j=1,3
!		 if (YWTatomcoeff(ind,j)==0D0) cycle
!		 calcatmdens=calcatmdens+YWTatomcoeff(ind,j)*exp(-r/YWTatomexp(ind,j))
!	 end do
!else
!	 call genatmraddens(ind,posarr,rhoarr,npt) !Extract spherically averaged radial density of corresponding element
!	 call lagintpol(posarr(1:npt),rhoarr(1:npt),npt,r,calcatmdens,der1r,der2r,1)
!end if
!end function
!!!---- Calculate promolecular density purely based on interpolation of radial density, the promolecular density obtained in this manner is quite accurate
!real*8 function calcprodens(x,y,z,indSTO)
!real*8 x,y,z
!integer indSTO
!calcprodens=0
!do iatm=1,ncenter
!	 calcprodens=calcprodens+calcatmdens(iatm,x,y,z,indSTO)
!end do
!end function

!!--------Output electron density at a point
real*8 function fdens(x,y,z)
real*8 x,y,z,wfnval(nmo)
call orbderv(1,1,nmo,x,y,z,wfnval)
fdens=0D0
do i=1,nmo
	fdens=fdens+MOocc(i)*wfnval(i)**2
end do
! Add in contribution of Electron density function
if (allocated(b_EDF)) then
	call EDFrho(1,x,y,z,EDFdens)
	fdens=fdens+EDFdens
end if
! if (fdens>0.5D0) fdens=0
!  write(*,"(' ----RPV----	   :',f20.10)") fdens
end function

!Calculate the Hartree-potential (for the Coulomb term J)
real*8 function eleesp(Cx,Cy,Cz)
implicit none
integer,parameter :: narrmax=396 !Enough for h-type GTF, 5+5=10. Alri(0:10,0:5,0:5)-->11*6*6=396
real*8 Cx,Cy,Cz,term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval
real*8 sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,addesp
real*8 Alri(narrmax),Amsj(narrmax),Antk(narrmax),Fn(0:10) !Enough for h-type GTF, 5+5=10
real*8 twoepsqPC,tl,tm,tn,espprivate,espexpcut
integer nu,imo,iprim,jprim,maxFn,maplri(narrmax),mapmsj(narrmax),mapntk(narrmax),tmpnuml,tmpnumm,tmpnumn
integer Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi
eleesp=0.0D0
espexpcut=log10(espprecutoff)*3
nthreads=getNThreads()
!$OMP parallel do private(Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi, &
!$OMP nu,imo,iprim,jprim,maxFn,maplri,mapmsj,mapntk,tmpnuml,tmpnumm,tmpnumn,&
!$OMP twoepsqPC,tl,tm,tn,Alri,Amsj,Antk,Fn,&
!$OMP sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,addesp,&
!$OMP term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval,espprivate) shared(eleesp) schedule(dynamic) NUM_THREADS(nthreads)
do iprim=1,nprims
	espprivate=0D0
	icen=b(iprim)%center
	Aexp=b(iprim)%exp
	Ax=a(icen)%x
	Ay=a(icen)%y
	Az=a(icen)%z
	Aix=type2ix(b(iprim)%functype)
	Aiy=type2iy(b(iprim)%functype)
	Aiz=type2iz(b(iprim)%functype)
	sumAi=Aix+Aiy+Aiz
	do jprim=iprim,nprims
		jcen=b(jprim)%center
		Bexp=b(jprim)%exp
		Bix=type2ix(b(jprim)%functype)
		Biy=type2iy(b(jprim)%functype)
		Biz=type2iz(b(jprim)%functype)
		Bx=a(jcen)%x
		By=a(jcen)%y
		Bz=a(jcen)%z
		sumBi=Bix+Biy+Biz
		ep=Aexp+Bexp
		Px=(Ax*Aexp+Bx*Bexp)/ep
		Py=(Ay*Aexp+By*Bexp)/ep
		Pz=(Az*Aexp+Bz*Bexp)/ep
		PAx=Px-Ax
		PAy=Py-Ay
		PAz=Pz-Az
		PBx=Px-Bx
		PBy=Py-By
		PBz=Pz-Bz
		sqAB=(Ax-Bx)**2+(Ay-By)**2+(Az-Bz)**2
		PCx=Px-Cx
		PCy=Py-Cy
		PCz=Pz-Cz
		sqPC=PCx*PCx+PCy*PCy+PCz*PCz

		tmpval=-Aexp*Bexp*sqAB/ep
		prefac=2*pi/ep*dexp(tmpval)
		
		if (-ep*sqPC>espexpcut) then
			expngaPC=dexp(-ep*sqPC)
		else
			expngaPC=0
		end if
		maxFn=sumAi+sumBi
		Fn(maxFn)=Fmch(maxFn,ep*sqPC,expngaPC)
		nu=maxFn
		twoepsqPC=2*ep*sqPC
		do while (nu>0)
			Fn(nu-1)=(expngaPC+twoepsqPC*Fn(nu))/(2*nu-1) !cook book p280
			nu=nu-1
		end do

		tmpnuml=0
		do l=0,Aix+Bix
			tl=1.0D0
			if (mod(l,2)==1) tl=-1.0D0
			fjtmp=fj(l,Aix,Bix,PAx,PBx)*tl*fact(l)
			do r=0,l/2.0D0
				do i=0,(l-2*r)/2.0D0
					tmpnuml=tmpnuml+1
					Alri(tmpnuml)=Afac(l,r,i,PCx,ep,fjtmp)
					maplri(tmpnuml)=l-2*r-i
				end do
			end do
		end do

		tmpnumm=0
		do m=0,Aiy+Biy
			tm=1.0D0
			if (mod(m,2)==1) tm=-1.0D0
			fjtmp=fj(m,Aiy,Biy,PAy,PBy)*tm*fact(m)
			do s=0,m/2.0D0
				do j=0,(m-2*s)/2.0D0
					tmpnumm=tmpnumm+1
					Amsj(tmpnumm)=Afac(m,s,j,PCy,ep,fjtmp)
					mapmsj(tmpnumm)=m-2*s-j
				end do
			end do
		end do

		tmpnumn=0
		do n=0,Aiz+Biz
			tn=1.0D0
			if (mod(n,2)==1) tn=-1.0D0
			fjtmp=fj(n,Aiz,Biz,PAz,PBz)*tn*fact(n)
			do t=0,n/2.0D0
				do k=0,(n-2*t)/2.0D0
					tmpnumn=tmpnumn+1
					Antk(tmpnumn)=Afac(n,t,k,PCz,ep,fjtmp)
					mapntk(tmpnumn)=n-2*t-k
				end do
			end do
		end do

		term=0.0D0
		!Now calc "term"=<psi(iprim)|1/r_Z|psi(jprim)>
		do l=1,tmpnuml
			do m=1,tmpnumm
				do n=1,tmpnumn
					term=term+Alri(l)*Amsj(m)*Antk(n)*Fn(maplri(l)+mapmsj(m)+mapntk(n))
				end do
			end do
		end do

		if (iprim/=jprim) term=2.0*term
		term=term*prefac
		addesp=0.0D0
		do imo=1,nmo
			addesp=addesp+MOocc(imo)*CO(imo,iprim)*CO(imo,jprim)
		end do
		espprivate=espprivate+addesp*term
	end do !end j primitive
!$OMP critical
	eleesp=eleesp+espprivate
!$OMP end critical
end do !end i primitive
!$OMP end parallel do
end function

!!!------------------------- Calculate total ESP
real*8 function totesp(x,y,z)
real*8 x,y,z
totesp=eleesp(x,y,z)+nucesp(x,y,z)
end function

!!!------------------------- Calculate Coulomb (J)
real*8 function calcJ(x,y,z)
real*8 x,y,z
real*8 wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo)
calcJ=0.5D0*eleesp(x,y,z)*fdens(x,y,z)
end function


 
real*8 function calcK(Cx,Cy,Cz)
implicit none
integer,parameter :: narrmax=396 !Enough for h-type GTF, 5+5=10. Alri(0:10,0:5,0:5)-->11*6*6=396
real*8 Cx,Cy,Cz,Aij,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval,wfnval(nmo)
real*8 sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp
real*8 Alri(narrmax),Amsj(narrmax),Antk(narrmax),Fn(0:10) !Enough for h-type GTF, 5+5=10
real*8 twoepsqPC,tl,tm,tn,esppeve1,esppeve2,esppeve3,esppeve4,espexpcut,Pik,Pjl,func
real*8 chik,kcen,kexp,kexpterm,kix,kiy,kiz,krr,ksftx,ksftx2,ksfty,ksfty2,ksftz,ksftz2
real*8 chil,lcen,lexp,lexpterm,lix,liy,liz,lrr,lsftx,lsftx2,lsfty,lsfty2,lsftz,lsftz2
integer nu,imo,iprim,jprim,kprim,lprim,maxFn,maplri(narrmax),mapmsj(narrmax),mapntk(narrmax),tmpnuml,tmpnumm,tmpnumn
integer Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi
!calcK=0.0D0
esppeve1=0.0D0
esppeve2=0.0D0
esppeve3=0.0D0
esppeve4=0.0D0
espexpcut=log10(espprecutoff)*3
nthreads=getNThreads()
!$OMP parallel do private(Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi, &
!$OMP nu,imo,iprim,jprim,kprim,lprim,maxFn,maplri,mapmsj,mapntk,tmpnuml,tmpnumm,tmpnumn,&
!$OMP twoepsqPC,tl,tm,tn,Alri,Amsj,Antk,Fn,&
!$OMP sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,&
!$OMP Aij,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval) shared(calcK) schedule(dynamic) NUM_THREADS(nthreads)

!attempt at implementing K: this is far from being complete, and it might be wrong
do kprim=1,nprims
    chik=0.0d0
	kcen=b(kprim)%center
	kexp=b(kprim)%exp

	kix=type2ix(b(kprim)%functype)
	kiy=type2iy(b(kprim)%functype)
	kiz=type2iz(b(kprim)%functype)
	
	ksftx=Cx-a(b(kcen)%center)%x
	ksfty=Cy-a(b(kcen)%center)%y
	ksftz=Cz-a(b(kcen)%center)%z
	ksftx2=ksftx*ksftx
	ksfty2=ksfty*ksfty
	ksftz2=ksftz*ksftz
	krr=ksftx2+ksfty2+ksftz2
	if (expcutoff>0.or.-kexp*krr>expcutoff) then
		kexpterm=exp(-kexp*krr)
	else
		kexpterm=0D0
	end if
	if (kexpterm==0D0) cycle
	!Calculate value for current GTF
	if (b(kprim)%functype==1) then !Some functype use manually optimized formula for cutting do
	chik=kexpterm
	else if (b(kprim)%functype==2) then
	chik=ksftx*kexpterm
	else if (b(kprim)%functype==3) then
	chik=ksfty*kexpterm
	else if (b(kprim)%functype==4) then
	chik=ksftz*kexpterm
	else if (b(kprim)%functype==5) then
	chik=ksftx2*kexpterm
	else if (b(kprim)%functype==6) then
	chik=ksfty2*kexpterm
	else if (b(kprim)%functype==7) then
	chik=ksftz2*kexpterm
	else if (b(kprim)%functype==8) then
	chik=ksftx*ksfty*kexpterm
	else if (b(kprim)%functype==9) then
	chik=ksftx*ksftz*kexpterm
	else if (b(kprim)%functype==10) then
	chik=ksfty*ksftz*kexpterm
	else !If above condition is not satisfied(Angular moment higher than f), the function w
	chik=ksftx**kix *ksfty**kiy *ksftz**kiz *kexpterm
	end if

    do lprim=1,nprims
        chil=0.0d0
    	lcen=b(lprim)%center
    	lexp=b(lprim)%exp
    
    	lix=type2ix(b(lprim)%functype)
    	liy=type2iy(b(lprim)%functype)
    	liz=type2iz(b(lprim)%functype)
    	
    	lsftx=Cx-a(b(lcen)%center)%x
    	lsfty=Cy-a(b(lcen)%center)%y
    	lsftz=Cz-a(b(lcen)%center)%z
    	lsftx2=lsftx*lsftx
    	lsfty2=lsfty*lsfty
    	lsftz2=lsftz*lsftz
    	lrr=lsftx2+lsfty2+lsftz2
    	if (expcutoff>0.or.-lexp*lrr>expcutoff) then
    		lexpterm=exp(-lexp*lrr)
    	else
    		lexpterm=0D0
    	end if
    
    	if (lexpterm==0D0) cycle
    	
    	!Calculate value for current GTF
    	if (b(lprim)%functype==1) then !Some functype use manually optimized formula for cutting do
    	chil=lexpterm
    	else if (b(lprim)%functype==2) then
    	chil=lsftx*lexpterm
    	else if (b(lprim)%functype==3) then
    	chil=lsfty*lexpterm
    	else if (b(lprim)%functype==4) then
    	chil=lsftz*lexpterm
    	else if (b(lprim)%functype==5) then
    	chil=lsftx2*lexpterm
    	else if (b(lprim)%functype==6) then
    	chil=lsfty2*lexpterm
    	else if (b(lprim)%functype==7) then
    	chil=lsftz2*lexpterm
    	else if (b(lprim)%functype==8) then
    	chil=lsftx*lsfty*lexpterm
    	else if (b(lprim)%functype==9) then
    	chil=lsftx*lsftz*lexpterm
    	else if (b(lprim)%functype==10) then
    	chil=lsfty*lsftz*lexpterm
    	else !If above condition is not satisfied(Angular moment higher than f), the function w
    	chil=lsftx**lix *lsfty**liy *lsftz**liz *lexpterm
    	end if
    
        do iprim=1,nprims
        	icen=b(iprim)%center
        	Aexp=b(iprim)%exp
        	Ax=a(icen)%x
        	Ay=a(icen)%y
        	Az=a(icen)%z
        	Aix=type2ix(b(iprim)%functype)
        	Aiy=type2iy(b(iprim)%functype)
        	Aiz=type2iz(b(iprim)%functype)
        	sumAi=Aix+Aiy+Aiz
        	do jprim=iprim,nprims
        		jcen=b(jprim)%center
        		Bexp=b(jprim)%exp
        		Bix=type2ix(b(jprim)%functype)
        		Biy=type2iy(b(jprim)%functype)
        		Biz=type2iz(b(jprim)%functype)
        		Bx=a(jcen)%x
        		By=a(jcen)%y
        		Bz=a(jcen)%z
        		sumBi=Bix+Biy+Biz
        
        		ep=Aexp+Bexp
        		Px=(Ax*Aexp+Bx*Bexp)/ep
        		Py=(Ay*Aexp+By*Bexp)/ep
        		Pz=(Az*Aexp+Bz*Bexp)/ep
        		PAx=Px-Ax
        		PAy=Py-Ay
        		PAz=Pz-Az
        		PBx=Px-Bx
        		PBy=Py-By
        		PBz=Pz-Bz
        		sqAB=(Ax-Bx)**2+(Ay-By)**2+(Az-Bz)**2
        		PCx=Px-Cx
        		PCy=Py-Cy
        		PCz=Pz-Cz
        		sqPC=PCx*PCx+PCy*PCy+PCz*PCz
        
        		tmpval=-Aexp*Bexp*sqAB/ep
        		prefac=2*pi/ep*dexp(tmpval)
        		
        		if (-ep*sqPC>espexpcut) then
        			expngaPC=dexp(-ep*sqPC)
        		else
        			expngaPC=0
        		end if
        		maxFn=sumAi+sumBi
        		Fn(maxFn)=Fmch(maxFn,ep*sqPC,expngaPC)
        		nu=maxFn
        		twoepsqPC=2*ep*sqPC
        		do while (nu>0)
        			Fn(nu-1)=(expngaPC+twoepsqPC*Fn(nu))/(2*nu-1) !cook book p280
        			nu=nu-1
        		end do
        
        		tmpnuml=0
        		do l=0,Aix+Bix
        			tl=1.0D0
        			if (mod(l,2)==1) tl=-1.0D0
        			fjtmp=fj(l,Aix,Bix,PAx,PBx)*tl*fact(l)
        			do r=0,l/2.0D0
        				do i=0,(l-2*r)/2.0D0
        					tmpnuml=tmpnuml+1
        					Alri(tmpnuml)=Afac(l,r,i,PCx,ep,fjtmp)
        					maplri(tmpnuml)=l-2*r-i
        				end do
        			end do
        		end do
        
        		tmpnumm=0
        		do m=0,Aiy+Biy
        			tm=1.0D0
        			if (mod(m,2)==1) tm=-1.0D0
        			fjtmp=fj(m,Aiy,Biy,PAy,PBy)*tm*fact(m)
        			do s=0,m/2.0D0
        				do j=0,(m-2*s)/2.0D0
        					tmpnumm=tmpnumm+1
        					Amsj(tmpnumm)=Afac(m,s,j,PCy,ep,fjtmp)
        					mapmsj(tmpnumm)=m-2*s-j
        				end do
        			end do
        		end do
        
        		tmpnumn=0
        		do n=0,Aiz+Biz
        			tn=1.0D0
        			if (mod(n,2)==1) tn=-1.0D0
        			fjtmp=fj(n,Aiz,Biz,PAz,PBz)*tn*fact(n)
        			do t=0,n/2.0D0
        				do k=0,(n-2*t)/2.0D0
        					tmpnumn=tmpnumn+1
        					Antk(tmpnumn)=Afac(n,t,k,PCz,ep,fjtmp)
        					mapntk(tmpnumn)=n-2*t-k
        				end do
        			end do
        		end do
        
        		Aij=0.0D0
        		!Now calc "Aij"=<psi(iprim)|1/r_Z|psi(jprim)>
        		do l=1,tmpnuml
        			do m=1,tmpnumm
        				do n=1,tmpnumn
        					Aij=Aij+Alri(l)*Amsj(m)*Antk(n)*Fn(maplri(l)+mapmsj(m)+mapntk(n))
        				end do
        			end do
        		end do
		
		
		        if (iprim/=jprim) Aij=2.0*Aij
        		Aij=Aij*prefac
                Pik=0d0
                Pjl=0d0
        		do imo=1,nmo
        			Pik=Pik+MOocc(imo)*CO(imo,iprim)*CO(imo,kprim)
                    Pjl=Pjl+MOocc(imo)*CO(imo,jprim)*CO(imo,lprim)
        		end do        
        ! this is where we calculate the value of eq. 30 of Imamura (2007)		
        !		if (Pij*Pkl*chik*chil*Aij>espexpcut) then
                esppeve1=esppeve1-Pik*Pjl*chik*chil*Aij
        !        end if
        	end do !end j primitive
        	
!        	calcK=calcK+esppeve1
!        	esppeve2=esppeve2+esppeve1	
!        	esppeve1=0d0
        end do !end i primitive

        calcK=calcK+esppeve1
!        	esppeve3=esppeve3+esppeve2
!        	esppeve2=0d0
        end do !end k primitive
!        esppeve4=esppeve4+esppeve3
!        esppeve3=0d0



!  write(*,"(' ----RPV----	   :',f20.5,f10.5,f10.5,f10.5,f10.5)") calcK
end do !end l primitive
!$OMP end parallel do
!calcK=calcK+esppeve1
!esppeve1=0d0	
calcK=calcK*fdens(Cx,Cy,Cz)*0.125D0

end function



















!!!------------------------- Calculate A-factor at Cook book p245 ---page 227 actually
real*8 function Afac(l,r,i,PC,gamma,fjtmp)
integer l,r,i,ti
real*8 gamma,PC,comp,PCterm,fjtmp
ti=1.0D0
if (mod(i,2)==1) ti=-1.0D0 !faster than ti=(-1)**i
PCterm=1.0D0
if (l-2*r-2*i/=0) PCterm=PC**(l-2*r-2*i)
comp=ti*PCterm*(0.25D0/gamma)**(r+i) / ( fact(r)*fact(i)*fact(l-2*r-2*i) )
Afac=fjtmp*comp
! comp=(-1)**i*fact(l)*PC**(l-2*r-2*i)*(1/(4*gamma))**(r+i) / ( fact(r)*fact(i)*fact(l-2*r-2*i) )
! Afac=(-1)**l * fj(l,l1,l2,PA,PB)*comp
end function

!!!------------------------- Calculate fj at Cook book p237	 ---page 219 actually
real*8 function fj(j,l,m,aa,bb)
real*8 aa,bb,pre,at,bt
integer j,l,m,k,imin,imax
imax=min(j,l)
imin=max(0,j-m)
fj=0.0D0
do k=imin,imax
	pre=fact(l)/fact(l-k)/fact(k) * fact(m)/fact(m-j+k)/fact(j-k)
	at=1.0D0
	bt=1.0D0
	if (l-k/=0) at=aa**(l-k)  !This determine helps to improve efficient
	if (m+k-j/=0) bt=bb**(m+k-j)
	fj=fj+pre*at*bt
end do
end function

!!!---- Calculate int('t^(2*m)*exp(-x*t^2)','t',0,1) see Cook book p281 for detail ---page 260, actually
real*8 function Fmch(m,x,expnx)
! expnx is input parameter, value should be exp(-x), because calculate the value is time-consuming and in
! other place this value also need be calculate, so not recalculate in this subroutine
IMPLICIT none
integer m,i
real*8 x,expnx,a,b,term,partsum,APPROX,xd,FIMULT,NOTRMS,eps,fiprop
eps=1.0D-8	!convergence precision
Fmch=0D0
if (x<=10) then
	if (expnx==0D0) RETURN
	A=m+0.5D0
	term=1.0D0/A
	partsum=term
	DO I=2,50
		A=A+1.0D0
		term=term*X/A
		partsum=partsum+term
		if ( term/partsum < eps) THEN
		   Fmch = 0.5D0*partsum*expnx
		   RETURN
		END IF
	END DO
	write(*,*) "Error: Fmch didn't converge"
else !x is big, use suitable method for solve this situation
	A=M
	B=A+0.5D0
	A=A-0.5D0
	XD=1.0D0/X
	APPROX=0.88622692D0*(dsqrt(XD)*XD**m)
	DO I=1,m
		B=B-1.0D0
		APPROX=APPROX*B
	END DO
	FIMULT=0.5D0*expnx*XD
	partsum=0.D0
	IF (FIMULT==0.0D0) THEN
		Fmch=APPROX-FIMULT*partsum
		return
	ELSE
		FIPROP=FIMULT/APPROX
		term=1.0d0
		partsum=term
		NOTRMS=X
		NOTRMS=NOTRMS+M
		DO I=2,NOTRMS
		   term=term*A*XD
		   partsum=partsum+term
		   IF (dabs(term*FIPROP/partsum)<eps)  THEN
			  Fmch=APPROX-FIMULT*partsum
			  RETURN
		   END IF
		   A=A-1.0D0
		END DO
		write(*,*) "Didn't converge"
	END IF
end if
end function

!!!------------- Output Exchange-correlation density, correlation hole and correlation factor
!rfx,rfy,rfz is reference point (commonly use refx,refy,refz in defvar module), namely r1
!x,y,z in the argument is the coordinate of r2
!Calculate which function is controlled by "pairfunctype" in settings.ini, correlation type is determined by "paircorrtype" in settings.ini

!    pairfunctype==1  Correlation hole for alpha
!    pairfunctype==2  Correlation hole for beta
!    pairfunctype==4  Correlation factor for alpha
!    pairfunctype==5  Correlation factor for beta
!    pairfunctype==7  Exc.-corr. density for alpha
!    pairfunctype==8  Exc.-corr. density for beta
!    pairfunctype==10 Pair density for alpha
!    pairfunctype==11 Pair density for beta
!    pairfunctype==12 Pair density for all electrons

real*8 function pairfunc(rfx,rfy,rfz,x,y,z,pairfunctype)
real*8 rfx,rfy,rfz,x,y,z,orbvalr1(nmo),orbvalr2(nmo)
integer pairfunctype
call orbderv(1,1,nmo,rfx,rfy,rfz,orbvalr1)
call orbderv(1,1,nmo,x,y,z,orbvalr2)
!Calculate alpha and beta density at r1 and r2
adensr1=0.0D0
bdensr1=0.0D0
adensr2=0.0D0
bdensr2=0.0D0
do i=1,nmo
    if (MOtype(i)==0) then
        adensr1=adensr1+MOocc(i)/2*orbvalr1(i)**2
        adensr2=adensr2+MOocc(i)/2*orbvalr2(i)**2
        bdensr1=bdensr1+MOocc(i)/2*orbvalr1(i)**2
        bdensr2=bdensr2+MOocc(i)/2*orbvalr2(i)**2
    else if (MOtype(i)==1) then
        adensr1=adensr1+MOocc(i)*orbvalr1(i)**2
        adensr2=adensr2+MOocc(i)*orbvalr2(i)**2
    else if (MOtype(i)==2) then
        bdensr1=bdensr1+MOocc(i)*orbvalr1(i)**2
        bdensr2=bdensr2+MOocc(i)*orbvalr2(i)**2
    end if
end do
totdensr1=adensr1+bdensr1
totdensr2=adensr2+bdensr2

ntime=1
if (pairfunctype==12) ntime=2 !Will need both alpha and beta information (aXCdens and bXCdens), so process twice
do itime=1,ntime
    !Calculate exchange-correlation density first, and then calculate correlation hole and correlation factor
    !For RHF/ROHF, we calculate them as if present system is open-shell
    if (pairfunctype==1.or.pairfunctype==4.or.pairfunctype==7.or.pairfunctype==10.or.(pairfunctype==12.and.itime==1)) then !Set start and end index of alpha orbitals
        !Cycle alpha orbitals first to obtain aXCdens
        istart=1
        if (wfntype==0.or.wfntype==2.or.wfntype==3) then !RHF,ROHF,R-post-HF
            iend=nmo
        else if (wfntype==1.or.wfntype==4) then !UHF, U-post-HF
            do iend=nmo,1,-1
                if (MOtype(iend)==1) exit
            end do
        end if
    else if (pairfunctype==2.or.pairfunctype==5.or.pairfunctype==8.or.pairfunctype==11.or.(pairfunctype==12.and.itime==2)) then !Set start and end index of beta orbitals
        if (wfntype==0.or.wfntype==3) then !RHF,R-post-HF
            istart=1
            iend=nmo
        else if (wfntype==2) then !ROHF
            istart=1
            do iend=1,nmo
                if (MOtype(iend)==1) exit
            end do
            iend=iend-1
        else if (wfntype==1.or.wfntype==4) then !UHF, U-post-HF
            do istart=1,nmo
                if (MOtype(istart)==2) exit
            end do
            iend=nmo
            if (nint(nbelec)==0) iend=0 !less than istart, so below cycle will be skipped
        end if
    end if

    XCtmp=0D0 !Really X+C density
    Xtmp=0D0 !Only X density
    Ctmp=0D0 !Only C density
    do i=istart,iend
        occi=MOocc(i)
        if (MOtype(i)==0) occi=occi/2D0 !Split close-shell orbital to spin orbital
        do j=istart,iend
            occj=MOocc(j)
            if (MOtype(j)==0) occj=occj/2D0
            tmpmul=orbvalr1(i)*orbvalr2(j)*orbvalr1(j)*orbvalr2(i)
            XCtmp=XCtmp-dsqrt(occi*occj)*tmpmul
            Xtmp=Xtmp-occi*occj*tmpmul
            Ctmp=Ctmp+(occi*occj-dsqrt(occi*occj))*tmpmul
        end do
    end do
        
    if (pairfunctype==1.or.pairfunctype==4.or.pairfunctype==7.or.pairfunctype==10.or.(pairfunctype==12.and.itime==1)) then
        if (paircorrtype==1) aXCdens=Xtmp
        if (paircorrtype==2) aXCdens=Ctmp
        if (paircorrtype==3) aXCdens=XCtmp
        acorrhole=aXCdens/adensr1
        acorrfac=acorrhole/adensr2
        if (pairfunctype==1) pairfunc=acorrhole
        if (pairfunctype==4) pairfunc=acorrfac
        if (pairfunctype==7) pairfunc=aXCdens
        if (pairfunctype==10) pairfunc=adensr1*adensr2+aXCdens
    else if (pairfunctype==2.or.pairfunctype==5.or.pairfunctype==8.or.pairfunctype==11.or.(pairfunctype==12.and.itime==2)) then
        if (paircorrtype==1) bXCdens=Xtmp
        if (paircorrtype==2) bXCdens=Ctmp
        if (paircorrtype==3) bXCdens=XCtmp
        bcorrhole=bXCdens/bdensr1
        bcorrfac=bcorrhole/bdensr2
        if (pairfunctype==2) pairfunc=bcorrhole
        if (pairfunctype==5) pairfunc=bcorrfac
        if (pairfunctype==8) pairfunc=bXCdens
        if (pairfunctype==11) pairfunc=bdensr1*bdensr2+bXCdens
    end if
end do
if (pairfunctype==12) pairfunc=adensr1*(adensr2+bdensr2)+aXCdens +bdensr1*(adensr2+bdensr2)+bXCdens
end function




real*8 function funcJ(x1,y1,z1,x2,y2,z2)
real*8 x1,y1,z1,x2,y2,z2,orbvalr1(nmo),orbvalr2(nmo)
funcJ=0.0d0
or12=0d0
r12=dsqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
if (r12/=0) or12=1.d0/r12
call orbderv(1,1,nmo,x1,y1,z1,orbvalr1)
call orbderv(1,1,nmo,x2,y2,z2,orbvalr2)
do i=1,nmo
    occi=MOocc(i)
    if (MOtype(i)==0) occi=occi/2D0 !Split close-shell orbital to spin orbital
    do j=1,nmo
        occj=MOocc(j)
        if (MOtype(j)==0) occj=occj/2D0
          funcJ=funcJ+2.d0*occi*occj*orbvalr1(i)*orbvalr1(i)*or12*orbvalr2(j)*orbvalr2(j)
     end do
 end do
end function



real*8 function funcK(x1,y1,z1,x2,y2,z2)
real*8 x1,y1,z1,x2,y2,z2,orbvalr1(nmo),orbvalr2(nmo)
funcK=0.0d0
or12=0d0
r12=dsqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
if (r12/=0) or12=1.d0/r12
call orbderv(1,1,nmo,x1,y1,z1,orbvalr1)
call orbderv(1,1,nmo,x2,y2,z2,orbvalr2)
do i=1,nmo
    occi=MOocc(i)
    if (MOtype(i)==0) occi=occi/2D0 !Split close-shell orbital to spin orbital
    do j=1,nmo
        occj=MOocc(j)
        if (MOtype(j)==0) occj=occj/2D0
          funcK=funcK-occi*occj*orbvalr1(i)*orbvalr1(j)*or12*orbvalr2(j)*orbvalr2(i)
     end do
 end do
end function


real*8 function func3p5(x1,y1,z1,x2,y2,z2)
real*8 x1,y1,z1,x2,y2,z2,orbvalr1(nmo),orbvalr2(nmo)
func3p5=0.0d0
or12=0d0
r12=dsqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
if (r12/=0) or12=1.d0/r12
call orbderv(1,1,nmo,x1,y1,z1,orbvalr1)
call orbderv(1,1,nmo,x2,y2,z2,orbvalr2)
do i=1,nmo
    occi=MOocc(i)
    if (MOtype(i)==0) occi=occi/2D0 !Split close-shell orbital to spin orbital
    do j=1,nmo
        occj=MOocc(j)
        if (MOtype(j)==0) occj=occj/2D0
          func3p5=func3p5-occi*occj*orbvalr1(i)*orbvalr1(j)*or12*orbvalr2(j)*orbvalr2(i)
     end do
 end do
end function






end module
