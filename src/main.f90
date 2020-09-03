program dftiog
use defvar
use util
use function
implicit real*8(a-h,o-z)
character nowdate*20,nowtime*20,inpstring*80,c200tmp*200,c200tmp2*200,c2000tmp*2000,outcubfile*200,selectyn,lovername*80,settingpath*200
real*8 :: inx,iny,inz,tmpvec(3),intgv,edisp,egcp,ebas
integer :: iprintfunc=1 !The default function whose gradient and Hessian will be outputted at a point by main function 1
integer,allocatable :: tmparrint(:)
integer walltime1,walltime2
real*8,allocatable :: d1add(:,:),d1min(:,:),d2add(:,:),d2min(:,:),d1addtmp(:,:),d1mintmp(:,:),d2addtmp(:,:),d2mintmp(:,:) !Store temporary data for drawing gradient map
real*8,allocatable :: planemat_cust(:,:) !For storing temporary data of doing custom map
real*8,allocatable :: planemat_bk(:,:) !Used to backup plane data
real*8,allocatable :: tmpmat(:,:),tmparr(:)
real*8 ,allocatable ::xyz (:,:)

!tst
integer, external   :: integrand, integrandJ, integrandK
real*8              :: integral(1), error(1), prob(1)
integer             :: nregions, fail
integer*8           :: neval
!endtst

call getarg(1,filename)
call walltime(iwalltime1)
CALL CPU_TIME(time_begin)
if (trim(filename)=="") then !Haven't defined filename variable
	inquire(file=filename,exist=alive)
		write(*,*) "No input file specified, exit program..."
		read(*,*)
		stop
end if

!call loadinput(filename)
!
!http://www.patorjk.com/software/taag/#p=display&f=Tiles&t=DFT%3AMCI
write(*,*) 
write(*,*) 
write(*,*) "     [.....    [........[... [......   [..       [..    [..   [.. "
write(*,*) "     [..   [.. [..           [..       [. [..   [... [..   [..[.. "
write(*,*) "     [..    [..[..           [..       [.. [.. [ [..[..       [.. "
write(*,*) "     [..    [..[......       [..       [..  [..  [..[..       [.. "
write(*,*) "     [..    [..[..           [..    [..[..   [.  [..[..       [.. "
write(*,*) "     [..   [.. [..           [..       [..       [.. [..   [..[.. "
write(*,*) "     [.....    [..           [..    [..[..       [..   [....  [.. "
write(*,*)                                                             
write(*,*) "       +======================================================+   "
write(*,*) "       |  Density Functional Theory: Monte Carlo Integration  |   "
write(*,*) "       |  Version 0.2 (2019.alpha)                            |   "
write(*,*) "       |  RPV Research                                        |   "
write(*,*) "       |  Peverati Lab at Florida Tech                        |   "
write(*,*) "       +------------------------------------------------------+   "
call date_and_time(nowdate,nowtime)
write(*,"('        |       Current date: ',a,'-',a,'-',a,'  Time: ',a,':',a,':',a,'       |')") &
		nowdate(1:4),nowdate(5:6),nowdate(7:8),nowtime(1:2),nowtime(3:4),nowtime(5:6)
write(*,*) "       +======================================================+   "
write(*,*)
call readinfile(filename,0)
!!-- Backup various information of first loaded (meanwhile unmodified) molecule
allocate(a_org(ncenter))
allocate(b_org(nprims))
allocate(CO_org(nmo,nprims))
allocate(MOocc_org(nmo))
a_org=a
b_org=b
CO_org=CO
MOocc_org=MOocc
nprims_org=nprims
nmo_org=nmo
ncenter_org=ncenter

!!-- Initialize fragment
nfragatmnum=ncenter !Default fragment is the whole molecule
nfragatmnumbackup=ncenter
allocate(fragatm(nfragatmnum),fragatmbackup(nfragatmnum))
forall (i=1:nfragatmnum) fragatm(i)=i
forall (i=1:nfragatmnum) fragatmbackup(i)=i
ifragcontri=0

!!-- Call some routines only once
if (ncenter>5000) write(*,"(a)") " Warning: There are a lot of atoms in your system. Generating the distance matrix. It might take a while..."
call gendistmat !Generate distance matrix

!-- Show related molecular information
if (ifiletype/=0.and.ifiletype/=8) then
	call showformula
	totmass=sum(atmwei(a%index))
	write(*,"(' Molecular weight:',f16.5)") totmass
end if
write(*,*)


if (doINT.eq.2 .or. doINT.eq.3) then
    write(*,"('---=< calling MC integration for xc >=--- ')")
    ndim=3
    ncomp=1
    if (MClvl.eq.1) then
        write(*,"(' Using Vegas algorithm')")
    	call vegas(ndim, ncomp, integrand, userdata,  &
         epsrel, epsabs, flags, seed,                 &
         mineval, maxeval, nstart, nincrease, nbatch, &
         gridno, "",                           &
         neval, fail, integral, error, prob)        
    else if (MClvl.eq.2) then
        write(*,"(' Using Suave algorithm')")
        call suave(ndim, ncomp, integrand, userdata,               & 
                   epsrel, epsabs, flags, seed, mineval, maxeval,  & 
                   nnew, flatness, "",                             & 
                   nregions, neval, fail, integral, error, prob)                    
    else if (MClvl.eq.3) then
        write(*,"(' Using Divonne algorithm')")
        call lldivonne(ndim, ncomp, integrand, userdata, & 
          epsrel, epsabs, flags, seed, mineval, maxeval, & 
          key1, key2, key3, maxpass,                     & 
          border, maxchisq, mindeviation,                & 
          ngiven, ldxgiven, 0., nextra, 0.,              & 
          "",                                            & 
          nregions, neval, fail, integral, error, prob)    
    else if (MClvl.eq.4) then      
        write(*,"(' Using Cuhre algorithm')")          
        call cuhre(ndim, ncomp, integrand, userdata, &    
                   epsrel, epsabs, flags,            &    
                   mineval, maxeval, key,            &    
                   "",                               &    
                   nregions, neval, fail, integral, error, prob)                
      end if
                    
      write(*,*) ' Results OF Monte Carlo integration:'
      write(*,"('  MCprc    =          ',ES12.0)") epsabs
      write(*,*) ' nregions =          ', nregions
      write(*,*) ' neval    = ',neval
      write(*,*) ' fail     =          ',     fail
      if(fail.gt.0)write(*,*) ' *** INTEGRATION FAILED TO CONVERGE *** '
      write(*,*)
      write(* ,30) integral(1), error, prob
      30 format(' Exc=',f17.8,'   error=',f14.8,'   p=',f10.3)
      write(*,*)
    CALL CPU_TIME(time_MC)
    call walltime(iwalltimeINT)
    write(*,"('---=< MC done. Total job time: CPU time',f12.2,'s, wall clock time',i10,'s >=--- ',/)") time_MC-time_begin,iwalltimeINT-iwalltime1
end if


!! Integrating J and K
if (doJMC.eq.1) then
    write(*,"('---=< calling MC integration for J >=--- ')")
    ndim2=6
    ncomp=1
    write(*,"(' Using Divonne algorithm')")
    call lldivonne(ndim2, ncomp, integrandJ, userdata, & 
      epsrel, epsabs, flags, seed, mineval, maxeval, & 
      key1, key2, key3, maxpass,                     & 
      border, maxchisq, mindeviation,                & 
      ngiven, ldxgiven, 0., nextra, 0.,              & 
      "",                                            & 
      nregions, neval, fail, integral, error, prob)    
    write(*,*) ' Results OF Monte Carlo integration:'
    write(*,"('  MCprc    =          ',ES12.0)") epsabs
    write(*,*) ' nregions =          ', nregions
    write(*,*) ' neval    = ',neval
    write(*,*) ' fail     =          ',     fail
    if(fail.gt.0)write(*,*) ' *** INTEGRATION FAILED TO CONVERGE *** '
    write(*,*)
    write(* ,10) integral(1), error, prob
    10 format(' Coulomb=',f17.8,'   error=',f14.8,'   p=',f10.3)
    write(*,*)
    CALL CPU_TIME(time_MCK)
    call walltime(iwalltimeINTK)
    write(*,"('---=< J done. Total job time: CPU time',f12.2,'s, wall clock time',i10,'s >=--- ',/)") time_MCK-time_MC,iwalltimeINTK-iwalltimeINT
end if

if (doK.eq.1) then
    write(*,"('---=< calling MC integration for K >=--- ')")
    ndim2=6
    ncomp=1
    write(*,"(' Using Divonne algorithm')")
    call lldivonne(ndim2, ncomp, integrandK, userdata, & 
      epsrel, epsabs, flags, seed, mineval, maxeval, & 
      key1, key2, key3, maxpass,                     & 
      border, maxchisq, mindeviation,                & 
      ngiven, ldxgiven, 0., nextra, 0.,              & 
      "",                                            & 
      nregions, neval, fail, integral, error, prob)    

    write(*,*) ' Results OF Monte Carlo integration:'
    write(*,"('  MCprc    =          ',ES12.0)") epsabs
    write(*,*) ' nregions =          ', nregions
    write(*,*) ' neval    = ',neval
    write(*,*) ' fail     =          ',     fail
    if(fail.gt.0)write(*,*) ' *** INTEGRATION FAILED TO CONVERGE *** '
    write(*,*)
    write(* ,20) integral(1), error, prob
    20 format(' HF Exchange=',f17.8,'   error=',f14.8,'   p=',f10.3)
    write(*,*)
    CALL CPU_TIME(time_MCK)
    call walltime(iwalltimeINTK)
    write(*,"('---=< K done. Total job time: CPU time',f12.2,'s, wall clock time',i10,'s >=--- ',/)") time_MCK-time_MC,iwalltimeINTK-iwalltimeINT
end if


!!--Call main routine to start the integration on the grid
if (doINT.eq.1 .or. doINT.eq.3) then
    write(*,"('---=< calling grid integration for xc >=--- ')")
    call intfunc(intgv)
!    call sensan(intgv)
end if
!call sys1eprop !Show some system 1e properties, it works with Cartesian basis functions only

!!--Call to ppot for DFT-D and other Grimme corrections
if (doGrimme.ne.0) then
	write(*,"('---=< entering ppot program for Grimme corrections >=--- ')")
	call ppot(edisp,egcp,ebas)
end if

if ((doJ.ne.0).and.(doGrimme.ne.0)) then
	write(*,"('---=< summing up some final total energies >=---')")
	write(*,"(' Electronic+nuclear energy (from IoG):',f20.10)") intgv
	write(*,"(' -D corrections                      :',f20.10)") edisp
	write(*,"(' gCP corrections                     :',f20.10)") egcp
	write(*,"(' Basis set incompleteness corrections:',f20.10)") ebas
	write(*,"('                                      --------------------')")
	write(*,"('   +--> Total ENERGY with corrections:',f20.10' <--+',/)") intgv+edisp+egcp+ebas
end if

!!--Call to DFT-C
if (doDFTc.ne.0) then
    allocate(xyz(4,ncenter))
    do i=1,ncenter
        xyz(1,i)=a(i)%index
        xyz(2,i)=a(i)%x
        xyz(3,i)=a(i)%y
        xyz(4,i)=a(i)%z
    end do
    call DFTC(ncenter,xyz,edftc,1)
end if

if ((doJ.ne.0).and.(doDFTc.ne.0)) then
	write(*,"('---=< summing up some final total energies >=---')")
	write(*,"(' Electronic+nuclear energy (from IoG):',f20.10)") intgv
	write(*,"(' -C corrections                      :',f20.10)") edftc
	write(*,"('                                      --------------------')")
	write(*,"('   +--> Total DFT-C ENERGY           :',f20.10' <--+',/)") intgv+edftc
end if

CALL CPU_TIME(time_end)
call walltime(iwalltime2)
write(*,"('---=< ALL GOOD! Total job time: CPU time',f12.2,'s, wall clock time',i10,'s >=--- ',/)") time_end-time_begin,iwalltime2-iwalltime1

end program








