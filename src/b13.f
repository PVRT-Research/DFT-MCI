      SUBROUTINE b13func(F,RA,RB,D1RA,D1RB,RLapA,RLapB,
     $           TA,TB,EX_HF_DENSA,EX_HF_DENSB,FUNCPAR)
c    ******************************************************************
c    *                                                                *
c    *  evaluates the non-dynamic correlation energy by the           *
c    *  real-space post-HF model B13 of Becke with RI and fully SCF.  *
c    *                                                                *
c    *  reference: A. D. Becke, J. Chem. Phys. 2013, 138, 074109.     *
c    *                                                                *
c    *  OUTPUT:                                                       *
c    *     F    - Functional values                                   *
c    *                                                                *
c    *  INPUT:                                                        *
c    *     RA,B   - Spin densities                                    *
c    *     D1RA,B   - Spin densitiy gradients                         *
c    *     TA,B   - Spin kinetic energy densities                     *
c    *     DLapA,B - Laplacian of the spin densities                  *
c    *     EX_HF_DENSA,B - exact HF-like exchange energy density      *
c    *     FUNCPAR: 1 for B13 parameters                              *
c    *              2 for B13SC (strong correlation) parameters       *
c    *              3 for B13SC for Transition Metals (see ref. below)*
c    *                                                                *
c    *                        ---WARNING---                           *
c    * For calculations of energy differences using B13SC we need a   *
c    * correction for the BR exchange-hole normalization, as          *
c    * explained in Appendix C of B13 paper. However, such correction *
c    * is very unsatisfying (and the biggest disappointment of the    *
c    * B13 functional) and it requires Hirshfeld's atomic partition   *
c    * scheme, and IS NOT IMPLEMENTED HERE!!!                         *
c    *                                                                *
c    ******************************************************************
c
        IMPLICIT REAL*8(A-H,O-Z)
        INTEGER FUNCPAR
        parameter(third=1.0d0/3.0d0,zero=0.0D0)
        parameter(dsmall=1.0D-15)
        parameter(dsmallb=1.0D-12)
        parameter(dsmall2=1.0D-14)
        parameter(detol=1.0D-08,dtol=1.0D-08)
        parameter(dtol2=1.0D-08)
        parameter(Tolbig=1.0d+15)
        parameter(delt=0.050d0)  ! smoothening of fa, fb
        parameter(smoth=115.0D0) ! unifies the definition of fcor
        parameter(smoth2=120.0D0)
        
        Tol = 1.0D-08
        smth1 = 1.0d0 - delt
        smth2 = 1.0d0 + delt
C PARAMETERS
        IF (FUNCPAR.EQ.1) THEN
C       B13 Parameters (JCP 2013, 138, 074109)
        aSop  = 0.59D0
        aSpar = 0.61D0
        aDop  = 0.66D0
        aDpar = 0.61D0

		aSC2  = 0.0D0
		aSC3  = 0.0D0
c
c-RPV- This is a flag to eventually turn on BR corrections (not yet implemented)
        DOBRCORR = 0
c
        ELSEIF (FUNCPAR.EQ.2) THEN
C       B13 strongC Parameters (JCP 2013, 138, 074109)
        aSop  = 0.534D0
        aSpar = 0.746D0
        aDop  = 0.635D0
        aDpar = 0.536D0
c
		aSC2  = 0.526D0
		aSC3  = 0.0D0
c
        DOBRCORR = 1
c
        ELSEIF (FUNCPAR.EQ.3) THEN
C       B13 strongC Parameters tuned for Transition Metals (JCP 2013, 138, 161101)
        aSop  = 0.552D0
        aSpar = 0.844D0
        aDop  = 0.640D0
        aDpar = 0.559D0
c
        aSC2  = 0.825D0
        aSC3  = -0.380D0
c
        DOBRCORR = 1
c
        ENDIF

        fc1 = zero
        fc2 = zero
        fcor = zero
        UCA = zero
        UCB = zero
        DSA =dsmall
        DSB =dsmall
        DDA = zero
        DDB = zero
        ASA = zero
        ASSA = zero
        AASA = zero
        ASB = zero
        ASSB = zero
        BBSB = zero


        EFNA = 1.0d0
        EFN2A = 1.0d0
        EFNB = 1.0d0
        EFN2B = 1.0d0
        UX1= -zero
        UX2= -zero
        gra=zero
        grb=zero
        tauA =zero
        tauB =zero
        DLapA =zero
        DLapB =zero

        VX1= -dsmall
        VX2= -dsmall
        SM1A = zero
        SM2A = zero
        SM1B = zero
        SM2B = zero
        ttgg = zero
        ttg1 = zero
        ttg2 = zero
        ttdB = zero
        ttdA = zero

        gc1 = zero
        gcc1 = zero
        ggg = zero
        ggp = zero
        gc2 = zero
        gcc2 = zero

        qqx = 0.0d0
        eey1=0.0d0
        eey2=0.0d0

        qqx = 0.0d0
        smtx = zero
        tx = 0.0d0
        qqfa = 0.0d0
        qqfb = 0.0d0
        xxx = zero
        yy1 = zero
        yy2 = zero
C-RPV-
        uSpar = zero
        uDpar = zero
        uSop = zero
        uDop = zero
        uCtot = zero
        Fdop = zero
        Fdpar = zero
        Fsop = zero
        Fspar = zero
        DESC = zero

        D1 = 0.0d0
        D2 = 0.0d0
        DF1 = 0.0d0
        DF2 = 0.0d0
        OD1 = 0.0d0
        OD2 = 0.0d0

        uSpar = zero
        uDpar = zero
        uSop = zero
        uDop = zero
        uCtot = zero
        Fdop = zero
        Fdpar = zero
        Fsop = zero
        Fspar = zero
        DESC = zero

        ZAA=zero
        ZBB=zero
        ZAB=zero
        ZFA=zero
        ZFB=zero
        uDA=zero
		uDB=zero
		fssa=zero
		fssb=zero

        fc1=zero
        fc2=zero
        fcor = zero
        UCA = zero
        UCB = zero
        DSA =dsmall
        DSB =dsmall
        DDA = zero
        DDB = zero
        ASA = zero
        ASSA = zero
        ASB = zero
        ASSB = zero

        EFNA = 1.0d0
        EFN2A = 1.0d0
        EFNB = 1.0d0
        EFN2B = 1.0d0
        UX1= -zero
        UX2= -zero
        gra=zero
        grb=zero
        tauA =zero
        tauB =zero

        SM1A = zero
        SM2A = zero
        SM1B = zero
        SM2B = zero
        ttgg = zero
        ttg1 = zero
        ttg2 = zero
        ttdB = zero
        ttdA = zero

        gc1 = zero
        gcc1 = zero
        ggg = zero
        ggp = zero
        gc2 = zero
        gcc2 = zero

        qqx = 0.0d0
        eey1=0.0d0
        eey2=0.0d0

        qqx = 0.0d0
        smtx = zero
        tx = 0.0d0
        qqfa = 0.0d0
        qqfb = 0.0d0
        xxx = zero
        yy1 = zero
        yy2 = zero
        fp1=zero
        fp2=zero
c
        D1 = RA
        D2 = RB

        UX1 = -dsmall2
        UX2 = -dsmall2
        VX1=EX_HF_DENSA
        VX2=EX_HF_DENSB

        IF((VX1.ge.dsmall).or.(VX2.ge.dsmall)) GO TO 666
        
        IF((D1.gt.detol).OR.(D2.gt.detol)) then
            IF(D1.gt.detol) then
                UX1 = VX1/D1
                gra = D1RA
                tauA = TA
                DLapA = RLapA
c brsc is the routine doing the main RSC numerics and it's in beck_rsc.F
                   CALL brsc(D1,DF1,gra,DLapA,tauA,VX1,UX1,EFNA,EFN2A
     $      ,SM1A,SM2A,DSA,ndrvA)
            ELSE  ! D1
                SM1A = zero
                SM2A = zero
                DSA = dsmall
                gra = 0.0d0
                tauA = zero
                DLapA = zero
                UX1 = -dsmall
                EFNA = 1.0d0
                EFN2A = 1.0d0
            ENDIF ! D1, VX1
c
            IF(D2.gt.detol) then
                UX2 = VX2/D2
                grb = D1RB
                tauB = TB
                DLapB = RLapB
                CALL brsc(D2,DF2,grb,DLapB,tauB,VX2,UX2,EFNB,EFN2B
     $      ,SM1B,SM2B,DSB,ndrvB)
c
            ELSE  !  D2
                SM1B= zero
                SM2B= zero
                DSB= dsmall
                grb = zero
                tauB = zero
                DLapB = zero
                UX2 = -dsmall
                EFNB = 1.0d0
                EFN2B = 1.0d0
            ENDIF  ! VX2
c
            fcor = 0.0d0
            ffcor =0.0d0
            FSop = 0.0d0
            FDop = 0.0d0
            IF((ndrvA.eq.0).or.(ndrvB.eq.0)) GO TO 777
c      If((NA.gt.0).and.(NB.gt.0)) then
            if((D1.gt.detol.and.ndrvA.eq.2).and.(D2.gt.detol
     $      .and.ndrvB.eq.2)) then
c
                IF(EFN2B.gt.dtol2) gcc1 = (1.0d0-EFN2A)/EFN2B
                IF(EFN2A.gt.dtol2) gcc2 = (1.0d0-EFN2B)/EFN2A
c           
                IF(EFNB.gt.dtol2) then
                    gc1 = (1.0d0-EFNA)/EFNB
                    ggg= gc1-1.0d0-delt
                    IF((gc1.gt.smth1).and.(gc1.lt.smth2)) then
                        fc1 = 1.0d0-ggg*ggg/(4.0D0*delt)
                        dfc1g1 = -ggg/(2.0D0*delt)
                    elseif (gc1.ge.smth2) then
                       fc1 = 1.0d0
                       dfc1g1 = 0.0d0
                    else
                        fc1 = gc1
                        dfc1g1 = 1.0D0
                    endif   ! gc1
                endif
c
                IF(EFNA.gt.dtol2) then
                    gc2 = (1.0d0-EFNB)/EFNA
c          gc2 = 1.0d0/EFNA-EFNB/EFNA
                    ggp= gc2-1.0d0-delt
                    IF((gc2.gt.smth1).and.(gc2.lt.smth2)) then
                       fc2 = 1.0d0-ggp*ggp/(4.0D0*delt)
                       dfc2g2 = -ggp/(2.0D0*delt)
                    elseif (gc2.ge.smth2) then
                        fc2 = 1.0d0
                        dfc2g2 = 0.0d0
                    else
                        fc2 = gc2
                        dfc2g2 = 1.0D0
                    endif   ! gc2
                endif
c
                tx = 1.0d0
                qqx = 0.50d0
                smtx = 0.0d0
                xxx=zero
                if((abs(fc1).gt.detol).or.(abs(fc2).gt.detol)) then
                    fp1 = fc1*fc1
                    fp2 = fc2*fc2
                    xxx = (fc1-fc2)/(fp1+fp2)
                    smtx = smoth*xxx
                    if(abs(smtx).lt.99.0d0) then
                        tx = exp(smtx)
                        qqx = 1.0d0/(tx+1.0d0)
                    elseif(smtx.lt.-99.0d0) then
                        tx= zero
                        qqx= 1.0d0
                    else
                        tx= zero
                        qqx= zero
                    endif  !smtx
                    fcor = fc1*qqx + fc2*(1.0d0-qqx)
                    if(gc1.le.-dsmall) gc1=0.0d0
                    if(gcc1.le.-dsmall) gcc1=0.0d0
                    if(gc2.le.-dsmall) gc2=0.0d0
                    if(gcc2.le.-dsmall) gcc2=0.0d0
                    ffcor = min(gcc1,gcc2,1.0d0)
                else  !  fc1,fc2
                   xxx=zero
                   qqx = 0.0d0
                   tx = 0.0d0
                   smtx=zero
                   fcor = zero
                endif  !  fc1,fc2
c
                ttdA = zero
                ttdB = zero
                IF(D2.gt.dtol) then
                    ttdB = UX2*D1/D2
                endif
                IF(D1.gt.dtol) then
                    ttdA = UX1*D2/D1
                endif
                ttg2=D1*UX2
                ttg1=D2*UX1
                ttgg = ttg1+ttg2
                uSop = 0.50d0*fcor*ttgg
                FSop=aSop*uSop
C-RPV start FDOP
                IF((UX1.lt.(-dtol)).and.(UX2.lt.(-dtol)))
     $        ZAB = -0.630D0*(EFNA/UX1+EFNB/UX2)
                ZSF = (ZAB**3)/(1.0D0+ZAB)
                uDop =-0.80d0*(1.0d0-fcor)*D1*D2*ZSF
                FDop= aDop*uDop
C-RPV- end FDOP
            else  ! NA
                gc1=zero
                gcc1 = zero
                fc1 = zero
                dfc1g1 = zero
                gc2=zero
                gcc2 = zero
                fc2 = zero
                dfc2g2 = zero
                xxx = zero
                tx = 0.0d0
                qqx = 0.0d0
                smtx = zero
                fcor = zero
                tx = 0.0d0
                FSop = 0.0d0
                FDop = 0.0d0
                ttdA=0.0d0
                ttdB=0.0d0
            endif   ! NA
cccccccccccc
 777  CONTINUE
            IF((ndrvA.eq.0).and.(ndrvB.eq.0)) GO TO 666
            yy1 = zero
            smtx2 = zero
            eey1 = 0.0d0
            ty1 = 0.0d0
            yy2 = zero
            smtx3 = zero
            eey2=0.0d0
            ty2 = 0.0d0
c
            ASSA =0.0d0
            A2SA = 0.0d0
            AASA = 0.0d0
            IF((D1.gt.detol).and.(ndrvA.eq.2)) then
                IF(abs(SM2A).gt.dtol) then
                    ASA=(1.0d0-EFNA-fcor*EFNB)/SM2A
                    A2SA=(1.0d0-EFN2A-ffcor*EFN2B)/SM2A
                ENDIF
                DDA = DSA/(3.0d0*D1)
                yy1 = ASA-DDA
                smtx2 = smoth2*yy1
c           ty1 = exp(smtx2)
                if(abs(smtx2).lt.99.0d0) then
                    ty1 = exp(smtx2)
                    eey1 = 1.0d0/(ty1+1.0d0)
                elseif(smtx2.lt.-99.0d0) then
                    ty1= zero
                    eey1 = 1.0d0
                else
                    ty1= zero
                    eey1 = zero
                endif
                ASSA = ASA*eey1+DDA*(1.0d0-eey1)
                AASA = min(A2SA,DDA)
            endif  ! D1
ccccccccccccccccccccc
            ASSB=0.0d0
            A2SB = 0.0d0
            BBSB = 0.0d0
            IF((D2.gt.detol).and.(ndrvB.eq.2)) then
                IF(abs(SM2B).gt.dtol) then
                    ASB=(1.0d0-EFNB-fcor*EFNA)/SM2B
                    A2SB=(1.0d0-EFN2B-ffcor*EFN2A)/SM2B
                endif
                DDB = DSB/(3.0d0*D2)
                yy2 = ASB-DDB
                smtx3 = smoth2*yy2
c           ty2 = exp(smtx3)
                if(abs(smtx3).lt.99.0d0) then
                    ty2 = exp(smtx3)
                    eey2 = 1.0d0/(ty2+1.0d0)
                elseif(smtx3.lt.-99.0d0) then
                    ty2= zero
                    eey2 = 1.0d0
                else
                    ty2= zero
                    eey2 = zero
                endif
                ASSB = ASB*eey2+DDB*(1.0d0-eey2)
                BBSB = min(A2SB,DDB)
            endif   ! D2
ccccccccc
            UCA=0.0d0
            UCB=0.0d0
            if(ndrvA.eq.2)  UCA = - ASSA*SM1A
            if(ndrvB.eq.2)  UCB = - ASSB*SM1B
            uSpar = (D1*UCA+D2*UCB)*0.50d0
c        else  !  NMO
c          uSpar= zero
c        endif
ccccccccccccccc
            FSpar = aSpar*uSpar
C-RPV- start FDPAR
            IF(DDA.gt.dtol) fssa=ASSA/DDA
            IF(DDB.gt.dtol) fssb=ASSB/DDB
            IF(UX1.lt.(-dtol)) ZAA = -0.880d0*2.0d0*EFNA/UX1
            IF(UX2.lt.(-dtol)) ZBB = -0.880d0*2.0d0*EFNB/UX2
            ZFA=(ZAA**5)/(1.0D0+ZAA*0.5D0)
            ZFB=(ZBB**5)/(1.0D0+ZBB*0.5D0)
            uDA=-0.005D0*(1.0D0-fssa)*D1*DSA*ZFA
            uDB=-0.005D0*(1.0D0-fssb)*D2*DSB*ZFB
            uDpar = uDA + uDB
            FDpar = aDpar*uDpar
C-RPV- end FDPAR  
c
C-RPV- start B13 strongC
		    uCtot = uDop + uDpar + uSop + uSpar
            xfac=zero
		    IF(uCtot.gt.dtol) xfac=(uSop+uSpar)/uCtot
		    DESC= aSC2*xfac**2*uCtot + aSC3*xfac**3*uCtot
C-RPV end
            F = FSop + FSpar + FDop + FDpar + DESC
c
!       ENDIF !rho.gt.tol
c Begin building the SCF potential:

        ENDIF  ! RA,B

 666    CONTINUE
      RETURN
      END
      
      
      
      
      SUBROUTINE brsc(DN,DFN,GR,DLAP,tau,VX,UX,EFN,EFN2,SM1,SM2,DS,ndrv)
c
      IMPLICIT REAL*8(A-H,O-Z)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                              c
c      DN is the electron density for the input spin direction c
c      (i.e. either D1, or D2, but not the sum of the two);    c
c                                                              c
c      GR is the dot product of the  density gradient with     c
c      iself (per sin direction), NOT the grad.modulus;        c
c                                                              c
c      tau is the electron kinetic energy density for the      c
c      given spin direction                                    c
c                                                              c
c      DLAP is the density Laplacian                           c
c                                                              c
c      VX is the exact Slater exchange potential multiplied    c
c      by DN (for the given input spin direction)              c
c                                                              c
c     OUTPUT                                                   c
c     UX is the effective X Slater potential                   c
c                                                              c
c     EFN is the effective X hole normalizations               c
c                                                              c
c     EFN required to evaluate the functional derivativres     c
c     later on                                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       parameter(zero=0.0d0,one=1.0d0,two=2.0d0,four=4.0d0)
       parameter(five=5.0d0,six=6.0d0,three=3.0d0)
       parameter(third=1.0d0/3.0d0)
       parameter(twthr=2.d0/3.d0, trihf=3.0d0/2.0d0)
       parameter(pi=3.141592653589793d0)
c
       parameter(aa1=0.9301841768238374D0,aa2=0.5485569431916153D0)
       parameter(ss1=0.9650036832623640d0)
       parameter(ss2=0.93018417682383740d0)
       parameter(ss3=0.54855694319161530d0)
       parameter(ss4=1.5158242410120670d0)
c
       parameter(c0 = -5.9685280129072020D0 )
       parameter(c1 = -2.1837477426038480D0 )
       parameter(c2 = -4.9858864412437560D0 )
       parameter(c3 = -1.1341612120636830D0 )
       parameter(c4 = -1.6921426426199750D0 )
       parameter(c5 =  0.57089595383468940D0)
c
       parameter(b0 = -5.9685280130660880D0 )
       parameter(b1 = -2.030780232084790D0  )
       parameter(b2 = -4.6796750480012860D0 )
       parameter(b3 = -1.1188490577541170D0 )
       parameter(b4 = -1.8087055034029230D0 )
       parameter(b5 =  0.59226482161752910D0)
c
       parameter(d1 = 10.850190580831960d0  )
       parameter(d2 = 423.51564269756920d0  )
       parameter(d3 = 45.317325658682150d0  )
       parameter(d4 = 9.6414424525059270d0  )
       parameter(d5 = 1.9494258703373230d0  )
       parameter(d6 = 2.4210661522846020d0  )
       parameter(d7 = 2.1892374225965940d0  )
c
       parameter(e1 = 613.92073109477340d0     )
       parameter(e2 = 60.413933318777930d0     )
       parameter(e3 = 13.001906863306460d0     )
       parameter(e4 = 3.4757233899623390d0     )
       parameter(e5 = 2.5384493464526640d0     )
       parameter(e6 = 2.5561187910239140d0     )
       parameter(ee0 = 3.360517076724290D0     )
       parameter(ee1 = 4.623703278485152D0     )
       parameter(ee2 = 2.688949840405010D0     )
       parameter(ee3 = 0.6007166968496472D0    )
       parameter(ee4 = 0.03922040006408070D0   )
       parameter(ee5 = 0.0005438465669613952D0 )
       parameter(ee6 = 0.00000078437439010087D0)
c
       parameter(q0 =   0.9129908719446399D0  )
       parameter(q1 =   3.655262558462426D0   )
       parameter(q2 =   0.1801828684494572D0  )
       parameter(q3 =  -3.062938667772561D0   )
       parameter(q4 =  -1.173405187745653D0   )
       parameter(q5 =  -1.662674088158794D0   )
       parameter(q6 =   0.6859613559654089D0  )
       parameter(q7 =   0.06595477584967326D0 )
       parameter(q8 =  -0.03038609318852905D0 )
       parameter(q9 =  -0.00000000000000077D0 )
c
       parameter(s0 =   8.382230306773296D0         )
       parameter(s1 =   19.60257290836815D0         )
       parameter(s2 =   19.71894106502952D0         )
       parameter(s3 =   10.77146542063455D0         )
       parameter(s4 =   3.415698370189622D0         )
       parameter(s5 =   0.5813153924917321D0        )
       parameter(s6 =   0.05426084060061605D0       )
       parameter(s7 =   0.002299629631716270D0      )
       parameter(s8 =   0.00005119354330427682D0    )
       parameter(s9 =   0.000000322977561012273D0   )
       parameter(s10 =  0.000000001405232963383258D0)
c
       parameter(a1 =  1.23767D0)
       parameter(a2 =  9.37587D0)
       parameter(a3 = 19.4777D0)
       parameter(a4 = 13.6816D0)
       parameter(a5 =  0.078655D0)
       parameter(a6 = 54.7264D0)
       parameter(a7 = 58.4331D0)
       parameter(a8 = 18.75174D0)
       parameter(a9 =  1.23767D0)
c
       parameter(dlam=1.0D0)
       parameter(AXC= 1.24070098179880D0)
       parameter(Tolbig=1.0d+10)
       parameter(dsmall=1.0D-13)
       parameter(dtol2=1.0D-08)
       parameter(dtol4=1.0D-06)
       parameter(dtol3=1.0D-08)
       parameter(detol0=1.0D-08)
       parameter(detol=1.0D-08,dtol=1.0D-08)
       parameter(delt=0.070d0)  ! smoothening of EFN

       REAL*8 xx1,xx2,yy,acc,X_y,EFN,DN,GR,DLAP,tau,SM1,SM2,DS
       INTEGER iter
ccccccccccccfor the numerical version onlyccccc
c      external FF
c      external rtsafe
cccccccccccccccccccccccccccccccccccccccccccccc
       Tol = 1.0D-08
c
           smth1 = 2.0d0 - delt
           smth2 = 2.0d0 + delt
c Initializations start
          dfdr=zero
          dfdt=zero
          dfdl=zero
          dfdVx = zero
          yy = dsmall
          X_y= dsmall
          g_y=zero
          P12_y=zero
          Q12_y=zero
          dxdy=zero
          dUdx=zero
          dUdr=zero
          dydn=zero
          dydg=zero
          dydt=zero
          dydl=zero
          dydVx=zero
          dQdn=zero
          dQdt=zero
          dQdl=zero
          dQdg=zero
          dNdr=zero
          dNdt=zero
          dNdl=zero
          dNdg=zero
          dNdVx=zero
          dM1dr = zero
          dM1dt = zero
          dM1dl = zero
          dM1dg = zero
          dM1dVx = zero
          dM2dr = zero
          dM2dt = zero
          dM2dl = zero
          dM2dg = zero
          dM2dVx = zero
          dg3dn = zero
          dg3dt = zero
          dg3dg = zero
          EFN = 1.0d0
          EFN1 = 1.0d0
          EFN2 = 1.0d0
          gc1=zero
          SM1 = 0.0d0
          SM2 = 0.0d0
          t1=zero
          t2=zero
          t3=zero
          t4=zero
          t5=zero
          t6=zero
          t7=zero
          t8=zero
          t9=zero
          t10=zero
          t11=zero
          t12=zero
          t13=zero
          t14=zero
          t15=zero
          t16=zero
          t17=zero
          tt9=zero
          ts14=zero
          ts15=zero
          ts16=zero
          ts23=zero
          tq1=zero
          tp6=zero
          tp7=zero
          t01=zero
          t02=zero
          tn1=zero
          tn2=zero
          tn3=zero
          tn4=zero
          tn5=zero
          tu1=zero
          tq1 = zero
          ndrv=2
c Initializations end
c
c    First calculate X_y and the energy density
c
           DN13 = DN**third
           DN43= DN13*DN
           DN53= DN*(DN13*DN13)
           pi13=pi**third
           pi23 = pi13*pi13
c      IF(tau.gt.dsmall) then
           if(tau.lt.dsmall) tau=dsmall
c      IF(DN.gt.detol) then
           t1 = four*pi*DN
           tq1 = t1*DN
           DS = tau-GR/(DN*four)
           QS = DLAP/six -DS/three
c  To keep in mind that VX = Ux*DN:
          IF(tq1.gt.dtol2) then
             yy = -three*QS*UX/tq1
          endif
          if(abs(yy).lt.dsmall) ndrv =0
c            yy = -three*QS*UX/tq1
c       if(yy.gt.-Tolbig.and.yy.lt.Tolbig) then
cccccccccccfor the numerical version onlyccccccccccccccc
c           xx1=dtol
c           xx2 = 100.0d0
c        IF(yy.le.(-dsmall)) then
c           xx1 = detol
c           xx2 = 2.00000d0
c        elseif ((yy.ge.detol).and.(yy.lt.1000.0d0)) then
c           xx1 = 2.0000000d0
c           xx2 = 10.50d0
c        elseif ((yy.ge.1000.0d0).and.(yy.lt.Tolbig)) then
c           xx1 = 9.0d0
c           xx2 = 100.0d0
c        endif
c           acc = 1.0D-8
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           cc = two*pi13*DN43
           cc1=twthr*pi13*DN13
           cc2=two*pi13*DN13
cccccc set up for numerical evaluation of X_y cccccccccc
c          X_y = rtsafe(FF, yy, xx1, xx2, acc)
c            GO TO 115
cccccc end of set up for numerical eval. of X_y ccccccc
c        IF(yy.le.(-detol)) then
         IF(yy.le.(-dsmall)) then
c        Region I : -infinity <= y <= 0
             t3 = atan(aa1*yy + aa2)
             t5 = atan(aa2)
             g_y = two*(two*t3 + pi)/(two*t5 + pi)
             t2 = yy*yy
             t4 = t2*yy
             t6 = t2*t2
             t8 = t6*yy
             tt9 = b0 +b1*yy +b2*t2 +b3*t4 +b4*t6 +b5*t8
             P12_y = dsmall
           IF(abs(tt9).gt.detol)
     $       P12_y = (c0+c1*yy+c2*t2+c3*t4+c4*t6+c5*t8)/tt9
             X_y=g_y*P12_y
c
             t5 = atan(aa2)
             tp1 = aa1*aa1
             tp7 = aa2*aa2
             dxx = 4.0d0*aa1/(one+tp1*t2+two*aa1*yy*aa2 + tp7)/
     $             (two*t5 + pi)
             ts14 = t6*t6
             ts15 = ts14*yy
             ts23 = tt9*tt9
             Q12_y = zero
           IF(ts23.gt.detol)
     $  Q12_y = (q0 +q1*yy +q2*t2 +q3*t4 +q4*t6 +q5*t8 +q6*t6*t2
     $        + q7*t6*t4 +q8*ts14 +q9*ts15)/ts23
             dxdy= dxx*P12_y + g_y*Q12_y
c
         elseif((yy.ge.dsmall).AND.(yy.le.1000.0D0)) then
c          Region II : 0 < y < 1001
             t2 = yy*yy
             t4 = t2*yy
             t6 = t2*t2
             t8 = t6*yy
             t5 = t4*t4
             t9 = ee0 +ee1*yy +ee2*t2 +ee3*t4 +ee4*t6 +ee5*t8 +ee6*t5
             t10 = t2+e5*yy+e6
             X_y=dsmall
         IF(abs(t10).gt.detol)
     $      X_y = d1*(yy+d2)*(yy+d3)*(yy+d4)*(yy+d5)*(t2+d6*yy+d7)
     $         / (yy+e1)/(yy+e2)/(yy+e3)/(yy+e4)/t10
c
             ts14 = t6*t6
             ts15 = ts14*yy
             ts16 = t8*t8
             ts23 = t9*t9
              dxdy=zero
           IF(ts23.gt.detol)
     $       dxdy = (s0+s1*yy+s2*t2+s3*t4+s4*t6+s5*t8+s6*t6*t2
     $            + s7*t6*t4+s8*ts14+s9*ts15+s10*ts16)/ts23
c
         elseif(yy.gt.1000.0D0) then
c            Region III : 1000 <= y
             t1 = log(yy)
             t2 = one/t1
c            t4 = t1**(one+t2)
             t22 = (t1+one)/t1
             t4 = t1**t22
             t6 = t1*t1
             t12 = t6*t6
        X_y = t4+a1/t1+a2/t6-a3/t6/t1+a4/t12-a5
             tp6 = log(t1)
             t7 = t1**t2
             dxdy1 = t7*(one - tp6/t1 + t2)/yy
             t8 = t2*t2/yy
             dxdy2 = t8*(-a9 - a8/t1 + a7/t6 - a6/t6/t1)
             dxdy = dxdy1+dxdy2
        endif
            IF(abs(yy).lt.dsmall) X_y=2.000010d0
c   end of calculating X(y)
c 115     CONTINUE   ! uncomment only for the numerical
c            twthr=2/3, trihf= 3/2
             t01 = sqrt(two)
             t02 = sqrt(three)
             tn0 = (twthr**trihf)*pi
             tn1 = DN*DN
             tn2 = sqrt(DN)
             tn4 = exp(X_y)
             tn5 = exp(-X_y)
             tu1 = X_y*X_y
             ts12 = 0.0d0
             tn10 = dsmall
             ts10 = dsmall
             xxxq = X_y*QS
c       IF((abs(X_y-2.0d0).gt.detol).and.X_y.gt.detol) then
        IF(X_y.gt.detol) then
             ts12 = X_y+4.0d0/X_y-(1.0d0+4.0d0/X_y)*tn5
          IF(abs(xxxq).gt.detol)
     $       tn10 = (X_y - two)*two/(xxxq*three)
           if(tn10.gt.dsmall) then
             ts10 = tn10*DN/4.0d0
c            tn11 = sqrt(abs(tn10))
             tn11 = sqrt(tn10)
c            ts11 = sqrt(abs(ts10))
             ts11 = sqrt(ts10)
             ALF3 = tn10*tn11
             EFN1 = ALF3*pi*tn4*DN*DN*sqrt(DN)
             EFN2 = EFN1
             IF(EFN2.gt.1.0d0) EFN2=1.0d0
           ggg= EFN1-2.0d0-delt
         IF((EFN1.gt.smth1).and.(EFN1.lt.smth2)) then
           EFN = 2.0d0-ggg*ggg/(4.0D0*delt)
           dNg1 = -ggg/(2.0D0*delt)
        elseif (EFN1.ge.smth2) then
           EFN = smth2
           dNg1 = 0.0d0
        else
            EFN = EFN1
           dNg1 = 1.0d0
         endif   ! EFN1
             SM1 = EFN*ts11*ts12
             SM2 = EFN*ts10*(tu1+12.0d0)
         else  !tn10
           EFN = 1.00d0
           EFN2=1.0d0
           SM1 = 0.0d0
           SM2 = 0.0d0
        endif   ! tn10
c
c   Now prepare for the SCF potential:
c       IF(DN.gt.dtol2)  then
            t3 = DN*DN
            t4 = 9.0d0*QS+GR/(four*DN)
            dydn= t4*UX/(four*DN)/(pi*DN)/DN
            dQdn = -GR/(four*DN)/(three*DN)
            dydg = -UX/(DN*pi)/(four*DN)/(four*DN)
            dQdg = one/(DN*12.0D0)
            dydl = -UX/(pi*DN*two)/(four*DN)
            dydt = UX/(pi*DN)/(4.0D0*DN)
            dQdl = one/six
            dydVx = -three*QS/(four*DN)/(pi*DN)/DN
            dQdt = - third
            IF(DN.gt.dtol4) then
            dg3dn = -tau/DN/(three*DN)+GR/DN/(DN*six)/DN
            dg3dt = one/(three*DN)
            dg3dg = -one/(four*DN)/(three*DN)
            else
              dg3dn = 0.0d0
              dg3dt = 0.0d0
              dg3dg = 0.0d0
            endif
            dNdn = dNg1*5.0d0*EFN1/(two*DN)
            tt12 = (X_y-two)/(X_y*QS)
            dNdx = zero
         if(tt12.gt.dsmall) then
            ts15 = sqrt(tt12)
c           ts15 = sqrt(abs(tt12))
           ts1 = sqrt(two)
            ts2 = sqrt(three)
            ts4 = DN*DN
            ts5 = sqrt(DN)
            ts9 = exp(X_y)
            ts17 = X_y*X_y
            dNdx = dNg1*(2.0d0/9.0d0)*ts1*ts2*pi*ts5*ts4*ts9*
     $      ts15*(ts17-2.0D0*X_y+3.0d0)/X_y/(X_y*QS)
          endif
            dNdQ = -dNg1*EFN1*three/(two*QS)
c
           t1 = DN*DN
           t6 = exp(X_y)
           t7 = X_y*X_y
           t12 = t7*t7
           t14 = QS*QS
           dM1dn = SM1*three/DN
           dM1dQ = -SM1*two/QS
           dM1dx = (two/9.0d0)*pi*t1*DN*(X_y-2.0d0)
     $     *(-t6*t7*X_y-12.0d0*t6*X_y+t6*t12+6.0d0*t6*t7+24.0d0*t6
     $     - 24.0d0)/t7/X_y/(QS*X_y)/QS
c
            t3 = sqrt(6.0d0)
            t4 = X_y - 2.0d0
            t5 = t4*t4
            t8 = X_y*X_y
            t10 = QS*QS
            tt13 = DN*t4/(X_y*QS)
            dM2dx = zero
           if(tt13.gt.(dtol)) then
            t18 = sqrt(tt13)
            t19 = exp(X_y)
            t23 = t8*t8
          dM2dx = pi*t1*DN*t3*t4/t8/X_y/QS
     $    *(t18*t19)/QS*(13.0d0*t8+60.0d0+t23-24.0d0*X_y)/27.0d0
        endif
           dM2dQ = -SM2*5.0d0/(two*QS)
           dM2dn=SM2*7.0d0/(two*DN)
c       else  !  dtol2
c         dydn = zero
c          dQdn = zero
c          dydg = zero
c          dQdg = zero
c          dydl = zero
c          dydt = zero
c          dQdl = zero
c          dydVx = zero
c          dQdt = zero
c          dg3dn = zero
c          dg3dt = zero
c          dg3dg = zero
c          dNdn = zero
c          dNdx = zero
c          dNdQ = zero
c          dM1dn=zero
c          dM1dx=zero
c          dM1dQ=zero
c          dM2dn=zero
c          dM2dx=zero
c          dM2dQ=zero
c      endif   ! dtol2
c    here dfdt == dNdt, etc.
         dNdr = dNdn + dNdx*dxdy*dydn +dNdQ*dQdn
         dNdt = dNdx*dxdy*dydt + dNdQ*dQdt
         dNdl = dNdx*dxdy*dydl+ dNdQ*dQdl
         dNdg = dNdx*dxdy*dydg + dNdQ*dQdg
         dNdVx = dNdx*dxdy*dydVx
         dM1dr = dM1dn + dM1dx*dxdy*dydn + dM1dQ*dQdn
         dM2dr = dM2dn + dM2dx*dxdy*dydn +dM2dQ*dQdn
         dM1dt = dM1dx*dxdy*dydt + dM1dQ*dQdt
         dM2dt = dM2dx*dxdy*dydt + dM2dQ*dQdt
         dM1dg = dM1dx*dxdy*dydg + dM1dQ*dQdg
         dM2dg = dM2dx*dxdy*dydg + dM2dQ*dQdg
         dM1dl = dM1dx*dxdy*dydl + dM1dQ*dQdl
         dM2dl = dM2dx*dxdy*dydl + dM2dQ*dQdl
         dM1dVx = dM1dx*dxdy*dydVx
         dM2dVx = dM2dx*dxdy*dydVx
       else   ! abs(X_y-2.0d0)
          dNdr = zero
          dNdt = zero
          dNdl = zero
          dNdg = zero
         dNdVx = zero
         dM1dr = zero
         dM2dr = zero
         dM1dt = zero
         dM2dt = zero
         dM1dg = zero
         dM2dg = zero
         dM1dl = zero
         dM2dl = zero
        dM1dVx = zero
        dM2dVx = zero
          EFN = 1.0d0
       endif  ! X_y
c      endif  ! Tolbig yy
c     ENDIF  ! DN and  tau
      RETURN
      END
      