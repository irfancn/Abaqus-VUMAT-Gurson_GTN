c ======================================================================
c subroutine for GTN
c All rights of reproduction or distribution in any form are reserved.
c By Irfan Habeeb CN (PhD, Technion - IIT)
c ======================================================================
      subroutine vumat(
C Read only -
     1     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     7     stressNew, stateNew, enerInternNew, enerInelasNew)
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1     charLength(nblock), strainInc(nblock,ndir+nshr),
     2     relSpinInc(nblock,nshr), tempOld(nblock),
     3     stretchOld(nblock,ndir+nshr),
     4     defgradOld(nblock,ndir+nshr+nshr),
     5     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6     stateOld(nblock,nstatev), enerInternOld(nblock),
     7     enerInelasOld(nblock), tempNew(nblock),
     8     stretchNew(nblock,ndir+nshr),
     9     defgradNew(nblock,ndir+nshr+nshr),
     1     fieldNew(nblock,nfieldv),
     2     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3     enerInternNew(nblock), enerInelasNew(nblock)

      character*80 cmname
      integer k, k1, k2, iter, Niter, iterC
C inputs
      real*8 e, xnu, dens, sigy0, epbar0, eprate0, mpw, npw, bg, 
     1 cc, alpha, Cp, chi, temp0, tempIn, f0, q1, q2, fc, ff, 
     2 fndev, fnamp, fnmean, fsdev, fsamp, fsmean, kw, tol
C scalar parameters
      real*8 mu, alamda, bulk, Gur, GurN, epbar, ep, dep,
     1 eprate, eprateN, sigy, sighyd, sigeqv, sigyT, sigyN, dGdep, 
     2 trInc, fT, fs, fN, fNs, sigdotp, epbarN,
     3 tempT, tempN, pi, sq2pi, arg, Astrn, Astrs, Bstrs, sigHrate,
     4 isnuc, sigeff, fdot, J3, Wnh, devdotp, stnx, stnz, stnxz, stn1, 
     5 stn2, stnrat, dWork, dPwork
C Tensor/vectors
      real*8 L(6,6), sOld(6), sNew(6), sigdev(6), sT(6), 
     1 np(6), dsig(6), Pep(6)

      e = props(1)                  ! Young's modulus
      xnu = props(2)                ! poisson's ratio
      dens = props(3)               ! density
      sigy0 = props(4)              ! yield stress
      eprate0 = props(5)            ! ref. plast. strain rate
      npw = props(6)                ! n, plst. strain expon.
      mpw = props(7)                ! m, plst. strain rate expon.
      bg = props(8)                 ! temp. soft. coef.
      cc = props(9)                 ! temp. soft. exponent. 
      alpha = props(10)             ! temp. increment plast. strain
      Cp = props(11)                ! spec. heat
      chi = props(12)               ! work - heat transfer
      temp0 = props(13)             ! ref. temp
      tempIn = props(14)            ! init. temp
      f0 = props(15)                ! initial porosity
      q1 = props(16)                ! q1, GTN
      q2 = props(17)                ! q2, GTN
      fc = props(18)                ! critical f
      ff = props(19)                ! final f
      fndev = props(20)             ! porosity distri. strain
      fnamp = props(21)
      fnmean = props(22)
      fsdev = props(23)             ! porosity distri. stress
      fsamp = props(24)
      fsmean = props(25)
      kw = props(26)                ! NH - shear strain factor
      tol = props(27)               ! tolerance for the N-R iteration

      epbar0 = sigy0/e              ! ref. plast. strain

      mu = e/(2.d0*(1.d0+xnu))
      alamda = e*xnu/((1.d0 + xnu) * (1.d0 - 2.d0*xnu))
      bulk = e/(3.d0*(1.d0 - 2.d0*xnu))

C stiffness matrix
      L = 0.d0
      do k1 = 1, 3
        do k2 = 1, 3
          L(k1, k2) = alamda
        end do 
        L(k1, k1) = alamda + 2.d0*mu
        L(k1+3, k1+3) = 2.d0*mu
      end do 
C -------------------- simulation first step & later -------------------
      do 30 k = 1, nblock
      if (stateOld(k, 1) .eq. 0.d0) then
        go to 10
      else 
        go to 20
      end if 
C----------------------------- initial state ---------------------------
 10   trInc = sum(strainInc(k, 1:3))
      do k1 = 1, 3
        stressNew(k, k1) = alamda*trInc + 2.d0*mu*strainInc(k, k1)
        stressNew(k, k1+3) = 2.d0*mu*strainInc(k, k1+3)
      end do 
      sigeff = sigy0 + sum( stressOld(k, 1:3) )/3.d0
      stnx = stateOld(k, 11) + strainInc(k, 1)
      stnz = stateOld(k, 12) + strainInc(k, 3)
      stnxz = stateOld(k, 13) + strainInc(k, 5)
      stn1 = 0.50*(stnx + stnz) + 0.5d0*sqrt((stnx - stnz)**2.d0 
     1  + 4.d0 * stnxz**2.d0)
      stn2 = 0.50*(stnx + stnz) - 0.5d0*sqrt((stnx - stnz)**2.d0 
     1  + 4.d0 * stnxz**2.d0)
      stnrat = stn2/stn1
      if (stn1 .eq. 0.d0) stnrat = 0.d0

      stateNew(k, 1) = 1.d0           ! initiation check
      stateNew(k, 2) = 0.d0           ! plastic strain
      stateNew(k, 3) = 0.d0           ! plastic strain rate 
      stateNew(k, 4) = sigy0          ! yield stress
      stateNew(k, 5) = f0             ! porosity
      stateNew(k, 6) = -1.d0          ! Gurson potential
      stateNew(k, 7) = 0.d0           ! initial increment of pl. strain 
      stateNew(k, 8) = 0.d0           ! number of iterations
      stateNew(k, 9) = sigeff         ! eff. max. stress
      stateNew(k, 10) = 1.d0          ! element damage
      stateNew(k, 11) = stnx          ! total strain, exx
      stateNew(k, 12) = stnz          ! total strain, eyy
      stateNew(k, 13) = stnxz         ! total strain, exy
      stateNew(k, 14) = stn1          ! principal strain, stn1 >= stn2
      stateNew(k, 15) = stn2          ! principal strain 2
      stateNew(k, 16) = 0.d0          ! principal strain rate 1
      stateNew(k, 17) = 0.d0          ! principal strain rate 2
      stateNew(k, 18) = 0.d0          ! triaxiality
      !tempNew(k) = tempIn
      go to 30
C-------------------------- 2nd step and later -------------------------
 20   epbar = stateOld(k, 2)
      eprate = 0.d0
      sigyT = max(sigy0, stateOld(k, 4))
      fs = stateOld(k, 5)
      tempT = tempOld(k)
      sOld(1:6) = stressOld(k, 1:6)
      trInc = sum(strainInc(k, 1:3))

      do k1 = 1, 3
        dsig(k1) = alamda*trInc + 2.d0*mu*strainInc(k, k1)
        dsig(k1+3) = 2.d0*mu*strainInc(k, k1+3)
      end do 
      sT(1:6) = sOld(1:6) + dsig(1:6)
      sigeff = sigT + sum( sOld(1:3) )/3.d0

      call fnsigcomp(sighyd, sigeqv, sigdev, sT)
      call fnfs(fT, fs, fc, ff, q1)

C values
      np = 0.d0
      sNew = sT
      epbarN = epbar
      eprateN = 0.d0
      sigyN = sigT
      fN = fs
      tempN = tempT
      GurN = stateOld(k, 6)

      Niter = 40.                  ! max number of N-R iterations
      pi = 2.d0 * asin(1.d0)
      sq2pi = sqrt(pi)
      iter = 0
      iterC = 0
      Pep = 0.d0

C plastic flow vector
      np = 3.d0 * sigdev/(sigyT*sigyT)
      do k1 = 1, 3
        np(k1)= np(k1)+ 3.d0*fs*q1*q2* sinh(1.5d0*sighyd*q2/sigyT)/sigyT
      end do 
      np = np/sqrt(dot_product(np, np))
      sigdotp = dot_product(sOld, np) + dot_product(sOld(4:6), np(4:6))
      devdotp = dot_product(sigdev,np)+dot_product(sigdev(4:6),np(4:6))
      do k1 = 1, 6
        do k2 = 1, 6
          Pep(k1) = Pep(k1) + L(k1, k2) * np(k2)
        end do 
      end do 

      if (dot_product(sOld,sOld) .lt. 1.d-10) go to 27

      ep = 0.d0
      epmin = 0.d0
      epmax = 1.d-5 !max(1.d-5, 3.d0*1.d-5)

C NR - iteration starts
      do while (iter .lt. Niter) 
        iter = iter + 1
        if (iter .gt. Niter-1.) then
          print*, 'too many iterations, iter = ', iter
          print*, 'Gurson = ', GurN
          print*, 'Porosity, epmax ', fN, epmax
          call XPLB_EXIT
        end if
        
        epbarN = epbar + ep
        eprateN = ep/dt
        tempN = tempT
        call fnsigy(sigyN, sigy0, epbarN, epbar0, eprateN, eprate0, mpw,
     1    npw, bg, cc, tempN, temp0)
        if (sigyN .lt. sigyT) sigyN = sigyT
        sigyrate = (sigyN - sigyT)/dt

C porosity distribution
        Astrn = 0.d0
        Astrs = 0.d0
        Bstrs = 0.d0
        isnuc = 0.d0
        arg = (epbarN - fnmean)/fndev
        if (arg .gt. 10.d0) arg = 10.d0
        Astrn = fnamp*exp(-0.5d0*( arg )**2)/ (fndev*sq2pi)

        if (sigeff .gt. stateOld(k, 9)) isnuc = 1.d0
      
        if (isnuc .lt. 1.d0) go to 26
        arg = (sigyT + sighyd - fsmean)/fsdev
        if (arg .gt. 10.d0) arg = 10.d0
        Bstrs = fsamp*exp(-0.5d0*( arg )**2)/ (fsdev*sq2pi)
 26     Astrs = Bstrs

C Nahson-Hachinson parameter for shear
        J3 = sigdev(1)*sigdev(2)*sigdev(3) + 2.d0*sigdev(4)*sigdev(5)*
     1    sigdev(6) - sigdev(1)*sigdev(4)*sigdev(4) - sigdev(2)*
     2    sigdev(5)*sigdev(5) - sigdev(3)*sigdev(6)*sigdev(6)
        Wnh = 1.d0 - (13.5d0 * J3 / sigeqv**3.d0)**2.d0

C updating the parameters
        fdot = (1.d0 - fs)*ep*sum(np(1:3))/dt + Astrn * eprateN
     1    + Bstrs * sigyrate + kw*fs*Wnh*devdotp*eprateN/sigeqv
        if (fdot .lt. 0.d0) fdot = 0.d0
        fN = fs + fdot*dt
        sNew = sT - Pep*ep * (1.d0 - fN)*sigyN/sigdotp
        sigHrate = sum( sNew(1:3) - sOld(1:3) )/dt
        if (sigHrate .lt. 0.d0) sigHrate = 0.d0

        call fnfs(fNs, fN, fc, ff, q1)
        call fnsigcomp(sighyd, sigeqv, sigdev, sNew)
        call fngur(GurN, sighyd, sigeqv, sigyN, fNs, q1, q2)
        if (abs(GurN) .lt. tol) exit

C elastic criteria, ep = 0 & Gur < 0
        if ((ep .eq. 0.d0) .and. (GurN .lt. 0.d0)) go to 27

C rearranging the margins and new plastic strain increm.
        if ((GurN .ge. 0.d0) .and. (ep .ge. epmin)) epmin = ep
        if ((GurN .lt. 0.d0) .and. (ep .lt. epmax)) epmax = ep
        ep = 0.5d0 * (epmax + epmin)
        if ((iter .eq. 1) .and. (iterC .eq. 0)) ep = stateOld(k, 7)

C extending the plastic strain maxima.
        if ((iter .gt. Niter-2) .and. ((epmax-ep) .lt. 1.d-8)) then
          epmin = epmax
          epmax = 10.d0 * epmax
          ep = 0.5d0 * (epmin + epmax)
          iter = 0
          iterC = iterC + 1
          if (epmax .gt. 1.d-2) then
            print*, 'epmax, Gur', epmax, GurN
            print*, 'reduce the strain rate'
            call XPLB_EXIT
          end if
        end if

      end do 
      
C evaluating the min/max. principal strain (planar)
 27   fN = fN  + Bstrs * sigHrate/3.d0 
      stnx = stateOld(k, 11) + strainInc(k, 1)
      stnz = stateOld(k, 12) + strainInc(k, 3)
      stnxz = stateOld(k, 13) + strainInc(k, 5)
      stn1 = 0.50*(stnx + stnz) + 0.5d0*sqrt((stnx - stnz)**2.d0 
     1  + 4.d0 * stnxz**2.d0)
      stn2 = 0.50*(stnx + stnz) - 0.5d0*sqrt((stnx - stnz)**2.d0 
     1  + 4.d0 * stnxz**2.d0)
      stnrat = stn2/stn1
      if (stn1 .eq. 0.d0) stnrat = 0.d0

C work, plastic work and temp
      dWork = dot_product( 0.5d0*(sOld(1:6) + sNew(1:6)),
     1  strainInc(k, 1:6))
      dPwork = 0.5d0 * ep * sigeqv
      enerInternNew(k) = enerInternOld(k) + dWork/dens
      enerInelasNew(k) = enerInelasOld(k) + dPwork/dens
      !tempNew(k) = tempN + chi * sigdotp * ep / (dens*Cp)

C updating state variables
      stressNew(k, 1:6) = sNew(1:6)
      stateNew(k, 1) = 1.d0
      stateNew(k, 2) = epbarN
      stateNew(k, 3) = eprateN
      stateNew(k, 4) = sigyN
      stateNew(k, 5) = fN
      stateNew(k, 6) = GurN
      stateNew(k, 7) = ep
      stateNew(k, 8) = iter + Niter*iterC
      stateNew(k, 9) = max(sigeff, stateOld(k, 9))
      stateNew(k, 10) = 1.d0
      stateNew(k, 11) = stnx
      stateNew(k, 12) = stnz
      stateNew(k, 13) = stnxz
      stateNew(k, 14) = stn1
      stateNew(k, 15) = stn2
      stateNew(k, 16) = (stn1 - stateOld(k, 14))/dt
      stateNew(k, 17) = (stn2 - stateOld(k, 15))/dt
      stateNew(k, 18) = -sighyd/sigeqv
      if (fN .ge. ff) stateNew(k, 10) = 0.d0          ! element deletion

 30   continue
      return
      end

C=========================== Equiv stress ==============================
      subroutine fnsigcomp(sighyd, sigeqv, sigdev, s)
      include 'vaba_param.inc'
      real*8 sighyd, sigeqv, sigdev(6), s(6)
      integer k

      sighyd = sum(s(1:3))/3.d0
      do k = 1, 3
        sigdev(k) = s(k) - sighyd
        sigdev(k+3) = s(k+3)
      end do 

      sigeqv = sqrt(3.d0/2.d0 *(sigdev(1)**2.d0 + sigdev(2)**2.d0 +
     1  sigdev(3)**2.d0 + 2.d0*sigdev(4)**2.d0 + 2.d0*sigdev(5)**2.d0 +
     2  2.d0*sigdev(6)**2.d0 ))
      end subroutine
C======================== Gurson yield function ========================
      subroutine fngur(Gur, sighyd, sigeqv, sigy, fs, q1, q2)
      include 'vaba_param.inc'
      real*8 Gur, sighyd, sigeqv, sigy, fs, q1, q2

      Gur = (sigeqv/sigy)**2 + 2.d0*q1*fs*cosh(1.5d0*q2*sighyd/sigy) - 
     1 (1.d0 + (q1*fs)**2)
      end subroutine
C============================== Porosity ===============================
      subroutine fnfs(f, fs, fc, ff, q1)
      include 'vaba_param.inc'
      real*8 f, fs, fc, ff, q1

      if (fs .le. fc) then
        f = fs
      else
        f = fc + (1.d0/q1 - fc)*(fs - fc)/(ff - fc)
      end if 
      if (f .gt. ff) f = ff
      end subroutine
C============================ Yield stress =============================
      subroutine fnsigy(sigy, sigy0, epbar, epbar0, eprate, eprate0,
     1  mpw, npw, bg, cc, temp, temp0)
      include 'vaba_param.inc'
      real*8 sigy, sigy0, epbar, epbar0, eprate, eprate0, mpw, npw,
     1  bg, cc, temp, temp0, eprsig

      eprsig = eprate/eprate0
      if (eprate .le. eprate0) eprsig = 1.d0

      sigy = sigy0*(eprsig**mpw)*((1.d0 + epbar/epbar0)**npw)*
     1 (1.d0 +bg*exp(-cc*(temp0-273.d0))*(exp(-cc*(temp-temp0)) - 1.d0))
      end subroutine
