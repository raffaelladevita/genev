

c----------------------------------------------------------------------
c                    radiate1
c----------------------------------------------------------------------
c  SIGMA_TRUE(es,ep,theta) 
c    - the function calculating non-rad xsec.
c  w,q2,omega,theta,phi  
c    - the values for measured electrons
c  w_vtx,q2_vtx,omega_vtx,theta_vtx,phi_vtx  - 
c    - the values in the (eeg) vertex
c  e_radgam, theta_radgam, phi_radgam 
c    - parameters og radiated photon (if e_radgam=0 ->soft photon)
c  mf_min, Ttarg,Twi,Twf,Ztarg,Zw, delta 
c    - ...
c
      SUBROUTINE radiate1(SIGMA_TRUE, xtotal,xpart, mfin,
     &          es,    ep,    w,    q2,    omega,    theta,    phi,
     &          es_vtx,ep_vtx,w_vtx,q2_vtx,omega_vtx,theta_vtx,phi_vtx,
     &          e_radgam, theta_radgam, phi_radgam,
     &          mf_min, Ttarg,Twi,Twf,Ztarg,Zw, delta)
      IMPLICIT none
      REAL*8    SIGMA_TRUE,sigma_rad
      EXTERNAL  SIGMA_TRUE
      REAL*4    xtotal,xpart
      REAL*4    mfin,mp,pii
      REAL*4    es,    ep,     w,    q2,    omega,    theta,    phi
      REAL*4    es_vtx,ep_vtx, w_vtx,q2_vtx,omega_vtx,theta_vtx,phi_vtx
      REAL*4    e_radgam, theta_radgam, phi_radgam
      REAL*4    mf_min, Ttarg,Twi,Twf,Ztarg,Zw, delta
      REAL*4    rnum, ar(1000), elow,ewid,eran
      INTEGER*4 i

      REAL*8             x_radtot,x_radsoft,x_radhardini,x_radhardfin
      COMMON/com_radxsec/x_radtot,x_radsoft,x_radhardini,x_radhardfin
      REAL*8                r_radsoft
      COMMON/com_radsoft/   r_radsoft
      INTEGER*4             np_radhardini,nn_radhardini
      REAL*8                int_radhardini,e_radhardini,r_radhardini
      COMMON/com_radhardini/np_radhardini,nn_radhardini,int_radhardini,
     &                      e_radhardini(1000),r_radhardini(1000)
      INTEGER*4             np_radhardfin,nn_radhardfin
      REAL*8                int_radhardfin,e_radhardfin,r_radhardfin
      COMMON/com_radhardfin/np_radhardfin,nn_radhardfin,int_radhardfin,
     &                      e_radhardfin(1000),r_radhardfin(1000)


cc      REAL*4     x_tot,x_soft,x_hardini,x_hardfin
cc      REAL*4     x_tot,x_soft,x_hardini,x_hardfin
cc      INTEGER*4  np_rhoini
cc      INTEGER*4  np_rhofin
cc      REAL*4     fcsoft
cc      REAL*4     omegalowini,omegawidini,rhoini(1000)
cc      REAL*4     omegalowfin,omegawidfin,rhofin(1000)


c ----- calculate radiative cross-section -----
c ----- RHOs are saved in commons         -----
      mp  = 0.93821
      pii = 3.1415926536
      xtotal = sigma_rad(SIGMA_TRUE, es*1.d+0,ep*1.d+0,theta*1.d+0,
     &                 mf_min*1.d+0, Ttarg*1.d+0,Twi*1.d+0,Twf*1.d+0,
     &                 Ztarg*1.d+0,Zw*1.d+0, delta*1.d+0)
      IF(xtotal.eq.0.) THEN
      xtotal = 0.
      xpart  = 0.
      e_radgam     = 0.
      theta_radgam = 0. 
      phi_radgam   = 0.
      RETURN
      ENDIF
      
      CALL ranlux(rnum,1)       !soft
         

ccc        rnum=0.5
ccc        print *,' rnum=',rnum

c ----- SOFT photon emission -----
      IF    ( rnum*1.d+0 .le. x_radsoft/x_radtot )THEN
        xpart   = x_radsoft
        e_radgam     = 0.d+0
        theta_radgam = 0.d+0
        phi_radgam   = 0.d+0
        es_vtx    = es
        ep_vtx    = ep
        theta_vtx = theta
        phi_vtx   = phi
c        print *,' SOF es_vtx,ep_vtx,e_radgam==',es_vtx,ep_vtx,e_radgam
c ----- HARD photon emission by the initial electron -----
      ELSEIF( rnum*1.d+0 .le. (x_radsoft+x_radhardini)/x_radtot )THEN
        do i=1,np_radhardini-1
        ar(i) = 0.5*(r_radhardini(i)+r_radhardini(i+1))
        if(ar(i).lt.0.) ar(i) = 0.
        enddo
        ewid = e_radhardini(2) - e_radhardini(1)
        elow = e_radhardini(1)
        CALL RNHPRE(ar,np_radhardini-1)
        CALL RNHRAN(ar,np_radhardini-1,elow,ewid,eran)
        CALL ranlux(rnum,1)
        xpart = x_radhardini
        e_radgam     = eran
        theta_radgam = 0. 
        phi_radgam   = 0. 
        es_vtx    = es-e_radgam
        ep_vtx    = ep
        theta_vtx = theta
        phi_vtx   = phi
c        print *,' INI es_vtx,ep_vtx,e_radgam==',es_vtx,ep_vtx,e_radgam
c ----- HARD photon emission by the final electron -----
      ELSE
        do i=1,np_radhardfin-1
        ar(i) = 0.5*(r_radhardfin(i)+r_radhardfin(i+1))
        if(ar(i).lt.0.) ar(i) = 0.
        enddo
        ewid = e_radhardfin(2) - e_radhardfin(1)
        elow = e_radhardfin(1)
        CALL RNHPRE(ar,np_radhardfin-1)
        CALL RNHRAN(ar,np_radhardfin-1,elow,ewid,eran)
        CALL ranlux(rnum,1)
        xpart = x_radhardfin
        e_radgam     = eran
        theta_radgam = theta
        phi_radgam   = phi
        es_vtx    = es
        ep_vtx    = ep+e_radgam
        theta_vtx = theta
        phi_vtx   = phi
cc        print *,' FIN es_vtx,ep_vtx,e_radgam==',es_vtx,ep_vtx,e_radgam
      ENDIF


c ----- calculate W_vtx,Q2_vtx,omega_vtx from Es_vtx,Ep_vtx,theta,phi
      theta_vtx = theta_vtx
      phi_vtx   = phi_vtx
      omega_vtx = Es_vtx - Ep_vtx
      q2_vtx = 4.*Es_vtx*Ep_vtx*(sin(theta_vtx/2.))**2
      w_vtx  = sqrt( max( (mp**2 + 2.*mp*omega_vtx - q2_vtx), (0.) ) )
      

      RETURN
      END


c----------------------------------------------------------------------
c                     sigma rad corrected (MoTsai) 
c----------------------------------------------------------------------
c
c  FUNCTION sigma_rad(SIGMA_TRUE, es,ep,theta, mf_min, 
c &                   Ttarg,Twi,Twf,Ztarg,Zw, delta)
c  calculates radiative cross-section from the noncorrected
c  (Born) cross-section. The ONLY corrections for the CONTINUOUS
c  spectra are applied.
c  ref: Mo.Tsai. Rev. Mod. Phys. V,41, N.1, P.205, (1969)
c  
c  SIGMA_TRUE - the function SIGMA_TRUE(es,ep,theta)
c               calculating not radiative cross-section
c  es,ep      - Es,Ep - initial and final electron energies
c  theta      - scattering angle 
c  mf_min     - The minimum of the final state mass (exept electron)
c  Ttarg,Twi,Twf - the length of the Target, initial window, and
c                  final window in units of the rad. len.
c  Ztarg,Zw   - Z of target and i/f windows 
c  delta      - Parameter Delta ~ 0.010GeV
c
      FUNCTION sigma_rad(SIGMA_TRUE, es,ep,theta,
     &                   mf_min, Ttarg,Twi,Twf,Ztarg,Zw, delta)
      IMPLICIT none
      REAL*8    sigma_rad,s1,s2,s3
      REAL*8    SIGMA_TRUE, es,ep,theta
      REAL*8    mf_min, Ttarg,Twi,Twf,Ztarg,Zw, delta
      REAL*8    FACTOR_RADSOFT,TAIL_RADHARDINI,TAIL_RADHARDFIN
      INTEGER*4 MODA
      EXTERNAL  SIGMA_TRUE,SIGMA_33_zw
      REAL*8               SIGMA_33_zw

      REAL*8             x_radtot,x_radsoft,x_radhardini,x_radhardfin
      COMMON/com_radxsec/x_radtot,x_radsoft,x_radhardini,x_radhardfin
     &                      

      MODA=1

      s1 = SIGMA_TRUE(es,ep,theta)
     &    *factor_radsoft(es,ep,theta,
     &                    Ttarg,Twi,Twf,Ztarg,Zw, delta)
      s2 = tail_radhardini(SIGMA_TRUE, es,ep,theta,
     &                     mf_min, Ttarg,Twi,Twf,Ztarg,Zw, delta)
      s3 = tail_radhardfin(SIGMA_TRUE, es,ep,theta,
     &                     mf_min, Ttarg,Twi,Twf,Ztarg,Zw, delta)
      x_radsoft    = s1
      x_radhardini = s2
      x_radhardfin = s3
      x_radtot     = s1+s2+s3
      sigma_rad    = s1+s2+s3
      
      RETURN
      END


c----------------------------------------------------------------------
c                   tail_radhardini  (MoTsai(IV.1))
c                   tail_radhardfin
c----------------------------------------------------------------------
      FUNCTION tail_radhardini(SIGMA_TRUE, es,ep,theta, mf_min, 
     &                         Ttarg,Twi,Twf,Ztarg,Zw, delta)
      IMPLICIT none
      REAL*8    tail_radhardini
      REAL*8    SIGMA_TRUE, es,ep,theta
      REAL*8    mf_min, Ttarg,Twi,Twf,Ztarg,Zw, delta
      EXTERNAL  SIGMA_TRUE
      INTEGER*4 nn,i
      REAL*8    omegas, omegas_min, omegas_max
      INTEGER*4             np_radhardini,nn_radhardini
      REAL*8                int_radhardini,e_radhardini,r_radhardini
      COMMON/com_radhardini/np_radhardini,nn_radhardini,int_radhardini,
     &                      e_radhardini(1000),r_radhardini(1000)
      REAL*8 DSIMPS
      REAL*8 OMEGA_MAXINI,RHOS_RADHARD

      nn = 801
      np_radhardini = nn
      omegas_min = delta
      omegas_max = OMEGA_MAXINI(es,ep,theta, mf_min)
      IF(omegas_max.gt.omegas_min) THEN
      DO i=1,nn
        omegas  = omegas_min+((omegas_max-omegas_min)*(i-1))/(nn-1)
        e_radhardini(i) = omegas
        r_radhardini(i) = RHOS_RADHARD(es,ep,theta, omegas,
     &                      mf_min, Ttarg,Twi,Twf,Ztarg,Zw, delta)
     &                   *sigma_true(es-omegas,ep,theta)
      ENDDO
      tail_radhardini = DSIMPS(r_radhardini, omegas_min,omegas_max,nn-1)
      int_radhardini  = tail_radhardini
      ELSE
      tail_radhardini = 0.d+0
      int_radhardini  = tail_radhardini
      ENDIF

      RETURN
      END


      FUNCTION tail_radhardfin(SIGMA_TRUE, es,ep,theta, mf_min, 
     &                         Ttarg,Twi,Twf,Ztarg,Zw, delta)
      IMPLICIT none
      REAL*8   tail_radhardfin
      REAL*8   SIGMA_TRUE, es,ep,theta
      REAL*8   mf_min, Ttarg,Twi,Twf,Ztarg,Zw, delta
      EXTERNAL SIGMA_TRUE
      INTEGER*4 nn,i
      REAL*8    omegap, omegap_min, omegap_max
      INTEGER*4             np_radhardfin,nn_radhardfin
      REAL*8                int_radhardfin,e_radhardfin,r_radhardfin
      COMMON/com_radhardfin/np_radhardfin,nn_radhardfin,int_radhardfin,
     &                      e_radhardfin(1000),r_radhardfin(1000)
      REAL*8 DSIMPS
      REAL*8 OMEGA_MAXFIN,RHOP_RADHARD

      nn = 801
      np_radhardfin = nn
      omegap_min = delta
      omegap_max = OMEGA_MAXFIN(es,ep,theta, mf_min)
      IF(omegap_max.gt.omegap_min) THEN
      DO i=1,nn
        omegap = omegap_min+((omegap_max-omegap_min)*(i-1))/(nn-1)
        e_radhardfin(i) = omegap
        r_radhardfin(i) = RHOP_RADHARD(es,ep,theta, omegap,
     &                      mf_min, Ttarg,Twi,Twf,Ztarg,Zw, delta)
     &                   *sigma_true(es,ep+omegap,theta)
      ENDDO
      tail_radhardfin = DSIMPS(r_radhardfin, omegap_min,omegap_max,nn-1)
      ELSE      
      tail_radhardfin=0.d+0
      ENDIF

      RETURN
      END


c----------------------------------------------------------------------
c                   omega_maxini 
c                   omega_maxfin 
c----------------------------------------------------------------------
      FUNCTION OMEGA_MAXINI(es,ep,theta, mf_min)
      IMPLICIT none
      REAL*8   OMEGA_MAXINI,es,ep,theta, mf_min, es_min,mp
      mp  = 0.938
      es_min = (mf_min**2 + 2.*mp*mf_min + 2.*mp*ep) / 
     &         (2.*mp - ep*(2.*(1.-cos(theta))))
      omega_maxini = es - es_min
      RETURN
      END
      FUNCTION OMEGA_MAXFIN(es,ep,theta, mf_min)
      IMPLICIT none
      REAL*8   OMEGA_MAXFIN,es,ep,theta, mf_min, ep_max,mp
      mp  = 0.938
      ep_max = (2.*mp*es - 2.*mp*mf_min - mf_min**2) /
     &         (2.*mp + es*(2.*(1.-cos(theta))))
      omega_maxfin = ep_max - ep
      RETURN
      END


c----------------------------------------------------------------------
c                   rhos_radhard 
c                   rhop_radhard 
c----------------------------------------------------------------------
      FUNCTION RHOS_RADHARD(es,ep,theta, omegas,
     &                      mf_min, Ttarg,Twi,Twf,Ztarg,Zw, delta)
      IMPLICIT none
      INTEGER*4 Isave
      REAL*8   rhos_radhard
      REAL*8   es,ep,theta,sp
      REAL*8   mf_min, Ttarg,Twi,Twf,Ztarg,Zw, Btarg,Bw, delta
      REAL*8    me,me2,mp,alfa,pi,alfadpi
      REAL*8    tr, c0,cs, fs,fp,  omegas, esprime,xs,ts
      SAVE   me,me2,alfa,pi,alfadpi,Btarg,Bw
      REAL*8 BZFUNMOTS
      DATA Isave/0/

      IF(Isave.eq.0) THEN
      Isave=1
      me  = 0.000511
      me2 = me**2
      mp  = 0.938
      alfa= 1.0/137.035
      pi  = 3.1415926535898
      alfadpi = alfa/pi
      Btarg = BZFUNMOTS(Ztarg)
      Bw    = BZFUNMOTS(Zw)
      ENDIF
        esprime = es - omegas
        sp  = Es*Ep*(1.-cos(theta))
        c0  = log(2*sp/me2)
        cs = (Bw*Twi+0.5*Btarg*Ttarg)
        tr  = (alfa/pi/Btarg) * (log(2*sp/me2) - 1.)
        fs = Btarg*tr + Bw*Twi + 0.5*Btarg*Ttarg
        fp = Btarg*tr + Bw*Twf + 0.5*Btarg*Ttarg
        xs = esprime/es
        ts = alfadpi*( 0.5*(1+xs**2)*c0 - xs )
        rhos_radhard = (1./(es-esprime))
     &       * ( ts + cs*(xs+0.75*(1-xs)**2) )
     &       * ( (log((1./xs)))**fs )
     &       * (delta/Ep)**(0.5*fp)

      RETURN
      END


      FUNCTION RHOP_RADHARD(es,ep,theta, omegap,
     &                      mf_min, Ttarg,Twi,Twf,Ztarg,Zw, delta)
      IMPLICIT none
      INTEGER*4 Isave
      REAL*8   rhop_radhard
      REAL*8   es,ep,theta,sp
      REAL*8   mf_min, Ttarg,Twi,Twf,Ztarg,Zw, Btarg,Bw, delta
      REAL*8    me,me2,mp,alfa,pi,alfadpi
      REAL*8    tr, c0,cp, fs,fp,  omegap, epprime,xp,tp
      SAVE   me,me2,alfa,pi,alfadpi,Btarg,Bw
      REAL*8 BZFUNMOTS
      DATA Isave/0/

      IF(Isave.eq.0) THEN
      Isave=1
      me  = 0.000511
      me2 = me**2
      mp  = 0.938
      alfa= 1.0/137.035
      pi  = 3.1415926535898
      alfadpi = alfa/pi
      Btarg = BZFUNMOTS(Ztarg)
      Bw    = BZFUNMOTS(Zw)
      ENDIF
        epprime = ep + omegap
        sp  = Es*Ep*(1.-cos(theta))
        c0  = log(2*sp/me2)
        cp = (Bw*Twf+0.5*Btarg*Ttarg)
        tr  = (alfa/pi/Btarg) * (log(2*sp/me2) - 1.)
        fs = Btarg*tr + Bw*Twi + 0.5*Btarg*Ttarg
        fp = Btarg*tr + Bw*Twf + 0.5*Btarg*Ttarg
        xp = ep/epprime
        tp = alfadpi*( 0.5*(1+xp**2)*c0 - xp )
        rhop_radhard = (1./(epprime-ep))
     &       * ( tp + cp*(xp+0.75*(1-xp)**2) )
     &       * ( (log((1./xp)))**fp )
     &       * (delta/Es)**(0.5*fs)

      RETURN
      END



c----------------------------------------------------------------------
c                     factor_radsoft(MoTsai )
c----------------------------------------------------------------------
      FUNCTION factor_radsoft(Es,Ep,theta, 
     &                        Ttarg,Twi,Twf,Ztarg,Zw, delta)
      IMPLICIT none
      REAL*8 factor_radsoft, Es,Ep,theta, sp
      REAL*8 Ttarg,Twi,Twf,Ztarg,Zw,Btarg,Bw, delta
      REAL*8 me,me2,alfa,pi,delta_t,delta_r
      REAL*8 SPENCEFU,BZFUNMOTS

      REAL*8             r_radsoft
      COMMON/com_radsoft/r_radsoft

      me  =0.000511
      me2 =me**2
      alfa=1.0/137.035
      pi  =3.1415926535898
      sp= Es*Ep*(1. - cos(theta))

      Btarg = BZFUNMOTS(Ztarg)
      Bw    = BZFUNMOTS(Zw)
c     MoTsai pg. 214
      delta_t = (Bw*Twi+0.5*Btarg*Ttarg)*log(Es/delta) +
     &          (Bw*Twf+0.5*Btarg*Ttarg)*log(Ep/delta)
      delta_t = -delta_t

c     MoTsai pg. 214
      delta_r = 28.d+0/9.d+0 - 
     &          (13.d+0/6.d+0)*log(2.*sp/me2) +
     &          (log(es/delta)+log(ep/delta))*(log(2*sp/me2)-1.) -
     &          spencefu( (ep-es)/ep ) -
     &          spencefu( (es-ep)/es )
      delta_r = delta_r * (-(alfa/pi))

      r_radsoft      = exp( delta_t + delta_r )
      factor_radsoft = r_radsoft

      RETURN
      END


c----------------------------------------------------------------------
c                           BZFUNMOTS
c----------------------------------------------------------------------
      FUNCTION BZFUNMOTS(Z)
      IMPLICIT none
      REAL*8  z,xi, BZFUNMOTS
      xi=(log(1440.d+0)-2.0*log(z)/3.0) / (log(183.d+0)-log(z)/3.0)
      BZFUNMOTS = (4.d+0/3.d+0) *
     &            ( 1. + ((z+1.)/(z+xi))/(log(183.d+0)-log(z)/3.)/9. )
      RETURN
      END


c----------------------------------------------------------------------
c                           spence function
c----------------------------------------------------------------------
      FUNCTION spencefu(x)
      IMPLICIT none
      REAL*8 spencefu,x,eps,dgauss,null,pii,spence1
      REAL*8   spencefunoi
      EXTERNAL spencefunoi
      pii=3.1415926535898
      spence1 = (pii**2)/6.d+0
      null=1.D-8
      eps=1.D-8
      IF    (x.ge.null.and.x.le.1) THEN
        spencefu=           dgauss(spencefunoi,null,    x,eps)
      ELSEIF(x.ge.1) THEN
        spencefu= spence1 + dgauss(spencefunoi,1.d+0,   x,eps)
      ELSEIF(x.le.-1.*null) THEN
        spencefu=          -dgauss(spencefunoi,x,-1.*null,eps)
      ELSE
        spencefu= 0.d+0
      ENDIF
      RETURN         
      END  
      FUNCTION spencefunoi(x)
      IMPLICIT none
      REAL*8 spencefunoi,x
      spencefunoi = -log(abs(1.0d+0-x))/x
      RETURN         
      END
cccccccccccccccccccccccc TEST SPENCE FUNC cccccccccccccccccccccccc
c      do i=1,1000
c        x=0. + 1.*i/1000.
c        f1=spencefu(x)
c        f2=spencefu(-x)   
c        write(11,*) x,f1,f2
c      enddo
c      do i=1,1000
c        x=1. + 99.*i/1000.
c        f1=spencefu(x)
c        f2=spencefu(-x)   
c        write(11,*) x,f1,f2
c      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c----------------------------------------------------------------------
c                     sigma_33 (MoTsai(III.5-9)) 
c----------------------------------------------------------------------
      FUNCTION sigma_33(es,ep,theta)
      IMPLICIT none
      REAL*8 sigma_33,sigma_33_zw, es,ep, theta, thetad2
      REAL*8 alfa,pi,me,mp,mp2,m33,m332,mpi,mpi2, mim
      REAL*8 q2,sq2,mf,mf2
      REAL*8 qq2,qq2star,eistar,c32,p2star,pstar,x,gmf,g1,g2,gg,ff
      LOGICAL zerow

      zerow=.false.
      GOTO 1001
      ENTRY sigma_33_zw(es,ep,theta)
      zerow=.true.
 1001 CONTINUE

      alfa = 7.2974d-3
      pi   = 3.1415926535898
      me  = 0.000511d+0
      mp  = 0.938d+0
      mp2 = mp**2
      m33 = 1.236d+0
      m332= m33**2
      mpi = 0.1396d+0
      mpi2= mpi**2
      thetad2 = theta/2.
      q2   = -( 4.*Es*Ep*(sin(thetad2))**2 )
      if(q2.gt.0.) sigma_33=0.d+0
      if(q2.gt.0.) RETURN
      sq2  = sqrt(-q2)
      mf2  = mp**2 + 2.*mp*(Es-Ep) + q2
      if(mf2.lt.0.) sigma_33=0.d+0
      if(mf2.lt.0.) RETURN
      if(mf2.lt.0.) print *,' sigma_true: mf2 < 0'
      if(mf2.lt.0.) print *,' Es,Ep,thetad=',Es,Ep,theta
      if(mf2.lt.0.) print *,' q2=',q2
      if(mf2.lt.0.) stop
      mf   = sqrt(mf2)

      qq2     = ((mf2 - q2 - mp2)/(2.*mp))**2 - q2 
      qq2star = mp2*qq2/mf2
      eistar  = (mf2+mp2-q2)/(2.*mf)
      c32     = (2.05**2)*exp((-6.3)*sq2)*(1.+9.*sq2) / mp2
      p2star  = ( (mf2-mp2+mpi2)/(2.*mf) )**2 - mpi2
      if(p2star.lt.0.)  sigma_33    = 0.d+0
      if(p2star.lt.0.)  sigma_33_zw = 0.d+0
      if(p2star.lt.0.)  RETURN
      if(p2star.lt.0.)  print *,' sigma_true: p2star < 0'
      if(p2star.lt.0.)  stop
      pstar   = sqrt(p2star)
      x = (0.85*pstar/mpi)
      gmf     = 0.1293*( x**3 / (1.+x**2) )
      gmf     = 0.0020*( x**3 / (1.+x**2) )

      IF( zerow ) THEN
        mim=sqrt((es-ep+mp)**2 - (es**2+ep**2-2.*es*ep*cos(theta)))
        if(abs(mim-m33).gt.1.d-6) print *,' sigma_33_zw: bad mim=',mim
        if(abs(mim-m33).gt.1.d-6) stop
        gg = (4.d+0/3.d+0)*m33*(eistar+mp)*qq2star*c32
        ff = gg*((-q2)/qq2)/mp2
        sigma_33 = (alfa**2*ep**2/(q2**2))
     &           * ( 1. / ( 1. + (Es/mp)*(1.-cos(theta)) ) )
     &           * ((cos(thetad2))**2)
     &           * ( ff + (2./mp2)*((tan(thetad2))**2)*gg )
        sigma_33 = sigma_33*389.379
        sigma_33_zw = sigma_33
        RETURN
      ENDIF

      g1 = (gmf*m33*mf/pi) / ((mf2-m332)**2 + gmf**2*m332)
      g1 = g1 * qq2star*2.*c32
      g1 = g1 * (eistar + mp)/(3.*mp)
      g2 = g1*(-q2)/qq2

      ff = (2./mp)*g2
      gg = 2.*mp*g1

      sigma_33 = (2.*alfa**2*ep**2*mp/(q2**2)) * ((cos(thetad2))**2)
     &           * ( ff + (2./mp2)*((tan(thetad2))**2)*gg )
      sigma_33 = sigma_33*389.379

      RETURN
      END





