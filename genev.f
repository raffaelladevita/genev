c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c Input files:
c  genev.inp		parameters and name of files
c  genev_chan.mtx	xsect differential for all channels
c  genev_rhodec.out	decaing pions distribution for rho0 decay
c                       in helicity rest frame
c
c Output files:
c  xxx.lund             name from input file; if needed
c  XXX.BOS	        BOS output file
c  genev.hbook	        HBOOK output file
c
c To be tested:
c
c Warnings:
c 1) Never checked for neutron target
c 2) .. 
c 3) ..
c
c Comments:
c It's running only on HP-UX and  SUN platforms
c On LINUX RAN intrinsic function is not defined
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   pgm modifie pour ecrire 2 fichiers bos (7.3.2003)
c   les photons ne sont ecrits que si e>0.01
c   acceptance individuelle calculee que pour e,p,pi+,pi-    
c   sous la ligne 99, modif possible pour le nombre d'evts Ok ou evts generes.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c  Values of target parameters (ie Ttarg, Twi, Twf, Ztarg, Zw) set for e1-6 target. 
c  (after "CROSS-SECTION (e,e',gamma_virt)").
c
c  Track from generated channel number now kept. Info put in MCTK bank-entry 11. 
c  (see after "Writing BOS file" : iw(imctk+11)=JCH)
c   
c  Subroutine GPhiLH to calculate the angle between hadronic and leptonic planes added.
c  Call to GPhiLH to get gphi added (after "Filling PLAB1"). 
c  For that purpose, add declarations : REAL Plab_Eout(3), Plab_Pout(3), res, gphi
c
c  Rho decay : 
c  - Call to READ_W_DEC_RHO subroutine commented (after "Reading rho decay distribution").
c  - Subroutine RHO_DECAY modified to calculate rho decay angular distributions from Schilling formula.
c    As a consequence, new arguments GPhi and EPS(epsilon) have been introduced (needed in the formula). 
c    New arguments Q2 and W have also been introduced to be used in the future 
c    in case the user want to vary spin density matrix elements as a function of these two variables.
c
c  Omega decay :
c  - Subroutine ROTATION added (private communication from Marco Battaglieri).
c  - Subroutine OMEGA_DECAY added 
c  -> idem as RHO_DECAY subroutine + theta_qp and p variables fixed according to Stevenson paper 
c                           (private communication from Marco Battaglieri).
c  - If JCH=19 : call to OMEGA_DECAY added 
c    + rotation of the 3 pions to come back in the lab frame.
c   (after "Omega decay using angular distribution function").
c  - Note that to get realistic omega decay (not flat !) you have to set FLAG_Decay_RHO variable 
c    in genev.inp file to 1 !!! 
c
c  Ludyvine Morand (March 2003)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


c+ Decomment this line to run on INFNGE
C	OPTIONS/ EXTEND_SOURCE
c-
c****************************************
	PROGRAM GENEV
c****************************************

	IMPLICIT NONE
	INCLUDE 'bcs.inc'
	REAL RRAN
	REAL CURRENT,ACC_P,ACC_PIP,ACC_PIM,ACC_EL,acc_had,acc_pro
	INTEGER SE_P,SE_PIP,SE_PIM,SE_EL,SPACODE
        REAL IMP,IMPQ,THP,PHP
        REAL BETD(4),RHO_CM(5),Q2MIN,Q2MAX
        INTEGER FLAGD
	INTEGER*4 IZ,IERR
        REAL FLAGC,E_test
        REAL BETM
	INTEGER MAXPAGES, NPARZ, NTOT,jj
	REAL MP,ME
	REAL TETMIN,TETMAX,OMEGAMIN,OMEGAMAX,EO,EO_BEAM,DUMMY,DEGRAD
	REAL RADDEG,F_O_MIN,F_O_MAX,WMIN,WMAX,FIEMIN,FIEMAX
	REAL PIGR,COSMIN,COSMAX,TETE,OMEGA
	REAL EE,PE,Q2,PSI,ARG,W,EPS,GAMMA_V,D_SIGMA,D_SIGMA_PART
        REAL E_TOTAL,X_TOTAL,Y_TOTAL,Z_TOTAL,CC
	REAL D_S,R_D_S,FE
	REAL PHI,UGV,VGV,WGV,UE1,VE1,WE1,UE2,VE2,WE2,UY,VY,WY,YYMOD
	REAL UX,VX,WX,BET,GAMMA,Y,XXMOD,TETA,U,FI,V
	REAL RMASS,PROD_SCAL
	INTEGER IPAWC,LRECL,N_TOT,ISTAT,NSTORY,FIFLAG,J,LUND,MAXCH
	INTEGER K, NVAR,INW,INTH,I,JCH,M,JK,IL,L,IT,IG,ICYCLE
	INTEGER NP_CODE,NP,I_Q2,I_W
	REAL DNUP(54)
        REAL PLAB1(5,8),PLAB2(5,8),P(4),PA1(3)
	REAL BETA(4)
	REAL PA(4),CD(3)
	REAL ECM,AMASS(18),KGENEV,PCM(5,18),WGT 
	INTEGER NCH,ICHAN(24)
	REAL SIGR(40,251),WR(251)
	REAL XSECT(14,50,0:18)
	CHARACTER*12 PARTIC(10)
	CHARACTER*12 TARGET
c	CHARACTER*12  FILELUND,FILEBOS
	CHARACTER*12  FILELUND
        CHARACTER*12 FILEBOS
	CHARACTER*16 FILEBOS_acc
	CHARACTER*8 NCODE(54),file
	REAL W_DEC_RHO(21,41)
	REAL rho_cm_r(3)
	REAL ptr(5),beta_c,gamma_c

	INTEGER PART_ID(40)
	INTEGER FLAG_DS_RHO,I_DUMMY,NP_FINAL
	real allchn(251),m_allchn,tagliow,taglioq2,tagliot,taglioo,costagliot
	real voltot,vol(4),wwmax,wwmin,qqmin,qqmax,smax,comin,comax,omin,omax
	real s_max(4),a,proba(0:4),f_allchn,sig
	INTEGER NW,IBID,M0,M8,M1,J2,MEVT,RUN,IHEAD,NR,IMCEV,IMCVX,IMCTK,
     %          ISTATUS1,NBANK,IBID_acc
	INTEGER J_J,K_K,ID_BOS
	REAL CHARGE,P_REST(5)
	INTEGER I_CM,flag_decay_rho
	CHARACTER*8 NCODE1(73)		
	REAL ENUP(100)
	INTEGER J_P,J_PP,J_PM,J_PZ
	REAL SECTOR,px_spool,py_spool
	CHARACTER*120 genev_parms, filename
	
	REAL TGRAD,TGLGT,TGOFFZ
	REAL vert_r,vert_phi,vert_x,vert_y,vert_z

        REAL eo_vtx,ee_vtx,w_vtx,q2_vtx,omega_vtx,tete_vtx,fe_vtx
        REAL eo_ini,ee_ini,w_ini,q2_ini,omega_ini,tete_ini,fe_ini
        REAL e_radgam, theta_radgam, phi_radgam
        REAL mf_min, mfin, Ttarg,Twi,Twf,Ztarg,Zw, delta

        REAL*8   SIGMA_EP_EET
        EXTERNAL SIGMA_EP_EET
        INTEGER*4 ModaRad,num_soft,num_hardini,num_hardfin

        REAL Plab_Eout(3), Plab_Pout(3), res, gphi

C--------------------------------------
	PARAMETER ( MAXPAGES=1000, NPARZ=10000, NTOT=100000 )
	PARAMETER ( genev_parms = 'GENOVA_PARMS')
C--------------------------------------
	COMMON/GENIN/NP,ECM,AMASS,KGENEV
	COMMON/GENOUT/PCM,WGT 
	COMMON/CHANNELS/NCH,ICHAN
        COMMON/SEZURTO/SIGR,WR
        COMMON/SIGMAS/XSECT
	COMMON/PAWC/IPAWC(MAXPAGES*128)
	COMMON/SIGMA_MAX/TETMIN,TETMAX,OMEGAMIN,OMEGAMAX
        COMMON/SIGMA_NEW/Q2MAX,Q2MIN,WMAX,WMIN
	COMMON/RHO/W_DEC_RHO
 	COMMON/RANDOM/IZ
        common/somma/allchn
	common/costanti/Eo_beam,Mp,me,pigr,degrad,raddeg
	COMMON/LUND/RHO_CM,DNUP,UE1,VE1,WE1,PE,EE,JCH,
     @              PART_ID

C--------------------------------------
	DATA DUMMY/0/
	DATA LRECL/8191/
c define buffersize for CWN
      	call hbset ('BSIZE',LRECL,IERR)	

	DATA NCODE  /
C    ELECTRON SDR = LAB
     +	             'PXE ', 'PYE' , 'PZE' , 'EE' , 'PE',
     +               'THE',  'PHE' ,
C    1-ST PARTICLE SDR=LAB
     +	             'PX1' , 'PY1' , 'PZ1' , 'E1' , 'P1',
     +	             'TH1',  'PH1',
C    2-ND PARTICLE SDR=LAB
     +	             'PX2' , 'PY2' , 'PZ2' , 'E2' , 'P2',
     +               'TH2',  'PH2',
C    3-RD PARTICLE SDR=LAB
     +	             'PX3 ', 'PY3' , 'PZ3' , 'E3' , 'P3',
     +               'TH3',  'PH3',
C    4-TH PARTICLE SDR=LAB
     +	             'PX4 ', 'PY4' , 'PZ4' , 'E4' , 'P4',
     +               'TH4',  'PH4',
C    5-TH PARTICLE SDR=LAB
     +	             'PX5 ', 'PY5' , 'PZ5' , 'E5' , 'P5',
     +               'TH5',  'PH5',
C    6-TH PARTICLE SDR=LAB
     +	             'PX6 ', 'PY6' , 'PZ6' , 'E6' , 'P6',
     +               'TH6',  'PH6',
C    Q2, OMEGA AND CHANNEL ID
     +		     'Q2'  , 'OMEGA', 'CHAN', 'W1P', 'W2P' 
     +               /

	DATA NCODE1 /
C Electron in LAB
     +	             'e_e', 'p_e' , 'th_e' , 'phi_e',

C Proton in LAB
     +	             'e_prot', 'p_prot' , 'th_prot' , 'phi_prot',
C Neutron in LAB
     +	             'e_neut', 'p_neut' , 'th_neut' , 'phi_neut',
C First pi+
     +	             'e_pp1', 'p_pp1' , 'th_pp1' , 'phi_pp1',
C First pi-
     +	             'e_pm1', 'p_pm1' , 'th_pm1' , 'phi_pm1',
C First pi0
     +	             'e_pz1', 'p_pz1' , 'th_pz1' , 'phi_pz1',
C Second pi+
     +	             'e_pp2', 'p_pp2' , 'th_pp2' , 'phi_pp2',
C Second pi-
     +	             'e_pm2', 'p_pm2' , 'th_pm2' , 'phi_pm2',
C Second pi0
     +	             'e_pz2', 'p_pz2' , 'th_pz2' , 'phi_pz2',
C Third pi+
     +	             'e_pp3', 'p_pp3' , 'th_pp3' , 'phi_pp3',
C Third pi-
     +	             'e_pm3', 'p_pm3' , 'th_pm3' , 'phi_pm3',
C Third pi0
     +	             'e_pz3', 'p_pz3' , 'th_pz3' , 'phi_pz3',
C Pion-Delta-Rho-Omega vars in CM
     +	             'e_cm', 'p_cm' , 'th_cm' , 'phi_cm',
C Pion decaing from Rho in REST frame
     +	             'e_rest', 'p_rest' , 'th_rest' , 'phi_rest',
C First k1
     +	             'e_k1', 'p_k1' , 'th_k1' , 'phi_k1',
C First k2
     +	             'e_k2', 'p_k2' , 'th_k2' , 'phi_k2',
C Q2, omega and channel Id 
     +		     'Q2'  , 'omega', 'chan', 'W_1part', 'W_2part', 
C Radiation Gamma in LAB
     +	             'e_g', 'p_g' , 'th_g' , 'phi_g' 
     + /

C--------------------------------------

	PIGR=ACOS(-1.)
	MP=0.93827231
	ME=0.51099906E-3
	N_TOT = 0
	DEGRAD=PIGR/180
	RADDEG=1/DEGRAD

c+ Calling some paw initialisation routines
	CALL HLIMIT(128*MAXPAGES)
	CALL HROPEN(1,'ESCA','genev.hbook','N',LRECL,ISTAT)
c-
c+ Calling some BOS initialization routines
c++ Decomment nexts line to compile on INFNGE
c	GOTO 911

        NW     = 11            ! number of colums in event bank
        IBID   = 12            ! BOS output device number
        IBID_acc  = 13
        M0     = 0             ! number of header bank
        M8     = 8             ! number of colums in header bank
        M1     = 1             ! number of rows in header bank
        J2     = 1             ! number of rows in event bank
        MEVT   = 1             ! event number        
        RUN    = 0             ! run number for GSIM 
        
        CALL BOS(IW,Nbcs)
        CALL BKFMT("HEAD","I")
        CALL BKFMT('MCEV','I')
        CALL BKFMT('MCTK','(6F,5I)')
        call bkfmt('MCVX','(4F,I)')     ! MC vertex parameters 
        CALL BLIST(IW,"E=","HEAD")
        call blist(iw,'E+','MCEV')
        call blist(iw,'E+','MCVX')
        CALL BLIST(IW,'E+','MCTK')
        
911 	CONTINUE
C--
c-
c+ Reading rho decay distribution
c	CALL READ_W_DEC_RHO 
c-

c+ Filling total xsection vector 
	CALL SEZTOT
c-
c+ Reading input parameters from external file (from current dir!)
c
	OPEN(UNIT=2,FILE='genev.inp',STATUS='OLD',FORM='FORMATTED')
c	READ(2,*) FILEBOS	 	! BOS file name
	READ(2,*) FILE	 	! BOS file name
	filebos(1:8)=file                                        ! modif dd 02.03
	filebos(9:12)='.bos'
	filebos_acc(1:8)=file
	filebos_acc(9:16)='_acc.bos'
	READ(2,*) NSTORY		! # EVENTI
	READ(2,*) EO_BEAM		! EO (GEV)
	READ(2,*) F_O_MIN		! OMEGA MINIMO COME FRAZIONE DI EO 
	READ(2,*) F_O_MAX		! OMEGA MAX COME FRAZIONE DI EO
	READ(2,*) TETMIN		! TETA_MIN ELETTRONE (DEG)
	READ(2,*) TETMAX		! TETA_MAX ELETTRONE (DEG)
        READ(2,*) I_Q2			! FLAG PER FISSARE Q2
        READ(2,*) Q2MIN                 ! VALORE INTORNO CUI FA L'ESTRAZIONE
        READ(2,*) Q2MAX                 ! VALORE INTORNO CUI FA L'ESTRAZIONE
        READ(2,*) I_W			! FLAG PER FISSARE W
        READ(2,*) WMIN                  ! VALORE MINIMO PER W
        READ(2,*) WMAX                  ! VALORE MASSIMO PER W 
	READ(2,*) FIFLAG		! SE 0 FI=RAND; SE 1 FI=(FI_MIN,FI_MAX)
	READ(2,*) FIEMIN		! FI_MIN E- (DEG)
	READ(2,*) FIEMAX		! FI-MAX E- (DEG)
	READ(2,*) NCH			! NUMERO CANALI
	READ(2,*) (ICHAN(J),J=1,NCH)	! ID DEI CANALI
	READ(2,*) LUND			! 0/1 FILE LUND NO/YES
	READ(2,*) FILELUND	        ! LUND file name
        READ(2,*) FLAG_DS_RHO           ! Flag ds Rho (0=ds/dt vera; 1=Uniform 2=flat 3=fixed angle)
        READ(2,*) FLAG_Decay_RHO        ! Flag ds Rho (0=decay flat 1=decay realistic)
	READ(2,*) CURRENT               ! Current of torus for PSEUDO-SPA
	READ(2,*) SPACODE               ! Code for extraction according to cuts
	READ(2,*) TGRAD                 ! Target radius
	READ(2,*) TGLGT                 ! Target length
	READ(2,*) TGOFFZ                ! Target offset in z
	READ(2,*) MODARAD		! Radiation moda (0=nonradiative, 1-radiative(no stagling) 2-radiative with stagling(e1 target)
	CLOSE(2)
	
	
	OPEN(UNIT=93,FILE='seed.dat',STATUS='OLD',FORM='FORMATTED')
	READ(93,*) IZ
	CLOSE(93)
	CALL RLUXGO(4,IZ,0,0)
	
c+ To ammeliorate see the gsim_bos.c routine
        CALL FPARM(
     >'OPEN UNIT=12 FILE='//FILEBOS//'  WRITE RECL=32760 '//
     >'ACTION=WRITE STATUS=NEW FORM=BINARY')
        CALL FPARM(
     >'OPEN UNIT=13 FILE='//FILEBOS_acc//'  WRITE RECL=32760 '//
     >'ACTION=WRITE STATUS=NEW FORM=BINARY')

C+ Making some conversions
	IZ=IZ*10E4+1 
c-

C+ Fixing n-tuple dim according max particle number
	MAXCH=ICHAN(1)
	DO K=1,NCH
	 IF (ICHAN(K).GT.MAXCH) MAXCH = ICHAN(K)
	ENDDO
c+ !!!! NTPLE 1  is wrong when channel from 25 are enabled
	IF(MAXCH.ge.28.and.MAXCH.le.29) THEN
		NVAR=54
		NCODE(50) =  'Q2'
		NCODE(51) =  'OMEGA'
		NCODE(52) =  'CHAN'
                NCODE(53) =  'W1P'
                NCODE(54) =  'W2P' 
	ENDIF
	IF(MAXCH.eq.27) THEN
		NVAR=47
		NCODE(43) =  'Q2'
		NCODE(44) =  'OMEGA'
		NCODE(45) =  'CHAN'
                NCODE(46) =  'W1P'
                NCODE(47) =  'W2P' 
	ENDIF
	IF(MAXCH.GE.19.AND.MAXCH.LE.26) THEN
		NVAR=40
		NCODE(36) =  'Q2'
		NCODE(37) =  'OMEGA'
		NCODE(38) =  'CHAN'
                NCODE(39) =  'W1P'
                NCODE(40) =  'W2P' 
	ENDIF
	IF(MAXCH.GE.5.AND.MAXCH.LE.18) THEN
		NVAR=33
		NCODE(29) =  'Q2'
		NCODE(30) =  'OMEGA'
		NCODE(31) =  'CHAN'
                NCODE(32) =  'W1P'
                NCODE(33) =  'W2P'  
	ENDIF
	IF(MAXCH.LE.4) THEN
		NVAR=26
		NCODE(22) =  'Q2'
		NCODE(23) =  'OMEGA'
		NCODE(24) =  'CHAN'
                NCODE(25) =  'W1P'
                NCODE(26) =  'W2P' 
	ENDIF

c-
c+ Calling another paw routine after known NVAR
	CALL HBOOKN(1,' ',NVAR,'//ESCA',NPARZ,NCODE )
	CALL HBOOKN(2,' ',73  ,'//ESCA',NPARZ,NCODE1)
c-

C	TETE E OMEGA SONO DETERMINATI DALLE RELAZIONI:
C
C	W^2=MP^2+2*MP*OMEGA-Q^2
C	Q^2=4*EO*(EO-OMEGA)*SIN^2(TETE/2)
C
C	PSI = ANGOLO POLARE DEL FOTONE VIRTUALE
C
C	COS(PSI)=(Q^2+2*EO*OMEGA)/(2*EO*SQRT(Q^2+OMEGA^2))
C
C	OMEGA=(W**2+Q2-MP**2)/(2*MP)
C	TETE=2*ASIN(SQRT(Q2/(4*EO*(EO-OMEGA))))

C+ Opening lund file
	IF(LUND.NE.0) THEN
          OPEN (UNIT=2,FILE=FILELUND,STATUS='UNKNOWN',
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
C	OPEN(UNIT=2,NAME=FILELUND,STATUS='UNKNOWN',FORM='FORMATTED')
 	ENDIF
c-
c  Get full path to the input file
	CALL revinm(genev_parms,'genev_chan.mtx',filename)
	
C+ Filling XSECT with angolar distribution for only 14 channels
          OPEN (UNIT=55,FILE=filename,STATUS='OLD')
          DO JCH=1,14 
            DO INW=1,50
              DO INTH=0,18
                READ(55,*) XSECT(JCH,INW,INTH)
              ENDDO
            ENDDO
          ENDDO
         CLOSE (UNIT=55)
c-
ccccccccccccccccccccccccccccccccccccccc


c       **************************
        call sezallchn
c       **************************
c	calculation the total cross sections of all channels


c------------------------------------------------------------------------------
c       calculation the max of the total cross section of all channels
c------------------------------------------------------------------------------
	do j=0,4
	   proba(j)=0
	enddo
	
	
	
	m_allchn=0.
        do j=1,251
         m_allchn = max(allchn(j),m_allchn)
        enddo

	 IF(I_W.EQ.1.AND.I_Q2.EQ.1) THEN
	  
	taglioW=0.5*(wmax-wmin)+wmin
	taglioq2=0.3*(q2max-q2min)+q2min

c	*************************************************************
	call massvol2(wmin,tagliow,q2min,taglioq2,s_max(1),vol(1),
     &	              m_allchn,modarad)
	call massvol2(wmin,tagliow,taglioq2,q2max,s_max(2),vol(2),
     &                m_allchn,modarad)
	call massvol2(tagliow,wmax,q2min,taglioq2,s_max(3),vol(3),
     &                m_allchn,modarad)
        call massvol2(tagliow,wmax,taglioq2,q2max,s_max(4),vol(4),
     &                m_allchn,modarad)
c	*************************************************************	

	 ELSE
       	omegamax=f_o_max*Eo_beam
	omegamin= f_o_min*Eo_beam
	cosmin=cos(tetmin*degrad)
	cosmax=cos(tetmax*degrad)

c       tagliot e' un val intermedio tra tetmin e tetmax 
c	taglioo e' un val intermedio tra omegamin e omegamax
c	servono per ottimizzare l'estrazione

	tagliot=0.37*(tetmax-tetmin)+tetmin
	taglioo=0.52*(omegamax-omegamin)+omegamin
	costagliot=cos(tagliot*degrad)

c       *************************************************************
        call massvol(tetmin,tagliot,omegamin,taglioo,s_max(1),vol(1),
     &	              m_allchn)
	call massvol(tetmin,tagliot,taglioo,omegamax,s_max(2),vol(2),
     &                m_allchn)
	call massvol(tagliot,tetmax,omegamin,taglioo,s_max(3),vol(3),
     &                m_allchn)
        call massvol(tagliot,tetmax,taglioo,omegamax,s_max(4),vol(4),
     &                m_allchn)
c       *************************************************************

	 ENDIF

c------------------------------------------------------------------------------
c       calcolo massimi nei 4 rettangoli
c------------------------------------------------------------------------------
        voltot=0.
	do j=1,4
	 voltot=voltot+vol(j)     
	enddo
ccccccccccccccccccccccccccccccccccccccc



c+ Starting loop over event number

	   acc_el=0.
	   acc_pro=0.
	   acc_pip=0.
	   acc_pim=0.
	   acc_had=0.
        num_soft=0
        num_hardini=0
        num_hardfin=0
	DO I=1,NSTORY
 99     CONTINUE
	N_TOT = N_TOT +1
c        if( mod(N_TOT,1000).eq.1 .or. N_TOT.le.3) then
c          print *,' NTOT,Istory==',N_TOT,I,num_soft,num_hardini,num_hardfin
c        endif
c        IF(N_TOT .gt. NSTORY) GOTO 9999
        EO = EO_beam

c++ Choosing between Q2,W or th_e,omega random extraction

	 IF(I_W.EQ.1.AND.I_Q2.EQ.1) THEN
ccccccccccccccccccccccccccccccccccccccc

c 	estrazione rettangolo
	a=rran()
	do j=1,4
	 proba(j)=vol(j)/voltot+proba(j-1)
	enddo
	if(a.le.proba(1)) then 
	 wwmax=tagliow
	 wwmin=wmin
	 qqmin=q2min
	 qqmax=taglioq2
	 smax=s_max(1)
	else if(a.le.proba(2)) then
	 wwmin=wmin
	 wwmax=tagliow
	 qqmin=taglioq2
	 qqmax=q2max
	 smax=s_max(2)
	else if(a.le.proba(3)) then
	 wwmin=tagliow
	 wwmax=wmax
	 qqmin=q2min
         qqmax=taglioq2
	 smax=s_max(3)
	else if(a.le.proba(4)) then
	 wwmin=tagliow
         wwmax=wmax
	 qqmin=taglioq2
         qqmax=q2max
	 smax=s_max(4)
	 endif
ccccccccccccccccccccccccccccccccccccccc
c++++ Random Q2 and W estraction and 
c     Scattered electron variables settings
          Q2=rran()*(QqMIN-QqMAX)+QqMAX
          W=rran()*(WwMIN-WwMAX)+WwMAX
	  OMEGA=(W**2+Q2-MP**2)/2./MP 
          EE=EO-OMEGA
	  PE=EE
	  ARG = (1-Q2/2./EE/EO)
	  IF(ARG.LE.-1.OR.ARG.GE.1)    GOTO 99!Regetting inconsistencies
	  IF(OMEGA.GT.EO.OR.OMEGA.LT.0)GOTO 99!Regetting inconsistencies
	  TETE = ACOS(ARG)    ! Theta of electron
	  ARG  = (Q2+2*EO*OMEGA)/(2*EO*SQRT(Q2+OMEGA**2))
	  IF(ARG.LE.-1.OR.ARG.GE.1)    GOTO 99!Regetting inconsistencies
	  FE=FIEMIN+(FIEMAX-FIEMIN)*rran()  ! Phi of electron
	  FE=FE*PIGR/180.
	  IF (FIFLAG.EQ.0) FE=2*PIGR*rran() ! Conditioned uniform extraction in 0-360
c----

	 ELSE
ccccccccccccccccccccccccccccccccccccccc
c 	estrazione rettangolo
	a=rran()
	do j=1,4
	 proba(j)=vol(j)/voltot+proba(j-1)
	enddo
	if(a.le.proba(1)) then 
	 comax=costagliot
	 comin=cosmin
	 omin=omegamin
	 omax=taglioo
	 smax=s_max(1)
	else if(a.le.proba(2)) then
	 comin=cosmin
	 comax=costagliot
	 omin=taglioo
	 omax=omegamax
	 smax=s_max(2)
	else if(a.le.proba(3)) then
	 comin=costagliot
	 comax=cosmax
	 omin=omegamin
         omax=taglioo
	 smax=s_max(3)
	else if(a.le.proba(4)) then
	 comin=costagliot
         comax=cosmax
	 omin=taglioo
         omax=omegamax
	 smax=s_max(4)
	 endif
ccccccccccccccccccccccccccccccccccccccc

c+++ Random theta_electron and omega estraction and 
c    Scattered electron variables settings
	  TETE=ACOS(rran()*(COMIN-COMAX)+COMAX)
	  OMEGA=rran()*(OMIN-OMAX)+OMAX
	  EE=EO-OMEGA
	  PE=SQRT(EE**2-ME**2)
	  Q2=4*EO*(EO-OMEGA)*(SIN(TETE/2))**2 
	  IF(I_Q2.EQ.1.and.(Q2.LT.Q2MIN.OR.Q2.GT.Q2MAX)) GOTO 99 !Reget
	  ARG=(MP**2+2*MP*OMEGA-Q2) ! W2 of the reaction
	  IF (ARG.LT.1.1) GOTO 99 !Regetting evts down 1-pion region
	  W=SQRT(ARG)
	  IF(I_W.EQ.1.and.(W.LT.WMIN.OR.W.GT.WMAX)) GOTO 99 !Regetting evts out W cut (if activated) 
	  FE=FIEMIN+(FIEMAX-FIEMIN)*rran()  ! Phi of electron
	  FE=FE*PIGR/180.
	  IF (FIFLAG.EQ.0) FE=2*PIGR*rran() ! Conditioned uniform extraction in 0-360

c----
	 ENDIF


c+++ Evaluating kinematical factors and xsection according
c+++ to the formula: d_sigma = GAMMA_V[SIGMA_T+EPS*SIGMA_L]
c+++ where SIGMA_T=(DIPOLE FIT) AND SIGMA_L=0
c+++   They will be recalculated after the cross-section calculation 
	  EPS = (1+2*(1+OMEGA**2/Q2)*(TAN(TETE/2.))**2)**(-1)
	  GAMMA_V = (1/137./2./PIGR**2)*(EE/EO)*((W**2-MP**2)/2./MP)/Q2
     @              /(1-EPS)

c+++ Starting calculation of full kinematics and
c+++ dynamics for survived events

c+++ CROSS-SECTION (e,e',gamma_virt)

ccc         ModaRad=2 
ccc         ! 0 -> nonradiative; 1-> radiative no stragling; 2 -> radiative
	 IF(I_W.EQ.1.AND.I_Q2.EQ.1) THEN
 	   d_sigma =f_allchn(w)* Gamma_V * (1+Q2/.7)**(-2)
     &             /(2*Eo*(Eo-omega)*Mp)*w
           IF(ModaRad.eq.0) THEN
           d_sigma =sigma_ep_eet(Eo*1.d+0, EE*1.d+0, TETE*1.d+0)
     &          /(2*Eo*(Eo-omega)*Mp)*w
           e_radgam= 0.
	   theta_radgam= 0.
	   phi_radgam= 0.
           ELSEIF(ModaRad.eq.1.or.ModaRad.eq.2) THEN
           mf_min = 0.1395679d+0
           if(ModaRad.eq.1) THEN
           Ttarg  = 0.0 
           Twi    = 0.0 
           Twf    = 0.0 
           else
           Ttarg  = 0.00577d+0
           Twi    = 0.00017d+0
           Twf    = 0.00083d+0
           endif
           Ztarg  = 1.d+0
           Zw     = 13.d+0
           delta  = 0.010d+0
c           mfin   = 0.93827231+0.1395679+0.1395679
           mfin   = 0.93827231+0.1395679
           CALL radiate1(sigma_ep_eet, d_sigma,d_sigma_part, mfin,
     &          Eo,    Ee,    w,    q2,    omega,    tete,    fe,
     &          Eo_vtx,Ee_vtx,w_vtx,q2_vtx,omega_vtx,tete_vtx,fe_vtx,
     &          e_radgam, theta_radgam, phi_radgam,
     &          mf_min, Ttarg,Twi,Twf,Ztarg,Zw, delta)
           d_sigma      = d_sigma      /(2*Eo*(Eo-omega)*Mp)*w
           d_sigma_part = d_sigma_part /(2*Eo*(Eo-omega)*Mp)*w
            Eo_ini=Eo
            Ee_ini=Ee
            Pe=Ee
            w_ini=w
            q2_ini=q2
            omega_ini=omega
            tete_ini=tete
            fe_ini=fe
              Eo=Eo_vtx
              Ee=Ee_vtx
              Pe=Ee
              w=w_vtx
              q2=q2_vtx
              omega=omega_vtx
              tete=tete
              fe=fe_vtx
           ELSE
            print *,' ModaRad is bad=',ModaRad 
            stop
           ENDIF
	 ELSE
          e_radgam= 0.
	  theta_radgam= 0.
	  phi_radgam= 0.
	  d_sigma =f_allchn(w)* Gamma_V * (1+Q2/.7)**(-2)
          d_sigma =sigma_ep_eet(Eo*1.d+0, EE*1.d+0, TETE*1.d+0)

	 ENDIF

ccccccccccccccccccccccccccccccccccccccccccccccc
          if(D_SIGMA.gt.Smax) then
            print *,' D_SIGMA.gt.Smax==',D_SIGMA,Smax
	    stop
          endif
          if(D_SIGMA.eq.0.) then
            goto 99
          endif
	  D_S = D_SIGMA/Smax
	  R_D_S=rran()

c++++ Regetting evts out dsigma distribution
      IF (R_D_S.GT.D_S) GOTO 99
c----
      if    (abs(Eo-Eo_ini).le.0.001.and.abs(Ee-Ee_ini).le.0.001)then
        num_soft=num_soft+1
      elseif(abs(Eo-Eo_ini).gt.0.001.and.abs(Ee-Ee_ini).le.0.001)then
        num_hardini=num_hardini+1
      elseif(abs(Eo-Eo_ini).le.0.001.and.abs(Ee-Ee_ini).gt.0.001)then
        num_hardfin=num_hardfin+1
      endif


c+++   Recalculated EPS,GAMMA_V 
	  EPS = (1+2*(1+OMEGA**2/Q2)*(TAN(TETE/2.))**2)**(-1)
	  GAMMA_V = (1/137./2./PIGR**2)*(EE/EO)*((W**2-MP**2)/2./MP)/Q2
     @              /(1-EPS)
c---
c+++ Angles of virtual photon angle (PHI)
	  ARG = (Q2+2*EO*OMEGA)/(2*EO*SQRT(Q2+OMEGA**2))
	  IF(ARG.LE.-1.OR.ARG.GE.1) GOTO 99
	  PSI=ACOS(ARG) ! Theta of virtual photon
          IF(FE.LT.PIGR) THEN
	    PHI=FE+PIGR ! Phi of virtual photon
          ELSE
	    PHI=FE-PIGR
          ENDIF
c---
c++++ Evaluating bETA in CM according to the formula:
c     BET= q/(MP+OMEGA)  where Q^2=q^2-OMEGA^2
	 BET=SQRT(Q2+OMEGA**2)/(MP+OMEGA)
	 GAMMA=1./SQRT(1.-BET**2)
         BETA(1)=0
         BETA(2)=0
         BETA(3)=BET
         BETA(4)=GAMMA
C++++ Changing BETA sign to move from CM to LAB frame 
  	 DO M=1,3
	  BETA(M)=-BETA(M)
	 ENDDO
c----
c++++ Virtual photon COS in frame=CLAS
	 UGV=SIN(PSI)*COS(PHI)
	 VGV=SIN(PSI)*SIN(PHI)
	 WGV=COS(PSI)
	 CD(1)=UGV
	 CD(2)=VGV
	 CD(3)=WGV
c++++ Scattered electron COS in: 1-> frame=CLAS;
c++++  2 -> Z along virtual photon
c++++      X in the plane containing scattered electron
c++++      Y perpendicular to the plane containing scat elec
	 UE1=SIN(TETE)*COS(FE)
	 VE1=SIN(TETE)*SIN(FE)
	 WE1=COS(TETE)
	 UE2=SIN(TETE+PSI)
	 VE2=0.
	 WE2=COS(TETE+PSI)
C++++ Using versors
C++++	Z' (Z-AX. IN REF 2) SAME VG (UGV,VGV,WGV)
C++++	Y' (Y-AX. IN REF. 2): Y'=VET(VG,VE')
	UY = VGV*WE1-WGV*VE1
	VY = WGV*UE1-UGV*WE1
	WY = UGV*VE1-VGV*UE1

	YYMOD=SQRT(UY**2+VY**2+WY**2)
	UY = UY/YYMOD
	VY = VY/YYMOD
	WY = WY/YYMOD
C++++	X' (X-AX. IN REF. 2): X'=VET(Y',Z')

	UX = VY*WGV-WY*VGV
	VX = WY*UGV-UY*WGV
	WX = UY*VGV-VY*UGV

	XXMOD=SQRT(UX**2+VX**2+WX**2)
	UX = UX/XXMOD
	VX = VX/XXMOD
	WX = WX/XXMOD
c----
c++++ Clearing DNUP vector
	 DO J=1,NVAR
	  DNUP(J)=0.
	 ENDDO  
c----
c++++ Calling final state (particles) routine where
c++++ NP= number of particles and theyr labels are in
c++++ (PARTIC(J), J=1,NP) vector and 
c++++ JCH is the selected channel
ccccccccccccccccccccccc

c+++ Choosing channel
	 CALL CANALE(W,JCH,SIG)
c---

c++++ Regetting dynamic inconsitency
        IF(JCH.EQ.0) GOTO 99
c----

ccccccccccccccccccccccc
	CALL PCLE(JCH,TARGET,NP,PARTIC)
        
ccc        print *,' JCH,NP=',JCH,NP
ccc        if(JCH.eq.15) print *,' --- 15 --- 111',NP
        

c----
C++++ Fixing energy in CM (ECM=W)
        ECM=W
c----
c++++ Assigning mass and id of final state particle
	DO J=1,NP
	 AMASS(J)=RMASS(PARTIC(J))
	 PART_ID(J)=NP_CODE(PARTIC(J))
	ENDDO
c----
c++++ Fixing new absolute variable to count number of final state parts
	NP_FINAL = NP
c----
c++++ Selecting procedure for 2 parts in final state
        IF(NP_FINAL.EQ.2) THEN
c+++++ Checking energy conservation for resonances
c+++++ production (rho, omega, phi...) in 2-part final state

          IMPQ=((W**2+AMASS(1)**2-AMASS(2)**2)**2-4*W**2*   
     $         AMASS(1)**2)/(4*W**2)

c++++++ Regetting events without energy conservation
          IF(IMPQ.LE.0) GOTO 99
c------
          IMP = SQRT(IMPQ)
c-----
c+++++ Special CM angular distribution for rho0  production (if enabled)	
	  IF(JCH.EQ.11.AND.FLAG_DS_RHO.ne.0) THEN
	  IF(FLAG_DS_RHO.EQ.1) THP=ACOS(-2.*rran()+1.) ! Uniform
          IF(FLAG_DS_RHO.EQ.2) THP=rran()*PIGR	        ! Flat distribution
          IF(FLAG_DS_RHO.EQ.3) thp=30.*degrad           ! fixed cm angle
c-----
  	  ELSE 
c+++++ Calling routine extracting angle according tabulated
c+++++ cross section
c+++++ The angle THP refers to:
c+++++ Ch  1-4  -> pion 
c+++++ Ch  5-10 -> pion (NOT DELTA!!!)
c+++++ Ch 11-14 -> rho
c+++++ Ch 19-20 -> omega (same as rho)
c+++++ Ch 25-26 -> phi   (same as rho)

           CALL ANG(JCH,W,THP)

ccc           print *,' 2================',JCH

c-----
	  ENDIF

c+++++ Inverting CM angle when a decaing particle is produced (THP=not decaing part)
c+++++ From now TH_MESON = PIGR-THP
         IF((JCH.GE.11.AND.JCH.LE.14).OR.JCH.EQ.19.OR.JCH.EQ.20) THEN
          THP=PIGR-THP
	  I_DUMMY = PART_ID(1)
c	  PART_ID(1) = PART_ID(2)
c	  PART_ID(2) = I_DUMMY
	  ENDIF
c+++++ Defining variables of first two particles in CM
c+++++ First particle has angle -THP then:
c+++++ Ch  1-4  ->   PCM(*,1)=nucleon  PCM(*,2)=pion
c+++++ and for decaing channels:
c+++++ Ch  5-10  ->  PCM(*,1)=delta    PCM(*,2)=pion
c+++++ Ch 11-14  ->  PCM(*,1)=rho      PCM(*,2)=nucleon
c+++++ Ch 19-20  ->  PCM(*,1)=omega    PCM(*,2)=nucleon
c+++++ Ch 25-26  ->  PCM(*,1)=omega    PCM(*,2)=nucleon
          PHP = 2 * PIGR * rran()
c++++++ Extracting correct phi distr for single pion channel
	if (JCH.le.4) call ang_distr_1pi(eps,thp,php)
c------	
     
          PCM(1,1)=IMP*SIN(THP)*COS(PHP)
          PCM(2,1)=IMP*SIN(THP)*SIN(PHP)
          PCM(3,1)=-IMP*COS(THP)
          PCM(1,2)=-PCM(1,1)
          PCM(2,2)=-PCM(2,1)
          PCM(3,2)=-PCM(3,1)
          PCM(4,1)=SQRT(IMPQ+AMASS(1)**2)
          PCM(4,2)=SQRT(IMPQ+AMASS(2)**2)
          PCM(5,1)=IMP
          PCM(5,2)=IMP
c++++++ Saving  CM variables (pion, delta, rho, omega)
c++++++ 1 pion channels     -> pion CM distribution
	IF(JCH.GE.1.AND.JCH.LE.4)   I_CM = 2
c++++++ delta-pi channel  -> delta CM distribution
c++++++ rho-N channel     -> rho CM distribution
c++++++ omega-N channel   -> omega CM distribution
c++++++ phi-N channel     -> phi CM distribution

	IF((JCH.GE.5.AND.JCH.LE.14).OR.(JCH.GE.19.AND.JCH.LE.20).OR.(JCH.GE.25.AND.JCH.LE.26))
     %         I_CM = 1
          DO L=1,4
           RHO_CM(L) = PCM(L,I_CM)
	   RHO_CM(5) = SQRT(PCM(4,I_CM)**2-PCM(5,I_CM)**2) ! this is the MASS not the mod(p)
          ENDDO
c------
c-----
c----    
	ELSE 

c++++ Calling GENBOD for more-than-3 parts in final state
c+ before to call genbod a test on energy is required
	   e_test=0.
	   do jj=1,np
	    E_test= E_test+ amass(jj)
	   enddo
	   if (ecm.lt.E_test) goto 99
30	   if(np.eq.3) then
	    call gen3n1bod
	   else
	    call genbod
	   endif
c----
        ENDIF 
c++++ Making Lorentz transformation for each particle
	 DO JK=1,NP
c+++++ Using vector P for each part
           DO IL=1,4
             P(IL)=PCM(IL,JK)
           ENDDO
c-----
C+++++ Lorentz boost in virtual gamma frame (result in PLAB2)
	   CALL GLOREN(BETA,P,PA)
           DO IL=1,4
             PLAB2(IL,JK)=PA(IL)
           ENDDO
           PLAB2(5,JK)=SQRT(PROD_SCAL(PA,PA))
c-----
C+++++ Rotation to come back in LAB frame (not used)
C	   CALL ROT(PA,CD)
c-----
C+++++ Rotation to come back in LAB frame using versors (results in PLAB1)
	  PA1(1)=PA(1)*UX+PA(2)*UY+PA(3)*UGV
	  PA1(2)=PA(1)*VX+PA(2)*VY+PA(3)*VGV
	  PA1(3)=PA(1)*WX+PA(2)*WY+PA(3)*WGV
	  PA(1)=PA1(1)
	  PA(2)=PA1(2)
	  PA(3)=PA1(3)
c+++++ Filling PLAB1 (final LAB vector) with variables of particles according with:
c+++++ Ch  1-4  ->   PLAB1(*,1)=nucleon PLAB1(*,2)=pion
c+++++ and for decaing channels:
c+++++ Ch  5-10  ->   PLAB1(*,1)=delta   PLAB1(*,2)=pion
c+++++ Ch 11-14  ->   PLAB1(*,1)=rho     PLAB1(*,2)=nucleon
c+++++ Ch 19-20  ->   PLAB1(*,1)=omega   PLAB1(*,2)=nucleon 
c+++++ Ch 25-26  ->   PLAB1(*,1)=phi     PLAB1(*,2)=nucleon 
c+++++ For phase space this is the final vector:
c+++++ Ch 15-18  ->   PLAB1(*,1)=nucleon PLAB1(*,2)=pion PLAB1(*,3)=pion
c+++++ Ch 21-24  ->   PLAB1(*,1)=nucleon PLAB1(*,2)=pion PLAB1(*,3)=pion PLAB1(*,4)=pion
          DO IL=1,4
             PLAB1(IL,JK)=PA(IL)
           ENDDO
           PLAB1(5,JK)=SQRT(PROD_SCAL(PA,PA))
c-----
	 ENDDO
c----

c++++ Calcul of GPhi = angle between leptonic and hadronic planes

	 Plab_Eout(1) = PE*UE1
	 Plab_Eout(2) = PE*VE1
	 Plab_Eout(3) = PE*WE1
	 
	 Plab_Pout(1) = PLAB1(1,2)
	 Plab_Pout(2) = PLAB1(2,2)
	 Plab_Pout(3) = PLAB1(3,2)
	 
	 call GPhiLH(EO_beam,Plab_Eout,Plab_Pout,res)
	 gphi = res

c++++ Staring routines ONLY for decaing particles
         FLAGD=0 ! Decay flag
         IF((JCH.GE.5.AND.JCH.LE.14).OR.(JCH.GE.19.AND.JCH.LE.20).OR.(JCH.GE.25.AND.JCH.LE.26)) THEN
c+++++ Assigning decaing particle mass to ECM
          ECM = AMASS(1)
c-----
C+++++ Increasing counter of final state parts according to differnt channels
	    NP_FINAL = NP_FINAL + 1 ! At least 3 parts in final state 
	    IF(JCH.EQ.19.or.JCH.EQ.20) NP_FINAL = NP_FINAL + 1 ! 4 parts in final state 
c-----	    
            FLAGD = 1
            Y=rran()
c+++++ Calling decaing routine to obtain mass and number of dauther partic
c+++++ (AMASS and NP are overwritten)
            CALL DECAD(JCH,NP,AMASS)      
c-----
c+++++ Fixing dauthers particles id
c++++++ Delta++
	    IF(JCH.EQ.5) THEN
	      PART_ID(1)=41
              PART_ID(2)=17
	      PART_ID(5-2)=NP_CODE(PARTIC(2)) ! moving ID of first pion in position 3             
	    END IF
c++++++ Delta+ decay: different decay according isospin coeff 
            IF(JCH.EQ.6.OR.JCH.EQ.8) THEN
	      PART_ID(5-2)=NP_CODE(PARTIC(2)) ! moving ID of first pion in position 3            
             IF(Y.LT.0.67) THEN
c+++++++ Delta+ -> p pi0
	      PART_ID(3-2)=41
	      PART_ID(4-2)=23
              FLAGC=+0.5
c-------
             ELSE
c+++++++ Delta+ -> n pi+
              PART_ID(3-2)=42
              PART_ID(4-2)=17
              FLAGC=-0.5
	     ENDIF
c-------
            ENDIF
c------
c++++++ Delta0 decay: different decay according isospin coeff 
            IF(JCH.EQ.7.OR.JCH.EQ.9)THEN
	      PART_ID(5-2)=NP_CODE(PARTIC(2)) ! moving ID of first pion in position 3             
c+++++++ Delta0 -> n pi0
             IF(Y.LT.0.67)THEN
              PART_ID(3-2)=42
              PART_ID(4-2)=23
              FLAGC=-0.5
c-------
             ELSE
c+++++++ Delta0 -> p pi-
              PART_ID(3-2)=41
              PART_ID(4-2)=-17
              FLAGC=+0.5
             ENDIF
c-------
	    ENDIF
c------
c++++++ Delta-
	    IF(JCH.EQ.10) THEN
	      PART_ID(3-2)=42
              PART_ID(4-2)=-17
	      PART_ID(5-2)=NP_CODE(PARTIC(2)) ! moving ID of first pion in position 3             
	    END IF
c++++++ Rho	    
	    IF(JCH.EQ.11) THEN
	      PART_ID(1)=41
	      PART_ID(2)=17
	      PART_ID(3)=-17
	    ENDIF
	    IF(JCH.EQ.12) THEN
	      PART_ID(1)=42
	      PART_ID(2)=17
	      PART_ID(3)=23
	    ENDIF
	    IF(JCH.EQ.13) THEN
	      PART_ID(1)=41	    
	      PART_ID(2)=23
	      PART_ID(3)=-17
	    ENDIF
	    IF(JCH.EQ.14) THEN
	      PART_ID(1)=42	    
	      PART_ID(2)=17
	      PART_ID(3)=-17
	    ENDIF
c++++++ Omega	    
	    IF(JCH.EQ.19) THEN
	      PART_ID(1)=41	    
	      PART_ID(2)=17
	      PART_ID(3)=-17
	      PART_ID(4)=23
	    ENDIF
	    IF(JCH.EQ.20) THEN
	      PART_ID(1)=42	    
	      PART_ID(2)=17
	      PART_ID(3)=-17
	      PART_ID(4)=23
	    ENDIF
c++++++ phi NOTE: no decad routine is called    
	    IF(JCH.EQ.25) THEN
             IF(Y.LT.0.491) then ! phi->k+ k-
	      PART_ID(1)=41	    
	      PART_ID(2)=18
	      PART_ID(3)=-18
	      np = 2
	      np_final = 3
	      amass(1) = 0.49367
	      amass(2) = 0.49367
	     elseif(Y.GT.0.491.and.y.lt.0.834) then ! phi->ks kl
	      PART_ID(1)=41	    
	      PART_ID(2)=37
	      PART_ID(3)=38
	      np = 2
	      np_final = 3
	      amass(1) = 0.49772
	      amass(2) = 0.49772
	     elseif(Y.GT.0.834) then ! phi->pi+pi-pi0 !!! IT'S WRONG but simple
	      PART_ID(1)=41	    
	      PART_ID(2)=17
	      PART_ID(3)=-17
	      PART_ID(4)=23
	      np_final = 4
	      np = 3
	      amass(1) =  0.1395679
	      amass(2) =  0.1395679
 	      amass(3) =  0.1349743
            endif
	    ENDIF
	    IF(JCH.EQ.26) THEN ! to be implemented
	      PART_ID(1)=42	    
	      PART_ID(2)=18
	      PART_ID(3)=-18
	    ENDIF

c------
c-----
c+++++ Rho0 decay using angular distribution function
	   IF(JCH.EQ.11) THEN
	    CALL RHO_DECAY(ECM,FLAG_Decay_RHO,Q2,W,GPhi,EPS)
C++++++ Saving REST FRAME decaing vars
               DO IL=1,5
                 P_REST(IL) = PCM(IL,1)
               ENDDO
c++++++ Rotating XY from hadronic to leptonic plane (rho prod angles = -thp php)
c+ Hadronic plane
c       Px_spool =   COS(-php)*PCM(1,1) + SIN(-php)*PCM(2,1) !pion 1
c       Py_spool =  -SIN(-php)*PCM(1,1) + COS(-php)*PCM(2,1)
c       PCM(1,1) = Px_spool 
c       PCM(2,1) = Py_spool 
c       Px_spool =   COS(-php)*PCM(1,2) + SIN(-php)*PCM(2,2) !pion 2
c       Py_spool =  -SIN(-php)*PCM(1,2) + COS(-php)*PCM(2,2)
c       PCM(1,2) = Px_spool 
c       PCM(2,2) = Py_spool 
c       rho_CM_r(1) =   COS(-php)*rho_CM(1) + SIN(-php)*rho_CM(2) !rho
c       rho_CM_r(2) =  -SIN(-php)*rho_CM(1) + COS(-php)*rho_CM(2)
c       rho_CM_r(3) =   rho_CM(3)

c+ hadronic cms z along rho

c+ tmp change of thp (I remember that thp is the proton angle in cm for decaying particles)
c+ boost along rho momentum
	BETA_C = sqrt(rho_cm(4)**2-rho_cm(5)**2)/rho_cm(4)
	GAMMA_C = 1/SQRT(1 - BETA_C**2)

	ptr(4) =  GAMMA_C*(pcm(4,1) + BETA_C*pcm(3,1))   !pion 1
        ptr(1) 	= pcm(1,1)
        ptr(2)	= pcm(2,1)
        ptr(3)  = GAMMA_C*(pcm(3,1) + BETA_C*pcm(4,1))
	ptr(5) = sqrt(ptr(1)**2+ ptr(2)**2+ ptr(3)**2)
	do  j= 1,5
         pcm(j,1)=ptr(j)
	enddo
	ptr(4) =  GAMMA_C*(pcm(4,2) + BETA_C*pcm(3,2))   !pion 2
        ptr(1) 	= pcm(1,2)
        ptr(2)	= pcm(2,2)
        ptr(3)  = GAMMA_C*(pcm(3,2) + BETA_C*pcm(4,2))
	ptr(5) = sqrt(ptr(1)**2+ ptr(2)**2+ ptr(3)**2)
	do  j= 1,5
         pcm(j,2)=ptr(j)
	enddo
c+	Rotation to the hadronic CMS 
c+	(z axis along q vector but x axis still on the hadronic plane)

       thp = pigr-thp
       Px_spool =   COS(-thp)*PCM(1,1) - SIN(-thp)*PCM(3,1) !pion 1
       Py_spool =   SIN(-thp)*PCM(1,1) + COS(-thp)*PCM(3,1)
       PCM(1,1) = Px_spool 
       PCM(3,1) = Py_spool 
       Px_spool =   COS(-thp)*PCM(1,2) - SIN(-thp)*PCM(3,2) !pion 2
       Py_spool =   SIN(-thp)*PCM(1,2) + COS(-thp)*PCM(3,2)
       PCM(1,2) = Px_spool 
       PCM(3,2) = Py_spool 
c       rho_cm_r(1) =   COS(-thp)*rho_CM(1) - SIN(-thp)*rho_CM(3) !rho
c       rho_CM_r(2) =   rho_CM(2)
c       rho_cm_r(3) =   SIN(-thp)*rho_CM(1) + COS(-thp)*rho_CM(3)
       thp=thp-pigr
c+ Hadronic plane
       Px_spool =   COS(-php)*PCM(1,1) + SIN(-php)*PCM(2,1) !pion 1
       Py_spool =  -SIN(-php)*PCM(1,1) + COS(-php)*PCM(2,1)
       PCM(1,1) = Px_spool 
       PCM(2,1) = Py_spool 
       Px_spool =   COS(-php)*PCM(1,2) + SIN(-php)*PCM(2,2) !pion 2
       Py_spool =  -SIN(-php)*PCM(1,2) + COS(-php)*PCM(2,2)
       PCM(1,2) = Px_spool 
       PCM(2,2) = Py_spool 
c       Px_spool =   COS(-php)*rho_CM_r(1) + SIN(-php)*rho_CM_r(2) !rho
c       Py_spool =  -SIN(-php)*rho_CM_r(1) + COS(-php)*rho_CM_r(2)
c       rho_cm_r(1) = Px_spool
c       rho_CM_r(2) = Py_spool

c++++++ Looping over two pions
            DO JK=1,2
c               DO IL=1,4
c                 P(IL) = PCM(IL,JK)
c               ENDDO
c             BETD(1)=-RHO_CM_r(1)/RHO_CM(4)
c             BETD(2)=-RHO_CM_r(2)/RHO_CM(4)
c             BETD(3)=-RHO_CM_r(3)/RHO_CM(4)
c             BETM=SQRT(BETD(1)**2+BETD(2)**2+BETD(3)**2)
c             BETD(4)=1./(SQRT(1.-BETM**2))
c              CALL GLOREN(BETD,P,PA)
              DO IL=1,4
                P(IL)=PCM(IL,JK)
              ENDDO
	     CALL GLOREN(BETA,P,PA)
              DO IL=1,4
                PLAB2(IL,JK)=PA(IL)
              ENDDO
           PLAB2(5,JK)=SQRT(PROD_SCAL(PA,PA))

	  PA1(1)=PA(1)*UX+PA(2)*UY+PA(3)*UGV
	  PA1(2)=PA(1)*VX+PA(2)*VY+PA(3)*VGV
	  PA1(3)=PA(1)*WX+PA(2)*WY+PA(3)*WGV
	  PA(1)=PA1(1)
	  PA(2)=PA1(2)
	  PA(3)=PA1(3)
c+++++++ Filling lab vector for two dauther pions and defining final configuration of PLAB1:
c+++++++ Ch 11  ->   PLAB1(*,3)=pi+     PLAB1(*,4)=pi-
            DO IL=1,4
             PLAB1(IL,JK+2)=PA(IL)
           ENDDO
           PLAB1(5,JK+2)=SQRT(PROD_SCAL(PA,PA))
            ENDDO
c-------
c------ Closing rho decay sub


c------
c-----
c+++++ Omega decay using angular distribution function
	  ELSE IF(JCH.EQ.19) THEN
	    CALL OMEGA_DECAY(ECM,FLAG_Decay_RHO,Q2,W,GPhi,EPS)
C++++++ Saving REST FRAME decaing vars
               DO IL=1,5
                 P_REST(IL) = PCM(IL,1)
               ENDDO
c++++++ Rotating XY from hadronic to leptonic plane (omega prod angles = -thp php)
c+ Hadronic plane
c       Px_spool =   COS(-php)*PCM(1,1) + SIN(-php)*PCM(2,1) !pion 1
c       Py_spool =  -SIN(-php)*PCM(1,1) + COS(-php)*PCM(2,1)
c       PCM(1,1) = Px_spool 
c       PCM(2,1) = Py_spool 
c       Px_spool =   COS(-php)*PCM(1,2) + SIN(-php)*PCM(2,2) !pion 2
c       Py_spool =  -SIN(-php)*PCM(1,2) + COS(-php)*PCM(2,2)
c       PCM(1,2) = Px_spool 
c       PCM(2,2) = Py_spool 
c       Px_spool =   COS(-php)*PCM(1,3) + SIN(-php)*PCM(2,3) !pion 3
c       Py_spool =  -SIN(-php)*PCM(1,3) + COS(-php)*PCM(2,3)
c       PCM(1,3) = Px_spool 
c       PCM(2,3) = Py_spool 
c       rho_CM_r(1) =   COS(-php)*rho_CM(1) + SIN(-php)*rho_CM(2) !rho
c       rho_CM_r(2) =  -SIN(-php)*rho_CM(1) + COS(-php)*rho_CM(2)
c       rho_CM_r(3) =   rho_CM(3)

c+ hadronic cms z along rho

c+ tmp change of thp (I remember that thp is the proton angle in cm for decaying particles)
c+ boost along rho momentum
	BETA_C = sqrt(rho_cm(4)**2-rho_cm(5)**2)/rho_cm(4)
	GAMMA_C = 1/SQRT(1 - BETA_C**2)

	ptr(4) =  GAMMA_C*(pcm(4,1) + BETA_C*pcm(3,1))   !pion 1
        ptr(1) 	= pcm(1,1)
        ptr(2)	= pcm(2,1)
        ptr(3)  = GAMMA_C*(pcm(3,1) + BETA_C*pcm(4,1))
	ptr(5) = sqrt(ptr(1)**2+ ptr(2)**2+ ptr(3)**2)
	do  j= 1,5
         pcm(j,1)=ptr(j)
	enddo
	ptr(4) =  GAMMA_C*(pcm(4,2) + BETA_C*pcm(3,2))   !pion 2
        ptr(1) 	= pcm(1,2)
        ptr(2)	= pcm(2,2)
        ptr(3)  = GAMMA_C*(pcm(3,2) + BETA_C*pcm(4,2))
	ptr(5) = sqrt(ptr(1)**2+ ptr(2)**2+ ptr(3)**2)
	do  j= 1,5
         pcm(j,2)=ptr(j)
	enddo
	ptr(4) =  GAMMA_C*(pcm(4,3) + BETA_C*pcm(3,3))   !pion 3
        ptr(1) 	= pcm(1,3)
        ptr(2)	= pcm(2,3)
        ptr(3)  = GAMMA_C*(pcm(3,3) + BETA_C*pcm(4,3))
	ptr(5) = sqrt(ptr(1)**2+ ptr(2)**2+ ptr(3)**2)
	do  j= 1,5
         pcm(j,3)=ptr(j)
	enddo
c+	Rotation to the hadronic CMS 
c+	(z axis along q vector but x axis still on the hadronic plane)

       thp = pigr-thp
       Px_spool =   COS(-thp)*PCM(1,1) - SIN(-thp)*PCM(3,1) !pion 1
       Py_spool =   SIN(-thp)*PCM(1,1) + COS(-thp)*PCM(3,1)
       PCM(1,1) = Px_spool 
       PCM(3,1) = Py_spool 
       Px_spool =   COS(-thp)*PCM(1,2) - SIN(-thp)*PCM(3,2) !pion 2
       Py_spool =   SIN(-thp)*PCM(1,2) + COS(-thp)*PCM(3,2)
       PCM(1,2) = Px_spool 
       PCM(3,2) = Py_spool
       Px_spool =   COS(-thp)*PCM(1,3) - SIN(-thp)*PCM(3,3) !pion 3
       Py_spool =   SIN(-thp)*PCM(1,3) + COS(-thp)*PCM(3,3)
       PCM(1,3) = Px_spool 
       PCM(3,3) = Py_spool
 
c       rho_cm_r(1) =   COS(-thp)*rho_CM(1) - SIN(-thp)*rho_CM(3) !rho
c       rho_CM_r(2) =   rho_CM(2)
c       rho_cm_r(3) =   SIN(-thp)*rho_CM(1) + COS(-thp)*rho_CM(3)
       thp=thp-pigr
c+ Hadronic plane
       Px_spool =   COS(-php)*PCM(1,1) + SIN(-php)*PCM(2,1) !pion 1
       Py_spool =  -SIN(-php)*PCM(1,1) + COS(-php)*PCM(2,1)
       PCM(1,1) = Px_spool 
       PCM(2,1) = Py_spool 
       Px_spool =   COS(-php)*PCM(1,2) + SIN(-php)*PCM(2,2) !pion 2
       Py_spool =  -SIN(-php)*PCM(1,2) + COS(-php)*PCM(2,2)
       PCM(1,2) = Px_spool 
       PCM(2,2) = Py_spool 
       Px_spool =   COS(-php)*PCM(1,3) + SIN(-php)*PCM(2,3) !pion 3
       Py_spool =  -SIN(-php)*PCM(1,3) + COS(-php)*PCM(2,3)
       PCM(1,3) = Px_spool 
       PCM(2,3) = Py_spool 
c       Px_spool =   COS(-php)*rho_CM_r(1) + SIN(-php)*rho_CM_r(2) !rho
c       Py_spool =  -SIN(-php)*rho_CM_r(1) + COS(-php)*rho_CM_r(2)
c       rho_cm_r(1) = Px_spool
c       rho_CM_r(2) = Py_spool

c++++++ Looping over THREE(!) pions
            DO JK=1,3
c               DO IL=1,4
c                 P(IL) = PCM(IL,JK)
c               ENDDO
c             BETD(1)=-RHO_CM_r(1)/RHO_CM(4)
c             BETD(2)=-RHO_CM_r(2)/RHO_CM(4)
c             BETD(3)=-RHO_CM_r(3)/RHO_CM(4)
c             BETM=SQRT(BETD(1)**2+BETD(2)**2+BETD(3)**2)
c             BETD(4)=1./(SQRT(1.-BETM**2))
c              CALL GLOREN(BETD,P,PA)
              DO IL=1,4
                P(IL)=PCM(IL,JK)
              ENDDO
	     CALL GLOREN(BETA,P,PA)
              DO IL=1,4
                PLAB2(IL,JK)=PA(IL)
              ENDDO
           PLAB2(5,JK)=SQRT(PROD_SCAL(PA,PA))

	  PA1(1)=PA(1)*UX+PA(2)*UY+PA(3)*UGV
	  PA1(2)=PA(1)*VX+PA(2)*VY+PA(3)*VGV
	  PA1(3)=PA(1)*WX+PA(2)*WY+PA(3)*WGV
	  PA(1)=PA1(1)
	  PA(2)=PA1(2)
	  PA(3)=PA1(3)
c+++++++ Filling lab vector for two dauther pions and defining final configuration of PLAB1:
c+++++++ Ch 19  ->   PLAB1(*,3)=pi+     PLAB1(*,4)=pi-     PLAB1(*,5)=pi0
            DO IL=1,4
             PLAB1(IL,JK+2)=PA(IL)
           ENDDO
           PLAB1(5,JK+2)=SQRT(PROD_SCAL(PA,PA))
            ENDDO
c-------
c------ Closing omega decay sub


	   ELSE
c+++++ Calling GENBOD for all others decay 
	   if(np.eq.3) then 
	    call gen3n1bod
	   else
	    call genbod
	   endif
c++++++ Making Lorentz transf
            DO JK=3,2+NP
               DO IL=1,4
                 P(IL) = PCM(IL,JK-2)
               ENDDO
            BETD(1)=-PLAB1(1,1)/PLAB1(4,1)
            BETD(2)=-PLAB1(2,1)/PLAB1(4,1)
            BETD(3)=-PLAB1(3,1)/PLAB1(4,1)
            BETM=SQRT(BETD(1)**2+BETD(2)**2+BETD(3)**2)
            BETD(4)=1./(SQRT(1.-BETM**2))
c++++++ To have decay distribution in CM is sufficient to comment 
c++++++ next trasformation 

c++++++ Lorents boost to move from CM of decaing part to LAB frame
              CALL GLOREN(BETD,P,PA)
              
c+++++++ Filling lab vector for two dauther pions and defining final configuration of PLAB1:
c+++++++ Ch  5-9   ->   PLAB1(*,3)=nucleon  PLAB1(*,4)=pion
c+++++++ Ch 12-14  ->   PLAB1(*,3)=pion     PLAB1(*,4)=pion
c+++++++ Ch 19-20  ->   PLAB1(*,3)=pion     PLAB1(*,4)=pion  PLAB1(*,5)=pion
c+++++++ Ch 25-26  ->   PLAB1(*,3)=k+     PLAB1(*,4)=k-  
c+++++++++++++ or  ->   PLAB1(*,3)=ks     PLAB1(*,4)=kl  
c+++++++++++++ or  ->   PLAB1(*,3)=pion     PLAB1(*,4)=pion  PLAB1(*,5)=pion

              DO IL=1,4
                PLAB1(IL,JK)=PA(IL)
              ENDDO
            PLAB1(5,JK)= SQRT(PROD_SCAL(PA,PA))
            ENDDO
c------
	   ENDIF
c----- Closing GENBOD loop
          ENDIF         
c---- Closing routines for decayng particles



c+++++ RESTORE initial and final electrons variavles before
c+++++ writing NTUPLEs and BOS
c+++++ angles of the electrons were not changed
          IF( ModaRad.eq.1.or.ModaRad.eq.2) THEN
            Eo=Eo_ini
            Ee=Ee_ini
            Pe=Ee
            w=w_ini
            q2=q2_ini
            omega=omega_ini
            tete=tete_ini
            fe=fe_ini
          ENDIF
c-----

ccc          if(JCH.eq.15) print *,' --- 15 --- 444',NP
ccc          print *,' NTUPLES'

C++++ Booking ntuple  General infos
	DNUP(NVAR-2)=JCH
	DNUP(NVAR-3)=OMEGA
	DNUP(NVAR-4)=Q2
c-----
c+++++ DNUP vector contains all particles in final state according the following scheme:
c+++++ DNUP(j*1)=p_x DNUP(j*2)=p_y DNUP(j*3)=p_z  DNUP(j*4)=energy DNUP(j*5)=momentum
c+++++ DNUP(j*6)=theta DNUP(j*7)=phi  beeing j the particles of final state 
c+++++ For all channels DNUP(1-7) are devoted to electron
c+++++ Channel number is conteined in JCH variable
c+++++ Number of final state particles is in NP_FINAL
c++++ Ch1 
c++      neutron  PLAB1(*,1) = DNUP(2*(1->7))    part_id(1) 	NP_FINAL = 2
c++      pi+      PLAB1(*,2) = DNUP(3*(1->7))    part_id(2)
c++++ Ch2
c++      proton   PLAB1(*,1) = DNUP(2*(1->7))    part_id(1) 	NP_FINAL = 2 
c++      pi0      PLAB1(*,2) = DNUP(3*(1->7))    part_id(2)
c++++ Ch 5
c++      delta++    PLAB1(*,1)   
c++      pi- (1st)  PLAB1(*,2)   same as PLAB1(*,5)
c++      proton     PLAB1(*,3) = DNUP(2*(1->7))    part_id(1)  	NP_FINAL = 3
c++      pi+        PLAB1(*,4) = DNUP(3*(1->7))    part_id(2)
c++      pi- (1st)  PLAB1(*,5) = DNUP(4*(1->7))    part_id(3)
c++++ Ch 6
c++      delta+     PLAB1(*,1)   
c++      pi0 (1st)  PLAB1(*,2)   same as PLAB1(*,5)
c++      nucleon    PLAB1(*,3) = DNUP(2*(1->7))    part_id(1)   NP_FINAL = 3
c++      pion       PLAB1(*,4) = DNUP(3*(1->7))    part_id(2)
c++      pi0 (1st)  PLAB1(*,5) = DNUP(4*(1->7))    part_id(3)
c++++ Ch 7
c++      delta0     PLAB1(*,1)   
c++      pi+ (1st)  PLAB1(*,2)   same as PLAB1(*,5)
c++      nucleon    PLAB1(*,3) = DNUP(2*(1->7))    part_id(1)   NP_FINAL = 3 
c++      pion       PLAB1(*,4) = DNUP(3*(1->7))    part_id(2)
c++      pi+ (1st)  PLAB1(*,5) = DNUP(4*(1->7))    part_id(3)
c++++ Ch 11
c++      rho0       PLAB1(*,1)   
c++      proton     PLAB1(*,2) = DNUP(2*(1->7))    part_id(1)    NP_FINAL = 3
c++      pi+        PLAB1(*,3) = DNUP(3*(1->7))    part_id(2)
c++      pi-        PLAB1(*,4) = DNUP(4*(1->7))    part_id(3)
c++++ Ch 12
c++      rho+       PLAB1(*,1)   
c++      neutron    PLAB1(*,2) = DNUP(2*(1->7))    part_id(1)    NP_FINAL = 3
c++      pi+        PLAB1(*,3) = DNUP(3*(1->7))    part_id(2)
c++      pi0        PLAB1(*,4) = DNUP(4*(1->7))    part_id(3)
c++++ Ch 15
c++      proton    PLAB1(*,1) = DNUP(2*(1->7))    part_id(1)    NP_FINAL = 3
c++      pi+       PLAB1(*,2) = DNUP(3*(1->7))    part_id(2)
c++      pi-       PLAB1(*,3) = DNUP(4*(1->7))    part_id(3)
c++++ Ch 16
c++      proton    PLAB1(*,1) = DNUP(2*(1->7))    part_id(1)    NP_FINAL = 3
c++      pi0       PLAB1(*,2) = DNUP(3*(1->7))    part_id(2)
c++      pi+       PLAB1(*,3) = DNUP(4*(1->7))    part_id(3)
c++++ Ch 19
c++      omega0     PLAB1(*,1)   
c++      proton     PLAB1(*,2) = DNUP(2*(1->7))    part_id(1)   NP_FINAL = 4
c++      pi+        PLAB1(*,3) = DNUP(3*(1->7))    part_id(2)
c++      pi-        PLAB1(*,4) = DNUP(4*(1->7))    part_id(3)
c++      pi0        PLAB1(*,5) = DNUP(4*(1->7))    part_id(4)
c++++ Ch 21	 
c++      proton    PLAB1(*,1) = DNUP(2*(1->7))    part_id(1)    NP_FINAL = 4
c++      pi+       PLAB1(*,2) = DNUP(3*(1->7))    part_id(2)
c++      pi-       PLAB1(*,3) = DNUP(4*(1->7))    part_id(3)
c++      pi0       PLAB1(*,4) = DNUP(4*(1->7))    part_id(4)
c++++ Ch 22	 
c++      neutron   PLAB1(*,1) = DNUP(2*(1->7))    part_id(1)    NP_FINAL = 4
c++      pi+       PLAB1(*,2) = DNUP(3*(1->7))    part_id(2)
c++      pi+       PLAB1(*,3) = DNUP(4*(1->7))    part_id(3)
c++      pi-       PLAB1(*,4) = DNUP(4*(1->7))    part_id(4)
c++++ Ch 25
c++      phi0       PLAB1(*,1)   
c++      proton     PLAB1(*,2) = DNUP(2*(1->7))    part_id(1) 
c+++       k+       PLAB1(*,3) = DNUP(3*(1->7))    part_id(2)  NP_FINAL = 3
c+++       k-       PLAB1(*,4) = DNUP(4*(1->7))    part_id(3)
c+++ or
c+++       kl       PLAB1(*,3) = DNUP(3*(1->7))    part_id(2)  NP_FINAL = 3
c+++       ks       PLAB1(*,4) = DNUP(4*(1->7))    part_id(3)
c+++ or
c+++       pi+       PLAB1(*,3) = DNUP(3*(1->7))    part_id(2) NP_FINAL = 4
c+++       pi-       PLAB1(*,4) = DNUP(4*(1->7))    part_id(3)
c+++       pi0       PLAB1(*,5) = DNUP(4*(1->7))    part_id(4)
c++++ Ch 27	 
c++      proton    PLAB1(*,1) = DNUP(2*(1->7))    part_id(1)    NP_FINAL = 5
c++      pi+       PLAB1(*,2) = DNUP(3*(1->7))    part_id(2)
c++      pi+       PLAB1(*,3) = DNUP(4*(1->7))    part_id(3)
c++      pi-       PLAB1(*,4) = DNUP(4*(1->7))    part_id(4)
c++      pi-       PLAB1(*,5) = DNUP(3*(1->7))    part_id(5)
c++++ Ch 28	 
c++      proton    PLAB1(*,1) = DNUP(2*(1->7))    part_id(1)    NP_FINAL = 6
c++      pi+       PLAB1(*,2) = DNUP(3*(1->7))    part_id(2)
c++      pi+       PLAB1(*,3) = DNUP(4*(1->7))    part_id(3)
c++      pi-       PLAB1(*,4) = DNUP(4*(1->7))    part_id(4)
c++      pi-       PLAB1(*,5) = DNUP(3*(1->7))    part_id(5)
c++      pi0       PLAB1(*,6) = DNUP(3*(1->7))    part_id(6)
c++++ Ch 29	 
c++      neutron   PLAB1(*,1) = DNUP(2*(1->7))    part_id(1)    NP_FINAL = 6
c++      pi+       PLAB1(*,2) = DNUP(3*(1->7))    part_id(2)
c++      pi+       PLAB1(*,3) = DNUP(4*(1->7))    part_id(3)
c++      pi+       PLAB1(*,4) = DNUP(4*(1->7))    part_id(4)
c++      pi-       PLAB1(*,5) = DNUP(3*(1->7))    part_id(5)
c++      pi-       PLAB1(*,6) = DNUP(3*(1->7))    part_id(6)
c-----
C+++++ Electron variables in first 7 values of DNUP 
	 DNUP(1)= PE*UE1
	 DNUP(2)= PE*VE1
	 DNUP(3)= PE*WE1
	 DNUP(4)= EE
	 DNUP(5)= PE
        DNUP(6)=ACOS(WE1)*RADDEG
        FI=ATAN(VE1/UE1)*RADDEG
        IF(UE1.LT.0.AND.VE1.GT.0) FI=FI+180
        IF(UE1.LT.0.AND.VE1.LT.0) FI=FI+180
        IF(UE1.GT.0.AND.VE1.LT.0) FI=FI+360
        DNUP(7)=FI
c-----
C+++++ Loop over NP particles (variables in LAB frame)
c++++++ No decaing particles (NP is the original one = nucleon+number of pions in final state)
        IF(FLAGD.EQ.0) THEN
	DO J=1,NP
	 K=J*7
c++++++	Filling DNUP 
	 DO L=1,5
	  DNUP(K+L)=PLAB1(L,J)
	 ENDDO
         IT=13+7*(J-1)
         IG=IT+1
         IF(DNUP(IT-1).NE.0.) THEN
           TETA=ACOS(DNUP(IT-3)/DNUP(IT-1))*RADDEG
           DNUP(IT)=TETA
           U=DNUP(IG-6)/DNUP(IG-2)
           V=DNUP(IG-5)/DNUP(IG-2)
           FI=ATAN(V/U)*RADDEG
           IF(U.LT.0.AND.V.GT.0) FI=FI+180
           IF(U.LT.0.AND.V.LT.0) FI=FI+180
           IF(U.GT.0.AND.V.LT.0) FI=FI+360
           DNUP(IG)=FI
         ENDIF    
	ENDDO
c------
       ELSE
c++++++ Decaing particles (NP refers to the number of decaing particle final state)
         DO J=1,1+NP
          K=J*7
          DO L=1,5
c+++++++ Mooving first pion (produced with delta) in last position PLAB1(*,5)=PLAB1(*,2)=pion     
          IF(JCH.GE.5.AND.JCH.LE.10) THEN
            IF(J.EQ.1) PLAB1(L,3+NP)=PLAB1(L,2)
            DNUP(K+L)=PLAB1(L,J+2)
          ELSE
             DNUP(K+L)=PLAB1(L,J+1)
C	   IF(J.EQ.2.OR.J.EQ.3) DNUP(K+L)=PCM(L,J-1) ! THE W_RHO DISTR IN CM
          ENDIF 
            ENDDO
c+++++++ Calculating theta and phi of each particle and filling DNUP(6*j) DNUP(7*j)           
	 IT=13+7*(J-1)
	 IG=IT+1
	 IF(DNUP(IT-1).NE.0.) THEN
	   ARG=DNUP(IT-3)/DNUP(IT-1)
	   IF(ARG.GT.1.) ARG=1.
           IF(ARG.LT.-1.) ARG=-1.
	   TETA=ACOS(ARG)*RADDEG
	   DNUP(IT)=TETA
	   U=DNUP(IG-6)/DNUP(IG-2)
	   V=DNUP(IG-5)/DNUP(IG-2)
	   FI=ATAN(V/U)*RADDEG
	   IF(U.LT.0.AND.V.GT.0) FI=FI+180
	   IF(U.LT.0.AND.V.LT.0) FI=FI+180
	   IF(U.GT.0.AND.V.LT.0) FI=FI+360
	   DNUP(IG)=FI
	 ENDIF
	ENDDO
c------
       ENDIF

c++++++ Changing NP according maximum number of particles in FS 
c++++++ (at this point NP=NP_FINAL only if 1 ch is enabled!!!!)
        NP=(NVAR-12)/7
C-------------------------------------------------------------------
C     DISTRIBUZIONE DEI PIONI RILEVATI AL VARIARE DI W DISTINGUENDO
C     IL CASO DI UN PIONE, DUE PIONI E TRE PIONI
C------------------------------------------------------------------- 
         IF(JCH.LE.4) THEN
            IF(DNUP(6).NE.185.AND.DNUP(20).NE.185) DNUP(NVAR-1)=W
         ELSEIF(JCH.GE.5.AND.JCH.LE.10) THEN
           IF((DNUP(20).NE.185.AND.DNUP(27).NE.185).AND.
     +       DNUP(6).NE.185) DNUP(NVAR)=W
           IF(((DNUP(20).EQ.185.AND.DNUP(27).NE.185).OR.(DNUP(20)
     +     .NE.185.AND.DNUP(27).EQ.185)).AND.
     +     DNUP(6).NE.185) DNUP(NVAR-1)=W
         ELSEIF(JCH.GE.11.AND.JCH.LE.18) THEN
           IF((DNUP(20).NE.185.AND.DNUP(27).NE.185).AND.DNUP(6).NE.185)
     +     DNUP(NVAR)=W
           IF(((DNUP(20).EQ.185.AND.DNUP(27).NE.185).OR.(DNUP(20)
     +     .NE.185.AND.DNUP(27).EQ.185)).AND.
     +      DNUP(6).NE.185) DNUP(NVAR-1)=W
         ELSEIF(JCH.GE.19) THEN
           IF(((DNUP(20).NE.185.AND.DNUP(27).NE.185).OR.(DNUP(27).NE.185
     +     .AND.DNUP(34).NE.185).OR.(DNUP(20).NE.185.AND.
     +     DNUP(34).NE.185).OR.(DNUP(20).NE.185.AND.DNUP(27)
     +     .NE.185.AND.DNUP(34).NE.185)).AND
     +     .DNUP(6).NE.185) DNUP(NVAR)=W
           IF(((DNUP(20).NE.185.AND.DNUP(27).EQ.185.AND.DNUP(34).EQ.185)
     +     .OR.(DNUP(20).EQ.185.AND.DNUP(27).NE.185.AND.DNUP(34).EQ.
     +     185).OR.(DNUP(20).EQ.185.AND.DNUP(27).EQ.185.AND.DNUP(34).
     +     NE.185)).AND.DNUP(6).NE.185) DNUP(NVAR-1)=W
         ELSE   
         CONTINUE
         ENDIF  
          
C+ Filling second ntuple ENUP according the following:
c+ ENUP(1-4)   -> electron
c+ ENUP(5-8)   -> proton
c+ ENUP(9-12)  -> neutron
c+ ENUP(13-16) -> 1pi+
c+ ENUP(17-20) -> 1pi-
c+ ENUP(17-24) -> 1pi0
c+ ENUP(25-28) -> 2pi+
c+ ENUP(29-32) -> 2pi-
c+ ENUP(33-36) -> 2pi0
c+ ENUP(37-40) -> 3pi+
c+ ENUP(41-44) -> 3pi-
c+ ENUP(45-48) -> 3pi0
C+ ENUP(49-52) -> CM vars of:
c+                            pion  ch 1-4
c+                            delta ch 5-10
c+                            rho   ch 11-14
c+                            omega ch 19-20
c+                            omega ch 25-26
c+ ENUP(53-56) -> REST FRAME vars of decaing:
c+                            pion in  rho   ch 11-14		    
C+ ENUP(57-60) -> 1k
c+ ENUP(61-64) -> 2k

        E_TOTAL=0.
        X_TOTAL=0.
        Y_TOTAL=0.
        Z_TOTAL=0.
	DO J = 1,100
	 ENUP(J) = -1000.
	ENDDO
	SECTOR = 0.
C+ Electron variables
         E_TOTAL=E_TOTAL+DNUP(4)
         X_TOTAL=X_TOTAL+DNUP(1)
         Y_TOTAL=Y_TOTAL+DNUP(2)
         Z_TOTAL=Z_TOTAL+DNUP(3)
	 ENUP(1) = DNUP(4)
	 ENUP(2) = DNUP(5)
	 ENUP(3) = DNUP(6)
	 ENUP(4) = DNUP(7)
c+ Loop over particles 
	 J_PP = 12
	 J_PM = 16
	 J_PZ = 20
	 
	 DO J = 1,NP_FINAL
	 IF(PART_ID(J).EQ.41) J_P = 4
	 IF(PART_ID(J).EQ.42) J_P = 8
	 IF(PART_ID(J).EQ.17) THEN
	   J_P  = J_PP
	   J_PP = J_PP + 12
	 ELSE IF(PART_ID(J).EQ.-17) THEN
	   J_P  = J_PM
	   J_PM = J_PM + 12
	 ELSE IF(PART_ID(J).EQ.23) THEN
	   J_P  = J_PZ
	   J_PZ = J_PZ + 12
	 ELSE IF(PART_ID(J).EQ.18.or.PART_ID(J).EQ.37) then
           J_P  = 56
	 ELSE IF(PART_ID(J).EQ.-18.or.PART_ID(J).EQ.38) then
           J_P  = 60
	 ENDIF
          if(DNUP(J*7+4).gt.0.) then
          E_TOTAL=E_TOTAL+DNUP(J*7+4)
          X_TOTAL=X_TOTAL+DNUP(J*7+1)
          Y_TOTAL=Y_TOTAL+DNUP(J*7+2)
          Z_TOTAL=Z_TOTAL+DNUP(J*7+3)
          endif
       	  ENUP(J_P+1) = DNUP(J*7+4)
       	  ENUP(J_P+2) = DNUP(J*7+5)
       	  ENUP(J_P+3) = DNUP(J*7+6)
       	  ENUP(J_P+4) = DNUP(J*7+7)


c+ Checking at least one part with PHI in sector 2
	IF(ENUP(J_P+4).GT.30.AND.ENUP(J_P+4).LT.90) SECTOR = 2.  	  

       	 ENDDO
	      	 
	IF(JCH.LE.14.OR.(JCH.GE.19.AND.JCH.LE.20).or.(JCH.GE.25.AND.JCH.LE.26)) THEN
	 ENUP(49) = RHO_CM(4)
	 ENUP(50) = SQRT(RHO_CM(4)**2-RHO_CM(5)**2)
          IF(ENUP(50).NE.0.) THEN
	   ENUP(51) = ACOS(RHO_CM(3)/ENUP(50))*RADDEG
           U=RHO_CM(1)/ENUP(50)
           V=RHO_CM(2)/ENUP(50)
           FI=ATAN(V/U)*RADDEG
           IF(U.LT.0.AND.V.GT.0) FI=FI+180.
           IF(U.LT.0.AND.V.LT.0) FI=FI+180.
           IF(U.GT.0.AND.V.LT.0) FI=FI+360.
           ENUP(52) = FI
          ENDIF    
	 ENDIF
       	 
	IF(JCH.EQ.11) THEN 
	 ENUP(53) = P_REST(4)
	 ENUP(54) = P_REST(5)
          IF(ENUP(54).NE.0.) THEN
	   ENUP(55) = ACOS(P_REST(3)/ENUP(54))*RADDEG
           U=P_REST(1)/ENUP(54)
           V=P_REST(2)/ENUP(54)
           FI=ATAN(V/U)*RADDEG
           IF(U.LT.0.AND.V.GT.0) FI=FI+180.
           IF(U.LT.0.AND.V.LT.0) FI=FI+180.
           IF(U.GT.0.AND.V.LT.0) FI=FI+360.
           ENUP(56) = FI
          ENDIF    
	 ENDIF
	ENUP(65) = Q2
	ENUP(66) = OMEGA
	ENUP(67) = JCH
	ENUP(68) = DNUP(NVAR-1)
	ENUP(69) = DNUP(NVAR)
c+ Radiative Gamma variables
        cc=acos(-1.0)/180.0
        E_TOTAL=E_TOTAL+e_radgam
        X_TOTAL=X_TOTAL+e_radgam*sin(theta_radgam*cc)*cos(phi_radgam*cc)
        Y_TOTAL=Y_TOTAL+e_radgam*sin(theta_radgam*cc)*sin(phi_radgam*cc)
        Z_TOTAL=Z_TOTAL+e_radgam*cos(theta_radgam*cc)
	ENUP(70) = e_radgam
	ENUP(71) = e_radgam
	ENUP(72) = theta_radgam
	ENUP(73) = phi_radgam

c++++ Check energy conservation

      
c      print *,' E=',abs(EO_beam + MP - E_TOTAL)
c      print *,' X=',abs(X_TOTAL) 
c      print *,' Y=',abs(Y_TOTAL) 
c      print *,' Z=',abs(EO_beam - Z_TOTAL) 
      IF( abs(EO_beam + MP - E_TOTAL) .gt. 0.001 ) THEN
        print *,' ENUP:EO_beam+MP,E_TOTAL=',EO_beam+MP,E_TOTAL,E_radgam
        print *,' Skipped'
        GOTO 99
      ENDIF


c+++++ Filling the ntuple

c          if(JCH.eq.15) print *,' --- 15 --- 555',NP
c          if(JCH.eq.15) print *,' --- 15 --- ENUP(67)',ENUP(67)
c          do j=1,67
c            print *,' NCODE1(j),ENUP(j)==',NCODE1(j),ENUP(j)
c          enddo
c          print *,' HFN'
c          pause
         
	 CALL HFN(1,DNUP)
	 CALL HFN(2,ENUP)
         


c-----
c----
c++++ Writing BOS file 
c+++++ Decomment next line to compile on INFNGE
c	GOTO 912
C+++++ 	Writing BOS header
        IHEAD = NBANK("HEAD",0,8,1)           ! BOS HEADER
        IW(iHEAD+1)=1
        IW(iHEAD+2)=RUN
        IW(iHEAD+3)=I
*        IW(iHEAD+4)=time()
        IW(iHEAD+4)=0
        IW(iHEAD+5)=-2
        IW(iHEAD+6)=0
        IW(iHEAD+7)=7
        IW(iHEAD+8)=1
c-----
C+++++ 	Writing	other MonteCarlo BOS banks	
	NR   = MEVT

        IMCEV = NBANK('MCEV',0,2,1)
        iw(iMCEV+1) = rran()*100000
        iw(iMCEV+2) = rran()*100000

        vert_r = rran()*  TGRAD ! Vertex radius
        vert_phi = rran()*  6.2832 ! Vertex phi
	vert_x = vert_r*cos(vert_phi)
	vert_y = vert_r*sin(vert_phi)
        vert_z = (rran()-.5)* TGLGT - TGOFFZ ! Target lenght + offset
	IMCVX = NBANK('MCVX',0,5,1)
	rw(iMCVX+1) = vert_x
	rw(iMCVX+2) = vert_y
	rw(iMCVX+3) = vert_z
	rw(iMCVX+4) = 0.
	iw(iMCVX+5) = 0
c-----

        IF( e_radgam.gt.0.010 ) THEN
        IMCTK = NBANK('MCTK',0,11,NP_FINAL+1+1)
        ELSE
        IMCTK = NBANK('MCTK',0,11,NP_FINAL+1)
        ENDIF
        E_TOTAL = 0. ! test the energy conservation

c+++++ Writing electron variables
           E_TOTAL = E_TOTAL+dnup(5)
           rw(imctk+1)=dnup(1)/dnup(5)
           rw(imctk+2)=dnup(2)/dnup(5)
           rw(imctk+3)=dnup(3)/dnup(5)
           rw(imctk+4)=dnup(5)
           rw(imctk+5)=me
           rw(imctk+6)=-1.
	   iw(imctk+7)=11                     ! e-
           iw(imctk+8)=0
           iw(imctk+9)=1
           iw(imctk+10)=0
           iw(imctk+11)=JCH

c ***********electron

	   call ElectronFiducialCut(DNUP(5),DNUP(6),DNUP(7) ,acc_el)
	   if(acc_el.eq.1.) iw(imctk+8)=1.
	   
c-----
c+++++ Looping and writing final state particles variables           
	   acc_pro=0.
	   acc_pip=0.
	   acc_pim=0.
	DO J_J = 1,NP_FINAL
	   K_K = J_J*7
	   IMCTK = IMCTK + 11
           E_TOTAL = E_TOTAL + DNUP(K_K+4)
           rw(imctk+1)=dnup(K_K+1)/dnup(K_K+5)
           rw(imctk+2)=dnup(K_K+2)/dnup(K_K+5)
           rw(imctk+3)=dnup(K_K+3)/dnup(K_K+5)
           rw(imctk+4)=dnup(K_K+5)
           rw(imctk+5)=SQRT(DNUP(K_K+4)**2-DNUP(K_K+5)**2)
           rw(imctk+6)=CHARGE(PART_ID(J_J))
           iw(imctk+7)=ID_BOS(PART_ID(J_J))
           iw(imctk+8)=0
           iw(imctk+9)=1
           iw(imctk+10)=0
           iw(imctk+11)=JCH

************DEFINING ACCEPTANCE FOR P P+ P-************


*** parameter in the bos bank acc (iw(imctk+8)): 

	   IF (ID_BOS(PART_ID(J_J)).EQ.2212) THEN      !proton
	     call HadronFiducialCut(DNUP(K_K+5),DNUP(K_K+6),DNUP(K_K+7) ,acc_had)
	     if(acc_had.eq.1.)then	 
	       acc_pro=1.
	       iw(imctk+8)=1.
	    endif
	  ENDIF
	   IF (ID_BOS(PART_ID(J_J)).EQ.211) THEN    !   pi+
	     call HadronFiducialCut(DNUP(K_K+5),DNUP(K_K+6),DNUP(K_K+7),acc_had)
	     if(acc_had.eq.1.)then
	       acc_pip=1.
	       iw(imctk+8)=1.
	    endif
	  ENDIF

	   IF (ID_BOS(PART_ID(J_J)).EQ.-211) THEN   !pi-
	     call HadronFiducialCut(DNUP(K_K+5),DNUP(K_K+6),DNUP(K_K+7),acc_had)
	     if(acc_had.eq.1.)then
	       acc_pim=1.
	       iw(imctk+8)=1.
	    endif
	  ENDIF



        ENDDO
c+++++ Writing final state Radiation Gamma       
        IF( e_radgam.gt.0.010 ) THEN
	   IMCTK = IMCTK + 11
           E_TOTAL = E_TOTAL+e_radgam
           rw(imctk+1)=sin(theta_radgam)*cos(phi_radgam)
           rw(imctk+2)=sin(theta_radgam)*sin(phi_radgam)
           rw(imctk+3)=cos(theta_radgam)
           rw(imctk+4)=e_radgam
           rw(imctk+5)=0.0
           rw(imctk+6)=0.0
           iw(imctk+7)=22
           iw(imctk+8)=0
           iw(imctk+9)=1
           iw(imctk+10)=0
           iw(imctk+11)=0
        ENDIF

c++++ Check energy conservation
      IF( abs(EO_beam + MP - E_TOTAL) .gt. 0.001 ) THEN
        print *,' BOS:EO_beam+MP,E_TOTAL=',EO_beam+MP,E_TOTAL,E_radgam
        print *,' Skipped'
        GOTO 99
      ENDIF

c******Writing file_acc.bos (all particles and their acc values)
c   *******  file used to calculate the geometrical acceptance table

	CALL FWBOS(IW,IBID_acc,'E',iSTATUS1)
c-----
c-----
* ACCEPTING REACTIONS 
** aucune condition 
	if(spacode.eq.0)goto 199

**  condition demandee
	if(spacode.eq.1)then
**test des differentes particules
c	   write(6,*)'acc',acc_el,acc_pro,acc_pip
	   if(acc_el.eq.1. .and. acc_pro.eq.1. .and. acc_pip.eq.1.)then
c	   if(acc_el.eq.1. .and. acc_pip.eq.1.)then
	      goto 199
	      else
	     goto 99
	   endif
	endif
***********
c  *****   writing bos file for accepted events

199	CALL FWBOS(IW,IBID,'E',iSTATUS1)
         CALL BDROP(IW,'E')
         CALL BGARB(IW)
912	 CONTINUE 
C-----        
c----
c++++ Creating LUND file
	 IF(LUND.EQ.1) CALL LUND_WRITE(PCM,NP_FINAL)
c---- Closing LUND if-loop
	ENDDO
 9999   CONTINUE
c-- Closing loop over number of events
c++ Closing paw ntuple
 	CALL HROUT(1,ICYCLE,' ')
	CALL HROUT(2,ICYCLE,' ') 	
	CALL HREND('ESCA')
c--
c++ Closing LUND file (if necessary)
	IF(LUND.EQ.1) CLOSE(2)
c--
c++ Closing BOS file
c++ Decomment next line to compile on INFNGE
c	GOTO 913
        CALL FWBOS(IW,IBID,'0',iSTATUS1)	
	print *,' end write status ',iSTATUS1
        CALL FCLOS()
913	CONTINUE	
c--
c++
c	WRITE(*,*) 'Total estractions:',N_TOT
	write(*,*) 'DONE!!!!'
	
	
	CALL RLUXUT(IZ)
	OPEN(UNIT=91,FILE='seed.dat',STATUS='OLD',FORM='FORMATTED')
	WRITE(91,*) IZ			! RANDOM NUMBER
	CLOSE(91)
	
	STOP


	END


C	*************************
        FUNCTION PROD_SCAL(A,B)
C	*************************
	IMPLICIT NONE
	
        REAL A(4),B(4), PROD_SCAL
        INTEGER J

	PROD_SCAL=0.
        DO J=1,3
          PROD_SCAL=PROD_SCAL+A(J)*B(J)
        END DO
        RETURN
        END


C	*************************
	FUNCTION RMASS(P)
C	*************************
	IMPLICIT NONE
	REAL RMASS,Q
	REAL RRAN
	CHARACTER*12 P
	INTEGER*4 IZ
	COMMON/RANDOM/IZ
	Q=rran()
	RMASS = -1.
	IF(P.EQ.'GAMMA')      RMASS = 0.0
	IF(P.EQ.'NEUTRINO')   RMASS = 0.0
	IF(P.EQ.'NULLA')      RMASS = 0.0
	IF(P.EQ.'ELECTRON')   RMASS = 0.511E-03
	IF(P.EQ.'POSITRON')   RMASS = 0.511E-03
	IF(P.EQ.'MU+')        RMASS = 0.10566
	IF(P.EQ.'MU-')        RMASS = 0.10566
	IF(P.EQ.'PROTON')     RMASS = 0.93827231
	IF(P.EQ.'NEUTRON')    RMASS = 0.93956563
	IF(P.EQ.'LAMDA')      RMASS = 1.116
	IF(P.EQ.'SIGMA+')     RMASS = 1.189
	IF(P.EQ.'SIGMA0')     RMASS = 1.192
	IF(P.EQ.'SIGMA-')     RMASS = 1.197
	IF(P.EQ.'CSI0')       RMASS = 1.3149
	IF(P.EQ.'CSI-')       RMASS = 1.3213
        IF(P.EQ.'DIBARYON')   RMASS =
     +  2.23 + Q*0.0165/2.355
	IF(P.EQ.'DELTA++') THEN
5          CALL GAUSS(0.0575,1.235,RMASS)
	  IF(RMASS.LE.1.08) GOTO 5 
	ENDIF
        IF(P.EQ.'DELTA+') THEN
6          CALL GAUSS(0.0575,1.235,RMASS)
           IF(RMASS.LE.1.08) GOTO 6 
	 ENDIF
	IF(P.EQ.'DELTA0')  THEN
7       CALL GAUSS(0.0575,1.235,RMASS)
           IF(RMASS.LE.1.08) GOTO 7   
	ENDIF
	IF(P.EQ.'DELTA-') THEN
8	CALL GAUSS(0.0575,1.235,RMASS)
           IF(RMASS.LE.1.08) GOTO 8
	ENDIF

	IF(P.EQ.'SIGMA+3/2')  RMASS = 1.382
	IF(P.EQ.'SIGMA03/2')  RMASS = 1.387
	IF(P.EQ.'SIGMA-3/2')  RMASS = 1.382
	IF(P.EQ.'CSI03/2')    RMASS = 1.532
	IF(P.EQ.'CSI-3/2')    RMASS = 1.535
	IF(P.EQ.'OMEGA-')     RMASS = 1.67245
	IF(P.EQ.'OMEGA') THEN
4        CALL GAUSS(0.0042,0.783,RMASS)
	 IF (RMASS.LT.0.415) GOTO 4  
	ENDIF
	IF(P.EQ.'PI+')        RMASS = 0.1395679
	IF(P.EQ.'PI0')        RMASS = 0.1349743
        IF(P.EQ.'PI-')        RMASS = 0.1395679
	IF(P.EQ.'KAPPA+')     RMASS = 0.49367
	IF(P.EQ.'KAPPA0')     RMASS = 0.49772
	IF(P.EQ.'KAPPA-')     RMASS = 0.49367
	IF(P.EQ.'KAPPAS')     RMASS = 0.49772
	IF(P.EQ.'KAPPAL')     RMASS = 0.49772
	IF(P.EQ.'AKAPPA0')    RMASS = 0.49772
	IF(P.EQ.'ETA')        RMASS = 0.5488
	IF(P.EQ.'ETAP')       RMASS = 0.95747
	IF(P.EQ.'RHO+') THEN
1	 CALL GAUSS(0.074,0.765,RMASS) 
	 IF (RMASS.LT.0.279) GOTO 1  
	ENDIF      
	IF(P.EQ.'RHO0') THEN
2        CALL GAUSS(0.074,0.77,RMASS)
         IF (RMASS.LT.0.279) GOTO 2  
	ENDIF
	IF(P.EQ.'RHO-')        THEN
3        CALL GAUSS(0.074,0.765,RMASS) 
         IF (RMASS.LT.0.279) GOTO 3  
        ENDIF      
	IF(P.EQ.'PHI')        RMASS = 1.020
	IF(P.EQ.'DEUTERON')   RMASS = 1.8755
	IF(P.EQ.'DEUTHEAVY')  RMASS = 1.878
	IF(P.EQ.'H3')         RMASS = 2.80875
	IF(P.EQ.'HE3')        RMASS = 2.80822
	IF(P.EQ.'HE4')        RMASS = 3.72715
	IF(P.EQ.'LI7')        RMASS = 6.53343
	IF(P.EQ.'BE9')        RMASS = 8.39225
	IF(P.EQ.'B11')        RMASS = 10.25191
	IF(P.EQ.'C12')        RMASS = 11.17418
C	IF(P.EQ.'C12STAR')    RMASS = 11.190
	IF(P.EQ.'O16')        RMASS = 14.89417
	IF(P.EQ.'F19')        RMASS = 17.69122
	IF(P.EQ.'MG24')       RMASS = 22.33506
	IF(P.EQ.'NA23')       RMASS = 24.40738
	IF(P.EQ.'AL27')       RMASS = 25.1257
	IF(P.EQ.'SI28')       RMASS = 26.05426
	IF(P.EQ.'K39')        RMASS = 36.28288
	IF(P.EQ.'CA40')       RMASS = 37.21355
	IF(P.EQ.'TI48')       RMASS = 44.64763
	RETURN
	END


C	*************************
	FUNCTION NP_CODE(P)
C	*************************
	IMPLICIT NONE
	INTEGER NP_CODE
	CHARACTER*12 P
	IF(P.EQ.'GAMMA')      NP_CODE=6
	IF(P.EQ.'NEUTRINO')   NP_CODE=0
	IF(P.EQ.'NULLA')      NP_CODE=0
	IF(P.EQ.'ELECTRON')   NP_CODE=7
17	IF(P.EQ.'POSITRON')   NP_CODE=-7
	IF(P.EQ.'MU+')        NP_CODE=-9
	IF(P.EQ.'MU-')        NP_CODE=9
	IF(P.EQ.'PROTON')     NP_CODE=41
	IF(P.EQ.'NEUTRON')    NP_CODE=42
	IF(P.EQ.'LAMDA')      NP_CODE=0
	IF(P.EQ.'SIGMA+')     NP_CODE=43
	IF(P.EQ.'SIGMA0')     NP_CODE=44
	IF(P.EQ.'SIGMA-')     NP_CODE=45
	IF(P.EQ.'CSI0')       NP_CODE=0
	IF(P.EQ.'CSI-')       NP_CODE=0
        IF(P.EQ.'DIBARYON')   NP_CODE=0
	IF(P.EQ.'DELTA')      NP_CODE=0
	IF(P.EQ.'DELTA++')    NP_CODE=61
	IF(P.EQ.'DELTA+')     NP_CODE=62
	IF(P.EQ.'DELTA0')     NP_CODE=63
	IF(P.EQ.'DELTA-')     NP_CODE=64
	IF(P.EQ.'SIGMA+3/2')  NP_CODE=0
	IF(P.EQ.'SIGMA03/2')  NP_CODE=0
	IF(P.EQ.'SIGMA-3/2')  NP_CODE=0
	IF(P.EQ.'CSI03/2')    NP_CODE=46
	IF(P.EQ.'CSI-3/2')    NP_CODE=47
	IF(P.EQ.'OMEGA-')     NP_CODE=0
	IF(P.EQ.'OMEGA')      NP_CODE=34
	IF(P.EQ.'PI+')        NP_CODE=17
	IF(P.EQ.'PI0')        NP_CODE=23
	IF(P.EQ.'PI-')        NP_CODE=-17
	IF(P.EQ.'KAPPA+')     NP_CODE=18
	IF(P.EQ.'KAPPA0')     NP_CODE=19
	IF(P.EQ.'KAPPA-')     NP_CODE=-18
	IF(P.EQ.'KAPPAS')     NP_CODE=37
	IF(P.EQ.'KAPPAL')     NP_CODE=38
	IF(P.EQ.'AKAPPA0')    NP_CODE=-19
	IF(P.EQ.'ETA')        NP_CODE=24
	IF(P.EQ.'ETAP')       NP_CODE=25
	IF(P.EQ.'RHO+')       NP_CODE=27
	IF(P.EQ.'RHO0')       NP_CODE=33
	IF(P.EQ.'RHO-')       NP_CODE=-27
	IF(P.EQ.'PHI')        NP_CODE=35
	IF(P.EQ.'DEUTERON')   NP_CODE=91
	IF(P.EQ.'DEUTHEAVY')  NP_CODE=0
	IF(P.EQ.'H3')         NP_CODE=92
	IF(P.EQ.'HE3')        NP_CODE=93
	IF(P.EQ.'HE4')        NP_CODE=94
	IF(P.EQ.'LI7')        NP_CODE=0
	IF(P.EQ.'BE9')        NP_CODE=0
	IF(P.EQ.'B11')        NP_CODE=0
	IF(P.EQ.'C12')        NP_CODE=0
C	IF(P.EQ.'C12STAR')    NP_CODE=0
	IF(P.EQ.'O16')        NP_CODE=0
	IF(P.EQ.'F19')        NP_CODE=0
	IF(P.EQ.'MG24')       NP_CODE=0
	IF(P.EQ.'NA23')       NP_CODE=0
	IF(P.EQ.'AL27')       NP_CODE=0
	IF(P.EQ.'SI28')       NP_CODE=0
	IF(P.EQ.'K39')        NP_CODE=0
	IF(P.EQ.'CA40')       NP_CODE=0
	IF(P.EQ.'TI48')       NP_CODE=0
	RETURN
	END

C*******************************   ROTATION   *********************************C
C------------------------------------------------------------------------------C
C	ROTATES A VECTOR FROM THE FRAME WITH THE Z AXIS ALONG THE  	       C
C	FLIGHT DIRECTION OF THE FLYING SYSTEM TO THE LAB SYSTEM                C
C									       C
C	U,V,W BEING THE DIR COS OF Z AXIS WITH RESPECT TO THE LAB              C
C									       C
C	THE OUTPUT VECTOR A(I) IS OVERWRITTEN ON THE INPUT VECTOR A(I)         C
C------------------------------------------------------------------------------C
C******************************************************************************C
	SUBROUTINE ROT(A,COSDIR)
C*******************************************************************************

	REAL A(4),B(4),DROT(3,3),COSDIR(3)
	REAL THETA,PHI

	PIGR=ACOS(-1.)
	U=COSDIR(1)
	V=COSDIR(2)
	W=COSDIR(3)

C**************************************************************
C	EVALUATION OF THETA, FI  FROM  U, V, W
C**************************************************************
	THETA = ACOS(W)
	IF(U.EQ.0.) THEN
	  PHI=PIGR*(2.-SIGN(1.,V))/2.
	ELSE
	  PHI =  ATAN(V/U)+PIGR*(1.-SIGN(1.,U))/2.
	  IF (PHI.LE.0.) PHI=PHI+2*PIGR	 
	ENDIF

C**************************************************************
C	EVALUATION OF THE COMPONENTS OF THE ROTATION MATRIX
C**************************************************************

	DROT(1,1) = - COS(THETA)*COS(PHI)
	DROT(1,2) = + COS(THETA)*SIN(PHI)
	DROT(1,3) = - SIN(THETA)

	DROT(2,1) = - SIN(PHI)
	DROT(2,2) = - COS(PHI)
	DROT(2,3) = + 0.

	DROT(3,1) = - SIN(THETA)*COS(PHI)
	DROT(3,2) = + SIN(THETA)*SIN(PHI)
	DROT(3,3) = + COS(THETA)

	DO L = 1,3
	 DO K = 1,3
	IF(ABS(DROT(L,K)-1.).LT.1.E-6) DROT(L,K) = 1.
	IF(ABS(DROT(L,K)+1.).LT.1.E-6) DROT(L,K) = -1.
	IF(ABS(DROT(L,K)).LT.1.E-4) DROT(L,K) = 0.
	 ENDDO
	ENDDO

C**************************************************************
C	EVALUATION OF THE SPACE COMPONENTS OF THE ROTATED VECTOR
C**************************************************************

	DO L = 1,3
	 B(L) = 0.
	 DO K = 1,3

	  B(L) = B(L) + A(K)*DROT(L,K)

	 ENDDO
	ENDDO

	DO L=1,3

	 A(L)=B(L)

	ENDDO

	RETURN

	END


	

C	*************************
	SUBROUTINE PCLE(JCH,TARGET,NP,PARTIC)
C	*************************
	IMPLICIT NONE

	INTEGER JCH,NP
	CHARACTER*12 PARTIC(10)
	CHARACTER*12 TARGET

	GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
     +  20,21,22,23,24,25,26,27,28,29) JCH

1	NP = 2
	TARGET    = 'PROTON'
	PARTIC(1) = 'NEUTRON'
	PARTIC(2) = 'PI+'
	RETURN

2	NP = 2
	TARGET    = 'PROTON'
	PARTIC(1) = 'PROTON'
	PARTIC(2) = 'PI0'
	RETURN

3	NP = 2
	TARGET    = 'NEUTRON'
	PARTIC(1) = 'PROTON'
	PARTIC(2) = 'PI-'
	RETURN

4	NP = 2
	TARGET    = 'NEUTRON'
	PARTIC(1) = 'NEUTRON'
	PARTIC(2) = 'PI0'
	RETURN

5	NP = 2
	TARGET    = 'PROTON'
	PARTIC(1) = 'DELTA++'
	PARTIC(2) = 'PI-'
	RETURN

6	NP = 2
	TARGET    = 'PROTON'
	PARTIC(1) = 'DELTA+'
	PARTIC(2) = 'PI0'
	RETURN

7	NP = 2
	TARGET    = 'PROTON'
	PARTIC(1) = 'DELTA0'
	PARTIC(2) = 'PI+'
	RETURN

8	NP = 2
	TARGET    = 'NEUTRON'
	PARTIC(1) = 'DELTA+'
	PARTIC(2) = 'PI-'
	RETURN

9	NP = 2
	TARGET    = 'NEUTRON'
	PARTIC(1) = 'DELTA0'
	PARTIC(2) = 'PI0'
	RETURN

10	NP = 2
	TARGET    = 'NEUTRON'
	PARTIC(1) = 'DELTA-'
	PARTIC(2) = 'PI+'
	RETURN

11	NP = 2
	TARGET    = 'PROTON'
	PARTIC(1) = 'RHO0'
	PARTIC(2) = 'PROTON'
	RETURN

12	NP = 2
	TARGET    = 'PROTON'
	PARTIC(1) = 'RHO+'
	PARTIC(2) = 'NEUTRON'
	RETURN


13      NP = 2
        TARGET    = 'NEUTRON'
        PARTIC(1) = 'RHO-'
        PARTIC(2) = 'PROTON'
        RETURN

14      NP = 2
        TARGET    = 'NEUTRON'
        PARTIC(1) = 'RHO0'
        PARTIC(2) = 'NEUTRON'
        RETURN
 
15      NP = 3
        TARGET    = 'PROTON'
        PARTIC(1) = 'PROTON'
        PARTIC(2) = 'PI+'
        PARTIC(3) = 'PI-'
        RETURN

16      NP = 3
        TARGET    = 'PROTON'
        PARTIC(1) = 'NEUTRON'
        PARTIC(2) = 'PI0'
        PARTIC(3) = 'PI+'
        RETURN

17      NP = 3
        TARGET    = 'NEUTRON'
        PARTIC(1) = 'NEUTRON'
        PARTIC(2) = 'PI+'
        PARTIC(3) = 'PI-'
        RETURN

18      NP = 3
        TARGET    = 'NEUTRON'
        PARTIC(1) = 'PROTON'
        PARTIC(2) = 'PI0'
        PARTIC(3) = 'PI-'
        RETURN

19      NP = 2
        TARGET    = 'PROTON'
        PARTIC(1) = 'OMEGA'
        PARTIC(2) = 'PROTON'
        RETURN

20      NP = 2
        TARGET    = 'NEUTRON'
        PARTIC(1) = 'OMEGA'
        PARTIC(2) = 'NEUTRON'
        RETURN

21      NP = 4
        TARGET    = 'PROTON'
        PARTIC(1) = 'PROTON' 
        PARTIC(2) = 'PI+'
        PARTIC(3) = 'PI-'
        PARTIC(4) = 'PI0'
        RETURN

22      NP = 4
        TARGET    = 'PROTON'
        PARTIC(1) = 'NEUTRON' 
        PARTIC(2) = 'PI+'
        PARTIC(3) = 'PI+'
        PARTIC(4) = 'PI-'
        RETURN

23      NP = 4
        TARGET    = 'NEUTRON'
        PARTIC(1) = 'NEUTRON'
        PARTIC(2) = 'PI+'
        PARTIC(3) = 'PI-'
        PARTIC(4) = 'PI0'
        RETURN 

24      NP = 4
        TARGET    = 'NEUTRON'
        PARTIC(1) = 'PROTON'
        PARTIC(2) = 'PI+'
        PARTIC(3) = 'PI-'
        PARTIC(4) = 'PI-'
        RETURN        

25	NP = 2
        TARGET    = 'PROTON'
        PARTIC(1) = 'PHI'
        PARTIC(2) = 'PROTON'
        RETURN        

26	NP = 2
        TARGET    = 'NEUTRON'
        PARTIC(1) = 'PHI'
        PARTIC(2) = 'NEUTRON'
        RETURN        
27      NP = 5
        TARGET    = 'PROTON'
        PARTIC(1) = 'PROTON'
        PARTIC(2) = 'PI+'
        PARTIC(3) = 'PI+'
        PARTIC(4) = 'PI-'
        PARTIC(5) = 'PI-'
        RETURN 
28	NP = 6
        TARGET    = 'PROTON'
        PARTIC(1) = 'PROTON'
        PARTIC(2) = 'PI+'
        PARTIC(3) = 'PI+'
        PARTIC(4) = 'PI-'
        PARTIC(5) = 'PI-'
        PARTIC(6) = 'PI0'
        RETURN 
29	NP = 6
        TARGET    = 'PROTON'
        PARTIC(1) = 'NEUTRON'
        PARTIC(2) = 'PI+'
        PARTIC(3) = 'PI+'
        PARTIC(4) = 'PI+'
        PARTIC(5) = 'PI-'
        PARTIC(6) = 'PI-'
        RETURN 
        END



c-
C	*************************
	SUBROUTINE SEZTOT
C	*************************



	IMPLICIT NONE
	REAL SIGR(40,251),WR(251)
	INTEGER I,J,K
	integer  iflagphi
        COMMON/SEZURTO/SIGR,WR
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C	TABELLA CANALE/REAZIONE
C
C         1) GAMMA  +  P  --->  N  +  PI+
C         2) GAMMA  +  P  --->  P  +  PI0
C         3) GAMMA  +  N  --->  P  +  PI-
C         4) GAMMA  +  N  --->  N  +  PI0
C         5) GAMMA  +  P  --->  DELTA++ +  PI-
C         6) GAMMA  +  P  --->  DELTA+  +  PI0
C         7) GAMMA  +  P  --->  DELTA0  +  PI+
C         8) GAMMA  +  N  --->  DELTA+  +  PI-
C         9) GAMMA  +  N  --->  DELTA0  +  PI0
C        10) GAMMA  +  N  --->  DELTA-  +  PI+
C        11) GAMMA  +  P  --->  RHO0  +  P
C        12) GAMMA  +  P  --->  RHO+  +  N
C        13) GAMMA  +  N  --->  RHO-  +  P
C        14) GAMMA  +  N  --->  RHO0  +  N
C        15) GAMMA  +  P  --->  P  +  PI   +  PI-  
C        16) GAMMA  +  P  --->  N  +  PI0  +  PI+   
C        17) GAMMA  +  N  --->  N  +  PI+  +  PI- 
C        18) GAMMA  +  N  --->  P  +  PI0  +  PI- 
C        19) GAMMA  +  P  --->  OMEGA  +  P
C        20) GAMMA  +  N  --->  OMEGA  +  N
C        21) GAMMA  +  P  --->  P  +  PI+  +  PI-  +  PI0  
C        22) GAMMA  +  P  --->  N  +  PI+  +  PI+  +  PI-
C        23) GAMMA  +  N  --->  N  +  PI+  +  PI-  +  PI0 
C        24) GAMMA  +  N  --->  P  +  PI+  +  PI-  +  PI- 
C        25) GAMMA  +  P  --->  PHI  +  P
C        26) GAMMA  +  N  --->  PHI  +  N
c        27) gamma  +  p  --->  pi+  +  pi+  +  pi-  +  pi-  +  p                   
c        28) gamma  +  p  --->  pi+  +  pi+  +  pi-  +  pi-  +  pi0  +  p                
c        29) gamma  +  p  --->  pi+  +  pi+  +  pi+  +  pi-  +  pi-  +  n                

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	DATA ((SIGR(I,J),J=1,251),I=1,1) / 0.00,   0.00,   0.00,   0.00,   0.00,
     &	           0.00,   0.00,   0.00,  18.16,  36.33,  54.49,  72.66, 
     &            90.82, 108.98, 127.15, 145.31, 163.48, 181.64, 200.88,
     &           220.41, 234.94, 239.35, 231.26, 213.19, 190.59, 168.22,
     &		 148.61, 132.51, 119.68, 109.55, 101.58,  95.34,  90.52,
     &		  86.87,  84.21,  82.40,  81.28,  80.76,  80.71,  81.07,
     &		  81.82,  82.98,  84.61,  86.78,  89.57,  92.97,  96.79,
     &		 100.51, 103.21, 103.70, 100.93,  94.73,  85.99,  76.31,
     &		  67.33,  60.14,  55.06,  51.88,  50.22,  49.73,  50.13,
     &		  51.16,  52.60,  54.15,  55.47,  56.17,  55.87,  54.28,
     &		  51.30,  47.18,  42.44,  37.64,  33.23,  29.43,  26.30,
     &		  23.78,  21.78,  20.20,  18.93,  17.02,  17.13,  16.51,
     &		  16.07,  15.75,  15.55,  15.44,  15.38,  15.33,  15.24,
     &		  15.10,  14.87,  14.56,  14.16,  13.68,  13.15,  12.59,
     &		  12.02,  11.48,  10.96,  10.49,  10.06,   9.89,   9.71,
     &		   9.35,   9.05,   8.88,   8.71,   8.53,   8.47,   8.32,
     &		   8.24,   8.06,   7.86,   7.71,   7.45,   7.23,   7.15,
     &		   7.08,   6.98,   6.91,   6.83,   6.72,   6.51,   6.45,
     &		   6.37,   6.29,   6.20,   6.16,   6.04,   5.95,   5.81,
     &		   5.70,   5.59,   5.53,   5.49,   5.40,   5.37,   5.32,
     &		   5.24,   5.20,   5.13,   5.07,   4.99,   4.93,   4.82,
c     &		   4.69,   4.32,   3.96,   3.60,   3.23,   2.86,   2.48,
c     &		   2.11,   1.74,   1.37,   1.00,   1.13,   1.26,   1.38,
c     &		   1.51,   1.64,   1.77,   1.90,   1.98,   2.07,   2.15,
c     &		   2.24,   2.32,   2.40,   2.49,   2.52,   2.55,   2.58,
c     &		   2.61,   2.64,   2.67,   2.70,   2.67,   2.64,   2.61,
c     &		   2.58,   2.55,   2.52,   2.49,   2.42,   2.36,   2.30,
c     &		   2.23,   2.17,   2.10,   2.04,   1.97,   1.91,   1.85,
c     &		   1.78,   1.72,   1.65,   1.59,   1.52,   1.46,   1.40,
c     &		   1.33,   1.27,   1.20,   1.14,   1.07,   1.01,   0.95,
c     &		   0.88,   0.82,   0.75,   0.69,   0.62,   0.56,   0.50,
c     &		   0.43,   0.37,   0.30,   0.24,   0.17,   0.11,   0.00,
     &		   4.73,   4.64,   4.55,   4.46,   4.36,   4.28,   4.20,
     &		   4.13,   4.06,   3.99,   3.92,   3.86,   3.80,   3.74,
     &		   3.68,   3.62,   3.56,   3.50,   3.44,   3.39,   3.34,
     &		   3.29,   3.24,   3.19,   3.14,   3.09,   3.04,   3.00,
     &		   2.96,   2.92,   2.88,   2.84,   2.80,   2.76,   2.72,     
     &		   2.69,   2.66,   2.63,   2.60,   2.57,   2.54,   2.51,
     &		   2.48,   2.45,   2.43,   2.41,   2.39,   2.37,   2.35,
     &		   2.33,   2.31,   2.29,   2.27,   2.25,   2.23,   2.21,
     &		   2.20,   2.19,   2.18,   2.17,   2.15,   2.14,   2.13,
     &		   2.12,   2.11,   2.10,   2.09,   2.08,   2.07,   2.06,
     &		   2.05,   2.04,   2.03,   2.02,   2.01,   2.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00/    
	DATA ((SIGR(I,J),J=1,251),I=2,2) / 0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,  10.90,  21.80,  32.69,  43.59,
     &		  54.49,  65.39,  76.29,  87.19,  98.09, 108.98, 172.62,
     &		 218.89, 262.70, 293.39, 302.38, 289.31, 261.40, 227.60,
     &		 194.21, 164.31, 138.81, 117.59, 100.09,  85.70,  73.86,
     &		  64.11,  56.07,  49.42,  43.93,  39.41,  35.66,  32.58,
     &		  30.09,  28.16,  26.79,  26.02,  25.93,  26.63,  28.23,
     &		  30.71,  33.85,  37.06,  39.57,  40.76,  40.39,  38.63, 
     &		  35.87,  32.67,  29.55,  26.90,  24.91,  23.57,  22.83,
     &		  22.60,  22.78,  23.26,  23.90,  24.52,  24.93,  24.92,
     &		  24.36,  23.22,  21.62,  19.79,  17.95,  16.25,  14.79,
     &		  13.58,  12.61,  11.84,  11.26,  10.85,  10.59,  10.46,
     &		  10.44,  10.52,  10.68,  10.90,  11.16,  11.43,  11.66,
     &		  11.83,  11.90,  11.84,  11.63,  11.30,  10.85,  10.34,
     &		  10.28,  10.23,  10.19,  10.15,  10.09,  10.05,  10.00,
     &		   9.98,   9.95,   9.89,   9.80,   9.74,   9.71,   9.67,
     &		   9.65,   9.61,   9.58,   9.52,   9.47,   9.45,   9.40,
     &		   9.36,   9.33,   9.29,   9.22,   9.17,   9.12,   9.07, 
     &		   9.00,   8.97,   8.95,   8.90,   8.85,   8.81,   8.78,
     &		   8.75,   8.19,   8.06,   7.93,   7.81,   7.68,   7.55,
     &		   7.42,   7.30,   7.17,   6.72,   6.28,   5.83,   5.38,
     &		   4.94,   4.49,   4.04,   3.60,   3.37,   3.14,   2.91,
     &		   2.68,   2.46,   2.23,   2.00,   2.03,   2.06,   2.08,
     &		   2.11,   2.14,   2.17,   2.20,   2.24,   2.28,   2.32,
     &		   2.36,   2.41,   2.45,   2.49,   2.51,   2.52,   2.54,
     &		   2.55,   2.57,   2.58,   2.60,   2.57,   2.54,   2.51,
     &		   2.48,   2.46,   2.43,   2.40,   2.38,   2.36,   2.34,
     &		   2.31,   2.29,   2.27,   2.25,   2.23,   2.21,   2.18,
     &		   2.16,   2.14,   2.12,   2.10,   2.08,   2.06,   2.04,
     &		   2.01,   1.99,   1.97,   1.95,   1.93,   1.91,   1.88,
     &		   1.86,   1.84,   1.82,   1.80,   1.78,   1.76,   1.74,
     &		   1.71,   1.69,   1.67,   1.65,   1.63,   1.61,   1.58,
     &		   1.56,   1.54,   1.52,   1.50,   1.48,   1.46,   1.43,
     &		   1.41,   1.39,   1.37,   1.35,   1.33,   1.31,   1.28,
     &		   1.26,   1.24,   1.22,   1.20,   1.18,   1.16,   1.13,
     &		   1.11,   1.09,   1.07,   1.05,   1.03,   1.01,   0.98, 
     &		   0.96/
	DATA ((SIGR(I,J),J=1,251),I=3,3) / 0.00,   0.00,   0.00,   0.00,   0.00,
     &	           0.00,   0.00,   0.00,  18.16,  36.33,  54.49,  72.65,
     &		  90.82, 108.98, 127.15, 145.31, 163.48, 181.64, 255.09,
     &		 278.94, 295.07, 297.65, 284.80, 260.68, 232.43, 205.75,
     &		 183.33, 165.62, 152.02, 141.58, 133.44, 126.97, 121.71,
     &		 117.33, 113.62, 110.46, 107.78, 105.56, 103.82, 102.55,
     &		 101.53, 100.80, 100.43, 100.46, 100.88, 101.58, 102.25,
     &		 102.33, 100.99,  97.33,  90.78,  81.78,  71.50,  61.40,
     &		  52.63,  45.75,  40.75,  37.34,  35.12,  33.82,  33.15,
     &		  32.82,  32.63,  32.38,  31.92,  31.15,  30.06,  28.72,
     &		  27.33,  25.96,  24.66,  23.44,  22.30,  21.26,  20.29,
     &		  19.39,  18.55,  17.79,  17.09,  16.49,  16.01,  14.16,
     &		  14.08,  14.00,  13.91,  13.83,  13.74,  13.65,  13.56,
     &		  13.47,  13.38,  13.29,  13.19,  13.10,  13.00,  12.91,
     &		  12.83,  12.75,  12.66,  12.58,  12.50,  12.38,  12.23,
     &		  12.18,  12.10,  12.07,  12.00,  11.96,  11.92,  11.86,
     &		  11.75,  11.68,  11.63,  11.55,  11.46,  11.37,  11.31,
     &		  11.25,  11.19,  11.08,  10.97,  10.84,  10.78,  10.69,
     &		  10.60,  10.29,   9.97,   9.66,   9.34,   9.03,   8.72,
     &	           8.40,   8.09,   7.89,   7.69,   7.49,   7.29,   7.10,
     &		   6.90,   6.70,   6.50,   6.14,   5.77,   5.41,   5.05,
     &		   4.69,   4.32,   3.96,   3.60,   3.23,   2.86,   2.49,
     &		   2.11,   1.74,   1.37,   1.00,   1.13,   1.26,   1.38,
     &		   1.51,   1.64,   1.77,   1.90,   1.98,   2.07,   2.15,
     &		   2.24,   2.32,   2.40,   2.49,   2.52,   2.55,   2.58,
     &		   2.61,   2.64,   2.70,   2.70,   2.67,   2.64,   2.61,
     &		   2.58,   2.55,   2.52,   2.49,   2.43,   2.36,   2.30,
     &		   2.23,   2.17,   2.10,   2.04,   1.97,   1.91,   1.85,
     &		   1.78,   1.72,   1.65,   1.59,   1.53,   1.46,   1.40,
     &		   1.33,   1.27,   1.20,   1.14,   1.07,   1.01,   0.95,
     &		   0.88,   0.82,   0.75,   0.69,   0.63,   0.56,   0.50,
     &		   0.43,   0.37,   0.30,   0.24,   0.18,   0.11,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00/
	DATA ((SIGR(I,J),J=1,251),I=4,4) / 0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   8.54,  17.07,  25.61,  34.15,
     &		  42.68,  51.22,  59.76,  68.30,  76.83,  85.37, 155.30,
     &		 200.65, 246.11, 281.58, 297.90, 293.09, 272.74, 244.88,
     &		 215.65, 188.35, 164.22, 143.25, 125.13, 109.52,  96.07,
     &		  84.49,  74.53,  65.97,  58.66,  52.45,  47.20,  42.80,
     &		  39.06,  35.96,  33.49,  31.66,  30.51,  30.07,  30.35,
     &		  31.30,  32.68,  34.04,  34.87,  34.87,  34.01,  32.43, 
     &		  30.37,  28.09,  25.83,  23.75,  21.90,  20.25,  18.74,
     &		  17.30,  15.86,  14.42,  12.97,  11.56,  10.26,   9.13,
     &		   8.19,   7.43,   6.82,   6.33,   5.95,   5.66,   5.48,
     &		   5.38,   5.34,   5.35,   5.42,   5.56,   5.81,   6.17,
     &		   6.66,   7.31,   8.11,   9.05,  10.14,  11.33,  12.58,
     &	          12.03,  11.68,  11.32,  10.97,  10.61,  10.25,   9.88,
     &		   9.82,   9.78,   9.72,   9.69,   9.65,   9.57,   9.53,
     &		   9.49,   9.42,   9.48,   9.49,   9.53,   9.61,   9.65,
     &		   9.67,   9.70,   9.73,   9.76,   9.75,   9.72,   9.69,
     &		   9.67,   9.65,   9.59,   9.53,   9.49,   9.45,   9.30,
     &		   8.89,   8.49,   8.10,   7.71,   7.31,   6.92,   6.52,
     &		   6.13,   5.73,   5.64,   5.55,   5.46,   5.38,   5.29,
     &		   5.20,   5.11,   5.02,   4.71,   4.39,   4.08,   3.77,
     &		   3.46,   3.14,   2.83,   2.52,   2.36,   2.20,   2.04,
     &		   1.88,   1.72,   1.56,   1.40,   1.42,   1.44,   1.46,
     &		   1.48,   1.50,   1.52,   1.54,   1.57,   1.60,   1.63,
     &		   1.66,   1.68,   1.71,   1.74,   1.75,   1.76,   1.77,
     &		   1.79,   1.80,   1.81,   1.82,   1.80,   1.78,   1.76,
     &		   1.74,   1.72,   1.70,   1.68,   1.66,   1.65,   1.63,
     &		   1.62,   1.60,   1.59,   1.57,   1.56,   1.54,   1.53,
     &		   1.51,   1.50,   1.48,   1.47,   1.45,   1.44,   1.42,
     &		   1.41,   1.39,   1.38,   1.36,   1.35,   1.33,   1.32,
     &		   1.30,   1.29,   1.27,   1.26,   1.24,   1.23,   1.21,
     &		   1.20,   1.18,   1.17,   1.15,   1.14,   1.12,   1.11,
     &		   1.09,   1.08,   1.06,   1.05,   1.03,   1.02,   1.00,
     &		   0.99,   0.97,   0.96,   0.94,   0.93,   0.91,   0.90,
     &		   0.88,   0.87,   0.85,   0.84,   0.82,   0.81,   0.79,
     &		   0.78,   0.76,   0.75,   0.73,   0.72,   0.70,   0.69,
     &		   0.67/
	DATA ((SIGR(I,J),J=1,251),I=5,5) / 0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.08,   0.19,   1.02,   1.73,
     &		   2.88,   3.95,   5.36,   6.82,   8.30,  10.37,  12.34,
     &		  15.33,  18.01,  21.88,  25.52,  30.35,  34.96,  40.28,
     &		  45.27,  50.18,  54.53,  58.24,  61.17,  63.19,  64.28,
     &		  64.36,  63.45,  61.67,  59.27,  56.74,  54.62,  53.38,
     &		  53.07,  53.50,  54.23,  55.00,  55.53,  55.83,  55.82,
     &		  55.56,  55.05,  54.35,  53.51,  52.52,  51.45,  47.50,
     &		  42.00,  39.20,  35.50,  33.00,  30.20,  28.50,  27.10,
     &		  26.30,  24.20,  22.30,  21.40,  20.60,  20.00,  19.60,
     &		  19.00,  18.00,  17.50,  17.00,  16.20,  15.80,  15.30,
     &		  14.60,  13.50,  12.80,  11.90,  11.00,  10.40,  10.00,
     &		   9.60,   9.40,   9.10,   8.80,   8.70,   8.50,   8.40,
     &		   8.40,   8.20,   8.00,   7.50,   7.00,   6.70,   6.30,
     &		   6.00,   5.80,   5.60,   5.40,   5.20,   5.00,   4.90,
     &		   4.90,   4.80,   4.70,   4.60,   4.60,   4.50,   4.40,
     &		   4.40,   4.30,   4.20,   4.20,   4.20,   4.10,   4.00,
     &		   4.00,   3.90,   3.90,   3.80,   3.80,   3.70,   3.70,
     &		   3.60,   3.60,   3.50,   3.50,   3.50,   3.50,   3.40,
     &		   3.40,   3.40,   3.40,   3.30,   3.30,   3.30,   3.30,
     &		   3.30,   3.30,   3.30,   3.20,   3.20,   3.20,   3.20,
     &		   3.20,   3.20,   3.20,   3.20,   3.15,   3.10,   3.06,
     &		   3.01,   2.97,   2.93,   2.89,   2.84,   2.80,   2.76,
     &		   2.73,   2.69,   2.65,   2.61,   2.58,   2.54,   2.51,
     &		   2.48,   2.44,   2.41,   2.38,   2.35,   2.32,   2.29,
     &		   2.26,   2.23,   2.20,   2.18,   2.15,   2.12,   2.10,
     &		   2.07,   2.05,   2.02,   2.00,   1.97,   1.95,   1.93,
     &		   1.90,   1.88,   1.86,   1.84,   1.82,   1.80,   1.78,
     &		   1.76,   1.74,   1.72,   1.70,   1.68,   1.66,   1.64,
     &		   1.63,   1.61,   1.59,   1.57,   1.56,   1.54,   1.52,
     &		   1.51,   1.49,   1.48,   1.46,   1.45,   1.43,   1.42,
     &		   1.40,   1.39,   1.38,   1.36,   1.35,   1.34,   1.32,
     &		   1.31,   1.30,   1.29,   1.27,   1.26,   1.25,   1.24,
     &		   1.23,   1.21,   1.20,   1.19,   1.18,   1.17,   1.16,
     &		   1.15/
	DATA ((SIGR(I,J),J=1,251),I=6,6) / 0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.01,   0.01,   0.02,   0.03,   0.04,   0.05,   0.06,
     &		   0.08,   0.11,   0.14,   0.19,   0.25,   0.33,   0.43,
     &		   0.55,   0.70,   0.89,   1.12,   1.40,   1.74,   2.16,
     &		   2.66,   3.23,   3.86,   4.48,   4.98,   5.25,   5.23,
     &		   4.94,   4.47,   3.93,   3.40,   2.93,   2.53,   2.21,
     &		   1.95,   1.76,   1.63,   1.54,   1.51,   1.51,   1.54,
     &		   1.58,   1.61,   1.61,   1.59,   1.53,   1.45,   1.37,
     &		   1.27,   1.18,   1.10,   1.02,   0.95,   0.89,   0.83,
     &		   0.78,   0.73,   0.69,   0.66,   0.62,   0.59,   0.57,
     &		   0.54,   0.52,   0.50,   0.48,   0.46,   0.44,   0.43,
     &		   0.41,   0.40,   0.39,   0.38,   0.37,   0.36,   0.35,
     &		   0.34,   0.33,   0.32,   0.31,   0.31,   0.30,   0.29,
     &		   0.29,   0.28,   0.28,   0.27,   0.26,   0.26,   0.26, 
     &		   0.25,   0.25,   0.24,   0.24,   0.23,   0.23,   0.23,
     &		   0.22,   0.22,   0.22,   0.21,   0.21,   0.21,   0.21,
     &		   0.20,   0.20,   0.20,   0.20,   0.19,   0.19,   0.19,
     &		   0.19,   0.18,   0.18,   0.18,   0.18,   0.18,   0.18,
     &		   0.17,   0.17,   0.17,   0.17,   0.17,   0.17,   0.16,
     &		   0.16,   0.16,   0.16,   0.16,   0.15,   0.14,   0.14,
     &		   0.13,   0.13,   0.13,   0.12,   0.12,   0.12,   0.11,
     &		   0.11,   0.11,   0.10,   0.10,   0.10,   0.09,   0.09,
     &		   0.09,   0.09,   0.09,   0.08,   0.08,   0.08,   0.08,
     &		   0.07,   0.07,   0.07,   0.07,   0.07,   0.07,   0.06,
     &		   0.06,   0.06,   0.06,   0.06,   0.06,   0.06,   0.05,
     &		   0.05,   0.05,   0.05,   0.05,   0.05,   0.05,   0.05,
     &		   0.05,   0.04,   0.04,   0.04,   0.04,   0.04,   0.04,
     &		   0.04,   0.04,   0.04,   0.04,   0.04,   0.04,   0.03,
     &		   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,
     &		   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,
     &		   0.03,   0.03,   0.02,   0.02,   0.02,   0.02,   0.02,
     &		   0.02,   0.02,   0.02,   0.02,   0.02,   0.02,   0.02,
     &		   0.02,   0.02,   0.02,   0.02,   0.02,   0.02,   0.02,
     &		   0.02/
	DATA ((SIGR(I,J),J=1,251),I=7,7) /251*1./
	DATA ((SIGR(I,J),J=1,251),I=8,8) / 0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.10,   0.23,   0.36,   0.44,   0.57,
     &		   0.71,   0.98,   1.36,   1.98,   2.74,   3.20,   4.10,
     &		   5.40,   6.30,   7.32,   8.50,   9.80,  11.60,  13.40,
     &		  15.71,  16.80,  17.60,  20.30,  22.10,  24.41,  25.80,
     &		  26.70,  28.50,  29.70,  31.15,  31.16,  31.17,  31.18,
     &		  31.19,  31.19,  30.00,  29.60,  27.00,  26.00,  25.58,
     &		  24.60,  23.50,  23.00,  22.00,  21.31,  20.80,  20.20,
     &		  19.50,  18.80,  18.08,  17.40,  16.90,  16.30,  15.90,
     &		  15.44,  14.80,  14.00,  13.80,  13.30,  13.20,  12.80,
     &		  12.20,  12.00,  11.70,  11.34,  11.00,  10.50,  10.30,
     &		  10.00,   9.80,   9.40,   9.10,   8.90,   8.70,   8.53,
     &		   8.30,   8.10,   7.90,   7.60,   7.47,   7.00,   6.90,
     &		   6.75,   6.66,   6.58,   6.35,   6.10,   6.00,   5.91,
     &		   5.83,   5.60,   5.50,   5.30,   5.20,   5.20,   5.10,
     &		   5.00,   4.80,   4.70,   4.66,   4.50,   4.40,   4.30,
     &		   4.22,   4.19,   4.10,   4.00,   3.90,   3.82,   3.79,
     &		   3.65,   3.52,   3.50,   3.47,   3.44,   3.38,   3.32,
     &		   3.28,   3.20,   3.14,   3.08,   3.02,   2.97,   2.92,
     &		   2.87,   2.83,   2.80,   2.75,   2.70,   2.64,   2.60,
     &		   2.55,   2.50,   2.47,   2.43,   2.40,   2.37,   2.33,
     &		   2.28,   2.24,   2.22,   2.18,   2.14,   2.10,   2.08,
     &		   2.06,   2.40,   2.02,   2.00,   1.93,   1.90,   1.88,
     &		   1.85,   1.82,   1.80,   1.78,   1.73,   1.70,   1.69,
     &		   1.68,   1.66,   1.64,   1.62,   1.60,   1.57,   1.53,
     &		   1.51,   1.49,   1.48,   1.47,   1.45,   1.43,   1.41,
     &		   1.39,   1.38,   1.36,   1.34,   1.32,   1.31,   1.30,
     &		   1.25,   1.20,   1.15,   1.10,   1.05,   1.00,   1.00,
     &		   1.00,   1.00,   1.00,   1.00,   1.00,   1.00,   1.00,
     &		   1.00,   1.00,   1.00,   1.00,   1.00,   1.00,   1.00,
     &		   1.00,   1.00,   1.00,   1.00,   1.00,   1.00,   1.00,
     &		   1.00,   1.00,   1.00,   1.00,   1.00,   1.00,   1.00,
     &		   1.00,   1.00,   1.00,   1.00,   1.00,   1.00,   1.00,
     &		   1.00,   1.00,   1.00,   1.00,   1.00,   1.00,   1.00,
     &		   1.00/
	DATA ((SIGR(I,J),J=1,251),I=9,9) / 0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.01,   0.02,   0.02,
     &		   0.03,   0.04,   0.05,   0.08,   0.11,   0.14,   0.17,
     &		   0.19,   0.24,   0.35,   0.44,   0.57,   0.62,   0.78,
     &		   0.89,   1.02,   1.31,   1.59,   1.59,   1.59,   1.59,
     &		   1.59,   1.59,   1.45,   1.21,   1.04,   0.89,   0.78,
     &		   0.67,   0.63,   0.60,   0.56,   0.53,   0.54,   0.54,
     &		   0.55,   0.56,   0.56,   0.52,   0.49,   0.47,   0.46,
     &		   0.44,   0.42,   0.40,   0.37,   0.33,   0.31,   0.29,
     &		   0.27,   0.26,   0.24,   0.23,   0.22,   0.21,   0.20,
     &		   0.19,   0.18,   0.18,   0.17,   0.17,   0.16,   0.15,
     &		   0.15,   0.14,   0.14,   0.13,   0.12,   0.12,   0.12,
     &		   0.11,   0.11,   0.11,   0.10,   0.10,   0.10,   0.10,
     &		   0.10,   0.10,   0.09,   0.09,   0.09,   0.09,   0.08,
     &		   0.08,   0.08,   0.08,   0.08,   0.08,   0.08,   0.08,
     &		   0.08,   0.07,   0.07,   0.07,   0.07,   0.07,   0.07,
     &		   0.07,   0.07,   0.07,   0.07,   0.06,   0.06,   0.06,
     &		   0.06,   0.06,   0.06,   0.06,   0.06,   0.06,   0.06,
     &		   0.06,   0.06,   0.06,   0.06,   0.06,   0.05,   0.05,
     &		   0.05,   0.05,   0.05,   0.05,   0.05,   0.05,   0.05,
     &		   0.05,   0.05,   0.05,   0.05,   0.05,   0.05,   0.05,
     &		   0.05,   0.05,   0.05,   0.05,   0.04,   0.04,   0.04,  
     &		   0.04,   0.04,   0.04,   0.04,   0.04,   0.04,   0.04,
     &		   0.04,   0.04,   0.04,   0.04,   0.04,   0.04,   0.04,
     &		   0.04,   0.04,   0.04,   0.04,   0.04,   0.04,   0.04,
     &		   0.04,   0.04,   0.04,   0.04,   0.04,   0.04,   0.04,
     &		   0.04,   0.04,   0.04,   0.04,   0.04,   0.04,   0.04,
     &		   0.04,   0.04,   0.03,   0.03,   0.03,   0.03,   0.03,
     &		   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,
     &		   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,
     &		   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,
     &		   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,
     &		   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,
     &		   0.03/
	DATA ((SIGR(I,J),J=1,251),I=10,10) /251*1./
	DATA ((SIGR(I,J),J=1,251),I=11,11) /0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &	           0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &             0.00,   0.00,   0.00,   0.00,   0.23,   0.23,   0.24,
     &	           0.26,   0.28,   0.29,   0.31,   0.33,   0.35,   0.38,
     &	           0.82,   1.38,   1.79,   2.18,   2.38,   2.58,   3.02,
     &		   3.50,   4.03,   4.57,   5.04,   5.36,   5.42,   5.24,
     &             5.87,   5.74,   5.44,   5.13,   4.86,   4.64,   4.50,
     &	           4.42,   4.41,   4.47,   7.13,   8.59,   9.74,  10.71,
     &	          11.48,  11.75,  12.03,  12.37,  12.56,  12.78,  14.32,
     &		  15.00,  15.50,  15.89,  16.22,  16.30,  16.50,  16.75,
     &            17.35,  18.14,  18.67,  19.08,  19.45,  19.77,  20.06,
     &            20.34,  20.71,  20.98,  21.22,  21.44,  21.65,  21.73,
     &            21.83,  22.09,  22.34,  22.50,  22.77,  23.02,  23.25,
     &            23.48,  23.71,  23.91,  24.10,  24.30,  24.43,  24.50,  
     &		  25.60,  24.75,  24.88,  24.99,  25.08,  25.15,  25.21,     
     &		  25.25,  25.28,  25.29,  25.29,  25.28,  25.25,  25.22,
     &		  25.18,  25.12,  25.06,  25.00,  24.93,  24.90,  24.85,     
     &		  24.77,  24.70,  24.60,  24.50,  24.40,  24.30,  24.20,
     &		  24.05,  24.00,  23.95,  23.80,  23.70,  23.60,  23.50,     
     &		  23.40,  23.25,  23.13,  23.00,  22.80,  22.60,  22.40,
     &		  22.20,  22.00,  21.70,  21.40,  21.10,  20.80,  20.30,     
     &		  20.00,  19.70,  19.50,  19.30,  19.10,  18.80,  18.70,
     &		  18.60,  18.50,  18.40,  18.20,  18.00,  17.50,  17.00,
     &		  16.50,  16.29,  16.27,  16.25,  16.23,  16.20,  16.18,
     &		  16.16,  16.14,  16.12,  16.10,  16.07,  16.04,  16.01,
     &		  15.98,  15.95,  15.92,  15.89,  15.86,  15.83,  15.80,
     &		  15.77,  15.74,  15.71,  15.68,  15.65,  15.62,  15.59,
     &		  15.56,  15.53,  15.50,  15.47,  15.44,  15.41,  15.38, 
     &		  15.35,  15.32,  15.29,  15.26,  15.23,  15.20,  15.17, 
     &		  15.14,  15.11,  15.08,  15.05,  15.02,  14.99,  14.96,
     &		  14.93,  14.90,  14.87,  14.84,  14.81,  14.78,  14.75,
     &		  14.72,  14.69,  14.66,  14.63,  14.60,  14.57,  14.54,
     &		  14.51,  14.48,  14.45,  14.42,  14.39,  14.36,  14.33,
     &		  14.30,  14.27,  14.24,  14.22,  14.19,  14.16,  14.13,
     &		  14.10/

	DATA ((SIGR(I,J),J=1,251),I=12,12) /0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.03,
     &		   0.07,   0.09,   0.11,   0.15,   0.18,   0.24,   0.31,
     &		   0.38,   0.49,   0.55,   0.67,   0.73,   0.89,   1.00,
     &		   1.32,   1.80,   2.20,   2.70,   3.40,   3.79,   4.00,
     &		   5.10,   6.50,   8.30,  10.73,  11.50,  11.90,  12.50,
     &		  12.90,  13.42,  13.00,  12.40,  11.30,  10.20,   9.13,
     &		   9.05,   9.00,   8.89,   8.65,   8.56,   8.60,   8.90,
     &		   9.30,   9.50,  10.00,   9.98,   9.92,   9.87,   9.73,
     &		   9.65,   9.46,   9.32,   9.11,   8.85,   8.67,   8.36,
     &		   8.28,   8.16,   8.07,   7.94,   7.85,   7.81,   7.75,
     &		   7.64,   7.50,   7.45,   7.38,   7.33,   7.30,   7.25,
     &		   7.21,   7.19,   7.17,   7.15,   7.13,   7.12,   7.11,
     &		   7.10,   7.09,   7.09,   7.09,   7.09,   7.09,   7.10,
     &		   7.09,   7.10,   7.10,   7.11,   7.12,   7.13,   7.14,
     &		   7.14,   7.15,   7.16,   7.17,   7.18,   7.19,   7.20,
     &		   7.21,   7.22,   7.22,   7.23,   7.24,   7.25,   7.26,
     &		   7.26,   7.27,   7.27,   7.28,   7.29,   7.29,   7.30,
     &		   7.30,   7.30,   7.31,   7.31,   7.31,   7.31,   7.32,
     &		   7.32,   7.32,   7.32,   7.32,   7.32,   7.32,   7.31,
     &		   7.31,   7.30,   7.30,   7.30,   7.29,   7.28,   7.28,
     &		   7.28,   7.27,   7.27,   7.26,   7.25,   7.24,   7.23,
     &		   7.22,   7.21,   7.20,   7.19,   7.18,   7.17,   7.15,
     &		   7.14,   7.13,   7.12,   7.10,   7.09,   7.07,   7.06,
     &		   7.05,   7.04,   7.03,   7.02,   7.00,   6.97,   6.95,
     &		   6.93,   6.91,   6.90,   6.89,   6.87,   6.87,   6.85,
     &		   6.83,   6.80,   6.78,   6.76,   6.74,   6.72,   6.70,
     &		   6.68,   6.66,   6.64,   6.62,   6.60,   6.58,   6.56,
     &		   6.54,   6.52,   6.50,   6.48,   6.46,   6.44,   6.42,
     &		   6.40,   6.38,   6.36,   6.34,   6.32,   6.29,   6.27,
     &		   6.25,   6.23,   6.20,   6.18,   6.15,   6.13,   6.11,
     &		   6.09,   6.07,   6.05,   6.03,   6.01,   6.00,   5.96,
     &		   5.93,   5.91,   5.89,   5.87,   5.85,   5.83,   5.81,
     &		   5.79,   5.75,   5.73,   5.70,   5.67,   5.65,   5.64,
     &		   5.62/
	DATA ((SIGR(I,J),J=1,251),I=13,13) /251*1./
	DATA ((SIGR(I,J),J=1,251),I=14,14) /251*1./
	DATA ((SIGR(I,J),J=1,251),I=15,15) /0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.05,   0.10,   0.20,   0.30,   0.80,
     &		   1.00,   1.50,   1.80,   2.30,   2.80,   3.10,   3.50,
     &		   4.10,   4.90,   5.40,   6.20,   7.30,   8.30,   8.90,
     &		   9.40,  10.50,  11.40,  12.50,  13.10,  13.60,  14.10,
     &		  14.70,  15.00,  15.50,  16.10,  16.80,  17.00,  17.50,
     &		  18.50,  19.50,  20.00,  20.70,  21.30,  21.90,  22.20,
     &		  22.80,  22.90,  23.00,  23.20,  23.40,  23.40,  23.30,
     &		  23.20,  23.10,  22.90,  22.80,  22.50,  22.00,  21.20,
     &		  20.00,  18.90,  17.70,  16.70,  15.60,  14.90,  14.30,
     &		  14.00,  13.20,  12.80,  12.50,  12.20,  12.00,  11.80,
     &		  11.40,  11.20,  10.80,  10.30,  10.00,   9.80,   9.60,
     &		   9.50,   9.40,   9.30,   9.20,   9.10,   9.00,   8.90,
     &		   8.80,   8.70,   8.50,   8.40,   8.30,   8.20,   8.00,
     &		   7.90,   7.70,   7.50,   7.30,   7.20,   7.15,   7.10,
     &		   7.03,   7.00,   6.96,   6.92,   6.89,   6.85,   6.83,
     &		   6.80,   6.72,   6.70,   6.65,   6.62,   6.59,   6.57,
     &		   6.55,   6.53,   6.51,   6.48,   6.47,   6.45,   6.39,
     &		   6.32,   6.25,   6.11,   6.06,   5.79,   5.64,   5.71,
     &		   5.46,   5.47,   5.33,   5.38,   5.44,   5.38,   5.40,
     &		   5.32,   5.33,   5.29,   5.25,   5.25,   5.25,   5.26,
     &		   5.19,   5.17,   5.17,   5.20,   5.09,   5.12,   5.07,
     &		   4.95,   4.77,   4.65,   4.55,   4.52,   4.43,   4.40,
     &		   4.38,   4.28,   4.29,   4.25,   4.36,   4.24,   4.01,
     &		   3.79,   3.65,   3.55,   3.43,   3.44,   3.44,   3.32,
     &		   3.18,   3.28,   3.12,   3.11,   3.25,   3.04,   2.94,
     &		   2.78,   2.73,   2.67,   2.55,   2.43,   2.33,   2.28,
     &		   2.20,   2.30,   2.22,   2.21,   2.23,   2.22,   2.22,
     &		   2.19,   2.17,   2.19,   2.18,   2.17,   2.17,   2.18,
     &		   2.16,   2.15,   2.15,   2.16,   2.14,   2.13,   2.11,
     &		   2.11,   2.10,   2.10,   2.10,   2.10,   2.11,   2.12,
     &		   2.10,   2.08,   2.07,   2.09,   2.08,   2.08,   2.09,
     &		   2.07,   2.06,   2.05,   2.04,   2.03,   2.02,   2.01,
     &		   2.01/
	DATA ((SIGR(I,J),J=1,251),I=16,16) /251*1./
	DATA ((SIGR(I,J),J=1,251),I=17,17) /251*1./
	DATA ((SIGR(I,J),J=1,251),I=18,18) /251*1./
	DATA ((SIGR(I,J),J=1,251),I=19,19) /0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.30,   0.70,   1.00,
     &		   1.50,   1.80,   2.00,   2.40,   2.76,   3.02,   3.50,
     &		   4.10,   4.60,   4.87,   5.20,   5.65,   6.03,   6.54,
     &		   6.70,   6.80,   6.92,   7.00,   7.15,   7.22,   7.32,
     &		   7.40,   7.50,   7.52,   7.50,   7.48,   7.48,   7.47,
     &		   7.47,   7.46,   7.45,   7.44,   7.42,   7.40,   7.41,
     &		   7.42,   7.41,   7.40,   7.40,   7.37,   7.35,   7.33,
     &		   7.30,   7.31,   7.30,   7.30,   7.30,   7.30,   7.25,
     &		   7.21,   7.18,   7.02,   6.99,   6.97,   6.94,   6.90,
     &		   6.85,   6.82,   6.81,   6.76,   6.72,   6.68,   6.60,
     &		   6.50,   6.40,   6.38,   6.33,   6.28,   6.23,   6.12,
     &		   5.90,   5.80,   5.70,   5.60,   5.50,   5.42,   5.38,
     &		   5.30,   5.25,   5.20,   5.18,   5.13,   5.08,   5.03,
     &		   4.98,   4.92,   4.87,   4.80,   4.78,   4.73,   4.68,
     &		   4.63,   4.58,   4.54,   4.50,   4.45,   4.40,   4.35,
     &		   4.30,   4.28,   4.23,   4.20,   4.18,   4.15,   4.13,
     &		   4.08,   4.05,   4.03,   4.00,   3.98,   3.95,   3.93,
     &		   3.89,   3.86,   3.83,   3.80,   3.75,   3.72,   3.70,
     &		   3.66,   3.63,   3.60,   3.58,   3.56,   3.54,   3.53,
     &		   3.52,   3.51,   3.50,   3.48,   3.46,   3.44,   3.42,
     &		   3.41,   3.40,   3.41,   3.40,   3.40,   3.40,   3.40,
     &		   3.40,   3.38,   3.36,   3.34,   3.32,   3.31,   3.30,
     &		   3.28,   3.26,   3.24,   3.23,   3.21,   3.20,   3.19,
     &		   3.18,   3.17,   3.16,   3.15,   3.14,   3.13,   3.12,
     &		   3.11,   3.10,   3.09,   3.08,   3.07,   3.06,   3.05,
     &		   3.04,   3.03,   3.02,   3.01,   3.00,   2.99,   2.98,
     &		   2.97,   2.96,   2.95,   2.94,   2.93,   2.92,   2.91,
     &		   2.90/
	DATA ((SIGR(I,J),J=1,251),I=20,20) /251*1./
	DATA ((SIGR(I,J),J=1,251),I=21,21) /0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.04,   0.09,   0.12,
     &		   0.18,   0.42,   0.54,   0.72,   0.90,   1.08,   1.20,
     &		   1.50,   1.80,   2.10,   2.40,   2.58,   2.82,   3.00,
     &		   3.30,   3.60,   3.78,   4.02,   4.20,   4.32,   4.50,
     &		   4.74,   4.92,   5.10,   5.28,   5.52,   5.70,   6.00,
     &		   6.30,   6.60,   6.78,   7.02,   7.20,   7.44,   7.80,
     &		   8.10,   8.46,   8.94,   9.60,  10.02,  10.31,  10.63,
     &		  10.95,  11.12,  11.43,  11.75,  12.00,  12.18,  12.42,
     &		  12.66,  12.96,  13.20,  13.38,  13.62,  13.86,  14.10,
     &		  14.28,  14.88,  15.30,  15.66,  15.90,  16.32,  16.56,
     &		  16.80,  16.56,  16.86,  17.10,  17.76,  17.88,  18.00,
     &		  18.24,  18.30,  18.42,  18.54,  18.60,  18.78,  18.90,
     &		  19.02,  19.74,  19.74,  19.79,  19.80,  19.80,  19.80,
     &		  19.86,  19.86,  19.86,  19.86,  19.86,  19.80,  19.74,
     &		  19.56,  19.44,  19.32,  19.20,  19.20,  19.20,  19.20,
     &		  19.20,  19.20,  19.05,  18.78,  18.66,  18.60,  18.60,
     &		  18.30,  18.00,  17.88,  17.82,  17.70,  17.64,  17.58,
     &		  17.52,  17.46,  17.40,  17.40,  17.34,  17.28,  17.22,
     &		  17.22,  17.19,  17.16,  17.10,  16.92,  16.80,  16.62,
     &		  16.38,  16.20,  16.02,  15.90,  15.87,  15.84,  15.82,
     &		  15.80,  15.79,  15.78,  15.77,  15.76,  15.75,  15.73,
     &		  15.73,  15.72,  15.72,  15.70,  15.69,  15.67,  15.64,
     &		  15.61,  15.60,  15.58,  15.55,  15.54,  15.51,  15.48,
     &		  15.45,  15.42,  15.39,  15.36,  15.33,  15.30,  15.27,
     &		  15.24,  15.21,  15.18,  15.15,  15.12,  15.09,  15.06,
     &		  15.03,  15.00,  14.97,  14.94,  14.91,  14.88,  14.85,
     &		  14.82,  14.79,  14.76,  14.73,  14.70,  14.67,  14.64,
     &		  14.61,  14.58,  14.55,  14.52,  14.49,  14.46,  14.43,
     &		  14.40,  14.37,  14.34,  14.31,  14.28,  14.25,  14.22,
     &		  14.19/
	DATA ((SIGR(I,J),J=1,251),I=22,22) /0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.50,   1.00,   1.50,   2.00,   2.50,   3.00,
     &		   4.00,   5.00,   6.00,   8.00,   8.50,   9.00,   9.70,
     &		  10.20,  10.80,  11.50,  12.20,  13.00,  13.60,  14.30,
     &		  14.80,  15.20,  15.70,  16.00,  16.50,  16.80,  16.90,
     &		  17.00,  17.30,  17.70,  18.00,  18.30,  18.70,  19.20,
     &		  19.50,  20.00,  20.40,  20.50,  20.60,  20.70,  20.80,
     &		  20.90,  21.00,  21.50,  21.70,  22.00,  22.10,  22.10,
     &		  22.10,  22.20,  22.20,  22.20,  22.10,  22.10,  22.10, 
     &		  22.10,  22.10,  22.00,  22.00,  22.00,  21.98,  21.98,
     &		  21.97,  21.96,  21.95,  21.95,  21.94,  21.94,  21.93,
     &		  21.92,  21.91,  21.90,  21.89,  21.89,  21.89,  21.88,
     &		  21.87,  21.87,  21.86,  21.85,  21.84,  21.84,  21.83,
     &		  21.81,  21.79,  21.78,  21.76,  21.74,  21.72,  21.70,
     &		  21.65,  21.60,  21.55,  21.50,  21.40,  21.30,  21.20,
     &		  21.10,  21.00,  20.80,  20.70,  20.60,  20.50,  20.40,
     &		  20.20,  20.10,  20.00,  19.90,  19.85,  19.80,  19.60,
     &		  19.50,  19.30,  19.10,  19.00,  18.80,  18.70,  18.50,
     &		  18.40,  18.30,  18.20,  18.10,  18.00,  17.90,  17.80,
     &		  17.70,  17.60,  17.50,  17.40,  17.30,  17.20,  17.00,
     &		  16.80,  16.60,  16.40,  16.20,  16.00,  15.80,  15.60,
     &		  15.40,  15.20,  15.00,  14.80,  14.60,  14.40,  14.20,
     &		  14.00,  13.80,  13.60,  13.40,  13.20,  13.00,  12.80,
     &		  12.60,  12.40,  12.20,  12.00,  11.80,  11.60,  11.40,
     &		  11.20,  11.00,  10.80,  10.60,  10.40,  10.20,  10.00,
     &		   9.80,   9.60,   9.40,   9.20,   9.00,   8.80,   8.60,
     &		   8.50,   8.40,   8.30,   8.20,   8.10,   8.00,   7.90,
     &		   7.80,   7.70,   7.60,   7.50,   7.40,   7.30,   7.20,
     &		   7.10,   7.00,   6.90,   6.80,   6.70,   6.60,   6.50,
     &		   6.40/
	DATA ((SIGR(I,J),J=1,251),I=23,23) /251*1./
	DATA ((SIGR(I,J),J=1,251),I=24,24) /251*1./
        data (sigr(25,j),j=1,251) /
     @   0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     @   0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     @   0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     @   0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     @   0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     @   0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     @   0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     @   0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     @   0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     @   0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     @   0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     @   0.00,  0.00,  0.00,  0.00,  0.00,  4.12,  4.22,  4.32,
     @   4.42,  4.52,  4.61,  4.71,  4.81,  4.91,  5.01,  5.11,
     @   5.20,  5.30,  5.40,  5.49,  5.59,  5.69,  5.78,  5.87,
     @   5.97,  6.06,  6.16,  6.25,  6.34,  6.44,  6.53,  6.62,
     @   6.71,  6.80,  6.89,  6.98,  7.07,  7.16,  7.25,  7.34,
     @   7.42,  7.51,  7.60,  7.68,  7.77,  7.86,  7.94,  8.03,
     @   8.11,  8.19,  8.28,  8.36,  8.44,  8.52,  8.60,  8.69,
     @   8.77,  8.85,  8.92,  9.00,  9.08,  9.16,  9.24,  9.31,
     @   9.39,  9.47,  9.54,  9.62,  9.69,  9.76,  9.84,  9.91,
     @   9.98, 10.06, 10.13, 10.20, 10.27, 10.34, 10.41, 10.48,
     @  10.54, 10.61, 10.68, 10.75, 10.81, 10.88, 10.94, 11.01,
     @  11.07, 11.14, 11.20, 11.26, 11.33, 11.39, 11.45, 11.51,
     @  11.57, 11.63, 11.69, 11.75, 11.81, 11.86, 11.92, 11.98,
     @  12.03, 12.09, 12.15, 12.20, 12.25, 12.31, 12.36, 12.42,
     @  12.47, 12.52, 12.57, 12.62, 12.67, 12.72, 12.77, 12.82,
     @  12.87, 12.92, 12.97, 13.01, 13.06, 13.11, 13.15, 13.20,
     @  13.24, 13.29, 13.33, 13.38, 13.42, 13.46, 13.50, 13.55,
     @  13.59, 13.63, 13.67, 13.71, 13.75, 13.79, 13.83, 13.87,
     @  13.91, 13.94, 13.98, 14.02, 14.05, 14.09, 14.13, 14.16,
     @  14.20, 14.23, 14.26, 14.30, 14.33, 14.37, 14.40, 14.43,
     @  14.46, 14.49, 14.52 /
        data (sigr(26,j),j=1,251) /251*1/

	DATA ((SIGR(I,J),J=1,251),I=27,27) /0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.01,   0.01,
     &		   0.01,   0.01,   0.01,   0.02,   0.03,   0.03,   0.03,
     &		   0.03,   0.04,   0.04,   0.04,   0.04,   0.04,   0.05,
     &		   0.05,   0.05,   0.05,   0.05,   0.05,   0.05,   0.05,
     &		   0.06,   0.06,   0.06,   0.06,   0.06,   0.06,   0.06,
     &		   0.06,   0.07,   0.07,   0.07,   0.07,   0.07,   0.08,
     &		   0.08,   0.09,   0.10,   0.11,   0.12,   0.13,   0.14,
     &		   0.16,   0.18,   0.20,   0.23,   0.26,   0.29,   0.30,
     &		   0.35,   0.40,   0.45,   0.50,   0.55,   0.60,   0.70,
     &		   0.72,   0.75,   0.80,   0.85,   0.90,   0.93,   0.95,
     &		   0.97,   0.99,   1.00,   1.05,   1.07,   1.07,   1.09,
     &	           1.09,   1.08,   1.09,   1.10,   1.11,   1.12,   1.14,
     &		   1.15,   1.16,   1.17,   1.18,   1.19,   1.20,   1.22,
     &		   1.28,   1.31,   1.35,   1.41,   1.55,   1.62,   1.70,
     &		   1.72,   1.75,   1.78,   1.82,   1.90,   1.95,   1.98,
     &		   2.00,   2.06,   2.13,   2.19,   2.26,   2.33,   2.40,
     &		   2.46,   2.48,   2.55,   2.63,   2.71,   2.79,   2.85,
     &		   2.92,   2.95,   3.06,   3.09,   3.14,   3.21,   3.27,
     &		   3.31,   3.36,   3.41,   3.43,   3.50,   3.55,   3.62,
     &		   3.74,   3.81,   3.85,   3.86,   3.92,   3.93,   3.98,
     &		   4.02,   4.08,   4.11,   4.16,   4.20,   4.25,   4.27,
     &		   4.33,   4.37,   4.44,   4.48,   4.52,   4.57,   4.61,
     &		   4.66,   4.66,   4.69,   4.72,   4.77,   4.82,   4.85,
     &		   4.90,   4.93,   4.96,   4.99,   5.00,   5.03,   5.10,
     &		   5.17,   5.19,   5.21,   5.25,   5.31,   5.37,   5.41,
     &		   5.45,   5.51,   5.53,   5.57,   5.58,   5.60,   5.61,
     &		   5.62,   5.63,   5.65,   5.68,   5.69,   5.71,   5.70,
     &		   5.68,   5.67,   5.71,   5.72,   5.73,   5.73,   5.71,
     &		   5.72,   5.71,   5.70,   5.69,   5.67,   5.67,   5.66,
     &		   5.64,   5.63,   5.63,   5.64,   5.65,   5.65,   5.67,
     &		   5.64,   5.65,   5.65,   5.66,   5.64,   5.64,   5.65,
     &		   5.63,   5.62,   5.61,   5.61,   5.61,   5.62,   5.60,
     &		   5.60/
	DATA ((SIGR(I,J),J=1,251),I=28,28) /0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.01,   0.02,   0.02,   0.03,   0.03,   0.04,   0.05,
     &		   0.05,   0.06,   0.07,   0.08,   0.09,   0.09,   0.10,
     &		   0.11,   0.12,   0.13,   0.15,   0.17,   0.20,   0.22,
     &		   0.23,   0.25,   0.26,   0.27,   0.29,   0.30,   0.32,
     &		   0.33,   0.35,   0.37,   0.40,   0.42,   0.45,   0.48,
     &		   0.52,   0.53,   0.55,   0.57,   0.59,   0.62,   0.65,
     &	           0.67,   0.69,   0.73,   0.75,   0.78,   0.81,   0.82,
     &		   0.85,   0.88,   0.91,   0.94,   0.97,   0.98,   1.02,
     &		   1.11,   1.17,   1.22,   1.28,   1.34,   1.48,   1.55,
     &		   1.67,   1.72,   1.78,   1.82,   1.84,   1.93,   2.05,
     &		   2.18,   2.26,   2.37,   2.55,   2.65,   2.71,   2.74,
     &		   2.78,   2.83,   2.90,   2.95,   3.01,   3.25,   3.33,
     &		   3.41,   3.57,   3.81,   3.92,   4.35,   4.67,   4.88,
     &		   4.98,   5.00,   5.05,   5.12,   5.35,   5.46,   5.67,
     &		   5.74,   5.91,   6.03,   6.15,   6.22,   6.38,   6.45,
     &		   6.51,   6.58,   6.64,   6.69,   6.73,   6.84,   6.94,
     &		   6.97,   7.26,   7.58,   7.89,   8.03,   8.24,   8.33,
     &		   8.45,   8.61,   8.74,   8.91,   9.05,   9.55,   9.89,
     &		  10.35,  10.77,  11.02,  11.23,  11.18,  11.15,  11.13,
     &		  11.10,  11.08,  11.08,  11.09,  11.03,  11.00,  11.00,
     &		  10.99,  10.99,  10.97,  10.98,  10.96,  10.95,  10.84,
     &		  10.95,  10.95,  10.94,  10.93,  10.93,  10.92,  10.92,
     &		  10.91,  10.91,  10.92,  10.90,  10.88,  10.87,  10.86,
     &		  10.85,  10.84,  10.85,  10.86,  10.84,  10.83,  10.83,
     &		  10.82,  10.83,  10.85,  10.83,  10.81,  10.81, 10.80,
     &		  10.81,  10.80,  10.78,  10.79,  10.80,  10.77, 10.78,
     &		  10.76,  10.77,  10.77,  10.77,  10.76,  10.74, 10.75,
     &		  10.78/

       DATA((SIGR(I,J),J=1,251),I=29,29) /0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
     &		   0.00,   0.01,   0.03,   0.03,   0.02,   0.05,   0.05,
     &		   0.06,   0.07,   0.07,   0.07,   0.09,   0.13,   0.13,
     &		   0.15,   0.16,   0.16,   0.18,   0.16,   0.18,   0.16,
     &		   0.19,   0.20,   0.21,   0.22,   0.21,   0.22,   0.23,
     &		   0.23,   0.25,   0.23,   0.22,   0.23,   0.26,   0.26,
     &		   0.27,   0.28,   0.30,   0.31,   0.32,   0.34,   0.35,
     &		   0.32,   0.34,   0.37,   0.39,   0.38,   0.41,   0.40,
     &		   0.40,   0.38,   0.42,   0.45,   0.47,   0.49,   0.50,
     &		   0.51,   0.53,   0.54,   0.56,   0.58,   0.59,   0.61,
     &		   0.63,   0.64,   0.67,   0.69,   0.71,   0.72,   0.74,
     &		   0.73,   0.72,   0.75,   0.73,   0.72,   0.75,   0.77,
     &	           0.79,   0.80,   0.82,   0.84,   0.85,   0.87,   0.89,
     &		   0.91,   0.93,   0.94,   0.95,   0.96,   0.97,   0.99,
     &		   1.00,   1.03,   1.07,   1.13,   1.17,   1.23,   1.36,
     &		   1.39,   1.41,   1.47,   1.56,   1.62,   1.71,   1.75,
     &		   1.81,   1.83,   1.85,   1.90,   1.87,   1.89,   1.92,
     &		   1.94,   1.98,   2.00,   2.03,   2.07,   2.14,   2.16,
     &		   2.15,   2.19,   2.23,   2.27,   2.28,   2.29,   2.30,
     &		   2.31,   2.33,   2.34,   2.35,   2.43,   2.48,   2.53,
     &		   2.59,   2.63,   2.71,   2.79,   2.83,   2.88,   2.94,
     &		   2.98,   3.02,   3.17,   3.25,   3.34,   3.55,   3.59,
     &		   3.63,   3.65,   3.74,   3.83,   3.89,   3.92,   4.03,
     &		   4.07,   4.11,   4.15,   4.17,   4.31,   4.39,   4.47,
     &		   4.58,   4.70,   4.85,   4.91,   4.97,   5.01,   5.21,
     &		   5.22,   5.29,   5.34,   5.39,   5.47,   5.58,   5.64,
     &		   5.68,   5.78,   5.82,   5.88,   5.95,   6.03,   6.11,
     &		   6.18,   6.25,   6.28,   6.49,   6.59,   6.74,   6.89,
     &		   7.05,   7.23,   7.49,   7.58,   7.69,   7.84,   7.91,
     &		   7.97,   7.93,   7.88,   7.87,   7.87,   7.91,   7.90,
     &		   7.88,   7.87,   7.86,   7.85,   7.86,   7.86,   7.84,
     &		   7.83,   7.81,   7.88,   7.86,   7.85,   7.84,   7.83,
     &		   7.83,   7.82,   7.82,   7.81,   7.81,   7.80,   7.81,
     &		   7.81/
        data (sigr(30,j),j=1,19) /19*0.0/
        data (sigr(30,j),j=20,251) /232*1.5/

        data (sigr(31,j),j=1,19) /19*0.0/
        data (sigr(31,j),j=20,251) /232*1./
        data iflagphi/0/


c*** riaggiusta la sezione d'urto della fi
        if (iflagphi.eq.0) then
        do j = 1,251
         sigr(25,j) = sigr(25,j)/30.
        enddo
        iflagphi=1
        endif


	 DO K=1,251
		SIGR(7,K)  = SIGR(5,K)*1/3
		SIGR(10,K) = SIGR(8,K)*2
		SIGR(13,K) = SIGR(12,K)
		SIGR(14,K) = SIGR(12,K)*2/3
		SIGR(16,K) = SIGR(15,K)/2
		SIGR(17,K) = SIGR(15,K)
		SIGR(18,K) = SIGR(15,K)/2
		SIGR(20,K) = SIGR(19,K)
		SIGR(23,K) = SIGR(21,K)
		SIGR(24,K) = SIGR(22,K)*5/3
		SIGR(26,K) = SIGR(25,K)
	ENDDO

	WR(1)=1.0
	DO J=2,251
	 WR(J)=WR(J-1)+0.01
	ENDDO


	RETURN
	END

C	*************************
	SUBROUTINE CANALE(W,JCH,SIG)
C	*************************
	IMPLICIT NONE
	REAL SIGMA(24),W,Q,TOT,SIG
	REAL PROB(0:24),SIGR(40,251),WR(251)
	INTEGER ICHAN(24),JCH,NCH,J,K
	INTEGER*4 IZ,JJ
	REAL RRAN
        COMMON/SEZURTO/SIGR,WR
        COMMON/CHANNELS/NCH,ICHAN

	COMMON/RANDOM/IZ

C	NCH = # OF CHANNELS ACCOUNTED
C	ICHAN(I) = TYPE OF ACCOUNTED CHANNELS
C	E.G. :
C	NCH = 4  AND  ICHAN=2,6,7,11  --> IT MEANS:
C	ONLY CHANNELS # 2,6,7 AND 11 ARE TAKEN INTO ACCOUNT
C
C
C	INITIALIZING VARIABLES..
	TOT = 0.
	PROB(0) = 0.

C
C	EVALUATING TOTAL CROSS SECTIOS AT W VALUE
C	WITH AN INTERPOLATION PROCEDURE......
	DO J=1,250
	 IF(W.GT.WR(J).AND.W.LE.WR(J+1)) then
          jj=j
	  GOTO 1
	 endif
	ENDDO
1	CONTINUE
	DO K=1,NCH
	 SIGMA(K) = SIGR(ICHAN(K),JJ) +
     &	 (W-WR(JJ))/(WR(JJ+1)-WR(JJ))*
     &	 (SIGR(ICHAN(K),JJ+1)-SIGR(ICHAN(K),JJ))
	ENDDO

C	EVALUATING TOTAL CROSS SECTION:
	DO K = 1, NCH
		TOT = TOT + SIGMA(K)
	ENDDO
        IF(TOT.NE.0) THEN
C
C	EVALUATING THE RELATIVE "WEIGHT" OF THE CHANNELS.......
	DO J = 1,NCH
		PROB(J) = PROB(J-1) + SIGMA(J)
	ENDDO
C
C	NORMALIZING......
	DO J = 1,NCH
		PROB(J) = PROB(J)/TOT
	ENDDO
C
C	EXTRACTING RANDOM NUMBER.....
	Q = rran()
C
C	CHOOSING THE CHANNEL....
	DO J=1,NCH
         IF (Q.GT.PROB(J-1).AND.Q.LE.PROB(J)) then
          JJ=J
          GO TO 2
         endif
	ENDDO

2	JCH = ICHAN(JJ)
	SIG = TOT 
        ELSE
        JCH = 0
        ENDIF


	RETURN
	END

C**************************************
        SUBROUTINE DECAD(JCH,NP,AMASS)
C*************************************

	 IMPLICIT NONE
         INTEGER NP,JCH
         REAL AMASS(18)

         IF(JCH.EQ.5)THEN
         NP=2
         AMASS(1)=0.93827231
         AMASS(2)=0.1395679
         ENDIF
         IF(JCH.EQ.6.OR.JCH.EQ.8)then
         NP=2
         AMASS(1)=0.93827231
         AMASS(2)=0.1349743
         ENDIF
         IF(JCH.EQ.7.OR.JCH.EQ.9)then
         NP=2
         AMASS(1)=0.93827231
         AMASS(2)=0.1395679
         ENDIF
         IF(JCH.EQ.10)then
         NP=2
         AMASS(1)=0.93956563
         AMASS(2)=0.1395679
         ENDIF

         IF(JCH.GE.11.AND.JCH.LE.14) THEN
         NP = 2
         AMASS(1)=0.1395679
         AMASS(2)=0.1395679
         ENDIF
         IF(JCH.GE.19.AND.JCH.LE.20) THEN
         NP = 3
         AMASS(1)=0.1395679
         AMASS(2)=0.1395679
         AMASS(3)=0.1349743
         ENDIF
         RETURN
         END

C**********************************
         SUBROUTINE ANG(JCH,W,THP)
C**********************************
	 IMPLICIT NONE
         REAL XSECT(14,50,0:18)
         REAL RRAN
         INTEGER JCHI,JCH,I,J,II
	 INTEGER*4 IZ	 
         REAL WI,WF,DW,THR,SIMAX,W,THP,SR,WR
	 REAL CONV,RJ,RJ1
	 REAL DECOS,ANGM

         COMMON/SIGMAS/XSECT
 	COMMON/RANDOM/IZ
C+ First 14 channels have the right angular cross section
         IF(JCH.LE.14)           JCHI=JCH
c++ Omega and PHI have the  same angular distribution of rho
         IF(JCH.EQ.19.OR.JCH.EQ.20.or.JCH.EQ.25.OR.JCH.EQ.26) JCHI=11
         WI=1.
         WF=3.5
         DW=0.05
         DO I=1,50
          IF(W.GE.(WI+I*DW).AND.W.LT.(WI+(I+1)*DW)) then
	   II=I
           GOTO 10
	  endif
         ENDDO
10       WR=WI+II*DW
C------------------------------------------------------------
C        CALCOLO SIMAX
C------------------------------------------------------------

         SIMAX=0 
            CONV = 3.1415/19.
         DO J=0,18
	    RJ=FLOAT(J)
	    RJ1=(RJ+1.)			
	    DECOS=COS((CONV*RJ))-COS((CONV*RJ1))
	    DECOS=ABS(DECOS)
	    ANGM=RJ*CONV+CONV/2.
C	    XSECT(JCHI,II,J)=1. 
            SIMAX=MAX(SIMAX,DECOS*XSECT(JCHI,II,J))
         ENDDO
20           THR=rran()*(0.-180.)+180.
	     J=INT((THR/180.*19.))
	     if(J.gt.18) J=18
            SR=rran()*SIMAX
	    RJ=FLOAT(J)
            RJ1=(RJ+1.)                 
	    DECOS=COS((CONV*RJ))-COS((CONV*RJ1))
            DECOS=ABS(DECOS)
            ANGM=RJ*CONV+CONV/2.
         IF(SR.LE.DECOS*XSECT(JCHI,II,J)) THEN
         THP=THR*3.1415/180.

         ELSE             
         GOTO 20
         ENDIF
         RETURN
         END
      
C***************************************************

	SUBROUTINE READ_W_DEC_RHO
C***************************************************
C--------------------------------------------------
C	OPENING THE PROBABILITY DISTRIBUTION FILE AND
C	NORMALIZTION OF W_DEC_RHO
C--------------------------------------------------
	
	IMPLICIT NONE
	
	REAL W_MAX,W_DEC_RHO(21,41),TH,PH,THDUMMY,PHDUMMY,
     $        W_DEC(21,41)

	INTEGER I,J
	CHARACTER *120 genev_parms,filename
	parameter (genev_parms  = 'GENOVA_PARMS')
	COMMON/RHO/W_DEC_RHO
c get full path to data file	
	CALL revinm(genev_parms,'genev_rhodec.dat',filename)
	OPEN(UNIT=62,FILE=filename,STATUS='OLD',FORM='FORMATTED')
	W_MAX = 0.
	DO I=1,21
	 DO J = 1,41
	 TH = FLOAT(I-1)*3.141592653589793/20.
	 PH = FLOAT(J-1)*3.141592653589793*2./40.
C	 W_DEC_RHO(I,J) = SIN(PH)*SIN(TH)
	 READ(62,*) THDUMMY,PHDUMMY,W_DEC(I,J)
	 W_DEC_RHO(I,J) = W_DEC(I,J)*SIN(TH)
	 W_MAX = MAX(W_MAX,W_DEC_RHO(I,J))
	 ENDDO
	ENDDO
	DO I=1,21
	 DO J = 1,41
	  W_DEC_RHO(I,J) = W_DEC_RHO(I,J)/W_MAX
	 ENDDO
	ENDDO
C	CLOSE(62)
	RETURN
	END

C*************************************************
	SUBROUTINE RHO_DECAY(WRO,FLAG,parm_Q2,parm_W,parm_GPhi,parm_EPS)
C*************************************************

	IMPLICIT NONE
	REAL   parm_Q2,parm_W,parm_GPhi,parm_EPS,
     @         rat, deno,
     @         rho00_0,rho00_4, Rerho10_0,Rerho10_4, rho1m1_0,rho1m1_4,
     @         rho11_1,rho00_1,Rerho10_1,rho1m1_1,Imrho10_2,Imrho1m1_2,
     @         rho11_5,rho00_5,Rerho10_5,rho1m1_5,Imrho10_6,Imrho1m1_6, 
     @         R00_04,ReR10_04,R1m1_04,
     @         R11_1,R00_1,ReR10_1,R1m1_1,ImR10_2,ImR1m1_2,
     @         R11_5,R00_5,ReR10_5,R1m1_5,ImR10_6,ImR1m1_6,
     @         decay_distr_max,decay_distr,pirad,
     @         thpi, phpi,
     @         RRAN,PCM(5,18),WGT,WRO,WRD,MP,MPI
	INTEGER FLAG,jthpi,jphpi 
	DATA MP,MPI/0.93827200,0.1395702/				
	COMMON/GENOUT/PCM,WGT 	
	
	pirad=3.14159


c define spin density matrix elements
c	rat = 0.75*(parm_Q2)**(1.09)
c	deno = 1 + parm_EPS*rat

c++ Definition from rhoXY_Z = spin density matrix elements in photoproduction
c	R00_04 =  (rho00_0 + parm_EPS*rat*rho00_4)  / deno
c	ReR10_04 =  (Rerho10_0 + parm_EPS*rat*Rerho10_4)  / deno
c	R1m1_04 = (rho1m1_0 + parm_EPS*rat*rho1m1_4)/ deno

c	R11_1 = rho11_1 / deno
c	R00_1 = rho00_1 / deno
c	ReR10_1 = Rerho10_1 / deno
c	R1m1_1 = rho1m1_1 / deno
c	ImR10_2 = Imrho10_2 / deno
c	ImR1m1_2 = Imrho1m1_2 / deno

c       R11_5 = sqrt(rat)*rho11_5 / deno
c	R00_5 = sqrt(rat)*rho00_5 / deno
c	ReR10_5 = sqrt(rat)*Rerho10_5 / deno
c	R1m1_5 = sqrt(rat)*rho1m1_5 / deno
c	ImR10_6 = sqrt(rat)*Imrho10_6 / deno
c	ImR1m1_6 = sqrt(rat)*Imrho1m1_6 / deno

c++ Values from Joos et al, NP B113 (1976) for Q2 = 1.05 and W = 2-2.8
	R00_04 =  0.47
	ReR10_04 = 0.16
        R1m1_04 = -0.06
	R11_1 = -0.16
	R00_1 = 0.26
	ReR10_1 = 0.07
	R1m1_1 = 0.0
	ImR10_2 = -0.01
	ImR1m1_2 = -0.17
	R11_5 = 0.00
	R00_5 = -0.04
	ReR10_5 = -0.06
	R1m1_5 = -0.04
	ImR10_6 = -0.08
	ImR1m1_6 = 0.09

c find maximum of decay distribution (Nota : neglect epsilon and gphi dependence here !!!!)
c cf equation 92 of Schilling and Wolf, NP B61 (1973)
	decay_distr_max=-1.
	do jthpi=0,30
	   do jphpi=0,90
	     thpi=jthpi*3.1416/30.
	     phpi=jphpi*6.2832/90.
	     decay_distr=3./4./pirad*(
     @                         0.5 * (1-R00_04)
     @                        +0.5 * (3*R00_04-1)*cos(thpi)**2
     @                        -1.414 * ReR10_04 * sin(2*thpi) * cos(phpi)
     @                        -R1m1_04 * sin(thpi)**2 * cos(2*phpi) ! end W^0
     @                     -parm_EPS*cos(2*parm_GPhi)*(
     @                         R11_1*sin(thpi)**2 
     @                        +R00_1*cos(thpi)**2
     @                        -1.414*ReR10_1*sin(2*thpi)*cos(phpi)
     @                        -R1m1_1*sin(thpi)**2*cos(2*phpi))   ! end W^1
     @                     -parm_EPS*sin(2*parm_GPhi)*(
     @                         1.414*ImR10_2*sin(2*thpi)*sin(phpi)
     @                        +ImR1m1_2*sin(thpi)**2*sin(2*phpi)) ! end W^2
     @                     +sqrt(2*parm_EPS*(1+parm_EPS))*cos(parm_GPhi)*(
     @                         R11_5*sin(thpi)**2
     @                        +R00_5*cos(thpi)**2
     @                        -1.414*ReR10_5*sin(2*thpi)*cos(phpi)
     @                        -R1m1_5*sin(thpi)**2*cos(2*phpi))   ! end W^3  
     @                     +sqrt(2*parm_EPS*(1+parm_EPS))*sin(parm_GPhi)*(
     @                         1.414*ImR10_6*sin(2*thpi)*sin(phpi)
     @                        +ImR1m1_6*sin(thpi)**2*sin(2*phpi))) ! end W^4
     @              *sin(thpi)
	     decay_distr_max=max(decay_distr_max,decay_distr)
	  enddo
	enddo
	decay_distr_max=decay_distr_max*1.5 ! safe factor (to not have fine thpi phpi binning)

c randum theta pion1 and phi pion1 in the rho rest frame  	
10	THPI = pirad*rran()
	PHPI = pirad*2.*rran()

c evaluate decay distribution
	decay_distr=3./4./pirad*(
     @                         0.5 * (1-R00_04)
     @                        +0.5 * (3*R00_04-1)*cos(thpi)**2
     @                        -1.414 * ReR10_04 * sin(2*thpi) * cos(phpi)
     @                        -R1m1_04 * sin(thpi)**2 * cos(2*phpi) ! end W^0
     @                     -parm_EPS*cos(2*parm_GPhi)*(
     @                         R11_1*sin(thpi)**2 
     @                        +R00_1*cos(thpi)**2
     @                        -1.414*ReR10_1*sin(2*thpi)*cos(phpi)
     @                        -R1m1_1*sin(thpi)**2*cos(2*phpi))   ! end W^1
     @                     -parm_EPS*sin(2*parm_GPhi)*(
     @                         1.414*ImR10_2*sin(2*thpi)*sin(phpi)
     @                        +ImR1m1_2*sin(thpi)**2*sin(2*phpi)) ! end W^2
     @                     +sqrt(2*parm_EPS*(1+parm_EPS))*cos(parm_GPhi)*(
     @                         R11_5*sin(thpi)**2
     @                        +R00_5*cos(thpi)**2
     @                        -1.414*ReR10_5*sin(2*thpi)*cos(phpi)
     @                        -R1m1_5*sin(thpi)**2*cos(2*phpi))   ! end W^3  
     @                     +sqrt(2*parm_EPS*(1+parm_EPS))*sin(parm_GPhi)*(
     @                         1.414*ImR10_6*sin(2*thpi)*sin(phpi)
     @                        +ImR1m1_6*sin(thpi)**2*sin(2*phpi))) ! end W^4
     @              *sin(thpi)

c if flag decay rho = 0 -> flat decay distrib
	IF(FLAG.EQ.0) THEN
	 decay_distr = decay_distr_max
	ENDIF

c selection of the extracted event in the decaying particle rest frame
	 WRD = rran()
	 IF((WRD*decay_distr_max).GT.(decay_distr)) GOTO 10


c	 open(unit=10,file='rho_dec.out',status='unknown')
c	 write(10,*) thpi, cos(thpi), phpi, decay_distr


c Pion1 and Pion2 variables in the rho rest frame
  	 PCM(5,1) = SQRT(0.25*WRO**2 - MPI**2)
	 PCM(5,2) = PCM(5,1)

	 PCM(1,1) = PCM(5,1)*SIN(THPI)*COS(PHPI)
	 PCM(2,1) = PCM(5,1)*SIN(THPI)*SIN(PHPI)
	 PCM(3,1) = PCM(5,1)*COS(THPI)
c	 PCM(4,1) = SQRT(PCM(1,1)**2 + PCM(2,1)**2 + PCM(3,1)**2  + MPI**2)
	 PCM(4,1) = 0.5*WRO
	 PCM(1,2) = - PCM(1,1)
	 PCM(2,2) = - PCM(2,1)
	 PCM(3,2) = - PCM(3,1)
	 PCM(4,2) =   PCM(4,1)
	
	
	RETURN
	END



C*********************************
	subroutine omega_decay(WOM,FLAG,parm_Q2,parm_W,parm_GPhi,parm_EPS)
C*********************************

	implicit none
	real   parm_Q2,parm_W,parm_GPhi,parm_EPS,
     @         rat, deno,
     @         rho00_0,rho00_4, Rerho10_0,Rerho10_4, rho1m1_0,rho1m1_4,
     @         rho11_1,rho00_1,Rerho10_1,rho1m1_1,Imrho10_2,Imrho1m1_2,
     @         rho11_5,rho00_5,Rerho10_5,rho1m1_5,Imrho10_6,Imrho1m1_6, 
     @         R00_04,ReR10_04,R1m1_04,
     @         R11_1,R00_1,ReR10_1,R1m1_1,ImR10_2,ImR1m1_2,
     @         R11_5,R00_5,ReR10_5,R1m1_5,ImR10_6,ImR1m1_6,
     @	       decay_distr_max,decay_distr,pirad,
     @         thN, phN,
     @         RRAN,PCM(5,18),WGT,WOM,WRD
	real   mpip,mpi0,mpi,m3,Qav,mu,q2,p2,A2,
     @         yx,p,q,Tz,Tm,Tp,y,x,th_dipi,phi_all,
     @         p0a(3),ppa(3),pma(3),p0b(3),ppb(3),pmb(3),
     @         p0c(3),ppc(3),pmc(3),p0int(3),ppint(3),pmint(3),
     @         test(5)
	integer flag,j,jthN,jphN
	COMMON/GENOUT/PCM,WGT 

	pirad=3.14159

c Extracting perpendicular of (pi+ pi-) plane
c++ Defining matrix element values (from ABBHHM and SLAC data) and deacy_distr where:
c++ th = theta pi+.vector.pi- plane
c++ ph = theta pi+.vector.pi- plane
c++++++++++++++++++++++++++++++++++

c define spin density matrix elements
c	rat = 0.75*(parm_Q2)**(1.09)
c	deno = 1 + parm_EPS*rat

c++ Definition from rhoXY_Z = spin density matrix elements in photoproduction
c	R00_04 =  (rho00_0 + parm_EPS*rat*rho00_4)  / deno
c	ReR10_04 =  (Rerho10_0 + parm_EPS*rat*Rerho10_4)  / deno
c	R1m1_04 = (rho1m1_0 + parm_EPS*rat*rho1m1_4)/ deno

c	R11_1 = rho11_1 / deno
c	R00_1 = rho00_1 / deno
c	ReR10_1 = Rerho10_1 / deno
c	R1m1_1 = rho1m1_1 / deno
c	ImR10_2 = Imrho10_2 / deno
c	ImR1m1_2 = Imrho1m1_2 / deno

c       R11_5 = sqrt(rat)*rho11_5 / deno
c	R00_5 = sqrt(rat)*rho00_5 / deno
c	ReR10_5 = sqrt(rat)*Rerho10_5 / deno
c	R1m1_5 = sqrt(rat)*rho1m1_5 / deno
c	ImR10_6 = sqrt(rat)*Imrho10_6 / deno
c	ImR1m1_6 = sqrt(rat)*Imrho1m1_6 / deno

c++ Values from Joos et al, NP B113 (1976) for Q2 = 1.05 and W = 2-2.8
	R00_04 =  0.47
	ReR10_04 = 0.16
        R1m1_04 = -0.06
	R11_1 = -0.16
	R00_1 = 0.26
	ReR10_1 = 0.07
	R1m1_1 = 0.0
	ImR10_2 = -0.01
	ImR1m1_2 = -0.17
	R11_5 = 0.00
	R00_5 = -0.04
	ReR10_5 = -0.06
	R1m1_5 = -0.04
	ImR10_6 = -0.08
	ImR1m1_6 = 0.09

c find maximum of decay distribution (Nota : neglect epsilon and gphi dependence here !!!!)
c cf equation 92 of Schilling and Wolf, NP B61 (1973)
	decay_distr_max=-1.
	do jthN=0,30
	   do jphN=0,90
	     thN=jthN*3.1416/30.
	     phN=jphN*6.2832/90.
	     decay_distr=3./4./pirad*(
     @                         0.5 * (1-R00_04)
     @                        +0.5 * (3*R00_04-1)*cos(thN)**2
     @                        -1.414 * ReR10_04 * sin(2*thN) * cos(phN)
     @                        -R1m1_04 * sin(thN)**2 * cos(2*phN) ! end W^0
     @                     -parm_EPS*cos(2*parm_GPhi)*(
     @                         R11_1*sin(thN)**2 
     @                        +R00_1*cos(thN)**2
     @                        -1.414*ReR10_1*sin(2*thN)*cos(phN)
     @                        -R1m1_1*sin(thN)**2*cos(2*phN))   ! end W^1
     @                     -parm_EPS*sin(2*parm_GPhi)*(
     @                         1.414*ImR10_2*sin(2*thN)*sin(phN)
     @                        +ImR1m1_2*sin(thN)**2*sin(2*phN)) ! end W^2
     @                     +sqrt(2*parm_EPS*(1+parm_EPS))*cos(parm_GPhi)*(
     @                         R11_5*sin(thN)**2
     @                        +R00_5*cos(thN)**2
     @                        -1.414*ReR10_5*sin(2*thN)*cos(phN)
     @                        -R1m1_5*sin(thN)**2*cos(2*phN))   ! end W^3  
     @                     +sqrt(2*parm_EPS*(1+parm_EPS))*sin(parm_GPhi)*(
     @                         1.414*ImR10_6*sin(2*thN)*sin(phN)
     @                        +ImR1m1_6*sin(thN)**2*sin(2*phN))) ! end W^4
     @              *sin(thN)
	     decay_distr_max=max(decay_distr_max,decay_distr)
	   enddo
	enddo
	decay_distr_max=decay_distr_max*1.5 ! safe factor (to not have fine th ph binning)
	
c randum thetaN and phiN in the omega rest frame 
 765	thN = pirad*rran()
	phN = 2.*pirad*rran()

c evaluate decay distribution
	decay_distr=3./4./pirad*(
     @                         0.5 * (1-R00_04)
     @                        +0.5 * (3*R00_04-1)*cos(thN)**2
     @                        -1.414 * ReR10_04 * sin(2*thN) * cos(phN)
     @                        -R1m1_04 * sin(thN)**2 * cos(2*phN) ! end W^0
     @                     -parm_EPS*cos(2*parm_GPhi)*(
     @                         R11_1*sin(thN)**2 
     @                        +R00_1*cos(thN)**2
     @                        -1.414*ReR10_1*sin(2*thN)*cos(phN)
     @                        -R1m1_1*sin(thN)**2*cos(2*phN))   ! end W^1
     @                     -parm_EPS*sin(2*parm_GPhi)*(
     @                         1.414*ImR10_2*sin(2*thN)*sin(phN)
     @                        +ImR1m1_2*sin(thN)**2*sin(2*phN)) ! end W^2
     @                     +sqrt(2*parm_EPS*(1+parm_EPS))*cos(parm_GPhi)*(
     @                         R11_5*sin(thN)**2
     @                        +R00_5*cos(thN)**2
     @                        -1.414*ReR10_5*sin(2*thN)*cos(phN)
     @                        -R1m1_5*sin(thN)**2*cos(2*phN))   ! end W^3  
     @                     +sqrt(2*parm_EPS*(1+parm_EPS))*sin(parm_GPhi)*(
     @                         1.414*ImR10_6*sin(2*thN)*sin(phN)
     @                        +ImR1m1_6*sin(thN)**2*sin(2*phN))) ! end W^4
     @              *sin(thN)
 
c if flag decay rho = 0 -> flat decay distrib
	if(FLAG.eq.0) then
	 decay_distr = decay_distr_max
	endif

c selection of the extracted event in the decaying particle rest frame
	WRD = rran()
	if((WRD*decay_distr_max).gt.(decay_distr)) goto 765


c Calculating the decay matrix element and the momentum according to omega quantum number (1-)
c++ In the omega rest frame, in the 3pi plane we have pion-dipion system
c++ vec(q) = dipion CMS
c++ vec(p) = remaining pion
c++ mu=vec(q)vec(p)/q/p (angle between p/q)
c++ therefore:
c++   pi0    vec(p)
c++   pi+   -vec(p)/2+vec(q)
c++   pi-   -vec(p)/2-vec(q)
c++ m3 = 3 pion mass system (M_omega) ad Q the available energy
c++   pi0    vec(p)
c++   pi+   -vec(p)/2+vec(q)
c++   pi-   -vec(p)/2-vec(q)
c++ A2 is the sqared matrix element assuming 
c++ A = (vec(pi0)extvec(pi+)) + (vec(pi+)extvec(pi-)) + (vec(pi-)extvec(pi0))
c++ compatible with JP=1- (omega quantum numbers) and minimal Lagrangian (l and L =1)
c++ cf Stevenson et al, PR 125 (1962) 687 and 2008.
c+++++++++++++++++++++++++++

	mpip=0.13957
	mpi0=0.13498
	mpi = (2.*mpip + mpi0)/3.
	m3 = WOM/mpi
	Qav = (m3 - 3.)*mpi 

c randum mu and q2 and extract p2 and A2
 976	mu = rran()*2.+(-1.)
	q2 = rran()*((m3-1.)**2/4.-1.)
	p2 = ((m3+1.)**2-4.*(q2+1))*((m3-1)**2-4.*(q2+1.))/4./m3**2
	A2 = 9.*p2*q2*(1-mu**2)

c selection of the extracted event according to A2
	yx=rran()*50. ! Note : 50 was found numerically with WOM=MOM+10MeV to normalize the extraction
	if (yx.gt.A2) goto 976

	p=0
	if(p2.gt.0) p=sqrt(p2)
	q=0
	if(q2.gt.0) q=sqrt(q2)

c calculating kinetic energy for the 3 pions (Tz->pi0, Tm->pi+ Tp->pi-)
c	Tz=(sqrt(p2+1.)-1.)*mpi
c	Tm=(sqrt(p2/4+q2-p*q*mu+1.)-1.)*mpi
c	Tp=(sqrt(p2/4+q2+p*q*mu+1.)-1.)*mpi

c defining normalized variables
c	y=Tz/Qav
c	x=(Tm-Tp)/sqrt(3.)/Qav

c defining thus deduced pions variables in the plane (label = a)
	th_dipi=acos(mu)
	p0a(1) = p*mpi
	p0a(2) = 0.
	p0a(3) = 0.
	ppa(1) = (-p/2+q*cos(th_dipi))*mpi
	ppa(2) = (q*sin(th_dipi))*mpi
	ppa(3) = 0.
	pma(1) = -p0a(1)-ppa(1)
	pma(2) = -ppa(2)
	pma(3) = 0.

c randum phi_all = azimutal angle defining the x and y axis position on this plane
c (indeed the 3-pions configuration can rotate by an overall phi angle in this plane)
	phi_all=2*pirad*rran()

c defining pions variables after rotation around phi_all (label = b)
	call ROTATION(0.,phi_all,p0a,p0b)
	call ROTATION(0.,phi_all,ppa,ppb)
	call ROTATION(0.,phi_all,pma,pmb)

c defining pions variables after rotation according to the decay matrix (label = c)
c+++ 1) rotating the vectors lying on the plane by theta around X(in the hel frame X=x hcm) axis
	call ROTATION(thN,0.,p0b,p0int)
	call ROTATION(thN,0.,ppb,ppint)
	call ROTATION(thN,0.,pmb,pmint)
c+++ 2) rotating the vectors lying on the plane by phi around Z(in the hel frame Z=omega dir in hcm) axis
	call ROTATION(0.,phN,p0int,p0c)
	call ROTATION(0.,phN,ppint,ppc)
	call ROTATION(0.,phN,pmint,pmc)


c	open(unit=10,file='omeg_dec.out',status='unknown')
c	write(10,*) thN, cos(thN), phN, decay_distr


c filling output variables : Pion1, Pion2 and Pion3 variables in the omega rest frame
	do j=1,3
	   pcm(j,1)=ppc(j)
	   pcm(j,2)=pmc(j)
c	   pcm(j,3)=p0c(j)
	   pcm(j,3)=-ppc(j)-pmc(j) ! to be sure that the energy is conserved, we get the pi0 as a funcion of pi+/pi-
	enddo
	pcm(4,1)=sqrt(ppc(1)**2+ppc(2)**2+ppc(3)**2+mpip**2)
	pcm(5,1)=sqrt(ppc(1)**2+ppc(2)**2+ppc(3)**2)
	pcm(4,2)=sqrt(pmc(1)**2+pmc(2)**2+pmc(3)**2+mpip**2)
	pcm(5,2)=sqrt(pmc(1)**2+pmc(2)**2+pmc(3)**2)
c	pcm(4,3)=sqrt(p0c(1)**2+p0c(2)**2+p0c(3)**2+mpi0**2)
	pcm(4,3)=WOM-pcm(4,1)-pcm(4,2) ! to be sure that the energy is conserved, we get the pi0 as a funcion of pi+/pi-
c	pcm(5,3)=sqrt(p0c(1)**2+p0c(2)**2+p0c(3)**2)
	pcm(5,3)=sqrt(pcm(1,3)**2+pcm(2,3)**2+pcm(3,3)**2) ! to be sure that the energy is conserved, we get the pi0 as a funcion of pi+/pi-

	do j=1,3
	   test(j)=pcm(j,1)+pcm(j,2)+pcm(j,3)
c       print *, test(j)
	enddo 
	test(4)=pcm(4,1)+pcm(4,2)+pcm(4,3)-WOM
c       print *, test(4)
	
	
	return
	end







C*********************************
	SUBROUTINE GAUSS(S,AM,V)
C********************************
	IMPLICIT NONE
	REAL RRAN
	REAL A,AM,V,S
	INTEGER*4 IZ,I
	COMMON/RANDOM/IZ
	A=0.
	DO 50 I=1,12
50	A=A+rran()
	V=(A-6.0)*S+AM
	END    


c        ********************
         subroutine sezallchn
c        ********************
	 IMPLICIT NONE 
	 real sigr(40,251),wr(251),allchn(251)
	 integer j,i,NCH,ICHAN(24)
	 real XSECT(14,50,0:18)
c-------------------------------
         common/sezurto/sigr,wr
         common/somma/allchn
	 COMMON/CHANNELS/NCH,ICHAN
        COMMON/SIGMAS/XSECT
c-------------------------------
c         open(unit=10,file='allchn.out',status='unknown')
c         open(unit=11,file='diffxsec.out',status='unknown')

	 allchn(1)=0
	 do i=1,251
           do j=1,NCH
  	    allchn(i)=allchn(i)+sigr(ichan(j),i)
	   enddo
c	  write(10,*) wr(i), allchn(i)
	 enddo
c	 do i=1,50 ! loop on energy
c	     do j=0,18 ! loop on angle
c	      write(11,*) j,i,xsect(1,i,j) ! n p+c
c	      write(11,*) j,i,xsect(2,i,j) ! p pi0
c	      write(11,*) j,i,xsect(5,i,j) ! delta++ pi-
c	      write(11,*) j,i,xsect(11,i,j) ! rho0  p
c	    enddo
c  	 enddo
15       format(f10.4)
	 return
	 end


c       ********************
        function f_allchn(s)
c       ********************
	 IMPLICIT NONE 
	real sigr(40,251),wr(251),allchn(251),f_allchn,s
	integer j
c------------------------------
        common/sezurto/sigr,wr
        common/somma/allchn
c------------------------------
        f_allchn=0
	do j=1,250
        if(s.ge.wr(j).and.s.le.wr(j+1))
     &   f_allchn=(allchn(j+1)-allchn(j))/0.01*s+
     &   (allchn(j)*wr(j+1)-allchn(j+1)*wr(j))/
     &   0.01
	enddo
        return
	end
       
c	*********************************************************
        subroutine massvol(tetmin,tetmax,omegamin,omegamax,fmax
     &                       ,vol,m_allchn)	
c       *********************************************************
	 IMPLICIT NONE 
        real mp2,Eo,Mp,me,pigr,degrad,raddeg,m_allchn
	real TETMIN,TETMAX,DTET,DOMEGA,omegamin,omegamax,fmax,vol
	real arg,Ee,pe,Q2,d_sigma,W,eps,gamma_v,tete,omega
c------------------------------
	common/costanti/Eo,Mp,me,pigr,degrad,raddeg
c------------------------------

        fmax= 0.
        mp2=mp**2

        TETMIN = TETMIN*DEGRAD
        TETMAX = TETMAX*DEGRAD

c       Set Delta_TETA (deg) and Delta_OMEGA (GeV)
        DTET =(tetmax-tetmin)/100
        DOMEGA = (omegamax-omegamin)/100
                          
	DO Tete = TETMIN, TETMAX, DTET
         DO Omega = OMEGAMIN, OMEGAMAX, DOMEGA

          Ee=Eo-omega
          pe=sqrt(Ee**2-me**2)
          Q2=4.*Eo*(Eo-omega)*(sin(tete/2.))**2
          arg=(Mp**2+2.*Mp*omega-Q2)
          if (arg.lt.1.1) then
	   d_sigma=0.
          else
	   W=sqrt(arg)
           eps = (1.+2.*(1+omega**2/Q2)*(TAN(TETe/2.))**2)**(-1)
           Gamma_V = (1/137./2./pigr**2)*(Ee/Eo)*((W**2-MP**2)/2./MP)/Q2
     &              /(1-EPS)
c          Inserting in d_sigma=Gamma_V*[SIGMA_T+eps*SIGMA_L]
c          SIGMA_T=(Dipole Fit) and SIGMA_L=0
           d_sigma = Gamma_V * m_allchn * (1+Q2/.7)**(-2)
          endif
	  fmax = max(d_sigma,fmax)
	 enddo
        enddo
	vol=fmax*(-cos(tetmax)+cos(tetmin))*(omegamax-omegamin)
          tetmin=tetmin*raddeg
          tetmax=tetmax*raddeg
        return
        end

c	*********************************************************
        subroutine massvol2(wmin,wmax,q2min,q2max,fmax
     &                       ,vol,m_allchn,modarad)	
c       *********************************************************	
	IMPLICIT NONE
        real mp2,Eo,Mp,me,pigr,degrad,raddeg,m_allchn,w2
	real q2max,dw,dq2,costete,wmin,wmax,q2min,fmax,vol
	real Ee,pe,Q2,d_sigma,W,eps,gamma_v,tete,omega      
	real mf_min,mfin,Ttarg,Twi,Twf,Ztarg,Zw, delta
	real d_sigma_part,fe
        real eo_vtx,ee_vtx,w_vtx,q2_vtx,omega_vtx,tete_vtx,fe_vtx
        real eo_ini,ee_ini,w_ini,q2_ini,omega_ini,tete_ini,fe_ini
	real e_radgam, theta_radgam, phi_radgam
	integer modarad,k
	real rcas
        REAL*8   SIGMA_EP_EET
        EXTERNAL SIGMA_EP_EET
	REAL*8    sigma_rad

c------------------------------------
	common/costanti/Eo,Mp,me,pigr,degrad,raddeg
c------------------------------------

        fmax= 0.
        mp2=mp**2

C       Set Delta_W (GeV) and Delta_q2 (GeV2)
        DW = (wmax-wmin)/100
        DQ2 = (q2max-q2min)/100
       DO w = WMIN, WMAX, DW
	  DO Q2 = Q2MIN, Q2MAX, DQ2
	     omega=(w**2+Q2-Mp**2)/2./Mp
	     if(omega.ge.Eo.or.omega.lt.0)then
		d_sigma=0
	     else
		costete=1-Q2/(2.*Eo*(Eo-omega))
		w2=w**2
		if(costete.ge.1.or.costete.le.-1.or.w2.lt.1.1) then
c     &     or.w.lt.1.1) then
		   d_sigma=0
		else
		   tete=acos(costete)
		   Ee=Eo-omega
		   pe=sqrt(Ee**2-me**2)
		   eps = (1.+2.*(1+omega**2/Q2)*(TAN(TETe/2.))**2)**(-1)
		   Gamma_V = (1/137./2./pigr**2)*(Ee/Eo)*((W**2-MP**2)/2./MP)
     &                              /Q2/(1-EPS)
c           Inserting in d_sigma=Gamma_V*[SIGMA_T+eps*SIGMA_L]
c           SIGMA_T=(Dipole Fit) and SIGMA_L=0
************************************************************************************
		   d_sigma = Gamma_V * m_allchn * (1+Q2/.7)**(-2)  
     1	                             /(2.*Eo*(Eo-omega)*Mp)*w   
************************************************************************************
c		   IF(ModaRad.eq.0) THEN
c		      d_sigma = Gamma_V * m_allchn * (1+Q2/.7)**(-2)  
c     1	                                 /(2.*Eo*(Eo-omega)*Mp)*w   
c     		        
c		   ELSEIF(ModaRad.eq.1.or.ModaRad.eq.2) THEN
c		      mf_min = 0.1395679d+0
c		      if(ModaRad.eq.1) THEN
c			 Ttarg  = 0.0 
c			 Twi    = 0.0 
c			 Twf    = 0.0 
c		      else
c			 Ttarg  = 0.00426d+0
c			 Twi    = 0.00087d+0
c			 Twf    = 0.00105d+0
c		      endif
c		      Ztarg  = 1.d+0
c		      Zw     = 29.d+0
c		      delta  = 0.010d+0
cc		      mfin   = 0.93827231+0.1395679+0.1395679
c		      mfin   = 0.93827231+0.1395679
c		      d_sigma =sigma_ep_eet(Eo*1.d+0, EE*1.d+0, TETE*1.d+0)
c     &                         /(2*Eo*(Eo-omega)*Mp)*w
c		      d_sigma =sigma_rad(sigma_ep_eet, eo*1.d+0,ee*1.d+0,tete*1.d+0,
c     &                          mf_min*1.d+0, Ttarg*1.d+0,Twi*1.d+0,Twf*1.d+0,
c     &                          Ztarg*1.d+0,Zw*1.d+0, delta*1.d+0)
c		      d_sigma      = d_sigma      /(2*Eo*(Eo-omega)*Mp)*w
c		   ELSE
c		      print *,' ModaRad is bad=',ModaRad 
c		      stop
c		   ENDIF
************************************************************************************
		endif
	     endif
	     fmax = max(d_sigma,fmax)
	  enddo
        enddo
c safety factor for maximum corss section
c	fmax=fmax*1.2    
	vol=fmax*(wmax-wmin)*(q2max-q2min)
        return
        end



C*******************************************
	SUBROUTINE LUND_WRITE(PCM_N,NP_N)
C*******************************************
	IMPLICIT NONE
	REAL DUMMY,PCM_N(5,18),RHO_CM(5),DNUP(54)
        REAL UE1,VE1,WE1,PE,EE
	REAL Eo,Mp,me,pigr,degrad,raddeg
	INTEGER JCH,K,J,NP_N,PART_ID(40)

C--------------------------
	COMMON/LUND/RHO_CM,DNUP,UE1,VE1,WE1,PE,EE,JCH,
     @              PART_ID
	common/costanti/Eo,Mp,me,pigr,degrad,raddeg
C--------------------------
c+++++ Header for Rho channel
	   IF(JCH.EQ.11) THEN
C++++++ Adding pions variables in CM
            WRITE (2,ERR=100) NP_N+3,DUMMY,DUMMY,EO,DUMMY,JCH*1.
C	               TYPE*, NP_N+3,DUMMY,DUMMY,EO,DUMMY,JCH*1.
              WRITE (2,ERR=100)
     +        0,17,
     +        PCM_N(1,1),PCM_N(2,1),PCM_N(3,1),PCM_N(4,1),0.1395679
              WRITE (2,ERR=100)
     +        0,33,
     +        RHO_CM(1),RHO_CM(2),RHO_CM(3),RHO_CM(4),RHO_CM(5)
C     +        PCM_N(1,2),PCM_N(2,2),PCM_N(3,2),PCM_N(4,2),0.1395679
C	TYPE*,
C     +        0,17,
C     +        PCM_N(1,1),PCM_N(2,1),PCM_N(3,1),PCM_N(4,1),0.1395679
C	TYPE*,
C     +        0,-17,
C     +        PCM_N(1,2),PCM_N(2,2),PCM_N(3,2),PCM_N(4,2),0.1395679
c------
	   ENDIF
c-----
C+++++ Header for others channels
            IF(JCH.NE.11) THEN 
                WRITE (2,ERR=100) NP_N+1,DUMMY,DUMMY,EO,DUMMY,JCH*1.
C	                   TYPE*, NP_N+1,DUMMY,DUMMY,EO,DUMMY,JCH*1.
	    ENDIF
c-----
c+++++ Writing particles
            DO J=1,NP_N
	       K=(J*7)
              WRITE (2,ERR=100)
     +        0,PART_ID(J),
     +        DNUP(1+K),DNUP(2+K),DNUP(3+K),DNUP(4+K),
     +                                  SQRT(DNUP(4+K)**2-DNUP(5+K)**2)
C	      TYPE*,
C     +        0,PART_ID(J),
C     +        DNUP(1+K),DNUP(2+K),DNUP(3+K),DNUP(4+K),
C     +                                  SQRT(DNUP(4+K)**2-DNUP(5+K)**2)
            END DO
c------
c++++++ Writing electron informations
              WRITE (2,ERR=100)
     +        0,7,
     +        PE*UE1,PE*VE1,PE*WE1,EE,ME

C	      TYPE*,
C     +        0,7,
C     +        PE*UE1,PE*VE1,PE*WE1,EE,ME

c---- Closing LUND if-loop
	 RETURN
100     WRITE (6,*)'ERROR IN WRITING LUND FILE'
	 END


C**********************************************
	FUNCTION CHARGE(I)
C*********************************************
	IMPLICIT NONE
	REAL CHARGE
	INTEGER I
	CHARGE=-1000.
	IF(I.EQ.41)  CHARGE=1.
	IF(I.EQ.42)  CHARGE=0.
	IF(I.EQ.-17) CHARGE=-1.
	IF(I.EQ.17)  CHARGE=1.
	IF(I.EQ.23)  CHARGE=0.
	if(I.eq.18)  CHARGE=1.
	if(I.eq.-18)  CHARGE=-1.
	if(I.eq.37)  CHARGE=0.
	if(I.eq.38)  CHARGE=0.
	RETURN 
	END
C**********************************************
	FUNCTION ID_BOS(I)
C*********************************************
	IMPLICIT NONE
	INTEGER ID_BOS,I
	ID_BOS=-99999999
	IF(I.EQ.41)  ID_BOS =  2212
	IF(I.EQ.42)  ID_BOS =  2112  
	IF(I.EQ.-17) ID_BOS = -211
	IF(I.EQ.17)  ID_BOS =  211
	IF(I.EQ.23)  ID_BOS =  111
	if(I.eq.18)   ID_BOS = 321
	if(I.eq.-18)  ID_BOS =-321   
	if(I.eq.37)   ID_BOS = 310
	if(I.eq.38)   ID_BOS = 130

	RETURN 
	END



C*******************************************
	SUBROUTINE ANG_DISTR_1PI(EPSILON,TH,PH)
C*******************************************
c+
c+ NZ = random init value                        INPUT
c+ EPSILON = EPSILON                             INPUT (future use)
c+ TH = CM pion angle                            INPUT (future use)
c+ PH = PH pion according                        OUTPUT
c+       a + b* cos(ph) + c*cos(2*ph)
c-------------------------------------------

	IMPLICIT NONE
        REAL RRAN
	REAL EPSILON, TH,PH

	REAL pi,a,b,c,fmax,f,phmax, fc
	data pi/3.1415/


	a = 139.7
	b =  17.4
	c = -54.5
	
        phmax = acos(-b/4./c)
        fmax = a + b* cos(phmax) + c*cos(2*phmax)	
101     continue       
	 ph = 2*pi*rran()
	 f = (a + b* cos(ph) + c*cos(2*ph))/fmax
	 fc = rran()
	if(fc.gt.f) goto 101
	 
	RETURN

	 END

C*******************************************
	SUBROUTINE SETTORE(ANGLE,SETT)
C*******************************************
c+
c+ ANGLE                  INPUT
c+ SETT                   OUTPUT
c+    
c-------------------------------------------



	REAL ANGLE
	INTEGER SETT
        if (angle.lt.30.or.angle.gt.330)   sett = 1
        if (angle.gt.30.and.angle.lt.90)   sett = 2
        if (angle.gt.90.and.angle.lt.150)  sett = 3
        if (angle.gt.150.and.angle.lt.210) sett = 4
        if (angle.gt.210.and.angle.lt.270) sett = 5
        if (angle.gt.270.and.angle.lt.330) sett = 6

	return
	end


C*******************************************
	SUBROUTINE PSEUDO_SPA(part_t,p,theta,phi,t_current,acp,sect,acc)
C*******************************************
c-      
c-  Pseudo Single Particle acceptance function from Volker Burkert.
c-
c-  Inputs:
c-          part_t    - INTEGER: 0=electron, 1=positive hadron, -1=negative hadron.
c-          p         - REAL:    particle momentum GeV/c.
c-          theta     - REAL:    particle polar angle in degrees (0 to 180.).
c-          phi       - REAL:    particle azimuthal angle in degrees in the sector (-30. to 30.)
c-          t_current - REAL:    torus current in A.
c-	        acp       - REAL:    fraction of phi acceptance; acp = 1. gives the widest range.
c-          sect      - INTEGER: sector number.    
c-
c-  Output:
c-          acc       - REAL:    acceptance - equals 0. if particle is not in the fiducial region, 
c-                               equals the solid angle weight if it is in the fiducial region.
c
      IMPLICIT NONE
      REAL p,theta,phi,t_current,acp    ! Input parameter
      INTEGER part_t,sect                 ! 
c
      REAL acc                          ! Output parameter
c
      REAL t_max
      parameter (t_max = 3375.)
      REAL phi0_el, phi0_nh, phi0_ph
      parameter (phi0_el = 30.)
      parameter (phi0_nh = 25.)
      parameter (phi0_ph = 25.)    ! Changed for better discription.
      REAL theta0_el, theta0_nh, theta0_ph
      parameter (theta0_el = 15.5)
      parameter (theta0_nh = 10.)
      parameter (theta0_ph = 10.)
      REAL thetas_el, thetas_nh, thetas_ph
      parameter (thetas_el = 15.)
      parameter (thetas_nh = 15.)
      parameter (thetas_ph = 25.)
      REAL p_shift, cel_ex, pel_ex, pim_ex
      parameter (p_shift = 0.05)
      parameter (pel_ex = 0.333)
      parameter (pim_ex = 0.5)
      parameter (cel_ex = 0.35)
      REAL  ch_ex,theta_cut
      parameter (theta_cut = 75.)
      parameter (ch_ex = 0.1)
      INTEGER electron,pos_hadron, neg_hadron
      parameter (electron = 0)
      parameter (pos_hadron = 1)
      parameter (neg_hadron = -1)
c- New, for very forward electrons -
      REAL theta_nip,cel_ex_nip,phi0_el_nip
      parameter (theta_nip = 2.)
      parameter (cel_ex_nip = 0.5)
      parameter (phi0_el_nip = 20.)
c
      REAL pi,d2r,u_acc
      parameter (pi = 3.141592653589793)
      parameter (d2r = 0.0174533)
      parameter (u_acc = 0.20944)
      INTEGER err_count
      data err_count/0/
      REAL theta_min, delta_phi, exp, pnorm,dtheta
c-      
      Acc=0.0
      pnorm=p*t_max/t_current
      IF(part_t.EQ.electron) THEN
        if(sect.eq.2)then
          dtheta=6./pnorm
          if(theta.gt.18.+dtheta.and.theta.lt.20.+dtheta)return
        elseif(sect.eq.5)then
          dtheta=6./(p*t_max/t_current)
          if(theta.gt.19.+dtheta.and.theta.lt.21.5+dtheta)return
        endif
        theta_min = theta0_el+thetas_el/(p*t_max/t_current+p_shift)
        if(theta.gt.theta_min.and.theta.lt.50.)then
          exp = cel_ex*(p*t_max/t_current)**pel_ex
          delta_phi = acp*phi0_el*sin((theta-theta_min)*d2r)**exp
          if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
c- New, for very forward electrons -
        elseif(theta.gt.theta_min-theta_nip.and.theta.lt.50.)then
          exp = cel_ex_nip*(p*t_max/t_current)**pel_ex
          delta_phi = acp*phi0_el_nip
     *   *sin((theta-theta_min+theta_nip)*d2r)**exp
          if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
        endif
      ELSEIF(part_t.EQ.pos_hadron) THEN
        if(sect.eq.3)then
	if(pnorm.gt.0.825+1.5*(theta/35.)**6.and.
     *	(pnorm.lt.1.0125+1.35*(theta/32.)**6))return
	if(pnorm.gt.0.15+ 0.75*(theta/65.)**8.and.
     *	(pnorm.lt.0.15+0.75*(theta/50.)**6))return
        elseif(sect.eq.4)then
        if(pnorm.lt.theta/5..and.pnorm.gt.theta/7.)return
        if(pnorm.gt.0.7+1.5*(theta/41.5)**8.and.
     *  pnorm.lt.0.85+1.5*(theta/35.5)**8)return
        elseif(sect.eq.5)then
	if(pnorm.gt.1.+2.*(theta/20.)**2.and.
     *	(pnorm.lt.1.2+1.8*(theta/15.)**2))return
        if(pnorm.gt.0.5+0.7*(theta/40.)**8.and.
     *  pnorm.lt.0.7+2.*(theta/40.)**8)return
        elseif(sect.eq.6)then
         if(pnorm.gt.0.2+0.5*(theta/80.)**6.and.
     *    pnorm.lt.0.25+0.5*(theta/70.)**6)return
        endif
        theta_min = theta0_ph+thetas_ph*(1.-p*t_max/t_current/5.)**24
        if(theta.gt.theta_min)then
          exp=(p*t_max/t_current/5.)**(1./8.)
          delta_phi = acp*phi0_ph*cos((theta-theta_cut)*d2r)**exp
          if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
        endif
      ELSEIF(part_t.EQ.neg_hadron) THEN
        theta_min = theta0_nh+thetas_nh/(p*t_max/t_current+p_shift)
        if(theta.gt.theta_min.and.theta.lt.130.)then
          exp = ch_ex*(p*t_max/t_current)**pim_ex
          delta_phi = acp*phi0_nh*sin((theta-theta_min)*d2r)**exp
          if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
        endif
      ELSE
      	err_count=err_count+1
      	if(err_count.le.10)WRITE(6,*)'Illegal particle type'
      		ACC=0.0
      ENDIF
      RETURN 
      END


C*******************************************
	SUBROUTINE PSEUDO_SPA_MORE1(part_t,p,theta,phi,t_current,acp,sect,acc)
C*******************************************
c-      
c-                     A WIDER version (no holes) of
c-  Pseudo Single Particle acceptance function from Volker Burkert.
c-
c-  Inputs:
c-          part_t    - INTEGER: 0=electron, 1=positive hadron, -1=negative hadron.
c-          p         - REAL:    particle momentum GeV/c.
c-          theta     - REAL:    particle polar angle in degrees (0 to 180.).
c-          phi       - REAL:    particle azimuthal angle in degrees in the sector (-30. to 30.)
c-          t_current - REAL:    torus current in A.
c-	        acp       - REAL:    fraction of phi acceptance; acp = 1. gives the widest range.
c-          sect      - INTEGER: sector number.    
c-
c-  Output:
c-          acc       - REAL:    acceptance - equals 0. if particle is not in the fiducial region, 
c-                               equals the solid angle weight if it is in the fiducial region.
c
      IMPLICIT NONE
      REAL p,theta,phi,t_current,acp    ! Input parameter
      INTEGER part_t,sect                 ! 
c
      REAL acc                          ! Output parameter
c
      REAL t_max
      parameter (t_max = 3375.)
      REAL phi0_el, phi0_nh, phi0_ph
      parameter (phi0_el = 30.)
      parameter (phi0_nh = 25.)
      parameter (phi0_ph = 25.)    ! Changed for better discription.
      REAL theta0_el, theta0_nh, theta0_ph
      parameter (theta0_el = 15.5)
      parameter (theta0_nh = 10.)
      parameter (theta0_ph = 10.)
      REAL thetas_el, thetas_nh, thetas_ph
      parameter (thetas_el = 15.)
      parameter (thetas_nh = 15.)
      parameter (thetas_ph = 25.)
      REAL p_shift, cel_ex, pel_ex, pim_ex
      parameter (p_shift = 0.05)
      parameter (pel_ex = 0.333)
      parameter (pim_ex = 0.5)
      parameter (cel_ex = 0.35)
      REAL  ch_ex,theta_cut
      parameter (theta_cut = 75.)
      parameter (ch_ex = 0.1)
      INTEGER electron,pos_hadron, neg_hadron
      parameter (electron = 0)
      parameter (pos_hadron = 1)
      parameter (neg_hadron = -1)
c- New, for very forward electrons -
      REAL theta_nip,cel_ex_nip,phi0_el_nip
      parameter (theta_nip = 2.)
      parameter (cel_ex_nip = 0.5)
      parameter (phi0_el_nip = 20.)
c
      REAL pi,d2r,u_acc
      parameter (pi = 3.141592653589793)
      parameter (d2r = 0.0174533)
      parameter (u_acc = 0.20944)
      INTEGER err_count
      data err_count/0/
      REAL theta_min, theta_max, delta_phi, exp, pnorm,dtheta
      REAL delta_angles
c-

      delta_angles = 1.5
      
      Acc=0.0
      pnorm=p*t_max/t_current
      IF(part_t.EQ.electron) THEN
        theta_min = theta0_el+thetas_el/(p*t_max/t_current+p_shift)
        theta_min = theta_min - delta_angles
        theta_max = 180. ! 50.
        if(sect.eq.2)then
          dtheta=6./pnorm
          if(theta.gt.18.+dtheta+delta_angles.and.
     &       theta.lt.20.+dtheta-delta_angles)return
        elseif(sect.eq.5)then
          dtheta=6./(p*t_max/t_current)
          if(theta.gt.19.+dtheta+delta_angles.and.
     &       theta.lt.21.5+dtheta-delta_angles)return
        endif
        if(theta.gt.theta_min.and.theta.lt.theta_max)then
          exp = cel_ex*(p*t_max/t_current)**pel_ex
          delta_phi = acp*phi0_el*sin((theta-theta_min)*d2r)**exp
          delta_phi = delta_phi + delta_angles
          if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
c- New, for very forward electrons -
        elseif(theta.gt.theta_min-theta_nip.and.theta.lt.theta_max)then
          exp = cel_ex_nip*(p*t_max/t_current)**pel_ex
          delta_phi = acp*phi0_el_nip
     &   *sin((theta-theta_min+theta_nip)*d2r)**exp   
          delta_phi = delta_phi + delta_angles
          if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
        endif
      ELSEIF(part_t.EQ.pos_hadron) THEN
        theta_min = theta0_ph+thetas_ph*(1.-p*t_max/t_current/5.)**24
        theta_min = theta_min - delta_angles
        theta_max = 180. ! 180.
        if(theta.gt.theta_min.AND.theta.lt.90)then
          exp=(p*t_max/t_current/5.)**(1./8.)
          delta_phi = acp*phi0_ph*cos((theta-theta_cut)*d2r)**exp
          delta_phi = delta_phi + delta_angles
          if(abs(phi).lt.delta_phi)Acc=1
        elseif(theta.gt.90)then
	  Acc=1
        endif
      ELSEIF(part_t.EQ.neg_hadron) THEN
        theta_min = theta0_nh+thetas_nh/(p*t_max/t_current+p_shift)
        theta_min = theta_min - delta_angles
        theta_max = 130. + delta_angles
        if(theta.gt.theta_min.and.theta.lt.theta_max)then
          exp = ch_ex*(p*t_max/t_current)**pim_ex
          delta_phi = acp*phi0_nh*sin((theta-theta_min)*d2r)**exp
          delta_phi = delta_phi + delta_angles
          if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
        endif
      ELSE
      	err_count=err_count+1
      	if(err_count.le.10)WRITE(6,*)'Illegal particle type'
      		ACC=0.0
      ENDIF
      RETURN 
      END



C***********************

	FUNCTION RRAN()
	REAL RRAN
	CALL RANLUX(RRAN,1)
	RETURN
	END



c======================================================================
c     GEN3N1BOD 
c     input/output through commons (see GENBOD)
c     NP - must be 3 - the num of particles -- input
c     TECM = W -- input 
c     AMASS(18)(1,2,3) masses -- input
c     PCM(1,i),PCM(2,i),PCM(3,i)-px,py,pz momenta for i-th part.--output
c     orther parameters in /GENIN / and /GENOUT/ are not used
c======================================================================

      SUBROUTINE gen3n1bod
      IMPLICIT none
      INTEGER NP,KGENEV
      REAL    TECM,AMASS,PCM,WT
      COMMON /GENIN /NP,TECM,AMASS(18),KGENEV
      COMMON /GENOUT/PCM(5,18),WT

      INTEGER nc
      REAL    R,m1,m2,m3
      REAL    pii,w,s12,s23,s12min,s12max,s23min,s23max
      REAL    cteta,fi,psi,teta
      REAL    q1(0:3),q2(0:3),q3(0:3)

      pii = acos(-1.)

      IF(NP.ne.3) THEN
        print *,' sub GEN3N1BOD: NP.ne.3=',NP
        stop
      ENDIF
      w=TECM      
      m1=AMASS(1)
      m2=AMASS(2)
      m3=AMASS(3)
      IF(w.lt.m1+m2+m3) THEN
        print *,' sub GEN3N1BOD: W is too small=',W
        stop
      ENDIF

 1001 continue
      s12min=(m1+m2)**2
      s12max=(W-m3)**2
      s23min=(m2+m3)**2
      s23max=(W-m1)**2
      CALL RANLUX(R,1)
      s12=s12min + R*(s12max-s12min)
      CALL RANLUX(R,1)
      s23=s23min + R*(s23max-s23min)
      CALL RANLUX(R,1)
      cteta=-1. + R*(2.)
      teta=acos(cteta)
      CALL RANLUX(R,1)
      fi = 0.  + R*(2.*pii-0.)
      CALL RANLUX(R,1)
      psi= 0. +  R*(2.*pii-0.)

      CALL kkk3genb3n1(q1,q2,q3,m1,m2,m3,w,s12,s23,teta,fi,psi, nc)
      IF(nc.ne.0) GOTO 1001

      PCM(1,1)=q1(1)
      PCM(2,1)=q1(2)
      PCM(3,1)=q1(3)
      PCM(5,1)=sqrt(pcm(1,1)**2+pcm(2,1)**2+pcm(3,1)**2)
      PCM(4,1)=sqrt(pcm(5,1)**2+m1**2)

      PCM(1,2)=q2(1)
      PCM(2,2)=q2(2)
      PCM(3,2)=q2(3)
      PCM(5,2)=sqrt(pcm(1,2)**2+pcm(2,2)**2+pcm(3,2)**2)
      PCM(4,2)=sqrt(pcm(5,2)**2+m2**2)


      PCM(1,3)=q3(1)
      PCM(2,3)=q3(2)
      PCM(3,3)=q3(3)
      PCM(5,3)=sqrt(pcm(1,3)**2+pcm(2,3)**2+pcm(3,3)**2)
      PCM(4,3)=sqrt(pcm(5,3)**2+m3**2)
  
      RETURN
      END


c------------------------------KKK3GENB--------------------------------
c  for:
c           w -> 1(pi-) + 2(pi+) + 3(Delta) (in CMS)
c  calculates:
c    4-momenta of the particles
c  starting from: w,s12,s23,teta1,fi1,fi23
c  ww      -                                               -- input
c  p1(0:3),p2(0:3),p3(0:3) - 4-momenta                     -- output
c  m1,m2,m3 - masses
c  nc - exit code  0 - OK                                  -- output
c                  1 - out of allowed kin. region
c  look for the comments in kinema3
c----------------------------------------------------------------------

      SUBROUTINE kkk3genb3n1(p1,p2,p3,m1,m2,m3,
     &                       ww,s12,s23,teta1,fi1,fi23, nc)
      IMPLICIT none
      REAL*4 ww,s12,s23
      REAL*4 p1  (0:3),p2  (0:3),p3  (0:3),p1_m,  p2_m,  p3_m
      REAL*4 p1_1(0:3),p2_1(0:3),p3_1(0:3),p1_1_m,p2_1_m,p3_1_m
      REAL*4 TH1,PH1, TH2,PH2, TH3,PH3
      REAL*4 m1,m2,m3,teta1,fi1,fi23,m1_2,m2_2,m3_2
      INTEGER*4 nc
      REAL*4 zero,quasizero,pii, s,s31, g, a,b
      REAL*4 G_BYCKLING4,x,y,z,u,v,w
      

c ----- G-function of Byckling, Byckling, p.89 -----
      G_BYCKLING4(x,y,z,u,v,w) =
     &   (x**2)*y+x*(y**2) + (z**2)*u+z*(u**2) + (v**2)*w+v*(w**2)
     &   + x*z*w + x*u*v + y*z*v + y*u*w
     &   - x*y*(z+u+v+w) - z*u*(x+y+v+w) - v*w*(x+y+z+u)


      zero=0.0d+0
      pii = acos(-1.d+0)
      nc=0
      s = ww**2
      s31 = s - s12 - s23 + m1**2 + m2**2 + m3**2
      m1_2=m1**2
      m2_2=m2**2
      m3_2=m3**2

c ----- -----
      G = G_BYCKLING4(s12,s23,s,m2_2,m1_2,m3_2)
      IF(G.gt.zero) THEN
        nc=1
        RETURN
      ENDIF

c----------------------------------------
c     Moduli of the momenta and energies
c----------------------------------------
      p1(0) = (s+m1_2-s23)/2./ww
      p2(0) = (s+m2_2-s31)/2./ww
      p3(0) = (s+m3_2-s12)/2./ww
      p1_m = sqrt(p1(0)**2-m1_2)
      p2_m = sqrt(p2(0)**2-m2_2)
      p3_m = sqrt(p3(0)**2-m3_2)

c------------------------------------------------------------------
c     Angles and three momenta of the particles in the frame 
c     with the Z-axes along the p1, X-axes perp. to the Z situated
c     on the a-b-1-plane
c------------------------------------------------------------------
      TH1 = 0.d+0
      PH1 = 0.d+0
      a = ((m1_2+m2_2 + 2.*p1(0)*p2(0) - s12)/2./p1_m/p2_m)
      IF(abs(a).gt.1.)THEN
        nc=1
        RETURN
      ENDIF
      TH2 = acos(a)
      PH2 = fi23
      a = ((m1_2+m3_2 + 2.*p1(0)*p3(0) - s31)/2./p1_m/p3_m)
      IF(abs(a).gt.1.)THEN
        nc=1
        RETURN
      ENDIF
      TH3 = acos(a)
      PH3 = fi23 + pii
      IF(PH3 .gt. 2.*pii) PH3 = PH3 - 2.*pii

c------------------------------------
c     Three momenta of the particles 
c------------------------------------

      p1(1) = 0.0
      p1(2) = 0.0
      p1(3) = p1_m
      p2(1) = p2_m*SIN(TH2)*COS(PH2)
      p2(2) = p2_m*SIN(TH2)*SIN(PH2)
      p2(3) = p2_m*COS(TH2)
      p3(1) = p3_m*SIN(TH3)*COS(PH3)
      p3(2) = p3_m*SIN(TH3)*SIN(PH3)
      p3(3) = p3_m*COS(TH3)

c----------------------------------------------------------------------
c     Rotation to the hadronic CMS 
c     (z axis along q vector but x axis is still on the hadronic plane)
c----------------------------------------------------------------------

      p1_1(1) =  COS(teta1)*p1(1) + SIN(teta1)*p1(3)
      p1_1(2) =  p1(2)
      p1_1(3) = (-SIN(teta1))*p1(1) + COS(teta1)*p1(3)
      p1_1(0) = p1(0)
      p1_1_m  = p1_m
      p2_1(1) =  COS(teta1)*p2(1) + SIN(teta1)*p2(3)
      p2_1(2) =  p2(2)
      p2_1(3) = (-SIN(teta1))*p2(1) + COS(teta1)*p2(3)
      p2_1(0) = p2(0)
      p2_1_m  = p2_m
      p3_1(1) = - p1_1(1) - p2_1(1)
      p3_1(2) = - p1_1(2) - p2_1(2)
      p3_1(3) = - p1_1(3) - p2_1(3)
      p3_1(0) = p3(0)
      p3_1_m  = p3_m

c----------------------------------------------------------
c     Rotation to the LAB plane (z axis is still along the 
c     virtual photon but x axis on the ee' plane)
c----------------------------------------------------------

      p1(1) = COS(fi1)*p1_1(1) - SIN(fi1)*p1_1(2)
      p1(2) = SIN(fi1)*p1_1(1) + COS(fi1)*p1_1(2)
      p1(3) = p1_1(3)
      p1(0) = p1_1(0)
      p1_m  = p1_1_m
      p2(1) = COS(fi1)*p2_1(1) - SIN(fi1)*p2_1(2)
      p2(2) = SIN(fi1)*p2_1(1) + COS(fi1)*p2_1(2)
      p2(3) = p2_1(3)
      p2(0) = p2_1(0)
      p2_m  = p2_1_m 
      p3(1) = - p1(1) - p2(1)
      p3(2) = - p1(2) - p2(2)
      p3(3) = - p1(3) - p2(3)
      p3(0) = p3_1(0)
      p3_m  = p3_1_m 

c--------------------------------------------------
c     Check for Kinematics 
c--------------------------------------------------

c ----- 4-momentum conservation:  -----
      quasizero=(ww)*1.d-5
      IF(
     &   abs(p1(1)+p2(1)+p3(1)).gt.quasizero .or.
     &   abs(p1(2)+p2(2)+p3(2)).gt.quasizero .or.
     &   abs(p1(3)+p2(3)+p3(3)).gt.quasizero .or.
     &   abs(p1(0)**2-p1(1)**2-p1(2)**2-p1(3)**2-m1**2).gt.quasizero.or.
     &   abs(p2(0)**2-p2(1)**2-p2(2)**2-p2(3)**2-m2**2).gt.quasizero.or.
     &   abs(p3(0)**2-p3(1)**2-p3(2)**2-p3(3)**2-m3**2).gt.quasizero
     &  ) THEN
        print *,' sub. KKK3: kinematics is bad'
        print *,' p1(1)+p2(1)+p3(1)=',p1(1)+p2(1)+p3(1)
        print *,' p1(2)+p2(2)+p3(2)=',p1(2)+p2(2)+p3(2)
        print *,' p1(3)+p2(3)+p3(3)=',p1(3)+p2(3)+p3(3)
        print *,' m1,3v_p1**2=',m1,
     &        abs(p1(0)**2-p1(1)**2-p1(2)**2-p1(3)**2-m1**2)
        print *,' m2,3v_p2**2=',m2,
     &        abs(p2(0)**2-p2(1)**2-p2(2)**2-p2(3)**2-m2**2)
        print *,' m3,3v_p3**2=',m3,
     &        abs(p3(0)**2-p3(1)**2-p3(2)**2-p3(3)**2-m3**2)
        STOP
      ENDIF
      IF(abs(p1(0)+p2(0)+p3(0)-ww).gt.quasizero ) THEN
        print *,' sub. KKK3: kinematics is bad'
        print *,' ww=',ww
        print *,' p1(0)+p2(0)+p3(0)=',p1(0)+p2(0)+p3(0)
        STOP
      ENDIF

c ----- check s12,s23 -----
      quasizero = (ww)*1.d-5
      a = (p1(0)+p2(0))**2
     &  - (p1(1)+p2(1))**2 - (p1(2)+p2(2))**2 - (p1(3)+p2(3))**2
      b = (p2(0)+p3(0))**2
     &  - (p2(1)+p3(1))**2 - (p2(2)+p3(2))**2 - (p2(3)+p3(3))**2
      IF(abs(a-s12).gt.quasizero .or. abs(b-s23).gt.quasizero) THEN
        print *,' ...sub. KKK3 kinematics is bad'
        print *,' s12_calc,s12=',a,s12
        print *,' s23_calc,s23=',b,s23
        STOP
      ENDIF

      RETURN
      END
      


c----------------------------------------------------------------------
c  CROSS-SECTION 
c----------------------------------------------------------------------

    
      FUNCTION sigma_ep_eet(Es,Ep,Theta)
      IMPLICIT none
      REAL*8   SIGMA_EP_EET,Es,Ep,Theta
      REAL*4   F_ALLCHN, W,omega,q2,eps,gamma_v, MP,PII

      MP=0.93827231
      PII = 3.141592653589793

      omega = Es - Ep                                              
      q2    = 4.*Es*Ep*(sin(theta/2.))**2         
      w     = MP**2 + 2.*MP*omega - q2                    
      if( w.le.0. ) then
        sigma_ep_eet = 0.d+0
        RETURN
      endif
      w     = sqrt(w)
      EPS = (1.+2.*(1+omega**2/q2)*(tan(theta/2.))**2)**(-1)
      GAMMA_V = (1/137./2./PII**2)*(Ep/Es)*((w**2-MP**2)/2./MP)/q2 
     @          /(1-EPS)
      SIGMA_EP_EET = f_allchn(w)* Gamma_V * (1+q2/.7)**(-2)     
      RETURN
      END


c*********************************************************************************************
c
c   Calcul of GPhi = angle between leptonic and hadronic planes              
c
c*********************************************************************************************

	Subroutine GPhiLH(Ebeam,QVlab_Eout,QVlab_Pout,res)
	real*4 Ebeam,QVlab_Eout(3),QVlab_Pout(3)
	real  res
	real*4 kx,ky,kz, pxgv,pygv,pzgv,pgv, pxp,pyp,pzp,pp
	real*4 costhp,thp,costhq,thq,sinthp,sinthq,prodvecx,prodvecy,prodvecz,prodsca
	real*4 cosgphi,singphi

**    incoming electron (k)
	kx=0.
	ky=0.
	kz=Ebeam
**    virtual photon (q)
	pxgv= -QVlab_Eout(1)
	pygv= -QVlab_Eout(2)
	pzgv= Ebeam-QVlab_Eout(3)
	pgv=sqrt(pxgv**2+pygv**2+pzgv**2)
**    outgoing proton 
	pxp = QVlab_Pout(1)
	pyp = QVlab_Pout(2)
	pzp = QVlab_Pout(3)
	pp = sqrt(pxp**2+pyp**2+pzp**2)	

** 
	costhp=(pxgv*pxp + pygv*pyp + pzgv*pzp)/(pgv*pp)
	thp=acos(costhp)
	costhq=(pzgv)/(pgv)
	thq=acos(costhq)
	sinthp=sin(thp)
	sinthq=sin(thq)

**    calcul of sin(GPhi) [prod vect q^k]
	prodvecx=pygv*kz-pzgv*ky 
	prodvecy=-pxgv*Ebeam + pzgv*kx
	prodvecz=pxgv*ky-pygv*kx
	prodsca=prodvecx*pxp+prodvecy*pyp+prodvecz*pzp
	singphi=-prodsca/(Ebeam*pgv*pp*sinthp*sinthq)
*     calcul of cos(Gphi)
	prodsca=Ebeam*pzp
	cosgphi=(costhq*costhp-prodsca/(Ebeam*pp))/(sinthq*sinthp)

	if(singphi.lt.-1.)singphi=-1.
	if(singphi.gt.1.)singphi=1.
	if(cosgphi.gt.1.)cosgphi=1.
	if(cosgphi.lt.-1)cosgphi=-1.
	
	if(singphi.lt.0)then
	  res=360.-acos(cosgphi)*180./acos(-1.)
	else
	  res=acos(cosgphi)*180./acos(-1.)
	endif
	
	RETURN
	END


C*******************************************************************************
	SUBROUTINE ROTATION
     @  (thx,phx,v,vout)
C*******************************************************************************
	implicit none
	REAL V(3),VPRIME(3),DROT(3,3)
	REAL V1,V2,V3,VPRIME1,VPRIME2,VPRIME3
	REAL THETA_E,PHI_E1
	real vout(3),thx,phx
	integer l,k

C**************************************************************
C	Evaluation of the components of the rotation matrix
C**************************************************************

	DROT(1,1) = COS(THX)*COS(PHX)
	DROT(1,2) = -COS(THX)*SIN(PHX)
	DROT(1,3) =   SIN(THX)

	DROT(2,1) = SIN(PHX)
	DROT(2,2) = COS(PHX)
	DROT(2,3) = 0.

	DROT(3,1) = -SIN(THX)*COS(PHX)
	DROT(3,2) = -SIN(THX)*SIN(PHX)
	DROT(3,3) = COS(THX)




	DO L = 1,3
	 DO K = 1,3
	IF(ABS(DROT(L,K)-1.).LT.1.E-6) DROT(L,K) = 1.
	IF(ABS(DROT(L,K)+1.).LT.1.E-6) DROT(L,K) = -1.
	IF(ABS(DROT(L,K)).LT.1.E-4) DROT(L,K) = 0.
	 ENDDO
	ENDDO


C**************************************************************
C	Evaluation of the components of the rotated vector
C**************************************************************

	DO L = 1,3
	VPRIME(L) = 0.
	 DO K = 1,3
	 VPRIME(L) = VPRIME(L) + V(K)*DROT(L,K)
	 ENDDO
	ENDDO

C--------------------------------------------------------
C	Fills the output values with the rotated vector
C--------------------------------------------------------

	do l=1,3
	   vout(l)=VPRIME(l)
	enddo

	RETURN

	END

c*********************************************************************************************
c
c     Fiducial cuts for electrons (same function for data)         
c
c*********************************************************************************************

	Subroutine ElectronFiducialCut(pmom,theta,phi,acc_el)
	real*4 pmom,theta,phi
	real  acc_el
	real*4 C1,C2,C3,C4,factor,pshift,theta_cut,expon,dphi
	integer sett
c	C1 = 12.5
	C1 = 10.
	C2 = 18.55
	C3 = 0.25
c	C4 = 25.
	C4 = 29.
	factor = 0.416667
	pshift = 0.14
	theta_cut = 0.
	expon = 0.
	dphi = 0. 
	acc_el=0.
	p = pmom
	call settore(phi,sett)
	phic=phi-((sett-1)*60.)
	if(sett.eq.1 .and. phi.ge.330)phic=phic-360.
c ******* Fiducial cut theta vs. p
	theta_cut = C1 + C2/(p+pshift)

c******** Fiducial cut phi vs. p
	expon = C3*p**factor
	dphi = C4*(  ( sin((theta-theta_cut)*3.14159/180.) )) **expon     

	if(abs(phic).le.dphi .and.  theta.ge.theta_cut)acc_el=1.

	RETURN
	END

c*********************************************************************************************
c
c   Fiducial cuts for hadrons              
c
c*********************************************************************************************

	Subroutine HadronFiducialCut(pmom,theta,phi,acc_had)
	real*4 pmom,theta,phi,a0p(6),a0m(6),a1p(6),a1m(6),a2p(6),a2m(6),
     +                   a3p(6),a3m(6)
	real acc_had,phic,phi_min,phi_max 
	integer s
	phi_min=0.
	phi_max=0.
	phic=0.
	acc_had=0.
c	data a0p /24.,24.,23.,23.5,24.5,24.5/
c	data a0m /25.,26.,26.,25.5,27.,26./
c	data a1p /0.22,0.23,0.20,0.20,0.22,0.22/
c	data a1m /0.22,0.22,0.22,0.22,0.16,0.16/
c	data a2p /8.,8.,8.,8.,8.,8./
c	data a2m /8.,8.,8.,8.,8.,8./
	data a0p /25.8,25.5,24.7,24.7,26.4,26.4/
	data a0m /26.,26.4,26.9,26.9,28.6,28.6/
	data a1p /0.22,0.24,0.2,0.20,0.22,0.22/
	data a1m /0.22,0.22,0.22,0.22,0.17,0.17/
	data a2p /7.,7.,7.,7.,7.,7./
	data a2m /7.,7.,7.,7.,7.,7./
	data a3p /1.,1.,1.,1.,1.,1./
	data a3m /1.,1.,1.,1.,1.,1./

                call settore(phi,s)
	phic=phi-((s-1)*60.)
	if(s.eq.1 .and. phi.ge.330)phic=phic-360.
c        Fiducial cut phi vs. theta
	phi_min = -a0m(s)*(1-exp(-a1m(s)*(theta-a2m(s))))+a3m(s)
	phi_max = a0p(s)*(1-exp(-a1p(s)*(theta-a2p(s))))+a3p(s)
c	write(6,*)'phi',phi_min,phi_max,phic,s,theta
	if(theta.ge.7. .and. phic.ge.phi_min  .and.  phic.le.phi_max)acc_had=1.
	RETURN
	END
