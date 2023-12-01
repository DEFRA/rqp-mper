*     common blocks for M-PER.FOR ----------------------------------------------
*     maximum number of data-points for a non-parametric distribution ----------
      parameter ( MPRF = 200 )
      parameter ( NPD = 9 )

      integer biometal
      integer rfpd,rcpd,rdpd,efpd,edpd,ecpd,tppd,ttpd,capd ! distributions -----
      integer SITEBIAS,bioERRORS,RQP
      real npval,npcum
      
      double precision udm,uds,tdm,tds,tem,tes,ddm,dds,ubm,ubs,tbm,tbs,
     &tbs3,tbs4,ubs3,wbs3,zbs3,zbs4,rbm,rbs,cbm,cbs,
     &tb9T,tb9Tlc,tb9Tuc,tbc9T,tbc9Tlc,tbc9Tuc,
     &td9T,td9Tlc,td9Tuc,tdc9T,tdc9Tlc,tdc9Tuc   

      double precision ! constants for copper ----------------------------------
     &CU433,CU432,CU431,CU430,CU423,CU422,CU421,CU420,
     &CU413,CU412,CU411,CU410,CU403,CU402,CU401,CU400,
     &CU333,CU332,CU331,CU330,CU323,CU322,CU321,CU320,
     &CU313,CU312,CU311,CU310,CU303,CU302,CU301,CU300,
     &CU233,CU232,CU231,CU230,CU223,CU222,CU221,CU220,
     &CU213,CU212,CU211,CU210,CU203,CU202,CU201,CU200,
     &CU133,CU132,CU131,CU130,CU123,CU122,CU121,CU120,
     &CU113,CU112,CU111,CU110,CU103,CU102,CU101,CU100,
     &CU033,CU032,CU031,CU030,CU023,CU022,CU021,CU020,
     &CU013,CU012,CU011,CU010,CU003,CU002,CU001,CU000

      double precision ! constants for nickel ----------------------------------
     &xni300,xni210,xni201,xni200,xni120,xni111,xni110,xni102, 
     &xni101,xni100,xni030,xni021,xni020,xni012,xni011,xni010, 
     &xni003,xni002,xni001,xni000

      character *30  name of discharge,name of river,name of metal
      character *15  name of metal2

      common /bt/ name of discharge,name of river,name of metal,
     &            name of metal2

      common /vi/ biometal,iteration,IBIG,ishot,
     &ndatsets,numba,erra,nsenum,
     &Ktargit,KPERC,K50,
     &small number,bioERRORS, ! wowowowowowowo
     &bioshotCUm,bioshotCUs1,bioshotCUs2,
     &bioshotZNm,bioshotZNs1,bioshotZNs2,
     &bioshotNIm,bioshotNIs1,bioshotNIs2,
     &bioshotMNm,bioshotMNs1,bioshotMNs2,
     &bioshotPBm,bioshotPBs1,bioshotPBs2,jgo49,
     &bioshotm,bioshots1,bioshots2
     
      common /vr1/ 
     &rfpd,rft,rfs,rfns,       ! upstream river flow
     &rcpd,rct,     ! upstream dissolved metal
     &rc9T,rc9Tlc,rc9Tuc,    

     &rdm,rds,rdpd,rdt,rdns,CDFC1,     ! upstream dissolved organic carbon
     &efpd,eft,efns,     ! discharge flow
     &edm,eds,edpd,edt,edns,CDfc5,     ! discharge quality: dissolved organic carbon
     &nodoc, ! if nodoc = 1 there is no dissolved organic carbon for the discharge
     &ecpd,ect, ! discharge quality 
     &tppd,pht,CPFC1,     ! downstream river: pH  
     &capd,cat,CCFC1,     ! downstream river: calcium
     &CCPC7                            ! correlation between calcium and pH

      common /mr1/                     ! master set of data - for use in sensitivity tests
     &rfmw,rf5w,rfsw,                  ! upstream river flow
     &rcmw,rcsw,rcmlcw,rcmucw,COFC1w,  ! upstream dissolved metal
     &rdmw,rdsw,rdllw,rdulw,CDFC1w,    ! upstream dissolved organic carbon
     &efmw,efsw,COFf2w,                ! discharge flow
     &edmw,edsw,edllw,edulw,CDfc5w,    ! discharge quality: DOC  
     &ecmw,ecsw,COfc5w,ecvw,           ! discharge dissolved metal 
     &phmw,phsw,phllw,phulw,CPFC1w,    ! downstream river: pH  
     &camw,casw,callw,caulw,CCFC1w,    ! downstream river: calcium
     &CCPC7w,                          ! correlation between calcium and pH
     &gofcl,gofclw  ! confidence factor for the PNEC calculation (site-based error)

      common /vrd/ 
     &tdm,tds, ! calculated downstream dissolved metal quality (5)
      
     &tb9T,tb9Tlc,tb9Tuc,
     &tbc9T,tbc9Tlc,tbc9Tuc,
     &td9T,td9Tlc,td9Tuc,   
     &tdc9T,tdc9Tlc,tdc9Tuc,   

     &ddm,dds, ! calculated dissolved organic carbon in downstream river (4)
     &tbm,tbs,tbs3,tbs4,zbs3,zbs4, ! d/s bioavailable metal (11)
     &ubm,ubs,ubs3,wbs3, ! calculated upstream bioavailable metal (12)
     &rbm,rbs, ! calculated mean sampling rate (14) ============================
     &cbm,cbs  ! calculated mean sampling rate (15) ============================

      common /vr3/targit,Covx2,
     &grdm,grds,gedm,geds,
     &gphm,gphs,gcam,gcas,
     &cmaster,bmaster,umaster,tmaster,dmaster,tdmr,
     &cumaster,bumaster,
     &e95c,e99c,e995c,e999c,
     &e95lc,e99lc,e995lc,e999lc,e95uc,e99uc,e995uc,e999uc,
     &tbmr,tbmu,phll,phul,call,caul,edll,edul,rdll,rdul
      
      common /vr2/ ! correlation coefficients ----------------
     &e1,e2,
     &fca1,fca2,fcp1,fcp2,pfails,
     &convirge,Purcentile,attampt
      
      common /curr2/ 
     &tdmc,tdsc,tdc05,tdc90,tdc95,tdc99,
     &tbmc,tbsc,tbc05,tbc90,tbc95,tbc99

      common /a/ CR(15,5000) ! arrays for shots ---------------------
      common /a1/ npdp(NPD)

*     non-parametric distributions ---------------------------------------------
      common /npar/ npval(NPD,mprf),npcum(NPD,mprf)

*     random normal deviates ---------------------------------------------------
      common /b/ RR5,RR6,RR9,RR7M,RR8M, ! wowowowowo
     &R6,R9,R7M,R8M ! wowowowowowowowowo
      common /cy/IR1,IR2,IR3,IR4,IR5,IR6,IR7,IR8,IR9,
     &JR1,JR2,JR3,JR4,JR5,JR6,JR7,JR8,JR9

*     data for an iteration ----------------------------------------------------
      common /i/ trial one,result one,trial two,result two
*     data for sensitivity tests -----------------------------------------------
      common /s/ jrun,krun,kkrun,NSEN,SITEBIAS,RQP
  
      common /zn/ PNECzn,mettal, ! constants for zinc
     &znA,znB,znC,znD,znE,znF,znG,znH,znI,znJ,znK,znL

      common /cu/ ! constants for copper ---------------------------------------
     &CU433,CU432,CU431,CU430,CU423,CU422,CU421,CU420,
     &CU413,CU412,CU411,CU410,CU403,CU402,CU401,CU400,
     &CU333,CU332,CU331,CU330,CU323,CU322,CU321,CU320,
     &CU313,CU312,CU311,CU310,CU303,CU302,CU301,CU300,
     &CU233,CU232,CU231,CU230,CU223,CU222,CU221,CU220,
     &CU213,CU212,CU211,CU210,CU203,CU202,CU201,CU200,
     &CU133,CU132,CU131,CU130,CU123,CU122,CU121,CU120,
     &CU113,CU112,CU111,CU110,CU103,CU102,CU101,CU100,
     &CU033,CU032,CU031,CU030,CU023,CU022,CU021,CU020,
     &CU013,CU012,CU011,CU010,CU003,CU002,CU001,CU000

      common /ni/ ! constants for nickel ---------------------------------------
     &xni300,xni210,xni201,xni200,xni120,xni111,xni110,xni102, 
     &xni101,xni100,xni030,xni021,xni020,xni012,xni011,xni010, 
     &xni003,xni002,xni001,xni000

      common /mn/ ! constants for manganese ------------------------------------
     &xa(9),xb(9),xc(9),xd(9),xe(9),xf(9),xg(9),xh(9),xi(9),xj(9)
      
      common /er/ ! combination of errors --------------------------------------
     &cph1,cph2,cca1,cca2,ced1,ced2,crd1,crd2,cgf1,cgf2,
     &uph1,uph2,uca1,uca2,urd1,urd2,ugf1,ugf2,
     &bph1,bph2,bca1,bca2,bed1,bed2,brd1,brd2,bgf1,bgf2,
     &bph3,bph4,bca3,bca4,bed3,bed4,brd3,brd4,bgf3,bgf4,
     &bph1x,bph2x,bca1x,bca2x,bed1x,bed2x,brd1x,brd2x,bgf1x,bgf2x

      common /gof/ tezt(5)
      common /allow/ icha
      
      common /sens/ bsens(50), csens(50)
      
      common /csv/ rfmlc,rfmuc,rf5lc,rf5uc,rcmlc,rcmuc,rcslc,rcsuc,
     &rc95,rc95lc,rc95uc,rc05lc,rc05uc,
     &phmlc,phmuc,camlc,camuc,rdmlc,rdmuc,
     &phslc,phsuc,caslc,casuc,rdslc,rdsuc,
     &ubmlc,ubmuc,ubslc,ubsuc,ub05lc,ub05uc,ub95,ub95lc,ub95uc,
     &ub9T,ub9Tlc,ub9Tuc,
     &ubmlcr,ubmucr,ubslcr,ubsucr,ub05lcr,ub05ucr,ub95lcr,ub95ucr,
     &efmlc,efmuc,efslc,efsuc,ef95lc,ef95uc,
     &ecmlc,ecmuc,ecslc,ecsuc,
     &edmlc,edmuc,edslc,edsuc,
      
     &ddmlc,ddmuc,ddslc,ddsuc,dd95lc,dd95uc,
      
     &tdmclc,tdsclc,tdc05lc,tdc95lc,
     &tdmcuc,tdscuc,tdc05uc,tdc95uc,
      
     &tbmclc,tbsclc,tbc05lc,tbc95lc,
     &tbmcuc,tbscuc,tbc05uc,tbc95uc, 
      
     &tdmlc,tdslc,td05lc,td95lc,td95,
     &tdmuc,tdsuc,td05uc,td95uc,
      
     &tbmlc,tbslc,tb05lc,tb95lc,tb95,
     &tbmuc,tbsuc,tb05uc,tb95uc,
      
     &temlc,teslc,te80lc,te95lc,te99lc,te995lc,te999lc, 
     &temuc,tesuc,te80uc,te95uc,te99uc,te995uc,te999uc 
      
      common /gener/ FMS(5000)