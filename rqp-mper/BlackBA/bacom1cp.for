*     common blocks for BLACK_AM.FOR
      parameter ( NAR = 16 )
      integer RR,RS,RA1,RA,RB1,RB,forw,sens,type
      integer free,sensa ! bio-available

      character *01 ANS
      character *150 rffile,rcfile,effile,ecfile
      
*     maximum number of data-points for a non-parametric distribution ----------
      parameter ( MPRF = 200 )
      parameter ( NPD = 9 )

      integer biometal
      integer rfpd,rcpd,rdpd,efpd,edpd,ecpd,phpd,capd ! distributions -----
      integer SITEBIAS
      real npval,npcum
      
      double precision udm,uds,tdm,tds,tem,tes,ddm,dds,ubm,ubs,tbm,tbs,
     &tbs3,tbs4,ubs3,wbs3,zbs3,zbs4,rbm,rbs,cbm,cbs,
     &tdc9T,tdc9Tlc,tdc9Tuc   

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
      
      double precision vbm,vbs

      character *30  name of discharge,name of river,name of metal
      character *15  name of metal2 ! ==========================================

      common /v1/ RR,RS,RA1,RA, ! random number seeds
     &NS,CFAC,ITER,IDIL,KDIL,mettal,nsamp,ans
      common /v2/rfm,rf5,reg,add,fsh,ifsh,iregq,rgm,rgs,
     &rcm,rcml,rcmu,rcs,rcsl,rcsu,rcq,rcql,rcqu,rcns,
     &forw,sens,type,efm,efs,udm,uds,sensa,
      
     &phm,phml,phmu,phs,phsl,phsu,phns,cam,caml,camu,
     &cas,casl,casu,cans,
     &udocm,udocml,udocmu,udocs,udocsl,udocsu,udocns,
     &edocm,edocml,edocmu,edocs,edocsl,edocsu,edocns,
     &ddocm,ddocml,ddocmu,ddocs,ddocsl,ddocsu,ddocns,
      
     &rcp,ecp,php,cap,udocp,edocp,
     &tdns,tbns,
     &tbxx,tbxxl,tbxxu,
      
     &nerr,tecsplus,
      
     &ecm,ecml,ecmu,ecs,ecsl,ecsu,ecv,ecns,targ,
     &ec95,ec95l,ec95u,ec99,ec99l,ec99u,ec995,ec995l,ec995u,
     &tcm,tcml,tcmu,tcs,tcsl,tcsu,tcxx,tcxxl,tcxxu,
     &tcx95,tcx99,tcns,
     &free,
     &tecm,tecml,tecmu,tecs,tecsl,tecsu,te80,te80l,te80u,
     &te95,te95l,te95u,
     &te99,te99l,te99u,te995,te995l,te995u,te999,te999l,te999u,
     &grfm,grfs,grgm,grgs,gefm,gefs,gecm,gecs,grcm,grcs,
     &msg1,msg2,msg3,msg4,
     &tecmx,tecsx,

     &COFC1,COFf2,COFc3,COCf4,COfc5,COCc6, ! correlation coefficients
     &capf,cacf,cacp,cafd,caed,

     &b1,b2,c1,c2,c3,d1,d2,d3,d4,p1,p2,f1,f2,z1,z2,y1,y2,

     &iset,XPER,!EN,RN,QN,
     &IHR,IMIN,IDAY,IMON,IYR,
     &Fini,Small,Big,newrun,more,one,
     &KTG,K80,K05,K90,K95,K99,K995,K999, ! array elements for percentiles
     &RR1,RR2,RR3,RR4,RR7,RR8,RRA,RRT,RRD,RRP,RRQ, ! random normal deviates
     &R1,R2,R3,R4,R5,R7,R8,
     &BM(8),BS(8),
     &tem1,out1,tem2,out2 ! data for an iteration

*     Stores for initial variables ...
      common /k1/ XRR,XRS,XRA,XRB,NSX,XPERX,
     &xrfm,xrf5,xreg,xadd,iregqx,xrgm,xrgs,xfsh,ifshx,xrcm,xrcs,
     &freex,
     &xphm,xphs,xphns,xcam,xcas,xcans,
     &xudocm,xudocs,xudocns,xedocm,xedocs,xedocns,
      
     &xgrgm,xgrgs,xefm,xefs,xecm,xecs,xecv,xtarg,
     &xgrcm,xgrcs,xgecm,xgecs,xgphm,xgphs,xgcam,xgcas, 
     &xgudocm,xgudocs,xgedocm,xgedocs,

     &rcmT,ecmT,phmT,camT,udocmT,edocmT,

     &xCOFC1,xCOFf2,xCOFc3,xCOCf4,xCOfc5,xCOCc6,
     &xcapf,xcacf,xcacp,xcafd,xcaed,
     &x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,
     &xa1,xa2,xa3,xa4,xa5,xa6,xa7,xa8,xa9,xa10,xa11,xa12,xa13,xa14,
     &xa15,xa16,xa17

      common /a/ C(NAR,5001) ! arrays for shots
      common /a/ TMC(8,5001) ! arrays for shots

*     non-parametric distributions
      common /npar1/ rfnpvl(1000), rfnpfd(1000), nprf, nprfx,
     &ecnpvl(1000), ecnpfd(1000), nprc, nprcx,
     &npef, npefx, npec, npecx, efnpvl(1000), efnpfd(1000),
     &rcnpvl(1000), rcnpfd(1000), npmax,
     &rffile,rcfile,effile,ecfile
*    &rfnxvl(1000), rfnxfd(1000),ecnxvl(1000), ecnxfd(1000),
*    &efnxvl(1000), efnxfd(1000),rcnxvl(1000), rcnxfd(1000),

      common /bt/ name of discharge,name of river,name of metal,
     &            name of metal2

      common /vi/ biometal,iteration,IBIG,
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
     &rcpd,rct,                ! upstream dissolved metal

     &rdpd,rdt,rdns,           ! upstream dissolved organic carbon
     &efpd,eft,efns,           ! discharge flow
     &edpd,edt,edns,CDfc5,     ! discharge quality: dissolved organic carbon
     &nodoc,                   ! if nodoc = 1 there is no dissolved organic carbon for the discharge
     &ecpd,ect,                ! discharge quality 
     &phpd,pht,                ! downstream river: pH  
     &capd,cat                 ! downstream river: calcium

      common /mr1/                    
     &gofcl   ! confidence factor for the PNEC calculation (site-based error)

      common /vrd/ 
     &tdm,tds, ! calculated downstream dissolved metal quality (5)
     &tdc9T,tdc9Tlc,tdc9Tuc,   
     &ddm,dds, ! calculated dissolved organic carbon in downstream river (4)
     &tbm,tbs,tbs3,tbs4,zbs3,zbs4, ! d/s bioavailable metal (11)
     &ubm,ubs,ubs3,wbs3, ! calculated upstream bioavailable metal (12)
     &vbm,vbs, ! vbs3, u/s biometal using downstream DOC (16)  
     &rbm,rbs, ! calculated mean sampling rate (14) ============================
     &cbm,cbs  ! calculated mean sampling rate (15) ============================

      common /vr3/targit,Covx2,
     &gudocm,gudocs,gedocm,gedocs,
     &gphm,gphs,gcam,gcas,
     &cmaster,
     &tbmr ! ,phll,phul,call,caul,edll,edul,rdll,rdul ##########
      
      common /vr2/ ! correlation coefficients ----------------
     &e1,e2,fca1,fca2,fcp1,fcp2,pfails,
     &Purcentile
      
      common /curr2/ 
     &tdmc,tdsc,tdc05,tdc90,tdc95,tdc99,
     &tbmc,tbsc,tbc05,tbc90,tbc95,tbc99

      common /a/ CR(15,5000) ! arrays for shots ---------------------
      common /a1/ npdp(NPD)

*     non-parametric distributions ---------------------------------------------
      common /npar/ npval(NPD,mprf),npcum(NPD,mprf)

*     random normal deviates ---------------------------------------------------
      common /b/ RR5,RR6,RR9,R6,R9 
      common /cy/IR1,IR2,IR3,IR4,IR5,IR6,IR7,IR8,IR9,
     &JR1,JR2,JR3,JR4,JR5,JR6,JR7,JR8,JR9
     &KR1,KR2,KR3,KR4,KR5,KR6,KR7,KR8,
     &LR1,LR2,LR3,LR4,LR5,LR6,LR7,LR8

*     data for sensitivity tests -----------------------------------------------
      common /s/ kkrun,SITEBIAS,nsampsim,iss,nbias,inkbias
  
      common /zn/ ! constants for zinc
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
      
      common /csv/ rfmlc,rfmuc,rf5lc,rf5uc,
     &rc95,rc95lc,rc95uc,rc05lc,rc05uc,
     &phmlc,phmuc,camlc,camuc,ddocmlc,ddocmuc,
     &phslc,phsuc,caslc,casuc,ddocslc,ddocsuc,
     &ubml,ubmu,ubsl,ubsu,ub05l,ub05u,ub95,ub95l,ub95u,
     &ub9T,ub9Tlc,ub9Tuc,ubxx,ubxxl,ubxxu,
     &ubmlcr,ubmucr,ubslcr,ubsucr,ub05lcr,ub05ucr,ub95lcr,ub95ucr,
     &ubsplus,vbsplus,
     &vbml,vbmu,vbsl,vbsu,vb05l,vb05u,vb95,vb95l,vb95u,
     &vb9T,vb9Tlc,vb9Tuc,vbxx,vbxxl,vbxxu,
     &vbmlcr,vbmucr,vbslcr,vbsucr,vb05lcr,vb05ucr,vb95lcr,vb95ucr,
     &efmlc,efmuc,efslc,efsuc,ef95lc,ef95uc,
     &ecmlc,ecmuc,ecslc,ecsuc,
     &edocmlc,edocmuc,edocslc,edocsuc,
      
     &ddmlc,ddmuc,ddslc,ddsuc,dd95lc,dd95uc,
      
     &tdml,tdmu,tdsl,tdsu,tdxx,tdxxl,tdxxu,
      
     &tdmclc,tdsclc,tdc05lc,tdc95lc, ! ###############################
     &tdmcuc,tdscuc,tdc05uc,tdc95uc,
      
     &tbml,tbsl,tb05l,tb95l,
     &tbmu,tbsu,tb05u,tb95u, 
      
     &td05lc,td95lc,td95,
     &td05uc,td95uc,
      
     &tb05lc,tb95lc,tb95,
     &tb05uc,tb95uc,
     &tbmplus,tbsplus,
      
     &temlc,teslc,te80lc,te95lc,te99lc,te995lc,te999lc, 
     &temuc,tesuc,te80uc,te95uc,te99uc,te995uc,te999uc 
