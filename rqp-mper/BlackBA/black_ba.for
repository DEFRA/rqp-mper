*     Fortran subroutine for Mass Balance Calculations -------------------------      
*     Called as a Dynamic Link Library (DLL) from Visual Basic interface -------
*     Used for Bio-available (BA) model ----------------------------------------
*     Tony Warn - November 2017 ------------------------------------------------
*     **************************************************************************

      subroutine MPER(vrfm,vrf5,vrcm,vrcs,vrcns,vrcml,vrcmu,
     &vrcsl,vrcsu,vrcq,vrcql,vrcqu,
     &vefm,vefs,vecm,vecs,vecns,vecml,vecmu,
     &vecsl,vecsu,vec95,vec95l,
     &vec95u,vec99,vec99l,vec99u,vec995,vec995l,vec995u,
     &vforw,vtarg,vxper,vtype,
     &vcofc1,vcoff2,vcofc3,vcocf4,vcofc5,vcocc6,     
     &vsens,vmsg1,vmsg2,vmsg3, 
     &vtcm,vtcs,vtcns,vtcml,vtcmu,vtcsl,vtcsu, 
     &vtcx,vtcxl,vtcxu,vtcxx,vtcxxl,vtcxxu, 
     &vtecm,vtecs,vtecml,vtecmu,vtecsl,vtecsu,vte95,vte95l,vte95u, 
     &vte99,vte99l,vte99u,vte995,vte995l,vte995u,
     &vx1,vx2,vx3,vx4,vx5,vx6,vx7,vx8,vx9,vx10,vx11,
     &vphm,vphs,vphns,vphml,vphmu,vphsl,vphsu,
     &vcam,vcas,vcans,vcaml,vcamu,vcasl,vcasu,
     &vudocm,vudocs,vudocns,vudocml,vudocmu,vudocsl,vudocsu,
     &vedocm,vedocs,vedocns,vedocml,vedocmu,vedocsl,vedocsu,
     &vddocm,vddocs,vddocns,vddocml,vddocmu,vddocsl,vddocsu,
     &vcapf,vcacf,vcacp,vcafd,vcaed,
     &vtbm,vtbml,vtbmu,vtbs,vtbsl,vtbsu,vtbxx,vtbxxl,vtbxxu,
     &vubm,vubml,vubmu,vubs,vubsl,vubsu,vubxx,vubxxl,vubxxu,
     &vudm,vuds,vudq,
     &vsensa,vfree,
     &vmsg4,vxa1,vxa2,vxa3,vxa4,vxa5,vxa6,vxa7,vxa8,vxa9,vxa10,vxa11,
     &vxa12,vxa13,vxa14,vxa17,vmettal)
*    &vadd,vreg,vrgm,vrgs,viregq,
*    &vrffile,vrcfile,veffile,vecfile,
      

*     ==========================================================================
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL,REFERENCE, ALIAS : "MPER" :: MPER

      integer vforw,viregq,vmsg1,vmsg2,vmsg3,vsens,vmsg4,vtype,vfree
      integer vrcns,vecns,vphns,vcans,vudocns,vedocns,vsensa,vmettal
      character*150 vrcfile,vrffile,veffile,vecfile,vidfile
      character *1 zed

      include 'bacom1cp.for'
      dimension xsampsim(5,5000)
      
      !nbias to be 0,1,0
      
      msampsim = 0 ! number of sets of samples to be simulated SSSSSSSSSSSSSSSSS
      newway = 0 ! use river quality confidence limits #########################
      NERR = 1 ! use conf.limits to assess errors linked to pH, Ca etc =========

      newrun = 1 ! this is a new run
      nbias = 0 ! a value of 0 triggers calculation of bias (BM and BS) 
      free = vfree
      if ( free .eq. 0 ) NERR = 0
      mettal = vmettal
      call check type of metal
      call write the heading
      
      if ( biometal .eq. 2 ) call metal data
      call set up data (vrfm,vrf5,vrcm,vrcs,vrcns,vefm,vefs,vecm,vecs,
     &vecns,vforw,vtarg,vxper,vtype,vcofc1,vcoff2,vcofc3,vcocf4,vcofc5,
     &vcocc6,vsens,vphm,vphs,vphns,vcam,vcas,vcans,vudocm,vudocs,
     &vudocns,vedocm,vedocs,vedocns,vcapf,vcacf,vcacp,vcafd,vcaed,
     &vsensa,vfree)
      
      if ( efm .lt. 0.00001 ) then
      call zero flow (vsens,vsensa,vforw,vrcml,vrcmu,vrcsl,vrcsu,
     &vrcq,vrcql,vrcqu,vtcm,vtcml,vtcmu,vtcs,vtcsl,vtcsu,
     &vtcq,vtcql,vtcqu,vtcns,
     &vecm,vecml,vecmu,vecs,vecsl,vecsu,
     &vec95,vec95l,vec95u,vec99,vec99l,vec99u,
     &vec995,vec995l,vec995u,vtecm,vtecml,vtecmu,
     &vtecs,vtecsl,vtecsu,vte95,vte95l,vte95u,vte99,vte99l,vte99u,
     &vte995,vte995l,vte995u)
      return
      endif

      call set the starting values of the data   
      call do the calculations (0) ! do the calculations
      
*     investigate the effect of rates of sampling SSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
*     set the random number starters SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      LR1 = KR1 ! river flow
      LR3 = KR3 ! discharge flow
      LR2 = KR2 ! river quality for dissolved metal 
      LR4 = KR4 ! discharge quality for dissolved metal
      LR5 = KR5 ! discharge quality for dissolved organic carbon 
      LR6 = KR6 ! dissolved organic carbon in the upstream river
      LR7 = KR7 ! pH in downstream river
      LR8 = KR8 ! calcium in downstream river
*     investigate the effect of rates of sampling SSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      nsampsim = 0
      if ( msampsim .gt. 0 ) then 
      nbias = 1 ! suppress a re-calculation of bias (BM and BS) 
      nsampsim = msampsim
      rcmT = 0.0 ! initialise the summaries of results of all sets of samples
      rcsT = 0.0
      ecmT = 0.0
      ecsT = 0.0
      phmT = 0.0
      phsT = 0.0
      camT = 0.0
      casT = 0.0
      udocmT = 0.0
      udocsT = 0.0
      edocmT = 0.0
      edocsT = 0.0
      tbmT = 0.0
      tbsT = 0.0
      tbmlT = 0.0
      tbmuT = 0.0
      tdmT = 0.0
      tdcvT = 0.0
      tdmlT = 0.0
      tdmuT = 0.0
      tecmT = 0.0
      tecsT = 0.0
      
      trcp = 0.0
      tecp = 0.0
      tphp = 0.0
      tcap = 0.0
      tudocp = 0.0
      tedocp = 0.0

      tecseT = 0.0 ! standard deviation from calculated means SSSSSSSSSSSSSSSSSS
      teccvT = 0.0 ! coefficient of variation SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      do iss = 1, nsampsim ! perform the simulations SSSSSSSSSSSSSSSSSSSSSSSSSSS
      do jss = 1, 4 ! initalise the stores for results of this simulation SSSSSS
      xsampsim(jss,iss) = 0.0
      enddo
      write(27,4000)iss
 4000 format(/100('=')/'Run number 27',i5/100('-'))
      write(3,4010)iss
 4010 format(/68('=')/'Run number',i5)
      call reset the starting data
      call take the samples ! select the sets of samples SSSSSSSSSSSSSSSSSSSSSSS
      
      rcmT = rcmT + rcm ! accumulate values for the global mean SSSSSSSSSSSSSSSS
      rcsT = rcsT + rcs ! and the global standard deviation SSSSSSSSSSSSSSSSSSSS
      ecmT = ecmT + ecm
      ecsT = ecsT + ecs
      phmT = phmT + phm
      phsT = phsT + phs
      camT = camT + cam ! calcium
      casT = casT + cas
      udocmT = udocmT + udocm
      udocsT = udocsT + udocs
      edocmT = edocmT + edocm
      edocsT = edocsT + edocs
      
      trcp = trcp + rcp
      tecp = tecp + ecp
      tphp = tphp + php
      tcap = tcap + cap
      tudocp = tudocp + udocp
      tedocp = tedocp + edocp

      tbmT = tbmT + tbm ! downstream bio-metal SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      tbsT = tbsT + tbs ! downstream bio-metal SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      tbmlT = tbmlT + tbml
      tbmuT = tbmuT + tbmu
      tdmT = tdmT + tdm ! downstream total metal SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      tdcvT = tdcvT + tds/tdm
      tdmlT = tdmlT + tdml
      tdmuT = tdmuT + tdmu
      tecmT = tecmT + tecm
      tecsT = tecsT + tecs
      tecseT = tecseT + tecm*tecm
      teccvT = teccvT + tecs/tecm
      ecv = ecs / ecm

      call do the calculations (0) ! do the calculations SSSSSSSSSSSSSSSSSSSSSSS
      write(27,5031)tdm,tds
 5031 format(2f8.2,' dissolved metal d/s of the discharge')
      write(27,5131)tbm,tbs
 5131 format(2f8.2,' bio-available d/s of the discharge')
      write(27,5132)
 5132 format(100('-'))
      write(27,5133)tecm,tecs,tecs/tecm
 5133 format(2f8.2,' required mean and standard deviation quality for ',
     &'the discharge',12x,3f9.4)
      write(27,5142)
 5142 format(100('='))
      
      xsampsim(1,iss) = tecm ! store the required discharge quality SSSSSSSSSSSS
      xsampsim(2,iss) = tecs
      xsampsim(3,iss) = tbml ! store the downstream bio-metal SSSSSSSSSSSSSSSSSS
      xsampsim(4,iss) = tbmu
      xsampsim(5,iss) = tbs
      enddo ! end of this simulation SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

*     all simulations of sampling have been completed SSSSSSSSSSSSSSSSSSSSSSSSSS
      rcmT = rcmT/float(nsampsim) ! compute the mean from all sets
      rcsT = rcsT/float(nsampsim) ! and the standard deviation
      ecmT = ecmT/float(nsampsim)
      ecsT = ecsT/float(nsampsim)
      phmT = phmT/float(nsampsim)
      phsT = phsT/float(nsampsim)
      camT = camT/float(nsampsim)
      casT = casT/float(nsampsim)
      udocmT = udocmT/float(nsampsim)
      udocsT = udocsT/float(nsampsim)
      edocmT = edocmT/float(nsampsim)
      edocsT = edocsT/float(nsampsim)
      tbmT = tbmT/float(nsampsim)
      tbsT = tbsT/float(nsampsim)
      tbmlT = tbmlT/float(nsampsim)
      tbmuT = tbmuT/float(nsampsim)
      tdmT = tdmT/float(nsampsim)
      tdcvT = tdcvT/float(nsampsim)
      tdmlT = tdmlT/float(nsampsim) ! the average of the lower confidence limit
      tdmuT = tdmuT/float(nsampsim)
      
      trcp = trcp/float(nsampsim)
      tecp = tecp/float(nsampsim)
      tphp = tphp/float(nsampsim)
      tcap = tcap/float(nsampsim)
      tudocp = tudocp/float(nsampsim)
      tedocp = tedocp/float(nsampsim)

      XXX = tecseT - tecmT*tecmT/FLOAT(nsampsim)
      if (XXX .le. Small) then 
      tecseT = 0.0
      else
      tecseT = SQRT (XXX/FLOAT(nsampsim-1))
      endif
      tecmT = tecmT/float(nsampsim)
      tecsT = tecsT/float(nsampsim)
      teccvT = teccvT/float(nsampsim)

      do 3 I = nsampsim, 2, -1 ! sequence the results SSSSSSSSSSSSSSSSSSSSSSSS
      do 6 J = 1, I - 1
      if (xsampsim(1,I) .gt. xsampsim(1,J)) goto 6
      CC = xsampsim(1,I)
      xsampsim(1,I) = xsampsim(1,J)
      xsampsim(1,J) = CC
      CC = xsampsim(2,I)
      xsampsim(2,I) = xsampsim(2,J)
      xsampsim(2,J) = CC
    6 continue
      do 26 J = 1, I - 1
      if (xsampsim(4,I) .gt. xsampsim(4,J)) goto 26
      CC = xsampsim(4,I)
      xsampsim(4,I) = xsampsim(4,J)
      xsampsim(4,J) = CC
      CC = xsampsim(3,I)
      xsampsim(3,I) = xsampsim(3,J)
      xsampsim(3,J) = CC
   26 continue
      do 36 J = 1, I - 1
      if (xsampsim(5,I) .gt. xsampsim(5,J)) goto 36
      CC = xsampsim(5,I)
      xsampsim(5,I) = xsampsim(5,J)
      xsampsim(5,J) = CC
   36 continue
    3 continue ! sequence the results SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

      write(7,6010)nsampsim,rcmT,xrcm,100.0*(rcmT-xrcm)/xrcm,
     &int(rcns)
 6010 format(/80('=')/'Summary of results of',i6,' simulations'/80('-')/
     &12x,'calculated     input         change .........'/80('-')/
     &'  mean rcm =',2f10.4,' ... ',f8.2,' %',18x,i7,' samples')
      write(7,6834)rcsT,xrcs,100.0*(rcsT-xrcs)/xrcs,trcp
 6834 format('  sdev rcs =',2f10.4,' ... ',2(f8.2,' %')/80('-'))
      write(7,6012)phmT,xphm,100.0*(phmT-xphm)/xphm,int(phns)
 6012 format('  mean phm =',2f10.4,' ... ',f8.2,' %',18x,i7,' samples')
      write(7,6394)phsT,xphs,100.0*(phsT-xphs)/xphs,tphp
 6394 format('  sdev phs =',2f10.4,' ... ',2(f8.2,' %')/80('-'))
      write(7,6013)camT,xcam,100.0*(camT-xcam)/xcam,int(cans)
 6013 format('  mean cam =',2f10.4,' ... ',f8.2,' %',18x,i7,' samples')
      write(7,6884)casT,xcas,100.0*(casT-xcas)/xcas,tcap
 6884 format('  sdev cas =',2f10.4,' ... ',2(f8.2,' %')/80('-'))
      write(7,6014)udocmT,xudocm,100.0*(udocmT-xudocm)/xudocm,
     &int(udocns)
 6014 format('mean udocm =',2f10.4,' ... ',f8.2,' %',18x,i7,' samples')
      write(7,6814)udocsT,xudocs,100.0*(udocsT-xudocs)/xudocs,tudocp
 6814 format('sdev udocs =',2f10.4,' ... ',2(f8.2,' %')/80('-'))
      write(7,6015)edocmT,xedocm,100.0*(edocmT-xedocm)/xedocm,
     &int(edocns)
 6015 format('mean edocm =',2f10.4,' ... ',f8.2,' %',18x,i7,' samples')
      write(7,6815)edocsT,xedocs,100.0*(edocsT-xedocs)/xedocs,tedocp
 6815 format('sdev edocs =',2f10.4,' ... ',2(f8.2,' %')/80('-'))
      write(7,6011)ecmT,xecm,100.0*(ecmT-xecm)/xecm,int(ecns)
 6011 format('  mean ecm =',2f10.4,' ... ',f8.2,' %',18x,i7,' samples')
      write(7,6311)ecsT,xecs,100.0*(ecsT-xecs)/xecs,tecp
 6311 format('  sdev ecs =',2f10.4,' ... ',2(f8.2,' %')/80('=')/)
      
      write(7,6046)tdmT,tdcvT,tdmlT,tdmuT
 6046 format(/80('=')/
     &'  mean tdm =',f10.4,5x,'mean cov =',f10.4,7x,
     &'(',f6.2,' -',f6.2,')  mean CLs'/80('=')/)
   
      i1 = 0.05 * nsampsim ! identify 5% biggest SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      i2 = 0.95 * nsampsim ! identify 95% biggest SSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      i3 = 0.50 * nsampsim ! identify the median SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

      write(7,7218)tbmT,xsampsim(3,i3),xsampsim(4,i3)
 7218 format(/80('=')/'downstream bio-quality ...'/80('-')/
     &'  mean tbm =',f10.4,5x,f9.4,'  -',f9.4,' ... (90% range ',
     &'by ranking) '/80('-'))
      !write(7,6016)tbmlT,tbmuT
 6016 format(25x,f9.4,'  -',f9.4,
     &' ... (average confidence limits)')
      write(7,7258)tbsT,xsampsim(5,i1),xsampsim(5,i2)
 7258 format(
     &'  mean tbs =',f10.4,5x,f9.4,'  -',f9.4,' ... (90% range ',
     &'by ranking) '/80('='))

      write(7,6217)tecmT,tecsT,int(ecns) ! SSSSSSSSSSSSS required discharge mean
 6217 format(/80('=')/'required discharge mean (tecm)'/80('-')/
     &' mean tecm =',f10.4,9x,'standard dev =',
     &f10.4,i17,' samples'/80('-'))
      T0=TDIST1(ecns-1.0,0.95)
      tecml=amax1(0.0,tecmT-tecseT*T0)
      tecsu=tecmT+tecseT*T0
      !write(7,6018)tecml,tecsu
 6018 format(
     &80('-')/11x,f9.4,'  -',f9.4,' ... (90% range calculated ',
     &'using standard error)')
      write(7,6218)xsampsim(1,i3),xsampsim(1,i1),xsampsim(1,i2)
 6218 format(' mean tecm =',2f10.4,'  -',f9.4,' ... (',
     &'calculated by ranking)'/80('=')) ! SSSSSSSS required discharge mean

      nsampsim = 0 ! sampling calculations finished: reset controls
      nbias = 1 ! a value of 0 triggers calculation of bias (BM and BS)
      iss = 0
      
      write(17,9482)(xsampsim(1,I),i=1,msampsim)
 9482 format(20f7.2)
      write(17,*)' ' ! biosumm.sm1
      write(17,9482)(xsampsim(5,I),i=1,msampsim)
      write(17,*)' '
      write(17,9482)(xsampsim(3,I),i=1,msampsim)

      call reset the starting data
      call set the starting values of the data       
      call do the calculations (0) ! now do the basic calculations

      endif ! sampling calculations finished SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      
      write(7,7412)tecml,tecmu,tecm
 7412 format(//80('=')/'with basic sampling errors ...'/80('-')/
     &11x,f9.2,'  -',f9.2,' ... (from MPER, dis mean discharge:',
     &f11.2,')')
      write(7,7482)tdml,tdmu,tdm
 7482 format(
     &11x,f9.2,'  -',f9.2,' ... (from MPER, dis mean d/s river:',
     &f11.2,')')
      write(7,7489)tbml,tbmu,tbm
 7489 format(
     &11x,f9.2,'  -',f9.2,' ... (from MPER, bio mean d/s river:',
     &f11.2,')')
      
      diff1 = tecmu - tecml
      write(3,7422)tecml,tecmu,tecm
      write(6,7422)tecml,tecmu,tecm
 7422 format(95('=')/
     &'with basic sampling errors ...'/95('-')/
     &11x,f9.3,'  -',f9.3,' ... (from MPER, dis mean discharge:',
     &f11.3,')')
      write(3,7432)tdml,tdmu,tdm
      write(6,7432)tdml,tdmu,tdm
 7432 format(
     &11x,f9.3,'  -',f9.3,' ... (from MPER, dis mean d/s river:',
     &f11.3,')')
      write(3,7589)tbml,tbmu,tbm
      write(6,7589)tbml,tbmu,tbm
 7589 format(
     &11x,f9.3,'  -',f9.3,' ... (from MPER, bio mean d/s river:',
     &f11.3,')'/95('-'))
      
      if ( nerr .eq. 1 .and. free .eq. 1) then ! ===============================
      call add errors from pH etc
      
      diff2 = tecmu - tecml
      
      perchange = 100.0 * (diff2 - diff1) / diff1

      write(3,7411)tecml,tecmu,tecm
      write(7,7411)tecml,tecmu,tecm
 7411 format(80('-')/'with extra sampling errors from pH etc ...'/
     &80('-')/
     &11x,f9.2,'  -',f9.2,' ... (from MPER, dis mean discharge:',
     &f11.2,')')
      write(3,7582)tdml,tdmu,tdm
      write(7,7582)tdml,tdmu,tdm
 7582 format(
     &11x,f9.2,'  -',f9.2,' ... (from MPER, dis mean d/s river:',
     &f11.2,')')
      write(3,4589)tbml,tbmu,tbm
      write(7,4589)tbml,tbmu,tbm
 4589 format(
     &11x,f9.2,'  -',f9.2,' ... (from MPER, bio mean d/s river:',
     &f11.2,')'/80('='))
      write(3,2589)perchange
 2589 format('Percentage change =',f6.2,'%'/80('=')//)
      
      endif ! if ( nerr .eq. 1 .and free .eq. 1) ===============================
      
      vtecml = tecml
      vtecmu = tecmu
      vtecsl = tecsl
      vtecsu = tecsu

*     use confidence limits on downstream river quality to calculate upper *****
*     confidence limit on the required effluent quality ************************
      if ( newway .eq. 1 ) then ! **********************************************
      targitlo = tbml ! first target
      targithi = tbmu ! second target
      call set the starting values of the data   
      targ = targitlo
      call do the calculations (0) ! newway - not used
      tecmlo = tecm ! lower confidence limit on effluent quality
      tecslo = tecs ! lower confidence limit on effluent quality
      call set the starting values of the data 
      targ = targithi
      call do the calculations (0) ! newway - not used
      tecmhi = tecm ! upper confidence limit on effluent quality
      tecshi = tecs ! lower confidence limit on effluent quality
      call set the starting values of the data 
      targ = xtarg
      call do the calculations (0) ! newway - not used
      vtecml = tecmlo
      vtecmu = tecmhi
      vtecsl = tecslo
      vtecsu = tecshi
      endif ! if ( newway .eq. 1 ) *********************************************
          
      vrfm = rfm ! store the outputs for transport back to the VB front-end
      vrf5 = rf5
      vrcm = rcm
      vrcml = rcml
      vrcmu = rcmu
      vrcs = rcs
      vrcsl = rcsl
      vrcsu = rcsu
      vrcq = rcq
      vrcql = rcql
      vrcqu = rcqu
      vefm = efm
      vefs = efs
      vecm=ecm
      vecml = ecml
      vecmu = ecmu
      vecs = ecs
      vecsl = ecsl
      vecsu = ecsu
      vec95 = ec95
      vec95l = ec95l
      vec95u = ec95u
      vec99 = ec99
      vec99l = ec99l
      vec99u = ec99u
      vec995 = ec995
      vec995l = ec995l
      vec995u = ec995u
      vmsg1 = msg1
      vmsg2 = msg2
      vmsg3 = msg3
      
      vtcm = tdm
      vtcml = tdml
      vtcmu = tdmu
      vtcs = tds
      vtcns = tdns
      vtcsl = tdsl
      vtcsu = tdsu
      vtcxx = tdxx
      vtcxxl = tdxxl
      vtcxxu = tdxxu
 
      vtecm = tecm
      vtecs = tecs
      vte95 = te95
      vte95l = te95l
      vte95u = te95u
      vte99 = te99
      vte99l = te99l
      vte99u = te99u
      vte995 = te995
      vte995l = te995l
      vte995u = te995u
      vte999 = te999
      vte999l = te999l
      vte999u = te999u

      vudm = udm
      vuds = uds

      vphml = phml
      vphmu = phmu
      vcaml = caml
      vcamu = camu
      vudocml = udocml
      vudocmu = udocmu
      vedocml = edocml
      vedocmu = edocmu
      vphsl = phsl
      vphsu = phsu
      vcasl = casl
      vcasu = casu
      vudocsl = udocsl
      vudocsu = udocsu
      vedocsl = edocsl
      vedocsu = edocsu

      vddocm = ddocm
      vddocml = ddocml
      vddocmu = ddocmu
      vddocs = ddocs
      vddocsl = ddocsl
      vddocsu = ddocsu
      vddocns = ddocns

      vcofc1 = cofc1
      vcoff2 = coff2
      vcofc3 = cofc3
      vcocf4 = cocf4
      vcofc5 = cofc5
      vcocc6 = cocc6
      vucq = ucq ! aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

      vubam = ubam ! aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
      vubaml = ubaml ! aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
      vubamu = ubamu ! aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
      vubas = ubas
      vubasl = ubasl
      vubasu = ubasu
      vubaxx = ubaxx ! aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
      
      vtbm = tbm
      vtbml = tbml
      vtbmu = tbmu
      vtbs = tbs
      vtbsl = tbsl
      vtbsu = tbsu
      vtbxx = tbxx
      vtbxxl = tbxxl
      vtbxxu = tbxxu
      
      vubm = ubm
      vubml = ubml
      vubmu = ubmu
      vubs = ubs
      vubsl = ubsl
      vubsu = ubsu
      vubxx = ubxx
      vubxxl = ubxxl
      vubxxu = ubxxu
      
      vubmp = ubmp
      vvbmp = vbmp

      if ( sens .eq. 1 ) then !  perform the sensitivity tests - backward
      call sens1
      endif ! if ( sens .eq. 1 ) perform the sensitivity tests - backward
      if ( sens .eq. 2 ) then !  perform the sensitivity tests - forward
      call sens2
      endif ! if ( sens .eq. 2 ) perform the sensitivity tests - forward
      
      if ( sensa .eq. 1 ) then
        call sensa1
      endif
      if ( sensa .eq. 2 ) then
          call sensa2
      endif ! aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

      vx1 = abs(x1)
      vx2 = abs(x2)
      vx3 = abs(x3)
      vx4 = abs(x4)
      vx5 = abs(x5)
      vx6 = abs(x6)
      vx7 = abs(x7)
      vx8 = abs(x8)
      vx9 = abs(x9)
      vx10 = abs(x10)
      vx11 = abs(x11)
      vx12 = abs(x11)
      vx13 = abs(x13)
      vx14 = abs(x14)
      vx15 = abs(x15)

      vxa1 = xa1 ! aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
      vxa2 = xa2
      vxa3 = xa3
      vxa4 = xa4
      vxa5 = xa5
      vxa6 = xa6
      vxa7 = xa7
      vxa8 = xa8
      vxa9 = xa9
      vxa10 = xa10
      vxa11 = xa11
      vxa12 = xa12
      vxa13 = xa13
      vxa14 = xa14
      vxa15 = xa15 ! aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
      vxa16 = xa16 ! aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
      vxa17 = xa17 ! aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

      close(03)
      close(06,status = "delete")
      close(09,status = "delete")
      close(07,status = "delete")
      close(17,status = "delete")
      close(27,status = "delete")

      return
      end



      
*     perform the calculation
      subroutine do the calculations (isens)
      include 'bacom1cp.for'

*     set discharge quality
      tecm = ecm ! current mean discharge quality
      tecs = ecs ! current standard deviation
      call convert correlation coefficients  ! set up correlation
      call compute logged summary statistics ! and initialise arrays

*     store the starting variables to allow resetting after repeat runs
*     which request a temporary change to the data 
      if ( newrun .eq. 1 .and. isens .eq. 0 ) then
          call store the starting data 
      endif

*     set the iteration counter for the backward calculation
      ITER = 0
    
*     compute the mean and standard deviation for logged variables
      call LogFCf

*     test for and execute the forward calculation
      if (forw .eq. 1) then ! ==================================================
      if (ecm .lt. Small) then
      gecm = -Big
      gecs = 0.0
      else
      gecm = ALOG( (ecm*ecm)/SQRT(ecm*ecm + ecs*ecs) )
      gecs = SQRT( ALOG(1. + (ecs*ecs)/(ecm*ecm)) )
      endif

      if ( isens .eq. 0 ) call massb setup ! forward
      call massb setup two ! forward using the BMs

      call massb (4) ! perform the mass balance ... for the forward calculation 
      goto 995
      endif ! if ( forw .eq. 1 ) do the forward calculation ====================

*     prepare for a backward calculation ---------------------------------------

  996 ITER = ITER + 1 ! increase the iteration counter - start next iteration
      
      if ( ITER .eq. 1 ) then
          if ( isens .eq. 0 ) call massb setup ! backward
          call massb setup two ! backward using the BMs
      endif
      call massb (3) ! perform the mass balance (backward)
      if ( iter .eq. 1 ) then
          tbmc = tbm
          tbsc = tbs
      endif

*     obtain resulting mean or percentile of downstream river quality ...
      if ( free .eq. 1 )then
      targp = C(11,KTG)
      targm = tbm
      else
      targp = C(1,KTG)
      targm = tcm
      endif
     
*     compare the target with the downstream river quality =====================
*     check whether the target is achievable =================================== 
      if ( IDIL .eq. 1 ) then ! ================================================
*     zero discharge quality was tried on the last iteration ===================
*     did this indicate that the target was achieveable? =======================
      if ( type .ne. 1 ) then
      if ( targp .gt. 1.0001 * targ ) goto 995
      else
      if ( targm .gt. 1.0001 * targ ) goto 995
      endif
*     it looks, after all, as though the target is achievable ==================
      IDIL = 0 ! reset the indicator ... =======================================
      endif ! if ( IDIL .eq. 1 ) ===============================================

*     test for convergence ( that the correct discharge quality has been found) 
      if ( type .ne. 1 ) then
      CONV = ABS( ( targp - targ) / targ )
      else
      CONV = ABS( ( targm - targ) / targ )
      endif 

*     finish if convergence is achieved or if the number of iterations exceeds 50 ...
      if (CONV .lt. CFAC .OR. ITER .gt. 50) goto 995

*     set up the data for the next iteration
      if ( free .eq. 1 ) then ! Select on un-ionised or dissolved metal ...
      call nextf (ichek)
      else
      call next (ichek)
      endif
      
      if (ichek .eq. 0) goto 996 ! proceed to the next iteration

  995 continue

      rcq = C(3,KTG)
      tcxx = C(1,KTG)
      te80 = C(2,K80)
      te95 = C(2,K95)
      te99 = C(2,K99)
      te995 = C(2,K995)
      te999 = c(2,K999)
      
      call CONAV1 (rcm,rcs,rcns,rcml,rcmu)
      !rcml = rcml! - udm + rdm
      !rcmu = rcmu! - udm + rdm
      call CONSDEV1 (rcs,rcns,rcsl,rcsu)
      !rcsl = rcs + rcsl - uds
      !rcsu = rcs + rcsu - uds
      call LNCL1 (rcm,rcs,rcns,xper,rcqxxx,rcql,rcqu,0) 
      rcql = rcq + rcql - rcqxxx
      rcqu = rcq + rcqu - rcqxxx
      
      call CONAV1 (ecm,ecs,ecns,ecml,ecmu)
      call CONSDEV1 (ecs,ecns,ecsl,ecsu)
      call LNCL1 (ecm,ecs,ecns,95.0,ec95,ec95l,ec95u,0) 
      call LNCL1 (ecm,ecs,ecns,99.0,ec99,ec99l,ec99u,0) 
      call LNCL1 (ecm,ecs,ecns,99.5,ec995,ec995l,ec995u,0) 
      
      call CONAV1 (edocm,edocs,edocns,edocml,edocmu)
      call CONSDEV1 (edocs,edocns,edocsl,edocsu)
      !call LNCL1 (ecm,ecs,ecns,95.0,ec95,ec95l,ec95u,0) 
      !call LNCL1 (ecm,ecs,ecns,99.0,ec99,ec99l,ec99u,0) 
      !call LNCL1 (ecm,ecs,ecns,99.5,ec995,ec995l,ec995u,0) 

      call CONAV1 (ddocm,ddocs,ddocns,ddocml,ddocmu)
      call CONSDEV1 (ddocs,ddocns,ddocsl,ddocsu)
      call LNCL1 (ddocm,ddocs,ddocns,xper,ddocxxx,ddocxxl,ddocxxu,0) 
      ddocxxl = ddocxx + ddocxxl - ddocxxx
      ddocxxu = ddocxx + ddocxxu - ddocxxx

      xubns = rcns
      xubm = ubm
      xubs = ubs
      xxubs = ubs
      if ( nerr .eq. 2 ) xxubs = ubsplus
      call CONAV1 (xubm,xxubs,xubns,ubml,ubmu)
      
      !write(03,7181)nerr,ubs,xubs,ubsplus,xxubs,ubml,ubmu
 7181 format(/50('*')/
     &'   nerr =',i8/
     &'    ubs =',f8.4/
     &'   xubs =',f8.4/
     &'ubsplus =',f8.4/
     &'  xxubs =',f8.4/
     &'   ubml =',f8.4/
     &'   ubmu =',f8.4/50('*')/)

      call CONSDEV2 (xubs,xxubs,xubns,ubsl,ubsu)
      call LNCL1 (xubm,xxubs,xubns,95.0,ub95,ub95l,ub95u,0) 

      !write(03,7185)ub95,ub95l,ub95u
 7185 format(/50('*')/
     &'     ub95 =',f8.4/
     &'    ub95l =',f8.4/
     &'    ub95u =',f8.4/
     &50('*'))

      ub95l = amax1(0.0, C(12,K95) + ub95l - ub95 )
      ub95u = C(12,K95) + ub95u - ub95 
      ub95 = C(12,K95)

      !write(03,7186)C(12,K95),ub95l,ub95u
 7186 format( 
     &'C(12,K95) =',f8.4/
     &'    ub95l =',f8.4/
     &'    ub95u =',f8.4/
     &50('*')/)

*     ==========================================================================
      xvbns = rcns
      xvbm = vbm
      xvbs = vbs
      xxvbs = vbs
      if ( nerr .eq. 2 ) xxvbs = vbsplus
      call CONAV1 (xvbm,xxvbs,xvbns,vbml,vbmu)
      
      !write(03,7581)nerr,vbs,xvbs,vbsplus,xxvbs,vbml,vbmu
 7581 format(/50('*')/
     &'   nerr =',i8/
     &'    vbs =',f8.4/
     &'   xvbs =',f8.4/
     &'vbsplus =',f8.4/
     &'  xxvbs =',f8.4/
     &'   vbml =',f8.4/
     &'   vbmu =',f8.4/50('*')/)

      call CONSDEV2 (xvbs,xxvbs,xvbns,vbsl,vbsu)
      call LNCL1 (xvbm,xxvbs,xvbns,95.0,vb95,vb95l,vb95u,0) 

      !write(03,7585)vb95,vb95l,vb95u
 7585 format(/50('*')/
     &'     vb95 =',f8.4/
     &'    vb95l =',f8.4/
     &'    vb95u =',f8.4/
     &50('*'))

      vb95l = amax1(0.0, C(16,K95) + vb95l - vb95 )
      vb95u = C(16,K95) + vb95u - vb95 
      vb95 = C(16,K95)

      !write(03,7586)C(16,K95),vb95l,vb95u
 7586 format( 
     &'C(16,K95) =',f8.4/
     &'    vb95l =',f8.4/
     &'    vb95u =',f8.4/
     &50('*')/)
*     ==========================================================================

      xtdns = tdns
      xtdns = C(14,NS+1)

      if ( forw .eq. 1 )xtdns = C(15,NS+1)
      tdns = xtdns
      xtdm = tdm
      xtds = tds
      call CONAV1 (xtdm,xtds,xtdns,tdml,tdmu)
      call CONSDEV1 (xtds,xtdns,tdsl,tdsu)

      !xtbns = tbns
      xtbns = C(14,NS+1)
      if ( forw .eq. 1 ) xtbns = C(15,NS+1)
      tbns = xtbns
      xtbm = tbm
      xtbs = tbs
      xxtbs = xtbs
      if ( nerr .eq. 2 ) xxtbs = tbsplus * xtbs
      call CONAV1 (xtbm,xxtbs,xtbns,tbml,tbmu)

      !write(3,5511)tbs,tbsplus,xxtbs,xtbns,tbml,tbmu
 5511 format(//50('-')/
     &' original sd (tbs) =',f10.4/
     &'           tbsplus =',f10.4/
     &'      new SD xxtbs =',f10.4/
     &'    samples(xtbns) =',f10.4/50('-')/
     &'             upper =',f10.4/
     &'             lower =',f10.4/50('-'))
      
      call CONSDEV2 (xtbs,xxtbs,xtbns,tbsl,tbsu)
      
      tecsX = tecs
      !if ( nerr .eq. 2 ) write(3,5521)tecsX,tecsplus,tecs,tecm
 5521 format(//50('-')/
     &'    tecsX =',f10.4/
     &' tecsplus =',f10.4/
     &'     tecs =',f10.4/
     &'     tecm =',f10.4/
     &50('-')/)
      if ( nerr .eq. 2 ) tecsX = tecsplus
    
      call CONAV1 (tecm,tecsX,ecns,tecml,tecmu)
      call CONSDEV2 (tecs,tecsX,ecns,tecsl,tecsu)
      
      te80xxx = 0.0
      te95xxx = 0.0
      te99xxx = 0.0
      te995xxx = 0.0
      te999xxx = 0.0
      
      te80l = 0.0
      te95l = 0.0
      te99l = 0.0
      te995l = 0.0
      te999l = 0.0
     
      if ( tecm .gt. 0.00001 ) then
      call LNCL1 (tecm,tecsX,ecns,80.0,te80xxx,te80l,te80u,0) 
      te80l = amax1(0.0,te80 + te80l - te80xxx)
      te80u = te80 + te80u - te80xxx
      call LNCL1 (tecm,tecsX,ecns,95.0,te95xxx,te95l,te95u,0) 
      te95l = amax1(0.0,te95 + te95l - te95xxx)
      te95u = te95 + te95u - te95xxx
      call LNCL1 (tecm,tecsX,ecns,99.0,te99xxx,te99l,te99u,0) 
      te99l = amax1(0.0,te99 + te99l - te99xxx)
      te99u = te99 + te99u - te99xxx
      call LNCL1 (tecm,tecsX,ecns,99.5,te995xxx,te995l,te995u,0) 
      te995l = amax1(0.0,te995 + te995l - te995xxx)
      te995u = te995 + te995u - te995xxx
      call LNCL1 (tecm,tecsX,ecns,99.9,te999xxx,te999l,te999u,0) 
      te999l = amax1(0.0,te999 + te999l - te999xxx)
      te999u = te999 + te999u - te999xxx
      endif

      call CONAV1 (phm,phs,phns,phml,phmu)
      call CONSDEV1 (phs,phns,phsl,phsu)
      call CONAV1 (cam,cas,cans,caml,camu)
      call CONSDEV1 (cas,cans,casl,casu)
      call CONAV1 (udocm,udocs,udocns,udocml,udocmu)
      call CONSDEV1 (udocs,udocns,udocsl,udocsu)

      if ( isens .eq. 0 ) then
      if ( nerr .eq. 0 .or. nerr .eq. 2 ) then
      call write out the input data
      call results ! write out the results
      call write downstream bioavailable metal (4,98) ! 222222222222222222222222
      endif
      endif

*     check for convergence failure      
      if ( IDIL .eq. 0 .and. CONV .ge. CFAC ) msg2 = 1

*     check for failure to meet the target ...
      if ( type .ne. 1 ) then
      if ( tcxx .gt. 1.0001 * targ ) msg3 = 1
      else
      if ( tcm .gt. 1.0001 * targ ) msg3 = 1
      endif

      return
      end




*     perform Mass Balance ...
*     calculate river quality downstream of discharge ...

      subroutine massb setup
      include 'bacom1cp.for'

*     accumulators for mean and standard deviations
      crfm = 0.0
      crf5 = 0.0
      crcm = 0.0
      crcs = 0.0
      cefm = 0.0
      cefs = 0.0
      cecm = 0.0
      cecs = 0.0
      cphm = 0.0
      cphs = 0.0
      ccam = 0.0
      ccas = 0.0
      cudocm = 0.0
      cudocs = 0.0 
      cedocm = 0.0 
      cedocs = 0.0  

*     set the random number starters -------------------------------------------
      IR1 = JR1 ! river flow
      IR3 = JR3 ! discharge flow
      IR2 = JR2 ! river quality for dissolved metal 
      IR4 = JR4 ! discharge quality for dissolved metal
      IR5 = JR5 ! discharge quality for dissolved organic carbon 
      IR6 = JR6 ! dissolved organic carbon in the upstream river
      IR7 = JR7 ! pH in downstream river
      IR8 = JR8 ! calcium in downstream river
      IR9 = JR9 ! goodness-of-fit in the bioavailability equation
      
      do 2 I = 1,NS ! loop around the Monte Carlo shots 
          
      call get a set of correlated random normal deviates
      call get the shot for river flow (I,RF)
      call get the shot for river quality (I,RC)
      call get the shot for discharge flow (I,EF)
      call get the shot for discharge quality (I,EC)
      call get the shot for upstream DOC (I,RD)
      call get the shot for discharge DOC (I,ED)
      call get the shot for downstream pH (I,PH)
      call get the shot for downstream calcium (I,CA)

      crfm = crfm + RF
      TMC(1,I) = RF
      crcm = crcm + RC
      crcs = crcs + RC * RC
      cefm = cefm + EF
      cefs = cefs + EF * EF
      cecm = cecm + EC
      cecs = cecs + EC * EC
      cphm = cphm + PH
      cphs = cphs + PH * PH
      ccam = ccam + CA
      ccas = ccas + CA * CA
      cudocm = cudocm + RD
      cudocs = cudocs + RD * RD
      cedocm = cedocm + ED
      cedocs = cedocs + ED * ED
      
    2 continue

*     order shots for upsteam river flow
      K = 1
      do 8 I = 1, NS-1
      do 11 J = I+1, NS
      if (TMC(K,I) .lt. TMC(K,J)) go to 11
      CC = TMC(K,I)
      TMC(K,I) = TMC(K,J)
      TMC(K,J) = CC
   11 continue
    8 continue

*     compute summary statistics for checking errors in Monte Carlo simulation
      crfm = crfm / FLOAT(NS)
      crf5 = TMC(1,k05)

      XXX = crcs - crcm*crcm/FLOAT(NS)
      if (XXX .le. Small) then 
      crcs = 0.0
      else
      crcs = SQRT (XXX/FLOAT(NS-1))
      endif
      crcm = crcm / FLOAT(NS)
      XXX = cefs - cefm*cefm/FLOAT(NS)
      if (XXX .le. Small) then 
      cefs = 0.0
      else
      cefs = SQRT (XXX/FLOAT(NS-1))
      endif
      cefm = cefm / FLOAT(NS)
      XXX = cecs - cecm*cecm/FLOAT(NS)
      if (XXX .le. Small) then 
      cecs = 0.0
      else
      cecs = SQRT (XXX/FLOAT(NS-1))
      endif
      cecm = cecm / FLOAT(NS)
      
      XXX = cphs - cphm*cphm/FLOAT(NS)
      if (XXX .le. Small) then 
      cphs = 0.0
      else
      cphs = SQRT (XXX/FLOAT(NS-1))
      endif
      cphm = cphm / FLOAT(NS)

      XXX = ccas - ccam*ccam/FLOAT(NS)
      if (XXX .le. Small) then 
      ccas = 0.0
      else
      ccas = SQRT (XXX/FLOAT(NS-1))
      endif
      ccam = ccam / FLOAT(NS)

      XXX = cudocs - cudocm*cudocm/FLOAT(NS)
      if (XXX .le. Small) then 
      cudocs = 0.0
      else
      cudocs = SQRT (XXX/FLOAT(NS-1))
      endif
      cudocm = cudocm / FLOAT(NS)

      XXX = cedocs - cedocm*cedocm/FLOAT(NS)
      if (XXX .le. Small) then 
      cedocs = 0.0
      else
      cedocs = SQRT (XXX/FLOAT(NS-1))
      endif
      cedocm = cedocm / FLOAT(NS)

      xx1 = 100.0 * (crfm/xrfm - 1.0)
      xx2 = 100.0 * (crcm/xrcm - 1.0)
      xx3 = 100.0 * (cefm/xefm - 1.0)
      xx4 = 100.0 * (cecm/xecm - 1.0)
      xx5 = 100.0 * (cphm/xphm - 1.0)
      xx6 = 100.0 * (ccam/xcam - 1.0)
      xx7 = 100.0 * (cudocm/xudocm - 1.0)
      xx8 = 100.0 * (cedocm/xedocm - 1.0)
      
      xy1 = 100.0 * (crf5/xrf5 - 1.0)
      xy2 = 100.0 * (crcs/xrcs - 1.0)
      xy3 = 100.0 * (cefs/xefs - 1.0)
      xy4 = 100.0 * (cecs/xecs - 1.0)
      xy5 = 100.0 * (cphs/xphs - 1.0)
      xy6 = 100.0 * (ccas/xcas - 1.0)
      xy7 = 100.0 * (cudocs/xudocs - 1.0)
      xy8 = 100.0 * (cedocs/xedocs - 1.0)

      if ( nbias .eq. 0 ) then ! calculate corrections in massb setup ----------
      BM(1) = xrfm / crfm
      BM(2) = xrcm / crcm
      BM(3) = xefm / cefm
      if ( forw .eq. 1 ) then
      BM(4) = xecm / cecm
      else
      BM(4) = tecm / cecm
      endif
      BM(5) = xphm / cphm
      BM(6) = xcam / ccam
      BM(7) = xudocm / cudocm
      BM(8) = xedocm / cedocm
      
      BS(1) = xrf5 / crf5
      BS(2) = xrcs / crcs
      BS(3) = xefs / cefs
      if ( forw .eq. 1 ) then
      BS(4) = xecs / cecs
      else
      BS(4) = tecs / cecs
      endif
      BS(5) = xphs / cphs
      BS(6) = xcas / ccas
      BS(7) = xudocs / cudocs
      BS(8) = xedocs / cedocs
      endif ! calculate corrections in massb setup -----------------------------
      
      write(9,7000)iss
 7000 format(/95('=')/'Checking the accuracy of Monte Carlo ',
     &'simulation'i6,'...............'/95('=')/
     &'     calc','     orig','     diff',6x,
     &'     calc','     orig','     diff',5x,'    BM','    BS'/
     &95('='))
      write(9,8000)crfm,xrfm,xx1,crf5,xrf5,xy1,BM(1),BS(1)
 8000 format(3f9.4,' %    ',3f9.4,' %',5x,2f6.3,'   RF')
      write(9,8001)crcm,xrcm,xx2,crcs,xrcs,xy2,BM(2),BS(2)
 8001 format(3f9.4,' %    ',3f9.4,' %',5x,2f6.3,'   RC')
      write(9,8002)cefm,xefm,xx3,cefs,xefs,xy3,BM(3),BS(3)
 8002 format(3f9.4,' %    ',3f9.4,' %',5x,2f6.3,'   EF')
      write(9,8003)cecm,xecm,xx4,cecs,xecs,xy4,BM(4),BS(4)
 8003 format(3f9.4,' %    ',3f9.4,' %',5x,2f6.3,'   EC')
      write(9,8004)cphm,xphm,xx5,cphs,xphs,xy5,BM(5),BS(5)
 8004 format(3f9.4,' %    ',3f9.4,' %',5x,2f6.3,'   PH')
      write(9,8005)ccam,xcam,xx6,ccas,xcas,xy6,BM(6),BS(6)
 8005 format(3f9.4,' %    ',3f9.4,' %',5x,2f6.3,'   CA')
      write(9,8006)cudocm,xudocm,xx7,cudocs,xudocs,xy7,BM(7),BS(7)
 8006 format(3f9.4,' %    ',3f9.4,' %',5x,2f6.3,'   UD')
      write(9,8007)cedocm,xedocm,xx8,cedocs,xedocs,xy8,BM(8),BS(8)
 8007 format(3f9.4,' %    ',3f9.4,' %',5x,2f6.3,'   ED')
      write(9,7001)
 7001 format(95('='))

      return 
      end

      
      
      
      subroutine massb setup two
      include 'bacom1cp.for'

*     accumulators for mean and standard deviations
      crfm = 0.0
      crf5 = 0.0
      crcm = 0.0
      crcs = 0.0
      cefm = 0.0
      cefs = 0.0
      cecm = 0.0
      cecs = 0.0
      cphm = 0.0
      cphs = 0.0
      ccam = 0.0
      ccas = 0.0
      cudocm = 0.0
      cudocs = 0.0 
      cedocm = 0.0 
      cedocs = 0.0  

*     set the random number starters -------------------------------------------
      IR1 = JR1 ! river flow
      IR3 = JR3 ! discharge flow
      IR2 = JR2 ! river quality for dissolved metal 
      IR4 = JR4 ! discharge quality for dissolved metal
      IR5 = JR5 ! discharge quality for dissolved organic carbon 
      IR6 = JR6 ! dissolved organic carbon in the upstream river
      IR7 = JR7 ! pH in downstream river
      IR8 = JR8 ! calcium in downstream river
      IR9 = JR9 ! goodness-of-fit in the bioavailability equation
      
      do 2 I = 1,NS ! loop around the Monte Carlo shots 
          
      call get a set of correlated random normal deviates
      call get the shot for river flow (I,RF)
      call get the shot for river quality (I,RC)
      call get the shot for discharge flow (I,EF)
      call get the shot for discharge quality (I,EC)
      call get the shot for upstream DOC (I,RD)
      call get the shot for discharge DOC (I,ED)
      call get the shot for downstream pH (I,PH)
      call get the shot for downstream calcium (I,CA)

      crfm = crfm + RF
      TMC(1,I) = RF
      crcm = crcm + RC
      crcs = crcs + RC * RC
      cefm = cefm + EF
      cefs = cefs + EF * EF
      cecm = cecm + EC
      cecs = cecs + EC * EC
      cphm = cphm + PH
      cphs = cphs + PH * PH
      ccam = ccam + CA
      ccas = ccas + CA * CA
      cudocm = cudocm + RD
      cudocs = cudocs + RD * RD
      cedocm = cedocm + ED
      cedocs = cedocs + ED * ED
      
    2 continue

*     order shots for upsteam river flow
      K = 1
      do 8 I = 1, NS-1
      do 11 J = I+1, NS
      if (TMC(K,I) .lt. TMC(K,J)) go to 11
      CC = TMC(K,I)
      TMC(K,I) = TMC(K,J)
      TMC(K,J) = CC
   11 continue
    8 continue

*     compute summary statistics for checking errors in Monte Carlo simulation
      crfm = crfm / FLOAT(NS)
      crf5 = TMC(1,k05)

      XXX = crcs - crcm*crcm/FLOAT(NS)
      if (XXX .le. Small) then 
      crcs = 0.0
      else
      crcs = SQRT (XXX/FLOAT(NS-1))
      endif
      crcm = crcm / FLOAT(NS)
      XXX = cefs - cefm*cefm/FLOAT(NS)
      if (XXX .le. Small) then 
      cefs = 0.0
      else
      cefs = SQRT (XXX/FLOAT(NS-1))
      endif
      cefm = cefm / FLOAT(NS)
      XXX = cecs - cecm*cecm/FLOAT(NS)
      if (XXX .le. Small) then 
      cecs = 0.0
      else
      cecs = SQRT (XXX/FLOAT(NS-1))
      endif
      cecm = cecm / FLOAT(NS)
      
      XXX = cphs - cphm*cphm/FLOAT(NS)
      if (XXX .le. Small) then 
      cphs = 0.0
      else
      cphs = SQRT (XXX/FLOAT(NS-1))
      endif
      cphm = cphm / FLOAT(NS)

      XXX = ccas - ccam*ccam/FLOAT(NS)
      if (XXX .le. Small) then 
      ccas = 0.0
      else
      ccas = SQRT (XXX/FLOAT(NS-1))
      endif
      ccam = ccam / FLOAT(NS)

      XXX = cudocs - cudocm*cudocm/FLOAT(NS)
      if (XXX .le. Small) then 
      cudocs = 0.0
      else
      cudocs = SQRT (XXX/FLOAT(NS-1))
      endif
      cudocm = cudocm / FLOAT(NS)

      XXX = cedocs - cedocm*cedocm/FLOAT(NS)
      if (XXX .le. Small) then 
      cedocs = 0.0
      else
      cedocs = SQRT (XXX/FLOAT(NS-1))
      endif
      cedocm = cedocm / FLOAT(NS)

      xx1 = 100.0 * (crfm/xrfm - 1.0)
      xx2 = 100.0 * (crcm/xrcm - 1.0)
      xx3 = 100.0 * (cefm/xefm - 1.0)
      xx4 = 100.0 * (cecm/xecm - 1.0)
      xx5 = 100.0 * (cphm/xphm - 1.0)
      xx6 = 100.0 * (ccam/xcam - 1.0)
      xx7 = 100.0 * (cudocm/xudocm - 1.0)
      xx8 = 100.0 * (cedocm/xedocm - 1.0)
      
      xy1 = 100.0 * (crf5/xrf5 - 1.0)
      xy2 = 100.0 * (crcs/xrcs - 1.0)
      xy3 = 100.0 * (cefs/xefs - 1.0)
      xy4 = 100.0 * (cecs/xecs - 1.0)
      xy5 = 100.0 * (cphs/xphs - 1.0)
      xy6 = 100.0 * (ccas/xcas - 1.0)
      xy7 = 100.0 * (cudocs/xudocs - 1.0)
      xy8 = 100.0 * (cedocs/xedocs - 1.0)
      
      write(9,7000)iss
 7000 format(/95('-')/'Checking the accuracy of Monte Carlo ',
     &'simulation',i6,' ...............'/95('-')/
     &'     calc','     orig','     diff',6x,
     &'     calc','     orig','     diff',5x,'    BM'/
     &95('-'))
      write(9,8000)crfm,xrfm,xx1,crf5,xrf5,xy1,BM(1)
 8000 format(3f9.4,' %    ',3f9.4,' %',5x,f6.3,'   RF')
      write(9,8001)crcm,xrcm,xx2,crcs,xrcs,xy2,BM(2)
 8001 format(3f9.4,' %    ',3f9.4,' %',5x,f6.3,'   RC')
      write(9,8002)cefm,xefm,xx3,cefs,xefs,xy3,BM(3)
 8002 format(3f9.4,' %    ',3f9.4,' %',5x,f6.3,'   EF')
      write(9,8003)cecm,xecm,xx4,cecs,xecs,xy4,BM(4)
 8003 format(3f9.4,' %    ',3f9.4,' %',5x,f6.3,'   EC')
      write(9,8004)cphm,xphm,xx5,cphs,xphs,xy5,BM(5)
 8004 format(3f9.4,' %    ',3f9.4,' %',5x,f6.3,'   PH')
      write(9,8005)ccam,xcam,xx6,ccas,xcas,xy6,BM(6)
 8005 format(3f9.4,' %    ',3f9.4,' %',5x,f6.3,'   CA')
      write(9,8006)cudocm,xudocm,xx7,cudocs,xudocs,xy7,BM(7)
 8006 format(3f9.4,' %    ',3f9.4,' %',5x,f6.3,'   UD')
      write(9,8007)cedocm,xedocm,xx8,cedocs,xedocs,xy8,BM(8)
 8007 format(3f9.4,' %    ',3f9.4,' %',5x,f6.3,'   ED')
      write(9,7001)
 7001 format(95('-'))

      return 
      end

   


      subroutine write out the input data
      include 'bacom1cp.for'

*     upstream river flow

      if ( forw .eq. 1 ) write(3,581)
  581 format(68('-')/'Calculation of the river quality d/s ',
     &'of the input discharge quality')

      if ( rfm .gt. Small ) then
      write(3,2)rfm,rf5
    2 format(68('-')/'Input data ... '/68('-')/
     &'Mean river flow upstream of discharge',F15.2/
     &'95-percent exceedence river flow     ',F15.2)
      endif

*     third parameter for Log-normal distribution
      if ( fsh .ne. 0.0 ) then
      write(3,8)fsh
    8 format(
     &'Shift parameter                      ',F15.2)
      endif

*     upstream river quality
      write(3,3)rcm,rcml,rcmu
    3 format(68('-')/
     &'Mean upstream dissolved metal',F23.2,2f8.2)
      if ( rcm .gt. Small ) then
      write(3,942)rcs,rcsl,rcsu,int(rcns)
  942 format(
     &'Standard deviation              ',F20.2,2f8.2/
     &'Number of samples               ',i20)
      endif

*     constant addition of river flow upstream of the discharge
      if ( add .ne. 0.0 ) then 
      write(3,6)add
    6 format(68('-')/
     &'Constant addition to river flow ',F20.2)
      endif

*     regulated flow upstream of the discharge
      if ( reg .ne. 0.0 ) then
      write(3,7)reg
    7 format(68('-')/
     &'Maintained river flow           ',F20.2)
      endif

*     quality of added water or regulation water
      if ( iregq .eq. 1 ) then
      write(3,10)rgm
   10 format(68('-')/
     &'Mean quality of added water     ',F20.2)
      write(3,11)rgs
   11 format(
     &'Corresponding standard deviation',F20.2)
      endif
      
      if ( free .eq. 1 ) then 
      write(3,3549)udocm,udocml,udocmu,udocs,udocsl,udocsu,int(udocns)
 3549 format(
     &'Mean upstream dissolved organic carbon',F14.2,2f8.2/
     &'Standard deviation                    ',F14.2,2f8.2/
     &'Number of samples                     ',i14)
*     pH downstream of the discharge
      !write(3,9349)phm,phml,phmu,phs,phsl,phsu,int(phns)
 9349 format(68('-')/
     &'Mean pH downstream of discharge       ',F14.2,2f8.2/
     &'Standard deviation                    ',F14.2,2f8.2/
     &'Number of samples                     ',i14)
*     calcium, Total Dissolved Solids and dissolved organic carbon
      !write(3,3349)cam,caml,camu,cas,casl,casu,int(cans)
 3349 format(
     &'Mean calcium downstream of discharge  ',F14.2,2f8.2/
     &'Standard deviation                    ',F14.2,2f8.2/
     &'Number of samples                     ',i14)
      endif

      return
      end




      subroutine results ! write out the results
      include 'bacom1cp.for'

      tcxx = c(1,KTG)
      tcx95 = c(1,K95)
      tcx99 = c(1,K99)

      xxx = float (int(xper)) - xper
      kkk = 0
      if ( xxx .lt. Small .and. xxx .gt. -Small ) kkk = int(xper)

*     check for illegal shift flow
*     if (ifsh .eq. 1) write(3,2)
    2 format(68('-')/
     &'###  The shift flow is too large a negative value ...'/
     &'###  The value has been reset to -0.99 of the '/
     &'###  95% exceedence river flow ...')
*     upstream river quality (if regulation used)
*     if (iregq .eq. 1) then
      write(3,958)udm,uds
  958 format(68('-')/
*    &'Effect of the added river flows ...'/80('-')/
     &'Mean upstream dissolved metal         ',F14.2/
     &'Standard deviation                    ',F14.2)
*     endif

      if ( type .eq. 0 ) then 
      if ( kkk .eq. 0 .and. rcm .gt. Small ) then
      write(3,819)xper,rcq,rcql,rcqu
  819 format(f4.1,'-percentile river quality',F21.2,2f8.2)
      else
      write(3,829)int(xper),rcq,rcql,rcqu
  829 format(i2,'-percentile river quality  ',F23.2,2f8.2)
      endif
      endif

      write(3,943)efm,efs ! discharge flow
  943 format(68('-')/
     &'Mean flow of discharge                ',F14.2/
     &'Standard deviation                    ',F14.2)

      write(3,974)ecm,ecml,ecmu,ecs,ecsl,ecsu,int(ecns) ! discharge quality
  974 format(68('-')/
     &'Mean quality of discharge             ',F14.2,2f8.2/
     &'Standard deviation                    ',F14.2,2f8.2/
     &'Number of samples                     ',i14)

      if ( free .eq. 1 ) then 
      write(3,944)edocm,edocml,edocmu,edocs,edocsl,edocsu,int(edocns) ! DOC in discharge quality
  944 format(
     &'Mean dissolved organic carbon in discharge',F10.2,2f8.2/
     &'Standard deviation                        ',F10.2,2f8.2/
     &'Number of samples                         ',i10)
*     pH downstream of the discharge
      write(3,9349)phm,phml,phmu,phs,phsl,phsu,int(phns)
 9349 format(68('-')/
     &'Mean pH downstream of discharge       ',F14.2,2f8.2/
     &'Standard deviation                    ',F14.2,2f8.2/
     &'Number of samples                     ',i14)
*     calcium downstream of the discharge
      write(3,3349)cam,caml,camu,cas,casl,casu,int(cans)
 3349 format(
     &'Mean calcium downstream of discharge  ',F14.2,2f8.2/
     &'Standard deviation                    ',F14.2,2f8.2/
     &'Number of samples                     ',i14)
      endif

      if ( free .eq. 1 ) then 
      write(3,3889)ubm,ubml,ubmu,ubs,ubsl,ubsu
 3889 format(68('-')/
     &'Mean upstream bioavailable metal      ',F14.4,2f8.4/
     &'Standard deviation                    ',F14.4,2f8.4)
      if ( type .eq. 1 ) then
      write(3,3882)ub95,ub95l,ub95u
 3882 format(
     &'95-percentile                         ',F14.4,2f8.4)
      else
      write(3,3819)int(xper),ubxx,ubxxl,ubxxu
 3819 format(
     &i2,'-percentile                         ',F14.4,2f8.4)
      endif
      endif

      if ( free .eq. 1 ) then 
      vbmp = 0.0
      ubmp = 0.0
      if ( targ .gt. 0.000001 ) then
      vbmp = 100.0 * (tbmc - vbm) / targ
      ubmp = 100.0 * (tbmc - ubm) / targ
      vbmp2 = 100.0 * (tbm - vbm) / targ
      ubmp2 = 100.0 * (tbm - ubm) / targ
      write(3,3489)vbm,tbmc,vbmp,ubmp,vbmp2,ubmp2 !,vbml,vbmu,vbs,vbsl,vbsu
 3489 format(68('=')/
     &'u/s bioavailable metal (with d/s DOC)  ',F13.4/
     &'d/s of current discharge (with d/s DOC)',F13.4/68('-')/
     &'(u/s - d/s) gap for current discharge ....'/68('-')/
     &'% of target (with d/s DOC for u/s metal)',F12.1,' %'/
     &'% of target (with u/s DOC for u/s metal)',F12.1,' %'/68('-')/
     &'(u/s - d/s) gap for modified discharge ....'/68('-')/
     &'% of target (with d/s DOC for u/s metal)',F12.1,' %'/
     &'% of target (with u/s DOC for u/s metal)',F12.1,' %')
      else
      write(3,3689)vbm,tbmc !,vbml,vbmu,vbs,vbsl,vbsu
 3689 format(68('=')/'u/s bioavailable metal (with ',
     &'d/s DOC) ',F14.4/
     &'d/s of current discharge (with d/s DOC)',F13.4)
      endif
      
    ! &'Mean upstream bioavailable metal      ',F14.4,2f8.4/
    !&'Standard deviation                    ',F14.4,2f8.4)
      if ( type .eq. 1 ) then
      !write(3,3482)ub95,vb95l,vb95u
 3482 format(
     &'95-percentile                         ',F14.4,2f8.4)
      else
      !write(3,3419)int(xper),vbxx,vbxxl,vbxxu
 3419 format(
     &i2,'-percentile                         ',F14.4,2f8.4)
      endif
      
      write(3,3559)ddocm,ddocml,ddocmu,ddocs,ddocsl,ddocsu,ddocns
 3559 format(68('=')/
     &'Mean downstream dissolved organic carbon',F12.2,2f8.2/
     &'Standard deviation                      ',F12.2,2f8.2/
     &'Number of samples                       ',f12.1)
      endif

      if ( free .eq. 0 ) then
      if ( forw .eq. 0 ) then ! ================================================
      if ( type .eq. 0 ) then ! target is a percentile -------------------------
      if ( kkk .eq. 0 ) then ! percentile point is not an integer --------------    
      write(3,125) xper,targ! River target (dissolved metal)
  125 format(68('-')/'River Quality Target (',F4.1,
     &'-percentile)',F14.2/68('-')/'Results ...')
      else ! if ( kkk .eq. 0 ) then
      write(3,425) kkk,targ! River target (dissolved metal)
  425 format(68('-')/
     &'River Quality Target (',i2,'-percentile)',f16.2/
     &68('-')/'Results ...')
      endif ! if ( kkk .eq. 0 ) then
      else ! then target is a maen ---------------------------------------------
      write(3,525) targ ! River target (dissolved metal)
  525 format(68('-')/10x,'River Quality Target (mean)',f24.2/
     &68('-')/'Results ...')
      endif ! if ( type .eq. 0 ) then
      endif ! if (forw .eq. 0) then ============================================

      write(3,908)tdm,tdml,tdmu,tds,tdsl,tdsu,td95,td95lc,td95uc,tdns 
  908 format(68('-')/
     &'Mean dissolved metal d/s of discharge   ',F12.2,2f8.2/
     &'Standard deviation                      ',F12.2,2f8.2/
     &'95-percentile dissolved metal           ',F12.2,2f8.2/
     &'Number of samples                       ',f12.1)

      else ! if (free .eq. 0) then
*     river target (bio-available metal)
      if ( type .eq. 0 ) then
      if ( kkk .eq. 0 ) then
      write(3,115) xper,targ
  115 format(68('-')/'River Quality Target (',F4.1,
     &'-percentile)',F12.4/68('-')/'Results ...')
      else ! if ( kkk .eq. 0 ) then
      write(3,415) kkk,targ
  415 format(68('-')/'River Quality Target (',i2,
     &'-percentile)',F16.4/68('-')/'Results ...')
      endif ! if ( kkk .eq. 0 ) then
      else ! if ( type .eq. 0 ) then
      write(3,515)targ
  515 format(68('-')/'River Quality Target (mean)',F25.4/
     &68('-')/'Results ...')
      endif ! if ( type .eq. 0 ) 
      endif !  ! if (free .eq. 0) 

*     bio-available metal downstream of the discharge
      if ( biometal .gt. 1 ) then ! ============================================
      write(3,908)tdm,tdml,tdmu,tds,tdsl,tdsu,td95,td95lc,td95uc,tdns 
      write(3,966)tbm,tbml,tbmu,tbs,tbsl,tbsu
  966 format(68('-')/
     &'Mean bioavailable metal d/s of discharge ',f11.4,2f8.4/
     &'Standard deviation                       ',f11.4,2f8.4)
      
      write(3,4969)C(11,K95),tb95lc,tb95uc,tbns
 4969 format(
     &'95-percentile                            ',f11.4,2f8.4/
     &'Number of samples                        ',f11.1)
      if ( type .eq. 0 ) then
      if ( kkk .eq. 0 ) then
      write(3,3939)xper,C(11,KTG),tbxxl,tbxxu
 3939 format(
     &f4.1,'percentile bio-available metal    ',F11.4,2f8.4)
      else
      write(3,4939)kkk,C(11,KTG),tbxxl,tbxxu
 4939 format(
     &i2,'-percentile bio-available metal   ',F16.4,2f8.4)
      endif
      endif ! if ( type .eq. 0 ) 
      endif ! if ( biometal .gt. 1 ) ===========================================
      
      if (efm .lt. Small) return

*     discharge quality needed to meet river target
      write(3,978)tecm,tecml,tecmu,tecs,tecsl,tecsu   
  978 format(68('-')/
     &'Mean dissolved metal in discharge       ',F12.2,2f8.2/
     &'Standard deviation                      ',F12.2,2f8.2)
      if ( nsamp .gt. 0 ) write(7,448)tecm,tecml,tecmu   
  448 format(
     &'Mean dissolved metal in discharge       ',F12.2,2f8.2)

*     discharge quality standards 
      if (C(2,K95) .gt. Small) then 
      write(3,168)C(2,K80),te80l,te80u
  168 format('80-percentile                    ',F19.2,2f8.2)
      write(3,968)C(2,K95),te95l,te95u
  968 format('95-percentile                    ',F19.2,2f8.2)
      write(3,948)C(2,K99),te99l,te99u
  948 format('99-percentile                    ',F19.2,2f8.2)
      write(3,868)C(2,K995),te995l,te995u
  868 format('99.5-percentile                  ',F19.2,2f8.2)
      write(3,768)C(2,K999),te999l,te999u,int(ecns)
  768 format('99.9-percentile                  ',F19.2,2f8.2/
     &       'Number of samples                ',i19)
      endif
      
      if ( nsampsim .eq. 0 ) then
      write(3,994)
  994 format(68('-'))
      write(3,822)coff2,cofc1,cofc5
  822 format('CORRELATION'/68('-')/
     &'River and discharge flow:       ',f20.4/
     &'River flow and river quality:   ',f20.4/
     &'Discharge flow and quality      ',f20.4)
      if ( free .eq. 1 ) then
      write(3,825)capf,cacf,cacp,cafd,caed
  825 format(
     &'River flow and pH:              ',f20.4/
     &'River flow and calcium:         ',f20.4/
     &'Calcium and pH:                 ',f20.4/
     &'River DOC and river flow:       ',f20.4/
     &'Discharge DOC and discharge flow:',f19.4)
      endif
      write(3,994)
      endif

      return
      end



*     Arrange high shots in increasing order ...
*     In order to obtain the percentiles ...

      subroutine sequence
      include 'bacom1cp.for'

      KL = min0 ( KTG, K90 )

      K = 3
      do 3 I = NS, KL, -1
      do 6 J = 1, I - 1
      if (C(K,I) .gt. C(K,J)) goto 6
      CC = C(K,I)
      C(K,I) = C(K,J)
      C(K,J) = CC
    6 continue
    3 continue

*     order river quality shots in downstream river   
      K = 1
      do 4 I = NS, KL, -1
      do 5 J = 1, I - 1
      if (C(K,I) .gt. C(K,J)) goto 5
      CC = C(K,I)
      C(K,I) = C(K,J)
      C(K,J) = CC
    5 continue
    4 continue

*     order quality shots in effluent
      K = 2
      do 7 I = NS, K95, -1
      do 10 J = 1, I - 1
	if (C(K,I) .gt. C(K,J)) goto 10
      CC = C(K,I)
      C(K,I) = C(K,J)
      C(K,J) = CC
   10 continue
    7 continue
    
      return
      end


      
      
*     compute the mean and standard deviation of logged variables   
      subroutine LogFCf
      include 'bacom1cp.for'

*     river flow
*     first check for non-parametric distributions
      if (nprf .eq. 0 ) then
      if (rfm .lt. Small) then
      grfm = - Big
      grfs =   0.0
      else

*     standard deviation in log domain
      grfs = SQRT( 2.7057 + 2.*ALOG((rfm + fsh)/(rf5 + fsh)) )
     &         - 1.6449

*     mean in log domain
      grfm = ALOG (rfm + fsh) - 0.5 * grfs * grfs
*     compute standard deviation of unlogged flows   
      rfs = rfm * SQRT ( EXP( grfs * grfs ) - 1.0 )
      endif
      
      endif

*     river quality
      if (rcm .lt. Small) then
      grcm = - Big
      grcs = 0.0
      else
*     mean in log domain
      grcm = ALOG( (rcm*rcm) / SQRT(rcm*rcm + rcs*rcs) )
*     standard deviation in log domain
      grcs = SQRT( ALOG(1.0 + (rcs*rcs)/(rcm*rcm)) )
      endif

      if (efm .lt. Small) then ! effluent flow
      gefm = - Big
      gefs =   0.0
      else
*     mean in log domain
      gefm = ALOG( (efm*efm) / SQRT(efm*efm + efs*efs) )
*     standard deviation in log domain
      gefs = SQRT( ALOG(1.0 + (efs*efs)/(efm*efm)) )
      endif

      return
      end




*     set starting value of effluent quality ...  

      subroutine steff
      include 'bacom1cp.for'

      IDIL = 0
      KDIL = 0
      
*     guess the mean river quality downstream of discharge 
      if (type .ne. 1) then
      trcm = 0.6 * targ
      else
      trcm = targ
      endif
      if ( free .eq. 1 ) trcm = 28.0 * 0.021 

*     guess the mean discharge quality
      tecm = amax1 ( 0.01 * targ, 
     &       (trcm * (rfm + efm)-(rcm * rfm))/efm)
*     standard deviation - preserve the coefficient of variation, ecv
      tecs = ecv * tecm

      return
      end



      subroutine fail ! Convergence failure detected
      include 'bacom1cp.for'

      PCFAC = 100.0 * CFAC

*     write(3,3)PCFAC
    3 format(8x,'###  Problem failed to converge to within',F5.2,
     &'% of the target'/
     &8x,'###  Check that the results make sense ....'/68('-'))

      return
      end




*     sensitivity analysis - backwards calculation
      subroutine sens1
      include 'bacom1cp.for'

      effmas = tecm 
      
      if ( tecm .lt. 0.001 ) then
      sens = 999
      write(3,4430)
 4430 format(//
     &'The sensitivity analysis was suppressed ...'/
     &'Effluent quality is zero ...'//)
      return
      endif

      if ( nprf .eq. 0 ) then
      rfm = rfm / 0.9
      call do the calculations (1)
      x1 = 100.0 * tecm / effmas - 100.0
      rfm = xrfm

      rf5 = rf5 / 0.9
      call do the calculations (1)
      x2 = 100.0 * tecm / effmas - 100.0
      rf5 = xrf5
      
      fsh = 0.5 * rf5
      call do the calculations (1)
      x12 = 100.0 * tecm / effmas - 100.0
      fsh = xfsh
      endif

      if ( npef .eq. 0 ) then
      efm = efm * 0.9
      call do the calculations (1)
      x3 = 100.0 * tecm / effmas - 100.0
      efm = xefm
      efs = efs / 0.9
      call do the calculations (1)
      x4 = 100.0 * tecm / effmas - 100.0
      efs = xefs
      endif

      if ( nprc .eq. 0 ) then
      rcm = 0.9 * rcm
      call do the calculations (1)
      x5 = 100.0 * tecm / effmas - 100.0
      rcm = xrcm
      rcs = 0.9 * rcs
      call do the calculations (1)
      x6 = 100.0 * tecm / effmas - 100.0
      rcs = xrcs
      endif

      if ( npec .eq. 0 ) then
      ecv = ecv / 0.9
      call do the calculations (1)
      x7 = 100.0 * tecm / effmas - 100.0
      ecv = xecs/xecm
      endif

      targ = targ / 0.9
      call do the calculations (1)
      x8 = 100.0 * tecm / effmas - 100.0
      targ = xtarg

      if (nprf+npef .eq. 0) then
      CoFf2=CoFf2-sign(0.1, CoFf2)
      call do the calculations (1)
      x9 = 100.0 * tecm / effmas - 100.0
      CoFf2=xCoFf2   
      endif

      if (nprf+nprc .eq. 0) then
      CoFc1=CoFc1-sign(0.1, CoFc1)
      call do the calculations (1)
      x10 = 100.0 * tecm / effmas - 100.0
      CoFc1=xCoFc1
      endif
     
      if (npef+npec .eq. 0) then
      Cofc5=Cofc5-sign(0.1, Cofc5)
      call do the calculations (1)
      x11 = 100.0 * tecm / effmas - 100.0
      Cofc5=xCofc5
      endif
      
*     reset value in case Sensa1 is called immediately after Sens1           
      tecm = effmas
	
      write(3,7742)x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12
 7742 format(
     &'Impact on the discharge standard of a 10% change in ...'/68('-')/
     &'Mean upstream river flow:                  ',f9.2,' %'/
     &'95-percentile river flow:                  ',f9.2,' %'/
     &'Mean discharge flow:                       ',f9.2,' %'/
     &'Standard deviation:                        ',f9.2,' %'/
     &'Mean upstream river quality:               ',f9.2,' %'/
     &'Standard deviation:                        ',f9.2,' %'/
     &'CoV of discharge quality:                  ',f9.2,' %'/
     &'River target:                              ',f9.2,' %'/
     &'Correlation of F on f (change of 0.1):     ',f9.2,' %'/
     &'Correlation of F on C (change of 0.1):     ',f9.2,' %'/
     &'Correlation of f on c (change of 0.1):     ',f9.2,' %'/
     &'River flow shift (half of q95):            ',f9.2,' %'/
     &68('-'))

      return
      end



      
      
*     sensitivity analysis - forwards calculation
      subroutine sens2
      include 'bacom1cp.for'

      call do the calculations (1)

      indax = NS+1
      tcxmas = c(1,NS+1)

      if ( tcxmas .lt. 0.001 ) then
      sens = 999
      write(3,4430)
 4430 format(//
     &'The sensitivity analysis was suppressed ...'/
     &'Effluent quality is zero ...'//)
      return
      endif

      if ( nprf .eq. 0 ) then
      rfm = rfm / 0.9          
      call do the calculations (1)
      x1 = 100.0 * c(1,indax) / tcxmas - 100.0
      rfm = xrfm

      rf5 = rf5 / 0.9
      call do the calculations (1)
      x2 = 100.0 * c(1,indax) / tcxmas - 100.0
      rf5 = xrf5

      fsh = 0.5 * rf5
      call do the calculations (1)
      x12 = 100.0 * c(1,indax) / tcxmas - 100.0
      fsh = xfsh
      endif

      if ( npef .eq. 0 ) then
      efm = 0.9 * efm
      call do the calculations (1)
      x3 = 100.0 * c(1,indax) / tcxmas - 100.0
      efm = xefm

      efs = efs / 0.9
      call do the calculations (1)
      x4 = 100.0 * c(1,indax) / tcxmas - 100.0
      efs = xefs
      endif

      if ( nprc .eq. 0 ) then
      rcm = 0.9 * rcm
      call do the calculations (1)
      x5 = 100.0 * c(1,indax) / tcxmas - 100.0
      rcm = xrcm
      
      rcs = 0.9 * rcs
      x6 = 100.0 * c(1,indax) / tcxmas - 100.0
      rcs = xrcs
      endif

*     mean discharge quality
      if ( npec .eq. 0 ) then
      ecm = 0.9 * ecm 
      call do the calculations (1)
      x7 = 100.0 * c(1,indax) / tcxmas - 100.0
      ecm = xecm
    
*     standard deviation for discharge quality
      ecs = ecs / 0.9
      call do the calculations (1)
   	x8 = 100.0 * c(1,indax) / tcxmas - 100.0
      ecs = xecs
      endif

      if (nprf+npef .eq. 0) then
      CoFf2=CoFf2-sign(0.1, CoFf2)
      call do the calculations (1)
      x9  =  100.0 * c(1,indax) / tcxmas - 100.0
      CoFf2=xCoFf2   
      endif
        
      if (nprf+nprc .eq. 0) then
      CoFc1=CoFc1-sign(0.1, CoFc1)
      call do the calculations (1)
      x10 = 100.0 * c(1,indax) / tcxmas - 100.0
      CoFc1=xCoFc1
      endif
     
      if (npef+npec .eq. 0) then
      Cofc5=Cofc5-sign(0.1, Cofc5)
      call do the calculations (1)
      x11 = 100.0 * c(1,indax) / tcxmas - 100.0
      Cofc5=xCofc5
      endif

*     reset values in case Sensa2 called after Sens2               
      c(1,K90) = effmas
	
      write(3,7742)x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12
 7742 format('Impact on the mean river quality of a 10% ',
     &' change in ...'/68('-')/
     &'Mean upstream river flow (F):                 ',f10.2,' %'/
     &'95-percentile river flow:                     ',f10.2,' %'/
     &'Mean discharge flow (f):                      ',f10.2,' %'/
     &'Standard deviation:                           ',f10.2,' %'/
     &'Mean upstream river quality (C):              ',f10.2,' %'/
     &'Standard deviation:                           ',f10.2,' %'/
     &'Mean discharge quality (c):                   ',f10.2,' %'/
     &'Standard deviation:                           ',f10.2,' %'/
     &'Correlation of F on f (change of 0.1):        ',f10.2,' %'/
     &'Correlation of F on C (change to 0.1):        ',f10.2,' %'/
     &'Correlation of f on c (change to 0.1):        ',f10.2,' %'/
     &'River flow distribution (shift = half of q95):',f10.2,' %'/
     &68('-'))

      return
      end



      subroutine next (ichek)
      include 'bacom1cp.for'

      ichek = 0

      if ( ITER .eq. 1 ) then ! set up the data for the first trial
*     set the values to be used to calculate data for the next trial
      if ( type .ne. 1 ) then
      tem1 = rcm ! effluent quality
      else
      tem1 = 0.5*rcm
      endif
      out1 = 0.0

*     store the results of the first trial
      tem2 = tecm
      if ( type .ne. 1 ) then
      out2 = C(1,KTG)
      else
      out2 = tdm
      endif 
      if (tem1 .eq. tem2) tem1 = 0.956 * tem1
      goto 891
      endif

      if (ITER .eq. 2 ) then ! second iteration
      tem1 = tecm ! store the results of the last trial
      if ( type .ne. 1 ) then
      out1 = C(1,KTG) 
      else
      out1 = tdm
      endif     
      goto 891
      endif

*     Iterations after the first two ...
*     Check which of the two retained trials is best ...

      if (ABS ( out1 - targ ) .lt. ABS ( out2 - targ ) ) then

*     over-write second set
      tem2 = tecm
      if ( type .ne. 1 ) then
      out2 = C(1,KTG)
      else
      out2 = tdm
      endif    
      else

*     over-write the first set
      tem1 = tecm
      if ( type .ne. 1 ) then
      out1 = C(1,KTG) 
      else
      out1 = tdm
      endif    
      endif

*     interpolate between latest trial and retained trial
  891 denom = out2 - out1
      if (denom. gt. Small .or. denom .lt. -Small) then
      tecm = tem1 + (targ-out1) * (tem2-tem1) / denom
      tecs = ecv * tecm
      else
      ichek = 1
      return
      endif

*     check for zero discharge quality
      if ( tecm .gt. 0.0 ) return
*     discharge quality is zero. Replace with fraction of previous trials.
      tecm = 0.1 * AMIN1( tem1, tem2 ) ! #############################
      tecs = ecv * tecm
*     check again for zero ...
      if (tecm .gt. 0.0) return ! #####################################
*     is this the second time that a zero has been obtained ?
*     if so, a fault has been detected ...

      if (KDIL .gt. 0) then
      msg1 = 1

*     set up a calculation to check whether target can be achieved at all
      else
      IDIL = 1
      KDIL = 1
      tecm = 0.0
      tecs = 0.0
      return
      endif    

      return
      end


      subroutine set the starting values of the data 
      include 'bacom1cp.for'

*     precision factor - iterative scheme stops when successive
*     trials are within 100*CFAC %
      CFAC = 0.00001

      one = 0.9999 ! test parameter for correlation coefficients

      Small = 1.0e-20
      Big = 1.0e5    
      
      do i = 1, 8
      BM(i) = 1.0
      BS(i) = 1.0
      enddo

      RR = -5119
      RR = -1858
      
      JR1 = -5119 ! random number starter - river flow
      JR3 = -1849 + 9 ! discharge flow
      JR2 = -7849 + 6 ! river quality for dissolved metal 
      JR4 = -5149 + 1 ! discharge quality for dissolved metal
      JR5 = -7765 ! discharge quality for dissolved organic carbon 
      JR6 = -8849 ! dissolved organic carbon in the upstream river
      JR7 = -1483 ! pH in downstream river
      JR8 = -5849 ! calcium in downstream river
      JR9 = -3849 ! errors in the bioavailable equation ========================

      KR1 = -5119 ! random number statert - river flow
      KR3 = -1849 + 9 ! discharge flow
      KR2 = -7849 + 6 ! river quality for dissolved metal 
      KR4 = -5149 + 1 ! discharge quality for dissolved metal
      KR5 = -7765 ! discharge quality for dissolved organic carbon 
      KR6 = -8849 ! dissolved organic carbon in the upstream river
      KR7 = -1483 ! pH in downstream river
      KR8 = -5849 ! calcium in downstream river
      KR9 = -3849 ! errors in the bioavailable equation ========================

      KTG = INT(0.5 + 0.01 * XPER * NS)
      K05 = INT(0.5 + 0.05 * NS)
      K80 = INT(0.5 + 0.80 * NS)
      K90 = INT(0.5 + 0.90 * NS)
      K95 = INT(0.5 + 0.95 * NS)
      K98 = INT(0.5 + 0.98 * NS)
      K99 = INT(0.5 + 0.99 * NS)
      K995 = INT(0.5 + 0.995 * NS)
      K999 = INT(0.5 + 0.999 * NS)

      RST = -6411
      TIM =  TDIST1 (9999.0, (1.0-0.01*spc) )
      kevin = 0 ! iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

      rs = -2379
      npmax = 1000
      
      if (rgm .le. Small) then
      vrgs =  0.0
      rgs=0.0
      grgm = - Big
      grgs = 0.0     
      else
      grgm = ALOG( (rgm*rgm)/SQRT(rgm*rgm + rgs*rgs) )
      grgs = SQRT( ALOG( 1.0 + (rgs*rgs)/(rgm*rgm)) )
      endif

      nprc = 0
      nprf = 0
      npef = 0
      npec = 0
      npdi = 0 ! iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

      if (type .eq. 3) then ! upper-tier uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
      type = 0
      iutflag = 1
      else
      iutflag = 0
      endif ! uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu

      msg1 = 0
      msg2 = 0
      msg3 = 0
      msg4 = 0
      
      do j = 1, 5000
      do i = 1, 4
      C(i,j) = 0.0
      enddo
      do i = 1, 8
      TMC(i,j) = 0.0
      enddo
      enddo

      bioshotCUm  =  20.0 ! % bias in the bio-equation (20%) bbbbbbbbbbbbbbbbbbb
      bioshotNIm  =  20.0 ! bias in the bio-equation (20%) bbbbbbbbbbbbbbbbbbbbb
      bioshotZNm  =  25.0 ! bias in the bio-equation (25%) bbbbbbbbbbbbbbbbbbbbb
      bioshotMNm  =  15.0 ! bias in the bio-equation (85%) bbbbbbbbbbbbbbbbbbbbb
      bioshotPBm  =   0.0 ! bias in the bio-equation (85%) bbbbbbbbbbbbbbbbbbbbb
      
      bioshotm = 0.0
            
      bioshotCUm  =  20.0 ! % bias in the bio-equation (20%) bbbbbbbbbbbbbbbbbbb
      bioshotCUs1 =  17.0 ! % standard deviation from the bio-equation (85%) bbb
      bioshotCUs2 =  17.0 ! % standard deviation from the bio-equation (85%) bbb
      bioshotNIm  =  20.0 ! bias in the bio-equation (20%) bbbbbbbbbbbbbbbbbbbbb
      bioshotNIs1 =   7.4 ! % standard deviation from the bio-equation (37%) bbb
      bioshotNIs2 =   5.6 ! % standard deviation from the bio-equation (28%) bbb
      bioshotZNm  =  25.0 ! bias in the bio-equation (25%) bbbbbbbbbbbbbbbbbbbbb
      bioshotZNs1 =  6.25 ! % standard deviation from the bio-equation (25%) bbb
      bioshotZNs2 = 21.25 ! % standard deviation from the bio-equation (85%) bbb
      bioshotMNm  =  15.0 ! bias in the bio-equation (85%) bbbbbbbbbbbbbbbbbbbbb
      bioshotMNs1 =   8.0 ! % standard deviation from the bio-equation bbbbbbbbb
      bioshotMNs2 =   8.0 ! % standard deviation from the bio-equation bbbbbbbbb
      bioshotPBm  =   0.0 ! bias in the bio-equation (85%) bbbbbbbbbbbbbbbbbbbbb
      bioshotPBs1 =   0.0 ! % standard deviation from the bio-equation bbbbbbbbb
      bioshotPBs2 =   0.0 ! % standard deviation from the bio-equation bbbbbbbbb
      
      bioshots1 = 0.0
      bioshots2 = 0.0

      return 
      end

      


      subroutine Correl(nvalid)

*     Takes 6 correlation coefficients, checks their validity, scales
*     them if necessary to make them valid, and returns nvalid=1 if they
*     have been scaled.

*     'test' is a scaling factor and testhigh and testlow are the
*     max and min values of test when iterating using the bisection
*     method.

      include 'bacom1cp.for'

      sCoFC1=CoFC1
      sCoFf2=CoFf2
      sCoFc3=CoFc3
      sCoCf4=CoCf4
      sCofc5=Cofc5
      sCoCc6=CoCc6

      testhigh=1
      testlow=0
      test=testhigh

      call Check4(jvalid,nvalid)
      if (jvalid .eq. 0) then
      return
      endif

      test=0
      do while (testhigh-testlow .gt. 0.005)
      CoFC1=sCoFC1*test
      CoFf2=sCoFf2*test
      CoFc3=sCoFc3*test
      CoCf4=sCoCf4*test
      Cofc5=sCofc5*test
      CoCc6=sCoCc6*test
      call Check4(jvalid,nvalid)
      If (jvalid .eq. 0) then
      testlow=test
      else
      testhigh=test
      endif
      test=0.5*(testhigh+testlow)
      enddo

      test=testlow
      CoFC1=sCoFC1*test
      CoFf2=sCoFf2*test
      CoFc3=sCoFc3*test
      CoCf4=sCoCf4*test
      Cofc5=sCofc5*test
      CoCc6=sCoCc6*test
    
      call Check4(jvalid,nvalid)
      nvalid=1

      return
      end


      subroutine Check4(jvalid,nvalid)

*     Checks validity of 6 correlation coefficients for 4 variables.
*     Each iteration involves checking CCs in turn for 4 sets of 3 variables,
*     scaling any set of 3 CCs as necessary, and then checking for
*     valid values of c3, d3 and d4.

      include 'bacom1cp.for'

      nvalid=0

      call Check3(CoFC1,CoFf2,CoCf4,f)
      if ( f .lt. 1 ) then
      nvalid=1
      endif

      call Check3(CoFC1,CoFc3,CoCc6,f)
      if ( f .lt. 1 ) then
      nvalid=1
      endif

      call Check3(CoFf2,CoFc3,Cofc5,f)
      if ( f .lt. 1 ) then
      nvalid=1
      endif

      call Check3(CoCf4,Cofc5,CoCc6,f)
      if ( f .lt. 1 ) then
      nvalid=1
      endif

      b1=CoFC1
      b2=sqrt(1.-b1*b1)
      c1=CoFf2

      if (b2 .gt. 0.) then
      c2 = (CoCf4 - b1 * c1) / b2
      else
      c2 = 0.
      endif

      ck = 1. - c1 * c1 - c2 * c2
      if (ck .gt. 0.) then
      c3 = sqrt (ck)
      else
      c3 = 0.
      endif

      d1 = CoFc3

      if (b2 .gt. 0.) then
      d2=(CoCc6 - b1 * d1 ) /b2
      else
      d2 = 0.
      endif

      if (c3 .gt. 0.) then
      d3 = (Cofc5 - c1 * d1 - c2 * d2) / c3
      else
      d3 = 0.
      endif

      ck = 1. - d1 * d1 - d2 * d2 - d3 * d3
      if (ck .gt. 0.) then
      d4 = sqrt (ck)
      else
      d4 = 0.
      endif

      jvalid=0
      if ( abs(d3) .gt. 1.0 .or. ck .lt. -0.00001 ) jvalid=1

*     correlation between calcium and pH
      f1 = cacp
      sq = 1.0 - f1*f1
      if (sq .gt. 0.00001) then
      f2 = SQRT (sq)
      else
      f2 = 0.0
      endif

*     correlation between pH and flow
      p1 = capf
      sq = 1.0 - p1*p1
      if (sq .gt. 0.00001) then
      p2 = SQRT (sq)
      else
      p2 = 0.0
      endif

*     correlation between calcium and flow
      z1 = cacf
      sq = 1.0 - z1*z1
      if (sq .gt. 0.00001) then
      z2 = SQRT (sq)
      else
      z2 = 0.0
      endif

      return
      end


      subroutine Check3 ( corr1, corr2, corr3, f )

*     Checks the validity of the three correlation coefficients relating
*     to a set of three variables. f is a scaling factor and fhigh and flow
*     are the max and min values of f when iterating using the bisection
*     method.

      scorr1=corr1
      scorr2=corr2
      scorr3=corr3

      fhigh=1
      flow=0
      f=fhigh

      temp=1.0-corr1*corr1-corr2*corr2-corr3*corr3+2*corr1*corr2*corr3
      if (temp .ge.0.0) then
      return
      endif

      f=0
      do while (fhigh-flow .gt. 0.005)
      corr1=scorr1*f
      corr2=scorr2*f
      corr3=scorr3*f
      temp=1.0-corr1*corr1-corr2*corr2-corr3*corr3+2*corr1*corr2*corr3
      If (temp.ge.0.0) then
      flow=f
      else
      fhigh=f
      endif
      f=0.5*(fhigh+flow)
      enddo

      f=flow
      corr1=scorr1*f
      corr2=scorr2*f
      corr3=scorr3*f

      return
      end




      function GAS1(IDUM) ! ---------------------------- random number gerarator
      save ISET,GSET
      data ISET/0/
      if (IDUM .lt. 0) ISET = 0
      if (ISET .eq. 0) then
    1 continue  
      V1 = 2.*RAN1(IDUM)-1.0
	V2 = 2.*RAN1(IDUM)-1.0
	R = V1**2 + V2**2
	if (R .ge. 1.) goto 1
	FAC = SQRT(-2.*LOG(R)/R)
	GSET = V1*FAC
	GAS1 = V2*FAC
	ISET = 1
      else
	GAS1 = GSET
	ISET = 0
      endif
      return
      end ! -------------------------------------------- random number gerarator
      function GAS5(IDUM) ! ---------------------------- random number gerarator
      save iset,gset
      data ISET/0/
      if (IDUM .lt. 0) ISET = 0
      if (ISET .eq. 0) then
    1 V1 = 2.*RAN5(IDUM)-1.0
	V2 = 2.*RAN5(IDUM)-1.0
	R = V1**2 + V2**2
	if (R.GE.1.) goto 1
	FAC = SQRT(-2.*LOG(R)/R)
	GSET = V1*FAC
	GAS5 = V2*FAC
	ISET = 1
      else
	GAS5 = GSET
	ISET = 0
      endif
      return
      end ! -------------------------------------------- random number gerarator
      function GAS6(IDUM) ! ---------------------------- random number gerarator
      save iset,gset
      data ISET/0/
      if (IDUM .lt. 0) ISET = 0
      if (ISET .eq. 0) then
    1 V1 = 2.*RAN6(IDUM)-1.0
	  V2 = 2.*RAN6(IDUM)-1.0
	R = V1**2 + V2**2
	if (R.GE.1.) goto 1
	FAC = SQRT(-2.*LOG(R)/R)
	GSET = V1*FAC
	GAS6 = V2*FAC
	ISET = 1
      else
	GAS6 = GSET
	ISET = 0
      endif
      return
      end ! -------------------------------------------- random number gerarator
      function GAS7(IDUM) ! ---------------------------- random number gerarator
      save iset,gset
      data ISET/0/
      if (IDUM .lt. 0) ISET = 0
      if (ISET .eq. 0) then
    1 V1 = 2.*RAN7(IDUM)-1.0
	V2 = 2.*RAN7(IDUM)-1.0
	R = V1**2 + V2**2
	if (R.GE.1.) goto 1
	FAC = SQRT(-2.*LOG(R)/R)
	GSET = V1*FAC
	GAS7 = V2*FAC
	ISET = 1
      else
	GAS7 = GSET
	ISET = 0
      endif
      return
      end ! -------------------------------------------- random number gerarator

      function RAN1(IDUM) ! ============================ random number gerarator
      dimension R(97)
      parameter (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      parameter (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      parameter (M3=243000,IA3=4561,IC3=51349)
      save IX1,IX2,IX3
      if (IDUM .lt. 0 ) then
	IFF = 1
	IX1 = MOD (IC1-IDUM,M1)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX2 = MOD (IX1,M2)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX3 = MOD (IX1,M3)
	DO J = 1,97
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      enddo
	IDUM = 1
      endif
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      IX3 = MOD (IA3*IX3 + IC3,M3)
      JJ = 1 + (97*IX3)/M3
      J = iabs(JJ)
      if (J .gt. 97 .or. J .lt. 1)PAUSE
      RAN1 = R(J)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      return
      end ! ============================================ random number gerarator
      function RAN5(IDUM) ! ============================ random number gerarator
      dimension R(97)
      parameter (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      parameter (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      parameter (M3=243000,IA3=4561,IC3=51349)
      save IX1,IX2,IX3
      if (IDUM .lt. 0 ) then
	IFF = 1
	IX1 = MOD (IC1-IDUM,M1)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX2 = MOD (IX1,M2)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX3 = MOD (IX1,M3)
	DO J = 1,97
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      enddo
	IDUM = 1
      endif
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      IX3 = MOD (IA3*IX3 + IC3,M3)
      J = 1 + (97*IX3)/M3
      if (J .gt. 97.OR.J .lt. 1)PAUSE
      RAN5 = R(J)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      return
      end ! ============================================ random number gerarator
      function RAN6(IDUM) ! ============================ random number gerarator
      dimension R(97)
      parameter (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      parameter (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      parameter (M3=243000,IA3=4561,IC3=51349)
      save IX1,IX2,IX3
      if (IDUM .lt. 0 ) then
	IFF = 1
	IX1 = MOD (IC1-IDUM,M1)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX2 = MOD (IX1,M2)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX3 = MOD (IX1,M3)
	DO J = 1,97
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      enddo
	IDUM = 1
      endif
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      IX3 = MOD (IA3*IX3 + IC3,M3)
      J = 1 + (97*IX3)/M3
      if (J .gt. 97.OR.J .lt. 1)PAUSE
      RAN6 = R(J)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      return
      end ! ============================================ random number gerarator
      function RAN7(IDUM) ! ============================ random number gerarator
      dimension R(97)
      parameter (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      parameter (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      parameter (M3=243000,IA3=4561,IC3=51349)
      save IX1,IX2,IX3
      if (IDUM .lt. 0 ) then
	IFF = 1
	IX1 = MOD (IC1-IDUM,M1)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX2 = MOD (IX1,M2)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX3 = MOD (IX1,M3)
	DO J = 1,97
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      enddo
	IDUM = 1
      endif
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      IX3 = MOD (IA3*IX3 + IC3,M3)
      J = 1 + (97*IX3)/M3
      if (J .gt. 97.OR.J .lt. 1)PAUSE
      RAN7 = R(J)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      return
      end ! ============================================ random number gerarator

      
      
      function GAS1s(IDUM) ! ---------------------------- random number gerarator
      save ISET,GSET
      data ISET/0/
      if (IDUM .lt. 0) ISET = 0
      if (ISET .eq. 0) then
    1 continue  
      V1 = 2.*RAN1s(IDUM)-1.0
	V2 = 2.*RAN1s(IDUM)-1.0
	R = V1**2 + V2**2
	if (R .ge. 1.) goto 1
	FAC = SQRT(-2.*LOG(R)/R)
	GSET = V1*FAC
	GAS1s = V2*FAC
	ISET = 1
      else
	GAS1s = GSET
	ISET = 0
      endif
      return
      end ! -------------------------------------------- random number gerarator
      function GAS5s(IDUM) ! ---------------------------- random number gerarator
      save iset,gset
      data ISET/0/
      if (IDUM .lt. 0) ISET = 0
      if (ISET .eq. 0) then
    1 V1 = 2.*RAN5s(IDUM)-1.0
	V2 = 2.*RAN5s(IDUM)-1.0
	R = V1**2 + V2**2
	if (R.GE.1.) goto 1
	FAC = SQRT(-2.*LOG(R)/R)
	GSET = V1*FAC
	GAS5s = V2*FAC
	ISET = 1
      else
	GAS5s = GSET
	ISET = 0
      endif
      return
      end ! -------------------------------------------- random number gerarator
      function GAS6s(IDUM) ! ---------------------------- random number gerarator
      save iset,gset
      data ISET/0/
      if (IDUM .lt. 0) ISET = 0
      if (ISET .eq. 0) then
    1 V1 = 2.*RAN6s(IDUM)-1.0
      V2 = 2.*RAN6s(IDUM)-1.0
	R = V1**2 + V2**2
	if (R.GE.1.) goto 1
	FAC = SQRT(-2.*LOG(R)/R)
	GSET = V1*FAC
	GAS6s = V2*FAC
	ISET = 1
      else
	GAS6s = GSET
	ISET = 0
      endif
      return
      end ! -------------------------------------------- random number gerarator
      function GAS7s(IDUM) ! ---------------------------- random number gerarator
      save iset,gset
      data ISET/0/
      if (IDUM .lt. 0) ISET = 0
      if (ISET .eq. 0) then
    1 V1 = 2.*RAN7s(IDUM)-1.0
	V2 = 2.*RAN7s(IDUM)-1.0
	R = V1**2 + V2**2
	if (R.GE.1.) goto 1
	FAC = SQRT(-2.*LOG(R)/R)
	GSET = V1*FAC
	GAS7s = V2*FAC
	ISET = 1
      else
	GAS7s = GSET
	ISET = 0
      endif
      return
      end ! -------------------------------------------- random number gerarator

      function RAN1s(IDUM) ! ============================ random number gerarator
      dimension R(97)
      parameter (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      parameter (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      parameter (M3=243000,IA3=4561,IC3=51349)
      save IX1,IX2,IX3
      if (IDUM .lt. 0 ) then
	IFF = 1
	IX1 = MOD (IC1-IDUM,M1)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX2 = MOD (IX1,M2)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX3 = MOD (IX1,M3)
	DO J = 1,97
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      enddo
	IDUM = 1
      endif
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      IX3 = MOD (IA3*IX3 + IC3,M3)
      JJ = 1 + (97*IX3)/M3
      J = iabs(JJ)
      if (J .gt. 97 .or. J .lt. 1)PAUSE
      RAN1s = R(J)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      return
      end ! ============================================ random number gerarator
      function RAN5s(IDUM) ! ============================ random number gerarator
      dimension R(97)
      parameter (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      parameter (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      parameter (M3=243000,IA3=4561,IC3=51349)
      save IX1,IX2,IX3
      if (IDUM .lt. 0 ) then
	IFF = 1
	IX1 = MOD (IC1-IDUM,M1)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX2 = MOD (IX1,M2)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX3 = MOD (IX1,M3)
	DO J = 1,97
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      enddo
	IDUM = 1
      endif
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      IX3 = MOD (IA3*IX3 + IC3,M3)
      J = 1 + (97*IX3)/M3
      if (J .gt. 97.OR.J .lt. 1)PAUSE
      RAN5s = R(J)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      return
      end ! ============================================ random number gerarator
      function RAN6s(IDUM) ! ============================ random number gerarator
      dimension R(97)
      parameter (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      parameter (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      parameter (M3=243000,IA3=4561,IC3=51349)
      save IX1,IX2,IX3
      if (IDUM .lt. 0 ) then
	IFF = 1
	IX1 = MOD (IC1-IDUM,M1)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX2 = MOD (IX1,M2)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX3 = MOD (IX1,M3)
	DO J = 1,97
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      enddo
	IDUM = 1
      endif
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      IX3 = MOD (IA3*IX3 + IC3,M3)
      J = 1 + (97*IX3)/M3
      if (J .gt. 97.OR.J .lt. 1)PAUSE
      RAN6s = R(J)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      return
      end ! ============================================ random number gerarator
      function RAN7s(IDUM) ! ============================ random number gerarator
      dimension R(97)
      parameter (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      parameter (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      parameter (M3=243000,IA3=4561,IC3=51349)
      save IX1,IX2,IX3
      if (IDUM .lt. 0 ) then
	IFF = 1
	IX1 = MOD (IC1-IDUM,M1)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX2 = MOD (IX1,M2)
	IX1 = MOD (IA1*IX1 + IC1,M1)
	IX3 = MOD (IX1,M3)
	DO J = 1,97
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      enddo
	IDUM = 1
      endif
      IX1 = MOD (IA1*IX1 + IC1,M1)
      IX2 = MOD (IA2*IX2 + IC2,M2)
      IX3 = MOD (IA3*IX3 + IC3,M3)
      J = 1 + (97*IX3)/M3
      if (J .gt. 97.OR.J .lt. 1)PAUSE
      RAN7s = R(J)
      R(J) = (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
      return
      end ! ============================================ random number gerarator


*     read in the data for a non-parametric distribution of river flow

      subroutine rfnprd
      logical exists
      character*80 string

      include 'bacom1cp.for'

*     attach data file
 
      Inquire(FILE=rffile,EXIST=exists)
      if ( .NOT. exists) then
      msg4 = msg4 + 1   
      return
      endif

      open(10,FILE=rffile,STATUS='OLD') ! open the file ------------------------
      read (10, *, ERR=20) nprf ! read the file --------------------------------
      if(nprf.eq.-1)goto 20
      iflag=0
      goto 50

*     test for Aardvark data
   20 read(10,'(A1)',err=7500)ans
      i=1
   30 read(10,31,end=40)string
   31 format(A80)
      read(string,*,err=30)idum1,idum2,idum3,cptemp
      rfnpvl(i)=cptemp
      i=i+1
      goto 30
   40 close(10)
      nprf=i-1
      iflag=1

*     check the number of data points - are there too many ?      
   50 if (nprf .gt. npmax) then
      msg4 = msg4 + 3
      return
      endif

*     check the number of data points - are there too few ?      

      if (nprf .lt. 5) then
      msg4 = msg4 + 4
      return
      endif

      if (iflag .eq. 0) then
      backspace (10)
      read (10, *, ERR=7500) nprf, (rfnpvl(i),i=1 , nprf)
      rewind (10)
      close(10)
      endif
      
*     arrange the data in sequence

      do 11 i = 1,   nprf-1
      do 12 j = i+1, nprf
      if (rfnpvl(i) .ge. rfnpvl(j)) goto 12
      x = rfnpvl (j)
      rfnpvl (j) = rfnpvl (i)
      rfnpvl (i) = x
   12 continue 
   11 continue

*     compute cumulative frequencies and store them in cfd
      CUMUL = 1.0 / float (nprf + 1)
      do 10 i = 1, nprf
   10 rfnpfd(i) = float(i) * CUMUL

*     compute mean and 5-percentile
      rfm = 0.0
      do 13 i = 1, nprf
      rfm = rfm + rfnpvl (i)
   13 continue
      rfm = rfm / float (nprf)  
      R1 = 0.95
      call rfnpar (rf5)

      return

 7500 msg4=msg4+2
      return
      end


*     read in the data for a non-parametric distribution of river qual.
      subroutine rcnprd

      logical exists
      character*80 string

      include 'bacom1cp.for'

*     attach data file
      Inquire(FILE=rcfile,EXIST=exists)
      if ( .NOT. exists) then
      msg4 = msg4 + 10
      return
      endif

*     open the file
      open(10,FILE=rcfile,STATUS='OLD')

*     read the file
      read (10, *, ERR=20) nprc
      if(nprc.eq.-1)goto 20
      iflag=0
      goto 50

*     Test for Aardvark data
   20 read(10,'(A1)',err=7500)ans
*     read(10,'(A1)',err=7500)ans
      i=1
   30 read(10,31,end=40)string
   31 format(A80)
      read(string,*,err=30)idum1,idum2,idum3,cptemp
      rcnpvl(i)=cptemp
      i=i+1
      goto 30
   40 close(10)
      nprc=i-1
      iflag=1

*     check the number of data points - are there too many ?      
   50 if (nprc .gt. npmax) then
      msg4 = msg4 + 30
      return
      endif

*     check the number of data points - are there too few ?      
      if (nprc .lt. 5) then
      msg4 = msg4 + 40
      return
      endif

      if (iflag .eq. 0) then
      backspace (10)
      read (10, *, ERR=7500) nprc, (rcnpvl(i),i=1 , nprc)
      rewind (10)
      close(10)
      endif

*     arrange the data in sequence
      do 11 i = 1,   nprc-1
      do 12 j = i+1, nprc

      if (rcnpvl(i) .ge. rcnpvl(j)) goto 12
      x = rcnpvl (j)
      rcnpvl (j) = rcnpvl (i)
      rcnpvl (i) = x
   12 continue 
   11 continue

*     compute cumulative frequencies and store them in cfd

      CUMUL = 1.0 / float (nprc + 1)

      do 10 i = 1, nprc
   10 rcnpfd(i) = float(i) * CUMUL

*     compute mean and standard deviation

      rcm = 0.0
      do 13 i = 1, nprc
      rcm = rcm + rcnpvl (i)
   13 continue
      rcm = rcm / float (nprc)  

      rc2=0.0
      do 14 i = 1, nprc
     
*     sum of the x^2
      rc2 = rc2 + rcnpvl(i)**2
   14 continue
      rcs=SQRT((1/(float(nprc)-1))*(rc2-float(nprc)*rcm**2))  

      return

 7500 msg4 = msg4 + 20
      return

      end


*     read in the data for a non-parametric distribution of effluent f.
      subroutine efnprd

      logical exists
      character*80 string

      include 'bacom1cp.for'

*     attach data file
      Inquire(FILE=effile,EXIST=exists)
      if ( .NOT. exists) then
      msg4 = msg4 + 100
      return
      endif

*     open the file
      open(10,FILE=effile,STATUS='OLD')

*     read the file
      read (10, *, ERR=20) npef
      if(npef.eq.-1)goto 20
      iflag=0
      goto 50

*     test for Aardvark data
   20 read(10,'(A1)',err=7500)ans
      i=1
   30 read(10,31,end=40)string
   31 format(A80)
      read(string,*,err=30)idum1,idum2,idum3,cptemp
      efnpvl(i)=cptemp
      i=i+1
      goto 30
   40 close(10)
      npef=i-1
      iflag=1

*     check the number of data points - are there too many ?      
   50 if (npef .gt. npmax) then
      msg4 = msg4 + 300
      return
      endif

*     check the number of data points - are there too few ?      

      if (npef .lt. 5) then
      msg4 = msg4 + 400
      return
      endif

      if (iflag .eq. 0) then
      backspace (10)
      read (10, *, ERR=7500) npef, (efnpvl(i),i=1 , npef)
      rewind (10)
      close(10)
      endif

*     arrange the data in sequence

      do 11 i = 1,   npef-1
      do 12 j = i+1, npef

      if (efnpvl(i) .ge. efnpvl(j)) goto 12

      x = efnpvl (j)
      efnpvl (j) = efnpvl (i)
      efnpvl (i) = x

   12 continue 
   11 continue

*     compute cumulative frequencies and store them in cfd

      CUMUL = 1.0 / float (npef + 1)

      do 10 i = 1, npef

   10 efnpfd(i) = float(i) * CUMUL

*     compute mean and standard deviation

      efm = 0.0
      do 13 i = 1, npef
      efm = efm + efnpvl (i)
   13 continue
      efm = efm / float (npef)  

      ef2=0.0
      do 14 i = 1, npef
     
*     sum of the x^2
     
      ef2 = ef2 + efnpvl(i)**2
   14 continue
      efs=SQRT((1/(float(npef)-1))*(ef2-float(npef)*efm**2))  

      return

 7500 msg4=msg4+200
      return
      end


*     read in the data for a non-parametric distribution of effluent quality
      subroutine ecnprd

      Logical exists
      character*80 string
	              
      include 'bacom1cp.for'

*     attach data file
      Inquire(FILE=ecfile,EXIST=exists)
      if ( .NOT. exists) then
      msg4 =msg4 + 1000
      return
      endif

*     open the file
      open(10,FILE=ecfile,STATUS='OLD')

*     read the file
      read (10, *, ERR=20) npec
      if(npec.eq.-1)goto 20
      iflag=0
      goto 50

*     test for Aardvark data
   20 read(10,'(A1)',err=7500)ans
      i=1
   30 read(10,31,end=40)string
   31 format(A80)
      read(string,*,err=30)idum1,idum2,idum3,cptemp
      ecnpvl(i)=cptemp
      i=i+1
      goto 30
   40 close(10)
      npec=i-1
      iflag=1

*     check the number of data points - are there too many ?
     
   50 if (npec .gt. npmax) then
      msg4 = msg4 + 3000
      return
      endif

*     check the number of data points - are there too few ? 
     
      if (npec .lt. 5) then
      msg4 = msg4 + 4000
      return
      endif

      if (iflag .eq. 0) then
      backspace (10)
      read (10, *, ERR=7500) npec, (ecnpvl(i),i=1 , npec)
      rewind (10)
      close(10)
      endif

*     arrange the data in sequence
      do 11 i = 1,   npec-1
      do 12 j = i+1, npec
      if (ecnpvl(i) .ge. ecnpvl(j)) goto 12
      x = ecnpvl (j)
      ecnpvl (j) = ecnpvl (i)
      ecnpvl (i) = x
   12 continue 
   11 continue

*     compute cumulative frequencies and store them in cfd
      CUMUL = 1.0 / float (npec + 1)

      do 10 i = 1, npec
   10 ecnpfd(i) = float(i) * CUMUL

*     compute mean and standard deviation
      ecm = 0.0
      do 13 i = 1, npec
      ecm = ecm + ecnpvl (i)
   13 continue
      ecm = ecm / float (npec)  

      ec2=0.0
      do 14 i = 1, npec

*     sum of the x^2 
      ec2 = ec2 + ecnpvl(i)**2
   14 continue
      ecs=SQRT((1/(float(npec)-1))*(ec2-float(npec)*ecm**2))  
      ecv = ecs / ecm

      return

 7500 msg4=msg4+2000
      return

      end


*     provide a random sample from a non-parametric distribution 
      subroutine rfnpar (RFN) 

*     rfnpvl    is an array of class marks,                           
*     rfnpfd    is an array of cumulative frequencies,      
*     nprf      is the number of classes,
*     R1        is the uniform variate,                     
*     RFN       is the returned variate.                      

      include 'bacom1cp.for'

*     prepare to deal with values beyond the first and the last values
*     entered as data ....
*     we shall interpolate between the last point and a point half the last
*     interval beyond the last point

*     calculate the first interval
      CLINT1 = rfnpvl (2) - rfnpvl (1)

*     TAU is the minimum - the first value less half the interval
      TAU = rfnpvl (1) - CLINT1 / 2.

*     the last interval 
      CLINT2 = rfnpvl (nprf) - rfnpvl (nprf-1)

*     THETA is the maximum - last value plus half the interval
      THETA = rfnpvl (nprf) + CLINT2 / 2.

*     calculate the gradients of these tails of the distribution
      GRAD1 = ( rfnpvl(1)      - TAU  ) /       rfnpfd (1)
      GRAD2 = ( THETA - rfnpvl (nprf) ) / (1. - rfnpfd (nprf) )

*     locate this point on the cumulative frequency distribution
      if ( R1 .le. rfnpfd(1) ) then
      RFN = GRAD1 * R1 + TAU

      else if ( R1 .ge. rfnpfd(nprf) ) then
      RFN = GRAD2 * R1 + THETA - GRAD2
      else
      CL = float(nprf)
      L = 1
      U = nprf
   10 cxx = (U + float(L))/2.
      I = INT (cxx + 0.5)

      if (R1 .le. rfnpfd(I) .and. R1 .ge. rfnpfd(I - 1)) then

      Y1 = rfnpvl(I - 1)
      Y2 = rfnpvl(I)
      XX1 = rfnpfd(I - 1)
      XX2 = rfnpfd(I)
 
      RFN = (R1 - XX1) * (Y2 - Y1) / (XX2 - XX1) + Y1

      else

      if (R1 .GT. rfnpfd(I)) L = I
      if (R1 .LT. rfnpfd(I)) U = I
      goto 10
   
      endif
      endif

      RFN = amax1( 0.0, RFN )

      return
      end




*     provide a random sample from a non-parametric distribution for river quality
      subroutine rcnpar (RCN) 

*     rcnpvl    is an array of class marks,                           
*     rcnpfd    is an array of cumulative frequencies,      
*     nprc      is the number of classes,
*     R2        is the uniform variate,                     
*     RCN       is the returned variate.                      

      include 'bacom1cp.for'

*     prepare to deal with values beyond the first and the last values
*     entered as data ....
*     we shall interpolate between the last point and a point half the last
*     interval beyond the last point

*     calculate the first interval
      CLINT1 = rcnpvl (2) - rcnpvl (1)

*     TAU is the minimum - the first value less half the interval
      TAU = rcnpvl (1) - CLINT1 / 2.

*     the last interval 
      CLINT2 = rcnpvl (nprc) - rcnpvl (nprc-1)

*     THETA is the maximum - last value plus half the interval
      THETA = rcnpvl (nprc) + CLINT2 / 2.

*     calculate the gradients of these tails of the distribution
      GRAD1 = ( rcnpvl(1)      - TAU  ) /       rcnpfd (1)
      GRAD2 = ( THETA - rcnpvl (nprc) ) / (1. - rcnpfd (nprc) )

*     locate this point on the cumulative frequency distribution.
      if ( R2 .le. rcnpfd(1) ) then

      RCN = GRAD1 * R2 + TAU

      else if ( R2 .ge. rcnpfd(nprc) ) then
      RCN = GRAD2 * R2 + THETA - GRAD2
      else
      CL = float(nprc)
      L = 1
      U = nprc
   10 cxx = (U + float(L))/2.
      I = INT (cxx + 0.5)

      if (R2 .le. rcnpfd(I) .and. R2 .ge. rcnpfd(I - 1)) then

      Y1 = rcnpvl(I - 1)
      Y2 = rcnpvl(I)
      XX1 = rcnpfd(I - 1)
      XX2 = rcnpfd(I)
 
      XXN = (R2 - XX1) * (Y2 - Y1) / (XX2 - XX1) + Y1
      RCN = XXN
      else ! ###### only in AM
      if (R2 .GT. rcnpfd(I)) L = I
      if (R2 .LT. rcnpfd(I)) U = I
      goto 10
   
      endif
      endif

      RCN = amax1( 0.0, RCN )

      return
      end


*     provide a random sample from a non-parametric distribution for effluent flow                   
      subroutine efnpar (EFN) 

*     efnpvl    is an array of class marks,                           
*     efnpfd    is an array of cumulative frequencies,      
*     npef      is the number of classes,
*     R3        is the uniform variate,                     
*     EFN       is the returned variate.                      

      include 'bacom1cp.for'

*     prepare to deal with values beyond the first and the last values
*     entered as data ....
*     we shall interpolate between the last point and a point half the last
*     interval beyond the last point

      R3 = 1.0 - R3

*     calculate the first interval
      CLINT1 = efnpvl (2) - efnpvl (1)
*     TAU is the minimum - the first value less half the interval
      TAU = efnpvl (1) - CLINT1 / 2.

*     the last interval 
      CLINT2 = efnpvl (npef) - efnpvl (npef-1)
*     THETA is the maximum - last value plus half the interval
      THETA = efnpvl (npef) + CLINT2 / 2.

*     calculate the gradients of these tails of the distribution
      GRAD1 = ( efnpvl(1)      - TAU  ) /       efnpfd (1)
      GRAD2 = ( THETA - efnpvl (npef) ) / (1. - efnpfd (npef) )

*     locate this point on the cumulative frequency distribution.

      if ( R3 .le. efnpfd(1) ) then

      EFN = GRAD1 * R3 + TAU

      else if ( R3 .ge. efnpfd(npef) ) then

      XXN = GRAD2 * R3 + THETA - GRAD2
 
      else

      CL = float(npef)
      L = 1
      U = npef

   10 cxx = (U + float(L))/2.
      I = INT (cxx + 0.5)

      if (R3 .le. efnpfd(I) .and. R3 .ge. efnpfd(I - 1)) then

      Y1 = efnpvl(I - 1)
      Y2 = efnpvl(I)
      XX1 = efnpfd(I - 1)
      XX2 = efnpfd(I)
 
      EFN = (R3 - XX1) * (Y2 - Y1) / (XX2 - XX1) + Y1

      else

      if (R3 .GT. efnpfd(I)) L = I
      if (R3 .LT. efnpfd(I)) U = I
      goto 10
   
      endif
      endif

      EFN = amax1( 0.0, EFN )

      return
      end





*     provide a random sample from a non-parametric distribution for effluent quality                 

      subroutine ecnpar (ECN) 

*     ecnpvl    is an array of class marks,                           
*     ecnpfd    is an array of cumulative frequencies,      
*     npec      is the number of classes,
*     R4        is the uniform variate,                     
*     ECN       is the returned variate.                      

      include 'bacom1cp.for'

*     prepare to deal with values beyond the first and the last values
*     entered as data ....
*     we shall interpolate between the last point and a point half the last
*     interval beyond the last point

*     calculate the first interval
      CLINT1 = ecnpvl (2) - ecnpvl (1)

*     TAU is the minimum - the first value less half the interval
      TAU = ecnpvl (1) - CLINT1 / 2.

*     the last interval 
      CLINT2 = ecnpvl (npec) - ecnpvl (npec-1)

*     THETA is the maximum - last value plus half the interval
      THETA = ecnpvl (npec) + CLINT2 / 2.

*     calculate the gradients of these tails of the distribution
      GRAD1 = ( ecnpvl(1)      - TAU  ) /       ecnpfd (1)
      GRAD2 = ( THETA - ecnpvl (npec) ) / (1. - ecnpfd (npec) )

*     locate this point on the cumulative frequency distribution.
      if ( R4 .le. ecnpfd(1) ) then
      ECN = GRAD1 * R4 + TAU
      else if ( R4 .ge. ecnpfd(npec) ) then
      ECN = GRAD2 * R4 + THETA - GRAD2
      else

      CL = float(npec)
      L = 1
      U = npec
   10 cxx = (U + float(L))/2.
      I = INT (cxx + 0.5)

      if (R4 .le. ecnpfd(I) .and. R4 .ge. ecnpfd(I - 1)) then
      Y1 = ecnpvl(I - 1)
      Y2 = ecnpvl(I)
      XX1 = ecnpfd(I - 1)
      XX2 = ecnpfd(I)
      ECN = (R4 - XX1) * (Y2 - Y1) / (XX2 - XX1) + Y1

      else
      if (R4 .GT. ecnpfd(I)) L = I
      if (R4 .LT. ecnpfd(I)) U = I
      goto 10
   
      end if
      end if

      ECN = amax1( 0.0, ECN )

      return
      end



*     Store the values of all the variables input for the first calculation ...
*     The original values can then always replace values changed for repeat runs ...

      subroutine store the starting data 
      include 'bacom1cp.for'
   
      newrun = 0

      tecmx = ecm ! ##################################
      tecsx = ecs

      xgrgm = grgm
      xgrgs = grgs
      xrgm = rgm
      xrgs = rgs

      xcofc1 = cofc1
      xcoff2 = coff2
      xcofc5 = cofc5

      xcapf = capf
      xcacf = cacf
      xcacp = cacp
      xcafd = cafd
      xcaed = caed
      
      xrfm = rfm
      xrf5 = rf5
*     xreg = reg
*     xadd = add
      xfsh = fsh
*     ishhx = ifsh
      xrcm = rcm
      xrcs = rcs
      xefm = efm
      xefs = efs
      xecm = ecm
      xecs = ecs
      xecv = ecs / ecm
      xtarg = targ

*     freex = free

      xphm = phm
      xphs = phs
      xcam = cam
      xcas = cas
      xudocm = udocm
      xudocs = udocs
      xedocm = edocm
      xedocs = edocs

      xgrcm = grcm
      xgrcs = grcs
      xgecm = gecm
      xgecs = gecs
      xgphm = gphm
      xgphs = gphs
      xgcam = gcam
      xgcas = gcas
      xgudocm = gudocm
      xgudocs = gudocs
      xgedocm = gedocm
      xgedocs = gedocs
      
      return
      end

      subroutine reset the starting data
      include 'bacom1cp.for'

*     iregq = iregqX
      grgm = xgrgm
      grgs = xgrgs
      rgm = xrgm
      rgs = xrgs

      COFC1 = xCOFC1
      COFf2 = xCOFf2
      COFc3 = xCOFc3
      COCf4 = xCOCf4
      COfc5 = xCOfc5
      COCc6 = xCOCc6

      capf = xcapf
      cacf = xcacf
      cacp = xcacp
      cafd = xcafd
      caed = xcaed

      rfm = xrfm
      rf5 = xrf5
*     reg = xreg
*     add = xadd
      fsh = xfsh
*     ifsh = ifshx
      rcm = xrcm
      rcs = xrcs

      phm = xphm
      phs = xphs
      cam = xcam
      cas = xcas
      udocm = xudocm
      udocs = xudocs
      edocm = xedocm
      edocs = xedocs

      grcm = xgrcm
      grcs = xgrcs
      gecm = xgecm
      gecs = xgecs
      gphm = xgphm
      gphs = xgphs
      gcam = xgcam
      gcas = xgcas
      gudocm = xgudocm
      gudocs = xgudocs
      gedocm = xgedocm
      gedocs = xgedocs

      efm = xefm
      efs = xefs
      ecm = xecm
      ecs = xecs
      ecv = ecs / ecm

      targ = xtarg

*      XPER = XPERX
*      PER = PERX

      return
      end

      
      subroutine write the heading
      include 'bacom1cp.for'
      character *1 a1(2)
      character *2 bh,bmin,bday,bmon
            
      open(03,file = "Background Results\MPER.out")
      open(09,file = "Background Results\MPER.mon")
      open(07,file = "Background Results\biometal.sam")
      open(17,file = "Background Results\biosumm.sm1")
      open(27,file = "Background Results\biosamps.sm2")
      open(06,file = "Background Results\biometal.eff")

      call GETDAT ( IYR, IMON, IDAY ) ! get the date
      call GETTIM ( IHR, IMIN, ISEC, IHUN ) ! get the time

      write(bh,2)IHR
    2 format(i2)
      read(bh,1)a1(1),a1(2)
    1 format(2a1)
      if ( a1(1) .eq. ' ' ) a1(1) = '0'
      write(bh,1)a1(1),a1(2)
 
      write(bmin,2)IMIN
      read(bmin,1)a1(1),a1(2)
      if ( a1(1) .eq. ' ' ) a1(1) = '0'
      write(bmin,1)a1(1),a1(2)

      write(bday,2)iday
      read(bday,1)a1(1),a1(2)
      if ( a1(1) .eq. ' ' ) a1(1) = '0'
      write(bday,1)a1(1),a1(2)

      write(bmon,2)imon
      read(bmon,1)a1(1),a1(2)
      if ( a1(1) .eq. ' ' ) a1(1) = '0'
      write(bmon,1)a1(1),a1(2)
      
      write(3,3)bday,bmon,IYR,BH,bmin!,name of discharge,name of metal,
 !    &name of river
    3 format(68('=')/'Mass Balance Calculation ... BIO-METAL',14x,
     &'Date: ',a2,'/',a2,'/',I4/
     &'VERSION 4.4 (Tony Warn 27/12/16)',25x,'Time: ',a2,'.',a2)!68('-'))
!     &10x,'Discharge: ',a30/
!     &10x,'Pollutant: ',a30/
!     &10x,'River:     ',a30)

      return
      end

      
      subroutine set up data (vrfm,vrf5,vrcm,vrcs,vrcns,
     &vefm,vefs,vecm,vecs,vecns,vforw,vtarg,vxper,vtype,
     &vcofc1,vcoff2,vcofc3,vcocf4,vcofc5,vcocc6,vsens,
     &vphm,vphs,vphns,vcam,vcas,vcans,
     &vudocm,vudocs,vudocns,vedocm,vedocs,vedocns,
     &vcapf,vcacf,vcacp,vcafd,vcaed,vsensa,vfree)
      
      include 'bacom1cp.for'
      
      integer vforw,viregq,vtype,vfree,vsens
      integer vrcns,vecns,vphns,vcans,vudocns,vedocns,vsensa
      character*150 vrcfile,vrffile,veffile,vecfile,vidfile
     
      rfm = vrfm !
      rf5 = vrf5 !
      rcm = vrcm !
      rcs = vrcs !
      efm = vefm !
      efs = vefs !
      ecm = vecm !
      ecs = vecs !
      ecv = ecs/ecm !
      forw = vforw !
      targ = vtarg !
      xper = vxper !
      sens = vsens !
      add = vadd !
      reg = vreg !
      rgm = vrgm !
      rgs = vrgs !

      rcns = vrcns
      ecns = vecns
      tcns = (ecm*efm*ecns+rcm*rfm*rcns)/(ecm*efm+rcm*rfm)
     
      fsh = vfsh !
      iregq = viregq !
     
      cofc1 = vcofc1 !
      coff2 = vcoff2 !
      cofc5 = vcofc5 !

      type = vtype

      capf = vcapf
      cacf = vcacf
      cacp = vcacp
      cafd = vcafd
      caed = vcaed

      sensa = vsensa
      
      phm = vphm ! bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
      phs = vphs
      phns = vphns
      cam = vcam
      cas = vcas
      cans = vcans
      
      udocm = vudocm
      udocs = vudocs ! bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
      udocns = vudocns
      edocm = vedocm
      edocs = vedocs ! bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
      edocns = vedocns
      
      free = vfree 

      rfpd = 2
      efpd = 2
      rcpd = 2
      ecpd = 2
      rdpd = 2
      edpd = 2
      phpd = 1 ! normal distribution for pH
      capd = 2
      
      rft = 0.0
      rct = 0.0
      eft = 0.0
      ect = 0.0
      rdt = 0.0
      edt = 0.0
      pht = 0.0
      cat = 0.0

*     open(unit=4,file="FILES.TMP",status="UNKNOWN")
*	read(4,40)vrffile
*	read(4,40)vrcfile
*	read(4,40)veffile
*	read(4,40)vecfile
*	read(4,40)vidfile
*  40 format(A150)
*	close(4)
*     vtype = 0: MC or NP with a percentile target
*     vtype = 1: MC or NP with a mean target
*     vtype = 2: IN with a percentile target
*     vtype = 3: UT with a percentile target (convert to 0, with iutflag = 1)

      nprc = 0
      nprf = 0
      npef = 0
      npec = 0

      rffile=vrffile
      rcfile=vrcfile
      effile=veffile
      ecfile=vecfile

      RS = -2379
      RA = -7765
      RB = -8322 ! bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb

      NS = 2000
      
      nsamp = 0
      
      if ( NS .lt. 12 ) then ! nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn
      read(rffile(1:1),'(a1)')zed
      if ( zed .ne. " " ) then
      call rfnprd
      endif
      read(rcfile(1:1),'(a1)')zed
      if ( zed .ne. " " ) then
      call rcnprd
      endif
      read(effile(1:1),'(a1)')zed
      if ( zed .ne. " ") then
      call efnprd
      endif
      read(ecfile(1:1),'(a1)')zed
      if ( zed .ne. " ") then
      call ecnprd
      endif
      if (msg4 .gt. 0) then
      vmsg4 = msg4 
      return
      endif  
      endif ! nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn
      
      x1 = 0.0
      x2 = 0.0
      x3 = 0.0
      x4 = 0.0
      x5 = 0.0
      x6 = 0.0
      x7 = 0.0
      x8 = 0.0
      x9 = 0.0
      x10 = 0.0
      x11 = 0.0
      x12 = 0.0
      x13 = 0.0
      x14 = 0.0
      x15 = 0.0

      xa1 = 0.0
      xa2 = 0.0
      xa3 = 0.0
      xa4 = 0.0
      xa5 = 0.0
      xa6 = 0.0
      xa7 = 0.0
      xa8 = 0.0
      xa9 = 0.0
      xa10 = 0.0
      xa11 = 0.0
      xa12 = 0.0
      xa13 = 0.0
      xa14 = 0.0
      xa15 = 0.0
      xa16 = 0.0
      xa17 = 0.0

      return
      end

      
      subroutine zero flow (vsens,vsensa,vforw,vrcml,vrcmu,vrcsl,vrcsu,
     &vrcq,vrcql,vrcqu,vtcm,vtcml,vtcmu,vtcs,vtcsl,vtcsu,
     &vtcq,vtcql,vtcqu,vtcns,
     &vecm,vecml,vecmu,vecs,vecsl,vecsu,
     &vec95,vec95l,vec95u,vec99,vec99l,vec99u,
     &vec995,vec995l,vec995u,vtecm,vtecml,vtecmu,
     &vtecs,vtecsl,vtecsu,vte95,vte95l,vte95u,vte99,vte99l,vte99u,
     &vte995,vte995l,vte995u)   
      
      include 'bacom1cp.for'
      
      vsens = 0
      vsensa = 0
      vforw = 1
      call CONAV1 (rcm,rcs,rcns,rcml,rcmu)
      call CONSDEV1 (rcs,rcns,rcsl,rcsu)
      call LNCL1 (rcm,rcs,rcns,xper,rcq,rcql,rcqu,0) 
      
      vrcml = rcml
      vrcmu = rcmu
      vrcsl = rcsl
      vrcsu = rcsu
      vrcq = rcq
      vrcql = rcql
      vrcqu = rcqu
      
      call CONAV1 (ecm,ecs,ecns,ecml,ecmu)
      call CONSDEV1 (ecs,ecns,ecsl,ecsu)
      call LNCL1 (ecm,ecs,ecns,95.0,ec95,ec95l,ec95u,0) 
      call LNCL1 (ecm,ecs,ecns,99.0,ec99,ec99l,ec99u,0) 
      call LNCL1 (ecm,ecs,ecns,99.5,ec995,ec995l,ec995u,0) 
      
      call CONAV1 (edocm,edocs,edocns,edocml,edocmu)
      call CONSDEV1 (edocs,edocns,edocsl,edocsu)
      !call LNCL1 (ecm,ecs,ecns,95.0,ec95,ec95l,ec95u,0) 
      !call LNCL1 (ecm,ecs,ecns,99.0,ec99,ec99l,ec99u,0) 
      !call LNCL1 (ecm,ecs,ecns,99.5,ec995,ec995l,ec995u,0) 

      vtcm = rcm
      vtcml = rcml
      vtcmu = rcmu
      vtcs = rcs
      vtcsl = rcsl
      vtcsu = rcsu
      vtcq = rcq
      vtcql = rcql
      vtcqu = rcqu
      vtcns = rcns

      vecm = ecm
      vecml = ecml
      vecmu = ecmu
      vecs = ecs
      vecsl = ecsl
      vecsu = ecsu
      vec95 = ec95
      vec95l = ec95l
      vec95u = ec95u
      vec99 = ec99
      vec99l = ec99l
      vec99u = ec99u
      vec995 = ec995
      vec995l = ec995l
      vec995u = ec995u
      
      vtecm = ecm
      vtecml = ecml
      vtecmu = ecmu
      vtecs = ecs
      vtecsl = ecsl
      vtecsu = ecsu
      vte95 = ec95
      vte95l = ec95l
      vte95u = ec95u
      vte99 = ec99
      vte99l = ec99l
      vte99u = ec99u
      vte995 = ec995
      vte995l = ec995l
      vte995u = ec995u
      
      close(03)
      close(09,status = "delete")
      close(07,status = "delete")

      return
      end

      
      subroutine take the samples
      include 'bacom1cp.for'
      
      rcm = 0.0
      rcs = 0.0
      rcp = 0.0
      ecm = 0.0
      ecs = 0.0
      ecp = 0.0
      phm = 0.0
      phs = 0.0
      php = 0.0
      cam = 0.0
      cas = 0.0
      cap = 0.0
      udocm = 0.0
      udocs = 0.0
      udocp = 0.0
      edocm = 0.0
      edocs = 0.0
      edocp = 0.0
      
      KS = max0 (int(rcns),int(ecns))
      KS = max0 (KS,int(phns))
      KS = max0 (KS,int(cans))
      KS = max0 (KS,int(udocns))
      KS = max0 (KS,int(edocns))

      do I = 1, KS ! loop through the Monte Carlo shots ----------------------
 
      call get a set of correlated random normal deviates 2
      call get the shot for river quality 2 (I,RC)
      call get the shot for discharge quality 2 (I,EC)
      call get the shot for upstream DOC 2 (I,RD)
      call get the shot for discharge DOC 2 (I,ED)
      call get the shot for downstream pH 2 (I,PH)
      call get the shot for downstream calcium 2 (I,CA)

      if ( i .le. rcns ) then
      rcm = rcm + rc
      rcs = rcs + rc*rc
      endif
      if ( i .le. ecns ) then
      ecm = ecm + ec
      ecs = ecs + ec*ec
      endif
      
      if ( i .le. phns ) then
      phm = phm + ph
      phs = phs + ph*ph
      endif
      if ( i .le. cans ) then
      cam = cam + ca
      cas = cas + ca*ca
      endif
      if ( i .le. udocns ) then
      udocm = udocm + rd
      udocs = udocs + rd*rd
      endif
      if ( i .le. edocns ) then
      edocm = edocm + ed
      edocs = edocs + ed*ed
      endif
      
      if ( I .le. int(rcns)) then
      if ( I .eq. 0 ) write(7,5025)
 5025 format(/100('-')/
     &3x,6x,'RC'6x,'EC',6x,'pH',6x,'Ca',4x,'uDOC',4x,'eDOC'/
     &100('-'))
!      write(7,5032)I,RC,EC,PH,CA,RD,ED,r4,gecm,gecs,EC,ecm/float(i)
!      write(9,5032)I,RC,rcm,xrcm
 5032 format(i3,' rc ... ',6f8.2)
      endif

      enddo
           
*     compute mean and standard deviation --------------------------------------
      XXX = rcs - rcm*rcm/rcns
      if (XXX .le. Small) then 
      rcs = 0.0
      else
      rcs = SQRT (XXX/(rcns-1.0))
      endif
      rcm = rcm / rcns
      !rcs = 1.09 * rcs
*     --------------------------------------------------------------------------
      XXX = ecs - ecm*ecm/ecns
      if (XXX .le. Small) then 
      ecs = 0.0
      else
      ecs = SQRT (XXX/(ecns-1.0))
      endif
      ecm = ecm / ecns
      !ecs = 1.04 * ecs
*     --------------------------------------------------------------------------
      XXX = phs - phm*phm/phns
      if (XXX .le. Small) then 
      phs = 0.0
      else
      phs = SQRT (XXX/(phns-1.0))
      endif
      phm = phm / phns
      !phs = 1.012 * phs
*     --------------------------------------------------------------------------
      XXX = cas - cam*cam/cans
      if (XXX .le. Small) then 
      cas = 0.0
      else
      cas = SQRT (XXX/(cans-1.0))
      endif
      cam = cam / cans
      !cas = 1.02 * cas
*     --------------------------------------------------------------------------
      XXX = udocs - udocm*udocm/udocns
      if (XXX .le. Small) then 
      udocs = 0.0
      else
      udocs = SQRT (XXX/(udocns-1.0))
      endif
      udocm = udocm / udocns
      !udocs = 1.045 * udocs
*     --------------------------------------------------------------------------
      XXX = edocs - edocm*edocm/edocns
      if (XXX .le. Small) then 
      edocs = 0.0
      else
      edocs = SQRT (XXX/(edocns-1.0))
      endif
      edocm = edocm / edocns
      !edocs = 1.04 * edocs
*     compute mean and standard deviation --------------------------------------
    
      write(27,5010)rcm,rcs,xrcm,xrcs,rcs/rcm,xrcs/xrcm,
     &100.0*(rcs/rcm)/(xrcs/xrcm)-100.0,BM(2),BS(2)      
 5010 format(2f8.2,' --- ',2f8.2,2f10.4,'  rcm  ',f9.1,' %',2f7.2) 
      rcp = rcp + 100.0*(rcs/rcm)/(xrcs/xrcm)-100.0
      
      write(27,5015)ecm,ecs,xecm,xecs,ecs/ecm,xecs/xecm,
     &100.0*(ecs/ecm)/(xecs/xecm)-100.0
 5015 format(2f8.2,' --- ',2f8.2,2f10.4,'  ecm  ',f9.1,' %')
      ecp = ecp + 100.0*(ecs/ecm)/(xecs/xecm)-100.0
      write(27,5132)
     
      write(27,5011)phm,phs,xphm,xphs,phs/phm,xphs/xphm,
     &100.0*(phs/phm)/(xphs/xphm)-100.0 
 5011 format(2f8.2,' --- ',2f8.2,2f10.4,'  phm  ',f9.1,' %')
      php = php + 100.0*(phs/phm)/(xphs/xphm)-100.0
      
      write(27,5012)cam,cas,xcam,xcas,cas/cam,xcas/xcam,
     &100.0*(cas/cam)/(xcas/xcam)-100.0
 5012 format(2f8.2,' --- ',2f8.2,2f10.4,'  cam  ',f9.1,' %')
      cap = cap + 100.0*(cas/cam)/(xcas/xcam)-100.0
      
      write(27,5013)udocm,udocs,xudocm,xudocs,udocs/udocm,
     &xudocs/xudocm,100.0*(udocs/udocm)/(xudocs/xudocm)-100.0     
 5013 format(2f8.2,' --- ',2f8.2,2f10.4,'  udocm',f9.1,' %')
      udocp = udocp + 100.0*(udocs/udocm)/(xudocs/xudocm)-100.0
     
      write(27,5014)edocm,edocs,xedocm,xedocs,edocs/edocm,
     &xedocs/xedocm,100.0*(edocs/edocm)/(xedocs/xedocm)-100.0
 5014 format(2f8.2,' --- ',2f8.2,2f10.4,'  edocm',f9.1,' %')
      edocp = edocp + 100.0*(edocs/edocm)/(xedocs/xedocm)-100.0 

      write(27,5132)
 5132 format(100('-'))

      return
      end
      
   
      
*     --------------------------------------------------------------------------
*     eliminate bias 
      function BIASQ (V,BM,BS,RM)
      Xbias = ( BM * V - RM ) * BS + RM
      if ( Xbias .lt. (0.00001*RM) ) Xbias = 0.0 ! 0.00001*RM
      if ( V .gt. 0.0 ) then
      if ( Xbias .gt. (2.0 * V) ) Xbias = V
      endif
      BIASQ = Xbias
      return
      end
*     --------------------------------------------------------------------------
      function BIASF (V, BM, B5, R5 ) ! log normal
      if ( V .gt. 1.2 * R5 ) V = 1.0004 * BM * V
      if ( V .lt. 1.2 * R5 ) V = B5 * V
      BIASF = V
      return
      end
*     --------------------------------------------------------------------------
