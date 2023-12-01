*     MASS-BALANCE CALCULATION: MONTE-CARLO SIMULATION FOR METALS --------------
*     M-PER VERSION 4.4 ... 27/12/16 ... Written by: Tony Warn -----------------
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine massb (ic)
      include 'bacom1cp.for'
      
*     initialise accumulators for calculating mean and standard deviation ------
      udm = 0.0 ! mean dissolved metal in upstream river 
      uds = 0.0 ! standard deviation
      tem = 0.0 ! mean dissolved metal in the discharge ------------------------
      tes = 0.0 ! standard deviation
      tdm = 0.0 ! mean dissolved metal in the river downstream of the discharge
      tds = 0.0 ! standard deviation
      ddocm = 0.0 ! mean dissolved organic carbon in the downstream river
      ddocs = 0.0 ! standard deviation
      ddocns = 0.0 ! number of samples
      ubm = 0.0 ! mean bioavailabe metal in the upstream river
      ubs = 0.0 ! standard deviation
      rbm = 0.0 ! mean number of samples for metal in the downstream river =====
      rbs = 0.0 ! standard deviation
      cbm = 0.0 ! mean number of samples for metal in the downstream river =====
      cbs = 0.0 ! standard deviation
      tbm = 0.0 ! mean bioavailable metal in the downstream river  
      tbs = 0.0 ! standard deviation
      vbm = 0.0 ! mean bioavailabe metal in the upstream river (with d/s DOC)
      vbs = 0.0 ! standard deviation
      
      tbxx = 0.0

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

      if (tecm .lt. Small) then ! set up the log domain for discharge quality -- 
      gecm = - Big ! mean in the log domain
      gecs = 0.0 ! standard deviation in the log domain
      else ! set mean etc in the log domain
      gecm = ALOG( (tecm*tecm)/SQRT(tecm*tecm + tecs*tecs) ) ! mean (log domain)
      gecs = SQRT( ALOG(1.0 + (tecs*tecs)/(tecm*tecm)) ) ! standard deviation --
      endif ! set up mean and standard deviation in the log domain -------------

      do 2 I = 1, NS ! loop through the Monte Carlo shots ----------------------
      call get a set of correlated random normal deviates
      call get the shot for river flow (I,RF)
      call get the shot for river quality (I,RC)
      call get the shot for discharge flow (I,EF)
      call get the shot for discharge quality (I,EC)
      C(2,I) = EC ! store the initial dissolved metal in the discharge
      C(3,I) = RC ! store the upstream dissolved metal
      C(5,I) = RF ! store the upstream river flow
      C(6,I) = EF ! discharge flow
      if ( free .eq. 1 ) then
      call get the shot for upstream DOC (I,RD)
      call get the shot for discharge DOC (I,ED)
      call get the shot for downstream pH (I,PH)
      call get the shot for downstream calcium (I,CA)
      C(4,I) = RD ! initially set downstream DOC to the upstream DOC
      C(9,I) = ED ! discharge DOC
      C(10,I) = RD ! upstream DOC
      C(13,I) = CA ! upstream and downstream calcium
      endif

*     store the data for tests based on regression -----------------------------   
      CR(2,I) = RF ! store the upstream river flow
      CR(3,I) = RC ! upstream dissolved metal
      CR(4,I) = EF ! discharge flow
      CR(5,I) = EC ! discharge dissolved metal
      CR(6,I) = ED ! dissolved organic carbon in the discharge
      CR(7,I) = PH ! river pH
      CR(8,I) = CA ! river calcium
      CR(9,I) = RD ! upstream dissolved organic carbon -------------------------
      
      if (EF .lt. Small) then ! if the discharge flow shot is zero -------------
      C(1,I) = RC ! set downstream dissolved metal to the upstream value -------
      C(14,I) = rcns ! set downstream sampling rate to the upstream value ------
      C(15,I) = rcns ! set downstream sampling rate to the upstream value ------
      C(4,I) = RD ! set the downstream DOC to the upstream DOC -----------------
      goto 1 ! if the discharge flow is zero -----------------------------------
      endif ! if the discharge flow is zero ------------------------------------

      if (RF .lt. Small) then ! if the river flow shot is zero -----------------
      C(1,I) = EC ! set the downstream quality to the discharge quality --------
      C(14,I) = ecns ! set downstream sampling rate to the upstream value ------
      C(15,I) = ecns ! set downstream sampling rate to the upstream value ------
      C(4,I) = ED ! set d/s DOC to the discharge DOC -------
      goto 1 ! if the river flow is zero ---------------------------------------
      endif ! if the river flow is zero ----------------------------------------
      
*     mass balance for dissolved metal downstream of the discharge -------------
      C(1,I) = ( (RF*RC) + (EF*EC) ) / (RF + EF)
      
      C(14,I) = ((RF*RC*rcns) + (EF*EC*ecns)) / (RF*RC + EF*EC) !  sampling rate
      
      if ( i .lt. 11 ) then
      !write(03,6192)I,RF,RC,EF,EC,C(14,I)
 6192 format(i4,4f8.3,f12.4)
      endif

      if ( iter .eq. 1 .or. forw .eq. 1 ) then
      C(15,I) = ((RF*RC*rcns) + (EF*EC*ecns)) / (RF*RC + EF*EC) !  sampling rate
      endif
      
*     and for dissolved organic carbon downstream of the discharge -------------
      if ( free .eq. 1 ) then
      !if ( nodoc .eq. 0 ) C(4,I) = ( (RF*RD) + (EF*ED) ) / (RF + EF)
      C(4,I) = ( (RF*RD) + (EF*ED) ) / (RF + EF)
      endif
    1 continue ! continue the loop on shots ====================================

*     calculate the cut-off in the calculations for copper #####################
      if ( free .eq. 1 ) then
      if ( mettal .eq. 2 ) then ! do this only for copper CCCCCCCCCCCCCCCCCCCCCC
      ttca = CR(8,I) ! river calcium for this shot ! CCCCCCCCCCCCCCCCCCCCCCCCCCC
      ttph = CR(7,I) ! river pH for this shot ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      tdc = 4.0 ! set starting point for the DOC ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call bioav for copper (pnecb,biof,ttca,ttph,tdc,99,1) ! cut-off sums CCCCC
      cutoff = biof ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do jjj = 1,500 ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      tdc = tdc + 0.1 ! increase DOC by 0.1 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call bioav for copper (pnecb,biof,ttca,ttph,tdc,99,2) ! cut-off sums CCCCC
      if ( biof .gt. cutoff ) goto 3 ! identify biof less than cut-off CCCCCCCCC
      cutoff = biof ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      enddo ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    3 continue ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      cutoff = tdc ! set cut-off ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      else ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      cutoff = 999999.9 ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      endif ! for copper only ! ################################################
      
      CR(10,I) = C(4,I) ! store downstream DOC for use in regression ===========
      CR(01,I) = C(1,I) ! store dissolved metal for use in regression ==========

      if ( biometal .eq. 2 ) then ! +++++++++++++++ calculate bioavailable metal
      if ( mettal .eq. 1 ) then ! zinc ==================================== zinc
      call bioav for zinc (C(8,i),CA,pH,RD) ! bioavailable zinc u/s of discharge
      C(12,i) = C(8,i) * C(3,i) ! bioavailable zinc u/s of discharge
      call bioav for zinc (C(7,i),CA,pH,C(4,I)) ! bioavailable zinc u/s === zinc 
      C(11,i) = C(7,i) * C(1,i) ! bioavailable zinc d/s of discharge
      call bioav for zinc (C(8,i),CA,pH,C(4,I)) ! bioavailable zinc u/s === zinc 
      C(16,i) = C(8,i) * C(3,i) ! bioavailable zinc d/s of discharge
      endif ! zinc ======================================================== zinc
      
      if ( mettal .eq. 2 ) then ! copper ================================ copper
      RDXu = amin1 ( cutoff, RD ) ! apply the cut-off for u/s DOC --------------
      call bioav for copper (pnecb,C(8,i),CA,pH,RDXu,ic,1) ! u/s bio is C(8,i) =
      C(12,i) = C(8,i) * C(3,i) ! bioavailable copper u/s of discharge ---------
      RDXd = amin1 ( cutoff, C(4,I) ) ! apply the cut-off for d/s DOC ----------
      call bioav for copper (pnecb,C(7,i),CA,pH,RDXd,ic,2) ! d/s bio is C(7,i) =
      C(11,i) = C(7,i) * C(1,i) ! copper d/s of discharge ---------
      call bioav for copper (pnecb,C(8,i),CA,pH,RDXd,ic,2) ! d/s bio is C(7,i) =
      C(16,i) = C(8,i) * C(3,i) ! copper d/s of discharge ---------
      endif ! copper ==================================================== copper
      
      
      if ( mettal .eq. 3 ) then ! manganese ========================== manganese
      call bioav for manganese (C(8,i),CA,pH,RD,jj) ! manganese u/s of discharge
      C(12,i) = C(8,i) * C(3,i) ! =             = bio manganese u/s of discharge
      call bioav for manganese (C(7,i),CA,pH,C(4,I),jj) ! manganese d/s ========
      C(11,i) = C(7,i) * C(1,i) ! =============== bio manganese d/s of discharge
      call bioav for manganese (C(8,i),CA,pH,C(4,I),jj) ! manganese d/s ========
      C(16,i) = C(8,i) * C(3,i) ! =============== bio manganese d/s of discharge
      endif ! manganese ============================================== manganese
    
      if ( mettal .eq. 4 ) then ! nickel ================================ nickel
      call bioav for nickel (C(8,i),CA,pH,RD)
      C(12,i) = C(8,i) * C(3,i) ! bioavailable nickel u/s of discharge
      call bioav for nickel (C(7,i),CA,pH,C(4,I))
      C(11,i) = C(7,i) * C(1,i) ! bioavailable nickel u/s of discharge
      call bioav for nickel (C(8,i),CA,pH,C(4,I))
      C(16,i) = C(8,i) * C(3,i) ! bioavailable nickel u/s of discharge
      endif ! nickel ==================================================== nickel
      
      if ( mettal .eq. 5 ) then ! lead =====================================lead
      call bioav for lead (C(8,i),RD)
      C(12,i) = C(8,i) * C(3,i) ! bioavailable lead u/s of discharge
      call bioav for lead (C(7,i),C(4,I))
      C(11,i) = C(7,i) * C(1,i) ! bioavailable lead u/s of discharge
      call bioav for lead (C(8,i),C(4,I))
      C(16,i) = C(8,i) * C(3,i) ! bioavailable lead u/s of discharge
      endif ! lead =========================================================lead
      endif ! end of +++++++++++++++++++++++++++++++++ bioavailable calculations
      
      CR(11,I) = C(11,I) ! store bioavailable metal for use in regression ======
      endif ! if ( free .eq. 1 ) ===============================================
   
*     accumulate for the calculation of mean and standard deviation ------------
      tdm = tdm + C(1,I) ! mean dissolved metal in the downstream river --------
      tds = tds + C(1,I) * C(1,I) ! standard deviation in the downstream river -
      tem = tem + C(2,I) ! mean dissolved metal in the discharge ---------------
      tes = tes + C(2,I) * C(2,I) ! standard deviation in the discharge --------
      udm = udm + C(3,I) ! mean dissolved metal in the upstream river ----------
      uds = uds + C(3,I) * C(3,I) ! standard deviation in the upstream river ---
      
      rbm = rbm + C(14,I) ! mean number of samples for metal in the d/s river --
      rbs = rbs + C(14,I) * C(14,I) ! standard deviation for rbm --------------- 
      cbm = cbm + C(15,I) ! mean number of samples for metal in the d/s river --
      cbs = cbs + C(15,I) * C(15,I) ! standard deviation for cbm --------------- 
      
      if ( biometal .eq. 2 ) then ! for biovailable metal
      ddocm = ddocm + C(4,I) ! mean dissolved organic carbon in the d/s river --
      ddocs = ddocs + C(4,I) * C(4,I) ! standard deviation in the d/s river ----
      tbm = tbm + C(11,I) ! bioavailable metal downstream of the discharge -----
      tbs = tbs + C(11,I) * C(11,I) ! d/s standard deviation of bio-metal ------    
      ubm = ubm + C(12,I) ! mean bioavailable metal upstream of the discharge --
      ubs = ubs + C(12,I) * C(12,I) ! u/s standard deviation of bio-metal ------
      vbm = vbm + C(16,I) ! mean bioavailable metal upstream of the discharge --
      vbs = ubs + C(16,I) * C(16,I) ! u/s standard deviation of bio-metal ------
      endif ! for biovailable metal --------------------------------------------
    2 continue ! loop on the Monte Carlo shots ---------------------------------

*     compute the statistics ---------------------------------------------------
      call compute double precision summary statistics (tem,tes) ! diss. metal in discharge -
      call compute double precision summary statistics (tdm,tds) ! d/s dissolved metal
      call compute ordinary summary statistics (ddocm,ddocs) ! downstream DOC --
      call compute double precision summary statistics (udm,uds) ! u/s DOC ----- 
      goto 5555
*     compute mean and standard deviation of u/s DOC ---------------------------
      XXX = uds - udm*udm/FLOAT(NS)
      if (XXX .le. Small) then 
      uds = 0.0
      else
      uds = SQRT (XXX/FLOAT(NS-1))
      endif
      udm = udm / FLOAT(NS)
*     compute mean and standard deviation of u/s DOC ---------------------------
 5555 continue
      
*     calculate d/s sampling rate for DOC --------------------------------------
      ddocns =((udocns*RFM*udm)+(edocns*EFM*ddocm ))/(RFM*udm+EFM*ddocm)
*     calculate d/s sampling rate for dissolved metal --------------------------
      tdns = ((rcns*RFM*rcm)+(ecns*EFM*tem))/(RFM*rcm+EFM*tem )
*     set d/s sampling rate for bioavailable metal -----------------------------
      tbns = tdns

*     call compute double precision summary statistics (udm,uds) ! u/s DOC 
      call compute double precision summary statistics (tbm,tbs) ! d/s bio-metal
      call compute double precision summary statistics (ubm,ubs) ! u/s bio-metal 
      call compute double precision summary statistics (vbm,vbs) ! u/s bio-metal 
      
      call compute double precision summary statistics (rbm,rbs) !
      call compute double precision summary statistics (cbm,cbs) ! 
      
*     ############################################################ STORE RESULTS
      C(1,NS + 1) = tdm ! --- store mean dissolved metal downstream of discharge
      C(3,NS + 1) = udm ! ----------- mean dissolved metal upstream of discharge
      if ( ic .ne. 4 ) then
      C(14,NS + 1) = rbm ! -------------------- mean number of samples for metal
      endif
      if ( iter .eq. 1 .or. forw .eq. 1 ) then
      C(15,NS + 1) = cbm ! -------------------- mean number of samples for metal
      endif
      
      if ( biometal .eq. 2 ) then ! -------- if bioavailable metal is being used
      C(4,NS + 1) = ddocm ! - mean dissolved organic carbon d/s of the discharge
      C(11,NS + 1) = tbm ! ----- mean bioavailable metal downstream of discharge
      C(12,NS + 1) = ubm ! ------- mean bioavailable metal upstream of discharge
      endif ! if ( biometal .eq. 2 ) ------- if bioavailable metal is being used
      call sequence the calculated shots ! --------------- calculate percentiles
      

*     ############################################################ STORE RESULTS
      if ( ic .eq. 4 ) then ! --- for calculation of impact of current discharge
      tdmc = tdm ! -------------------- mean downstream of the current discharge
      call compute confidence limits d (0.0,2,tdm,tds,0.0,C(15,NS+1),
     &con1,con2) ! -- d/s of the current discharge
      tdmclc = con1 ! ---------------------- downstream of the current discharge
      tdmcuc = con2 ! ---------------------- downstream of the current discharge
      tdml = con1 ! ------------------------ downstream of the current discharge
      tdmu = con2 ! ------------------------ downstream of the current discharge
      tdsc = tds ! ------ standard deviation downstream of the current discharge
      tdc05 = C(01,K05) ! - 05-percentile dissolved d/s of the current discharge 
      tdc90 = C(01,K90) ! - 90-percentile dissolved d/s of the current discharge 
      tdc95 = C(01,K95) ! - 95-percentile dissolved d/s of the current discharge 
      call compute confidence limits d (95.0,2,tdm,tds,tdc95,C(15,NS+1),
     &con1,con2) ! d/s of the current discharge
      tdc99 = C(01,K99) ! - 99-percentile dissolved d/s of the current discharge
      tdc95lc = con1 ! --------------------- downstream of the current discharge
      tdc95uc = con2 ! --------------------- downstream of the current discharge
      call compute confidence limits d (Purcentile,2,tdm,tds,
     &C(01,KPERC), ! --------------- for target-percentile
     &C(15,NS+1),con1,con2) ! -------------- downstream of the current discharge
      tdc9T = C(01,KPERC) ! - target %ile dissolved d/s of the current discharge
      tdc9Tlc = con1
      tdc9Tuc = con2
      
      if ( biometal .eq. 2 ) then ! -------- if bioavailable metal is being used
      tbmc = tbm ! ------------ mean biometal downsteam of the current discharge 
      tbsc = tbs ! ---------- st.dev biometal downsteam of the current discharge 
      tbc05 = C(11,K05) ! -- 05-percentile biometal d/s of the current discharge 
      tbc90 = C(11,K90) ! -- 90-percentile biometal d/s of the current discharge 
      tbc95 = C(11,K95) ! -- 95-percentile biometal d/s of the current discharge 
      tbc99 = C(11,K99) ! -- 99-percentile biometal d/s of the current discharge 
      endif ! if ( biometal .eq. 2 ) ------- if bioavailable metal is being used
      endif ! if ( ic .eq. 4 ) -- for calculation of impact of current discharge
*     ############################################################ STORE RESULTS

*     ############################################################ STORE RESULTS
      !if ( ic .eq. 3 .or. ic .eq. 5 .or. ic .eq. 6 ) then ! --------------------
          
      td05 = C(01,K05) ! - 05-percentile dissolved d/s of the modified discharge 
      td90 = C(01,K90) ! - 90-percentile dissolved d/s of the modified discharge 
      td95 = C(01,K95) ! - 95-percentile dissolved d/s of the modified discharge 
      td99 = C(01,K99) ! - 99-percentile dissolved d/s of the modified discharge
     
      call compute confidence limits d (95.0,2,tdm,tds,C(1,K95),
     &C(14,K95),con1,con2) ! d/s of the modified discharge
      td95 = C(1,K95)
      td95lc = con1 ! --------------------- downstream of the modified discharge
      td95uc = con2 ! --------------------- downstream of the modified discharge

      if ( type .eq. 0 ) then
      call compute confidence limits d (Purcentile,2,tdm,tds,
     &C(01,KPERC),C(14,KPERC),con1,con2) ! ------- d/s of the modified discharge
      td9T = C(01,KPERC) ! - target %ile dissolved d/s of the modified discharge
      td9Tlc = con1
      td9Tuc = con2
      endif
      
      if ( biometal .eq. 2 ) then ! -------- if bioavailable metal is being used
      tb05 = C(11,K05) ! -- 05-percentile biometal d/s of the modified discharge 
      tb90 = C(11,K90) ! -- 90-percentile biometal d/s of the modified discharge 
      tb95 = C(11,K95) ! -- 95-percentile biometal d/s of the modified discharge 
      tb99 = C(11,K99) ! -- 99-percentile biometal d/s of the modified discharge
      tbxx = C(11,KTG) ! -- 99-percentile biometal d/s of the modified discharge
 
      ub05 = C(12,K05) ! -- 05-percentile biometal u/s of the discharge 
      ub90 = C(12,K90) ! -- 90-percentile biometal u/s of the discharge 
      ub95 = C(12,K95) ! -- 95-percentile biometal u/s of the discharge 
      ub99 = C(12,K99) ! -- 99-percentile biometal u/s of the discharge
      ubxx = C(12,KTG) ! -- 99-percentile biometal u/s of the discharge
      
      ub05 = C(16,K05) ! -- 05-percentile biometal u/s of the discharge 
      ub90 = C(16,K90) ! -- 90-percentile biometal u/s of the discharge 
      ub95 = C(16,K95) ! -- 95-percentile biometal u/s of the discharge 
      ub99 = C(16,K99) ! -- 99-percentile biometal u/s of the discharge
      ubxx = C(16,KTG) ! -- 99-percentile biometal u/s of the discharge

      !write(03,6188)ubm,ubs,ubsplus,rcns,ub95,ubxx
 6188 format(/40('-')/
     &'       ubm =',f10.4,' OK'/
     &'       ubs =',f10.4,' OK'/
     &'  ubs plus =',f10.4,' OK'/
     &'      rcns =',f10.4,' OK'/
     &'      ub95 =',f10.4/
     &'      ubxx =',f10.4/40('-')/)

      call compute confidence limits d(95.0,2,tbm,tbs*tbsplus,C(11,K95),
     &C(14,K95),con1,con2) ! d/s of the modified discharge
      tb95lc = con1 ! --------------------- downstream of the modified discharge
      tb95uc = con2 ! --------------------- downstream of the modified discharge
      call compute confidence limits d(xper,2,tbm,tbs*tbsplus,C(11,KTG),
     &C(14,K95),con1,con2) ! d/s of the modified discharge
      tbxxl = con1 ! --------------------- downstream of the modified discharge
      tbxxu = con2 ! --------------------- downstream of the modified discharge

      call compute confidence limits d(95.0,2,ubm,ubs,C(12,K95),
     &rcns,con1,con2) ! d/s of the  discharge
      tu95lc = con1 ! --------------------- downstream of the discharge
      tu95uc = con2 ! --------------------- downstream of the discharge
      call compute confidence limits d(xper,2,ubm,ubs,C(12,KTG),
     &rcns,con1,con2) ! d/s of the modified discharge
      ubxxl = con1 ! --------------------- downstream of the discharge
      ubxxu = con2 ! --------------------- downstream of the discharge
      
      call compute confidence limits d(95.0,2,vbm,vbs,C(16,K95),
     &rcns,con1,con2) ! d/s of the  discharge
      tv95lc = con1 ! --------------------- downstream of the discharge
      tv95uc = con2 ! --------------------- downstream of the discharge
      call compute confidence limits d(xper,2,vbm,vbs,C(16,KTG),
     &rcns,con1,con2) ! d/s of the modified discharge
      vbxxl = con1 ! --------------------- downstream of the discharge
      vbxxu = con2 ! --------------------- downstream of the discharge

      if ( type .eq. 0 ) then
      call compute confidence limits d (Purcentile,2,tbm,tbs,
     &C(11,KPERC),C(14,KPERC),con1,con2) ! ------- d/s of the modified discharge
      td9T = C(01,KPERC) ! - target %ile dissolved d/s of the modified discharge
      td9Tlc = con1 ! --------------------- downstream of the modified discharge
      td9Tuc = con2 ! --------------------- downstream of the modified discharge
      endif
      
      endif ! if ( biometal .eq. 2 ) ------- if bioavailable metal is being used
      !endif ! if ( ic .eq. 4 ) - for calculation of impact of modified discharge
*     ############################################################ STORE RESULTS

      return
      end


      subroutine get a set of correlated random normal deviates
      include 'bacom1cp.for'

      RR1 = GAS1(IR1) ! river flow
      RR2 = GAS2(IR2) ! upstream dissolved metal 
      RR3 = GAS3(IR3) ! discharge flow
      RR4 = GAS4(IR4) ! discharge - concentration of dissolved metal
      RR5 = GAS5(IR5) ! discharge - concentration of DOC
      RR6 = GAS6(IR6) ! dissolved organic carbon in the upstream river
      RR7 = GAS7(IR7) ! pH in the downstream river
      RR8 = GAS8(IR8) ! calcium in downstream river

      R1 = RR1             ! river flow
      R2 = b1*RR1 + b2*RR2 ! correlate river quality and river flow
      R3 = c1*RR1 + c2*RR3 ! correlate discharge flow and river flow
      R4 = d1*R3 + d2*RR4  ! correlate discharge quality and discharge flow
      
      R6 = f1*RR1 + f2*RR6 ! DOC in the u/s river and river flow
      R5 = e1*R3 + e2*RR5  ! DOC in the discharge and discharge flow
      
      R7 = p1*RR1 + p2*RR7     ! d/s pH - correlate with river flow ------------
      R8 = fca1*RR1 + fca2*RR8 ! d/s calcium - correlate with river flow -------
      R8 = fcp1*RR7 + fcp2*R8  ! d/s calcium - correlate with pH ---------------
      
      return
      end

      subroutine get a set of correlated random normal deviates 2
      include 'bacom1cp.for'

      RR1 = GAS1s(LR1) ! river flow
      RR2 = GAS2s(LR2) ! upstream dissolved metal 
      RR3 = GAS3s(LR3) ! discharge flow
      RR4 = GAS4s(LR4) ! discharge - concentration of dissolved metal
      RR5 = GAS5s(LR5) ! discharge - concentration of DOC
      RR6 = GAS6s(LR6) ! dissolved organic carbon in the upstream river
      RR7 = GAS7s(LR7) ! pH in the downstream river
      RR8 = GAS8s(LR8) ! calcium in downstream river

      R1 = RR1             ! river flow
      R2 = b1*RR1 + b2*RR2 ! correlate river quality and river flow
      R3 = c1*RR1 + c2*RR3 ! correlate discharge flow and river flow
      R4 = d1*R3 + d2*RR4  ! correlate discharge quality and discharge flow
      
      R6 = f1*RR1 + f2*RR6 ! DOC in the u/s river and river flow
      R5 = e1*R3 + e2*RR5  ! DOC in the discharge and discharge flow
      
      R7 = p1*RR1 + p2*RR7     ! d/s pH - correlate with river flow ------------
      R8 = fca1*RR1 + fca2*RR8 ! d/s calcium - correlate with river flow -------
      R8 = fcp1*RR7 + fcp2*R8  ! d/s calcium - correlate with pH ---------------

      return
      end

      
      
      subroutine convert correlation coefficients
      include 'bacom1cp.for'

      b1 = COFC1 ! upstream river flow and upstream river quality --------------
      sq = 1.0 - b1*b1
      if (sq .gt. 0.00001) then
      b2 = SQRT (sq)
      else
      b2 = 0.0
      endif ! upstream river flow and upstream river quality -------------------
      
      c1 = COFf2 ! upstream river flow on discharge flow -----------------------
      sq = 1.0 - c1*c1
      if (sq .gt. 0.00001) then
      c2 = SQRT (sq)
      else
      c2 = 0.0
      endif ! upstream river flow on discharge flow ----------------------------

      d1 = COfc5 ! discharge flow on discharge quality --------------------------
      sq = 1.0 - d1*d1
      if (sq .gt. 0.00001) then
      d2 = SQRT (sq)
      else
      d2 = 0.0
      endif ! upstream river flow on discharge flow ----------------------------

      p1 = capf ! correlation between river flow and pH -----------------------
      sq = 1.0 - p1*p1
      if (sq .gt. 0.00001) then
      p2 = SQRT (sq)
      else
      p2 = 0.0
      endif ! correlation between river flow and pH ----------------------------

      fca1 = cacf ! correlation between river flow and calcium ----------------
      sq = 1.0 - fca1*fca1
      if (sq .gt. 0.00001) then
      fca2 = SQRT (sq)
      else
      fca2 = 0.0
      endif ! correlation between river flow and calcium -----------------------
      
      fcp1 = cacp ! correlation between calcium and pH ------------------------
      sq = 1.0 - fcp1*fcp1
      if (sq .gt. 0.00001) then
      fcp2 = SQRT (sq)
      else
      fcp2 = 0.0
      endif ! correlation between calcium and pH -------------------------------

      f1 = cafd ! upstream river flow and upstream dissolved organic carbon --
      sq = 1.0 - f1*f1
      if (sq .gt. 0.00001) then
      f2 = SQRT (sq)
      else
      f2 = 0.0
      endif ! upstream river flow and upstream dissolved organic carbon --------

      e1 = caed ! discharge flow and dissolved organic carbon ------------------
      sq = 1.0 - e1*e1
      if (sq .gt. 0.00001) then
      e2 = SQRT (sq)
      else
      e2 = 0.0
      endif ! discharge flow and dissolved organic carbon ----------------------

      return 
      end


      subroutine set up the calculations
      include 'bacom1cp.for'
      
      inkbias = 0 ! correct the bias in the equations (test) ===================
      
      bioERRORS = 0 ! include the errors from the bioavailable equation SSSSSSSS 
           
      kkrun = 10 ! number of checks on errors ----------------------------------
      if ( SITEBIAS .eq. 1 ) kkrun = 12 ! do more when SITEBIAS equals 1 FFFFFFF

      NS = 2000  ! the number of Monte-Carlo shots -----------------------------
*     precision factor - the iterative scheme stops when successive trials are -
*     within 100 * CFAC % ------------------------------------------------------
      CFAC = 0.00001

*     array element for obtaining the percentiles from a ranked list -----------
      KTG = INT(0.5 + 0.01 * XPER * NS)
      K95 = INT(0.5 + 0.95 * NS) ! 95-percentile from a ranked list ------------
      K90 = INT(0.5 + 0.90 * NS) ! 90-percentile from a ranked list ------------
      K80 = INT(0.5 + 0.80 * NS) ! 80-percentile from a ranked list ------------
      K99 = INT(0.5 + 0.99 * NS) ! 99-percentile from a ranked list ------------
      K995 = INT(0.5 + 0.995 * NS) ! 99.5-percentile from a ranked list --------
      K999 = INT(0.5 + 0.999 * NS) ! 99.9-percentile from a ranked list --------
      K05 = INT(0.5 + 0.05 * NS) ! 5-percentile from a ranked list -------------
      K50 = INT(0.5 + 0.50 * NS) ! 50-percentile from a ranked list ------------

*     default values of correlation coefficients not read in with data ---------
      COFc3 = 0.0 ! river flow on discharge quality ----------------------------
      COCf4 = 0.0 ! river quality on discharge flow ----------------------------
      COCc6 = 0.0 ! river quality on discharge quality -------------------------
      
      small = 1.0e-15 ! a small number
      big = 1.0e5 ! a big number

      do im = 1, NPD
      npdp(im) = 0
      do i = 1, mprf
      npval(im,i) = 0.0
      npcum(im,i) = 0.0
      enddo
      enddo
      
      return      
      end

      
      
    
      
      
      
*     arrange shots in increasing order ----------------------------------------
*     in order to obtain the percentiles ---------------------------------------
      subroutine sequence the calculated shots
      include 'bacom1cp.for'
      
*     items stored in array C(K,I) ---------------------------------------------
*     1  ... the dissolved metal shots in the downstream river -----------------
*     2  ... the dissolved metal shots in the discharge ------------------------
*     3  ... the dissolved metal shots in the upstream river -------------------
*     4  ... the dissolved organic carbon shots in the downstream river --------
*     5  ... river flow --------------------------------------------------------
*     6  ... discharge flow ----------------------------------------------------
*     7  ... BIOF in in the downstream river -----------------------------------
*     8  ... BIOF in in the upstream river -------------------------------------
*     9  ... dissolved organic carbon in the discharge -------------------------
*     10 ... not used ----------------------------------------------------------
*     11 ... the bioavailable concentration in the downstream river ------------
*     12 ... the bioavailable concentration in the upstream river --------------
*     13 ... calcium in the river ----------------------------------------------
*     14 ... numbers of samples downstream of modified discharge ---------------
*     15 ... numbers of samples downstream of current discharge ----------------

      rmin = 99999.9
      rmax = -99999.0
      rmen = 0.0
      do i = 1, NS
      rmin = amin1 ( rmin, C(14,i) )    
      rmax = amax1 ( rmax, C(14,i) ) 
      rmen = rmen + C(14,i)
      enddo
      rmen = rmen / float(NS)
  
      cmin = 99999.9
      cmax = -99999.0
      cmen = 0.0
      do i = 1, NS
      cmin = amin1 ( cmin, C(15,i) )    
      cmax = amax1 ( cmax, C(15,i) ) 
      cmen = cmen + C(15,i)
      enddo
      cmen = cmen / float(NS)

      do K = 1, NAR-2
      do I = NS, 1, -1
      do J = 1, I - 1
      if (C(K,I) .lt. C(K,J)) then
      CC = C(K,I)
      C(K,I) = C(K,J)
      C(K,J) = CC
      
      if ( K .eq. 1 ) then
      CC = C(14,I)
      C(14,I) = C(14,J)
      C(14,J) = CC
      endif
      if ( K .eq. 1 ) then
      CC = C(15,I)
      C(15,I) = C(15,J)
      C(15,J) = CC
      endif
      
      endif
      enddo
      enddo 
      enddo
      
      return
      end

*     compute the mean and standard deviation of logged variables -------------- 
      subroutine compute logged summary statistics 
      include 'bacom1cp.for'

*     river flow ---------------------------------------------------------------
      if (rfm .lt. Small) then
      grfm = - Big
      grfs =   0.0
      else
*     standard deviation in log domain -----------------------------------------
      grfs = SQRT( 2.7057 + 2.0 * ALOG((rfm + rft)/(rf5 + rft)) )
     & - 1.6449
*     mean in the log domain ---------------------------------------------------
      grfm = ALOG (rfm + rft) - 0.5 * grfs * grfs
*     compute standard deviation of unlogged flows -----------------------------  
      rfs  = rfm * SQRT ( EXP( grfs * grfs ) - 1.0 )
      endif ! river flow -------------------------------------------------------

*     dissolved metal in the upstream river ------------------------------------
      if ( rcm .lt. Small ) then
      grcm = - Big
      grcs =   0.0
      else
*     mean in log domain -------------------------------------------------------
      grcm = ALOG( (rcm*rcm) / SQRT(rcm*rcm + rcs*rcs) )
*     standard deviation in log domain -----------------------------------------
      grcs = SQRT( ALOG(1.0 + (rcs*rcs)/(rcm*rcm)) )
      endif ! metal in the upstream river --------------------------------------

*     dissolved organic carbon in upstream river -------------------------------
      if ( udocm .lt. Small ) then
      gudocm = - Big
      gudocs =   0.0
      else
*     mean in log domain -------------------------------------------------------
      gudocm = ALOG( (udocm*udocm) / SQRT(udocm*udocm + udocs*udocs) )
*     standard deviation in log domain -----------------------------------------
      gudocs = SQRT( ALOG(1.0 + (udocs*udocs)/(udocm*udocm)) )
      endif ! dissolved organic carbon in upstream river -----------------------
      
*     discharge flow -----------------------------------------------------------
      if ( efm .lt. Small ) then
      gefm = - Big
      gefs = 0.0
      else
*     mean in log domain -------------------------------------------------------
      gefm = ALOG( (efm*efm) / SQRT(efm*efm + efs*efs) )
*     standard deviation in log domain -----------------------------------------
      gefs = SQRT( ALOG(1.0 + (efs*efs)/(efm*efm)) )
      endif ! discharge flow ---------------------------------------------------

*     compute the mean and standard deviation of metal in the discharge --------
      if ( ecm .lt. Small ) then
      gecm = - Big
      gecs =   0.0
      else
      gecm = ALOG( (ecm*ecm)/SQRT(ecm*ecm + ecs*ecs) )
      gecs = SQRT( ALOG(1. + (ecs*ecs)/(ecm*ecm)) )
      endif
      
*     compute the mean and standard deviation of dissolved organic carbon in ---
*     the discharge -------------------------------------------------------------
      if ( edocm .lt. Small ) then
      gedocm = - Big
      gedocs = 0.0
      else
      gedocm = ALOG( (edocm*edocm)/SQRT(edocm*edocm + edocs*edocs) )
      gedocs = SQRT( ALOG(1.0 + (edocs*edocs)/(edocm*edocm)) )
      endif
      
*     compute the mean and standard deviation of pH in the downstream river ----
      if ( phm .lt. Small ) then
      gphm = - Big
      gphs = 0.0
      else
      gphm = ALOG( (phm*phm)/SQRT(phm*phm + phs*phs) )
      gphs = SQRT( ALOG(1.0 + (phs*phs)/(phm*phm)) )
      endif
     
*     compute the mean and standard deviation of calcium in the downstream river
      if ( cam .lt. Small ) then
      gcam = - Big
      gcas = 0.0
      else
      gcam = ALOG( (cam*cam)/SQRT(cam*cam + cas*cas) )
      gcas = SQRT( ALOG(1.0 + (cas*cas)/(cam*cam)) )
      endif

      do is = 1, NS+1
      do k = 1, NAR-2
      C(k,is) = 0.0 ! initialise all the shots
      enddo
      enddo

      return
      end


      
      


      subroutine get the shot for river flow (IS,RF)
      include 'bacom1cp.for'
      
      RF = 0.0
      if ( rfpd .eq. 0 ) then ! contant river flow -----------------------------
      RF = amax1 ( small, rfm )
      return
      endif
      
      if ( rfpd .eq. 1 ) then ! normal distribution for river flow
      RF = amax1 ( small, ( R1 * rfs + rfm) )
      if ( nbias .eq. 1 ) RF = BIASF (RF,BM(1),BS(1),xrf5)
      return
      endif
      
      if ( rfpd .eq. 2 ) then ! 2-parameter log-normal river flow --------------
      RF = EXP (R1 * grfs + grfm)
      if ( nbias .eq. 1 ) RF = BIASF (RF,BM(1),BS(1),xrf5)
      return
      endif

      if ( rfpd .eq. 3 ) then ! 3-parameter log-normal river flow --------------
      RF = amax1 ( small, EXP (R1 * grfs + grfm) - rft )
      if ( nbias .eq. 1 ) RF = BIASF (RF,BM(1),BS(1),xrf5)
      return
      endif

*     non-parametric upstream river flow ---------------------------------------
      if ( rfpd .eq. 4 ) call get non parametric shot ( R1, EF, 1)

      return
      end
      
      subroutine get the shot for river quality (IS,RC)
      include 'bacom1cp.for'
      
      RC = 0.0
*     contant river quality ----------------------------------------------------
      if ( rcpd .eq. 0 ) then
      RC = amax1 ( small, rcm )
      return
      endif
      
*     normal upstream river quality --------------------------------------------
      if ( rcpd .eq. 1 ) then
      RC = amax1 ( 0.001*xrcm, ( R2 * xrcs + xrcm) )
      if ( nbias .eq. 1 ) RC = BIASQ (RC,BM(2),BS(2),xrcm)
      return
      endif
      
*     log-normal upstream river quality ----------------------------------------
      if ( rcpd .eq. 2 ) then
      RC = EXP (R2 * grcs + grcm)
      if ( nbias .eq. 1 ) RC = BIASQ (RC,BM(2),BS(2),xrcm)
      return
      endif

*     log-normal upstream river quality ----------------------------------------
      if ( rcpd .eq. 3 ) then
      RC = amax1 ( small, EXP (R2 * grcs + grcm) - rct )
      if ( nbias .eq. 1 ) RC = BIASQ (RC,BM(2),BS(2),xrcm)
      return
      endif

*     non-parametric upstream river quality ------------------------------------
      if ( rcpd .eq. 4 ) call get non parametric shot ( R2, RC, 2)

      return
      end

      subroutine get the shot for discharge flow (IS,EF)
      include 'bacom1cp.for'

      EF = 0.0
*     contant discharge flow ---------------------------------------------------
      if ( efpd .eq. 0 ) then
      EF = amax1 ( small, efm )
      return
      endif
      
*     normal discharge flow ----------------------------------------------------
      if ( efpd .eq. 1 ) then
      EF = amax1 ( 0.001*xefm, ( R3 * xefs + xefm) )
      if ( nbias .eq. 1 ) EF = BIASQ (EF,BM(3),BS(3),xefm)
      return
      endif
      
*     log-normal discharge flow ------------------------------------------------
      if ( efpd .eq. 2 ) then
      EF = EXP (R3 * gefs + gefm)
      if ( nbias .eq. 1 ) EF = BIASQ (EF,BM(3),BS(3),xefm)
      return
      endif

*     log-normal discharge flow ------------------------------------------------
      if ( efpd .eq. 3 ) then
      EF = amax1 ( small, EXP (R3 * gefs + gefm) - eft )
      if ( nbias .eq. 1 ) EF = BIASQ (EF,BM(3),BS(3),xefm)
      return
      endif

*     non-parametric discharge flow --------------------------------------------
      if ( efpd .eq. 4 ) call get non parametric shot ( R3, EF, 3)

      return
      end
      
      subroutine get the shot for discharge quality (IS,EC)
      include 'bacom1cp.for'
      
      EC = 0.0
*     contant discharge quality ------------------------------------------------
      if ( ecpd .eq. 0 ) then
      EC = amax1 ( small, ecm )
      return
      endif
      
*     normal discharge quality ----------------------------------------
      if ( ecpd .eq. 1 ) then
      EC = amax1 ( 0.001*xecm, ( R4 * xecs + xecm) )
      if ( nbias .eq. 1 ) EC = BIASQ (EC,BM(4),BS(4),xecm)
      return
      endif
      
*     log-normal discharge quality ------------------------------------
      if ( ecpd .eq. 2 ) then
      EC = EXP (R4 * gecs + gecm)
      if ( nbias .eq. 1 ) EC = BIASQ (EC,BM(4),BS(4),xecm)
      endif

*     log-normal discharge quality ------------------------------------
      if ( ecpd .eq. 3 ) then
      EC = amax1 ( small, EXP (R4 * gecs + gecm) - ect )
      if ( nbias .eq. 1 ) EC = BIASQ (EC,BM(4),BS(4),xecm)
      return
      endif

*     non-parametric discharge quality ---------------------------------
      if ( ecpd .eq. 4 ) call get non parametric shot ( R4, EC, 4)
 
      return
      end
      
      subroutine get the shot for upstream DOC (IS,RD)
      include 'bacom1cp.for'
      
      RD = 0.0
*     contant upstream dissolved organic carbon ---------------------------------
      if ( rdpd .eq. 0 ) then
      RD = amax1 ( small, udocm )
      return
      endif
      
*     normal upstream dissolved organic carbon ----------------------------------
      if ( rdpd .eq. 1 ) then
      RD = amax1 ( 0.001*xudocm, ( R6 * xudocs + xudocm) )
      if ( nbias .eq. 1 ) RD = BIASQ (RD,BM(7),BS(7),xudocm)
      return
      endif
      
*     log-normal upstream dissolved organic carbon ------------------------------
      if ( rdpd .eq. 2 ) then
      RD = EXP (R6 * gudocs + gudocm)
      if ( nbias .eq. 1 ) RD = BIASQ (RD,BM(7),BS(7),xudocm)
      return
      endif

*     log-normal upstream dissolved organic carbon ------------------------------
      if ( rdpd .eq. 3 ) then
      RD = amax1 ( small, EXP (R6 * gudocs + gudocm) - rdt )
      if ( nbias .eq. 1 ) RD = BIASQ (RD,BM(7),BS(7),xudocm)
      return
      endif

*     non-parametric upstream dissolved organic carbon --------------------------
      if ( rdpd .eq. 4 ) RD = EXP (R6 * gudocs + gudocm)

      return
      end
      
      subroutine get the shot for discharge DOC (IS,ED)
      include 'bacom1cp.for'
      
      ED = 0.0
*     contant discharge DOC ----------------------------------------------------
      if ( edpd .eq. 0 ) then
      ED = amax1 ( small, edocm )
      return
      endif
      
*     normal discharge DOC -----------------------------------------------------
      if ( edpd .eq. 1 ) then
      ED = amax1 ( 0.001*xedocm, ( R5 * xedocs + xedocm) )
      if ( nbias .eq. 1 ) ED = BIASQ (ED,BM(8),BS(8),xedocm)
      return
      endif
      
*     log-normal discharge DOC -------------------------------------------------
      if ( edpd .eq. 2 ) then
      ED = EXP (R5 * gedocs + gedocm)
      if ( nbias .eq. 1 ) ED = BIASQ (ED,BM(8),BS(8),xedocm)
      return
      endif

*     log-normal discharge DOC -------------------------------------------------
      if ( edpd .eq. 3 ) then
      ED = amax1 ( small, EXP (R5 * gedocs + gedocm) - edt )
      if ( nbias .eq. 1 ) ED = BIASQ (ED,BM(8),BS(8),xedocm)
      return
      endif

*     non-parametric discharge DOC ----------------------------------------------
      if ( edpd .eq. 4 ) then
      ED = EXP (R5 * gedocs + gedocm)
      return
      endif
      return
      end

      subroutine get the shot for downstream pH (IS,PH)
      include 'bacom1cp.for'

      PH = 0.0
*     contant downstream pH ----------------------------------------------------
      if ( phpd .eq. 0 ) then
      PH = AMAX1 (4.5, AMIN1(xphm, 10.0))
      return
      endif
      
*     normal downstream pH -----------------------------------------------------
      if ( phpd .eq. 1 ) then
      PH = AMAX1 (4.5, AMIN1(xphm + R7 * xphs, 10.0))
      if ( nbias .eq. 1 ) PH = BIASQ (PH,BM(5),BS(5),xphm)
      return
      endif
      
*     log-normal downstream pH -------------------------------------------------
      if ( phpd .eq. 2 ) then
      PH = AMAX1 (4.5, AMIN1(EXP (R7 * gphs + gphm), 10.0))
      if ( nbias .eq. 1 ) PH = BIASQ (PH,BM(5),BS(5),xphm)
      return
      endif

*     log-normal downstream pH -------------------------------------------------
      if ( phpd .eq. 3 ) then
      PH = AMAX1 (4.5, AMIN1(EXP (R7 * gphs + gphm) - pht, 10.0))
      if ( nbias .eq. 1 ) PH = BIASQ (PH,BM(5),BS(5),xphm)
      return
      endif
      
*     non-parametric downstream pH ---------------------------------------------
      if ( phpd .eq. 4 ) then
      PH = AMAX1 (4.5, AMIN1(xphm + R7 * xphs, 10.0))
      return
      endif
      return
      end

 
      subroutine get the shot for downstream calcium (IS,CA)
      include 'bacom1cp.for'
      
      CA = 0.0
*     contant downstream calcium -----------------------------------------------
      if ( capd .eq. 0 ) then
      CA = amax1 ( small, cam )
      return
      endif
      
*     normal downstream calcium -------------------------------------------------
      if ( capd .eq. 1 ) then
      CA = amax1 ( 0.001*xcam, ( R8 * xcas + xcam) )
      if ( nbias .eq. 1 ) CA = BIASQ (CA,BM(6),BS(6),xcam)
      return
      endif
      
*     log-normal downstream calcium --------------------------------------------
      if ( capd .eq. 2 ) then
      CA = EXP (R8 * gcas + gcam)
      if ( nbias .eq. 1 ) CA = BIASQ (CA,BM(6),BS(6),xcam)
      return
      endif

*     log-normal downstream calcium --------------------------------------------
      if ( capd .eq. 3 ) then
      CA = amax1 ( small, EXP (R8 * gcas + gcam) - cat )
      if ( nbias .eq. 1 ) CA = BIASQ (CA,BM(6),BS(6),xcam)
      return
      endif

*     non-parametric downstream calcium ----------------------------------------
      if ( capd .eq. 4 ) then
*     CA = EXP (R8 * gcas + gcam)
      return
      endif
      
      return
      end
      
      


      subroutine bioav for zinc (biof,Ca,pH,DOC)
      include 'bacom1cp.for'

      pnec = ((znA * pH + (znB)) * alog (Ca) 
     &     + (znC * pH + znD)) * DOC * DOC 
     &     + ((znE * pH + (znF)) * alog (Ca) 
     &     + (znG * pH + (znH))) * DOC 
     &     + ((znI * pH + znJ) * alog (Ca) 
     &     + (znK * pH + (znL)))

*     calculate the impact of the errors in the bioavailability equation =======
      pnecz = pnec ! store the value of pnec that excludes the errors ==========
      biof = 10.9 / pnec 
      biofz = biof ! store the value of biof that excludes the errors ==========

      RR9 = 0.0 ! initialise the random normal deviate SSSSSSSSSSSSSSSSSSSSSSSSS
      pnecERROR = 0.0 ! initialise the error SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      RR9 = GAS9(IR9) ! random deviate for error in the bio-equation ===========
      if ( RR9 .lt. 0.0 ) then
      pnecERROR = (bioshotm + RR9 * bioshots2) * 0.01 ! ========================
      else
      pnecERROR = (bioshotm + RR9 * bioshots1) * 0.01 ! ========================
      endif
      pnec = amax1 ( 0.0, pnec * (1.0 + pnecERROR) ) ! ========================
      
      if ( SITEBIAS .eq. 1 ) then ! insert FIXED error from the bio-equation FFF
      pnecst2 = gofcl * pnec + pnec ! estimate the change in pnec (site-based) F
      pnec = amax1 ( 0.0, pnecst2 ) ! new value of pnec FFFFFFFFFFFFFFFFFFFFFFFF 
      endif ! if ( SITEBIAS .eq. 1 ) insert FIXED errors in bio-equation FFFFFFF

      pnec = amax1 ( 10.9, pnec ) ! ============================================
      biof = 10.9 / pnec 

      return
      end ! zinc
      
      

      subroutine bioav for copper (pnecb,biof,Ca,pH,DOC,ic,iup)
      include 'bacom1cp.for'
      integer pHband

      pHv = pH
      phb = pH

      if ( pH .lt. 4.0 ) then
      pHband = 1
      pHv = 5.5
      goto 1000
      endif
      
      if ( pH .lt. 5.5 ) then
      pHband = 2
      pHv = 5.5
      goto 1000
      endif
      
      if ( pH .le. 8.5 ) then
      pHband = 3
      pHv = 5.5
      goto 1000
      endif
      
      if ( pH .gt. 8.5 ) then
      pHband = 4
      pHv = 8.5
      endif
 1000 continue 
      
      pnecb = ! PNEC basic
     & CU433*pHb**4.0*DOC**3.0*Ca**3.0
     &+CU432*pHb**4.0*DOC**3.0*Ca**2.0
     &+CU431*pHb**4.0*DOC**3.0*Ca**1.0
     &+CU430*pHb**4.0*DOC**3.0*Ca**0.0
     &+CU423*pHb**4.0*DOC**2.0*Ca**3.0
     &+CU422*pHb**4.0*DOC**2.0*Ca**2.0
     &+CU421*pHb**4.0*DOC**2.0*Ca**1.0
     &+CU420*pHb**4.0*DOC**2.0*Ca**0.0
     &+CU413*pHb**4.0*DOC**1.0*Ca**3.0
     &+CU412*pHb**4.0*DOC**1.0*Ca**2.0
     &+CU411*pHb**4.0*DOC**1.0*Ca**1.0
     &+CU410*pHb**4.0*DOC**1.0*Ca**0.0
     &+CU403*pHb**4.0*DOC**0.0*Ca**3.0
     &+CU402*pHb**4.0*DOC**0.0*Ca**2.0
     &+CU401*pHb**4.0*DOC**0.0*Ca**1.0
     &+CU400*pHb**4.0*DOC**0.0*Ca**0.0
     &+CU333*pHb**3.0*DOC**3.0*Ca**3.0
     &+CU332*pHb**3.0*DOC**3.0*Ca**2.0
     &+CU331*pHb**3.0*DOC**3.0*Ca**1.0
     &+CU330*pHb**3.0*DOC**3.0*Ca**0.0
     &+CU323*pHb**3.0*DOC**2.0*Ca**3.0
     &+CU322*pHb**3.0*DOC**2.0*Ca**2.0
     &+CU321*pHb**3.0*DOC**2.0*Ca**1.0
     &+CU320*pHb**3.0*DOC**2.0*Ca**0.0
     &+CU313*pHb**3.0*DOC**1.0*Ca**3.0
     &+CU312*pHb**3.0*DOC**1.0*Ca**2.0
     &+CU311*pHb**3.0*DOC**1.0*Ca**1.0
     &+CU310*pHb**3.0*DOC**1.0*Ca**0.0
     &+CU303*pHb**3.0*DOC**0.0*Ca**3.0
     &+CU302*pHb**3.0*DOC**0.0*Ca**2.0
     &+CU301*pHb**3.0*DOC**0.0*Ca**1.0
     &+CU300*pHb**3.0*DOC**0.0*Ca**0.0
     &+CU233*pHb**2.0*DOC**3.0*Ca**3.0
     &+CU232*pHb**2.0*DOC**3.0*Ca**2.0
     &+CU231*pHb**2.0*DOC**3.0*Ca**1.0
     &+CU230*pHb**2.0*DOC**3.0*Ca**0.0
     &+CU223*pHb**2.0*DOC**2.0*Ca**3.0
     &+CU222*pHb**2.0*DOC**2.0*Ca**2.0
     &+CU221*pHb**2.0*DOC**2.0*Ca**1.0
     &+CU220*pHb**2.0*DOC**2.0*Ca**0.0
     &+CU213*pHb**2.0*DOC**1.0*Ca**3.0
     &+CU212*pHb**2.0*DOC**1.0*Ca**2.0
     &+CU211*pHb**2.0*DOC**1.0*Ca**1.0
     &+CU210*pHb**2.0*DOC**1.0*Ca**0.0
     &+CU203*pHb**2.0*DOC**0.0*Ca**3.0
     &+CU202*pHb**2.0*DOC**0.0*Ca**2.0
     &+CU201*pHb**2.0*DOC**0.0*Ca**1.0
     &+CU200*pHb**2.0*DOC**0.0*Ca**0.0
     &+CU133*pHb**1.0*DOC**3.0*Ca**3.0
     &+CU132*pHb**1.0*DOC**3.0*Ca**2.0
     &+CU131*pHb**1.0*DOC**3.0*Ca**1.0
     &+CU130*pHb**1.0*DOC**3.0*Ca**0.0
     &+CU123*pHb**1.0*DOC**2.0*Ca**3.0
     &+CU122*pHb**1.0*DOC**2.0*Ca**2.0
     &+CU121*pHb**1.0*DOC**2.0*Ca**1.0
     &+CU120*pHb**1.0*DOC**2.0*Ca**0.0
     &+CU113*pHb**1.0*DOC**1.0*Ca**3.0
     &+CU112*pHb**1.0*DOC**1.0*Ca**2.0
     &+CU111*pHb**1.0*DOC**1.0*Ca**1.0
     &+CU110*pHb**1.0*DOC**1.0*Ca**0.0
     &+CU103*pHb**1.0*DOC**0.0*Ca**3.0
     &+CU102*pHb**1.0*DOC**0.0*Ca**2.0
     &+CU101*pHb**1.0*DOC**0.0*Ca**1.0
     &+CU100*pHb**1.0*DOC**0.0*Ca**0.0
     &+CU033*pHb**0.0*DOC**3.0*Ca**3.0
     &+CU032*pHb**0.0*DOC**3.0*Ca**2.0
     &+CU031*pHb**0.0*DOC**3.0*Ca**1.0
     &+CU030*pHb**0.0*DOC**3.0*Ca**0.0
     &+CU023*pHb**0.0*DOC**2.0*Ca**3.0
     &+CU022*pHb**0.0*DOC**2.0*Ca**2.0
     &+CU021*pHb**0.0*DOC**2.0*Ca**1.0
     &+CU020*pHb**0.0*DOC**2.0*Ca**0.0
     &+CU013*pHb**0.0*DOC**1.0*Ca**3.0
     &+CU012*pHb**0.0*DOC**1.0*Ca**2.0
     &+CU011*pHb**0.0*DOC**1.0*Ca**1.0
     &+CU010*pHb**0.0*DOC**1.0*Ca**0.0
     &+CU003*pHb**0.0*DOC**0.0*Ca**3.0
     &+CU002*pHb**0.0*DOC**0.0*Ca**2.0
     &+CU001*pHb**0.0*DOC**0.0*Ca**1.0
     &+CU000*pHb**0.0*DOC**0.0*Ca**0.0
      
      pnecbt = amax1 ( 1.0, pnecb ) ! PNEC basic and trimmed
      biofb = 1.0 / pnecbt ! BioF basic 
 
      pnecd = ! PNEC default
     & CU433*pHv**4.0*DOC**3.0*Ca**3.0
     &+CU432*pHv**4.0*DOC**3.0*Ca**2.0
     &+CU431*pHv**4.0*DOC**3.0*Ca**1.0
     &+CU430*pHv**4.0*DOC**3.0*Ca**0.0
     &+CU423*pHv**4.0*DOC**2.0*Ca**3.0
     &+CU422*pHv**4.0*DOC**2.0*Ca**2.0
     &+CU421*pHv**4.0*DOC**2.0*Ca**1.0
     &+CU420*pHv**4.0*DOC**2.0*Ca**0.0
     &+CU413*pHv**4.0*DOC**1.0*Ca**3.0
     &+CU412*pHv**4.0*DOC**1.0*Ca**2.0
     &+CU411*pHv**4.0*DOC**1.0*Ca**1.0
     &+CU410*pHv**4.0*DOC**1.0*Ca**0.0
     &+CU403*pHv**4.0*DOC**0.0*Ca**3.0
     &+CU402*pHv**4.0*DOC**0.0*Ca**2.0
     &+CU401*pHv**4.0*DOC**0.0*Ca**1.0
     &+CU400*pHv**4.0*DOC**0.0*Ca**0.0
     &+CU333*pHv**3.0*DOC**3.0*Ca**3.0
     &+CU332*pHv**3.0*DOC**3.0*Ca**2.0
     &+CU331*pHv**3.0*DOC**3.0*Ca**1.0
     &+CU330*pHv**3.0*DOC**3.0*Ca**0.0
     &+CU323*pHv**3.0*DOC**2.0*Ca**3.0
     &+CU322*pHv**3.0*DOC**2.0*Ca**2.0
     &+CU321*pHv**3.0*DOC**2.0*Ca**1.0
     &+CU320*pHv**3.0*DOC**2.0*Ca**0.0
     &+CU313*pHv**3.0*DOC**1.0*Ca**3.0
     &+CU312*pHv**3.0*DOC**1.0*Ca**2.0
     &+CU311*pHv**3.0*DOC**1.0*Ca**1.0
     &+CU310*pHv**3.0*DOC**1.0*Ca**0.0
     &+CU303*pHv**3.0*DOC**0.0*Ca**3.0
     &+CU302*pHv**3.0*DOC**0.0*Ca**2.0
     &+CU301*pHv**3.0*DOC**0.0*Ca**1.0
     &+CU300*pHv**3.0*DOC**0.0*Ca**0.0
     &+CU233*pHv**2.0*DOC**3.0*Ca**3.0
     &+CU232*pHv**2.0*DOC**3.0*Ca**2.0
     &+CU231*pHv**2.0*DOC**3.0*Ca**1.0
     &+CU230*pHv**2.0*DOC**3.0*Ca**0.0
     &+CU223*pHv**2.0*DOC**2.0*Ca**3.0
     &+CU222*pHv**2.0*DOC**2.0*Ca**2.0
     &+CU221*pHv**2.0*DOC**2.0*Ca**1.0
     &+CU220*pHv**2.0*DOC**2.0*Ca**0.0
     &+CU213*pHv**2.0*DOC**1.0*Ca**3.0
     &+CU212*pHv**2.0*DOC**1.0*Ca**2.0
     &+CU211*pHv**2.0*DOC**1.0*Ca**1.0
     &+CU210*pHv**2.0*DOC**1.0*Ca**0.0
     &+CU203*pHv**2.0*DOC**0.0*Ca**3.0
     &+CU202*pHv**2.0*DOC**0.0*Ca**2.0
     &+CU201*pHv**2.0*DOC**0.0*Ca**1.0
     &+CU200*pHv**2.0*DOC**0.0*Ca**0.0
     &+CU133*pHv**1.0*DOC**3.0*Ca**3.0
     &+CU132*pHv**1.0*DOC**3.0*Ca**2.0
     &+CU131*pHv**1.0*DOC**3.0*Ca**1.0
     &+CU130*pHv**1.0*DOC**3.0*Ca**0.0
     &+CU123*pHv**1.0*DOC**2.0*Ca**3.0
     &+CU122*pHv**1.0*DOC**2.0*Ca**2.0
     &+CU121*pHv**1.0*DOC**2.0*Ca**1.0
     &+CU120*pHv**1.0*DOC**2.0*Ca**0.0
     &+CU113*pHv**1.0*DOC**1.0*Ca**3.0
     &+CU112*pHv**1.0*DOC**1.0*Ca**2.0
     &+CU111*pHv**1.0*DOC**1.0*Ca**1.0
     &+CU110*pHv**1.0*DOC**1.0*Ca**0.0
     &+CU103*pHv**1.0*DOC**0.0*Ca**3.0
     &+CU102*pHv**1.0*DOC**0.0*Ca**2.0
     &+CU101*pHv**1.0*DOC**0.0*Ca**1.0
     &+CU100*pHv**1.0*DOC**0.0*Ca**0.0
     &+CU033*pHv**0.0*DOC**3.0*Ca**3.0
     &+CU032*pHv**0.0*DOC**3.0*Ca**2.0
     &+CU031*pHv**0.0*DOC**3.0*Ca**1.0
     &+CU030*pHv**0.0*DOC**3.0*Ca**0.0
     &+CU023*pHv**0.0*DOC**2.0*Ca**3.0
     &+CU022*pHv**0.0*DOC**2.0*Ca**2.0
     &+CU021*pHv**0.0*DOC**2.0*Ca**1.0
     &+CU020*pHv**0.0*DOC**2.0*Ca**0.0
     &+CU013*pHv**0.0*DOC**1.0*Ca**3.0
     &+CU012*pHv**0.0*DOC**1.0*Ca**2.0
     &+CU011*pHv**0.0*DOC**1.0*Ca**1.0
     &+CU010*pHv**0.0*DOC**1.0*Ca**0.0
     &+CU003*pHv**0.0*DOC**0.0*Ca**3.0
     &+CU002*pHv**0.0*DOC**0.0*Ca**2.0
     &+CU001*pHv**0.0*DOC**0.0*Ca**1.0
     &+CU000*pHv**0.0*DOC**0.0*Ca**0.0

      pnecdt = amax1 ( 1.0, pnecd )
      biofd = 1.0 / pnecdt ! Biof Default

*     BIOF extender calculations -----------------------------------------------
      if ( pHband .eq. 2 ) then ! pH is less than 5.5
      vaShape2 = ((1.0 - (1.0 - 2.718282**-(((1.0 - biofd) / 0.027) 
     & **0.35))) * 12.0 + 1.0) - (0.5 * ((1.0 - biofd)**4.0))
      BioFe = (1.0 -(1.0 - 2.718282**-(((pH - 4.0) / vaShape2)** 2.7)))
      PNECe = 1.0 /BioFe ! wowowowowowowowowowowowo
*     check whether the above "1.0" values are "targit" ------------------------

      if ( PNECe .le. 1.0 ) then ! trim
      PNECet = 1.0 ! PNEC extended and trimmed
      else
      PNECet = PNECe ! PNEC extended and trimmed
      endif
      endif ! BIOF extender calculations
    
      if ( pHband .eq. 1 ) then ! pH is less than 4.0 
      pnec = 1.0 ! PNEC generic
      biof = 1.0 ! PNEC generic
      endif
      if ( pHband .eq. 2 ) then ! pH is less than 5.5
      pnec = PNECet ! PNEC extended and trimmed
      biof = BioFe ! PNEC extended
      endif
      if ( pHband .eq. 3 ) then ! pH lies between 5.5 and 8.5
      pnec = pnecbt ! PNEC basic and trimmed
      biof = BioFb ! PNEC basic
      endif
      if ( pHband .eq. 4 ) then ! pH is greater than 8.5
      pnec = pnecdt ! PNEC default and trimmed
      biof = biofd ! PNEC default
      endif
      
*     calculate the impact of the errors in the bioavailability equation =======
      pnecz = pnec ! store the value of pnec that excludes the errors ==========
      biofz = biof ! store the value of biof that excludes the errors ==========

      RR9 = 0.0 ! initialise the random normal deviate SSSSSSSSSSSSSSSSSSSSSSSSS
      pnecERROR = 0.0 ! initialise the error SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      RR9 = GAS9(IR9) ! random deviate for error in the bio-equation SSSSSSSSSSS
      if ( pHband .gt. 1 ) then ! exclude this where the pH is less than 4 =====
      if ( RR9 .lt. 0.0 ) then
      pnecERROR = (bioshotm + RR9 * bioshots2) * 0.01 ! ========================
      else
      pnecERROR = (bioshotm + RR9 * bioshots1) * 0.01 ! ========================
      endif
      pnec = amax1 ( 0.0, pnec * (1.0 + pnecERROR) ) ! =========================   
      endif ! if ( pHband .gt. 1 ) =============================================

      if ( SITEBIAS .eq. 1 ) then ! insert FIXED error from the bio-equation FFF
      pnecst2 = gofcl * pnec + pnec ! estimate the change in pnec (site-based) F
      if ( pHband .gt. 1 ) then ! exclude this where pH is less than 4 FFFFFFFFF
      pnec = amax1 ( 0.0, pnecst2 ) ! new value of pnec FFFFFFFFFFFFFFFFFFFFFFFF 
      endif ! if ( pHband .gt. 1 ) FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
      endif ! if ( SITEBIAS .eq. 1 ) insert FIXED errors in bio-equation FFFFFFF

      if ( inkbias .eq. 1 ) then ! correct the bias in the equation ============
      pnec = (100.0 - bioshotCUm) * 0.01 * amax1 ( 1.0, pnec ) ! =============== 
      else
      pnec = amax1 ( 1.0, pnec ) ! ============================================= 
      endif ! if ( inkbias .eq. 1 ) ============================================
      biof = 1.0 / pnec ! new value of biof with the error =====================

      return
      end
      
      
      
      subroutine bioav for manganese (biof,Ca,pH,DOC,jj)
      include 'bacom1cp.for'

      if ( Ca .lt.   6.0 ) i = 1 ! minimum calcium is 6.0
      if ( Ca .ge.   6.0 .and. Ca .lt.  24.0 ) i = 2     
      if ( Ca .ge.  24.0 .and. Ca .lt.  70.0 ) i = 3     
      if ( Ca .ge.  70.0 .and. Ca .lt.  90.0 ) i = 4     
      if ( Ca .ge.  90.0 .and. Ca .lt. 105.0 ) i = 5     
      if ( Ca .ge. 105.0 .and. Ca .lt. 150.0 ) i = 6     
      if ( Ca .ge. 150.0 .and. Ca .lt. 170.0 ) i = 7     
      if ( Ca .ge. 170.0 .and. Ca .lt. 190.0 ) i = 8     
      if ( Ca .ge. 190.0 .and. Ca .lt. 200.0 ) i = 9 ! maximum is 200.0    

      Cat = Ca
      
      PNEC1 = xa(i) * DOC + xb(i)
      PNEC2 = xc(i) * pH * pH - xd(i) * pH + xe(i) 
      PNEC3 = xf(i) * Ca**xg(i) + xh(i)
      
      EC10 = amax1 ( 0.00001, PNEC1/PNEC2 + PNEC3 )
      BIOF1 = 0.426/EC10
      EC10 = amax1 ( 0.0001, ( 10.0**( (xi(i)*pH) + xj(i) ) ))	
      BIOF2 = 0.310/EC10
      
      BIOF3 = amax1 (BIOF1,BIOF2)
      biof = targit/BIOF3 ! Richard III

      EC10 = (( (xa(i) * DOC) + xb(i) ) / 
     &       (( xc(i)*(pH*pH))-(xd(i)*pH)+xe(i)))
     &       + (xf(i)*(Ca**xg(i)))+xh(i)
      EC10 = amax1 ( 0.0001, EC10 )	

      BIOF1 = 0.426/EC10
      jj = 1

*     algal equation ------------------------------------------------------------
      EC10 = amax1 ( 0.0001, ( 10.0**( (xi(i)*pH) + xj(i) ) ))	
      BIOF2 = 0.310/EC10
      if ( biof2 .gt. biof1 ) jj = 2
      
      biof = amax1(BIOF1,BIOF2)
      pnec = 0.31 / biof

*     calculate the impact of the errors in the bioavailability equation =======
      pnecz = pnec ! store the value of pnec that excludes the errors ==========
      biofz = biof ! store the value of biof that excludes the errors ==========

      RR9 = 0.0 ! initialise the random normal deviate SSSSSSSSSSSSSSSSSSSSSSSSS
      pnecERROR = 0.0 ! initialise the error SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      RR9 = GAS9(IR9) ! random deviate for error in the bio-equation ===========
      if ( RR9 .lt. 0.0 ) then
      pnecERROR = (bioshotm + RR9 * bioshots2) * 0.01 ! ========================
      else
      pnecERROR = (bioshotm + RR9 * bioshots1) * 0.01 ! ========================
      endif
      pnec = amax1 ( 0.0, pnec * (1.0 + pnecERROR) ) ! ======================== 
      
      if ( SITEBIAS .eq. 1 ) then ! insert FIXED error from the bio-equation FFF
      pnecst2 = gofcl * pnec + pnec ! estimate the change in pnec (site-based) F
      pnec = amax1 ( 0.0, pnecst2 ) ! new value of pnec FFFFFFFFFFFFFFFFFFFFFFFF 
      endif ! if ( SITEBIAS .eq. 1 ) insert FIXED errors in bio-equation FFFFFFF

      pnec = amax1 ( 0.31, pnec ) ! ======================== 
      biof = 0.31 / pnec 

      return
      end

      subroutine bioav for nickel (biof,Ca,pH,DOC)
      include 'bacom1cp.for'

      pnec = 
     &   xni300 * (pH**3.0)
     & + xni210 * (pH**2.0) * DOC
     & + xni201 * (pH**2.0) * Ca
     & + xni200 * (pH**2.0)
     & + xni120 * pH * (DOC**2.0)
     & + xni111 * pH * DOC * Ca
     & + xni110 * pH * DOC
     & + xni102 * pH * (Ca**2.0)
     & + xni101 * pH * Ca
     & + xni100 * pH
     & + xni030 * (DOC**3.0)
     & + xni021 * (DOC**2.0) * Ca
     & + xni020 * (DOC**2.0)
     & + xni012 * DOC * (Ca**2.0)
     & + xni011 * DOC * Ca
     & + xni010 * DOC
     & + xni003 * (Ca**3.0)
     & + xni002 * (Ca**2.0)
     & + xni001 * Ca
     & + xni000
      
      biof = 4.0 / pnec
      pnecz = pnec ! store the value of pnec that excludes the errors ==========
      biofz = biof ! store the value of biof that excludes the errors ==========

      RR9 = 0.0 ! initialise the random normal deviate SSSSSSSSSSSSSSSSSSSSSSSSS
      pnecERROR = 0.0 ! initialise the error SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      RR9 = GAS9(IR9) ! random deviate for error in the bio-equation ===========
      if ( RR9 .lt. 0.0 ) then
      pnecERROR = (bioshotm + RR9 * bioshots2) * 0.01 ! ========================
      else
      pnecERROR = (bioshotm + RR9 * bioshots1) * 0.01 ! ========================
      endif
      pnec = amax1 ( 0.0, pnec * (1.0 + pnecERROR) ) ! ========================= 
      
      if ( SITEBIAS .eq. 1 ) then ! insert FIXED error from the bio-equation FFF
      pnecst2 = gofcl * pnec + pnec ! estimate the change in pnec (site-based) F
      pnec = amax1 ( 0.0, pnecst2 ) ! new value of pnec FFFFFFFFFFFFFFFFFFFFFFFF 
      endif ! if ( SITEBIAS .eq. 1 ) insert FIXED errors in bio-equation FFFFFFF

      pnec = amax1 ( 4.0, pnec ) ! ========================= 
      biof = 4.0 / pnec
      
      return
      end

      
      
      subroutine bioav for lead (biof,DOC)
      include 'bacom1cp.for'

      EQS = 1.2
      Cons = 1.2

      if (DOC .lt.  1.0) then ! lower limit 
      pnec = EQS
      else
      pnec = EQS + Cons * (DOC - 1.0)
      endif
       
      if (DOC .gt. 20.0) then ! upper limit
      DOC = 20.0
      pnec = EQS + Cons * (DOC - 1.0)
      endif
      
      Biof = EQS / pnec
      
      return
      end

      
      subroutine get the shot for river flow 2 (IS,RF)
      include 'bacom1cp.for'
      
      RF = 0.0
      if ( rfpd .eq. 0 ) then ! contant river flow -----------------------------
      RF = amax1 ( small, rfm )
      return
      endif
      
      if ( rfpd .eq. 1 ) then ! normal distribution for river flow
      RF = amax1 ( small, ( R1 * rfs + rfm) )
      return
      endif
      
      if ( rfpd .eq. 2 ) then ! 2-parameter log-normal river flow --------------
      RF = EXP (R1 * grfs + grfm)
      return
      endif

      if ( rfpd .eq. 3 ) then ! 3-parameter log-normal river flow --------------
      RF = amax1 ( small, EXP (R1 * grfs + grfm) - rft )
      return
      endif

*     non-parametric upstream river flow ---------------------------------------
      if ( rfpd .eq. 4 ) call get non parametric shot ( R1, EF, 1)

      return
      end
      
      subroutine get the shot for river quality 2 (IS,RC)
      include 'bacom1cp.for'
      
      RC = 0.0
*     contant river quality ----------------------------------------------------
      if ( rcpd .eq. 0 ) then
      RC = amax1 ( small, rcm )
      return
      endif
      
*     normal upstream river quality --------------------------------------------
      if ( rcpd .eq. 1 ) then
      RC = amax1 ( 0.001*xrcm, ( R2 * xrcs + xrcm) )
      return
      endif
      
*     log-normal upstream river quality ----------------------------------------
      if ( rcpd .eq. 2 ) then
      RC = EXP (R2 * grcs + grcm)
      return
      endif

*     log-normal upstream river quality ----------------------------------------
      if ( rcpd .eq. 3 ) then
      RC = amax1 ( small, EXP (R2 * grcs + grcm) - rct )
      return
      endif

*     non-parametric upstream river quality ------------------------------------
      if ( rcpd .eq. 4 ) call get non parametric shot ( R2, RC, 2)

      return
      end

      subroutine get the shot for discharge flow 2 (IS,EF)
      include 'bacom1cp.for'

      EF = 0.0
*     contant discharge flow ---------------------------------------------------
      if ( efpd .eq. 0 ) then
      EF = amax1 ( small, efm )
      return
      endif
      
*     normal discharge flow ----------------------------------------------------
      if ( efpd .eq. 1 ) then
      EF = amax1 ( 0.001*xefm, ( R3 * xefs + xefm) )
      return
      endif
      
*     log-normal discharge flow ------------------------------------------------
      if ( efpd .eq. 2 ) then
      EF = EXP (R3 * gefs + gefm)
      return
      endif

*     log-normal discharge flow ------------------------------------------------
      if ( efpd .eq. 3 ) then
      EF = amax1 ( small, EXP (R3 * gefs + gefm) - eft )
      return
      endif

*     non-parametric discharge flow --------------------------------------------
      if ( efpd .eq. 4 ) call get non parametric shot ( R3, EF, 3)

      return
      end
      
      subroutine get the shot for discharge quality 2 (IS,EC)
      include 'bacom1cp.for'
      
      EC = 0.0
*     contant discharge quality ------------------------------------------------
      if ( ecpd .eq. 0 ) then
      EC = amax1 ( small, ecm )
      return
      endif
      
*     normal upstream discharge quality ----------------------------------------
      if ( ecpd .eq. 1 ) then
      EC = amax1 ( 0.001*xecm, ( R3 * xecs + xecm) )
      return
      endif
      
*     log-normal upstream discharge quality ------------------------------------
      if ( ecpd .eq. 2 ) then
      EC = EXP (R4 * gecs + gecm)
      return
      endif

*     log-normal upstream discharge quality ------------------------------------
      if ( ecpd .eq. 3 ) then
      EC = amax1 ( small, EXP (R4 * gecs + gecm) - ect )
      return
      endif

*     non-parametric upstream discharge quality ---------------------------------
      if ( ecpd .eq. 4 ) call get non parametric shot ( R4, EC, 4)
 
      return
      end
      
      subroutine get the shot for upstream DOC 2 (IS,RD)
      include 'bacom1cp.for'
      
      RD = 0.0
*     contant upstream dissolved organic carbon ---------------------------------
      if ( rdpd .eq. 0 ) then
      RD = amax1 ( small, udocm )
      return
      endif
      
*     normal upstream dissolved organic carbon ----------------------------------
      if ( rdpd .eq. 1 ) then
      RD = amax1 ( 0.001*xudocm, ( R6 * xudocs + xudocm) )
      return
      endif
      
*     log-normal upstream dissolved organic carbon ------------------------------
      if ( rdpd .eq. 2 ) then
      RD = EXP (R6 * gudocs + gudocm)
      return
      endif

*     log-normal upstream dissolved organic carbon ------------------------------
      if ( rdpd .eq. 3 ) then
      RD = amax1 ( small, EXP (R6 * gudocs + gudocm) - rdt )
      return
      endif

*     non-parametric upstream dissolved organic carbon --------------------------
      if ( rdpd .eq. 4 ) RD = EXP (R6 * gudocs + gudocm)

      return
      end
      
      subroutine get the shot for discharge DOC 2 (IS,ED)
      include 'bacom1cp.for'
      
      ED = 0.0
*     contant discharge DOC ----------------------------------------------------
      if ( edpd .eq. 0 ) then
      ED = amax1 ( small, edocm )
      return
      endif
      
*     normal discharge DOC -----------------------------------------------------
      if ( edpd .eq. 1 ) then
      ED = amax1 ( 0.001*xedocm, ( R5 * xedocs + xedocm) )
      return
      endif
      
*     log-normal discharge DOC -------------------------------------------------
      if ( edpd .eq. 2 ) then
      ED = EXP (R5 * gedocs + gedocm)
      return
      endif

*     log-normal discharge DOC -------------------------------------------------
      if ( edpd .eq. 3 ) then
      ED = amax1 ( small, EXP (R5 * gedocs + gedocm) - edt )
      return
      endif

*     non-parametric discharge DOC ----------------------------------------------
      if ( edpd .eq. 4 ) then
      ED = EXP (R5 * gedocs + gedocm)
      return
      endif
      return
      end

      subroutine get the shot for downstream pH 2 (IS,PH)
      include 'bacom1cp.for'

      PH = 0.0
*     contant downstream pH ----------------------------------------------------
      if ( phpd .eq. 0 ) then
      PH = AMAX1 (4.5, AMIN1(xphm, 10.0))
      return
      endif
      
*     normal downstream pH -----------------------------------------------------
      if ( phpd .eq. 1 ) then
      PH = AMAX1 (4.5, AMIN1(xphm + R7 * xphs, 10.0))
      return
      endif
      
*     log-normal downstream pH -------------------------------------------------
      if ( phpd .eq. 2 ) then
      PH = AMAX1 (4.5, AMIN1(EXP (R7 * gphs + gphm), 10.0))
      return
      endif

*     log-normal downstream pH -------------------------------------------------
      if ( phpd .eq. 3 ) then
      PH = AMAX1 (4.5, AMIN1(EXP (R7 * gphs + gphm) - pht, 10.0))
      return
      endif
      
*     non-parametric downstream pH ---------------------------------------------
      if ( phpd .eq. 4 ) then
      PH = AMAX1 (4.5, AMIN1(xphm + R7 * xphs, 10.0))
      return
      endif
      return
      end
      

 
      subroutine get the shot for downstream calcium  2 (IS,CA)
      include 'bacom1cp.for'
      
      CA = 0.0
*     contant downstream calcium -----------------------------------------------
      if ( capd .eq. 0 ) then
      CA = amax1 ( small, cam )
      return
      endif
      
*     normal downstream calcium -------------------------------------------------
      if ( capd .eq. 1 ) then
      CA = amax1 ( 0.001*xcam, ( R8 * xcas + xcam) )
      return
      endif
      
*     log-normal downstream calcium --------------------------------------------
      if ( capd .eq. 2 ) then
      CA = EXP (R8 * gcas + gcam)
      return
      endif

*     log-normal downstream calcium --------------------------------------------
      if ( capd .eq. 3 ) then
      CA = amax1 ( small, EXP (R8 * gcas + gcam) - cat )
      return
      endif

*     non-parametric downstream calcium ----------------------------------------
      if ( capd .eq. 4 ) then
*     CA = EXP (R8 * gcas + gcam)
      return
      endif
      
      return
      end
      
