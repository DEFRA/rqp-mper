*     set up the data for the next iteration ...
      subroutine nextf (ichek)
      include 'bacom1cp.for'

      ichek = 0

      if (ITER .eq. 1) then
*     store the results from the current discharge
      tem2 = tecm
      if ( type .ne. 1 ) then
      out2 = tcxx
      if ( free .eq. 1 ) out2 = tbxx
      else
      out2 = tcm
      if ( free .eq. 1 ) out2 = tbm
      endif      
*     set up the data for the first trial
      if ( out2 .lt. targ ) tecm = tecm * 1.2
      if ( out2 .gt. targ ) tecm = tecm * 0.8
      goto 892
      endif ! if (ITER .eq. 1) then

      if (ITER .eq. 2) then
*     store the results of the first trial
      tem1 = tecm
      if ( type .ne. 1 ) then
      out1 = tcxx
      if ( free .eq. 1 ) out1 = tbxx
      else
      out1 = tcm
      if ( free .eq. 1 ) out1 = tbm
      endif 
      goto 891
      endif

*     iterations after the first two ...
*     check which of the two retained trials is best ...
      !write(3,3000)Iter,out1,tem1,out2,tem2
 3000 format(68('=')/
     &'Iteration: ',i4,' ... out1 and tem1 =',2f10.3/
     &'           ',4x,' ... out2 and tem2 =',2f10.3/68('-'))

      if ( ABS ( out1 - targ ) .lt. ABS ( out2 - targ ) ) then
*     over-write second set
      if ( type .ne. 1 ) then
      out2 = tcxx
      if ( free .eq. 1 ) out2 = tbxx
      else
      out2 = tcm
      if ( free .eq. 1 ) out2 = tbm
      endif 
      tem2 = tecm
      else
*     over-write the first set
      if ( type .ne. 1 ) then
      out1 = tcxx
      if ( free .eq. 1 ) out1 = tbxx
      else
      out1 = tcm
      if ( free .eq. 1 ) out1 = tbm
      endif
      tem1 = tecm
      endif

  891 continue
*     interpolate between latest trial and retained trial
      denom = out2 - out1
      !write(03,3100)Iter,out1,tem1,out2,tem2,
    ! &(targ-out1),(tem2-tem1),
    ! &(targ-out1)*(tem2-tem1),denom
 3100 format(68('-')/
     &'Updated:   ',i4,' ... out1 and tem1 =',2f10.3/
     &'           ',4x,' ... out2 and tem2 =',2f10.3/
     &'           ',4x,' ...   targ - out1 =',f10.3/
     &'           ',4x,' ...   tem2 - tem1 =',f10.3/
     &'           ',4x,' ...          mult =',f10.3/
     &'           ',4x,' ...         denom =',f10.3/68('-'))

      if (denom. gt. Small .or. denom .lt. -Small) then
      tecm = tem1 + (targ-out1) * (tem2-tem1) / denom
      tecs = ecv * tecm
      !write(7,681)tecm
  681 format(15x,' ...          tecm =',f10.3/68('-'))
      
      else
	ichek = 1
      if ( tecm .lt. small ) then
          tecm = 0.0
          tecs = 0.0
      endif
	return
      endif
  892 continue

*     check for a zero discharge quality
      if (tecm .gt. 0.0) return
*     is this the second time that a zero has been obtained?
*     if so, a fault has been detected ...
      if (KDIL .gt. 0) then
      msg1 = 1
*     set up calculation to check whether target can be achieved at all
      else
	IDIL = 1
	KDIL = 1
	tecm = 0.0
	!tecs = 0.0
	return
      endif    

      return
      end



*     sensitivity tests on bio-available items ... the backwards calculations ....
      subroutine sensa1
      include 'bacom1cp.for'

      effmas = tecm 
      if ( tecm .gt. 0.0001 ) then

      xold = phm
      phm = phm / 0.98
      xphm = phm
      call do the calculations (1)
	xa1 = 100.0 * tecm / effmas - 100.0
      phm = xold
      xphm = phm

      phm = 0.98 * phm          
      xphm = phm
      call do the calculations (1)
      xa15 = 100.0 * tecm / effmas - 100.0
      phm = xold
      xphm = phm

      xold = phs
      phs = phs / 0.9
      xphs = phs
      call do the calculations (1)
      xa2 = 100.0 * tecm / effmas - 100.0
      phs = xold
      xphs = phs

      xold = cam
      cam = 0.9 * cam
      xcam = cam
      call do the calculations (1)
      xa3 = 100.0 * tecm / effmas - 100.0
      cam = xold
      xcam = cam
      
      xold = cams
      cas = 0.9 * cas
      xcas = cas
      call do the calculations (1)
      xa4 = 100.0 * tecm / effmas - 100.0
      cas= xold
      xcas = cas

      xold = udocm
      udocm = udocm / 0.9
      xudocm = udocm
      call do the calculations (1)
      xa5 = 100.0 * tecm / effmas - 100.0
      udocm = xold
      xudocm = udocm

      xold = udocs
      udocs = udocs / 0.9
      xudocs = udocs
      call do the calculations (1)
      xa6 = 100.0 * tecm / effmas - 100.0
      udocs = xold
      xudocs = udocs

      xold = edocm
      edocm = edocm / 0.9
      edocmx = edocm
      call do the calculations (1)
      xa7 = 100.0 * tecm / effmas - 100.0
      edocm = xold
      xedocm = edocm

      xold = edocs
      edocs = edocs / 0.9
      xedocs = edocs
      call do the calculations (1)
      xa8 = 100.0 * tecm / effmas - 100.0
      edocs = xold
      xedocs = edocs

      xold = capf
      capf = capf - sign(0.1, capf)
      call do the calculations (1)
      xa9 = 100.0 * tecm / effmas - 100.0
      capf = xold
      tecm = effmas

      xold = cacf
      cacf = cacf - sign(0.1, cacf)
      call do the calculations (1)
      xa10 = 100.0 * tecm / effmas - 100.0
      cacf = xold

      xold = cacp
      cacp = cacp - sign(0.1, cacp)
      call do the calculations (1)
      xa11 = 100.0 * tecm / effmas - 100.0
      cacp = xold

      xold = cafd
      cafd = cafd - sign(0.1, cafd)
      call do the calculations (1)
      xa12 = 100.0 * tecm / effmas - 100.0
      cafd = xold

      xold = caed
      caed = caed - sign(0.1, caed)
      call do the calculations (1)
      xa13 = 100.0 * tecm / effmas - 100.0
      caed = xold
      
      if ( mettal .eq. 1 ) bioshots1 = -bioshotZNs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 1 ) bioshots2 =  bioshotZNs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 2 ) bioshots1 = -bioshotCUs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 2 ) bioshots2 =  bioshotCUs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 3 ) bioshots1 = -bioshotMNs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 3 ) bioshots2 =  bioshotMNs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 4 ) bioshots1 = -bioshotNIs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 4 ) bioshots2 =  bioshotNIs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 5 ) bioshots1 = -bioshotPBs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 5 ) bioshots2 =  bioshotPBs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      call do the calculations (1)
      xa14 = 100.0 * tecm / effmas - 100.0 ! errors in the bio-equation 
      bioshots1 = 0.0
      bioshots2 = 0.0

      if ( mettal .eq. 1 ) bioshotm = bioshotZNm ! bias in the bio-equation SSSS
      if ( mettal .eq. 2 ) bioshotm = bioshotCUm ! bias in the bio-equation SSSS
      if ( mettal .eq. 3 ) bioshotm = bioshotMNm ! bias in the bio-equation SSSS
      if ( mettal .eq. 4 ) bioshotm = bioshotNIm ! bias in the bio-equation SSSS
      if ( mettal .eq. 5 ) bioshotm = bioshotPBm ! bias in the bio-equation SSSS
      call do the calculations (1)
      xa16 = 100.0 * tecm / effmas - 100.0 ! bias in the bio-equation
      bioshotm = 0.0

      if ( mettal .eq. 1 ) bioshotm = bioshotZNm ! bias in the bio-equation SSSS
      if ( mettal .eq. 2 ) bioshotm = bioshotCUm ! bias in the bio-equation SSSS
      if ( mettal .eq. 3 ) bioshotm = bioshotMNm ! bias in the bio-equation SSSS
      if ( mettal .eq. 4 ) bioshotm = bioshotNIm ! bias in the bio-equation SSSS
      if ( mettal .eq. 5 ) bioshotm = bioshotPBm ! bias in the bio-equation SSSS
      if ( mettal .eq. 1 ) bioshots1 = -bioshotZNs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 1 ) bioshots2 =  bioshotZNs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 2 ) bioshots1 = -bioshotCUs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 2 ) bioshots2 =  bioshotCUs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 3 ) bioshots1 = -bioshotMNs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 3 ) bioshots2 =  bioshotMNs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 4 ) bioshots1 = -bioshotNIs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 4 ) bioshots2 =  bioshotNIs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 5 ) bioshots1 = -bioshotPBs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 5 ) bioshots2 =  bioshotPBs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      call do the calculations (1)
      xa17 = 100.0 * tecm / effmas - 100.0 ! bias and errors in the bio-equation
      bioshotm = 0.0
      bioshots1 = 0.0
      bioshots2 = 0.0

      xa1 = amax1 (abs(xa1),abs(xa15))

      write(3,7742)xa1,xa2,xa3,xa4,xa5,xa6,xa7,xa8,xa9,xa10,xa11,xa12,
     &xa13,xa14,xa16,xa17
 7742 format(
     &'Mean pH (2%)                              ',f10.2,' %'/
     &'Standard deviation for pH:                ',f10.2,' %'/
     &'Mean calcium:                             ',f10.2,' %'/
     &'Standard deviation:                       ',f10.2,' %'/
     &'Mean u/s dissolved organic carbon:        ',f10.2,' %'/
     &'Standard deviation:                       ',f10.2,' %'/
     &'Mean discharge dissolved organic carbon:  ',f10.2,' %'/
     &'Standard deviation:                       ',f10.2,' %'/
     &68('-')/
     &'Change of 0.1 in correlation coefficient for:'/68('-')/
     &'river pH and river flow:                  ',f10.2,' %'/
     &'river calcium and river flow:             ',f10.2,' %'/
     &'river calcium and river pH:               ',f10.2,' %'/
     &'river DOC on river flow:                  ',f10.2,' %'/
     &'discharge DOC on discharge flow:          ',f10.2,' %'/
     &68('-')/
     &'errors in the equation:                   ',f10.2,' %'/
     &'bias in the equation:                     ',f10.2,' %'/
     &'bias and errors in the equation:          ',f10.2,' %'/
     &68('='))

      if ( abs(xa15) .gt. abs(xa1) ) xa1 = xa15 
      return

	else
      sensa = 999
      write(3,4430)
 4430 format(//
     &'The sensitivity analysis was suppressed ...'/
     &'Effluent quality is zero ...'//)
      return
      endif

      return
      end





*     sensitivity tests on bio-metal items ... the forwards calculations ....
      subroutine sensa2
      include 'bacom1cp.for'

      loop = 0

      if ( free .eq. 0 ) then
      dowmas = C(1, NS+1)
      else
      dowmas = C(11, NS+1)
      endif

      if ( dowmas .gt. 0.0001 ) then

      phm = phm / 0.975          
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa1 = 100.0 * dqual / dowmas - 100.0
      phm = xphm

      phm = 0.975 * phm          
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa15 = 100.0 * dqual / dowmas - 100.0
      phm = xphm

      phs = phs / 0.9
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa2 = 100.0 * dqual / dowmas - 100.0
      phs = xphs

      cam = 0.9 * cam
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa3 = 100.0 * dqual / dowmas - 100.0
      cam = xcam

      cas = 0.9 * cas
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa4 = 100.0 * dqual / dowmas - 100.0
      cas = xcas

      udocm = udocm / 0.9
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa5 = 100.0 * dqual / dowmas - 100.0
      udocm = xudocm

      udocs = udocs / 0.9
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa6 = 100.0 * dqual / dowmas - 100.0
      udocs = xudocs

      edocm = edocm / 0.9
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa7 = 100.0 * dqual / dowmas - 100.0
      edocm = xedocm

      edocs = edocs / 0.9
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa8 = 100.0 * dqual / dowmas - 100.0
      edocs = xedocs

      capf = capf - sign(0.1, capf)
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa9 = 100.0 * dqual / dowmas - 100.0
      capf = xcapf

      cacf = cacf - sign(0.1, cacf)
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa10 = 100.0 * dqual / dowmas - 100.0
      cacf = xcacf

      cacp= cacp - sign(0.1, cacp)
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa11 = 100.0 * dqual / dowmas - 100.0
      cacp = xcacp

      cafd = cafd - sign(0.1, cafd)
      cafd = cafd / 0.9
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa12 = 100.0 * dqual / dowmas - 100.0
      cafd = xcafd
 
      caed = caed - sign(0.1, caed)
      caed = caed / 0.9
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa13 = 100.0 * dqual / dowmas - 100.0
      caed = xcaed
      
      if ( mettal .eq. 1 ) bioshots1 = -bioshotZNs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 1 ) bioshots2 =  bioshotZNs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 2 ) bioshots1 = -bioshotCUs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 2 ) bioshots2 =  bioshotCUs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 3 ) bioshots1 = -bioshotMNs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 3 ) bioshots2 =  bioshotMNs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 4 ) bioshots1 = -bioshotNIs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 4 ) bioshots2 =  bioshotNIs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 5 ) bioshots1 = -bioshotPBs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 5 ) bioshots2 =  bioshotPBs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa14 = 100.0 * dqual / dowmas - 100.0 ! errors in the bio-equation 
      
      bioshots1 = 0.0
      bioshots2 = 0.0

      if ( mettal .eq. 1 ) bioshotm = bioshotZNm ! bias in the bio-equation SSSS
      if ( mettal .eq. 2 ) bioshotm = bioshotCUm ! bias in the bio-equation SSSS
      if ( mettal .eq. 3 ) bioshotm = bioshotMNm ! bias in the bio-equation SSSS
      if ( mettal .eq. 4 ) bioshotm = bioshotNIm ! bias in the bio-equation SSSS
      if ( mettal .eq. 5 ) bioshotm = bioshotPBm ! bias in the bio-equation SSSS
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa16 = 100.0 * dqual / dowmas - 100.0 ! bias in the bio-equation 

      bioshotm = 0.0

      if ( mettal .eq. 1 ) bioshotm = bioshotZNm ! bias in the bio-equation SSSS
      if ( mettal .eq. 2 ) bioshotm = bioshotCUm ! bias in the bio-equation SSSS
      if ( mettal .eq. 3 ) bioshotm = bioshotMNm ! bias in the bio-equation SSSS
      if ( mettal .eq. 4 ) bioshotm = bioshotNIm ! bias in the bio-equation SSSS
      if ( mettal .eq. 5 ) bioshotm = bioshotPBm ! bias in the bio-equation SSSS
      if ( mettal .eq. 1 ) bioshots1 = -bioshotZNs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 1 ) bioshots2 =  bioshotZNs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 2 ) bioshots1 = -bioshotCUs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 2 ) bioshots2 =  bioshotCUs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 3 ) bioshots1 = -bioshotMNs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 3 ) bioshots2 =  bioshotMNs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 4 ) bioshots1 = -bioshotNIs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 4 ) bioshots2 =  bioshotNIs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 5 ) bioshots1 = -bioshotPBs1 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      if ( mettal .eq. 5 ) bioshots2 =  bioshotPBs2 ! SSSSSSSSSSSSSSSSSSSSSSS 45
      call do the calculations (1)
      if ( free .eq. 0 ) then
           dqual = C(1, NS+1)
      else
           dqual = C(11, NS+1)
      endif
      xa17 = 100.0 * dqual / dowmas - 100.0 ! bias and errors in the bio-equation 

      bioshotm = 0.0
      bioshots1 = 0.0
      bioshots2 = 0.0

      if ( abs(xa15) .gt. abs(xa1) ) xa1 = xa15 

      write(3,7742)xa1,xa2,xa3,xa4,xa5,xa6,xa7,xa8,xa9,xa10,xa11,xa12,
     &xa13,xa14,xa16,xa17
 7742 format('Impact on the mean river quality of a 10%',
     &' change in ...'/68('-')/
     &'Mean pH (2.5% in this case)                   ',f10.2,' %'/
     &'Standard deviation:                           ',f10.2,' %'/
     &'Mean calcium:                                 ',f10.2,' %'/
     &'Standard deviation:                           ',f10.2,' %'/
     &'Mean u/s dissolved organic carbon:            ',f10.2,' %'/
     &'Standard deviation:                           ',f10.2,' %'/
     &'Mean discharge dissolved organic carbon:      ',f10.2,' %'/
     &'Standard deviation:                           ',f10.2,' %'/
     &68('-')/
     &'Change of 0.1 in correlation coefficient for:'/68('-')/
     &'river pH and river flow:                      ',f10.2,' %'/
     &'river calcium and river flow:                 ',f10.2,' %'/
     &'river calcium and river pH:                   ',f10.2,' %'/
     &'river DOC on river flow:                      ',f10.2,' %'/
     &'discharge DOC on discharge flow:              ',f10.2,' %'/
     &68('-')/
     &'errors in the equation:                       ',f10.2,' %'/
     &'bias in the equation:                         ',f10.2,' %'/
     &'bias and errors in the equation:              ',f10.2,' %'/
     &68('='))

      else
      sensa = 999
      endif

      return
      end

      
      


      subroutine add errors from pH etc
      include 'bacom1cp.for'
      
      tecmaster = tecm ! discharge quality =====================================
      tecsplus = tecs
      osterror = tecsplus / SQRoot( 117333, ecns ) ! original standard error ===
      
      tbmaster = tbm ! downstream bioavailble ==================================
      tbmcmaster = tbmc
      tbscmaster = tbsc
      tbsplus = 1.0
      tbnsc = C(15,NS+1)
      bsterror = tbscmaster / SQRoot( 117333, tbnsc) ! original standard error =
 
      ubmaster = ubm ! upstream bioavailble ====================================
      ubsplus = 1.0
      vbsplus = 1.0

      !write(03,8419)tecm,tecs,tecmaster
 8419 format(/50('+')/
     &'     tecm =',f12.4/ 
     &'     tecs =',f12.4/ 
     &'tecmaster =',f12.4/
     &50('+')/)

*     --------------------------------------------------------------------------
      write(6,186)
  186 format(//95('=')/'Calculate the effects of pH sampling:'/
     &95('-'))
      vphmxx = phm ! store the input mean pH
      vphsxx = phs
      phmlx = phml
      phmux = phmu 
      !write(3,97)phm,phs,phml,phmu
   97 format(2f10.3,5x,'(',f7.4,' to',f9.4,')'/95('-'))
      
      phm = phml ! lower confidence limit
      phs = phm * vphsxx/vphmxx
      xphmhold = xphm
      xphshold = xphs
      xphm = phm
      xphs = phs
      call do the calculations (1) ! lower confidence limit
      tecmphml = tecm ! lower confidence limit on discharge quality
      tbmphml = tbmc
      ubmphml = ubm
      
      phm = phmux ! upper confidence limit of pH
      phs = phm * vphsxx/vphmxx
      xphm = phm
      xphs = phs
      call do the calculations (1) ! upper confidence limit
      tecmphmu = tecm ! upper confidence limit on discharge quality
      tbmphmu = tbmc
      ubmphmu = ubm
      
      xphm = xphmhold ! restore original mean pH
      xphs = xphshold

      write(6,1)tecmphml,phmlx,tecmphmu,phmux
    1 format(
     &'The lower confidence limit produces a discharge quality of: ',
     &f10.3,3x,'from a pH of',f10.3/
     &'The upper confidence limit produces one of:',f27.3,3x,
     &'from a pH of',f10.3/95('-'))
      
      t = errx (0.95)
      slph = (tecmphml - tecmaster)/t ! standard error
      suph = (tecmphmu - tecmaster)/t ! standard error
      seph2 = 0.5 * (abs(slph)+abs(suph)) ! average standard error ------------- 
      
      !write(03,8479)tecmaster,tecmphml,tecmphmu,slph,
      !&suph,seph2
 8479 format(/50('+')/
     &'tecmaster =',f12.4/
     &' tecmphml =',f12.4/
     &' tecmphmu =',f12.4/
     &'     slph =',f12.4/
     &'     suph =',f12.4/
     &'    seph2 =',f12.4/
     &50('+')/)
      
 
      slph = (tbmphml - tbmcmaster)/t ! standard error
      suph = (tbmphmu - tbmcmaster)/t ! standard error
      sbph2 = 0.5 * (abs(slph)+abs(suph)) ! average standard error ------------- 

      ulph = (ubmphml - ubmaster)/t ! standard error
      uuph = (ubmphmu - ubmaster)/t ! standard error
      suph2 = 0.5 * (abs(ulph)+abs(uuph)) ! average standard error ------------- 

      phm = vphmxx
      phs = vphsxx

*     **************************************************************************
      write(6,3)osterror,tecsplus
    3 format('The original standard error of discharge quality was:',
     &f9.3,4x,'with sdev of:',f9.3/95('-'))
      suma = osterror*osterror + seph2*seph2 ! sum squares of standard errors --
      suma = SQRoot( 117333,suma) ! compounded standard error ------------------
      ubsa = suma * SQRoot( 117333, ecns ) ! compounded standard deviation -----
      write(6,14)seph2,ubsa ! write to file .EFF
   14 format(/'    The new standard error is: ',f11.4/
     &        'The new standard deviation is: ',f11.4/95('='))
      suma = bsterror*bsterror + sbph2*sbph2 ! sum squares of standard errors --
      suma = SQRoot( 117333,suma) ! compounded standard error ------------------
      bbsa = suma * SQRoot( 117333, rcns ) ! compounded standard deviation -----
      write(6,14)sbph2,bbsa ! write to file .EFF
   64 format(/'    The new standard error is: ',f11.4,' (sbph2)'95('='))
      !write(6,11)
   11 format(/95('=')/'RESTORE THE ORIGINAL DATA - ',
     &'prior to proceeding with calcium'/95('='))
*     **************************************************************************

*     --------------------------------------------------------------------------
      write(6,187)
  187 format(//95('=')/'Calculate the effects of calcium sampling:'/
     &95('-'))
      vcamxx = cam
      vcasxx = cas
      camlx = caml 
      camux = camu
      write(6,98)cam,cas,caml,camu
   98 format(2f10.3,5x,'(',f10.3,' to',f10.3,')'/95('-'))
      cam = caml ! lower confidence limit
      cas = cam * vcasxx/vcamxx
      
      xcamhold = xcam
      xcashold = xcas
      xcam = cam
      xcas = cas

      call do the calculations (1) ! lower confidence limit
      tecmcaml = tecm ! lower confidence limit
      tbmcaml = tbmc
      ubmcaml = ubm
      
      cam = camux ! upper confidence limit for calcium
      cas = cam * vcasxx/vcamxx   
      
      xcam = cam
      xcas = cas
      
      call do the calculations (1) ! upper confidence limit for calcium
      tecmcamu = tecm ! upper confidence limit
      tbmcamu = tbmc
      ubmcamu = ubm
      
      xcam = xcamhold
      xcas = xcashold
      
      write(6,1)tecmcaml,camlx,tecmcamu,camux
      
      t=errx (0.95)
      slca = (tecmcaml - tecmaster)/t ! standard error
      suca = (tecmcamu - tecmaster)/t ! standard error
      seca2 = 0.5 * (abs(slca)+abs(suca)) ! average standard error --------------     

      slca = (tbmcaml - tbmcmaster)/t ! standard error
      suca = (tbmcamu - tbmcmaster)/t ! standard error
      sbca2 = 0.5 * (abs(slca)+abs(suca)) ! average standard error ------------- 

      ulca = (ubmcaml - ubmaster)/t ! standard error
      uuca = (ubmcamu - ubmaster)/t ! standard error
      suca2 = 0.5 * (abs(ulca)+abs(uuca)) ! average standard error ------------- 

      cam = vcamxx
      cas = vcasxx
*     --------------------------------------------------------------------------
 

*     --------------------------------------------------------------------------
      write(6,188)
  188 format(//95('=')/'Calculate the effects of u/s DOC sampling:'/
     &95('-'))
      vudocmxx = udocm
      vudocsxx = udocs
      udocmlx = udocml 
      udocmux = udocmu 

      write(6,198)udocm,udocs,udocml,udocmu
  198 format(2f10.3,5x,'(',f10.3,' to',f10.3,')'/95('-'))

      udocm = udocml ! lower confidence limit
      udocs = udocm * vudocsxx/vudocmxx
      
      xudocmhold = xudocm
      xudocshold = xudocs
      xudocm = udocm
      xudocs = udocs

      call do the calculations (1) ! lower confidence limit
      tecmudocml = tecm ! lower confidence limit
      tbmudocml = tbmc ! lower confidence limit
      ubmudocml = ubm
      
      udocm = udocmux ! upper confidence limit
      udocs = udocm * vudocsxx/vudocmxx
      xudocm = udocm
      xudocs = udocs

      call do the calculations (1) ! upper confidence limit
      tecmudocmu = tecm ! upper confidence limit
      tbmudocmu = tbmc ! lower confidence limit
      ubmudocmu = ubm
      
      xudocm = xudocmhold
      xudocs = xudocshold
      
      write(6,1)tecmudocml,udocmlx,tecmudocmu,udocmux
      
      t=errx (0.95)
      slud = (tecmudocml - tecmaster)/t ! standard error
      suud = (tecmudocmu - tecmaster)/t ! standard error
      sudoc2 = 0.5 * (abs(slud)+abs(suud)) ! biggest standard error --------------     

      slud = (tbmudocml - tbmaster)/t ! standard error
      suud = (tbmudocmu - tbmaster)/t ! standard error
      sbudoc2 = 0.5 * (abs(slud)+abs(suud)) ! average standard error -------------    
     
      slud = (tbmudocml - tbmcmaster)/t ! standard error
      suud = (tbmudocmu - tbmcmaster)/t ! standard error
      sbud2 = 0.5 * (abs(slud)+abs(suud)) ! average standard error ------------- 

      ulud = (ubmudocml - ubmaster)/t ! standard error
      uuud = (ubmudocmu - ubmaster)/t ! standard error
      suud2 = 0.5 * (abs(ulud)+abs(uuud)) ! average standard error ------------- 

      udocm = vudocmxx
      udocs = vudocsxx

      !write(6,32)
   32 format(/95('=')/'RESTORE THE ORIGINAL DATA - ',
     &'proceed with discharge dissolved organic carbon ...')
*     --------------------------------------------------------------------------

      
*     --------------------------------------------------------------------------
      write(6,488)
  488 format(//95('=')/'Calculate the effects of discharge DOC ',
     &'sampling:'/95('-'))
      
      vedocmxx = edocm
      vedocsxx = edocs
      edocmlx = edocml 
      edocmux = edocmu 
      
      write(6,298)edocm,edocs,edocml,edocmu
  298 format(2f10.3,5x,'(',f10.3,' to',f10.3,')'/95('-'))

      edocm = edocml ! lower confidence limit
      edocs = edocm * vudocsxx/vudocmxx
      
      xedocmhold = xedocm
      xedocshold = xedocs
      xedocm = edocm
      xedocs = edocs

      call do the calculations (1) ! lower confidence limit
      tecmedocml = tecm ! lower confidence limit
      tbmedocml = tbmc ! lower confidence limit
      ubmedocml = ubm
      
      edocm = edocmux ! upper confidence limit
      edocs = edocm * vedocsxx/vedocmxx
      xedocm = edocm
      xedocs = edocs

      call do the calculations (1) ! upper confidence limit
      tecmedocmu = tecm ! upper confidence limit
      tbmedocmu = tbmc ! upper confidence limit
      ubmedocmu = ubm
     
      xedocm = xedocmhold
      xedocs = xedocshold
      
      write(6,1)tecmedocml,edocmlx,tecmedocmu,edocmux
      
      t=errx (0.95)
      sled = (tecmedocml - tecmaster)/t ! standard error
      sued = (tecmedocmu - tecmaster)/t ! standard error
      sedoc2 = 0.5 * (abs(sled)+abs(sued)) ! biggest standard error --------------     
      
      sled = (tbmedocml - tbmaster)/t ! standard error
      sued = (tbmedocmu - tbmaster)/t ! standard error
      sbedoc2 = 0.5 * (abs(sled)+abs(sued)) ! average standard error -------------    

      sled = (tbmedocml - tbmcmaster)/t ! standard error
      sued = (tbmedocmu - tbmcmaster)/t ! standard error
      sbed2 = 0.5 * (abs(sled)+abs(sued)) ! average standard error ------------- 

      edocm = vedocmxx
      edocs = vedocsxx
*     --------------------------------------------------------------------------

      write(6,861)
  861 format(//95('=')/'USE THE ADDED ERRORS FOR ALL FOUR ITEMS',
     &' - pH, Ca, u/s DOC and discharge DOC')

      sum1 = osterror*osterror + seph2*seph2
      sum1 = SQRoot( 117333,sum1) ! compounded standard error ------------------
      sum1a = sum1 * SQRoot( 117333, ecns ) ! biggest standard deviation -------
      sum2 = osterror*osterror + seca2*seca2 
      sum2 = SQRoot( 117333,sum2) ! compounded standard error ------------------
      sum2a = sum2 * SQRoot( 117333, ecns ) ! biggest standard deviation -------
      sum3 = osterror*osterror + sudoc2*sudoc2
      sum3 = SQRoot( 117333,sum3) ! compounded standard error ------------------
      sum3a = sum3 * SQRoot( 117333, ecns ) ! biggest standard deviation -------
      sum4 = osterror*osterror + sedoc2*sedoc2 
      sum4 = SQRoot( 117333,sum4) ! compounded standard error ------------------
      sum4a = sum4 * SQRoot( 117333, ecns ) ! biggest standard deviation -------

      write(6,456)osterror,tecsplus,
     &slph,suph,seph2,sum1,sum1a,
     &slca,suca,seca2,sum2,sum2a,
     &slud,suud,sudoc2,sum3,sum3a,
     &sled,sued,sedoc2,sum4,sum4a
  456 format(95('-')/'standard errors'/95('-')/
     &19x,'lower','    upper','  average','      sum',
     &' sd.error'
     &/95('-')/
     &'dissolved data:',27x,2f9.3/
     &'            pH:',5f9.3/
     &'       calcium:',5f9.3/
     &'       u/s DOC:',5f9.3/
     &' discharge DOC:',5f9.3/95('-'))
      
      sum = osterror*osterror + seph2*seph2 + seca2*seca2 + 
     &sudoc2*sudoc2 + sedoc2*sedoc2 
      sum = SQRoot( 117333,sum) ! compounded standard error --------------------
      tecsplus = sum * SQRoot( 117333, ecns ) ! standard deviation -------------
      write(6,556)sum,tecsplus
  556 format('           SUM:',27x,2f9.3,/95('='))
      !write(03,8831)tecm,tecs,osterror,seph2,seca2,sudoc2,sedoc2,
      !  &tecsplus,tecmaster
 8831 format(/50('+')/
     &'     tecm =',f12.4/ 
     &'     tecs =',f12.4/ 
     &' osterror =',f12.4/ 
     &'    seph2 =',f12.4/ 
     &'    seca2 =',f12.4/ 
     &'   sudoc2 =',f12.4/ 
     &'   sedoc2 =',f12.4/ 
     &' tecsplus =',f12.4/ 
     &'tecmaster =',f12.4/ 
     &50('+')/)

      sum = bsterror*bsterror + sbph2*sbph2 + sbca2*sbca2 + 
     &sbud2*sbud2 + sbed2*sbed2 
      sum = SQRoot( 117333,sum) ! compounded standard error --------------------
      tbsplus = sum * SQRoot( 117333, tbnsc ) ! standard deviation -------------

      usterror = ubs / SQRoot( 117333, rcns ) ! original standard error ==
      sum = usterror*usterror + suph2*suph2 + suca2*suca2 + 
     &suud2*suud2 
      sum = SQRoot( 117333,sum) ! compounded standard error --------------------
      ubsplus = sum * SQRoot( 117333, rcns ) ! standard deviation --------------

      vsterror = vbs / SQRoot( 117333, rcns ) ! original standard error ==
      sum = vsterror*vsterror + suph2*suph2 + suca2*suca2 + 
     &suud2*suud2 
      sum = SQRoot( 117333,sum) ! compounded standard error --------------------
      vbsplus = sum * SQRoot( 117333, rcns ) ! standard deviation --------------

      !write(03,4396)usterror,suph2,suca2,suud2,ubsplus
 4396 format(/50('*')/
     &'usterror = ',f10.5/
     &'   suph2 = ',f10.5/
     &'   suca2 = ',f10.5/
     &'   suud2 = ',f10.5/
     &' ubsplus = ',f10.5/50('*')/)

      !write(3,1456)bsterror,sbph2,sbca2,sbud2,sbed2,sum,tbsplus,
      !&tbscmaster,tbsplus/tbscmaster
 1456 format(/50('-')/'Standard Error ...'/50('-')/
     &'dissolved data:',f9.3,'  (bsterror)'/
     &'            pH:',f9.3/
     &'       calcium:',f9.3/'       u/s DOC:',f9.3/
     &' discharge DOC:',f9.3/50('-')/
     &'compounded standard error ...'/50('-')/
     &'           sum:',f9.3,'  (sum)'/50('=')/
     &'standard deviation ...'/50('=')/
     &'       tbsplus:',f9.3/
     &'   master tbcs:',f9.3/
     &'         ratio:',f9.3/50('=')/)
     
      tbsplus = tbsplus/tbscmaster
    
      !write(3,56)
   56 format(//68('=')/'USE THE ADDED ERRORS FOR ALL FOUR ITEMS',
     &' - pH, Ca, u/s DOC and discharge DOC')
      nerr = 2
      call set the starting values of the data    
      call do the calculations (0) ! do the calculations
      nerr = 1

      return
      end
      
      
      subroutine check type of metal
      include 'bacom1cp.for'

      biometal = 0
      if ( free .eq. 1 ) biometal = 2

      if ( mettal .eq. 1 ) then
      name of metal = 'Zinc ...'
      name of metal2 = '           Zinc'
      endif
      
      if ( mettal .eq. 2 ) then
      name of metal = 'Copper ...'
      name of metal2 = '         Copper'
      endif
     
      if ( mettal .eq. 3 ) then
      name of metal = 'Manganese ...'
      name of metal2 = '      Manganese'
      endif
     
      if ( mettal .eq. 4 ) then
      name of metal = 'Nickel ...'
      name of metal2 = '         Nickel'
      endif

      if ( mettal .eq. 5 ) then
      name of metal = 'Lead ...'
      name of metal2 = '           Lead'
      endif
     
      return
      end

      subroutine metal data
      include 'bacom1cp.for'
      
      if ( mettal .eq. 1 ) call zinc data
      if ( mettal .eq. 2 ) call copper data
      if ( mettal .eq. 3 ) call manganese data
      if ( mettal .eq. 4 ) call nickel data
      if ( mettal .eq. 5 ) call lead data

      return
      end