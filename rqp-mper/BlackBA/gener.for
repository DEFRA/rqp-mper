*     ==========================================================================
*     SIMCAT - Monte-Carlo simulation of water quality in a river ....
*     ==========================================================================
      


*     generate shots for river quality ... Log-Normal Distribution -------------
      subroutine generate log normal river quality
      include 'bacom1cp.for'

      GRCM = 0.0
      GRCS = 0.0
      RC3 = 0.0 ! shift --------------------------------------

      if ( RCS .lt. 1.0E-08 ) then
      do IS = 1, NS
      CMS(JP,IS) = RCM
      enddo
      goto 56
      endif

      RC3 = 0.0 ! shift quality ------------------------------------------------
      RM3 = (RCM + RC3) * (RCM + RC3)
      spcorrRFRC = cofc1 ! river flow on river quality --- QRAN
      GRCM=ALOG(RM3/SQRoot(1003,(RM3+RCS*RCS))) ! mean of logged data ----------
      GRCS=SQRoot(1004,ALOG(1.0+(RCS*RCS)/RM3)) ! standard deviation -----------

*     calculate Monte-Carlo sampling errors ------------------------------------
      call bias in log normal river quality 

      if ( BM(2) .gt. 1.0E-08 ) BM(2) = RCM/BM(2)
      if ( BS(2) .gt. 1.0E-08 ) BS(2) = RCS/BS(2)
      if ( MONQ .gt. 1 ) then
*	if ( nobigout .le. 0 ) write(01,12)BM(2),BS(2)
   12 format(77('-')/
     &'Corrections for Monte-Carlo sampling: log-normal quality'/
     &77('-')/'Mean          =',F8.4/
     &'95-percentile =',F8.4/77('-'))
      endif ! if ( MONQ .gt. -1 ) then
      
*     sample the distributions -------------------------------------------------
      do is = 1, NS
      js = JSET ( IQ, NS, is )
      qdln = QRAN function ( JP, js, is )
      CMS(JP,is) = Vlogn (qdln,GRCS,GRCM,RC3,BM(2),BS(2),RCM)
      enddo

   56 continue
      if ( MONQ .gt. 1 ) call write generated shots for river quality
      return
      end




*     generate quality shots for head of reach ------------- normal distribution
      subroutine generate normal river quality
      include 'bacom1cp.for'

*     mean and standard deviation ----------------------------------------------
      GRCM = 0.0
      GRCS = 0.0
      GRCM = rcm
      if ( GRCM .lt. 1.0E-8 ) goto 56
      GRCS = rcs

      if ( GRCS .lt. 1.0E-8 ) then
      do IS=1,NS
      CMS(JP,IS) = AMAX1(1.0e-8,GRCM)
      enddo
      goto 56
      endif

      RC3 = 0.0 !  shift quality -----------------------------------------------
      spcorrRFRC = cofc1 ! river flow on river quality -- QRAN

*     calculate Monte-Carlo sampling errors ------------------------------------
      call bias in normal river quality 

      if ( BM(2) .gt. 1.0E-08 ) BM(2) = GRCM/BM(2)
      if ( BS(2) .gt. 1.0E-08 ) BS(2) = GRCS/BS(2)

      if ( MONQ .gt. 1 ) then
	if ( nobigout .le. 0 ) write(01,12)BM(2),BS(2)
   12 format(77('-')/
     &'Corrections for Monte-Carlo sampling errors: normal quality'/
     &77('-')/'Mean               =',F8.3/
     &'Standard deviation =',F8.3/77('-'))
	endif

   56 continue

      return
      end






*     compute values of randomly selected variables ----------------------------
      function Vlogn (R,GS,GM,RC33,BM,BS,RCM)
      XVAL = EXP ( GM + GS * R ) - RC33
      XVlogn = ( BM * XVAL - RCM ) * BS + RCM
      if ( XVlogn .lt. (0.00001*RCM) ) XVlogn = 0.00001*RCM
      if ( Xval .gt. 0.0 ) then
      if ( XVlogn .gt. (2.0 * XVAL) ) XVlogn = XVAL
      endif
      Vlogn = Xvlogn
      return
      end


      function JSET ( jstart, jmax, ival )
*     jstart is set to IQ, jmax is usually NS, IVAL is element (IS)
      jj = jstart - 1 + ival

    2 continue
      if ( jj .le. jmax ) goto 1
      jj = jj - jmax
      goto 2
    1 continue

      jset = jj

      return
      end


*     get random normal deviate for river flow  --------------------------------
      function FRAN function (JS)
      include 'bacom1cp.for'
      FRAN function = FRAN (JS)
      return
      end

      
      
*     get random normal deviate for river quality ------------------------------
*     imposing the correct correlation -----------------------------------------
      function QRAN function ( JDET,JS,IS )
      include 'bacom1cp.for'
      itest = 0
*     check there is no special correlation ------------------------------------
      if ( cofc1 .lt.  0.001 .and. cofc1 .gt. -0.001 ) then
      xx = CRAN (JDET,JS)
      else
      RRS = FRAN (IS)
      xx = cofc1 * RRS + CRAN (JDET,JS) * SQRMB(115, 1.0 - cofc1*cofc1 )
      endif
      QRAN function = xx
      return
      end

