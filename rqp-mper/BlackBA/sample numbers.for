*     --------------------------------------------------------------------------
*     random number generators   
*     --------------------------------------------------------------------------
*     generate random samples from normal distribution for river flow ----------
*     --------------------------------------------------------------------------

*     --------------------------------------------------------------------------
*     generate random samples from normal distribution for river quality -------
*     --------------------------------------------------------------------------
      function GAS2(IDUM)
      save iset,gset
      data ISET /0/
      if ( IDUM .lt. 0 ) ISET=0
      if ( ISET .eq. 0 ) then
    1 continue 
      V1=2.*RAN2(IDUM)-1.
      V2=2.*RAN2(IDUM)-1.
      R=V1**2+V2**2
      if (R .ge. 1.) goto 1
      FAC=SQRoot(1010,-2.*LOG(R)/R)
      GSET=V1*FAC
      GAS2=V2*FAC
      ISET=1
      else
      GAS2=GSET
      ISET=0
      endIF
      return
      end
      
*     --------------------------------------------------------------------------
*     generate random samples from normal distribution for discharge flow ------
*     --------------------------------------------------------------------------
      function GAS3(IDUM)
      save iset,gset
      data ISET /0/
      if (IDUM .lt. 0) ISET=0
      if (ISET .eq. 0) then
    1 continue     
      
      V1=2.*RAN3(IDUM)-1.
      V2=2.*RAN3(IDUM)-1.
      R=V1**2+V2**2
      if (R .ge. 1.) goto 1
      FAC=SQRoot(1012,-2.*LOG(R)/R)
      GSET=V1*FAC
      GAS3=V2*FAC
      ISET=1
      else
      GAS3=GSET
      ISET=0
      endIF
      return
      end
*     --------------------------------------------------------------------------
*     generate random samples from normal distribution for discharge quality ---
*     --------------------------------------------------------------------------
      function GAS4(IDUM)
      save iset,gset
      data ISET /0/
      if ( IDUM .lt. 0 ) ISET=0
      if ( ISET .eq. 0 ) then
    1 continue     
      V1=2.*RAN4(IDUM)-1.
      V2=2.*RAN4(IDUM)-1.
      R=V1**2+V2**2
      if (R .ge. 1.) goto 1
      FAC=SQRoot(1013,-2.*LOG(R)/R)
      GSET=V1*FAC
      GAS4=V2*FAC
      ISET=1
      else
      GAS4=GSET
      ISET=0
      endIF
      return
      end


      function GAS8(IDUM)
      save iset,gset
      data ISET /0/
      if ( IDUM .lt. 0 ) ISET=0
      if ( ISET .eq. 0 ) then
    1 continue     
      V1=2.*RAN8(IDUM)-1.0
      V2=2.*RAN8(IDUM)-1.0
      R=V1**2+V2**2
      if (R .ge. 1.0) goto 1
      FAC=SQRoot(1014,-2.0*LOG(R)/R)
      GSET=V1*FAC
      GAS8=V2*FAC
      ISET=1
      else
      GAS8=GSET
      ISET=0
      endIF
      return
      end
      
      function GAS9(IDUM)
      save iset,gset
      data ISET /0/
      if ( IDUM .lt. 0 ) ISET=0
      if ( ISET .eq. 0 ) then
    1 continue     
      V1=2.*RAN9(IDUM)-1.0
      V2=2.*RAN9(IDUM)-1.0
      R=V1**2+V2**2
      if (R .ge. 1.0) goto 1
      FAC=SQRoot(1014,-2.0*LOG(R)/R)
      GSET=V1*FAC
      GAS9=V2*FAC
      ISET=1
      else
      GAS9=GSET
      ISET=0
      endIF
      return
      end

*     --------------------------------------------------------------------------
*     Random number generators -------------------------------------------------  
*     --------------------------------------------------------------------------
*     generate random numbers for river flow -----------------------------------
*     --------------------------------------------------------------------------
      function RAN2(IDUM)
      dimension R(97)
      parameter (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      parameter (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      parameter (M3=243000,IA3=4561,IC3=51349)
      save IX1,IX2,IX3
      if ( IDUM .lt. 0 ) then
      IX1=MOD(IC1-IDUM,M1)
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IX1,M2)
      IX1=MOD(IA1*IX1+IC1,M1)
      IX3=MOD(IX1,M3)
      do J=1,97
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      enddo
      IDUM=1
      endIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      JJ=1+(97*IX3)/M3
      J = iabs (JJ)
      if (J .gt. 97 .or. J .lt. 1) then
      write(3,*)'Stopped in random number generation 2',JJ
      stop
      endif
      RAN2=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      return
      end
*     --------------------------------------------------------------------------
*     generate random numbers for discharge flow -------------------------------
*     --------------------------------------------------------------------------
      function RAN3(IDUM)
      dimension R(97)
      parameter (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      parameter (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      parameter (M3=243000,IA3=4561,IC3=51349)
      save IX1,IX2,IX3
      if ( IDUM .lt. 0 ) then
      IX1=MOD(IC1-IDUM,M1)
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IX1,M2)
      IX1=MOD(IA1*IX1+IC1,M1)
      IX3=MOD(IX1,M3)
      do J=1,97
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      enddo
      IDUM=1
      endIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      JJ=1+(97*IX3)/M3
      J = iabs (JJ)
      if (J .gt. 97 .or. J .lt. 1) then
      write(3,*)'Stopped in random number generation 3',JJ
      stop
      endif
      RAN3=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      return
      end
*     --------------------------------------------------------------------------
*     generate random numbers for discharge quality ----------------------------
*     --------------------------------------------------------------------------
      function RAN4(IDUM)
      dimension R(97)
      parameter (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      parameter (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      parameter (M3=243000,IA3=4561,IC3=51349)
      save IX1,IX2,IX3
      if ( IDUM .lt. 0 ) then
      IX1=MOD(IC1-IDUM,M1)
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IX1,M2)
      IX1=MOD(IA1*IX1+IC1,M1)
      IX3=MOD(IX1,M3)
      do 11 J=1,97
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11    continue
      IDUM=1
      endIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      JJ=1+(97*IX3)/M3
      J = iabs (JJ)
      if (J .gt. 97 .or. J .lt. 1) then
      write(3,*)'Stopped in random number generation 4',JJ
      endif
      RAN4=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      return
      end


      function RAN8(IDUM)
      dimension R(97)
      parameter (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      parameter (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      parameter (M3=243000,IA3=4561,IC3=51349)
      save IX1,IX2,IX3
      if ( IDUM .lt. 0 ) then
      IX1=MOD(IC1-IDUM,M1)
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IX1,M2)
      IX1=MOD(IA1*IX1+IC1,M1)
      IX3=MOD(IX1,M3)
      do J=1,97
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      enddo
      IDUM=1
      endIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      JJ=1+(97*IX3)/M3
      J = iabs (JJ)
      if (J .gt. 97 .or. J .lt. 1) then
      write(3,*)'Stopped in random number generation 8',JJ
      endif
      RAN8=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      return
      end

      function RAN9(IDUM)
      dimension R(97)
      parameter (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      parameter (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      parameter (M3=243000,IA3=4561,IC3=51349)
      save IX1,IX2,IX3
      if ( IDUM .lt. 0 ) then
      IX1=MOD(IC1-IDUM,M1)
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IX1,M2)
      IX1=MOD(IA1*IX1+IC1,M1)
      IX3=MOD(IX1,M3)
      do J=1,97
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      enddo
      IDUM=1
      endIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      JJ=1+(97*IX3)/M3
      J = iabs (JJ)
      if (J .gt. 97 .or. J .lt. 1) then
      write(3,*)'Stopped in random number generation 8',JJ
      stop
      endif
      RAN9=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      return
      end

      
      subroutine compute confidence limits d (perc,idist,XMd,XSd,XP,
     &qualn,XL,XU)
      !include 'bacom1cp.for'
      double precision XMd,XSd
      
      XL = 0.0 ! lower confidence limit
      XU = 0.0 ! upper confidence limit
      if ( idist .eq. 0 ) return
      A = XMd  ! mean
      S = XSd  ! standard deviation
*     check for a zero calculated mean -----------------------------------------
      if ( A .lt. small .or. S .lt. small ) then	
      S = 0.0
	XL = A   ! set confidence limits as the mean
	XU = A
	return
      endif

*     calculations for the mean ------------------------------------------------
      if ( perc .gt. -0.00001 .and. perc .lt. 0.00001 ) then
      xqualn = amax1 ( 0.05, QUALN )
      t = TDIST1(xqualn-1.0,0.95)
      SEM = S / SQRoot( 107333, xqualn ) ! standard error
*     calculate the 95% confidence limits about the mean -----------------------
      XL = amax1 (0.0, (A - t * SEM) )   ! lower confidence limit
      XU = amax1 (0.0, (A + t * SEM) )   ! upper confidence limit
      return
      endif ! end of calculations for the mean ---------------------------------

*     calculations for the 50-percentile ---------------------------------------
      if ( perc .gt. 49.99999 .and. perc .lt. 50.00001 ) then
      t = errx ( 0.95 )
      xqualn = amax1 ( 0.05, QUALN )
      SEM = S / SQRoot( 107333, xqualn ) ! standard error
      XL = amax1 (0.0, ( XP - t * SEM) )
      XU = amax1 (0.0, ( XP + t * SEM) )
      return
      endif ! calculations for the 50-percentile -------------------------------

*     log-normal distribution and percentiles ----------------------------------
      if ( idist .eq. 2 .or. idist .eq. 3 .or. idist .eq. 4) then
      GM = alog ( (A*A) / SQRoot(10778,A*A+S*S) )
      GS = SQRoot(100775, alog ( 1.0 + (S*S) / (A*A) ) )
      if ( perc .ge. 50.1 .and. perc .lt. 99.9999 ) then
      t = errx ( 0.01 * perc )
*     calculate the percentile assuming a log-normal distribution --------------
      DP = exp ( GM + t * GS )
      TL = tshift ( QUALN, t, -1.0)
      TU = tshift ( QUALN, t,  1.0)
*     calculate the confidence limits ------------------------------------------
      XL = amax1 (0.0, XP - (DP - exp ( GM + TL * GS)) )
      XU = amax1 (0.0, XP + (exp ( GM + TU * GS) - DP) )
      if ( XL .gt. 999928.0 .and. XL .lt. 28.5 ) then
          write(03,*)'A  = ',A
          write(03,*)'N  = ',qualn
          write(03,*)'S  = ',S
          write(03,*)'GM = ',GM,GS
          write(03,*)'TU = ',TU,t
          write(03,*)'XP = ',XP
          write(03,*)'DP = ',DP
          write(03,*)'Xl = ',XL
      endif
      endif
*     low percentiles ----------------------------------------------------------
      if ( perc .gt. 0.0001 .and. perc .lt. 49.9999 ) then
      t = - errx ( 0.01 * perc )
      DP = exp ( GM - t * GS )
      TL = -tshift ( QUALN, t,  1.0)
      TU = -tshift ( QUALN, t, -1.0)
      XL = amax1 (0.0, XP - (DP - exp ( GM + TL * GS)) )
      XU = amax1 (0.0, XP + (exp ( GM + TU * GS) - DP) )
      endif
      return
      endif 

*     normal distributions -----------------------------------------------------
      if ( idist .eq. 1 ) then
      GM = A
      GS = S
      if ( perc .ge. 49.9999 .and. perc .lt. 99.9999 ) then
      t = errx ( 0.01 * perc )
*     calculate the percentile assuming a normal distribution ------------------
      DP = GM + t * GS
*     calculate confidence limits ----------------------------------------------
      TL = tshift ( QUALN, t, -1.0)
      TU = tshift ( QUALN, t,  1.0)
      XL = amax1 (0.0, XP - (DP - ( GM + TL * GS)) )
      XU = amax1 (0.0, XP + (( GM + TU * GS) - DP) )
      endif
*     low percentiles ----------------------------------------------------------
      if ( perc .gt. 0.0001 .and. perc .lt. 49.9999 ) then
      t = - errx ( 0.01 * perc )
      DP = GM - t * GS
      TL = -tshift ( QUALN, t,  1.0)
      TU = -tshift ( QUALN, t, -1.0)
      XL = amax1 (0.0, XP - (DP - ( GM + TL * GS)) )
      XU = amax1 (0.0, XP + (( GM + TU * GS) - DP) )
      
	endif
      endif 

      return
      end
*     GENERATE NON-PARAMETRIC DATA FOR FLOW AND QUALITY ------------------------
*     non-parametric distribution of river flow --------------------------------
      
*     provide a random sample from a non-parametric distribution ---------------
      subroutine get non parametric shot (tnvar,val,im)
      include 'bacom1cp.for'
      integer U
      
      np = npdp (im)
      
*     tnvar   is the random normal deviate -------------------------------------
*     np      is the number of data points entered -----------------------------
*     npval   is an array of the data points entered ---------------------------
*     npcum   is an array of cumulative frequencies that matches these ---------
*     unvar   is the uniform random number (lying between zeo and 1.0) ---------
*     val     is the returned value that corresponds to "tnvar" ----------------

*     convert 'tnvar' from uniform deviate to normal deviate "unvar" -----------
      unvar = tdist2 ( 999999.0, tnvar )
      
*     prepare to deal with values beyond the first and the last values entered -
*     as data ... we shall interpolate between the last point and a point half -
*     the last interval beyond the last point

*     calculate the first interval ---------------------------------------------
      CLINT1 = npval(im,2) - npval(im,1)

*     TAU is the minimum - the first value less half the interval --------------
      TAU = npval(im,1) - CLINT1 / 2.0

*     the last interval --------------------------------------------------------
      CLINT2 = npval(im,np) - npval(im,np-1)

*     THETA is the maximum - last value plus half the interval
      THETA = npval(im,np) + CLINT2 / 2.0

*     calculate the gradients of these tails of the distribution
      GRAD1 = ( npval(im,1) - TAU  ) /  npcum(im,1)
      GRAD2 = ( THETA - npval(im,np) ) / (1.0 - npcum(im,np) )

*     locate this point on the cumulative frequency distribution.
      if ( unvar .le. npcum(im,1) ) then
      val = GRAD1 * unvar + TAU
      else if ( unvar .ge. npcum(im,np) ) then
      val = GRAD2 * unvar + THETA - GRAD2
      else

      CL = float(np)
      LL = 1
      U = np
      i10 = 0

   10 i10 = i10 + 1

      if ( i10 .eq. mprf ) then
      write(3,*)'Error in random sampling ...'
      stop    
      endif

      cxx = 0.5 * float (U + LL)
      I = INT (cxx + 0.5)
      if ( unvar .le. npcum(im,I) .and. 
     &unvar .ge. npcum(im,I-1) ) then
      Y1 = npval(im,I-1)
      Y2 = npval(im,I)
      X1 = npcum(im,I-1)
      X2 = npcum(im,I)
      val = (unvar - X1) * (Y2 - Y1) / (X2 - X1) + Y1
      else
      if ( unvar .gt. npcum(im,I) ) LL = I
      if ( unvar .lt. npcum(im,I) ) U = I
      goto 10
      endif
      
      endif

      val = amax1 ( 0.0, val )

      return
      end

      
      
