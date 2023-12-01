*     approximation to cumulative Normal function...

      function CUMNOR(U)

      Y = SQRT(2/3.14159)*U*(1 + 0.044715*U**2)
      CUMNOR = EXP(2*Y)/(1 + EXP(2*Y))

      return
      end


*     approximation to Inverse Normal function...

      function fnor (CumP)

      t = SQRT(-LOG(4*CumP*(1.-CumP)))
      fnor = 1.238*t*(1 + 0.0262*t)
      IF(CumP.LT.0.5) fnor = -fnor

      return
      end



*     Cumulative probability and percentage point of Student's
*     T-Distribution.
*                  ... TWO-TAIL TEST ...
    
      FUNCTION TDIST1(EN,P)
      Q=P
      IF(P.LT.0.5)Q=1.0-P
      Q=2.0*Q-1.0
      TDIST1=TDIST(1,EN,Q)
      IF(P.LT.0.5)TDIST1=-TDIST1
      RETURN
      END

      FUNCTION TDIST2(EN,T)
      P=TDIST(0,EN,T)
      P=0.5*(1.0+P)
      IF(T.LT.0.0)P=1.0-P
      TDIST2=P
      RETURN
      END
 
      FUNCTION TDIST(IX,A,X)
      DATA A2/.5/
      A1=A/2.
      IF(IX.EQ.1)GO TO 20
      Y=X*X
      Y=Y/(A+Y)
      TDIST=BETA(0,Y,A2,A1)
      RETURN
   20 IF(ABS(X).LT.0.5) GO TO 40
      FX=BETA(1,1.-X,A1,A2)
      FX=AMAX1(FX,A*1.E-30)
      TDIST=SQRT((1./FX-1.)*A)
      RETURN
   40 FX=BETA(1,X,A2,A1)
      TDIST=SQRT((FX/(1.-FX))*A)
      RETURN
      END
    
    
    
    
*     INCOMPLETE BETA FUNCTION AND ITS INVERSE
*     ****************************************
    
      FUNCTION BETA(IND,X,A,B)
      IF(X.GT.0.0) GO TO 4
      BETA=0.0
      RETURN
    4 IF(X.LT.1.0) GO TO 6
      BETA=1.0
      RETURN
    6 CONTINUE
      CAB=CGAM(A+B)-CGAM(A)-CGAM(B)-.5*ALOG((A+B)*6.28318531)
      IF(IND)10,10,20
   10 EP=CAB+A*ALOG(X*(1.+B/A))+B*ALOG((1.-X)*(1.+A/B))
      IF(X-A/(A+B))12,12,14
   12 BETA=ZI(X,A,B)*EXP(EP+.5*ALOG(B/A))
      RETURN
   14 BETA=1.-ZI(1.-X,B,A)*EXP(EP+.5*ALOG(A/B))
      RETURN
   20 IF(X-.5)22,22,24
   22 QZ=ALOG(X)
      IGO=1
      AA=A
      BB=B
      GOTO 26
   24 QZ=ALOG(1.0-X)
      IGO=2
      AA=B
      BB=A
   26 XT=AA/(AA+BB)
      CABB=CAB+.5*ALOG(BB/AA)+AA*ALOG(1.+BB/AA)+BB*ALOG(1.+AA/BB)
      DO 40 NC=1,100
      ZZ=ZI(XT,AA,BB)
      QX=CABB+AA*ALOG(XT)+BB*ALOG(1.-XT)+ALOG(ZZ)
      XC=(QZ-QX)*(1.0-XT)*ZZ/AA
      XC=AMAX1(XC,-0.99)
      XC=AMIN1(XC,0.5/XT-0.5)
      XT=XT*(1.0+XC)
      IF(ABS(XC)-1.E-6)42,40,40
   40 CONTINUE
   42 GO TO(44,46),IGO
   44 BETA=XT
      RETURN
   46 BETA=1.-XT
      RETURN
      END
    
    
    
    
*     AUXILLARY SUBPROGRAMS FOR 'BETA'
*     ********************************
    
      FUNCTION ZI(X,A,B)
      FN=.7*(ALOG(15.+A+B))**2+AMAX1(X*(A+B)-A,0.)
      N=INT(FN)
      C=1.-(A+B)*X/(A+2.*FN)
      ZI=2./(C+SQRT(C**2-4.*FN*(FN-B)*X/(A+2.*FN)**2))
      DO 60 J=1,N
      FN=N+1-J
      A2N=A+2.*FN
      ZI=(A2N-2.)*(A2N-1.-FN*(FN-B)*X*ZI/A2N)
      ZI=1./(1.-(A+FN-1.)*(A+FN-1.+B)*X/ZI)
   60 CONTINUE
      RETURN
      END
    
      FUNCTION CGAM(A)
      AA=A
      CAC=0.0
      IF(A-2.)2,8,8
    2 IF(A-1.)4,6,6
    4 CAC=-2.+(A+.5)*ALOG(1.+1./A)+(A+1.5)*ALOG(1.+1./(A+1.))
      AA=A+2.
      GOTO 8
    6 CAC=-1.+(A+.5)*ALOG(1.+1./A)
      AA=A+1.0
    8 CA=2.269489/AA
      CA=0.52560647/(AA+1.0115231/(AA+1.5174737/(AA+CA)))
      CA=.083333333/(AA+.033333333/(AA+.25238095/(AA+CA)))
      CGAM=CA+CAC
      RETURN
      END



      subroutine LNCL (dcm,dcs,n,per,dcp,dcpl,dcpu,normal) ! LNCL

      en = n
      
      if ( dcs .lt. 0.00001 ) then
          dcp = dcm
          dcpl = dcm
          dcpl = dcm
          return
      endif
          
      if ( normal .eq. 0) then ! the distribution is log-normal
      
*     mean of logged variables
      GDCM=ALOG(DCM/(SQRT(1.0+DCS*DCS/(DCM*DCM))))
    
*     standard deviation of the logged variables
      GDCS=SQRT(ALOG(1.0+DCS*DCS/(DCM*DCM)))
    
      T0 = errx (0.01*per )
      DCP=EXP(GDCM+T0*GDCS) ! percentile assuming a log-normal distribution
      TL=TSHIFT(EN,T0,-1.0)
      TU=TSHIFT(EN,T0,1.0)
      DCPL=EXP(GDCM+TL*GDCS) ! lower confidence limit
      DCPU=EXP(GDCM+TU*GDCS) ! upper confidence limit
      else ! the distribution is normal
      GDCM=DCM ! mean
      GDCS=DCS ! standard deviation
      T0 = errx (0.01*per) 
      DCP=GDCM+T0*GDCS ! 95-percentile assuming a normal distribution
      TL=TSHIFT(EN,T0,-1.0)
      TU=TSHIFT(EN,T0,1.0)
      DCPL=GDCM+TL*GDCS ! lower confidence limit
      DCPU=GDCM+TU*GDCS ! lower confidence limit
      endif  
      
      return
      END

      
      subroutine LNCL1 (dcm,dcs,en,per,dcp,dcpl,dcpu,normal) ! LNCL1

      if ( dcm .lt. 0.0000001 ) then
          dcp = dcm
          dcpl = dcm
          dcpl = dcm
          return
      endif

      if ( dcs .lt. 0.00001 ) then
          dcp = dcm
          dcpl = dcm
          dcpl = dcm
          return
      endif

      if ( normal .eq. 0) then ! the distribution is log-normal
      
*     mean of logged variables
      GDCM=ALOG(DCM/(SQRT(1.0+DCS*DCS/(DCM*DCM))))
*     standard deviation of the logged variables
      GDCS=SQRT(ALOG(1.0+DCS*DCS/(DCM*DCM)))
    
      T0 = errx (0.01*PER )
      DCP=EXP(GDCM+T0*GDCS) ! percentile assuming a log-normal distribution
      TL=TSHIFT(EN,T0,-1.0)
      TU=TSHIFT(EN,T0,1.0)
      DCPL=EXP(GDCM+TL*GDCS) ! lower confidence limit
      DCPU=EXP(GDCM+TU*GDCS) ! upper confidence limit
      
      else ! the distribution is normal
      GDCM=DCM ! mean
      GDCS=DCS ! standard deviation
      T0 = errx (0.01*PER) 
      DCP=GDCM+T0*GDCS ! 95-percentile assuming a normal distribution
      TL=TSHIFT(EN,T0,-1.0)
      TU=TSHIFT(EN,T0,1.0)
      DCPL=amax1(0.0,GDCM+TL*GDCS) ! lower confidence limit
      DCPU=GDCM+TU*GDCS ! lower confidence limit
      endif  
      
      return
      END
     

*     normal Probability Function
      function Errx (x)
      p = 2.0 * x - 1.0
      if ( x .lt. 0.5 ) p = -p
      xx = errinv ( p ) * 1.41421
      Errx = sign ( xx, x - 0.5 )
      return
      end


*     inverse Error Function
      function errinv ( c )
      double precision y(7)

      data Y/0.00D+00,0.842700793,0.995322265,0.999977910,0.999999984,
     &       1.00D+00,1.00D+00/

      if ( c .ge. 1.0 ) then
        errinv = 6.0
        goto 99
      endif

      if ( c .le. 0.0 ) then
        errinv = 0.00001
        goto 99
      endif

      do 50 ii = 1, 7
      i = ii
      if ( Y(ii) .eq. c ) then
        errinv = float ( i - 1 )
        goto 99
      endif
      if ( Y(ii) .gt. c ) goto 70
   50 continue

   70 xc = i - 2

      do 90 k = 1, 20
      A = errf ( xc )
      temp = xc + (c - A) * ( 0.886226925 * exp ( xc**2 ) )
      B = errf ( temp ) 
      z = c - B
      xc = temp
      if ( z .lt. 1.0e-10 ) goto 80
   90 continue

   80 errinv = temp

   99 return
      end


*     Return a T-statistic from the Shifted T-distribution.
      FUNCTION TSHIFT(EN,KPROB,SIGN)
      REAL LAMBDA,KPROB
    
      EF=EN-1.0 ! degrees of freedom
      RN=SQRT(EN)
      EF2=2.0*EF
    
*     NON-CENTRALITY PARAMETER FOR SHIFTED T-DISTRIBUTION
    
      DELTA=RN*KPROB
    
*     VARIABLE FOR LOOK-UP TABLE FOR LAMBDA
    
      ETA=SIGN*DELTA/SQRT(EF2*(1.0+DELTA*DELTA/EF2))
      LAMBDA=TABLE(ETA,EF)
      TSHIFT=SIGN*(SIGN*DELTA+LAMBDA*SQRT(1.0+DELTA*DELTA/EF2
     #      -LAMBDA*LAMBDA/EF2))/(RN*(1.0-LAMBDA*LAMBDA/EF2))
      
      return
      end
      
      
      FUNCTION TABLE(X,Y)
      DIMENSION U(25),V(11),A(11,25)
      DATA V/2.,3.,4.,5.,6.,7.,8.,9.,16.,36.,144./
      DATA U/-1.00,-0.95,-0.90,-0.85,-0.80,-0.75,-0.70,-0.65,-0.60,
     #-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,
     #0.9,1.0/
      DATA A/1.4616,1.5039,1.5277,1.5431,1.5542,1.5625,1.5691,
     #       1.5744,1.5952,1.6141,1.6307,
     #       1.4648,1.5037,1.5262,1.5411,1.5519,1.5601,1.5666,
     #       1.5720,1.5930,1.6124,1.6298,
     #       1.4798,1.5108,1.5302,1.5435,1.5534,1.5610,1.5671,
     #       1.5721,1.5924,1.6116,1.6293,
     #       1.5047,1.5242,1.5388,1.5497,1.5580,1.5646,1.5700,
     #       1.5745,1.5931,1.6115,1.6290,
     #       1.5362,1.5422,1.5511,1.5588,1.5651,1.5704,1.5749,
     #       1.5786,1.5951,1.6121,1.6290,
     #       1.5707,1.5631,1.5658,1.5700,1.5742,1.5779,1.5813,
     #       1.5843,1.5980,1.6133,1.6293,
     #       1.6049,1.5853,1.5819,1.5826,1.5845,1.5867,1.5889,
     #       1.5910,1.6017,1.6149,1.6297,
     #       1.6364,1.6075,1.5986,1.5960,1.5957,1.5963,1.5973,
     #       1.5985,1.6060,1.6170,1.6304,
     #       1.6631,1.6285,1.6152,1.6096,1.6072,1.6063,1.6062,
     #       1.6065,1.6108,1.6194,1.6312,
     #       1.6984,1.6641,1.6456,1.6357,1.6300,1.6265,1.6244,
     #       1.6230,1.6213,1.6251,1.6333,
     #       1.7101,1.6888,1.6701,1.6582,1.6505,1.6453,1.6417,
     #       1.6390,1.6322,1.6313,1.6359,
     #       1.7045,1.7023,1.6874,1.6758,1.6675,1.6614,1.6569,
     #       1.6535,1.6437,1.6378,1.6387,
     #       1.6892,1.7065,1.6979,1.6883,1.6805,1.6744,1.6696,
     #       1.6658,1.6525,1.6442,1.6417,
     #       1.6700,1.7040,1.7025,1.6959,1.6894,1.6839,1.6793,
     #       1.6755,1.6611,1.6503,1.6447,
     #       1.6501,1.6970,1.7024,1.6994,1.6947,1.6902,1.6862,
     #       1.6828,1.6682,1.6558,1.6477,
     #       1.6312,1.6876,1.6991,1.6995,1.6969,1.6937,1.6905,
     #       1.6875,1.6739,1.6607,1.6504,
     #       1.6141,1.6769,1.6934,1.6970,1.6965,1.6947,1.6924,
     #       1.6901,1.6780,1.6647,1.6529,
     #       1.5990,1.6660,1.6864,1.6927,1.6941,1.6936,1.6923,
     #       1.6907,1.6807,1.6678,1.6551,
     #       1.5861,1.6553,1.6785,1.6872,1.6902,1.6909,1.6906,
     #       1.6897,1.6819,1.6699,1.6568,
     #       1.5751,1.6452,1.6705,1.6808,1.6853,1.6870,1.6875,
     #       1.6873,1.6818,1.6711,1.6581,
     #       1.5661,1.6361,1.6625,1.6741,1.6796,1.6823,1.6834,
     #       1.6838,1.6805,1.6714,1.6588,
     #       1.5590,1.6280,1.6549,1.6673,1.6736,1.6769,1.6786,
     #       1.6795,1.6782,1.6707,1.6590,
     #       1.5536,1.6210,1.6479,1.6607,1.6674,1.6712,1.6734,
     #       1.6746,1.6750,1.6692,1.6586,
     #       1.5497,1.6153,1.6416,1.6544,1.6613,1.6654,1.6678,
     #       1.6693,1.6710,1.6667,1.6576,
     #       1.5470,1.6160,1.6362,1.6487,1.6556,1.6597,1.6622,
     #       1.6638,1.6664,1.6635,1.6560/
    
      IU=1
      DO 10 I=1,25
      IF(U(I).GT.X)GOTO 11
      IU=I
   10 CONTINUE
   11 JU=IU+1
    
      IF(Y.GT.9.0)GOTO 12
      IV=INT(Y)-1
      GOTO 13
   12 IV=1
      DO 14 I=8,11
      IF(V(I).GT.Y)GOTO 17
      IV=I
   14 CONTINUE
    
   13 IF(IU.GT.1.AND.IU.LT.25)GOTO 16
      TABLE=A(IU,IV)
      RETURN
   16 P=(X-U(IU))/(U(JU)-U(IU))
      TABLE=(1.0-P)*A(IV,IU)+P*A(IV,JU)
      RETURN
    
   17 JV=IV+1
      A00=A(IV,IU)
      A01=A(IV,JU)
      A10=A(JV,IU)
      A11=A(JV,JU)
    
      P=(X-U(IU))/(U(JU)-U(IU))
      T=(12.0/SQRT(Y)-12.0/SQRT(V(IV)))/(12.0/SQRT(V(JV))-
     #12.0/SQRT(V(IV)))
      TABLE=(1.0-P)*(1.0-T)*A00+T*(1.0-P)*A10+T*P*A11+(1.0-T)*P*A01
    
      RETURN
      END

      
*     calculations for the errors function
      function errf (w)
      double precision a(25), b(30), f

      data a/16443152242714.D-13,-9049760497548.D-13,
     &643570883797.D-13,196418177368.D-13,-1244215694.D-13,
     &-9101941905.D-13,-1796219835.D-13,139836786.D-13,
     &164789417.D-13,39009267.D-13,-893145.D-13,-3747896.D-13,
     &1298818.D-13,136773.D-13,77107.D-13,46810.D-13,
     &11844.D-13,-5.D-13,-1384.D-13,-652.D-13,145.D-13,
     &10.D-13,24.D-13,11.D-13,2.D-13/

      m = 24
      x = abs ( w )

      if ( x - 0.01 ) 1,2,2
    1 xerr = 2.0 / ( 3.0 * 1.77245385 ) * x * ( 3.0 - x**2 )
      goto 6

    2 z = ( x - 1.0 ) / ( x + 1.0 )

      do 3 i = 1, 30
      b(i) = 0.0
    3 continue

      do 4 i = 1, m

      m1 = ( m + 1 ) - i
      b ( m1 ) = 2.0 * z * b( m1 + 1 )-b( m1 + 2 ) + a( m1 + 1 )
    4 continue

      f = -b(2) + z * b(1) + 0.5 * a(1)
      xerr = 1.0 - ( 1.0 / 1.77245385 ) * (exp (-(x**2))) * f

      if ( x - 0.01 ) 6,7,7
    6 cerr = 1.0 -xerr
      goto 5

    7 cerr = ( 1.0/1.77245385 ) * (exp(-(x**2))) * f
    5 if ( w ) 9,8,8
    8 errf = xerr
      goto 13

    9 errf = cerr
   13 return
      end

********************************************************************************

      
      subroutine CONAV(dcm,dcs,n1,dcml,dcmu)    
      E1=N1
      T0=TDIST1(E1-1.0,0.95)
      S1=SQRT(dcs*dcs/E1)
      dcml=amax1(0.0,dcm-S1*T0)
      dcmu=dcm+S1*T0
      return
      end
      
      subroutine CONAV1(dcm,dcs,e1,dcml,dcmu) 
      T0=TDIST1(E1-1.0,0.95)
      S1=SQRT(dcs*dcs/E1)
      dcml=amax1(0.0,dcm-S1*T0)
      dcmu=dcm+S1*T0
      return
      end

      subroutine confidence in standard deviation (std,Sd,sn,XL,XU)
      xn = sn - 1.0
      SES = Sd / sqrt ( 2.0 * xn )
      t = errx ( 0.95 )
      XL = amax1 (0.0, (std - t * SES) )
      XU = amax1 (0.0, (std + t * SES) )
	return
      end

      subroutine CONSDEV (std,n1,XL,XU)
      xn = n1 - 1.0
      SES = std / sqrt ( 2.0 * xn )
      t = errx ( 0.95 )
      XL = amax1 (0.0, (std - t * SES) )
      XU = amax1 (0.0, (std + t * SES) )
	return
      end
      
      subroutine CONSDEV1 (std,sn,XL,XU)
      xn = sn - 1.0
      SES = std / sqrt ( 2.0 * xn )
      t = errx ( 0.95 )
      XL = amax1 (0.0, (std - t * SES) )
      XU = amax1 (0.0, (std + t * SES) )
	return
      end
      
      subroutine CONSDEV2 (std,std2,sn,XL,XU)
      xn = sn - 1.0
      SES = std2 / sqrt ( 2.0 * xn )
      t = errx ( 0.95 )
      XL = amax1 (0.0, (std - t * SES) )
      XU = amax1 (0.0, (std + t * SES) )
	return
      end
    

      subroutine compute double precision summary statistics (XM,XS)
      include 'bacom1cp.for'
      double precision xm,xs,xxx,xns
      
      xns = NS

      XXX = XS-XM*XM/xns
      if (XXX .le. Small) then
      XS = 0.0
      else
      XS = SQRT (XXX/(xns-1.0)) ! standard deviation
      endif
      
      XM = XM / xns ! mean

      return
      end

      subroutine compute ordinary summary statistics (XM,XS)
      include 'bacom1cp.for'
      real xm,xs

      XXX = XS-XM*XM/FLOAT(NS)
      if (XXX .le. Small) then
      XS = 0.0
      else
      XS = SQRT (XXX/FLOAT(NS-1)) ! standard deviation
      endif
      XM = XM / FLOAT(NS) ! mean

      return
      end
     
      
*     compute square root and check for error ----------------------------------
      function SQRoot (Jdentity,X)
      if ( x .lt. 1.0E-20) X=1.0E-20
      if ( x .ge. 1.0e20 ) X=1.0E20  
      SQRoot = SQRT(X)
      return
      end

