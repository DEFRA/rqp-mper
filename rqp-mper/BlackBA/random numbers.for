*     --------------------------------------------------------------------------
*     random number generators   
*     --------------------------------------------------------------------------
*     generate random samples from normal distribution for river flow ----------
*     --------------------------------------------------------------------------

*     --------------------------------------------------------------------------
*     generate random samples from normal distribution for river quality -------
*     --------------------------------------------------------------------------
      function GAS2s(IDUM)
      save iset,gset
      data ISET /0/
      if ( IDUM .lt. 0 ) ISET=0
      if ( ISET .eq. 0 ) then
    1 continue 
      V1=2.*RAN2s(IDUM)-1.
      V2=2.*RAN2s(IDUM)-1.
      R=V1**2+V2**2
      if (R .ge. 1.) goto 1
      FAC=SQRoot(1010,-2.*LOG(R)/R)
      GSET=V1*FAC
      GAS2s=V2*FAC
      ISET=1
      else
      GAS2s=GSET
      ISET=0
      endIF
      return
      end
      
*     --------------------------------------------------------------------------
*     generate random samples from normal distribution for discharge flow ------
*     --------------------------------------------------------------------------
      function GAS3s(IDUM)
      save iset,gset
      data ISET /0/
      if (IDUM .lt. 0) ISET=0
      if (ISET .eq. 0) then
    1 continue     
      
      V1=2.*RAN3s(IDUM)-1.
      V2=2.*RAN3s(IDUM)-1.
      R=V1**2+V2**2
      if (R .ge. 1.) goto 1
      FAC=SQRoot(1012,-2.*LOG(R)/R)
      GSET=V1*FAC
      GAS3s=V2*FAC
      ISET=1
      else
      GAS3s=GSET
      ISET=0
      endIF
      return
      end
*     --------------------------------------------------------------------------
*     generate random samples from normal distribution for discharge quality ---
*     --------------------------------------------------------------------------
      function GAS4s(IDUM)
      save iset,gset
      data ISET /0/
      if ( IDUM .lt. 0 ) ISET=0
      if ( ISET .eq. 0 ) then
    1 continue     
      V1=2.*RAN4s(IDUM)-1.
      V2=2.*RAN4s(IDUM)-1.
      R=V1**2+V2**2
      if (R .ge. 1.) goto 1
      FAC=SQRoot(1013,-2.*LOG(R)/R)
      GSET=V1*FAC
      GAS4s=V2*FAC
      ISET=1
      else
      GAS4s=GSET
      ISET=0
      endIF
      return
      end


      function GAS8s(IDUM)
      save iset,gset
      data ISET /0/
      if ( IDUM .lt. 0 ) ISET=0
      if ( ISET .eq. 0 ) then
    1 continue     
      V1=2.*RAN8s(IDUM)-1.0
      V2=2.*RAN8s(IDUM)-1.0
      R=V1**2+V2**2
      if (R .ge. 1.0) goto 1
      FAC=SQRoot(1014,-2.0*LOG(R)/R)
      GSET=V1*FAC
      GAS8s=V2*FAC
      ISET=1
      else
      GAS8s=GSET
      ISET=0
      endIF
      return
      end
      
      function GAS9s(IDUM)
      save iset,gset
      data ISET /0/
      if ( IDUM .lt. 0 ) ISET=0
      if ( ISET .eq. 0 ) then
    1 continue     
      V1=2.*RAN9s(IDUM)-1.0
      V2=2.*RAN9s(IDUM)-1.0
      R=V1**2+V2**2
      if (R .ge. 1.0) goto 1
      FAC=SQRoot(1014,-2.0*LOG(R)/R)
      GSET=V1*FAC
      GAS9s=V2*FAC
      ISET=1
      else
      GAS9s=GSET
      ISET=0
      endIF
      return
      end

*     --------------------------------------------------------------------------
*     Random number generators -------------------------------------------------  
*     --------------------------------------------------------------------------
*     generate random numbers for river flow -----------------------------------
*     --------------------------------------------------------------------------
      function RAN2s(IDUM)
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
      write(3,*)'Stopped in random number generation 2s',JJ
      stop
      endif
      RAN2s=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      return
      end
*     --------------------------------------------------------------------------
*     generate random numbers for discharge flow -------------------------------
*     --------------------------------------------------------------------------
      function RAN3s(IDUM)
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
      write(3,*)'Stopped in random number generation 3s',JJ
      stop
      endif
      RAN3s=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      return
      end
*     --------------------------------------------------------------------------
*     generate random numbers for discharge quality ----------------------------
*     --------------------------------------------------------------------------
      function RAN4s(IDUM)
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
      write(3,*)'Stopped in random number generation 4s',JJ
      endif
      RAN4s=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      return
      end


      function RAN8s(IDUM)
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
      write(3,*)'Stopped in random number generation 8s',JJ
      endif
      RAN8s=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      return
      end

      function RAN9s(IDUM)
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
      write(3,*)'Stopped in random number generation 8s',JJ
      stop
      endif
      RAN9s=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      return
      end

      
