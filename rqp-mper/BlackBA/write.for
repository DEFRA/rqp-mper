      subroutine write downstream bioavailable metal (ic,iic)
      include 'bacom1cp.for'
      
      call sort out the errors (ic,iic) ! ######################################

      return
      end
      
     
      
      subroutine sort out the errors (ic,iic)
      include 'bacom1cp.for'
      
      if ( ic .eq. 4 .or. targit .lt. small ) then
      qualn = C(15,NS+1) ! number of samples for downstream quality ------------
      else
      qualn = C(14,NS+1) ! number of samples for downstream quality ------------
      endif

*     ##########################################################################
      if ( NERR .gt. 0 .and. ic .eq. 4 ) then ! sort upper and lower conf.limits 
      bph1 = bph1x ! mean bioavailable downstream of current discharge =========
      bph2 = bph2x ! mean bioavailable downstream of current discharge =========
      bca1 = bca1x ! mean bioavailable downstream of current discharge =========
      bca2 = bca2x ! mean bioavailable downstream of current discharge =========
      bed1 = bed1x ! mean bioavailable downstream of current discharge =========
      bed2 = bed2x ! mean bioavailable downstream of current discharge =========
      brd1 = brd1x ! mean bioavailable downstream of current discharge =========
      brd2 = brd2x ! mean bioavailable downstream of current discharge =========
      bgf1 = bgf1x ! mean bioavailable downstream of current discharge =========
      bgf2 = bgf2x ! mean bioavailable downstream of current discharge =========
          
      if ( bph1 .gt. bph2 ) then ! effect of pH on mean bioavailable d/s -------
      xxxx = bph1 ! swap the values --------------------------------------------
      bph1 = bph2
      bph2 = xxxx
      endif ! effect of pH
      if ( bca1 .gt. bca2 ) then ! effect of calcium
      xxxx = bca1 ! on mean bioavailable downstream
      bca1 = bca2
      bca2 = xxxx
      endif ! effect of calcium
      if ( bed1 .gt. bed2 ) then ! effect of discharge DOC
      xxxx = bed1 ! on mean bioavailable downstream
      bed1 = bed2
      bed2 = xxxx
      endif ! effect of discharge DOC
      if ( brd1 .gt. brd2 ) then ! effect of u/s river DOC
      xxxx = brd1 ! on mean bioavailable downstream
      brd1 = brd2
      brd2 = xxxx
      endif ! effect of u/s river DOC
      if ( bgf1 .gt. bgf2 ) then ! effect of PNEC
      xxxx = bgf1 ! on mean bioavailable downstream
      bgf1 = bgf2
      bgf2 = xxxx
      endif ! effect of PNEC
      
      if ( iic .ne. 98 ) then ! write the results wwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      write(3,276)tbm
  276 format(/146('=')/'Calculations for river biometal downstream ',
     &'of the current discharge'/146('=')/
     &'Effect on river quality of using confidence limits for pH, ',
     &'Calcium (Ca) and DOC (RD) ...',4x,
     &'... and the bio-equation (GF)'/146('=')/
     &'Resulting effect on the mean downstream concentration of ',
     &'biometal of ...',f9.2/146('='))
      write(3,376)bph1,bph2,bca1,bca2,bed1,bed2,brd1,brd2,bgf1,bgf2
  376 format('Confidence limits            ',
     &' pH: ',2f8.3,3x,' Ca: ',2f8.3,3x,
     &' ED: ',2f8.3,3x,' RD: ',2f8.3,3x,' GF: ',2f8.3)
      endif ! if ( iic .ne. 98 ) wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

      bph1z = (bph1 - tbmr) ! difference from master value 
      bph2z = (bph2 - tbmr) 
      bca1z = (bca1 - tbmr) 
      bca2z = (bca2 - tbmr) 
      bed1z = (bed1 - tbmr) 
      bed2z = (bed2 - tbmr) 
      brd1z = (brd1 - tbmr) 
      brd2z = (brd2 - tbmr) 
      bgf1z = (bgf1 - tbmr) 
      bgf2z = (bgf2 - tbmr) 
      
      if ( iic .ne. 98 ) then ! write the results wwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      write(3,476)bph1z,bph2z,bca1z,bca2z,bed1z,bed2z,brd1z,brd2z,
     &bgf1z,bgf2z
  476 format('Shift from mean of           ',
     &' pH: ',2f8.3,3x,' Ca: ',2f8.3,3x,
     &' ED: ',2f8.3,3x,' RD: ',2f8.3,3x,' GF: ',2f8.3)
      write(3,976)100.0*bph1z/tbmr,100.0*bph2z/tbmr,100.0*bca1z/tbmr,
     &100.0*bca2z/tbmr,100.0*bed1z/tbmr,100.0*bed2z/tbmr,
     &100.0*brd1z/tbmr,100.0*brd2z/tbmr,
     &100.0*bgf1z/tbmr,100.0*bgf2z/tbmr
  976 format('% change in mean             ',
     &' pH: ',2f8.3,3x,' Ca: ',2f8.3,3x,
     &' ED: ',2f8.3,3x,' RD: ',2f8.3,3x,' GF: ',2f8.3)
      endif ! if ( iic .ne. 98 ) wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      
      bph1z = (bph1z ) / tbmr
      bph2z = (bph2z ) / tbmr
      bca1z = (bca1z ) / tbmr
      bca2z = (bca2z ) / tbmr
      bed1z = (bed1z ) / tbmr
      bed2z = (bed2z ) / tbmr
      brd1z = (brd1z ) / tbmr
      brd2z = (brd2z ) / tbmr
      bgf1z = (bgf1z ) / tbmr
      bgf2z = (bgf2z ) / tbmr
      
      if ( iic .ne. 98 ) then ! write the results wwwwwwwwwwwwwwwwwwwwwwwwwwwwww
*     write(3,576)bph1z,bph2z,bca1z,bca2z,bed1z,bed2z,brd1z,brd2z,
*    &bgf1z,bgf2z
  576 format(
     &' pH: ',2f8.3,3x,' Ca: ',2f8.3,3x,
     &' ED: ',2f8.3,3x,' RD: ',2f8.3,3x,' GF: ',2f8.3,
     &' ... shift to mean ratio')
      endif ! if ( iic .ne. 98 ) wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

      t = errx ( 0.95 )
      bph1z = (bph1z / t) * tbmr
      bph2z = (bph2z / t) * tbmr
      bca1z = (bca1z / t) * tbmr
      bca2z = (bca2z / t) * tbmr
      bed1z = (bed1z / t) * tbmr
      bed2z = (bed2z / t) * tbmr
      brd1z = (brd1z / t) * tbmr
      brd2z = (brd2z / t) * tbmr
      bgf1z = (bgf1z / t) * tbmr
      bgf2z = (bgf2z / t) * tbmr
      
      if ( iic .ne. 98 ) then ! write the results wwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      write(3,676)bph1,bph2,bca1,bca2,bed1,bed2,brd1,brd2,bgf1,bgf2
  676 format('Standard errors              ',
     &' pH: ',2f8.3,3x,' Ca: ',2f8.3,3x,
     &' ED: ',2f8.3,3x,' RD: ',2f8.3,3x,' GF: ',2f8.3)
      endif ! if ( iic .ne. 98 ) wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww 
      
      bph1z = bph1z * bph1z
      bph2z = bph2z * bph2z
      bca1z = bca1z * bca1z
      bca2z = bca2z * bca2z
      bed1z = bed1z * bed1z
      bed2z = bed2z * bed2z
      brd1z = brd1z * brd1z
      brd2z = brd2z * brd2z
      bgf1z = bgf1z * bgf1z
      bgf2z = bgf2z * bgf2z
      
      if ( iic .ne. 98 ) then ! write out the results wwwwwwwwwwwwwwwwwwwwwwwwww
      write(3,674)bph1z,bph2z,bca1z,bca2z,bed1z,bed2z,brd1z,brd2z,
     &bgf1z,bgf2z
  674 format('Squared standard errors      ',
     &' pH: ',2f8.3,3x,' Ca: ',2f8.3,3x,
     &' ED: ',2f8.3,3x,' RD: ',2f8.3,3x,' GF: ',2f8.3/146('='))
      endif ! if ( iic .ne. 98 ) wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      
      SEM = tbs / SQRoot( 117333, qualn ) ! original standard error ------------
      zumt1 = SEM*SEM+bph1z+bca1z+bed1z+brd1z ! sum squares for pH etc ---------
      zumt2 = SEM*SEM+bph2z+bca2z+bed2z+brd2z ! sum squares for pH etc ---------
      zumt1 = SQRoot( 117333,zumt1) ! compounded standard error for pH etc -----
      zumt2 = SQRoot( 117333,zumt2) ! compounded standard error for pH etc -----
      zumt3 = amax1 (zumt1,zumt2) ! biggest compounded standard error for pH etc
      zbs4 = zumt3 * SQRoot( 117333, qualn ) ! st.dev with pH and Ca errors ----
      zcov3 = zbs4/tbmr ! biggest coefficient of variation for pH etc ----------

      sumt1 = SEM*SEM+bph1z+bca1z+bed1z+brd1z+bgf1z ! squares of standard errors
      sumt2 = SEM*SEM+bph2z+bca2z+bed2z+brd2z+bgf1z ! squares of standard errors
      sumt1 = SQRoot( 117333,sumt1) ! compounded standard error with equation --
      sumt2 = SQRoot( 117333,sumt2) ! compounded standard error with equation --
      sumt3 = amax1 (sumt1,sumt2) ! biggest compounded standard error ----------
      tbs4 = sumt3 * SQRoot( 117333, qualn ) ! st.dev with equation errors -----
      tcov3 = tbs4/tbmr !  biggest coefficient of variation with equation errors 
      
      if ( iic .ne. 98 ) then ! write the results wwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      write(3,4576)
 4576 format(/108('-')/'Calculations of bioavailable metal in the ',
     &'river downstream of the current discharge quality ...'/108('-'))
*     write(3,4433)
 4433 format('Introduce the errors from the bio-available ',
     &'equation ...'/108('-'))
*     write(3,4443)pfit2
*    &'Scaling factor for mean bio-available metal =',f10.4/108('-'))
      write(3,4422)
 4422 format(' [A] Based on the sampling rates for metal'/
     &' [B] as [A] plus the effects of sampling rates for pH, ',
     &'Calcium and DOC'/
     &' [C] as [B] plus the effects of errors in the biometric ',
     &'equation'/
     &108('-')/39x,'[A]',23x,'[B]',23x,'[C]'/108('-'))
      write(3,6426)sem,zumt3,sumt3
 6426 format('          Standard error = ',f15.2,11x,f15.2,11x,f15.2)
      write(3,4426)tbs,zbs4,tbs4
 4426 format('      Standard deviation = ',f15.2,11x,f15.2,11x,f15.2)
      write(3,4436)tbs/tbmr,zcov3,tcov3
 4436 format('Coefficient of variation = ',f15.2,11x,f15.2,11x,f15.2/
     &108('-'))
      endif ! if ( iic .ne. 98 ) wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      
      endif ! if ( NERR .gt. 0 .and. ic .eq. 4 )  sort out the confidence limits
*     ##########################################################################

*     ##########################################################################
      if ( NERR .gt. 0 ) then ! sort upper and lower conf.limits ===============
      if ( ic .eq. 3 .or. ic .eq. 5 .or. ic .eq. 6 ) then ! modified discharge = 
      if ( bph3 .gt. bph4 ) then ! effect of pH
      xxxx = bph1 ! on mean bioavailable downstream
      bph3 = bph4
      bph4 = xxxx
      endif ! effect of pH
      if ( bca3 .gt. bca4 ) then ! effect of calcium
      xxxx = bca3 ! on mean bioavailable downstream
      bca3 = bca4
      bca4 = xxxx
      endif ! effect of calcium
      if ( bed3 .gt. bed4 ) then ! effect of discharge DOC
      xxxx = bed3 ! on mean bioavailable downstream
      bed3 = bed4
      bed4 = xxxx
      endif ! effect of discharge DOC
      if ( brd3 .gt. brd4 ) then ! effect of u/s river DOC
      xxxx = brd3 ! on mean bioavailable downstream
      brd3 = brd4
      brd4 = xxxx
      endif ! effect of u/s river DOC
      if ( bgf3 .gt. bgf4 ) then ! effect of PNRC
      xxxx = bgf3 ! on mean bioavailable downstream
      bgf3 = bgf4
      bgf4 = xxxx
      endif ! effect of PNEC

      if ( iic .ne. 98 ) then ! write the results wwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      write(3,8276)bpm
 8276 format(/145('+')/'Calculations for river biometal downstream',
     &'of the modified discharge'/145('+')/
     &'Effect on river quality of using confidence limits for pH, ',
     &'Calcium (Ca) and DOC (RD) ...',26x,
     &'... and the bio-equation (GF)'/145('+')/
     &'Resulting effect on the mean downstream concentration of ',
     &'biometal of ...',f9.2/145('+'))
      write(3,8376)bph3,bph4,bca3,bca4,bed3,bed4,brd3,brd4,bgf3,bgf4
 8376 format('Confidence limits           ',
     &' pH: ',2f8.3,3x,' Ca: ',2f8.3,3x,
     &' ED: ',2f8.3,3x,' RD: ',2f8.3,3x,' GF: ',2f8.3)
      endif ! if ( iic .ne. 98 ) wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

      bph3x = (bph3 - tbm) ! difference from master value 
      bph4x = (bph4 - tbm) 
      bca3x = (bca3 - tbm) 
      bca4x = (bca4 - tbm) 
      bed3x = (bed3 - tbm) 
      bed4x = (bed4 - tbm) 
      brd3x = (brd3 - tbm) 
      brd4x = (brd4 - tbm) 
      bgf3x = (bgf3 - tbm) 
      bgf4x = (bgf4 - tbm) 
      
      if ( iic .ne. 98 ) then ! write the results wwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      write(3,8476)tbm,bph3x,bph4x,bca3x,bca4x,bed3x,bed4x,
     &brd3x,brd4x,bgf3x,bgf4x
 8476 format('Shift from mean of ',f9.2,
     &' pH: ',2f8.3,3x,' Ca: ',2f8.3,3x,
     &' ED: ',2f8.3,3x,' RD: ',2f8.3,3x,' GF: ',2f8.3)
      write(3,8976)100.0*bph3x/tbm,100.0*bph4x/tbm,100.0*bca3x/tbm,
     &100.0*bca4x/tbm,100.0*bed3x/tbm,100.0*bed4x/tbm,
     &100.0*brd3x/tbm,100.0*brd4x/tbm,100.0*bgf3x/tbm,100.0*bgf4x/tbm
 8976 format('% change in mean            ',
     &' pH: ',2f8.3,3x,' Ca: ',2f8.3,3x,
     &' ED: ',2f8.3,3x,' RD: ',2f8.3,3x,' GF: ',2f8.3)
      endif ! if ( iic .ne. 98 ) wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      
      bph3x = (bph3x) / tbm
      bph4x = (bph4x) / tbm
      bca3x = (bca3x) / tbm
      bca4x = (bca4x) / tbm
      bed3x = (bed3x) / tbm
      bed4x = (bed4x) / tbm
      brd3x = (brd3x) / tbm
      brd4x = (brd4x) / tbm
      bgf3x = (bgf3x) / tbm
      bgf4x = (bgf4x) / tbm
      
      if ( iic .ne. 98 ) then ! write the results wwwwwwwwwwwwwwwwwwwwwwwwwwwwww
*     write(3,8576)bph3x,bph4x,bca3x,bca4x,bed3x,bed4x,brd3x,brd4x,
*    &bgf3x,bgf4x
 8576 format(
     &' pH: ',2f8.3,3x,' Ca: ',2f8.3,3x,
     &' ED: ',2f8.3,3x,' RD: ',2f8.3,3x,' GF: ',2f8.3,
     &' ... shift to mean ratio')
      endif ! if ( iic .ne. 98 ) wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

      t = errx ( 0.95 )
      bph3x = (bph3x / t) * tbm
      bph4x = (bph4x / t) * tbm
      bca3x = (bca3x / t) * tbm
      bca4x = (bca4x / t) * tbm
      bed3x = (bed3x / t) * tbm
      bed4x = (bed4x / t) * tbm
      brd3x = (brd3x / t) * tbm
      brd4x = (brd4x / t) * tbm
      bgf3x = (bgf3x / t) * tbm
      bgf4x = (bgf4x / t) * tbm
      
      if ( iic .ne. 98 ) then ! write the results wwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      write(3,8676)bph3x,bph4x,bca3x,bca4x,bed3x,bed4x,brd3x,brd4x,
     &bgf3x,bgf4x
 8676 format('Standard errors             ',
     &' pH: ',2f8.3,3x,' Ca: ',2f8.3,3x,
     &' ED: ',2f8.3,3x,' RD: ',2f8.3,3x,' GF: ',2f8.3)
      endif ! if ( iic .ne. 98 ) wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

      bph3x = bph3x * bph3x
      bph4x = bph4x * bph4x
      bca3x = bca3x * bca3x
      bca4x = bca4x * bca4x
      bed3x = bed3x * bed3x
      bed4x = bed4x * bed4x
      brd3x = brd3x * brd3x
      brd4x = brd4x * brd4x
      bgf3x = bgf3x * bgf3x
      bgf4x = bgf4x * bgf4x
      
      if ( iic .ne. 98 ) then ! write the results wwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      write(3,8674)bph3x,bph4x,bca3x,bca4x,bed3x,bed4x,brd3x,brd4x,
     &bgf3x,bgf4x
 8674 format('Squared standard errors     ',
     &' pH: ',2f8.3,3x,' Ca: ',2f8.3,3x,
     &' ED: ',2f8.3,3x,' RD: ',2f8.3,3x,' GF: ',2f8.3/145('+'))
      endif ! if ( iic .ne. 98 ) wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      
      SEM = tbs / SQRoot( 117333, qualn )
      
      zumt1 = SEM*SEM+bph3x+bca3x+bed3x+brd3x ! sum squares of standard errors -
      zumt2 = SEM*SEM+bph4x+bca4x+bed4x+brd4x ! sum squares of standard errors -
      zumt1 = SQRoot( 117333,zumt1) ! calculate compounded standard error ------
      zumt2 = SQRoot( 117333,zumt2) ! calculate compounded standard error ------
      zumt3 = amax1 (zumt1,zumt2) ! biggest compounded standard error ----------
      zbs3 = zumt3 * SQRoot( 117333, qualn ) ! new st.dev d/s modified ---------
      zcov4 = zbs3/tbmr ! biggest coefficient of variation ---------------------

      sumt1 = SEM*SEM+bph3x+bca3x+bed3x+brd3x+bgf3x ! sum squares of errors ----
      sumt2 = SEM*SEM+bph4x+bca4x+bed4x+brd4x+bgf4x ! sum squares of errors ----
      sumt1 = SQRoot( 117333,sumt1) ! calculate compounded standard error ------
      sumt2 = SQRoot( 117333,sumt2) ! calculate compounded standard error ------
      sumt3 = amax1 (sumt1,sumt2) ! biggest compounded standard error ----------
      tbs3 = sumt3 * SQRoot( 117333, qualn ) ! new st.dev d/s modified ---------
      tcov4 = tbs3/tbmr ! biggest coefficient of variation ---------------------
      
      if ( iic .ne. 98 ) then ! write the results wwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      write(3,8586)
 8586 format(/108('+')/'Calculations of bioavailable metal in the ',
     &'river downstream of the modified discharge quality ...'/108('+'))
      write(3,8422)
 8422 format(' [A] Based on the sampling rates for metal'/
     &' [B] as [A] plus the effects of sampling rates for pH, ',
     &'Calcium and DOC'/
     &' [C] as [B] plus the effects of errors in the biometric '/
     &108('+')/39x,'[A]',23x,'[B]',23x,'[C]'/108('+'))
*     write(3,8442)pfit2
 8442 format('Average scaling factor =',f7.3/108('+'))
      write(3,8826)sem,zumt3,sumt3
 8826 format('          Standard error = ',f15.2,11x,f15.2,11x,f15.2)
      write(3,8426)tbs,zbs3,tbs3 ! d/s modified discharge =====================
 8426 format('      Standard deviation = ',f15.2,11x,f15.2,11x,f15.2)
      write(3,8436)tbs/tbmr,zcov4,tcov4
 8436 format('Coefficient of variation = ',f15.2,11x,f15.2,11x,f15.2/
     &108('+'))
      endif ! if ( iic .ne. 98 ) wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      
      endif ! if ( ic .eq. 3 .or. ic .eq. 5 .or. ic .eq. 6 ) modified ========== 
      endif ! if ( NERR .gt. 0 ) sort upper & lower conf.limits ================
*     ##########################################################################

      return
      end