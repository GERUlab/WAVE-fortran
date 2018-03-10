c     version 3.0
c	author: m. vanclooster 15/02/2002

c###################################################################################
      subroutine sol_print
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'gen.com'
      include  'sol.com'

      if (printing_time(nday))call sol_sum
      if (pr_sol_cum) call sol_cum
      call sol_tseries
      end



c###################################################################################
      subroutine sol_sum
c     in   : cberr_sol, cberr_sol_isd, conirsol, decsoli, decsolm, dleasa,
c            dleasa_isd, dsol, dsol_isd, dt, fsol, idinof, idirr, idplow,
c            imobsw, isd, ncs, ncsplow, ninof, nirr, nplo, nr_of_sol, nts,
c            out_sol, pleasa, pleasa_isd, pr_sol_sum, ressol, rleasa,
c            rleasa_isd, sflxsa, simnit, t, tcsink, tcsink_isd, tcsolio,
c            tcsolmo, tcsolo, tflsol, tsoli, tsoli_isd, x
c     out  : -
c     calls: date_wave, date_title, heading, open_file, small_title
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      integer sum_file(kt_sol),sp,zero
      logical first_time
      integer year,month,day
      character *(15) n_sol_name(3)
      data n_sol_name /'ureum','ammonium','nitrate'/
      data first_time /.true./
      data iteliof,itelplow,itelirr/1,1,1/

      if (.not.pr_sol_sum) return
      do sp = 1, nr_of_sol
		if (out_sol(sp)) then
c			output for solute species sp
			if (first_time) then
				zero = ichar('0')
c				initialisation of file with overview of main state variables
				call open_file(sum_file(sp),'sol_sum'//char(sp+zero)// 
     $			'.out')
				call heading (sum_file(sp),'sol_sum'//char(sp+zero)//
     $			'.out' ,
     $			'summary of main state variables for the solute')
				if (simnit.and. sp.le.3) then
					write(sum_file(sp),1010) n_sol_name(sp)
1010					format (/,'solute : ',a,/,
     $					  '++++++++++++++++',//)
				else
					write(sum_file(sp),1015) sp
1015					format (/,'solute : ',i3,/,
     $					  '++++++++++++++++',//)
				endif
			endif
c			at the end of the day t = t_begin_day +1
			call date_title (sum_file(sp),idint(t)-1)
			call small_title(sum_file(sp),'time step')
			write(sum_file(sp),1020) dt,nts
1020			format(
     $		'time step                                  : ',e10.3,' day',/,
     $		'number of time step                        : ',i10,/)

c			upper boundary						
			call small_title(sum_file(sp),
     $		'terms of the upper boundary condition (mg/m**2)')
			write (sum_file(sp), 1070) ressol(sp), tflsol(sp)
1070			format(
     $		'amount of solute in the reservoir          : ' ,f10.1,/,
     $		'cumulative solute influx                   : ' ,f10.1,/)

c			irrigation
			if (nirr.gt.0) then
1080				if(idint(t).gt.idirr(itelirr)) then
					itelirr = min0(itelirr+1,nirr+1)
					if (itelirr.le.nirr) goto 1080
				endif	  
				if (idint(t).gt.idirr(1)) then
					call date_wave(idirr(itelirr - 1), year,month,day)
					write(sum_file(sp),1100)year,month,day,
     $					conirsol(itelirr-1,sp)
				endif
			endif
1100			format (/,
     $		'irrigation',/,
     $		'- most recent irrigation was on            : ',i4,1x,i2,1x,i2,/,
     $		'- amount of solute in irrigation solution  : ',f10.1)

c			inorganic fertilizer
			if (ninof.gt.0) then
1105				if(idint(t).gt.idinof(iteliof)) then
					iteliof = min0(iteliof+1,ninof+1) 
					if (iteliof.le.ninof) goto 1105
				endif	 
				if (idint(t).gt.idinof(1)) then
					call date_wave(idinof(iteliof - 1), year,month,day)
					write(sum_file(sp),1110)year,month,day,
     $					fsol(iteliof-1,sp)
				endif
			endif
1110			format (/,
     $		'inorganic fertilisation',/,
     $		'- most recent fertilizer application on    : ',i4,1x,i2,1x,i2,/,
     $		'- amount of solute applied                 : ', f10.1)

c			plowing
			if(nplo.gt.0) then
1115				if(idint(t).gt.idplow(itelplow)) then
					itelplow = min0(itelplow+1,nplo+1)
					if (itelplow.le.nplo) goto 1115
				endif	  
				if (idint(t).gt.idplow(1)) then
					call date_wave(idplow(itelplow - 1), year,month,day)
					write(sum_file(sp),1120)year,month,day,
     $					ncsplow(itelplow-1)
				endif
			endif
1120			format (/,
     $		'plowing',/,
     $		'- most recent time of plowing              : ',i4,1x,i2,1x,i2,/,
     $		'- number of compartments affected          : ', i12,/,
     $		'- solute influx with plowing               : ', f10.1)

c			lower boundary
			call small_title(sum_file(sp),
     $			'terms of the lower boundary condition (mg/m**2)')
			write(sum_file(sp),1150)
1150			format(
     $		'                                             whole profile', 
     $		'    isd',/,
     $		'                                             -------------', 
     $		'  -------')
			write (sum_file(sp),1200) 
     $			pleasa(sp),pleasa_isd(sp),
     $            dleasa(sp),dleasa_isd(sp),
     $            rleasa(sp),rleasa_isd(sp)
1200			format(
     $		'integrated flux at the bottom                 ',/,
     $		'- positive (= upward) flux                  :  ',f9.1,2x,f9.1,/,
     $		'- negative (= downward) flux                :  ',f9.1,2x,f9.1,/,
     $		'                                             -------------', 
     $		'  -------',/
     $		,'                                              ',f10.1,1x,f10.1)

c			balance terms
			call small_title(sum_file(sp),
     $			'solute balance (mg/m**2)')
			write(sum_file(sp),1150)
			write (sum_file(sp),1300)
     $			tsoli(sp),tsoli_isd(sp),
     $			tsoli(sp)+dsol(sp),tsoli_isd(sp)+dsol_isd(sp),
     $			dsol(sp),dsol_isd(sp),
     $			tflsol(sp),tflsol(sp),
     $			plo_in(sp),plo_in(sp),
     $			tcsink(sp),tcsink_isd(sp),
     $			rleasa(sp),rleasa_isd(sp),
     $			cberr_sol(sp),cberr_sol_isd(sp)
1300			format(
     $		'initial profile total                       :  ',f9.1,2x,f9.1,/,
     $		'actual  profile total                       :  ',f9.1,2x,f9.1,/,
     $	    'change for the profile                      :  ',f9.1,2x,f9.1,/,
     $	    'solute inflow at the top                    :  ',f9.1,2x,f9.1,/,
     $		'solute influx with plowing                  :  ',f9.1,2x,f9.1,/,
     $		'sink term                                   :  ',f9.1,2x,f9.1,/,
     $	    'integrated flux at the bottom               :  ',f9.1,2x,f9.1,/,
     $	    'integrated error on the balance             :  ',f9.1,2x,f9.1,/)

c			profile terms			
			call small_title(sum_file(sp),'profile')
			if (imobsw) then
				write(sum_file(sp),1600) 
			else
				write(sum_file(sp),1650)
			endif
			do i=1,ncs
				if (imobsw) then
					write(sum_file(sp),1700) i,x(i),tcsolio(i,sp),
     $					tcsolmo(i,sp)
				else
					write(sum_file(sp),1750) i,x(i),tcsolo(i,sp)
				endif
			enddo 
1600			format(
     $		'comp    depth   tcsolio   tcsolmo  ',/,
     $		'        (mm)   (mg/m**2) (mg/m**2)')
1650			format(
     $		'comp  depth      tcsol ',/,
     $		'      (mm)   (mg/m**2)')
1700			format(i4,1x,f8.1,1x,2(f9.1,1x))
1750			format(i4,1x,f8.1,1x,1(f9.1,1x))
		endif
	enddo
      first_time = .false.
      return
      end

c###################################################################################
      subroutine sol_cum
c     in   : cberr_sol, cberr_sol_isd, dsol, dsol_isd, nr_of_sol, out_sol,
c            rleasa, rleasa_isd, simnit, t, tcsink, tcsink_isd, tflsol
c     out  : -
c     calls: date_wave, heading, open_file
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      integer cum_file(kt_sol)
      integer zero,sp
      integer year, month, day
      logical first_time
      character *(15) n_sol_name(3)
      data n_sol_name /'ureum-n','ammonium-n','nitrate-n'/
      data first_time /.true./

      do sp = 1, nr_of_sol
		if (out_sol(sp)) then
			if (first_time) then
				zero = ichar('0')
				call open_file(cum_file(sp), 'sol_cum'//char(sp+zero)//'.out')
				call heading (cum_file(sp),'sol_cum'//char(sp+zero)//'.out' ,
     $			'main cumulative terms of the solute balance')
				if (simnit.and. sp.le.3) then
					write(cum_file(sp),900) n_sol_name(sp)
				else
					write(cum_file(sp),910) sp
				endif
900				format (/,'solute : ',a,/,
     $			'++++++++++++++++',//)
910				format (/,'solute : ',i3,/,
     $			'++++++++++++++++',//)
				write(cum_file(sp),1010)
				write(cum_file(sp),1050)
				write(cum_file(sp),1060)
1010				format(
     $			' tflsol     = cumulative solute inflow (mg/m**2)',/,
     $			' rleasa     = cumulative netto flux at bottom of the profile
     $			(mg/m**2)', /,
     $			' sink       = cumulative sink (mg/m**2)',/,
     $			' dsol       = change in solute content for the entire profile
     $			(mg/m**2)',/,
     $			' cberr_sol  = error on the balance (mg/m**2)',/,
     $			' rleasa_isd = cumulative netto flux bottom of the profile ',
     $			' (mg/m**2)', /,
     $			' sink_isd   = cumulative sink (mg/m**2)',/,
     $			' dsol_isd   = change in solute content for the entire profile
     $			(mg/m**2)',/,
     $			' cberr_sol_isd= error on the balance (mg/m**2)',//)
1050				format(
     $			'year month day      tflsol    rleasa     sink        dsol  ',
     $			'cberr_sol  rleasa_isd  sink_isd  dsol_isd  cberr_sol_isd')
1060				format(
     $			'                  (mg/m**2) (mg/m**2)  (mg/m**2)  (mg/m**2)',
     $			'(mg/m**2)   (mg/m**2)  (mg/m**2) (mg/m**2)   (mg/m**2)')
			endif

			call date_wave (idint(t)-1, year, month,day)
			write(cum_file(sp),1100) year, month,day, tflsol(sp), rleasa(sp),
     $		tcsink(sp), dsol(sp), cberr_sol(sp), rleasa_isd(sp),  
     $		tcsink_isd(sp), dsol_isd(sp), cberr_sol_isd(sp)
1100			format(i4,2x,i2,3x,i2,1x,9(f10.1,1x))
		endif
10    enddo
      first_time = .false.
      return
      end

c###################################################################################
      subroutine sol_tseries
c     in   : flxs, from, ncs, nday, nr_of_points, nr_of_sol, out_sol, 
c            pr_sol_tseries, sflxsa, simnit, tcsolo, tinit, to, tsinki,
c            tsinkm
c     out  : -
c     calls: heading, heading_table, open_file, sum, wr_table
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      integer solc_file(kt_sol),sink_file(kt_sol)
      integer zero,sp,time,l_name(3)
      dimension out(kt_comps+1),solout(kt_comps),sinkout(kt_comps)
      logical first_time
      character *(15) n_sol_name(3)
      data n_sol_name /'ureu','nh4','no3'/
      data l_name /4,3,3/
      data first_time /.true./ 

      time = nday +idint(tinit) -1
      do sp = 1, nr_of_sol
		if (out_sol(sp)) then

c			file with values for solute concentration in selected compartments
			if (pr_sol_tseries(1)) then
			    if (first_time) then
					if (simnit.and.sp.lt.4) then
						call open_file(solc_file(sp),n_sol_name(sp)
     $					(1:l_name(sp))//'conc.out')
						call heading (solc_file(sp),n_sol_name(sp)
     $					(1:l_name(sp))//'conc.out',
     $					'solute concentration (mg/m**2) for solute: '
     $					//n_sol_name(sp)(1:l_name(sp))//'-n')
						call heading_table(solc_file(sp),from,to,
     $					nr_of_points,9)
					else
						zero = ichar('0')
						call open_file(solc_file(sp),
     $					'conc'//char(sp+zero)//'.out')
						call heading (solc_file(sp),
     $					'conc'//char(sp+zero)//'.out',
     $					'solute concentration (mg/m**2) for solute: '
     $					//char(sp+zero))
						call heading_table(solc_file(sp),from,to,
     $					nr_of_points,9)
					endif
				endif
				do i = 1, ncs
					solout(i) = tcsolo(i,sp)
				enddo
				call sum(solout,out,1.d0,from,to,nr_of_points)
				call wr_table (solc_file(sp),out,time,nr_of_points, 9,2)
			endif

c			file with values for solute sink term in selected compartments
			if (pr_sol_tseries(3)) then
				if (first_time) then
					if (simnit.and.sp.lt.4) then
						call open_file(sink_file(sp),n_sol_name(sp)
     $					(1:l_name(sp))//'sink.out')
						call heading (sink_file(sp),n_sol_name(sp)
     $					(1:l_name(sp))//'sink.out',
     $					'solute sink (mg/m**2) for solute : '// 
     $					n_sol_name(sp)(1:l_name(sp)))
						call heading_table(sink_file(sp),from,to,
     $					nr_of_points,9)
					else
						zero = ichar('0')
						call open_file(sink_file(sp),
     $					'sink'//char(sp+zero)//'.out')
						call heading (sink_file(sp),
     $					'sink'//char(sp+zero)//'.out',
     $					'solute sink (mg/m**2)for solute : '//char(sp+zero))
						call heading_table(sink_file(sp),from,to,
     $					nr_of_points,9)
					endif
				endif
				do i = 1, ncs
				    sinkout(i) = tsinki(i,sp) + tsinkm(i,sp)
				enddo
				call sum(sinkout,out,1.d0,from,to,nr_of_points)
				call wr_table (sink_file(sp),out,time,nr_of_points, 9,2)
			endif
		endif
      enddo
      first_time = .false.
      return
      end