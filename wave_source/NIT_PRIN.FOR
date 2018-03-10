c	version 3.0
c	author:. vanclooster 24/02/2002
c###################################################################################
      subroutine nit_print
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'gen.com'
      include  'nit.com'

      if (printing_time(nday))call nit_sum
      if (pr_nit_cum) call nit_cum
      call nit_tseries
      end

c###################################################################################
      subroutine nit_sum
c     in   : carborg, cberr_corg, cberr_corg_isd, cberr_norg, cberr_norg_isd,
c            cchumo, cclito, ccmano, cnhumo, cnlito, cnmano, dcorg,
c            dcorg_isd, dnorg, dnorg_isd, dt, idman, isd, ncs, nof, nts,
c            pr_nit_sum, ptcarborg, ptcarborg_isd, ptchum, ptchum_isd,
c            ptclit, ptclit_isd, ptcman, ptcman_isd, ptden, ptden_isd,
c            pthyd, pthyd_isd, ptmin, ptmin_isd, ptnhum, ptnhum_isd, ptnit,
c            ptnitorg, ptnitorg_isd, ptnit_isd, ptnlit, ptnlit_isd, ptnman,
c            ptnman_isd, ptscorg, ptscorg_isd, ptsnorg, ptsnorg_isd, ptup,
c            ptup_isd, ptvol, ptvol_isd, rnitorg, t, tcorgs, tcorgs_isd,
c            tflcorg, tflnorg, tnorgs, tnorgs_isd, x
c     out  : -
c     calls: date, date_title, denitrification, heading, hydrolysis, nitrification,
c            open_file, small_title, volatilisation
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'nit.com'
      integer sum_file
      integer year,month,day
      logical first_time
      data first_time /.true./
      data itelof/1/

      if (.not.pr_nit_sum) return
      if (first_time) then
		first_time = .false.
		call open_file(sum_file, 'nit_sum.out')
		call heading(sum_file, 'nit_sum.out',
     $		'summary of the main state variables for the nitrogen balance')
      endif
      call date_title (sum_file, idint(t)-1)
      call small_title(sum_file,'time step')
      write(sum_file,20) dt,nts
      call small_title
     $(sum_file,'terms of the upper boundary condition (mg/m**2)')
      write (sum_file, 30) tflcorg, tflnorg
      if (nof.ne.0) then
10		if(idint(t).gt.idman(itelof)) then
			itelof = min0(itelof+1,nof+1)
			if (itelof.le.nof) goto 10
		endif
          if (idint(t).gt.idman(1)) then
			call date_wave(idman(itelof - 1), year,month,day)
			write(sum_file,40) year,month,day,
     $			carborg(max(1,itelof-1)), rnitorg(max(1,itelof-1))
		endif
      endif
      call small_title(sum_file,'totals of the organic matter pool (mg/m**2)')
      write (sum_file,50)
      write (sum_file,60)
     $	ptclit,ptnlit,ptclit_isd,ptnlit_isd,
     $    ptcman,ptnman,ptcman_isd,ptnman_isd,
     $    ptchum,ptnhum,ptchum_isd,ptnhum_isd,
     $    ptcarborg,ptnitorg,ptcarborg_isd,ptnitorg_isd   
      call small_title(sum_file,'terms of the mineral nitrogen balance (mg/m2)')
      write(sum_file,70)
      write(sum_file,80)
     $	ptup(2),ptup_isd(2),
     $	ptup(3),ptup_isd(3),
     $	ptmin(2),ptmin_isd(2),
     $	ptmin(3),ptmin_isd(3),
     $	ptvol,ptvol_isd,
     $	ptnit,ptnit_isd,
     $	ptden,ptden_isd,
     $	pthyd,pthyd_isd
      call small_title(sum_file,'terms of the organic matter balance (mg/m**2)')
      write (sum_file,90) 
      write (sum_file,100) 
     $	ptcarborg, ptnitorg,ptcarborg_isd,ptnitorg_isd,
     $	tcorgs, tnorgs,tcorgs_isd,tnorgs_isd,
     $	dcorg, dnorg,dcorg_isd,dnorg_isd,
     $	tflcorg, tflnorg,tflcorg,tflnorg,
     $	ptscorg,ptsnorg,ptscorg_isd,ptsnorg_isd,
     $	cberr_corg, cberr_norg,cberr_corg_isd,cberr_norg_isd
      call small_title(sum_file,'organic matter profile') 
      write (sum_file,110)
      do i = 1, ncs
        write (sum_file, 120) i,x(i),cclito(i), cnlito(i), ccmano(i),
     $  cnmano(i), cchumo(i) , cnhumo(i)
	enddo 
20    format(
     $' time step                      : ',e10.3,' (day)',/,
     $' number of time step            : ',i10)
30    format(' total inflow ',/,
     $' - amount of carbon             : ',f10.1,/,
     $' - amount of nitrogen           : ',f10.1,/)
40    format(
     $' most recent manure application',/,
     $' - application on               : ',
     $i4,1x,i2,1x,i2/,
     $' - amount of carbon             :  ',f9.1 ,/,
     $' - amount of nitrogen           :  ',f9.1)
50    format(
     $'         xxxxxxwhole profilexxx ',
     $'   xxxxxxxxxxisdxxxxxxxxx',/,
     $'             carbon    nitrogen ',
     $'      carbon    nitrogen',/,
     $'         ----------  ---------- ',
     $'   ----------  ----------')
60    format('litter: ',f10.1,2x,f10.1,5x,f10.1,2x,f10.1,/,
     $       'manure: ',f10.1,2x,f10.1,5x,f10.1,2x,f10.1,/,
     $       'humus : ',f10.1,2x,f10.1,5x,f10.1,2x,f10.1,/,
     $       '        ----------  ----------',
     $'     ----------  ----------',/,
     $       'total : ',f10.1,2x,f10.1,5x,f10.1,2x,f10.1,/)

70    format(
     $'                                  ',
     $'whole profile      isd',/
     $'                                  ',
     $'-------------      ---')
80    format(
     $'uptake of ammonium-n            : ',f10.1,2x,f10.1,/,
     $'uptake of nitrate-n             : ',f10.1,2x,f10.1,/,
     $'amount of ammonium mineralised  : ',f10.1,2x,f10.1,/,
     $'amount of nitrate mineralised   : ',f10.1,2x,f10.1,/,
     $'amount of volatilisation        : ',f10.1,2x,f10.1,/,
     $'amount of nitrification         : ',f10.1,2x,f10.1,/,
     $'amount of denitrification       : ',f10.1,2x,f10.1,/,
     $'amount of hydrolysis            : ',f10.1,2x,f10.1,/)
90    format(
     $'         xxxxxxwhole profilexxx ',
     $'     xxxxxxxxxxisdxxxxxxxxx',/,
     $'           carbon    nitrogen ',
     $'         carbon    nitrogen',/,
     $'         ----------  ---------- ',
     $'     ----------  ----------')
100   format(
     $'actual : ', f10.1,2x,f10.1,5x,f10.1,2x,f10.1,/,
     $'initial: ', f10.1,2x,f10.1,5x,f10.1,2x,f10.1,/,
     $'change : ', f10.1,2x,f10.1,5x,f10.1,2x,f10.1,/,
     $'inflow : ', f10.1,2x,f10.1,5x,f10.1,2x,f10.1,/,
     $'sink   : ', f10.1,2x,f10.1,5x,f10.1,2x,f10.1,/,
     $'error  : ', f10.1,2x,f10.1,5x,f10.1,2x,f10.1,/)
110   format (
     $' comp depth   cclito  cnlito  ccmano  cnmano  cchumo',
     $'  cnhumo',/,
     $'      (mm)                  (mg/m**2)')
120   format (i3,1x,f8.1,1x,6(e7.1,1x))
      return
      end

c###################################################################################
      subroutine nit_cum
c     in   : cberr_corg, cberr_corg_isd, cberr_norg, cberr_norg_isd,
c            dcorg, dcorg_isd, dnorg, dnorg_isd, ptscorg, ptscorg_isd,
c            ptsnorg, ptsnorg_isd, t, tflcorg, tflnorg
c     out  : -
c     calls: date, heading, open_file
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'nit.com'
      integer cum_file
      logical first_time
      integer year, month, day
      data first_time /.true./

      if (first_time) then
		first_time = .false.
		call open_file(cum_file, 'nit_cum.out')
		call heading(cum_file, 'nit_cum.out',
     $		'cumulative terms of the organic pool balance ')
		write(cum_file,10)
		write(cum_file,20)
      endif
      call date_wave (idint(t)- 1, year, month,day)
      write(cum_file,30) year, month,day, tflnorg, tflcorg, 
     $ptsnorg, ptscorg, dnorg, dcorg, cberr_norg, cberr_corg,
     $ptsnorg_isd, ptscorg_isd, dnorg_isd, dcorg_isd, 
     $cberr_norg_isd, cberr_corg_isd
10    format(
     $	'tflnorg        = cumulative  organic nitrogen inflow (mg/m**2)'
     $	,/,
     $	'tflcorg        = cumulative  organic carbon inflow (mg/m**2)'
     $	,/,/,
     $	'terms for the whole profile:',/,
     $	'ptsnorg        = cumulative organic nitrogen sink (mg/m**2)',/,
     $	'ptscorg        = cumulative organic carbon sink (mg/m**2)',/,
     $	'dnorg          = change in organic nitrogen content ',
     $	'since start (mg/m**2)',/,
     $	'dcorg          = change in organic carbon content',
     $	' since start (mg/m**2)',/,
     $	'cberr_norg     = error on the balance for nitrogen (mg/m**2)'
     $	,/,
     $	'cberr_corg     = error on the balance for carbon (mg/m**2)',/
     $	,/,
     $	'terms for part of the profile between compartments 1 and isd:'
     $	,/,
     $	'ptsnorg_isd    = cumulative organic nitrogen sink (mg/m**2)'
     $	,/,
     $	'ptscorg_isd    = cumulative organic carbon sink (mg/m**2)',/,
     $	'dnorg_isd      = change in organic nitrogen content ',
     $	'since start (mg/m**2)',/,
     $	'dcorg_isd      = change in organic carbon content ',
     $	'since start (mg/m**2)',/,
     $	'cberr_norg_isd = error on the balance for nitrogen (mg/m**2)',
     $	/,
     $	'cberr_corg_isd = error on the balance for carbon (mg/m**2)'
     $	,//)
20    format('year month day (see table above for sequence: ',
     $	'top table is first column)')
30    format(i4,2x,i2,3x,i2,1x,14(1x,f8.0))
      return
      end

      subroutine nit_tseries
c###################################################################################
c     in   : cchumo, cclito, ccmano, cden, chyd, cmin, cnhumo, cnit, cnlito,
c            cnmano, cupt, cvol, from, ncs, nday, nr_of_points,
c            pr_nit_tseries, tinit, to
c     out  : -
c     calls: denitrification, heading, heading_table, hydrolysis, nitrification,
c            open_file, sum, volatilisation, wr_table
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'nit.com'
      dimension cminout(kt_comps),cuptout(kt_comps),out(kt_comps)
      integer nit_file(16)
      integer time
      logical first_time
      data first_time /.true./

      time = nday +idint(tinit)- 1
c     n in organic litter 
      if (pr_nit_tseries(1)) then
		if (first_time) then
			call open_file(nit_file(1),'orgnlit.out')
			call heading(nit_file(1),'orgnlit.out','litter nitrogen  (mg/m**2)')
			call heading_table(nit_file(1),from,to,nr_of_points,9)
		endif
		call sum (cnlito,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(1),out,time,nr_of_points, 9,2)
      endif
c     n in organic manure 
      if (pr_nit_tseries(2)) then
		if (first_time) then
			call open_file(nit_file(2),'orgnman.out')
			call heading(nit_file(2),'orgnman.out','manure nitrogen  (mg/m2)')
			call heading_table(nit_file(2),from,to,nr_of_points,9)
		endif
		call sum (cnmano,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(2),out,time,nr_of_points, 9,2)
      endif
c     n in organic humus 
      if (pr_nit_tseries(3)) then
		if (first_time) then
			call open_file(nit_file(3),'orgnhum.out')
			call heading(nit_file(3),'orgnhum.out','humus nitrogen  (mg/m2)')
			call heading_table(nit_file(3),from,to,nr_of_points,9)
		endif
		call sum (cnhumo,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(3),out,time,nr_of_points, 9,2)
      endif
c	c in organic litter
      if (pr_nit_tseries(4)) then
		if (first_time) then
			call open_file(nit_file(4),'orgclit.out')
			call heading(nit_file(4),'orgclit.out','litter carbon  (mg/m**2)')
			call heading_table(nit_file(4),from,to,nr_of_points,9)
		endif
		call sum (cclito,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(4),out,time,nr_of_points, 9,2)
      endif
c	c in organic manure
      if (pr_nit_tseries(5)) then
		if (first_time) then
			call open_file(nit_file(5),'orgcman.out')
			call heading(nit_file(5),'orgcman.out','manure carbon  (mg/m**2)')
			call heading_table(nit_file(5),from,to,nr_of_points,9)
		endif
		call sum (ccmano,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(5),out,time,nr_of_points, 9,2)
      endif
c	c in organic humus
      if (pr_nit_tseries(6)) then
		if (first_time) then
			call open_file(nit_file(6),'orgchum.out')
			call heading(nit_file(6),'orgchum.out','humus carbon  (mg/m**2)')
			call heading_table(nit_file(6),from,to,nr_of_points,9)
		endif
		call sum (cchumo,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(6),out,time,nr_of_points, 9,2)
      endif
c     nh4 uptake
      if (pr_nit_tseries(7)) then
		if (first_time) then
			call open_file(nit_file(7),'nh4upt.out')
			call heading(nit_file(7),'nh4upt.out',' nh4 uptake (mg/m2)')
			call heading_table(nit_file(7),from,to,nr_of_points,9)
		endif
		do i = 1, ncs
			cuptout(i) = cupt(i,2)
		enddo
		call sum (cuptout,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(7),out,time,nr_of_points, 9,2)
      endif
c	no3 uptake
      if (pr_nit_tseries(8)) then
		if (first_time) then
			call open_file(nit_file(8),'no3upt.out')
			call heading(nit_file(8),'no3upt.out', 'no3 uptake  (mg/m2)')
			call heading_table(nit_file(8),from,to,nr_of_points,9)
		endif
		do i = 1, ncs
			cuptout(i) = cupt(i,3)
		enddo
		call sum (cuptout,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(8),out,time, nr_of_points, 9,2)
      endif
c     denitrification 
      if (pr_nit_tseries(9)) then
		if (first_time) then
			call open_file(nit_file(9),'denitrif.out')
			call heading(nit_file(9),'denitrif.out','denitrification (mg/m2)')
			call heading_table(nit_file(9),from,to,nr_of_points,9)
		endif
		call sum (cden,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(9),out,time,nr_of_points, 9,2)
      endif
c	nitrification 
      if (pr_nit_tseries(10)) then
		if (first_time) then
			call open_file(nit_file(10),'nitrific.out')
			call heading(nit_file(10),'nitrific.out','nitrification  (mg/m2)')
			call heading_table(nit_file(10),from,to,nr_of_points,9)
		endif
		call sum (cnit,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(10),out,time,nr_of_points, 9,2)
      endif
c	ureum hydrolysis
      if (pr_nit_tseries(11)) then
		if (first_time) then
			call open_file(nit_file(11),'hydrolys.out')
			call heading(nit_file(11),'hydrolys.out','hydrolysis (mg/m2)')
			call heading_table(nit_file(11),from,to,nr_of_points,9)
		endif
		call sum (chyd,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(11),out,time,nr_of_points, 9,2)
      endif
c	volatilisation 
      if (pr_nit_tseries(12)) then
		if (first_time) then
			call open_file(nit_file(12),'volatil.out')
			call heading(nit_file(12),'volatil.out','volatilisation (mg/m2)')
			call heading_table(nit_file(12),from,to,nr_of_points,9)
		endif
		call sum (cvol,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(12),out,time,nr_of_points, 9,2)
      endif
c	mineralisation in nh4
      if (pr_nit_tseries(13)) then
		if (first_time) then
			call open_file(nit_file(13),'nh4min.out')
			call heading(nit_file(13),'nh4min.out','min. to nh4 (mg/m2)')
			call heading_table(nit_file(13),from,to,nr_of_points,9)
		endif
		do i = 1, ncs
			cminout(i) = cmin(i,2)
		enddo
		call sum (cminout,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(13),out,time,nr_of_points, 9,2)
      endif
c	mineralisation in nO3
      if (pr_nit_tseries(14)) then
		if (first_time) then
			call open_file(nit_file(14),'no3min.out')
			call heading(nit_file(14),'no3min.out','min. to no3 (mg/m2)')
			call heading_table(nit_file(14),from,to,nr_of_points,9)
		endif
		do i = 1, ncs
			cminout(i) = cmin(i,3)
		enddo
		call sum (cminout,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(14),out,time,nr_of_points, 9,2)
      endif
c	convective no3 uptake
      if (pr_nit_tseries(15)) then
		if (first_time) then
			call open_file(nit_file(15),'no3_cupt.out')
			call heading(nit_file(15),'no3_cupt.out',
     $			'conv. no3 uptake (mg/m2)')
			call heading_table(nit_file(15),from,to,nr_of_points,9)
		endif
		call sum (cunc_no3,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(15),out,time,nr_of_points, 9,4)
      endif
c	diffusive no3 uptake
      if (pr_nit_tseries(16)) then
		if (first_time) then
			call open_file(nit_file(16),'no3_dupt.out')
			call heading(nit_file(16),'no3_dupt.out',
     $			'diffus. no3 uptake (mg/m2)')
			call heading_table(nit_file(16),from,to,nr_of_points,9)
		endif
		call sum (cund_no3,out,1.0d0,from,to,nr_of_points)
		call wr_table (nit_file(16),out,time,nr_of_points, 9,4)
      endif
      first_time = .false.
      end