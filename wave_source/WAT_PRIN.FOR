c	version 3.0
c	author:m vanclooster 19/02/02
c###################################################################################
      subroutine wat_print
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'gen.com'
      include  'wat.com'

      if (printing_time(nday)) call wat_sum
      if (pr_wat_cum) call wat_cum
      if (pr_wat_et)  call wat_et
      call wat_tseries
      end

c###################################################################################
      subroutine wat_sum
c     in   : cberr_wat, cberr_wat_isd, cev_intc, cflbu, cflbup, cflsd, 
c            cflsdp, cinf, cirr, conduc, cprec, cpsev, cptra, csev, ctra, 
c            cum_intc, cwatupt_isd, drz, dt, flxs, gwl, ibboco, isd, 
c            itertot, ncomp, ncs, nday, nrepeat, nr_no_conv, nts, ph, 
c            pr_wat_sum, t, visd, visdi, vol, voli, wc, x
c     out  : -
c     calls: date_title, heading, open_file, small_title
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      integer sum_file,day_from,month_from,year_from
      integer day_to,month_to,year_to
      logical first_time
      data first_time /.true./
      data i_pond_old /1/

      if (.not.pr_wat_sum) return

c	initialisation of file with overview of main state variables
      if (first_time) then
		first_time = .false.
		call open_file(sum_file, 'wat_sum.out')
		call heading(sum_file, 'wat_sum.out',
     $		'summary of the main state variables for water')
      endif

c     output at the end of the day
      call date_title (sum_file,idint(t)-1)
      call small_title(sum_file,'time step and iterations')
      ratio=dfloat(itertot)/dfloat(nts)
      write(sum_file,10) dt,nts,itertot,ratio,nrepeat,nr_no_conv 
10    format(
     $' time step                                 : ',e10.3,' day',/,
     $' number of time step                       : ',i10,/,
     $' total number of newton-raphson iterations : ',i10,/,
     $' average number of iterations per time step: ',f10.3,/,
     $' number of times the time step was halved  : ',i10,/,
     $' number of times there was no convergence  : ',i10)
      call small_title(sum_file,'terms of the upper boundary condition (mm)')

      if (pond_max.gt.0.d0) then
		write (sum_file,20) pond
		if (i_pond.gt.0) write (sum_file, 30)
		do i = i_pond_old, i_pond
			t_day_from  = idint(pond_from(i))
			t_hour_from = dfloat(idint((pond_from(i)-t_day_from)*24.d0))
			t_min_from  =
     $		(((pond_from(i)-t_day_from)*24.d0) -t_hour_from)*60.d0
			call date_wave(idint(t_day_from), year_from, month_from, day_from)
			if (i.eq.i_pond.and.ponded)then
				t_day_to = t
				t_hour_to = 0.d0
				t_min_to = 0.d0
				pond_maximum(i) = dmax1(pond_maximum(i),pond)
			else
				t_day_to  = dfloat(idint(pond_to(i)))
				t_hour_to = dfloat(idint((pond_to(i)-t_day_to)*24))
				t_min_to  = (((pond_to(i)-t_day_to)*24) -t_hour_to)*60
			endif
			call date_wave(idint(t_day_to), year_to, month_to, day_to)
			write (sum_file,40) year_from, month_from, day_from, 
     $		idint(t_hour_from), t_min_from, year_to, month_to, day_to,
     $		idint(t_hour_to), t_min_to, pond_maximum(i)
		enddo
		i_pond_old = i_pond+1
		if (ponded) then
			pond_maximum(i_pond+1)= pond_maximum(i_pond)
			i_pond = i_pond + 1
			pond_from(i_pond) = t
		endif
20		format(' ponding',/,
     $		' - ponding depth at the end of this day  : ',f6.2)
30		format(' - ponding history since last output',/,
     $        '      from                to             max depth',/,
     $        '  year/mn/dy/hr/min   year/mn/dy/hr/min      ')
40		format(2x,i4,'/',i2,'/',i2,'/',i2,'/',f5.2,1x,i4,'/',i2,'/',
     $         i2,'/',i2,'/',f5.2,1x, f6.2)
      endif

      cpetr=cptra+cpsev
      cetr=ctra+csev + cev_intc
      cpinf=cprec+cirr-cum_intc
      cruno=cpinf-cinf-pond
      write(sum_file,50) cprec,cruno,cirr,cpinf,cinf,cum_intc,
     $	cptra, ctra, cev_intc, cpsev,csev, cpetr,cetr
50    format(
     $' precipitation                    : ',f7.2,/,
     $' runoff                           : ',f7.2,/,
     $' applied irrigation water         : ',f7.2,/,
     $' potential infiltration           : ',f7.2,/,
     $' actual infiltration              : ',f7.2,/,
     $' interception                     : ',f7.2,/,
     $ /,
     $' potential transpiration          : ',f7.2,/,
     $' actual transpiration             : ',f7.2,/,
     $' evaporation of intercepted water : ',f7.2,/,
     $' potential evaporation            : ',f7.2,/,
     $' actual evaporation               : ',f7.2,/,
     $' potential evapotranspiration     : ',f7.2,/,
     $' actual evapotranspiration        : ',f7.2,/)
      call small_title(sum_file,
     $	'terms of the bottom boundary condition')
      write(sum_file, 60)
60    format(
     $'                                   ','whole profile',
     $ 4x,'isd',/
     $'                                   ','-------------',
     $ 2x,'-------')

      cflbun=cflbu-cflbup
      cflsdn=cflsd-cflsdp
      write(sum_file,70) cflbup,cflsdp,cflbun,cflsdn,cflbu,cflsd
70    format(
     $' integrated flux at the bottom    ',/,
     $' - positive (= upward) flux       :  ',f9.1,2x,f9.1,/,
     $' - negative (= downward) flux     :  ',f9.1,2x,f9.1,/,
     $'                                   ','-------------', 2x,'-------',/,
     $'                                     ',f9.1,2x,f9.1,/)
      write(sum_file,80) drz(nday)
      if(ibboco.lt.4) write(sum_file,90) gwl(nday),ncomp
80    format(
     $	' actual rooting depth             : ',f10.1,' (mm)',/)
90    format(
     $	' actual groundwater level         : ',f10.1,' (mm)',/,
     $	' number of the first unsaturated node above ',/,
     $	'   the groundwater table          :       ',i5)

      call small_title(sum_file,'water balance ')
      write(sum_file, 60)
      write(sum_file,100)
     $	voli,visdi,
     $	vol,visd,
     $	vol - voli,visd-visdi,
     $	cinf-csev,cinf-csev,
     $	-ctra,-cwatupt_isd,
     $	cflbu,cflsd,
     $	cberr_wat,cberr_wat_isd
100   format(
     $	' water storage at the initial time:  ',f9.1,2x,f9.1,/,
     $	' actual water storage             :  ',f9.1,2x,f9.1,/,
     $	' change in water storage          :  ',f9.1,2x,f9.1,/,
     $	' integrated flux at the top       :  ',f9.1,2x,f9.1,/,
     $	' integrated plant water uptake    :  ',f9.1,2x,f9.1,/,
     $	' integrated flux at the bottom    :  ',f9.1,2x,f9.1,/,
     $	' integrated error                 :  ',f9.1,2x,f9.1)
      call small_title(sum_file, 'profile')
      write(sum_file,110) 
110   format( 
     $'comp.  depth    theta      pr.head      conduc   ',/,
     $'       (mm)   (m**3/m**3)    (cm)      (mm/day)  ' ,/)
      do i=1,ncs
		write(sum_file,120) i,x(i),wc(i),ph(i)/10.,conduc(i)
120		format(1x,i4,1x,f8.1,2x,f5.3,1x,e13.3,1x,e13.3)
      enddo
      end

c###################################################################################
      subroutine wat_cum
c     in   : cberr_wat, cberr_wat_isd, cflbu, cflbup, cflsd, cflsdp, cinf,
c            cirr, cprec, csev, ctra, cum_intc, isd, t, visd, visdi, vol,
c            voli
c     out  : -
c     calls: date, heading, open_file
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      integer cum_file
      integer year, month,day
      logical first_time
      data first_time /.true./

c     initialisation
      if (first_time) then
		first_time = .false.
		call open_file(cum_file, 'wat_cum.out')
		call heading(cum_file, 'wat_cum.out',
     $		'cumulative terms of the water balance')
		write(cum_file,10)
		write(cum_file,20)
		write(cum_file,30)
10		format( /,
     $		' cprec         = cumulative precipitation (mm)',/,
     $		' cruno         = cumulative runoff (mm)',/,
     $		' cirr          = cumulative amount of applied irrigation ',
     $		'water',
     $		' (mm)',/,
     $		' cum_intc      = cumulative interception (mm)',/,
     $		' ctra          = cumulative actual transpiration (mm)',/,
     $		' csev          = cumulative actual evaporation (mm)',/,
     $		' cflbp         = integrated positive (= upward) flux at the ',
     $		'bottom of the soil profile (mm)',/,
     $		' cflbn         = integrated negative (= downward) flux at the ',
     $		'bottom of the soil profile (mm)',/,
     $		' dif           = change of water storage in the ',
     $		'entire profile ',
     $		'(mm)',/,
     $		' cberr_wat     = error on the balance (mm)',/,/,
     $		' cflsdp        = integrated positive (= upward) flux at the ',
     $		'bottom of the isd profile (mm)',/,
     $		' cflsdn        = integrated negative (= downward) ',
     $		'flux at the ',
     $		'bottom of the isd profile (mm)',/,
     $		' dif_isd       = change of water storage in the ',
     $		'isd profile (mm)',/,
     $		' cberr_wat_isd = error on the balance in the ',
     $		'isd profile (mm)',//)
20		format(
     $		'year month day  cprec  cruno  cirr  cum_intc ctra '
     $  ,		' csev  cflxbp  cflxbn voldif cberr_wat cflsdp cflsdn'
     $		,' dif_isd cberr_wat_isd ')
30		format(
     $		'                 (mm)   (mm)   (mm)   (mm)   (mm)'
     $		,'  (mm)   (mm)   (mm)    (mm)    (mm)    (mm)    '
     $		,'(mm)    (mm)    (mm)')
      endif

      cpinf = cprec+cirr-cum_intc
      cruno = cpinf-cinf-pond
      cflbun = cflbu-cflbup
      dif = vol-voli
      cflsdn = cflsd - cflsdp
      dif_isd = visd - visdi
      
	call date_wave(idint(t)-1, year, month,day)
      write(cum_file,40) year, month, day,cprec,cruno,cirr,cum_intc,
     $ ctra,csev,cflbup,cflbun,dif,cberr_wat, cflsdp, cflsdn, dif_isd, 
     $ cberr_wat_isd
40    format(i4,2x,i2,3x, i2, 1x,8(f10.1,1x),1x,6(f10.1,2x)  )
      return
      end

c###################################################################################
      subroutine wat_et
c     in   : akc, cpsev, cptra, csev, ctra, gwl, ibboco, nday, rlai, t
c     out  : -
c     calls: date, heading, open_file
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      integer et_file
      integer year, month,day
      logical first_time
      data first_time/.true./ 
      if (first_time) then
		first_time = .false.
		call open_file(et_file,'wat_et.out')
		if(ibboco.le.3) then
			call heading(et_file,'wat_et.out',
     $		'evapotransp. rates and depth of the ground water table')
			write(et_file,10)
10			format(/,
     $		'ptra  : potential transpiration rate         (mm/day)',/,
     $		'tra   : actual transpiration rate            (mm/day)',/,
     $		'psev  : potential evaporation rate           (mm/day)',/,
     $		'sev   : actual evaporation rate              (mm/day)',/,
     $		'cptra : cumulative potential transpiration   (mm)',/,
     $		'ctra  : cumulative actual transpiration      (mm)',/,
     $		'cpsev : cumulative potential evaporation     (mm)',/,
     $		'csev  : cumulative actual evaporation        (mm)',/,
     $		'lai   : leaf area index                      (m**2/m**2)',/,
     $		'kc    : kc factor',/,
     $		'gwla  : groundwater level                    (mm)',//,
     $		'year month day      ptra       tra       psev      sev    '
     $		'  cptra      ctra    cpsev      csev       lai       kc     '
     $		'gwla',/,
     $		'                     (mm)      (mm)      (mm)      (mm)    ',
     $		'  (mm)      (mm)      (mm)      (mm)     (m2/m2)     (-)   ',
     $		' (mm)  ')
		else
			call heading(et_file,'wat_et.out',
     $		'evapotranspiration rates')
			write(et_file,20)
20			format(/,
     $		'ptra  : potential transpiration rate         (mm/day)',/,
     $		'tra   : actual transpiration rate            (mm/day)',/,
     $		'psev  : potential evaporation rate           (mm/day)',/,
     $		'sev   : actual evaporation rate              (mm/day)',/,
     $		'cptra : cumulative potential transpiration   (mm)',/,
     $		'ctra  : cumulative actual transpiration      (mm)',/,
     $		'cpsev : cumulative potential evaporation     (mm)',/,
     $		'csev  : cumulative actual evaporation        (mm)',/,
     $		'lai   : leaf area index                      (m**2/m**2)',/,
     $		'kc    : kc factor',/,
     $		'year month day      ptra       tra       psev      sev    '
     $		'  cptra      ctra    cpsev      csev       lai       kc     '
     $		,/,
     $		'                     (mm)      (mm)      (mm)      (mm)    ',
     $		'  (mm)      (mm)      (mm)      (mm)     (m2/m2)     (-)   ')
        endif
      endif

      x1 = cptra-cptra2
      x2 = ctra-ctra2
      x3 = cpsev-cpsev2
      x4 = csev-csev2
      cptra2 = cptra
      ctra2 = ctra
      cpsev2 = cpsev
      csev2 = csev
      call date_wave(idint(t)-1, year, month,day)
      if(ibboco.le.3) then
		write(et_file,30) year, month, day,x1,x2,x3,x4,cptra,
     $	ctra,cpsev,csev,rlai(nday),akc, gwl(nday)
      else
		write(et_file,30) year, month, day,x1,x2,x3,x4,cptra,
     $	ctra,cpsev,csev,rlai(nday), akc
      endif
30	format(i4,2x,i2,3x,i2,1x,12f10.5)
      end

c###################################################################################
      subroutine wat_tseries
c     in   : crtex, flxs, from, nday, nr_of_points, ph, pr_wat_tseries,
c            rtex, tinit, to, wc
c     out  : -
c     calls: avg, heading, heading_table, open_file, sum, wr_table
c###################################################################################
      implicit double precision (a-h,o-z)
      include 'constant'
      include 'gen.com'
      include 'wat.com'
      dimension out(kt_comps)
      integer ph_file,wc_file,crtex_file,time
      logical first_time
      data first_time /.true./

      time = nday +idint(tinit)-1

c     pressure head in selected compartments
      if (pr_wat_tseries(1)) then
		if (first_time) then
			call open_file(ph_file,'ph.out')
			call heading(ph_file,'ph.out',
     $		'average pressure head (cm)')
			call heading_table(ph_file,from,to,nr_of_points,8)
		endif
		call avg (ph,out,0.1d0,1,from,to,nr_of_points)
		call wr_table (ph_file,out,time,nr_of_points,8,1)
      endif

c     water contents in selected compartments
      if (pr_wat_tseries(2)) then
		if (first_time) then
			call open_file(wc_file,'wc.out')
			call heading(wc_file,'wc.out',
     $		'average moisture content (cm**3/cm**3)')
			call heading_table(wc_file,from,to,nr_of_points,6)
		endif
		call avg (wc,out,1.0d0,1,from,to,nr_of_points)
		call wr_table (wc_file,out,time,nr_of_points,6,4)
      endif

c     cumulative root extraction rates in selected compartments
      if (pr_wat_tseries(5)) then
		if (first_time) then
			call open_file(crtex_file,'crtex.out')
			call heading(crtex_file,'crtex.out',
     $		'total cumulative root extraction (mm)')
			call heading_table(crtex_file,from,to,nr_of_points,7)
		endif
		call sum (crtex,out,1.0d0,from,to,nr_of_points)
		call wr_table (crtex_file,out,time,nr_of_points, 7,2)
      endif
      first_time = .false.
      return
      end