c	version 3.0
c	author m. vanclooster 22/02/2002
c###################################################################################
      subroutine watdata
c     in   : dsp, dx, dxinter, isd, ncs, ndsim, simplant, x, yearin
c     out  : arel, arer, brel, brer, drz, 
c            dxn, eq_grwt, gprh, gwl, histor_file, hystlim,
c            ibbcct, ibboco, idvs, igiveph, ilinr, ismxnl, isucr, laict,
c            modelk, modelpf, nparkon, nparpf, npl, pond_max, p0, p1, p2h, p2l,
c            p3, parakon, parapf, ph, phsa, plant_date, pr_wat_cum,
c            pr_wat_et, pr_wat_histor, pr_wat_sum, pr_wat_tseries,
c            qdeep, rlai, rnam, rootct, smx, tb, te, visd,
c            visdi, vol, voli, wc
c     calls: calc_ph, calc_wc, check_et, date, day_nr, end_file, hold_line,
c            init_file, init_pfk, jdate, next_line, open_file, rd_date, rd_int,
c            rd_logic, rd_real, report_err, warn_err
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      character *100 message
      integer day,month
      integer year_for
      integer i, j
      double precision p1max, p2hmax,  p3max, p2lmax,rt_upt
      logical rd_ext, rd_logic,mo_than_1
      integer rd_int
      dimension test_min(kt_in_row-1),test_max(kt_in_row-1)
      dimension result (kt_in_row, 366)
      integer nr_par_pf(kt_modelpf), nr_par_kon(kt_modelkon) 
      data nr_par_pf/5,5,6,0/
      data nr_par_kon/3,2,2,2,2,0/
      call init_file('watdata.in')
c     water: parameters for each soil layer
      rd_ext =rd_logic ()
      if (rd_ext) call init_file('watpar.wp')
      multimodal = rd_logic()
      if (multimodal) then
		call init_n_p
		modelpf = 4
		modelk  = 6
      else
		modelpf= rd_int(1, kt_modelpf)
		modelk= rd_int(1, kt_modelkon)
		nparpf = nr_par_pf (modelpf)
		nparkon = nr_par_kon (modelk)
c		moisture retention parameters for each soil layer
		do i=1,npl
			call hold_line
			layer = rd_int(1, npl)
			if (layer.ne.i) then
				call report_err
     $				(' error in sequence of compartments')
			else 
c				residual water content
				parapf(i,1) = rd_real(0.d0, mx_real)
c				saturated water content > residual water content
				parapf(i,2) = rd_real(parapf(i,1)+mn_real, mx_real)
c				alpha (model 1)/ alphaw(model 2) / alphad(model 3)
				parapf(i,3) = rd_real(mn_real, mx_real)
c				n 
				parapf(i,4) = rd_real(mn_real, mx_real)
c				m
				parapf(i,5) = rd_real(mn_real, mx_real)
				if (modelpf.eq.3) then
c					alphaw(model 3)
					parapf(i,6) = rd_real (parapf(i,3),mx_real)
c					nw = nd
					parapf(i,7) = parapf(i,4)
c					mw = md
					parapf(i,8) = parapf(i,5)
				endif
			endif
			call next_line
		enddo
		call check_et
c		hysteresis model: minimum rel.press.head change at return point
		if (modelpf.eq.2.or.modelpf.eq.3) then
			hystlim = rd_real(0.00001d0,1.d0)
		endif
c		hydraulic conductivity parameters for each soil layer
		do i=1,npl
			call hold_line
			layer = rd_int(1, npl)
			if (layer.ne.i) then
				call report_err
     $				('error in sequence of compartments')
			else 
				do j=1,nparkon
					parakon(i,j) = rd_real(-mx_real, mx_real)
				enddo
			endif
			call next_line
		enddo
		call check_et
		if (modelk.eq.5) then 
			call warn_err
     $			('using m = 1 - 1/n for mrc for mualem model hcc')
			do i = 1,npl
				parapf(i,5) = (1.d0 - 1.d0/parapf(i,4))
				if (modelpf.eq.3) parapf(i,8) = parapf(i,5)
				parakon(i,3) = parapf(i,5)
			enddo
		endif
      endif
      if (rd_ext) call end_file 
c     water: upper boundary conditions
      rd_ext =rd_logic()
      if (rd_ext) call init_file('watubc.wu')
c     minimum allowed pressure head at the surface 
      phsa=- dabs(rd_real (-mx_ph, mx_ph))
c     ponding depth
      pond_max = rd_real (0.d0, mx_real)
      if (rd_ext) call end_file 

c     evaporation and transpiration
      rd_ext =rd_logic()
      if (rd_ext) call init_file('watcrop.wc')
      if (isucr) then
		idvs = rd_logic()
      else
		idvs = .false.
      endif
      nr_of_kc_values = rd_int (1,365)
      do i = 1, nr_of_kc_values
		call hold_line 
		if (idvs) then
			dvs_kc(i) = rd_real (0.d0, mx_real)
		else
			call rd_date (month,day)
			iday_kc(i) = jdate(year_for(month,day),month,day,yearin)
		endif
		if (i.gt.1) then
			if (idvs) then
				if (dvs_kc(i).lt.dvs_kc(i-1)) call report_err
     $				('dvs values for kc factors not ascending')
			else
				if (iday_kc(i).lt.iday_kc(i-1)) call report_err
     $				('dates for kc factors not ascending')
			endif   
		endif   
		kc(i) = rd_real (0.d0,mx_real)
		call next_line
      enddo
      call check_et

c	root extraction parameters 
      if (simplant) then
		call hold_line
		call rd_date (month, day)
		call next_line
		tb = dfloat(jdate(year_for(month,day),month,day,yearin))
		call hold_line
		call rd_date (month, day)
		call next_line
		te = dfloat(jdate(year_for(month,day),month,day,yearin))
		if (te.lt.tb) then 
			call report_err
     $		('roots start to become inactive after date at '//
     $		'which they reach maximum inactivity')
		endif
c		value of the rootlength below which the roots can no longer extract water
		rnam = - dabs(rd_real(-dsp,dsp))
c		lai and rooting depth
		if(.not.isucr) then
			test_min(1) = 0.d0
			test_max(1) = 15.d0
			call interp_tseries (result,1, 'lai',test_min, test_max,9)
			do i = 1, ndsim
				rlai(i) = result (1,i)
			enddo
			test_min(1) = -dsp
			test_max(1) = dsp
			call interp_tseries (result,1,'rooting depth',test_min,test_max,9)
			do i = 1,ndsim
				drz(i) = - dabs(result (1,i))
			enddo
		endif
c		water: sink term variables
		ilinr= rd_logic()
		
		p0 = rd_real (-mx_ph,mx_ph)
		p1max =  dabs(p0)
		p0  = -10.0d0 * dabs(p0)
		p1 = rd_real (-mx_ph,p1max)
		p2hmax =  dabs(p1)
		p1  = - 10.d0* dabs(p1)
		p1  =  dmin1(p1,(p0-0.00000000001d0))
		p2h = rd_real (-mx_ph,p2hmax)
		p3max =  dabs(p2h)
		p2lmax =  dabs(p2h)
		p2h = - 10.0d0 * dabs(p2h)
		p2l = rd_real (-mx_ph,p2lmax)
		p3max = dmax1( dabs(p2l), p3max)
		p2l = - 10.0d0 * dabs(p2l)
		p3 = rd_real (-mx_ph, p3max)
		p3  = - 10.d0* dabs(p3)
		p3  =  dmin1(p3,p2h-0.00000000001d0,p2l-0.00000000001d0)
		
		ismxnl = rd_logic()
		if (ismxnl) then
			arer = dabs(rd_real (0.0d0, mx_real))
			brer = - dabs(rd_real (-mx_real, mx_real))
			i = 1
10			if (i .le.ncs) then 
				rt_upt = arer+brer* dabs(x(i))
				if (rt_upt.ge.0.0d0) then
					rt_distr(i) = rt_upt
					i = i+1
					goto 10
				else
					write (message, 20) i,ncs
20					format('root water uptake was negative from compartment ',
     $				i3 , ' to ', i3	,' ; replaced bij 0')
					call warn_err(message)
					do j = i, ncs
						rt_distr(j)= 0.0d0
					enddo 
				endif
			endif
		else
			test_min(1) = 0.d0
			test_max(1) = mx_real
c			interp_profile (result, nr_rows, nr_cols, table_name
c                         ,test_min, test_max,nine)
			call interp_profile (result, ncs, 1, 'rt_distr(z)'
     $                         ,test_min, test_max,9)
			do i = 1,ncs
				rt_distr(i) = result (1,i)
			enddo
		endif
      endif
      if (rd_ext) call end_file 

c     water: bottom boundary condition
      rd_ext =rd_logic()
      if (rd_ext) call init_file('watbbc.wb')
      ibboco = 0
      mo_than_1 = .false.
c     bottom boundary condition #1: groundwater level 
      if (rd_logic()) then
		ibboco = 1 
		test_min(1) = -dsp
		test_max(1) = dsp
		call interp_tseries (result,1, 'gwl(t)',test_min, test_max,9)
		do i = 1,ndsim
			gwl(i) = - dabs(result (1,i))
		enddo
      endif
c     bottom boundary condition #2: flux from the saturated zone
      if (rd_logic()) then
		if (ibboco.eq.0) then
			ibboco = 2
		else
			mo_than_1 = .true.
		endif
		gwl(1) = -abs(rd_real (- dsp,dsp))
		test_min(1) = -mx_real
		test_max(1) = mx_real
		call interp_tseries (result,1, 'qdeep(t)',test_min, test_max,9)
		do  i = 1,ndsim
			qdeep(i) = result (1,i)
		enddo 
      endif
c     bottom boundary condition #3: flux from the saturated zone as a func. of gwl
      if (rd_logic()) then
		if (ibboco.eq.0) then
			ibboco = 3
		else
			mo_than_1 = .true.
		endif
		gwl(1) = - dabs(rd_real (- dsp,dsp))
		arel = rd_real (-mx_real, mx_real)
		brel = rd_real (-mx_real, mx_real)
      endif
c     bottom boundary condition #4: pressure head 
      if (rd_logic()) then
		if (ibboco.eq.0) then
			ibboco = 4 
		else
			mo_than_1 = .true.
		endif
		test_min(1) = -mx_ph
		test_max(1) = mx_ph
		call interp_tseries (result,1, 'gprh(t)'
     $		,test_min, test_max,9)
		do i = 1, ndsim
			gprh(i) = -10.d0 * dabs(result (1,i))
		enddo
      endif
c     bottom boundary condition #5: flux known
      if (rd_logic()) then
		if (ibboco.eq.0) then
			ibboco = 5 
		else
			mo_than_1 = .true.
		endif
      endif 
c     bottom boundary condition #6: free drainage 
      if (rd_logic()) then
		if (ibboco.eq.0) then
			ibboco = 6 
		else
			mo_than_1 = .true.
		endif
      endif 
c     bottom boundary condition #7: lysimeter condition
      if (rd_logic()) then
		if (ibboco.eq.0) then
			ibboco = 7 
		else
			mo_than_1 = .true.
		endif
      endif
	
c     check input bottom boundary condition
      if (ibboco.eq.0) then
		call report_err ('no bottom boundary was chosen')
      end if
      if (mo_than_1) then
        call report_err ('more than one choice for bottom boundary condition')
      end if
      if (rd_ext) call end_file

c     initial conditions
      rd_ext =rd_logic()
      if (rd_ext) call init_file('watinit.wi')
      eq_grwt = rd_logic() 
	
c	equilibrium groundwater table
      if (eq_grwt) then
		if (ibboco.le.3) then
			do i=1,ncs
				ph(i)=-x(i)+gwl(1)
			enddo
		else 
			call report_err
     $			('can not calculate ph profile in equilibrium with groundwater
     $			table: no ground water table')
		endif 
      else
		igiveph = rd_logic()
c		pressure head values given
		if (igiveph) then 
			test_min(1) = -mx_ph
			test_max(1) =  mx_ph
			call interp_profile (result, ncs, 1, 'initial ph(z)'
     $                   ,test_min, test_max,9)
			do i = 1, ncs
				ph(i) = 10.0d0 *  result(1,i)
			enddo
c		water content values given 
		else  
			test_min(1) = 0.d0
			test_max(1) =  1.d0
			call interp_profile (result, ncs, 1, 'initial wc(z)'
     $                   ,test_min, test_max,9)
			do i = 1, ncs
				wcr = calc_wc(-mx_real, i)
				wcs = calc_wc(0.d0, i)
				if (result(1,i).ge.wcs*1.1.or.result(1,i).le.wcr*0.9) then
					call report_err
     $					('wc(i)_initial larger than wcs or smaller than wcr')
					wc(i) = (wcr+wcs)/2.d0
				else
					wc(i) = result(1,i)
				end if
			enddo
		endif
      endif
      if (rd_ext) call end_file
c     water: printing control
c     iteration history (wat_hist.out)
      pr_wat_histor = rd_logic()
      if (pr_wat_histor) call open_file (histor_file, 'wat_hist.out')
c     summary table (wat_sum.out)
      pr_wat_sum = rd_logic()
c     cumulative terms of the water balance (wat_cum.out)
      pr_wat_cum = rd_logic()
c     evapotranspiration and level of the ground water table (wat_et.out)
      pr_wat_et = rd_logic()
c     pressure head (ph.out)
      pr_wat_tseries(1) = rd_logic()
c	water content (wc.out)
      pr_wat_tseries(2) = rd_logic()
c	root extraction (crtex.out)
      pr_wat_tseries(5) = rd_logic()
      call end_file

c	postprocessing
c	compartments under watertable: ph = 0
      if(ibboco.le.3) then
		dxn=0.d0
		i= 0
30		if(dxn.gt.gwl(1))then
			if(i.lt.ncs)then
				i = i+ 1
				dxn=dxn-dxinter(i)
				goto 30
			else
				call report_err
     $			('ground water table depth exceeds soil profile depth')
			endif
		endif
		do  j=i,ncs
			ph(j)=0.0d0
		enddo
      endif
c     initialise wc and ph
      if (igiveph.or.eq_grwt) then
		do i=1,ncs
			wc(i)=calc_wc(ph(i),i)
		enddo
      else
		do i=1,ncs
			ph(i)=calc_ph(wc(i),i)
		enddo
      endif
c     water storage
      vol = 0.0d0
      do i=1,ncs
		vol=vol+wc(i)*dx
	enddo
      voli=vol
      visd=0.d0
      do i=1,isd
		visd=visd+wc(i)*dx
	enddo 
	
      visdi = visd
      end