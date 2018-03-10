c	version 3.0
c	author: m.vanclooster (18/02/02)
c###################################################################################
      subroutine soldata
c     in   : dx, imobsw, isd, ncs, nla, npl, wcio, wcmo, wco, yearin
c     out  : acsolio, acsolmo, acsolo, conirsol, csolio, csolmo, csolo,
c            ddepsol, dsp, fsol, idinof, idirr, idplow, ncsplow, ninof,
c            nirr, nplo, nr_of_sol, out_sol, parasol, pr_sol_cum, pr_sol_sum,
c            pr_sol_tseries, rates, tcsolio, tcsolmo, tcsolo, tsoli,
c            tsoli_isd, wdepsol
c     calls: calc_wat_sol, check_et, date, end_file, hold_line, init_file, jdate,
c            next_line, rd_date, rd_int, rd_logic, rd_real, report_err
c###################################################################################
      implicit double precision (a-h, o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      integer applicat, day, month,year_for,sp,sp_rd, rd_int
      dimension ppdepth(kt_plo)
      dimension test_min(kt_in_row-1),test_max(kt_in_row-1)
      dimension result (kt_in_row, 366)
      logical rd_logic,rd_ext
      logical sol_layer(kt_layrs)

	do i = 1, kt_ninof
		do j = 1, kt_sol
			fsol(i,j) = 0.d0
		enddo
	enddo
      call init_file ('soldata.in')
c     solute: number of species
      nr_of_sol = rd_int(1,kt_sol)
c     immobile/mobile water ?      
      imobsw = rd_logic()
c     solute:  parameters for each soil layer
      rd_ext =rd_logic()
      if (rd_ext) call init_file('solpar.sp')
      do  i = 1, npl
		sol_layer (i) = .false.
      enddo
      do i=1,npl
		call hold_line
		layer = rd_int (1,npl)
		if (sol_layer(layer)) then
			call report_err
     $		(' more than one set of solute parameters for one layer')
		else
			sol_layer(layer) =  .true.
			if(imobsw) then 
				do j = 6,9
					parasol(j,i,1) = rd_real(0.d0,mx_real)
				enddo
			else
				parasol(6,i,1) = rd_real(0.d0,mx_real)
			endif 
		endif
		call next_line
      enddo
      call check_et

c     input a table of parameters per layer for each solute species
      do sp= 1, nr_of_sol
		do i = 1, npl
			sol_layer (i) = .false.
		enddo
		do i=1,npl
			call hold_line
			layer = rd_int (1,npl)
			if (sol_layer(layer)) then
				call report_err
     $			(' more than one set of solute parameters for one layer')
			else
				sol_layer(layer) =  .true.
				do j = 2,5 
					parasol(j,i,sp) = rd_real(0.d0,mx_real)
				enddo
			endif
			call next_line
		enddo
		call check_et
      enddo
      if (rd_ext) call end_file 

c	assign parasol parameters
      do i=1,npl
		parasol(1,i,1) = bulk_dens(i)
		do sp=2,nr_of_sol
			parasol(1,i,sp) = parasol(1,i,1)
			do j = 6,9
				parasol(j,i,sp)=parasol(j,i,1)
			enddo
		enddo
	enddo
      if (.not.imobsw) then
		do sp = 1,nr_of_sol
			do j = 1,npl
				parasol(7,j,sp) = 1.d0
				parasol(8,j,sp) = 0.d0
				parasol(9,j,sp) = 1.d0
			enddo
		enddo
      endif

c     solute : upper boundary conditions
      rd_ext =rd_logic()
      if (rd_ext) call init_file('solubc.su')

c	wet and dry deposition 
      wet_and_dry = rd_logic()
      if (wet_and_dry) then
		do sp = 1, nr_of_sol
			call hold_line
			sp_rd = rd_int (1,nr_of_sol)
			if (sp_rd.ne.sp) then
				call report_err
     $			('error in species numbers for dry and wet depositions')
			else
				wdepsol(sp) = rd_real(0.d0,mx_real)
				ddepsol(sp) = rd_real(0.d0,mx_real)
			endif
			call next_line
		enddo 
		call  check_et
      endif

c     inorganic fertiliser
      ninof = rd_int(0, kt_ninof)
      do i=1,ninof
		call hold_line
		applicat = rd_int (1,ninof)
		if (applicat.ne.i) then
			call report_err
     $		("error in sequence of fertilizer applications")
		else
			call rd_date (month, day)
			idinof(i) = jdate(year_for(month,day),month,day,yearin)
			do sp =1, nr_of_sol
				fsol(i,sp)= rd_real (0.d0,mx_real)
			enddo
		end if
		call next_line
      enddo
      if (ninof.gt. 0) call check_et

c     irrigation 
      nirr = rd_int(0, kt_nirr)
      do i=1,nirr
		call hold_line
		applicat = rd_int (1,nirr)
		if (applicat.ne.i) then
			call report_err ("error in sequence of irrigations")
		else
			call rd_date (month, day)
			idirr(i) = jdate(year_for(month,day),month,day,yearin)
			do sp = 1, nr_of_sol
				conirsol(i,sp)= rd_real (0.d0,mx_real)
			enddo
		endif
		call next_line
      enddo
      if (nirr.gt.0) call check_et

c     plowing
      nplo = rd_int (0,kt_plo)
      do i = 1, nplo
		call hold_line
		nr = rd_int (1, nplo)
		if (nr.ne.i) then 
			call report_err("error in sequence of plowing numbers")
		else
			call rd_date (month, day)
			idplow(i) = jdate(year_for(month,day),month,day,yearin)
			ppdepth(i)= rd_real(-dsp,dsp)
		endif
		call next_line
      enddo
      if (nplo.gt.0) call check_et
      if (rd_ext) call end_file
      do i =1, nplo
		depthplo = dabs(ppdepth(i))
		depth=0.d0
		ncsplow(i) = 0
80		if(depth.lt.depthplo) then
			depth=depth+dx
			ncsplow(i) = ncsplow(i) +1
			if(ncsplow(i).lt.ncs) goto 80
		endif
      enddo 

c     solute rate constants decomposition
      rd_ext =rd_logic()
      if (rd_ext) call init_file('solrates.sr')
      do sp = 1, nr_of_sol
		test_min(sp) = 0.d0
		test_max(sp) = mx_real
      enddo
      call interp_profile (result, ncs, nr_of_sol,
     $                     'solute decay rates',test_min, test_max,9)
      do i = 1,ncs
		do sp = 1, nr_of_sol
			rates (i,sp) =  result (sp,i)
		enddo 
      enddo
      if (rd_ext) call end_file

c     solute : initial values
      rd_ext =rd_logic()
      if (rd_ext) call init_file('solinit.si')
	do sp = 1, nr_of_sol
		test_min(sp) = 0.d0
		test_max(sp) = mx_real
      enddo
      call interp_profile (result, ncs, nr_of_sol,
     $                     'initial solute concentration',test_min, test_max,9)
      do i = 1,ncs
		do sp = 1, nr_of_sol
			tcsolo(i,sp) = result (sp,i)
		enddo
      enddo
      if (rd_ext) call end_file 

c     initialisation of solute profiles (sp)
      do sp = 1, nr_of_sol
		tsoli(sp) = 0.d0
		do i= 1, ncs
			tsoli(sp) = tsoli(sp) + tcsolo(i,sp)
			if (i.eq.isd) tsoli_isd(sp) = tsoli(sp)
		enddo
      enddo
      call calc_wat_sol
      volume=dx
      do sp = 1, nr_of_sol
		do i=1,ncs
			if(imobsw)then
				tcsolio(i,sp)=tcsolo(i,sp)*wcio(i)/wco(i)
				tcsolmo(i,sp)=tcsolo(i,sp)*wcmo(i)/wco(i)
				csolio(i,sp)=tcsolio(i,sp)/(volume*(parasol(2,nla(i),sp)*
     $			(1-parasol(9,nla(i),sp))
     $			*parasol(1,nla(i),sp)+wcio(i)))
				csolmo(i,sp)=tcsolmo(i,sp)/(volume*(parasol(2,nla(i),sp)*
     $			parasol(9,nla(i),sp)
     $			*parasol(1,nla(i),sp)+wcmo(i)))
				csolo(i,sp)=(csolio(i,sp)*wcio(i)+csolmo(i,sp)*wcmo(i))/
     $			wco(i)
				acsolio(i,sp)=parasol(2,nla(i),sp)*
     $			(1-parasol(9,nla(i),sp))*csolio(i,sp)
				acsolmo(i,sp)=parasol(2,nla(i),sp)*
     $			parasol(9,nla(i),sp)*csolmo(i,sp)
				acsolo(i,sp)=acsolio(i,sp)+acsolmo(i,sp)
			else
				tcsolio(i,sp)=tcsolo(i,sp)*wcio(i)/wco(i)
				tcsolmo(i,sp)=tcsolo(i,sp)*wcmo(i)/wco(i)
				csolo(i,sp)=tcsolo(i,sp)/(volume*(parasol(2,nla(i),sp)*
     $			parasol(1,nla(i),sp) +wco(i)))
				csolio(i,sp)=0.0d0
				csolmo(i,sp)=csolo(i,sp)
				acsolio(i,sp)=0.0d0
				acsolo(i,sp)=parasol(2,nla(i),sp)*csolo(i,sp)
				acsolmo(i,sp)=acsolo(i,sp)
			endif
		enddo 
	enddo

c     solute: printing control
      do sp = 1, nr_of_sol
		call hold_line
		sp_rd = rd_int(1,nr_of_sol)
		if (sp_rd.ne.sp) then
			call report_err
     $		('error in sequence of compartments')
		else
			out_sol(sp) = rd_logic()
		endif
		call next_line
	enddo
      call check_et
      pr_sol_sum = rd_logic()
      pr_sol_cum = rd_logic()
      pr_sol_tseries(1) = rd_logic()
      pr_sol_tseries(3) = rd_logic()
      call end_file
      end