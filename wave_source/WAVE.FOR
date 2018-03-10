c	Version 3.0
c	Author: m.vanclooster (22/02/2002)
c###################################################################################
      program wave
c     in   : dt, isucr, ndsim, simnit, simsol, simtemp, tinit
c     out  : nday, ndts, nts, t, tmd
c     calls: calcdt, calc_wat_sol, calc_wat_temp, climdata, cropdata, crop_print, date,
c            gendata, input_error, nitdata, nit_boco, nit_intgr, nit_print, nit_sink, 
c            prhead, soiltemp, soldata, soln, sol_boco, sol_intgr, sol_print , sol_sink, 
c            stop_simulation, sucros, tempdata, temp_print, watdata, wat_boco, wat_intgr, 
c            wr_climdata, wr_cropdata, wr_gendata, wr_nitdata, wr_soldata, wr_tempdata, 
c            wr_watdata
c###################################################################################
      implicit double precision (a-h,o-z)
      include 'constant'
      include 'gen.com'
      include 'wat.com'
      logical input_error
      integer year, month, day

c	prompt interface
	call program_version
c	read input
	call reading
	
c	start daily loop
      nts = 0
      t=tinit

      do nday =1, ndsim
		call date_wave (idint(t),year,month,day)
		write (*, 10) year, month, day
10		format ('+calculations for day: ',i4,' / ', i2,' / ', i2)
c		new day
		ndts=0
		tmd = t +1.0d0
		if (isucr)call sucros
		call wat_boco
		if (simsol)then
			call sol_boco
			if (simnit) then
				call nit_boco
			endif
		endif
20		nts=nts+1
		call calcdt
c		water module
		call prhead
		call wat_intgr
c		heat module 
		if (simtemp) then
			call calc_wat_temp
			call soiltemp 
		endif

c		solute module 
		if (simsol) then
			call calc_wat_sol
			call sol_sink
			if (simnit) call nit_sink
			call soln
			call sol_intgr
			if (simnit) call nit_intgr
		endif
		t=t+dt
		if (t.lt.tmd) goto 20  
c		end of time step loop
		call wat_print
		if(isucr)call crop_print
		if(simsol)call sol_print 
		if(simtemp)call temp_print
		if(simnit)call nit_print
30	enddo
c     end of day loop
      call prep_nxt_yr
      call stop_simulation ('programme finished execution normally')
      end