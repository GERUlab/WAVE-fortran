c######################################################################################
      subroutine gendata
c	version 3.0
c	author m. vanclooster 15/02/2002
c     in   : -
c     out  : devstop, dsp, dt, dthm, dtmax, dtmin, dx, dxinter, from,
c            harvest_date, iprct, isd, isucr, ncs, ncsl, ndsim, nla, npl,
c            npr, nr_of_points, plant_date, printing_time, range_err,
c            simnit, simplant, simsol, simtemp, tend, tinit, to, tprint,
c            x, yearin
c     calls: check_et, date, end_file, hold_line, init_file, jdate, next_line,
c            rd_date, rd_int, rd_logic, rd_real report_err, set_range_checking,
c            sum
c######################################################################################
      implicit double precision (a-h,o-z)
      include	'constant'
      include	'gen.com'
      integer	rd_int
      integer	i, j,nr_comps, jdate,day, month, year, year_for
      double precision rd_real
      logical  rd_logic

      call init_file ('gendata.in') 

c	simulation type
      simplant = rd_logic()
      if (simplant) isucr = rd_logic()
      simtemp = rd_logic()
      simsol = rd_logic()
      simnit = rd_logic()
      simsol = simsol.or.simnit
      simtemp = simtemp.or.simnit
      isucr = isucr.and.simplant
	
c	soil profile development
c	compartment size (mm) 
      dx = rd_real (mn_real,mx_real)
c	number of soil layers
      npl=rd_int(1 ,kt_layrs)
      nr_comps = kt_comps
      do i=1,npl
c     number of compartments for each soil layer
		call hold_line
		lyread = rd_int(1,npl)
		if (lyread.ne.i) then
			call report_err('error in the sequence of the compartments')
		else
			ncsl (i) = rd_int (1,kt_comps)
		end if
		call next_line
		nr_comps = nr_comps - ncsl(i)
	enddo

c     test if sum of compartments for each layer > total number of compartments available
      if (nr_comps.lt.0) then
		call report_err ('too many compartments')
      endif
      call check_et

	ncs =0
      do i=1,npl
		do  j = 1,ncsl(i)
			ncs = ncs +1
			nla(ncs) = i
		enddo 
      enddo

c     added for solute >> nla(ncs +1) needed
      nla(ncs+1) = nla(ncs)

c     initial calculation of the array containing
c     possible variable compartiment size, depth to the middle of each node 
      dxinter(1)=dx/2.d0

      x(1)=-dxinter(1)

      do i= 2, ncs
		dxinter(i)=dx
		x(i)=x(i-1)-dxinter(i)
      enddo 
      dxinter(ncs+1)=dx/2.d0
c     calculation of depth of the soil profile
      dsp=-x(ncs)+dx/2.

c	general soil parameters (bulk density, ...)
      if (simsol.or.simtemp) then
		do i=1,npl
			call hold_line
			lyread = rd_int(1,npl)
			if (lyread.ne.i) then
				call report_err
     $			('error in the sequence of the layers bulk density')
			else
				bulk_dens(i) = rd_real (mn_real, mx_real)
			endif
			call next_line
		enddo
      call check_et
      endif

c	simulation time variables
c     calculations ( year month day )
      call hold_line
      yearin = rd_int (1900, 3000)
      init_year=yearin
      call rd_date (month,day)
      init_month=month
      init_day=day
      call next_line
      tinit = dfloat(jdate(yearin,month,day,yearin))
      call hold_line
      year = rd_int (1900, 3000)
      call rd_date (month,day)
      call next_line
      mo_1_year = year.gt.yearin
      tend = dfloat(jdate(year,month,day,yearin))
      ndsim = idint(tend) - idint(tinit) + 1
      if (simplant) then
c       planting or emergence date 
        call hold_line
        call rd_date (month, day)
        call next_line
        plant_date = jdate(year_for(month,day),month,day,yearin)
c       harvesting date 
        call hold_line
        call rd_date (month, day)
        call next_line
        harvest_date = jdate(year_for(month,day),month,day,yearin)
      endif

c	parameters concerning the numerical solution
c     maximum time step (days) ( > 0.001 days)
      dtmax = rd_real(0.001d0 ,366.d0)
c     minimum time step (days) ( < 1 day )
      dtmin = rd_real(mn_real,dtmax)
c     maximum change of the moisture content duringone time step (m**3/m**3) (0.002<dthm<0.05)
      dthm = rd_real(0.002d0,0.05d0)
c     maximum allowed balance error for each compartiment(0.00001<stop<1) (m*3/m*3/day)
      devstop = rd_real(0.00001d0,1.d0)
      dt = 1.e-10

c	parameters concerning printing and plotting
c     check range errors (y/n)
      range_err = rd_logic()
      if (range_err) call set_range_checking
c     constant time increment for printing(y/no)
      iprct = rd_logic()
      if(iprct) then
c		time increment
		tprint =rd_int (1,366)
		do i= 1,ndsim,tprint 
			printing_time(i) = .true.
		enddo 
      else
c		printing dates 
		npr= rd_int (0,366)
		do i=1,npr
			call hold_line
			call rd_date (month, day)
			idate = jdate(year_for (month,day),month,day,yearin)
			call next_line
			if (idate.lt.idint(tinit).or.idate.gt.idint(tend)) then
				call report_err ('printing date outside simulation time interval')
			else
				printing_time(idate-idint(tinit)+1) = .true.
			endif
		enddo
		call check_et
      endif
c     integration profile (ids) 
      isd = rd_int (1, ncs)
c	upper and lower compartments
      nr_of_points = rd_int(0,ncs )
      do i=1,nr_of_points
		call hold_line
		from(i)=rd_int(1,ncs)
		to(i)=rd_int(from(i),ncs)
		call next_line
	enddo 
      call check_et
      call end_file
      
	return
      end