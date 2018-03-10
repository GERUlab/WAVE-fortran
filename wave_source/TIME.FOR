c	version 3.0
c	m. vanclooster 18/02/02
c###################################################################################
      integer function jdate(year,month,day,init_year)
c###################################################################################
      integer ndm(12), sum,year,month,day,init_year
      logical ok
      parameter (ok = .true.)
      data ndm/31,28,31,30,31,30,31,31,30,31,30,31/
      
	if(mod(init_year,4).eq.0) then
		init_sum=366
      else
		init_sum=365
      end if
      if(mod(year,4).eq.0) then
		ndm(2)=29
      else
		ndm(2)=28
      endif
      sum=0
      do i=1,month -1
		sum=sum+ndm(i)
      enddo
      if(year.eq.init_year) then
		jdate=sum+day
      else
		if(year.gt.init_year) then
			jdate=sum+init_sum+day
		else
			call report_err
     $			('year for calculation of the julian day number '//
     $		'before the initial year')
			call stop_simulation
     $			('programme stopped: check err_file')
		endif
      endif
      return
      end

c###################################################################################
      subroutine date_wave(j_date, year, month, day)
c###################################################################################
      implicit double precision (a-h, p-z)
      include   'constant'
      include   'gen.com'
      integer dy_nr, j_date, year, month, day, ndm(12),days_gone, dy_yearin
      data ndm/31,28,31,30,31,30,31,31,30,31,30,31/

      if (mod(yearin,4).eq.0) then
		dy_yearin = 366
      else
		dy_yearin = 365
      endif
      if (j_date. le. 0) then
		call report_err
     $		('julian day number for calculation date negative') 
      else if (j_date. gt. dy_yearin) then
		dy_nr = j_date - dy_yearin
		year = yearin + 1
      else 
		dy_nr = j_date
		year = yearin
      endif
      if(mod(year,4).eq.0) then
		ndm(2) = 29
      else
		ndm(2) = 28
      endif
      month = 1       
10    days_gone = ndm(month)
      if (dy_nr.gt. days_gone) then
		dy_nr =dy_nr - days_gone
		month = month + 1
		goto 10
      endif
      day = dy_nr 
      return
      end

c###################################################################################
      integer function day_nr (year,month,day,iyear)
c###################################################################################
      implicit double precision (a-h, p-z)
      include   'constant'
      include   'gen.com'
      integer year,month,day,iyear, jdate
      integer t_init, t_actual
      
	t_init = idint(tinit)
      t_actual = jdate (year,month,day,iyear)
      if (year.ge.iyear) then
		day_nr = t_actual - t_init +1
      else
		day_nr = t_actual - t_init 
      endif
      return
      end

c###################################################################################
      integer function day_nr_no_yr (month,day)
c###################################################################################
      implicit double precision (a-h, p-z)
      include   'constant'
      include   'gen.com'
      integer month,day, year_for, day_nr
      
	day_nr_no_yr = int(day_nr(year_for(month,day),month,day,yearin))
      return
      end    

c###################################################################################
      integer function year_for (month,day)
c###################################################################################
      implicit double precision (a-h, p-z)
      include   'constant'
      include   'gen.com'
      integer month,day
      logical first_time
      data first_time/.true./

      if (first_time) then
		first_time = .false.
		call date_wave (idint(tinit), init_year, init_month,init_day)
      endif
      if(mo_1_year.and.month.lt.init_month) then
		year_for = yearin+1
      else if(mo_1_year.and.month.eq.init_month.and.day.lt.init_day) then
		year_for = yearin+1
      else
		year_for = yearin
      endif
      return
      end    

c###################################################################################
      subroutine calcdt
c     in   : dt, dthm, dtmax, dtmin, dx, flxa, flxs, iter, ncs, ndts, phsa,
c            phsurf, rtex, t, tmd
c     out  : dt, ndts
c     calls: report_err, stop_simulation
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'

c     define maximum timestep in terms of iteration (only if maxiter>1 and if ndts>0)
      if((maxiter.gt.1).and.(ndts.gt.0))then
		fac=0.5d0+0.25d0*(7-iter)
		fac= dmax1(0.5d0,fac)
		dtmax1=dt*fac
      else
		dtmax1=dtmax
      endif
c     define maximum time step in terms of water content change 
      dtdt=1.0e-10
      if(ndts.le.1)then
		if((flxa.gt.0.).and.(phsurf.le.phsa))then
			flxsn=-(rtex(1)+(flxs(1)-flxs(2))/dx)
		elseif((flxa.lt.0.).and.(phsurf.gt.0.d0))then
			flxsn=-(rtex(1)+(flxs(1)-flxs(2))/dx)
		else
			flxsn=-(rtex(1)+( flxa  -flxs(2))/dx)
		endif
      else
		flxsn=-(rtex(1)+(flxs(1)-flxs(2))/dx)
      endif
      dtdt= dmax1(dtdt, dabs(flxsn))
      do i=2,ncs
		flxsn=-(rtex(i)+(flxs(i)-flxs(i+1))/dx)
		dtdt= dmax1(dtdt, dabs(flxsn))
      enddo
c     define maximum time step in terms of solute transport 
      if (simsol) then
		dtmax_sol = dt_sol(dtmax)
      else
		dtmax_sol = dtmax
      endif
c     define actual dt
      dt= dmax1(dtmin, dmin1(dthm/dtdt,dtmax,dtmax1,dtmax_sol))
      ndts=ndts+1
      if(ndts.gt.9999) then
		call report_err 
     $		('too many timesteps for one day: check iteration history')
		call stop_simulation
     $		('program stopped: check err_file')
      endif
      if(t+dt.ge.tmd)dt=tmd-t
      return
      end

c###################################################################################
      subroutine halve (no_convergence)
c     in   : dt, dtmin, histor_file, nrepeat, nr_no_conv  ???, pr_wat_histor
c     out  : dt, no_convergence, nrepeat, nr_no_conv
c     calls: -
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      logical no_convergence

      if(dt.gt.dtmin)then
		nrepeat=nrepeat+1
		if(pr_wat_histor)write(histor_file,10) maxiter
10		format ('timestep halved after ', i3, ' iterations')
		dt= dmax1((dt/2.),dtmin)
		no_convergence = .false.
      else
		nr_no_conv=nr_no_conv+1
		no_convergence = .true.
      endif
      return
      end