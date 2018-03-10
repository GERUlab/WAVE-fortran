c######################################################################################
      subroutine climdata
c######################################################################################

c	version 3.0
c	author m. vanclooster 15/02/2002
c	corrections: m javaux 03/04/2002 
c     in   : isucr, ndsim, simtemp, yearin
c     out  : et0, fin, hsh, prec, rirr, tmax, tmin, hsh
c     calls: check_et, date, day_nr, end_file, hold_line, init_file, next_line,
c            rd_date, rd_real, report_err

      implicit double precision (a-h,o-z)
      include 'constant'
      include 'gen.com'
      include 'clim.com'
      integer year,month,day,day_nr,dn,rd_int

      call init_file('climdata.in')

      do 10 i= 1, ndsim
        call hold_line
        year = rd_int ( 1900, 3000)
        call rd_date (month,day)
        dn = day_nr(year,month,day,yearin)
        if (i.ne.dn) then
		call report_err ('error in date sequence')     
        else
		et0(i)  = rd_real(0.d0 ,  1000.d0)
		prec(i) = rd_real(0.d0 ,  1000.d0)
		rirr(i) = rd_real(0.d0 ,  100.d0)
		fin(i)  = rd_real(0.d0 ,  100.d0)
		if (simtemp.or.isucr) then
			tmin(i) = rd_real(-100.d0 , 100.d0)
			tmax(i) = rd_real(tmin(i), 100.d0)
		endif
		if (isucr) hsh(i) = rd_real(0.d0, mx_real)
        endif
        call next_line
10    continue

      call check_et
      call end_file 

      end