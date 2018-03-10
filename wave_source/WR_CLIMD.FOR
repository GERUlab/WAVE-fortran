c	Version 3.0
c	Author: m. vanclooster 22/02/2002
c###################################################################################
      subroutine wr_climdata
c     in   : et0, fin, hsh, isucr, ndsim, prec, rirr, simnit, simtemp,
c            tinit, tmax, tmin
c     out  : -
c     calls: close_file, date, heading, open_file
c     correction 10.3.94 also output min and max temp when only temperature is modeled
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'clim.com'
      integer res_file
      integer day, month, year

      call open_file(res_file,'climdata.res')
      call heading(res_file,'climdata.res','climatological data input')
      if (isucr) then
		write (res_file,20)
      else if (simtemp) then
		write (res_file,30) 
      else
		write (res_file,40)
      endif
      do 10 i= 1, ndsim
		call date_wave (i+idint(tinit)-1, year, month, day)
		if (isucr) then
			write (res_file,50) year, month , day, et0(i),prec(i),
     $		rirr(i),fin(i), tmin(i),tmax(i),hsh(i)
		else if (simnit.or.simtemp) then
			write (res_file,60) year, month, day, et0(i),prec(i),
     $		rirr(i), fin(i), tmin(i),tmax(i)
		else
			write (res_file,70)year, month, day, et0(i),prec(i),
     $		rirr(i), fin(i)
		endif
10    enddo 
      call close_file(res_file)

20    format( 'year month day    et0     prec    rirr    fin     ',
     $'tmin   tmax     hsh' ,/,
     $        '                  (mm)    (mm)    (mm)    (mm)     (c)',
     $'   (c) (j/cm**2day)')
30    format('year month day     et0      prec    rirr   fin     ',
     $'tmin   tmax'
     $	  ,/,'                   (mm)     (mm)    (mm)   (mm)     ',
     $'(c)    (c)') 
40    format( 'year month day      et0      prec    rirr   fin'
     $     ,/,'                   (mm)     (mm)    (mm)   (mm)',/)
50    format ( i4,2x,i2,3x,i2,1x,7(f7.1,1x))
60    format ( i4,2x,i2,3x,i2,1x,6(f7.1,1x))
70    format ( i4,2x,i2,3x,i2,1x,4(f7.1,1x))
      end
