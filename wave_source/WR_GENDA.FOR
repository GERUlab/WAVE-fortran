c	version 3.0
c	author m. vanclooster 21/2/2002	
c###################################################################################
      subroutine wr_gendata
c     in   : devstop, dthm, dtmax, dtmin, dx, from, harvest_date,
c              iprct, isd, isucr, ncsl, ndsim, npl, npr, nr_of_points,
c            plant_date, printing_time, range_err, simnit, simplant,
c            simsol, simtemp, table, tend, tinit, to, tprint
c     out  : -
c     calls: big_title, close_file, date, heading, open_file, small_title
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      integer res_file
      integer day, month, year

      call open_file (res_file,'gendata.res')
      call heading (res_file,'gendata.res','general input for the model')

c     simulation type
      call big_title (res_file,'model switches')
      if (simplant) then
		write(res_file,10)plus
		if (isucr) then
			write(res_file,20)plus
		else
			write(res_file,20)neg
		endif
      else
		write(res_file,10)neg
      endif
      if (simtemp) then 
		write(res_file,30)plus
      else
		write(res_file,30)neg
      endif
      if (simsol)then 
		write(res_file,40)plus
      else
		write(res_file,40)neg
      endif
      if (simnit)then 
		write(res_file,50)plus
      else
		write(res_file,50)neg
      endif
10    format('plants     : ',a)
20    format('sucros     : ',a)
30    format('temperature: ',a)
40    format('solutes    : ',a)
50    format('nitrogen   : ',a)

c     write profile development
      call big_title(res_file,'profile development')
      write(res_file,60) dx,npl
      do i=1,npl
		write(res_file,70) i, ncsl(i)
	enddo 
60    format('compartment size            = ',f8.2,' mm',/
     $       'number of soil layers       = ',i5 )
70    format('layer',i3,' has: ' i3,' compartment(s)')
      call small_title(res_file,'bulk density (kg/l)')
      write (res_file,80)
80    format ('layer    1')
      do i=1,npl
		write (res_file,90) i,bulk_dens(i)
      enddo
90    format (i2,1x,f8.4)

c     simulation time variables
      call big_title(res_file,'simulation time variables')
      call date_wave (idint(tinit), year, month,day)
      write(res_file,100)year, month,day
      call date_wave (idint(tend), year, month,day)
      write(res_file,110)year, month,day
      if (simplant) then
		call date_wave (plant_date, year, month,day)
		write(res_file,120)year, month,day
		call date_wave (harvest_date, year, month,day)
		write(res_file,130)year, month,day
      endif
100   format('start of calculations: ',i4,1x,i2,1x,i2)
110   format('end of calculations  : ',i4,1x,i2,1x,i2)
120   format('planting date        : ',i4,1x,i2,1x,i2)
130   format('harvest date         : ',i4,1x,i2,1x,i2)

c     parameters concerning the numerical solution
      call big_title(res_file,'parameters concerning the numerical solution')
      write(res_file,140) dtmax,dtmin,dthm,devstop
140   format(
     $'maximum time step                   :',f10.2,' days',/,
     $'minimum time step                   :',f10.4,' days',/,
     $'maximum change of moisture content  :',f10.4,' m**3/m**3',/,
     $'maximum balance error               :',f10.5,
     $' m**3/m**3/day')

c     parameters concerning printing and plotting
      call big_title(res_file,'parameters concerning output')
      call small_title(res_file,'input testing')
      if (range_err) then
		write(res_file,150)
      else
		write(res_file,160)
      endif
150   format('input is tested for range errors')
160   format('input is not tested for range errors: ',
     $'only warnings are written to err_file')

c     plot dates
      call small_title(res_file,'output dates')
      if(iprct) then
		write(res_file,170) tprint
      else
		write (res_file, 180) npr
		do i = 1,ndsim
			if (printing_time(i)) then
				call date_wave (i+idint(tinit)-1, year, month,day)
				write (res_file,190) year, month, day
			endif
		enddo
      endif
170   format('the time increment between printing of ', /, 
     $'the summary table of the state variables is constant:',i4,
     $' day(s)')
180   format( 'number of dates output will be produced :',i3)
190   format(1x, i4, 1x, i2, 1x, i2)

c	isd
      write (res_file, 200) isd
200   format( 'the bottom compartment of the top layer, isd: ', i4)
      call small_title(res_file,'compartments for which there is output')
      write(res_file,210) 
      do i=1,nr_of_points
		write(res_file,220) from(i),to(i)
	enddo
210   format('  from   to') 
220   format(i5,1x,i5)
      call close_file(res_file)
      return
      end