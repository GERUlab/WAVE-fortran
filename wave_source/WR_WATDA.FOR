c	version 3.0
c	author: m. vanclooster 22/02/2002	
c###################################################################################
      subroutine wr_watdata
c     in   : arel, arer, brel, brer, drz, 
c            eq_grwt, gprh, gwl, hystlim, ibbcct, ibboco, 
c            idvs, igiveph, ilinr, ismxnl, isucr, laict, modelk, modelpf, 
c            ncs, ndsim, nparkon, nparpf, npl, pond_max, p0, p1, p2h, p2l, p3,
c            parakon, parapf, phsa, ph, pr_wat_cum, pr_wat_et, 
c            pr_wat_histor, pr_wat_sum, pr_wat_tseries, qdeep, rlai, 
c            rootct, simplant, rt_distr, tb, te, tinit, wc, x
c     out  : -
c     calls: big_title, close_file, date, heading, open_file, small_title
c###################################################################################
      implicit double precision (a-h,o-z)
      include 'constant'
      include 'gen.com'
      include 'wat.com'
      integer res_file
      integer day, month, year

      call open_file (res_file,'watdata.res')
      call heading (res_file,'watdata.res','input for modeling water transport')

c     water: parameters for each soil layer
      if (multimodal) then
		call wr_init_n_p (res_file)
      else
		call big_title (res_file,'parameters for each soil layer')
		write(res_file,10) modelpf
		write(res_file,20) modelk
		call small_title (res_file, 'parameters for mrc model')
		write(res_file,30) (n, n = 1, nparpf)
		do i=1,npl
			write(res_file,40) i,(parapf(i,j),j=1,nparpf)
		enddo
		call small_title (res_file,'parameters for conductivity model')
		write(res_file,50) (n, n = 1, nparkon)
		do  i=1,npl
			write(res_file,60) i,(parakon(i,j),j=1,nparkon)
		enddo
		if (modelpf.gt.1) write (res_file,70) hystlim
      endif
10    format('model type for the mrc                :',i5)
20    format('model type for the conductivity curve :',i5)
30    format(' layer',10(4x,i1,6x))
40    format(1x,i3,1x,10(e10.3,1x))
50    format(' layer',10(4x,i1,6x))
60    format(1x,i3,1x,10(e10.3,1x))
70    format('minimum (relative) pressure head change for return',
     $'point (hysteresis models) : ',f8.5)

c     water: upper boundary conditions
      call big_title(res_file,'upper boundary condition')
      write(res_file,80) phsa,pond_max
80    format(/,
     $'the minimum allowed pressure head at the soil surface    :',
     $e10.2,' (cm)',/,
     $'the maximum ponding depth',f7.2,' (mm)')

c     water: kc factors 
      call big_title (res_file,'sink (evaporation and plant water uptake)')
      call small_title(res_file,'kc values')
      if (idvs) then
		write (res_file,90)
      else
		write (res_file,100)
      endif
      do i = 1, nr_of_kc_values
		if (idvs) then
			write (res_file, 110) dvs_kc(i), kc(i)
		else
			call date_wave(iday_kc(i),year,month,day)
			write (res_file, 120) year, month, day, kc(i)
		endif
      enddo
90    format('kc as a function of dvs')
100   format('kc as a function of time')
110   format( f5.2, 1x, f7.3)
120   format(i4, 1x, i2, 1x, i2, 1x, f7.3)

      if (simplant) then
c		inactivation of the roots
		call small_title(res_file,'inactivation of the roots')
		call date_wave(idint(tb),year,month,day)
		write (res_file,130) year, month,day
		call date_wave(idint(te),year,month,day)
		write (res_file,140) year, month,day
c		lai and rooting depth
		if(.not.isucr) then
			call small_title(res_file,'lai (sucros is not used)')
			write (res_file,150) 
			do i = 1,ndsim
				call date_wave(i+idint(tinit)-1,year,month,day)
				write (res_file,160)year,month,day,rlai(i)
			enddo
			call small_title(res_file,'rooting depth (sucros is not used)')
			write (res_file,170)
			do i = 1,ndsim
				call date_wave(i+idint(tinit)-1,year,month,day)
				write (res_file,180)year,month,day,drz(i)
			enddo
		endif
c		water: water uptake
		call small_title(res_file,'water uptake reduction')
		if(ilinr) then
			write(res_file,190)
		else
			write(res_file,200)
		endif
		write(res_file,210) p0/10.d0,p1/10.d0,p2h/10.d0,p2l/10.d0,
     $	p3/10.d0
		call small_title(res_file,'root distribution per compartment')
		if(ismxnl) then
			write(res_file,220) brer,arer
		else
			write(res_file,230)
		endif
		write (res_file,240)
		do i = 1,ncs
			write (res_file,250)i,rt_distr(i)
		enddo 
      endif
130   format('date at which roots start to become inactive   : '
     $  ,i4,1x,i2,1x,i2)
140   format('date at which roots reach their max inactivity : '
     $  ,i4,1x,i2,1x,i2)
150   format('year month day  lai'/,
     $           '              (m**2/m**2)')
160   format (i4,2x,i2,2x,i2,1x,f8.2)
170   format('year month day  rooting depth',/,
     $    '                   (mm)')
180   format (i4,2x,i2,2x,i2,1x,f10.3)
190   format(
     $  'linear relation between reduction factor for water',
     $  ' uptake and pressure head')
200   format(
     $  'nonlinear relation between reduction factor for water',
     $  ' uptake and pressure head')
210   format('pressure head below which water uptake ...',/,
     $  'starts                                    : ',f12.1,' (cm)',/,
     $  'is optimal                                : ',f12.1,' (cm)',/,
     $  'is not optimal at high evaporative demand : ',f12.1,' (cm)',/,
     $  'is not optimal at low evaporative demand  : ',f12.1,' (cm)',/,
     $  'ceases (wilting point)                    : ',f12.1,' (cm)')
220	format(
     $  'root distribution is a linear function of depth'
     $  ,/,'rt_distr =',e10.2,' * dabs(depth) ',sp,e10.2,
     $  10x,'(1/day)')
230   format(
     $  'root distribution is specified for each compartment')
240   format('comp    rt_distr ',/,
     $  '      (1/day)')
250   format( i3,1x,f9.5)

c     water: bottom boundary condition
      call big_title(res_file,'bottom boundary condition')
      if(ibboco.eq.1) then
		write(res_file,260)
          do i = 1, ndsim 
			call date_wave(i+idint(tinit)-1,year,month,day)
			write(res_file,270) year, month , day, gwl(i)
		enddo
	elseif(ibboco.eq.2) then
		write(res_file,280) gwl(1)
          write(res_file,290)
          do  i = 1, ndsim 
			call date_wave(i+idint(tinit)-1,year,month,day)
			write(res_file,300) year, month , day, qdeep(i)
		enddo
      elseif(ibboco.eq.3) then
		write(res_file,310) arel,brel,gwl(1)
      elseif(ibboco.eq.4) then
		write(res_file,320)
          do i = 1, ndsim 
			call date_wave(i+idint(tinit)-1,year,month,day)
			write(res_file,330) year, month , day, gprh(i)/10.d0
		enddo
      elseif(ibboco.eq.5) then
		write(res_file,340)
      elseif(ibboco.eq.6) then
		write(res_file,350)
      elseif(ibboco.eq.7) then
		write(res_file,360)
      endif
260   format('the groundwater level: ',/,
     $'year mo day   level',/,
     $'               (mm)')
270   format(i4,1x,i2,1x,i2,1x,f10.1)
280   format('the initial groundwater level :',f10.2,' mm')
290   format('the flux of water through the bottom:',/,
     $	'year month day    flux',/,
     $	'                (mm/day)')
300   format(i4,1x,i2,1x,i2,1x,f10.1)
310   format('flux-ground water level relationship: ',/,1x,
     $	'flux =',e10.3,' *exp(',e10.3,' *abs(ground water level))  (mm/
     $	day)',/,
     $	' initial groundwater level:',f10.0,' mm')
320   format(' the pressure head of the bottom compartment: ',
     $	/,'year month day    pr.head',/,
     $	'                   (cm)')
330   format(i4,1x,i2,1x,i2,1x,f10.1)
340   format('zero flux at the bottom of the soil profile')
350   format('free drainage at the bottom of the soil profile')
360   format('lysimeter bottom boundary condition')

c     water: initial values
      call big_title(res_file,'initial values')
      if(eq_grwt)then
		write(res_file,370)
      else
		if(igiveph) then
			write(res_file,380)
		else
			write(res_file,390)
		endif
      endif
      call small_title(res_file,'initial wc / ph profile')
      write(res_file,400)
      write(res_file,410) (x(i),ph(i)/10.d0,wc(i),i=1,ncs)
370   format('the pressure head profile is calculated in equilibrium',
     $' with the groundwater level')
380   format('the pressure head profile is given')
390   format('the water content profile is given')
400   format(5x,'depth',9x,'ph',6x,'wc',/,
     $7x,'(mm)',7x,'(cm)',1x,'(m**3/m**3)')
410   format(1x,f10.1,f12.1,2x,f6.3)

c     water: printing control
      call big_title(res_file,'output')
      if(pr_wat_histor ) then
		write (res_file,420)plus
      else
		write (res_file,420)neg
      endif
      if(pr_wat_sum )then
		write (res_file,430)plus
      else
		write (res_file,430)neg
      endif
      if(pr_wat_cum)then
		write (res_file,440)plus
	else
		write (res_file,440)neg
      endif
      if(pr_wat_et )then
		write (res_file,450)plus
      else
		write (res_file,450)neg
      endif
      if(pr_wat_tseries(1) )then
		write (res_file,460)plus
      else
		write (res_file,460)neg
      endif
      if(pr_wat_tseries(2))then
		write (res_file,470)plus
      else
		write (res_file,470)neg
      endif
      if(pr_wat_tseries(5))then
		write (res_file,480)plus
      else
		write (res_file,480)neg
      endif
420   format
     $('the iteration history (wat_hist.out)                 : ',a)
430   format
     $('overview of main state variables (wat_sum.out)       : ',a)
440   format
     $('cumulative terms of the water balance (wat_cum.out)  : ',a)
450   format
     $('evapotransp. and gr. water table level               : ',a)
460   format
     $('the pressure head (ph.out)                           : ',a)
470   format
     $('the water content (wc.out)                           : ',a)
480   format
     $('the cumulative root extraction (crtex.out)           : ',a)
      call close_file (res_file)
      return
      end