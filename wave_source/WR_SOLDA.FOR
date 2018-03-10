c	version 3.0
c	author m. vanclooster 21/2/2002	
c###################################################################################
      subroutine wr_soldata
c     in   : conirsol, ddepsol, fsol, idinof, idirr, idplow, imobsw,
c            ncs, ncsplow, ninof, nirr, npl, nplo, nr_of_sol, out_sol,
c            parasol, pr_sol_cum, pr_sol_sum, pr_sol_tseries, rates, simnit,
c            tcsolo, wdepsol
c     out  : -
c     calls: big_title, close_file, date, heading, open_file, small_title
c###################################################################################
      implicit double precision (a-h, o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      integer sp,appl
      integer res_file
      integer day, month, year
      character *(15) n_sol_name(3)
      integer  l_sol_name(3)
      character *(20) sm_t
      data n_sol_name /'ureum   ','ammonium','nitrate '/
      data l_sol_name /5,8,7 /

      call open_file (res_file, 'soldata.res')
      call heading (res_file, 'soldata.res',
     $'input for modeling solute transport')

c     number of species
      call big_title(res_file,'general information')
      write(res_file,10)nr_of_sol
      if (imobsw) then
		write(res_file,20)
      else
		write(res_file,30)
      endif
10    format ('number of solutes: ',i1)
20    format ('mobile immobile water is considered')
30    format ('mobile immobile water is not considered')

c     solute parameters for each soil layer
      call big_title(res_file,'soil specific transport parameters')
      write (res_file,40)
      do i=1,npl
		write (res_file,50) i,(parasol(j,i,1),j = 6,9)
      enddo
      call big_title(res_file,'solute specific transport parameters')
      do sp= 1, nr_of_sol
		if (simnit .and. sp.lt.4) then
			write (sm_t,60) n_sol_name(sp)(1:l_sol_name(sp))
			l_nt = 8 + l_sol_name(sp)
		else
			write (sm_t,70) sp 
			l_nt = 9
		endif
		call small_title(res_file,sm_t(1:l_nt))
		write (res_file,40)
		do i=1,npl
			write (res_file,50) i,(parasol(j,i,sp),j = 2,5)
		enddo
      enddo 
40    format ('layer    1      2      3        4   ')
50    format (i2,1x,4(f8.4,1x))
60    format ('solute: ',a)
70    format ('solute: ',i1)

c     solute upper boundary conditions
      call big_title(res_file,'upper boundary conditions')

c	wet and dry deposition 
      call small_title(res_file,'dry and wet deposition')
      if (wet_and_dry) then
		write (res_file,80)
		do sp = 1,nr_of_sol
			write (res_file,90) sp,wdepsol(sp),ddepsol(sp)
		enddo
      else
		write(res_file,100)
      endif
80    format ('spec wdepsol ddepsol')
90    format (1x,i1,2x,f8.3,1x,f8.3)
100   format ('no wet and dry deposition input')

c     applications of inorganic fertiliser
      call small_title(res_file,'inorganic fertilisation')
      write (res_file,110) ninof
      if (ninof.gt.1) write (res_file,120)
      do appl = 1,ninof
		call date_wave (idinof(appl),year, month,day)
		write (res_file,130)year, month, day, (fsol(appl,sp),sp=1,nr_of_sol)
      enddo
110   format ('the number of applications of inorganic fertilizer: ',
     $   1x,i2)
120   format ('year month day', 1x, 'amounts for each solute ',
     $'species 1..nr of solutes',/,'              (mg/m**2)')
130   format (i4,2x,i2, 2x, i2, 1x,5(f10.3,1x))

c     irrigations
      call small_title(res_file,'irrigation')
      write (res_file,140) nirr
      if(nirr.lt.0) write (res_file,150)
      do appl = 1,nirr
		call date_wave (idirr(appl),year, month,day)
		write (res_file,160)year,month,day,(conirsol(appl,sp),sp=
     $		1,nr_of_sol)
      enddo
140   format('the number of irrigations                       : ',1x,i2)
150   format('ir_day solute amounts for 1..nr of solutes',/,'   (mg)??')                      
160   format(i4,2x,i2, 2x, i2, 1x,5(f8.3,1x))

c     plowing
      call small_title(res_file,'plowing')
      write (res_file,170) nplo
      if (nplo.gt.0) write (res_file,180)
      do i = 1, nplo
		call date_wave (idplow(i),year, month,day)
		write (res_file,190)year, month, day ,ncsplow(i)
      enddo
170   format ('the field is plowed ',i2,' times')
180   format ('plowing dates and number of compartments affected')
190   format (i4,2x,i2, 2x, i2, 1x,i3)

c     solute : transformation rates
      call big_title(res_file,'transformation rates(first order decay)')
      write (res_file,200)
      do i = 1, ncs
		write (res_file,210)i, (rates(i,sp),sp=1,nr_of_sol)
      enddo
200   format (
     $'comp   transformation rates for solutes 1 .. nr_of_sol',/,
     $'        (day-1)' )
210   format (i3,1x,8(f8.3,1x))

c     solute : initial values
      call big_title(res_file,'initial conditions')
      write (res_file,220)
      do i = 1, ncs
		write (res_file,230)i,(tcsolo(i,sp),sp= 1,nr_of_sol)
      enddo
220   format ('comp',1x,'sol. conc. for 1..nr of sol (first three n)' 
     $,/,     '           (mg/m**2)')
230   format (i3,1x,5(f8.3,1x))

c     solute: output control
      call big_title(res_file,'output')
      call small_title(res_file,'solutes for which output is produced')
      do sp = 1, nr_of_sol
		if( out_sol(sp)) then
			if (simnit.and. sp.le.3) then
				write(res_file,240) n_sol_name(sp),plus
			else
				write(res_file,250) sp,plus
			endif
		else
			if (simnit.and. sp.le.3) then
				write(res_file,240) n_sol_name(sp),neg
			else
				write(res_file,250) sp,neg
			endif
		endif
      enddo
      call small_title(res_file,'the kind of output')
      if (pr_sol_sum) then
		write (res_file,260)plus
      else
		write (res_file,260)neg
      endif
      if(pr_sol_cum) then
		write (res_file,270)plus
      else
		write (res_file,270)neg
      endif
      if(pr_sol_tseries(1))then
		write (res_file,280)plus
      else
		write (res_file,280)neg
      endif
      if(pr_sol_tseries(3))then
		write (res_file,290)plus
      else
		write (res_file,290)neg
      endif
240   format 
     $('output for solute ',a8,    '                         : ',a)
250   format 
     $('output for solute ',i2,'                            : ',a)
260   format
     $('overview of main state variables (sol_sum.out)      : ',a)
270   format  
     $('cumulative terms of the solute balance (sol_cum.out): ',a)
280   format  
     $('solute concentration (conc*.out)                    : ',a)
290   format  
     $('sink term (sink*.out)                               : ',a)
      call close_file(res_file)
      end