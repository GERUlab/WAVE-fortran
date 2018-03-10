c	version 3.0
c	author: m vanclooster 19/02/02
c###################################################################################
      subroutine wat_boco

c     in   : arel, brel, dt, dxinter, esa, fin, gwl, ibboco, ncs, nday, pond, rirr
c     out  : dstor, dxn, flxa, ncomp, ph, prec, qdeep, rinf, stor,flxa
c     calls: calgwl, etsplit, report_err, stop_simulation
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'clim.com'
      include   'wat.com'
      common /bo_cond/flxar,flxa1,nbot,nupboco

      if (ibboco.gt.1.and.ibboco.lt.4) call calgwl  

c     calculation of a first estimate for the flux through the soil surface 
      dstor =  dmin1( dmax1(fin(nday)-stor,0.d0),(prec(nday) + rirr(nday)))
      stor = stor + dstor 
      rinf = prec(nday) + rirr(nday) - dstor
      call etsplit
      flxa1 = esa-rinf
      flxa = flxa1 - pond/dt

c     determination of the bottom boundary condition   
c     the groundwater level is input  
      if (ibboco.eq.1) then
c		determine first node above groundwater level (ncomp) and depth of first saturated node (dxn)
		dxn=0.
		i = 1
		ncomp=ncs
10		dxn=dxn-dxinter(i)
		if (dxn.lt.gwl(nday))  then
			ncomp= i - 1
		else 
			if (i.lt.ncs) then
				i = i+1
				goto 10
			endif
		endif
		if(ncomp.eq.ncs) then
			call report_err('groundwater level reached bottom of soil')
			call stop_simulation
     $		('programme stopped: check err_file')
			return
		endif
c		determine pressure head at node (ncomp+1) assuming unit gradient
		ph(ncomp+1)=gwl(nday)-dxn
		return

c     flux - groundwater level relationship  
      else if (ibboco.eq.3) then
		qdeep(nday)=arel* dexp(brel* dabs(gwl(nday)))
		return
      else
		return
      endif
      end

c###################################################################################
      subroutine calc_bboco(ncbot,phbot,flxsbot,flxs_cond,free_drainage)
c     in   : conin, gprh, ibboco, ncomp, ncs, nday, ph, qdeep
c     out  : flxs_cond, flxsbot, free_drainage, ncbot, phbot
c     calls: report_err, stop_simulation
c###################################################################################
      implicit double precision (a-h, o-z)
      include   'constant'
      include   'gen.com'
      include   'clim.com'
      include   'wat.com'
      common /bo_cond/flxar,flxa1,nbot,nupboco
      logical flxs_cond,free_drainage

c     determine number of bottom compartment(ncbot) and bottom
c     boundary condition  nbot=1 : press. head (phbot) is given
c                         nbot=2 : flux (flxsbot) is given
      if(ibboco.eq.1)then
		nbot=1
		ncbot=ncomp
		phbot=ph(ncomp+1)
      elseif((ibboco.gt.1).and.(ibboco.lt.4))then
		nbot=2
		ncbot=ncomp
		flxsbot=qdeep(nday)
      elseif(ibboco.eq.4)then
		nbot=1
		ncbot=ncs
		phbot=gprh(nday)
      elseif(ibboco.eq.5)then
		nbot=2
		ncbot=ncs
		flxsbot=0.
      elseif(ibboco.eq.6)then
		nbot=2
		ncbot=ncs
		flxsbot=-conin(ncbot+1)
      elseif(ibboco.eq.7)then
		ncbot=ncs
		if(nbot.eq.1)then
			phbot=0.
		else
			flxsbot=-0.
		endif
      endif
      flxs_cond = nbot.eq.2
      free_drainage = ibboco.eq.6
      return
      end

c###################################################################################
      subroutine check_bboco(phbot,flxsbot,flxs_cond,ok_cb)
c     in   : flxs, ibboco, ncs, ph
c     out  : flxs_cond, flxsbot, ok_cb, phbot
c     calls: -
c###################################################################################
      implicit double precision (a-h, o-z)
      include   'constant'
      include   'gen.com'
      include   'clim.com'
      include   'wat.com'
      common /bo_cond/flxar,flxa1,nbot,nupboco

      logical ok_cb,flxs_cond
      ok_cb = .true.
c     change bottom boundary condition if necessary
c     lysimeter bottom boundary condition 
      if(ibboco.eq.7)then
		if((nbot.eq.1).and.(flxs(ncs+1).gt.0.))then
			nbot = 2
			flxs_cond = .true.
			flxsbot= 0.d0
			ok_cb = .false.
		elseif((nbot.eq.2).and.(ph(ncs).gt.5.e1))then
			nbot = 1
			flxs_cond = .false.
			phbot = 0.d0
			ok_cb = .false.
		endif
      endif
      return
      end

c###################################################################################
      subroutine calc_uboco(phsurf_calc,flxsbot, flxar_out,flxs_cond)
c     in   : dt, dx, flxa, histor_file, ncs, phsa, pond, pr_wat_histor,
c            rtex, vol
c     out  : flxs_cond, flxar_out, flxsbot, phsurf_calc
c     calls: calc_wc
c###################################################################################
      implicit double precision (a-h, o-z)
      include   'constant'
      include   'gen.com'
      include   'clim.com'
      include   'wat.com'
      common /bo_cond/flxar,flxa1,nbot,nupboco
      logical flxs_cond,first_time
      data first_time/.true./

      if (first_time) then
		first_time = .false.
		volmax = 0.d0
		do i = 1, ncs 
			volmax = volmax + calc_wc(0.d0,i)*dx
		enddo
      endif
c     determine upper boundary condition 
c                   nupboco=1 : press. head at surface = phsa
c                   nupboco=2 : flux condition (flxa)
c                   nupboco=3 : press. head at surface = 0
      if(ponded.and.pond.eq.0.d0) then
		ponded = .false.
          pond_to(i_pond)= t
      else if (.not.ponded.and.pond.gt.0.d0) then
          ponded = .true.
          i_pond = i_pond + 1
          pond_from(i_pond)= t
          pond_maximum(i_pond) = 0.d0
      endif
      if (ponded) then
          if (pond_maximum(i_pond).lt.pond) pond_maximum(i_pond)=pond
      endif
      flxa = flxa1 - pond/dt
      
	if(nupboco.eq.1)then
		if(flxa.gt.0.d0)then
			phsurf_calc=phsa
		else
			nupboco=2
		endif
      elseif(nupboco.eq.3)then
		if(flxa.gt.0.d0)then
			nupboco=2
		else
			phsurf_calc=pond
		endif
      endif
      flxs_cond = nupboco.eq.2
c     avoid completely saturated profile by limitting flxa ---> flxar
      flxar=flxa
      if(flxs_cond)then
		epra=0.d0
		do i=1,ncs
			epra=epra+rtex(i)*dx
		enddo 
		excess=(flxsbot-flxa-epra)*dt-volmax+vol+1.
		if(excess.gt.0.d0)then
			flxar=flxa+excess/dt
			if(pr_wat_histor)write(histor_file,*)
     $		' soil profile nearly saturated, flxa has been limited'
		endif
      endif
      flxar_out = flxar
      return
      end

c###################################################################################
      subroutine check_uboco(phsurf_calc,flxs_cond,ok_cu)
c     in   : conduc, dxinter, flxs, histor_file, ph, phsa, pond, pr_wat_histor
c     out  : flxs_cond, ok_cu, phsurf_calc
c     calls: calc_con
c###################################################################################
      implicit double precision (a-h, o-z)
      include   'constant'
      include   'gen.com'
      include   'clim.com'
      include   'wat.com'
      common /bo_cond/flxar,flxa1,nbot,nupboco
      logical ok_cu,flxs_cond
      logical first_time
      data first_time /.true./

      if (first_time) then
		first_time = .false.
		csurfdry = calc_con(phsa,1)
		csurfsat = calc_con (0.d0,1)
      endif
c     calculate maximum infiltration and evaporation flux
      if(nupboco.eq.2)then
		conins=sqrt(csurfsat*conduc(1))
		flxsmin=-conins*( (pond-ph(1))/dxinter(1) +1.)
		conins=sqrt(csurfdry*conduc(1))
		flxsmax=-conins*( (phsa-ph(1))/dxinter(1) +1.)
      endif
      if((nupboco.eq.1).and.(flxs(1).gt.flxar))then
		nupboco=2
		ok_cu = .false.
      elseif((nupboco.eq.3).and.(flxs(1).lt.flxar))then
		nupboco=2
		ok_cu = .false.
      elseif((nupboco.eq.2).and.(flxar.gt.flxsmax))then
		nupboco=1
		phsurf_calc=phsa
		ok_cu = .false.
      elseif((nupboco.eq.2).and.(flxar.lt.flxsmin))then
		nupboco=3
		phsurf_calc=pond
		ok_cu = .false.
      else
		ok_cu = .true.
      endif
      flxs_cond = nupboco.eq.2
      if (.not.ok_cu.and.pr_wat_histor)write(histor_file,*)
     $	' nupboco has been changed -->',nupboco
      end

c###################################################################################
      subroutine fix_uboco (phsurf_calc)
c     in   : flxs, histor_file, phsa, pond, pr_wat_histor, t
c     out  : phsurf_calc
c     calls: -
c###################################################################################
      implicit double precision (a-h, o-z)
      include   'constant'
      include   'gen.com'
      include   'clim.com'
      include   'wat.com'
      common /bo_cond/flxar,flxa1,nbot,nupboco

c     in case of two unsuccessful trials for top boundary condition :
c     switch to dirichlet condition, and repeat loop without checking
      if(flxar.gt.0.d0)then
		nupboco=1
		phsurf_calc=phsa
      else
		nupboco=3
		phsurf_calc=pond
      endif
      if(pr_wat_histor)then
		write(histor_file,*)
     $	'two loops for top b.c. not successful,',
     $	'nupboco has been changed -->',nupboco
		write(histor_file,10)nupboco,flxs(1),flxar,t
10		format(
     $	'*******************warning***************************',/,
     $	'nupboco fixed at: ',i1,' for 3rd iteration for top b.c.',/,
     $	'flxs(1) = ',f7.2,' flxa = ',f7.2,' time = ',f9.4,/,
     $	'*****************************************************')
      endif
      return
      end


c###################################################################################
      subroutine calgwl
c     in   : dxinter, ncs, nday, ph, phsurf
c     out  : gwl, ncomp
c     calls: report_err, stop_simulation
c###################################################################################
      implicit double precision (a-h,o-z)
      include 'constant'
      include 'gen.com'
      include 'clim.com'
      include 'wat.com'

c     locate first unsaturated compartment
      ncomp=ncs
10    if(ph(ncomp).ge.0.d0)then
		ncomp = ncomp - 1
		if (ncomp.ge.1) goto 10
      endif
      if(ncomp.eq.ncs) then
		call report_err ('ground water level reached bottom of profile')
		call stop_simulation
     $		('programme stopped: check err_file')
		return
      endif
c     if first node is saturated  
      if(ncomp.eq.0)then
		if(phsurf.gt.0.d0)then
			gwl(nday)=phsurf
		else
			gwl(nday)=dxinter(1)/(ph(1)-phsurf)*(phsurf)
		endif
c     else calculate depth of first unsaturated node and distance from 
c     that node to the groundwater level   
      else
		dunsat=0.d0
		do i = 1,ncomp
			dunsat=dunsat+dxinter(i)
		enddo
		if (ph(ncomp+1).ne.ph(ncomp)) then
			dist=dxinter(ncomp+1)/(ph(ncomp+1)-ph(ncomp))*(-ph(ncomp))
		else
			dist = 0.d0
		endif
		gwl(nday)=-(dunsat+dist)
      endif           
      return
      end

c###################################################################################
      block data watboco
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      common /bo_cond/flxar,flxa1,nbot,nupboco
      data nbot,nupboco / 2,2/
      end