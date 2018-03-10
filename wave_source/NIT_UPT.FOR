c	version 3.0
c	author:	m.vanclooster (26/02/2002) 
c###################################################################################
      subroutine nit_upt
c     in   : csolo, dt, dx, g, harvest_date, isucr, ncs, nla, parasol,
c            plant_date, rdo, rnmaxp, rorad, simplant, t
c     out  : rdens_left, rfract_tot, suni, sunm
c     calls: calc_nitreductgrow, calc_rdens, calc_rtex, calc_xncle, plant_weights, 
c            rdens_sucros, report_err, stop_simulation
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      include   'nit.com'
      common /local_nit_upt/ ancl, rlncl,rmncl
      integer sp
      logical first_time
      dimension hcsolo(kt_sol)
      dimension rdens(kt_comps)
      dimension rcsolo(kt_sol),diffus_rm(kt_comps,kt_sol)
      parameter ( pi=3.14159 )
      data rlncl/0.005/
      data rcsolo(3)/6.2d0/
      data rcsolo(2)/1.8d0/
      data anst/mn_real/
      data first_time/.true./

      volume=dx
      if(.not.simplant .or.
     $	idint(t).lt.plant_date.or.idint(t).ge.harvest_date) then
		do i=1,ncs
			do sp = 2,3
				sunm(i,sp)=0.0d0
				suni(i,sp)=0.0d0
			enddo
		enddo
		return
      else
		if (first_time) then
			first_time = .false.
			if (isucr) then
				anlv = sucr_anlv()
			else
				totdem = rnmaxp/100.d0
				g = g * (dfloat(harvest_date - plant_date)+1.d0)
				factor = rnmaxp/(((g/2.d0) - (g/3.d0))*(g**2))
			endif
		endif
      endif
      if (isucr) then
		call rdens_sucros(rdens)
      else
		call calc_rdens(rdens)
      endif
      if ( dint(t).eq.t.and.idint(t).eq.harvest_date-1) then
		rfract_tot = 0.d0
		do i= 1,ncs
			rfract_tot=rfract_tot+rdens(i)
			rdens_left(i) = rdens(i)
		enddo
		rfract_tot = dmax1(rfract_tot,mn_real)
      endif
      if (isucr) then
		xncle= calc_xncle()
		xncst=calc_xncst()
		xncrt=calc_xncrt()
		rmncl=0.5d0*xncle
		if (rmncl.eq.0.d0) then
			call report_err ('rmncl is zero')
			call stop_simulation ('programme stopped: check err_file')
		endif
		call plant_weights(wso,wlv,wrt,wst,wlvg)
		ancl=anlv/wlvg
		red = calc_nitreductgrow()
		fndef=1-sqrt(1-red*red)
		rndeml= dmax1(0.d0,(wlv*xncle-anlv))
		rndemst= dmax1(0.d0,(wst*xncst-anst))
		rndemrt= dmax1(0.d0,(wrt*xncrt-anrt))
		xncso=calc_xncso()
		rndemso= dmin1( dmax1(0.d0,(wso*xncso-anso))*fndef,
     $          rndeml+rndemst)
		totdem= dmax1(0.d0,rndeml+rndemst+rndemrt)*100.0d0
		if (totdem.eq.0.d0) then
			do sp = 2,3
				do i = 1, ncs 
					sunm(i,sp)=0.0d0
					suni(i,sp)=0.0d0
				enddo
			enddo
			return
		endif
      else
		if((rnmaxp-totup).lt.0.1d0) then
			do i=1,ncs
				do sp = 2,3
					sunm(i,sp)=0.0d0
					suni(i,sp)=0.0d0
				enddo
			enddo
			return
		endif
		time_since_pl = t-dfloat(plant_date)+1.d0
		totdem= dmax1(factor*time_since_pl*(g-time_since_pl),0.d0)
     $		   * dt
		if (rnmaxp .lt. totup + totdem ) totdem = rnmaxp - totup
      endif
      tconn=0.d0
      do sp=2,3
		do i=1,ncs
			hcsolo(sp) = dmax1(0.d0,(csolo(i,sp)-rcsolo(sp)))
			unc(i,sp) = calc_rtex(i)*hcsolo(sp)*wco(i)*dt*dx
			tconn=tconn+unc(i,sp)
		enddo
	enddo
      if (tconn.gt.totdem)then
		do sp=2,3
			do j=1, ncs
				unc(j,sp)=unc(j,sp)*totdem/tconn
				und(j,sp) =0.d0
			enddo
		enddo
		tconn=totdem
		tund=0.d0
      else
		tpdnup=totdem-tconn
		tund=0.d0
		do sp=2,3
			do i=1,ncs
				hcsolo(sp) = dmax1(0.d0,(csolo(i,sp)-rcsolo(sp)))
				diffus_rm(i,sp) = parasol(4,nla(i),sp)*
     $				dexp(parasol(5,nla(i),sp)*wco(i))*parasol(3,nla(i),sp)/
     $				wco(i)
				und(i,sp)=rdens(i)*rorad*2*pi*diffus_rm(i,sp)*
     $				hcsolo(sp)*wco(i)*dx*dt/(rdo*100000.d0)
				tund=tund+und(i,sp)
			enddo
		enddo
		if (tund.gt.tpdnup) then
			do sp = 2,3
				do j=1, ncs
					und(j,sp)=und(j,sp)*tpdnup/tund
				enddo
			enddo
			tund=tpdnup
		endif
      endif
c     calculate the total uptake (mg/m**2)
      totvr=(tconn+tund)
      totup=totup+totvr
      if(isucr) then
		rnupl=totvr*rndeml/(totdem+0.000000001d0)
		rnupst=totvr*rndemst/(totdem+0.0000000001d0)
		rnuprt=totvr*rndemrt/(totdem+0.0000000001d0)
c		calculate new concentrations in the different parts
		anlv=anlv+(rnupl-rndemso*anlv/(anlv+anst))
		anst=anst+(rnupst-rndemso*anst/(anlv+anst))
		anrt=anrt+rnuprt
		anso=anso+rndemso
      endif
c     calculate the uptake rates
      do sp=2,3
		do i=1,ncs
			sunm(i,sp)=(unc(i,sp)+und(i,sp))*parasol(7,nla(i),sp)/(dx*dt)
			suni(i,sp)=(unc(i,sp)+und(i,sp))*(1-parasol(7,nla(i),sp))/
     $		(dx*dt)
		enddo
	enddo
      return
      end  

c###################################################################################
      subroutine calc_rdens(rdens)
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'gen.com'
      include  'nit.com'
      dimension rdens(kt_comps)
      dimension fractions(kt_comps)

      call calc_rt_distr(fractions)
      i = 0
10    i = i + 1
      if(fractions(i) .eq.0.) goto 10
	rdens(i) = w0_rdens
      do j = i+1, ncs
		rdens(j)=w0_rdens* fractions(j)/fractions(i)
      enddo
	return
      end