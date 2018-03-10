c	version 3.0
c	author: m vanclooster 19/02/02
c###################################################################################
      subroutine wat_intgr
c     in   : dstor, dt, dx, epa, esa, ev_intc, flxs, irz, isd, ncs, nday, 
c            prec, rinf, rirr, rtex, t, transp_act, visd, visdi, vol, voli, 
c            wc
c     out  : cberr_wat, cberr_wat_isd, cev_intc, cflbu, cflbup, cflsd, 
c            cflsdp, cinf, cirr, cprec, cpsev, cptra, crtex, csev, ctra, 
c            ctra1, cum_intc, cwatupt_isd
c     calls: -
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'clim.com'
      include   'wat.com'

      if( dint(t)-t.eq.0.d0) transp_act=0.0d0

c     profile integration 
      cprec=cprec+prec(nday)*dt
      cirr =cirr +rirr(nday)*dt
      cpsev=cpsev+esa*dt
      cptra=cptra+epa*dt
      cev_intc = cev_intc +ev_intc*dt
      cum_intc = cum_intc + dstor *dt
c     actual infiltration and evaporation  
      if(flxs(1).gt.-1.0e-10) then
		csev=csev+flxs(1)*dt
		if(rinf.gt.1.0e-10) then
			csev=csev+rinf*dt
			cinf=cinf+rinf*dt
		endif
      else
		cinf=cinf-flxs(1)*dt+esa*dt
		csev=csev+esa*dt
      endif
c     actual transpiration  
      ctra1=ctra
      do i=1,ncs
		ctra=ctra+rtex(i)*dx*dt
		crtex(i)=crtex(i)+rtex(i)*dx*dt
      enddo
      transp_act=transp_act+ctra-ctra1
c     drainage 
      cflbup=cflbup+ dmax1(0.0d0,flxs(ncs+1))*dt
      cflbu=cflbu+flxs(ncs+1)*dt
      vol=0.0
      do i=1,ncs
		vol=vol+wc(i)*dx
      enddo
c     balance error
      cberr_wat = (vol - voli) -(cflbu + cinf - csev - ctra)

c	subprofile (isd) integration 
c	water uptake
      do i = 1,isd
		cwatupt_isd = cwatupt_isd + rtex(i)*dx*dt
	enddo
c	drainage 
      cflsdp=cflsdp+ dmax1(0.0d0,flxs(isd+1))*dt
      cflsd=cflsd+flxs(isd+1)*dt
c     storage 
      visd =0.d0
      do i=1,isd
		visd=visd+wc(i)*dx
      enddo
c	balance error 
      cberr_wat_isd = (visd -visdi)-(cflsd+cinf-csev-cwatupt_isd)
      return
      end