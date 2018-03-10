c	version 3.0
c	author: m. vanclooster 21/02/2002
c###################################################################################
      subroutine rer
c     in   : drz, dt, dx, epa_min_intc, ilinr, isucr, ncs, nday, p0, p1,
c            p2h, p2l, p3, ph, rnam, rt_distr, simplant, t, tb, te, wc, x
c     out  : irz, rtex
c     calls: calc_wc, root_sucros
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'clim.com'
      include   'wat.com'
      dimension wcwp(kt_comps)
      dimension alpha_red(kt_comps)
      dimension rtex1(kt_comps)
      dimension fr_pot_et(kt_comps)
      double precision rest

      if (.not.simplant) return
      if (isucr)  drz(nday) =  root_sucros()
      if( dint(t)-t.eq.0.d0) call calc_dim_roots(ira,irz,rna)
      if( dint(t)-t.eq.0.d0) call calc_rt_distr(fr_pot_et)
      do i = 1,irz
		rtex1(i) = fr_pot_et(i)* epa_min_intc/dx
      enddo
      if (ira.gt.ncs) then
		do i = 1, ncs
			rtex(i) = 0.d0
		enddo
      else
		rest = 0
		do i = 1, ira - 1
			rest = rest + rtex1(i)
			rtex(i) = 0.d0
		enddo
		p2=p2h
		if(epa_min_intc.lt.0.1d0) then
			p2=p2l
		else
			if(epa_min_intc.lt.0.5d0)
     $			p2=p2h+((0.5d0-epa_min_intc)/0.4d0)*(p2l-p2h)
		endif
		do i=ira,irz
			if(ph(i).gt.p1) then
				if (p0.gt.p1) then
					alpha_red(i)=dmax1((p0-ph(i))/(p0-p1),0.0d0)
				else
					alpha_red(i) = 0.d0
				endif
			endif
			if(ph(i).lt.p2) then
				if (ph(i).gt.p3) then
					if (p3.lt.p2) then
						if(ilinr) then
							alpha_red(i)=dmax1((p3-ph(i))/(p3-p2),0.0d0)
						else
							alpha_red(i)=10**((p2-ph(i))/p3)
						endif
					else
						alpha_red(i) = 0.d0
					endif
				else
					alpha_red(i) = 0.d0
				endif
			else
				alpha_red(i) = 1.d0
			endif
		enddo

c		root extraction rate is also limited to avoid extraction below
c		wilting point (wcwp)
		do i = ira, irz
			wcwp(i)=calc_wc(p3,i)
		enddo
		do i = ira, irz
			rtex(i)= dmax1(0.d0,
     $                 dmin1(alpha_red(i)*(rtex1(i)+rest),
     $                         (wc(i)-wcwp(i))/dt) )
			rest = dmax1(rest-rtex(i) + rtex1(i)-rtex(i),0.d0)
		enddo

c		rtex(i)= 0 for compartments below the root zone
		do i = irz+1, ncs
			rtex(i) = 0.d0
		enddo
      endif
      return
      end

c###################################################################################
      subroutine etsplit
c     in   : et0, harvest_date,
c            idvs, isucr, nday, plant_date, rlai, simplant, t, 
c     out  : akc, dvs_act, epa, epa_min_intc, esa, ev_intc, stor
c     calls: dvs_sucros, rlai_sucros
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      include   'clim.com'
      data kc_num /1/
	
c	calculation of kc factor
      if (idvs) then
		dvs_act = dvs_sucros()
		if (kc_num.eq.1. and. dvs_act.le.dvs_kc(1)) then
			akc = kc(kc_num)
		else if (kc_num.eq.nr_of_kc_values. and. 
     $                      dvs_act.ge.dvs_kc(kc_num)) then
			akc = kc(kc_num)
		else
			if (dvs_act.gt.dvs_kc(kc_num)) kc_num = kc_num + 1
			akc = kc(kc_num - 1) + 
     $			(kc(kc_num)     -     kc(kc_num-1)) * 
     $			(dvs_act        - dvs_kc(kc_num-1))   /
     $			(dvs_kc(kc_num) - dvs_kc(kc_num-1))
		endif
	else
		if (kc_num.eq.1. and. t.le.dfloat(iday_kc(1))) then
			akc = kc(kc_num)
		else if (kc_num.eq.nr_of_kc_values. and. 
     $                 t.ge.dfloat(iday_kc(kc_num))) then
			akc = kc(kc_num)
		else
			if (t.gt.dfloat(iday_kc(kc_num))) kc_num = kc_num + 1
			akc = kc(kc_num - 1) + 
     $			(kc(kc_num)     -     kc(kc_num-1)) *
     $			(t              - dfloat(iday_kc(kc_num-1))) /
     $			(dfloat(iday_kc(kc_num))-dfloat(iday_kc(kc_num-1)))
		endif
      endif

c	splitting transpiration-evaporation 
	if(isucr) rlai(nday) = rlai_sucros()
      et0a=akc*et0(nday)
      esa= dexp(-0.6d0* rlai(nday))*et0a
      epa=et0a-esa

c     evaporation of intercepted water
      ev_intc =  dmin1 (stor , epa)
      epa_min_intc = epa - ev_intc
      stor = stor - ev_intc
      end

c###################################################################################
      double precision function calc_wat_drz(day)
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'wat.com'
      integer day

      calc_wat_drz = drz(day)
      return
      end

c###################################################################################
      double precision function calc_rtex(nr_of_comp)
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'wat.com'

      calc_rtex = rtex(nr_of_comp)
      return
      end

c###################################################################################
      double precision function calc_watreductgrow ()
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'

      if (idint(t).lt.plant_date.or.idint(t).gt.harvest_date) then
		calc_watreductgrow=0.d0
		return
      endif
      if (epa.eq.0.d0) then
		calc_watreductgrow = 0.d0
      else
		calc_watreductgrow =  dmin1(((transp_act+ev_intc)/epa),1.d0)
      endif
      end


c###################################################################################
      subroutine calc_dim_roots(ira,irz,rna)
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
c     ira =  compartment for which the bottom is below the first compartment
c            in which there is extraction
c     irz =  compartment for which the lower boundary is less or equal to drz
c     rna = depth above which roots are inactive

      if(t.ge.tb) then
		if (te.eq.tb) then
			rna = dmax1(rnam,drz(nday))
		else
			rna= dmax1(rnam* dmin1((t-tb)/(te-tb),1.d0),drz(nday))
		endif
      else
		rna=0.0
      endif

      ira = 1
10    if((x(ira)-dx/2.).ge.rna)then
		ira=ira + 1
		if (ira.le.ncs) goto 10
      endif

      irz = ira
      if (ira.le.ncs) then
c		determine compartment for which the lower boundary is
c		less or equal to drz (irz)
20		if((x(irz)-dx/2.).gt.drz(nday))then
			irz = irz +1
			if (irz .le.ncs) goto 20
		endif
      endif
      irz = min(irz, ncs)
      return
      end

c###################################################################################
      subroutine calc_rt_distr(fractions)
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      dimension fractions(*)

      call calc_dim_roots(ira,irz,rna)
      sum_of_rt_distr = 0.d0
      if (ira.ne.irz) then
		sum_of_rt_distr = rt_distr(ira)* (rna - (x(ira)-dx/2.))/dx
		do i = ira+1,irz -1
			sum_of_rt_distr = sum_of_rt_distr + rt_distr(i)
		enddo
c		for last compartment correct for partial uptake
		sum_of_rt_distr = sum_of_rt_distr+rt_distr(irz)*
     $                    ((dx/2.)+x(irz)-drz(nday))/dx
      else
		sum_of_rt_distr = (rna - drz(nday))/dx
      endif

      if (sum_of_rt_distr.gt.0.d0) then
		do i = 1, ira-1
			fractions(i) = 0.d0
		enddo
		do i = ira,irz
			fractions(i) = rt_distr(i)/sum_of_rt_distr
		enddo
c		for first & last compartment correct for partial uptake
		if (ira.ne.irz) then
			fractions(ira) = (rt_distr(ira)*(rna - (x(ira)-dx/2.))/dx)
     $         /sum_of_rt_distr
			fractions(irz) = (rt_distr(irz)*((dx/2.)+x(irz)-drz(nday))
     $                    /dx)/sum_of_rt_distr
		else
			fractions(ira) = 1.d0
		endif
      else
		do i = 1, irz
			fractions(i) = 0.d0
		enddo
      endif
      return
      end