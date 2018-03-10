c	Version 3.0
c	Author: m. vanclooster 15/02/2002
c###################################################################################
      subroutine sol_boco
c     in   : acsolio, acsolmo, conirsol, csolio, csolmo, csolo, ddepsol,
c            dx, fsol, idinof, idirr, idplow, imobsw, ncsplow, nday,
c            ninof, nirr, nla, nplo, nr_of_sol, parasol, prec, solinfl,
c            t, tcsolio, tcsolmo, tcsolo, wcio, wcmo, wco, wdepsol
c     out  : acsolib, acsolio, acsolmb, acsolmo, acsolo, csolb, csolib,
c            csolio, csolmb, csolmo, csolo, ressol, solsur, tcsolio,
c            tcsolmo, tcsolo
c     calls: -
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'clim.com'
      include   'sol.com'
      logical application_passed, application_day 
      integer dt_sol_count, sp 
      common/dt_sol_boco/application_day,application_passed,dt_sol_count
      data iteliof,itelplow,itelirr/1,1,1/

      if(application_day)application_passed=.true.
      application_day=.false.
      volume=dx
      do sp = 1,nr_of_sol
c		wet deposition
		ressol(sp) = ressol(sp)+prec(nday)*wdepsol(sp)
c		dry deposition
		ressol(sp) = ressol(sp) + ddepsol(sp)
c		irrigation
		if (nirr.gt.0) then
			if(idint(t) .eq. (idirr(itelirr))) then
				ressol(sp) = ressol(sp)+conirsol(itelirr,sp)
				if(sp.eq.nr_of_sol) then 
					itelirr = itelirr+1
					application_day=.true.
				endif
			endif
		endif
c		inorganic fertilizer
		if (ninof.gt.0) then
			if(idint(t) .eq. (idinof(iteliof))) then
				ressol(sp) = ressol(sp)+fsol(iteliof,sp)
				if(sp.eq.nr_of_sol) then 
					iteliof=iteliof+1
					application_day=.true.
				endif
			endif
		endif
c		depletion during last day
		ressol(sp)= dmax1(ressol(sp),mn_real)-dmax1(solinfl(sp),mn_real)
c		plowing
		if(nplo.gt.0) then
			if(idint(t) .gt. idplow(itelplow)) then
				itelplow = itelplow +1
				if (itelplow.gt.nplo) itelplow = nplo          
			endif
			if(idint(t) .eq. idplow(itelplow)) then
c			calculate the total amount of sol in reservoir and plow layer
c			calculate total amount of sp in profile   
				pltsol  = 0.d0
				do i=1,ncsplow(itelplow)
					pltsol = pltsol+tcsolo(i,sp)
				enddo 
				pltsol=pltsol+ dmax1(ressol(sp),0.0d0)
				plo_in(sp) = plo_in(sp) + ressol(sp)
				ressol(sp)   = 0.d0
c			recalculate concentrations in plow layer
				do i=1,ncsplow(itelplow)
					csolb(i,sp)=csolo(i,sp)
					csolmb(i,sp)=csolmo(i,sp)
					csolib(i,sp)=csolio(i,sp)
					acsolmb(i,sp)=acsolmo(i,sp)
					acsolib(i,sp)=acsolio(i,sp)
					tcsolo(i,sp) = pltsol/ncsplow(itelplow)
					tcsolio(i,sp)=tcsolo(i,sp)*wcio(i)/wco(i)
					tcsolmo(i,sp)=tcsolo(i,sp)*wcmo(i)/wco(i)
					if(imobsw) then
						csolio(i,sp)=tcsolio(i,sp)/
     $					(volume*(parasol(2,nla(i),sp)*
     $					(1-parasol(9,nla(i),sp)) *
     $					parasol(1,nla(i),sp)+wcio(i)))
					else
						csolio(i,sp) = 0.0d0
					endif
					csolmo(i,sp)=tcsolmo(i,sp)/(volume*(parasol(2,nla(i),sp)*
     $				parasol(9,nla(i),sp) *parasol(1,nla(i),sp)+wcmo(i)))
					csolo(i,sp)=(csolio(i,sp)*wcio(i)+csolmo(i,sp)*wcmo(i))/
     $				wco(i)
					if(imobsw) then
						acsolio(i,sp)=parasol(2,nla(i),sp)*
     $					(1-parasol(9,nla(i),sp))*csolio(i,sp)
					else
						acsolio(i,sp) = 0.0d0
					endif
					acsolmo(i,sp)=parasol(2,nla(i),sp)*
     $				parasol(9,nla(i),sp)*csolmo(i,sp)
					acsolo(i,sp)=acsolio(i,sp)+acsolmo(i,sp)
				enddo
			endif
		endif
	enddo 
c     calculate solsur
      do sp = 1,nr_of_sol
		if(prec(nday).gt.0.d0.or.rirr(nday).gt.0.d0)then
			solsur(sp)  = ressol(sp)
		else
			solsur(sp) = 0.d0
		endif
	enddo 
      return
      end	 

c###################################################################################
      double precision function dt_sol(dtmax_in)
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      logical application_day,application_passed
      integer dt_sol_count, sp
      common/dt_sol_boco/application_day,application_passed,dt_sol_count
 
      dt_sol = dtmax_in
      do sp = 1, nr_of_sol
		if(solsur(sp).gt.0.d0)dt_sol= dmin1(dt_sol,0.0099d0)
      enddo
      if(application_passed)then
		dt_sol =dtmin
		dt_sol_count=dt_sol_count+1
		if(dt_sol_count.eq.10) then
			application_passed=.false.
			dt_sol_count=0
		endif
      endif
      return
      end