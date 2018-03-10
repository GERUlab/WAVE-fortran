c	Version 3.0
c	author: M. vanclooster 25/02/2002
c###################################################################################
      subroutine nit_sink
c     in   : csolio, csolmo, decsoli, decsolm, dt, imobsw, ncs, nla,
c            nr_of_sol, parasol, wcio, wcmo
c     out  : decsoli, decsolm, sinki, sinkm
c     calls: nit_upt, seqtrans, sol_sink
c###################################################################################
      implicit double precision (a-h,o-z)
      call nit_upt
      call seqtrans 
      call nit_correct_sink
      end

      subroutine nit_correct_sink
c###################################################################################
c     in   : csolio, csolmo, decsoli, decsolm, dt, imobsw, ncs, nla,
c            nr_of_sol, parasol, wcio, wcmo, rhuream,rnitrim,sunm,
c            rhureai,rnitrii,suni
c     out  : sinki, sinkm,rhuream,rnitrim,sunm,rhureai,rnitrii,suni
c     calls:
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      include   'nit.com'
      integer sp

      do sp = 1, nr_of_sol
c		the mobile sinkterm
		do i=1,ncs-1
			sinkm(i,sp)=decsolm(i,sp)
     $			-parasol(8,nla(i),sp)*(csolmo(i,sp)-csolio(i,sp))
			supply=csolmo(i,sp)*(wcmo(i)+parasol(2,nla(i),sp)*
     $		parasol(1,nla(i),sp)*parasol(9,nla(i),sp))/dt
			if(sinkm(i,sp).lt.0.0d0) then
				demand=-sinkm(i,sp)
				if(demand.gt.supply) then
					deficit=supply-demand
					sinkm(i,sp)= dmin1(sinkm(i,sp)-deficit,0.d0)
					if(sp.eq.1) rhuream(i)= dmax1(rhuream(i)+deficit,0.0d0)
					if(sp.eq.2) rnitrim(i)= dmax1(rnitrim(i)+deficit,0.0d0)
					if(sp.eq.3) sunm(i,3)= dmax1(sunm(i,3)+deficit,0.0d0)
				endif
			endif
c		the immobile sink term
			if(imobsw) then
				sinki(i,sp)=decsoli(i,sp)
				supply=csolio(i,sp)*(wcio(i)+parasol(2,nla(i),sp)*
     $				parasol(1,nla(i),sp)* (1-parasol(9,nla(i),sp)))/dt
				if(sinki(i,sp).lt.0.0d0) then
					demand=-sinki(i,sp)
					if(demand.gt.supply) then
						deficit=supply-demand
						sinki(i,sp)= dmin1(sinki(i,sp)-deficit,0.d0)
						if(sp.eq.1) rhureai(i)=dmax1(rhureai(i)+deficit,0.0d0)
						if(sp.eq.2) rnitrii(i)=dmax1(rnitrii(i)+deficit,0.0d0)
						if(sp.eq.3) suni(i,3)=dmax1(suni(i,3)+deficit,0.0d0)
					endif
				endif
			endif
		enddo
	enddo
      return
      end

c###################################################################################
      subroutine seqtrans
c     in   : ncs, rcarbmin, rdeniti, rdenitm, rhureai, rhuream, rmini,
c            rminm, rnitrii, rnitrim, rvoli, rvolm, suni, sunm 
c     out  : deccorg, decnorg, decsoli, decsolm,
c     calls: denitrification, hydrolysis, miner_immob, nitrification, rate_reduction,
c            sol_sink, volatilisation
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      include   'nit.com'
      common/local_nit_sink/
     $		reddenit(kt_comps), redmoi(kt_comps), redtemp(kt_comps)
      dimension temp(kt_comps)

      call what_is_soiltemp(temp)
      do i= 1,ncs
c		calculation of the temperature reduction factor
		redtemp(i)=3**((temp(i)-16.0d0)/10.)
c		calculation of the different moisture threshold values
		thetaw=calc_wc(-150000.d0,i)
		thetas=calc_wc(0.0d0,i)
		thetah=calc_wc(-100.d0,i)
		thetal=calc_wc(-1000.d0,i)
c		calculation of the moisture reduction factors
		if(thetas.gt.wco(i).and.wco(i).ge.thetah) then
			redmoi(i)=0.5+ 0.5*((thetas-wco(i))/(thetas-thetah))
		elseif(thetah.gt.wco(i).and.wco(i).ge.thetal) then
			redmoi(i)=1.d0
		elseif(thetal.gt.wco(i).and.wco(i).ge.thetaw) then
			redmoi(i)=((wco(i)-thetaw)/(thetal-thetaw))	
		elseif(wco(i).ge.thetas) then
			redmoi(i)=0.5d0
		elseif(wco(i).lt.thetaw) then
			redmoi(i)=0.0d0
		endif
		if (thetas .eq.0.d0) then
			call report_err ('saturated moisture content is zero')
			call stop_simulation ('programme stopped : check err_file')
		endif
		rsat=wco(i)/thetas
		if ((rsat-0.80d0).lt.0.0001) then
			reddenit(i)=0.d0
		else
			reddenit(i)=((rsat-0.80d0)/.20d0)**2
		endif
      enddo
      do i=1,ncs-1
		rkhyd_act=redtemp(i)*redmoi(i)*rkhyd(i)
		rhuream(i)=rkhyd_act*csolmo(i,1)*wco(i)*parasol(7,nla(i),1)
		rhureai(i)=rkhyd_act*csolio(i,1)*wco(i)*(1-parasol(7,nla(i),1))
      enddo
      do i=1,ncs-1
		rknitri_act=redtemp(i)*redmoi(i)*rknitri(i)
		rnitrim(i)=rknitri_act*csolmo(i,2)*wco(i)*parasol(7,nla(i),2)
		rnitrii(i)=rknitri_act*csolio(i,2)*wco(i)*
     $		(1-parasol(7,nla(i),2))
      enddo
      do i=1,ncs-1
		rkdenit_act=redtemp(i)*reddenit(i)*rkdenit(i)
		rdenitm(i)=rkdenit_act*csolmo(i,3)*wco(i)*parasol(7,nla(i),3)
		rdeniti(i)=rkdenit_act*csolio(i,3)*wco(i)*
     $		(1-parasol(7,nla(i),3))
      enddo
      do  i=1,ncs-1
		rkvol_act=redtemp(i)*redmoi(i)*rkvol(i)
		rvolm(i)=rkvol_act*csolmo(i,2)*wco(i)*parasol(7,nla(i),2)
		rvoli(i)=rkvol_act*csolio(i,2)*wco(i)*(1-parasol(7,nla(i),2))
      enddo
      call miner_immob
c     new values for decsol* 
c	urea
      do i=1,ncs-1
		decsolm(i,1)=decsolm(i,1)-rhuream(i)
		decsoli(i,1)=decsoli(i,1)-rhureai(i)
      enddo
c	ammonia
      do i=1,ncs-1
		decsolm(i,2)=decsolm(i,2)+rhuream(i)-sunm(i,2)-rnitrim(i)+
     $		rminm(i,2)-rvolm(i)
		decsoli(i,2)=decsoli(i,2)+rhureai(i)-suni(i,2)-rnitrii(i)+
     $		rmini(i,2)-rvoli(i)
      enddo
c     nitrate
      do i=1,ncs-1
		decsolm(i,3)= decsolm(i,3)-sunm(i,3)+rnitrim(i)-rdenitm(i)+
     $		rminm(i,3)
		decsoli(i,3)= decsoli(i,3)-suni(i,3)+rnitrii(i)-rdeniti(i)+
     $	rmini(i,3)
      enddo
c     organic nitrogen
      do i = 1,ncs-1
		decnorg(i)=-rminm(i,2)-rminm(i,3)-rmini(i,2)-rmini(i,3)
      enddo
c     organic carbon
      do i = 1,ncs-1
		deccorg(i)=-rcarbmin(i)
      enddo
      end


c###################################################################################
      subroutine miner_immob  
c     in   : cchumo, cclito, ccmano, cnhumo, cnlito, cnmano, csolo, dt,
c            dx, fe, fh, ncs, rkhum, rklit, rkman, ro, wco
c     out  : cchumo, cclito, ccmano, cnhumo, cnlito, cnmano, rcarbmin,
c            rmini, rminm
c     calls: report_err, stop_simulation
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      include   'nit.com'
      common/local_nit_sink/
     $reddenit(kt_comps), redmoi(kt_comps), redtemp(kt_comps)
      dimension cco2o(kt_comps)
      data cco2o/kt_comps*0.d0/

      volume=dx
c     pxl= ratio defining mineralisation or immobilisation for litter
c     pxm= ratio defining mineralisation or immobilisation for liiter
c          > 0 = mineralisation
c          < 0 = immobilisation
      do i=1,ncs
		red_fact = redmoi(i) * redtemp(i)
		pxl=(cnlito(i)/cclito(i)-(fe/ro))
		pxm=(cnmano(i)/ccmano(i)-(fe/ro))
		eff_rate_lit = rklit(i)*red_fact
		eff_rate_man = rkman(i)*red_fact
		eff_rate_hum = rkhum(i)*red_fact
		check_mim = ( eff_rate_hum*cnhumo(i)+
     $		pxm*eff_rate_man*ccmano(i)+
     $		pxl*eff_rate_lit*cclito(i))*dt
		if(check_mim .ge. 0.d0) then
			change_nh4 = check_mim
			change_no3 = 0.d0
			eff_rate_lit = rklit(i)*red_fact
			eff_rate_man = rkman(i)*red_fact
		else
			if((csolo(i,2)-1.8d0)*wco(i)*dx .ge. (-check_mim)) then
				change_nh4 = check_mim
				change_no3 = 0.d0
				eff_rate_lit = rklit(i)*red_fact
				eff_rate_man = rkman(i)*red_fact
			else
				change_nh4 = -(csolo(i,2)-1.8d0)*wco(i)*dx
				chno3  = (csolo(i,2)-1.8d0)*wco(i)*dx+check_mim
				if((csolo(i,3)-6.2d0)*wco(i)*dx .ge. (-chno3)) then
					change_no3 = chno3
					eff_rate_lit = rklit(i)*red_fact
					eff_rate_man = rkman(i)*red_fact
				else
					change_no3 = -(csolo(i,3)-6.2d0)*wco(i)*dx
					shortage_min_n   = (csolo(i,3)-6.2d0)*wco(i)*dx+chno3
					if(pxm .lt. 0.d0  .and. pxl .ge. 0.d0) then
						vimmb    = pxm*eff_rate_man*ccmano(i)*dt
						vimma    = vimmb-shortage_min_n   
						eff_rate_lit = rklit(i)*red_fact
						eff_rate_man = vimma/(pxm*ccmano(i)*dt)
					elseif(pxl .lt. 0.d0  .and.pxm .ge. 0.d0) then
						vimmb    = pxl*eff_rate_lit*cclito(i)*dt
						vimma    = vimmb-shortage_min_n   
						eff_rate_lit = vimma/(pxl*cclito(i)*dt)
						eff_rate_man = rkman(i)*red_fact
					elseif(pxl .lt. 0.d0 .and. pxm .lt. 0.d0) then
						vimmb    = pxl*eff_rate_lit*cclito(i)*dt
						vimma    = vimmb-shortage_min_n   
						if(vimma .gt. 0.d0) then
							vimmb    = pxm*eff_rate_man*ccmano(i)*dt
							vimma    = vimmb+vimma
							eff_rate_lit = 0.d0
							eff_rate_man = vimma/(pxm*ccmano(i)*dt)
						else
							call report_err('troubles in rminimm')
							call stop_simulation 
     $						('programme stopped : check err_file')
						endif
					else
						call report_err('troubles in rminimm')
						call stop_simulation 
     $						('programme stopped : check err_file')
					endif
				endif
			endif
		endif
c		the gain/losses of the nitrogen litter pool
c		through mineralisation/immobilisation and nitrogen humification
		cnlito(i)=cnlito(i)+(-pxl*eff_rate_lit-fe*fh/ro*
     $		eff_rate_lit)*cclito(i)*dt
c		the gain/losses of the nitrogen litter pool
c		through mineralisation/immobilisation and nitrogen humification
		cnmano(i)=cnmano(i)+(-pxm*eff_rate_man-fe*fh/ro*
     $		eff_rate_man)*ccmano(i)*dt
c		the production of co2
		cco2o(i)=cco2o(i)+(1-fe)*(eff_rate_man*ccmano(i)+eff_rate_lit*
     $	cclito(i))*dt+eff_rate_hum*cchumo(i)*dt
c		the production of carbon humus
		cchumo(i)=cchumo(i)+fe*fh*(eff_rate_man*ccmano(i)+
     $		eff_rate_lit* cclito(i))*dt-eff_rate_hum*cchumo(i)*dt
c		the changes in the n-humus pool through mineralisation and n-litter humification
		cnhumo(i)=cnhumo(i)+(fe*fh/ro*(eff_rate_lit*cclito(i)+
     $		eff_rate_man*ccmano(i))-eff_rate_hum*cnhumo(i))*dt
c		the change in the manure pool
		ccmano(i)=ccmano(i)+(-fe*fh-(1-fe))*eff_rate_man*ccmano(i)*dt
c		the change in the litter pool
		cclito(i)=cclito(i)+(-fe*fh-(1-fe))*eff_rate_lit*cclito(i)*dt
c		retain the rate of change of the ammonia and nitrate species
c		to be incorporated in the sequential transformation process  (mg/(day*liter))
		rminm(i,2)=change_nh4/(dt*volume)
		rminm(i,3)=change_no3/(dt*volume)
		rmini(i,2)=0.0d0
		rmini(i,3)=0.0d0
		rcarbmin(i)=(1-fe)*(eff_rate_man*ccmano(i)+
     $		eff_rate_lit*cclito(i))/volume + eff_rate_hum*cchumo(i)/volume
      enddo
      return
      end