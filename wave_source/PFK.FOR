c	version 3.0
c	author m. vanclooster 15/02/2002

c######################################################################################
      double precision function calc_wc (suc_in, i_comp)
c     in   : i_comp, khys, modelpf, nla, parapf, suc_in, sucreturn, wcreturn
c     out  : calc_wc
c     calls: check_turn_point, init_pfk, report_err, stop_simulation
c######################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      integer kontro , prob_n
      double precision wclim,suclim,dmcaplim
      common /wc_dmc_prob/wclim(kt_layrs),suclim(kt_layrs),
     $dmcaplim(kt_layrs), kontro(kt_comps), prob_n(kt_layrs) 

      wcgenuch(suc,wcr,wcs,alpha,rn,rm) = 
     $wcr+(wcs-wcr)/((1.d0+(alpha*suc)**rn)**rm)
      call init_pfk
      n_layr =nla(i_comp)
      if (suc_in.ge.-mn_ph.and.modelpf.ne.4) then 
		calc_wc = parapf(n_layr,2)
		return
      endif
      suc =  dmin1( dabs( dmin1(suc_in, -mn_ph)),mx_ph)/10.d0

      goto (100, 200, 300, 400, 90,190,290),modelpf+(prob_n(n_layr)*kt_modelpf)
      call report_err  ('modelpf >4 or < 0')
      call stop_simulation('programme stopped: check err_file')
90    if (suc.lt.suclim(n_layr)) then
		calc_wc = wclim(n_layr) + dmcaplim(n_layr)* (suclim(n_layr)-suc)
	return
      endif

c     van genuchten model with 5 parameters
100   wcr=parapf(n_layr,1)
      wcs=parapf(n_layr,2)
      alpha=parapf(n_layr,3)
      rn=parapf(n_layr,4)
      rm=parapf(n_layr,5)
      calc_wc=wcgenuch(suc,wcr,wcs,alpha,rn,rm)
      return         
190   if (suc.lt.suclim(n_layr)) then
		calc_wc = wclim(n_layr) + dmcaplim(n_layr)*(suclim(n_layr)-suc)
	return
      endif

c     universal hysteresismodel (derived from mualem model 2)
200   call check_turn_point
      wcr=parapf(n_layr,1)
      wcs=parapf(n_layr,2)
      alphaw=parapf(n_layr,3)
      rnw=parapf(n_layr,4)
      rmw=parapf(n_layr,5)
c     control if the main loop is reached (kontro=1)
      kontro(i_comp)=0
      if(khys(i_comp).eq.+1)then
c		wetting                 
		v1=wcgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
		v2=wcgenuch(sucreturn(i_comp),wcr,wcs,alphaw,rnw,rmw)
		calc_wc=wcreturn(i_comp)+((v1-v2)*(wcs-v2)/wcs-wcr) +wcr
c		check if the wc is not outside the main loop
		if (calc_wc.lt.v1)then
			kontro(i_comp)=1
			calc_wc = v1
		endif
      elseif(khys(i_comp).eq.-1)then    
c		drying
		v1=wcgenuch(sucreturn(i_comp),wcr,wcs,alphaw,rnw,rmw)
		v2=wcgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
		calc_wc=wcreturn(i_comp)-((v1-v2)*(wcs-v2)/wcs-wcr)+wcr         
		wc_main_dry = ((2*wcs)-v2-wcr)*((v2-wcr)/(wcs-wcr))+wcr
		if(calc_wc.gt.wc_main_dry)then
			kontro(i_comp)=1
			calc_wc=wc_main_dry
		endif
      elseif(khys(i_comp).eq.-2)then
c		main drying
		v1=wcgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
		calc_wc=(2.d0*wcs -v1-wcr)*((v1-wcr)/(wcs-wcr))+wcr
      endif
      return
290   if (suc.lt.suclim(n_layr)) then
	calc_wc = wclim(n_layr) + 
     $  dmcaplim(n_layr)* (suclim(n_layr)-suc)
	return
      endif

c     hysteresismodel mualem 2
300   call check_turn_point
      wcr=parapf(n_layr,1)
      wcs=parapf(n_layr,2)
      alphad=parapf(n_layr,3)
      rnd=parapf(n_layr,4)
      rmd=parapf(n_layr,5)
      alphaw=parapf(n_layr,6)
      rnw=parapf(n_layr,7)
      rmw=parapf(n_layr,8)
      kontro(i_comp)=0
c	wetting                 
      if(khys(i_comp).eq.+1)then
		v1=wcgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
		v2=wcgenuch(sucreturn(i_comp),wcr,wcs,alphaw,rnw,rmw)
		v3=wcgenuch(sucreturn(i_comp),wcr,wcs,alphad,rnd,rmd)
		if( dabs(wcs-v2).lt.0.00001d0)then
			calc_wc= v1
			kontro(i_comp)=1
		else
			calc_wc=wcreturn(i_comp)+((v1-v2)*((wcs-v3)/(wcs-v2))) 
			if(calc_wc.lt.v1)then
				calc_wc=v1
				kontro(i_comp)=1
		endif
		endif
c	drying
      elseif(khys(i_comp).eq.-1)then
		v1=wcgenuch(sucreturn(i_comp),wcr,wcs,alphaw,rnw,rmw)
		v2=wcgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
		v3=wcgenuch(suc,wcr,wcs,alphad,rnd,rmd)
		if( dabs(wcs-v2).lt.0.00001d0)then
			calc_wc= v3
			kontro(i_comp)=1
		else
			calc_wc=wcreturn(i_comp)-((v1-v2)*((wcs-v3)/(wcs-v2)))
			if(calc_wc.gt.v3)then
				calc_wc=v3
				kontro(i_comp)=1
			endif
		endif
c       main drying
      elseif(khys(i_comp).eq.-2)then
		calc_wc=wcgenuch(suc,wcr,wcs,alphad,rnd,rmd)
      endif
      return
 	
c     multi-modal retention curve    
400   calc_wc= calc_wc_mm(suc_in,i_comp)
      return
      end


c######################################################################################
      double precision function calc_ph (wc_in, i_comp)
c     in   : i_comp, modelpf, nla, parapf, wc_in
c     out  : calc_ph
c     calls: calc_dmc, calc_wc, report_err, stop_simulation,
c######################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      double precision wclim,suclim,dmcaplim
      dimension table (100,kt_layrs)
      integer kontro , prob_n 
      logical f_time_ph_muniv
      common /wc_dmc_prob/wclim(kt_layrs),suclim(kt_layrs),
     $dmcaplim(kt_layrs), kontro(kt_comps), prob_n (kt_layrs)
      data f_time_ph_muniv/.true./

      phgenuch(wcin,wcr,wcs,alpha,rn, rm)=
     $((((wcin-wcr)/(wcs-wcr))**(-1.d0/rm)-1.)**(1./rn))/alpha

      call init_pfk
      n_layr =nla(i_comp)
      wcr=parapf(n_layr,1)
      wcs=parapf(n_layr,2)
      wcont =  dmin1( dmax1(wc_in,wcr+0.0001d0),wcs)
      goto (100, 200, 100, 400, 90,190,90),modelpf + (prob_n(n_layr)*kt_modelpf)
      call report_err ('modelpf >4 or < 0')
      call stop_simulation ('programme stopped: check err_file')
c     van genuchten model with 5 parameters
90    if (wcont.gt.wclim(n_layr)) then
		if (wcont.ge.wcs) then
			calc_ph = 0.d0
		else
			calc_ph = - dmax1( mx_ph/10.d0, suclim(n_layr) +
     $		(wclim(n_layr)-wcont)/dmcaplim(n_layr)  )*10.d0
		endif
		return
      endif
100   wcr=parapf(n_layr,1)
      wcs=parapf(n_layr,2)
      alpha=parapf(n_layr,3)
      rn=parapf(n_layr,4)
      rm=parapf(n_layr,5)
      if (wcont.gt.wcs) then 
		calc_ph = 0.d0
      else
		calc_ph= - phgenuch(wcont,wcr,wcs,alpha,rn,rm)*10.0d0
      endif
      return         

c     universal hysteresismodel derived from mualem model 2
190   if (wcont.gt.wclim(n_layr)) then
		if (wcont.gt.wcs) then
			calc_ph = 0.d0
		else
			calc_ph =  - dmax1( mx_ph/10.d0,  suclim(n_layr) +
     $		(wclim(n_layr)-wcont)/dmcaplim(n_layr))*10.d0
		endif
	return
      endif
200   if (f_time_ph_muniv) then
		f_time_ph_muniv = .false.
		komp_wc_call = 1
		do layer = 1, npl
			wcr=parapf(layer,1)
			wcs=parapf(layer,2)
205			if (nla(komp_wc_call).ne.layer) then
				komp_wc_call = komp_wc_call + 1
				goto 205
			endif
			do pfvalue = -5.d0, 5.d0, 0.1d0
				press_head = 10.d0** pfvalue
				thetaw_ph = calc_wc (-10.d0*press_head, komp_wc_call)
				table(idint((pfvalue+5.d0)*10),layer) =
     $			(2.d0*wcs - thetaw_ph)*(thetaw_ph - wcr)/(wcs-wcr)
			enddo
		enddo
      endif
      if (wc_in.gt. table(1,n_layr)) then
		wcs=parapf(n_layr,2)
		calc_pf = -5.d0 + (log10(mn_ph) + 5.d0)*
     $	  (dmin1(wc_in,wcs)-table(1,n_layr))/(wcs-table(1,n_layr))
		calc_ph = -10.d0 * 10**calc_pf
      else if (wc_in.lt.table(100,n_layr)) then
		wcr=parapf(n_layr,1)
		calc_pf = 5.d0 +(log10(mx_ph) - 5.d0)*
     $	  ( dmax1(wcr,wc_in)-table(100,n_layr))/(wcr-table(100,n_layr))
		calc_ph = -10.d0 * 10.d0** calc_pf
      else
c     conduct a binary search
		i_l  = 100
		i_h  = 1
210		i = (i_l+i_h)/2
		if (wc_in.gt.table(i,n_layr)) then
			i_l = i
		else
			i_h = i
		endif
		if (i_l.gt.i_h+1) goto 210
		calc_ind = dfloat(i_l)-(wc_in-table(i_l,n_layr))
     $       /(table(i_h,n_layr)-table(i_l,n_layr))
		calc_ph = -10.d0 * 10.d0** ((calc_ind/10.d0)-5.d0)
      endif
      return

c     multiporosity
400   calc_ph = calc_ph_mm(wc_in,i_comp)
      return
      end

c######################################################################################
      double precision function calc_dmc (suc_in, i_comp)
c     in   : i_comp, khys, modelpf, nla, parapf, suc_in, sucreturn, wcreturn
c     out  : calc_dmc
c     calls: check_turn_point, report_err, stop_simulation,
c######################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      logical first_time(kt_layrs)
      double precision wclim,suclim,dmcaplim,succrit(kt_layrs)
      integer kontro , prob_n 
      common /wc_dmc_prob/wclim(kt_layrs),suclim(kt_layrs),
     $dmcaplim(kt_layrs), kontro(kt_comps), prob_n (kt_layrs)
      data first_time/kt_layrs*.true./

      wcgenuch(suc,wcr,wcs,alpha,rn,rm)=
     $wcr+(wcs-wcr)/((1.d0+(alpha*suc)**rn)**rm)

      cgenuch(suc,wcr,wcs,alpha,rn,rm)= 
     $(wcs-wcr)*rn*rm*alpha*((1.d0+((alpha*suc)**rn))
     $**(-rm-1.d0))*((alpha*suc)**(rn-1.d0))

      call init_pfk
      n_layr =nla(i_comp)
      if (first_time(n_layr)) then
		first_time(n_layr) = .false.
		if (.not.multimodal) then
			wcr=parapf(n_layr,1)
			wcs=parapf(n_layr,2)
			alpha=parapf(n_layr,3)
			rn=parapf(n_layr,4)
			rm=parapf(n_layr,5)
			secrit= (wcs -0.00001 -wcr)/(wcs-wcr)
			succrit(n_layr) = (((1.d0 - secrit)/secrit)**(1/rn))/alpha
		else
			succrit(n_layr) = -mn_ph/10.d0
		endif
      endif
      suc = dmin1( dabs( dmin1(suc_in,-mn_ph,-10.d0*succrit(n_layr))),
     $     mx_ph)/10.d0
      goto (100, 200, 300, 400,90,190,290), modelpf+(prob_n(n_layr)*kt_modelpf)
      call report_err ('modelpf>4 or modelpf<0')
      call stop_simulation ('programme stopped: check err_file')
c     van genuchten model with 5 parameters
90    if (suc.lt.suclim(n_layr)) then
		calc_dmc = dmcaplim(n_layr)/10.d0
	return
      endif
100   wcr=parapf(n_layr,1)
      wcs=parapf(n_layr,2)
      alpha=parapf(n_layr,3)
      rn=parapf(n_layr,4)
      rm=parapf(n_layr,5)
      calc_dmc=cgenuch(suc,wcr,wcs,alpha,rn,rm)/10.d0
      return
c     universeel hysteresismodel (derived from mualem 2)
190   if (suc.lt.suclim(n_layr)) then
		calc_dmc = dmcaplim(n_layr)/10.d0
	return
      endif
200   call check_turn_point
      wcr=parapf(n_layr,1)
      wcs=parapf(n_layr,2)
      alphaw=parapf(n_layr,3)
      rnw=parapf(n_layr,4)
      rmw=parapf(n_layr,5)
      if(kontro(i_comp).eq.1)then
		if(khys(i_comp).eq.+1)then
c         wetting                 
		   calc_dmc=cgenuch(suc,wcr,wcs,alphaw,rnw,rmw)/10.d0
		else
c		drying
			c=cgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
			v1=wcgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
			calc_dmc=(2.d0*c*((wcs-v1)/(wcs-wcr)))/10.0d0
		endif 
      else
		if(khys(i_comp).eq.+1)then
c         wetting                 
			c=cgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
			v1=wcgenuch(sucreturn(i_comp),wcr,wcs,alphaw,rnw,rmw)
			calc_dmc=(c*(wcs-v1)/(wcs- wcr))/10.0d0
		else if (khys(i_comp).eq.-1)then
c         drying
			c=cgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
			v1=wcgenuch(sucreturn(i_comp),wcr,wcs,alphaw,rnw,rmw)
			v2=wcgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
			calc_dmc=((c/(wcs-wcr))*(v1+wcs-(2.d0*v2)))/10.0d0
		else
c         main drying
			c=cgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
			v1=wcgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
			calc_dmc=(2.d0*c*((wcs-v1)/(wcs-wcr)))/10.0d0
		endif
      endif 
      return
c     hysteresismodel mualem ii
290   if (suc.lt.suclim(n_layr)) then
		calc_dmc = dmcaplim(n_layr)/10.d0
	return
      endif
300   call check_turn_point
      wcr=parapf(n_layr,1)
      wcs=parapf(n_layr,2)
      alphad=parapf(n_layr,3)
      rnd=parapf(n_layr,4)
      rmd=parapf(n_layr,5)
      alphaw=parapf(n_layr,6)
      rnw=parapf(n_layr,7)
      rmw=parapf(n_layr,8)
      if(kontro(i_comp).eq.1)then
		if(khys(i_comp).eq.+1)then
c			main wetting                 
			calc_dmc=cgenuch(suc,wcr,wcs,alphaw,rnw,rmw)/10.0d0
		else
c			main drying
			calc_dmc=cgenuch(suc,wcr,wcs,alphad,rnd,rmd)/10.0d0
		endif 
      else
		if(khys(i_comp).eq.+1)then
c			scanning wetting
			c=cgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
			v1=wcgenuch(sucreturn(i_comp),wcr,wcs,alphad,rnd,rmd)
			v2=wcgenuch(sucreturn(i_comp),wcr,wcs,alphaw,rnw,rmw)
			calc_dmc=(c*((wcs-v1)/(wcs-v2)))/10.0d0
		elseif (khys(i_comp).eq.-1) then
c			scanning drying
			c1=cgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
			c2=cgenuch(suc,wcr,wcs,alphad,rnd,rmd)
			v1=wcgenuch(suc,wcr,wcs,alphad,rnd,rmd)
			v2=wcgenuch(suc,wcr,wcs,alphaw,rnw,rmw)
			v3=wcgenuch(sucreturn(i_comp),wcr,wcs,alphaw,rnw,rmw)
			t1=c1*(wcs-v1)/(wcs-v2)
			t2=c2*(v3-v2)/(wcs-v2)
			t3=c1*((wcs-v2)**(-2))*(v3-v2)*(wcs-v1)
			calc_dmc=(t1+t2-t3)/10.0d0
		else
c			main drying
			calc_dmc=cgenuch(suc,wcr,wcs,alphad,rnd,rmd)/10.0d0
		endif
      endif
      return
c     multiporosity
400   calc_dmc =  calc_dmc_mm(suc_in, i_comp)
      return
      end

c######################################################################################
      double precision function calc_con (suc_in, i_comp)
c     in   : i_comp, modelk, nla, parakon, parapf, suc_in
c     out  : calc_con
c     calls: calc_wc, report_err, stop_simulation,
c######################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      n_layr =nla(i_comp)
      suc =  dmin1( dabs( dmin1(suc_in,-mn_ph)),mx_ph)/10.d0
      goto (100, 200, 300, 400,500,600), modelk 
      call report_err ('incorrect choice of conductivity model')
      call stop_simulation('programme stopped: check err_file')
c     gardner k(h) model with 3 parameters
100   rksat=parakon(n_layr,1)
      b=parakon(n_layr,2)
      sn=parakon(n_layr,3)
      calc_con=rksat/(1.+(suc*b)**sn)*10.d0
      return
c     transformed gardner k(h) model with 3 parameters
200   calc_con = parakon(n_layr,1)* dexp(-parakon(n_layr,2)*suc)*10.d0
      return
c     gilham k(theta) model with 2 parameters
300   wcont=calc_wc(-suc*10.d0,i_comp)
      a=parakon(n_layr,1)
      sn=parakon(n_layr,2)
      calc_con=a*wcont**sn*10.d0
      return
c     brooks&corey k(se) model with 2 parameters
400   wcont=calc_wc(-suc*10.d0,i_comp)
      wcr=parapf(n_layr,1)
      wcs=parapf(n_layr,2)
      se=(wcont-wcr)/(wcs-wcr)
      rksat=parakon(n_layr,1)
      rlambda=parakon(n_layr,2)
      calc_con=rksat*se**((2.d0+3.d0*rlambda)/rlambda)*10.d0
      return
c     mualem k(se) with restriction m=1-1/n on h(theta)
500   wcont=calc_wc(-suc*10.d0,i_comp)
      wcr=parapf(n_layr,1)
      wcs=parapf(n_layr,2)
      se=(wcont-wcr)/(wcs-wcr)
      rksat=parakon(n_layr,1)
      xexp=parakon(n_layr,2)
      rmkon= parakon(n_layr,3)
      calc_con=rksat*((se)**xexp)*((1.d0-(1.d0-
     $(se**(1.d0/rmkon)))**rmkon)**2.d0)*10.d0
      return
c     multimodaal
600   calc_con = calc_con_mm(suc_in,i_comp)
      return
      end

c######################################################################################
      subroutine init_pfk 
c     in   : modelpf, ncs, npl, parapf
c     out  : ihist, khys
c     calls: report_err, stop_simulation
c######################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      double precision wclim,suclim,dmcaplim
      integer kontro , prob_n
      common /wc_dmc_prob/wclim(kt_layrs),suclim(kt_layrs),
     $dmcaplim(kt_layrs), kontro(kt_comps), prob_n (kt_layrs)
      logical first_time
      data first_time /.true./

      cgenuch(suc,wcr,wcs,alpha,rn,rm)= 
     $(wcs-wcr)*rn*rm*alpha*((1.d0+((alpha*suc)**rn))
     $**(-rm-1.d0))*((alpha*suc)**(rn-1.d0))
c	linear approximation of the moisture  retention curve near saturation 
      if (.not.first_time) return
      first_time = .false.
	goto(100,200,300,400), modelpf
c	van genuchten 5 parameters 
100   do i= 1,npl
		wcr=parapf(i,1)
		wcs=parapf(i,2)
		alpha=parapf(i,3)
		rn=parapf(i,4)
		rm=parapf(i,5)
		if(rn.lt.1.d0)then
			wclim(i)=wcs-0.01d0
			selim=(wclim(i)-wcr)/(wcs-wcr)
			suclim(i)= dmax1(((((1.d0/selim)**(1.d0/rm))-1.d0)
     $		**(1.d0/rn))/alpha,mn_ph/10.d0)
			dmcaplim(i)=
     $		cgenuch(suclim(i),wcr,wcs,alpha,rn,rm)
			prob_n(i) = 1
		endif
	enddo 
	return 
c	universal hysteresis model derived from mualem 2
200	do i= 1,npl
		wcr=parapf(i,1)
		wcs=parapf(i,2)
		alpha=parapf(i,3)
		rn=parapf(i,4)
		rm=parapf(i,5)
		if(rn.lt.1.d0)then
			wclim(i)=wcs-0.01d0
			selim=(wclim(i)-wcr)/(wcs-wcr)
			suclim(i)= dmax1(((((1.d0/selim)**(1.d0/rm))-1.d0)
     $		**(1.d0/rn))/alpha,mn_ph/10.d0)
			dmcaplim(i)=
     $		cgenuch(suclim(i),wcr,wcs,alpha,rn,rm)
			prob_n(i) = 1
		endif
	enddo 
      do  i = 1, ncs
		khys(i) = - 2
		ihist(i) = 0
	enddo 
      return
c	hysteresis model mualem 2
300   do i= 1,npl
		rn1=parapf(i,4)
		rn2=parapf(i,7)
		prob_n(i) = 0
		if(rn1.eq.rn2)then
			wcr=parapf(i,1)
			wcs=parapf(i,2)
			alpha=parapf(i,3)
			rn=parapf(i,4)
			rm=parapf(i,5)
			if(rn.lt.1.d0)then
				wclim(i)=wcs-0.01d0
				selim=(wclim(i)-wcr)/(wcs-wcr)
				suclim(i)= dmax1(((((1/selim)**(1/rm))-1.)
     $			**(1/rn))/alpha,mn_ph/10.d0)
				dmcaplim(i)=
     $			cgenuch(suclim(i),wcr,wcs,alpha,rn,rm)
				prob_n(i) = 1
			endif
		elseif (rn1.lt.1.d0.or.rn2.lt.1.d0) then
			write (*,*) 'the n parameters for the mrc are different'
		endif
	enddo
      do  i = 1, ncs
		khys(i) = - 2
		ihist(i) = 0
	enddo 
c	multimodal
400   return
      end

c######################################################################################
      subroutine check_turn_point
c     in   : hystlim, khys, modelpf, ncs, ph, ph1, t, wc
c     out  : ihist, khys, sucreturn, wcreturn
c     calls: init_pfk, report_err, stop_simulation
c######################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      logical first_time 
      dimension suchist(kt_comps,10), wchist(kt_comps,10)
      data first_time/.true./

      if (first_time) then
		first_time = .false.
		call init_pfk
      endif 
c     return points are only checked once per timestep
      if (t.ne.time_before) then
		time_before = t
      else
		return
      endif
      goto (100, 200, 200, 400), modelpf 

100   return
c hysteresis models
200   do i=1,ncs
		suc =  dmin1( dabs( dmin1(ph(i),-mn_ph)),mx_ph)/10.d0
		rkhys=dfloat(max0(khys(i),-1))
		if(ph(i).gt.-mn_ph) then
			reldif = 0.d0
		else
			reldif=(ph1(i)-ph(i))/ dabs(ph(i))
		endif
		if(reldif/rkhys.gt.hystlim)then
			if(ihist(i).lt.20)then
				ihist(i)=ihist(i)+1
				if(ihist(i).ge.1)then
					suchist(i,ihist(i))=suc
					wchist(i,ihist(i))=wc(i)
					wcreturn(i)=wc(i)
					sucreturn(i)=suc
				endif
			else
				call report_err('hysteresis history list too large: '//
     $			'take larger value for hystlim')
				call stop_simulation ('programme stopped: check err_file')
			endif
c			change direction index
c			-2  :  main drying
c			-1  :  scanning drying
c			+1  :  wetting
			if(khys(i).eq.1)then
				if(ihist(i).eq.0)then
					khys(i)=-2
				else
					khys(i)=-1
				endif
			else
				khys(i)=1
			endif
		endif
c		check if scanning loop has been closed
		if(ihist(i).gt.2)then
c         in case of tertiary (or higher order) scanning curve
			if((khys(i).eq.-1).and.(suc.gt.suchist(i,ihist(i)-1)))then
c				if drying and press. head more negative than at 
c				last but one return point, remove information
c				of the closed scanning loop
				ihist(i)=ihist(i)-2
				wcreturn(i)=wchist(i,ihist(i))
				sucreturn(i)=suchist(i,ihist(i))
			elseif((khys(i).eq.1).and.(suc.lt.suchist(i,ihist(i)-1)))then
c				if wetting and press. head less negative than at 
c				last but one return point, remove information
c				of the closed scanning loop
				ihist(i)=ihist(i)-2
				wcreturn(i)=wchist(i,ihist(i))
				sucreturn(i)=suchist(i,ihist(i))
			endif  
		elseif(ihist(i).eq.2)then
c			in case of secondary drying curve, check if
c			main drying curve has been reached
			if(suc.gt.suchist(i,1))then
c				remove information of the closed scanning loop
c				and set khys = -2 (= main drying curve)
				ihist(i)=0
				khys(i)=-2
			endif
		elseif(ihist(i).eq.1)then
c			in case of primary wetting curve, check if
c			saturation has been reached
			if(ph(i).ge.-mn_ph) ihist(i)=-1
		endif
	enddo
400   return
      end
      
c######################################################################################
      double precision function calc_ph_mm (suc_in, i_comp)
c######################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      stop 'initial profile should be pressure head profile'
      suc_in = suc_in
      i_comp = i_comp
      calc_ph_mm = 0.d0
      return
      end

c######################################################################################
      double precision function calc_wc_mm (suc_in, i_comp)
c######################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'multi.com'

      lyr = nla(i_comp)
      if (suc_in.ge.-mn_ph) then 
		calc_wc_mm = wcs(lyr)
		return
      endif
      suc =  dmin1( dabs( dmin1(suc_in, -mn_ph)),mx_ph)/10.d0
      se = 0.d0
      do i = 1, nr_poros(lyr)
		se = se + w(i,lyr)*(1+(alpha(i,lyr)*suc)**rn(i,lyr))**
     >    (-rm(i,lyr))
      enddo  
      calc_wc_mm = wcr(lyr)+(wcs(lyr)-wcr(lyr))*se
      return
      end

c######################################################################################
      double precision function calc_dmc_mm (suc_in,i_comp)
c######################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'multi.com'
      dimension dse(nr_por_class)
      
      lyr = nla(i_comp)
      suc =  dmin1( dabs( dmin1(suc_in,-mn_ph)),mx_ph)/10.d0
      do i = 1, nr_poros(lyr)
		dse(i) = 
     $	w(i,lyr)*rn(i,lyr)*rm(i,lyr)*alpha(i,lyr)*
     $	((1.d0+((alpha(i,lyr)*suc)**rn(i,lyr)))
     $	**(-rm(i,lyr)-1.d0))*((alpha(i,lyr)*suc)**
     $	(rn(i,lyr)-1.d0))
      enddo  
      dse_tot = 0.d0
      do i = 1, nr_poros(lyr)
		dse_tot = dse_tot + dse(i)
      enddo
      calc_dmc_mm = (wcs(lyr)-wcr(lyr))*dse_tot/10.d0
      return
      end

c######################################################################################
      double precision function calc_con_mm (suc_in,i_comp)
c######################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'multi.com'
      
      lyr = nla(i_comp)
      suc =  dmin1( dabs( dmin1(suc_in, -mn_ph)),mx_ph)/10.d0
      se = 0.d0
      teller = 0.d0
      do i = 1, nr_poros(lyr)
		se = se + w(i,lyr)*(1+(alpha(i,lyr)*suc)**rn(i,lyr))**
     >    (-rm(i,lyr))
		si = (1+(alpha(i,lyr)*suc)**rn(i,lyr))**(-1)
		teller = teller + 
     >	w(i,lyr)*alpha(i,lyr)*rm(i,lyr)*
     >	betai(p(i,lyr),q(i,lyr),si)*beta(p(i,lyr),q(i,lyr))
      enddo
      calc_con_mm = csat(lyr)*(se**xexp(lyr))*((teller/total(lyr))
     >           **2)*10.d0
      return
      end

c######################################################################################
      subroutine init_n_p
c######################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'multi.com'
      integer rd_int

      do nl = 1, npl
		call hold_line
		nr_poros(nl) = rd_int(1, nr_por_class)
		wcr(nl) = rd_real (0.d0, 1.d0)
		wcs(nl) = rd_real (wcr, 1.d0)
		xexp(nl) = rd_real (-100.d0, 100.d0)
		csat(nl) = rd_real (0.d0, mx_real)
		call next_line
		total(nl) = 0.d0
		do i = 1, nr_poros(nl)
			call hold_line
			w(i,nl)     =  rd_real (0.d0,1.d0)
			alpha(i,nl) =  rd_real (mn_real, 10.d0)
			rn(i,nl)    =  rd_real (1.d0  , mx_real)
			rm(i,nl)    =  rd_real (0.d0  , mx_real)
			call next_line
			p(i,nl) = rm(i,nl) +(1/rn(i,nl))
			q(i,nl) = 1 - (1/rn(i,nl))
			total(nl) = total(nl)+ w(i,nl)*alpha(i,nl)*rm(i,nl)*
     >        beta(p(i,nl),q(i,nl))
		enddo
      enddo
      call check_et 
      end

c######################################################################################
      subroutine wr_init_n_p(res_file)
c######################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'multi.com'
      integer	   res_file

      do nl = 1, npl 
		write (res_file, 100) nl, nr_poros(nl), wcr(nl),wcs(nl),
     >		xexp(nl), csat(nl)
		do i = 1, nr_poros(nl)
			write (res_file, 110) i, w(i,nl), alpha(i,nl),
     >		rn(i,nl), rm(i,nl)
		enddo
      enddo
100   format ('layer: ',i1,/,
     >        '.........',/,
     >        'nr of porosities: ',i2,/,
     >        '  wcr       wcs       lexp     csat',/,
     >        '(cm3cm-3) (cm3cm-3)    (-)   (cmday-1)',/,
     >         f7.5,3x,f7.5,3x, f7.3,1x,f10.3,/,
     >        'por     w    alpha      rn       rm',/,
     >        '             (cm-1)            ')
110   format (i3,1x,f7.5,1x,f7.5,1x, f8.5,1x, f8.5)
      end

c######################################################################################
      subroutine param_one_layer(res_file,nl)
c######################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'multi.com'
      integer res_file

      write (res_file, 100) nl, nr_poros(nl), wcr(nl),wcs(nl),
     >  xexp(nl), csat(nl)
      do i = 1, nr_poros(nl)
		write (res_file, 110) i, w(i,nl), alpha(i,nl),
     >		rn(i,nl), rm(i,nl)
      enddo
      write (res_file, 120) 
100   format ('parameters for layer: ',i1,/,
     >        '-----------------------',/,
     >        'nr of porosities: ',i2,/,
     >        '  wcr       wcs       lexp     csat',/,
     >        '(cm3cm-3) (cm3cm-3)    (-)   (cmday-1)',/,
     >         f7.5,3x,f7.5,3x, f7.3,1x,f10.3,/,
     >        'por     w    alpha      rn       rm',/,
     >        '             (cm-1)            ')
110   format (i3,1x,f7.5,1x,f7.5,1x, f8.5,1x, f8.5)
120   format ('-------------------------------------',/)
      end


c######################################################################################
      double precision function betai (a,b,x)
c######################################################################################
      implicit double precision (a-h,o-z)
      if(x.lt.0.d0.or.x.gt.1.d0) stop 'bad argument x in betai'
      if(x.eq.0.d0.or.x.eq.1.d0) then
		bt = 0.d0
      else
		bt = dexp(gammln(a+b)-gammln(a)-gammln(b)+
     >	a*dlog(x)+b*dlog(1.d0-x))
      endif
      if(x.lt.(a+1.d0)/(a+b+2.d0)) then
		betai  = bt*betacf(a,b,x)/a
      else
		betai  = 1.d0-bt*betacf(b,a,1.d0-x)/b
      endif
      return
      end

c######################################################################################
      double precision function betacf (a,b,x)
c######################################################################################
      implicit double precision (a-h,o-z)
      parameter (itmax = 100, eps = 0.0000003d0)
      am = 1.d0
      bm = 1.d0
      az = 1.d0
      qab = a+b
      qap = a+1.d0
      qam = a-1.d0
      bz = 1.d0-qab*x/qap
      do m = 1,itmax
		em = m
		tem = em+em
		d = em*(b-m)*x/((qam+tem)*(a+tem))
		ap = az + d*am
		bp = bz + d*bm
		d = -(a+em)*(qab+em)*x/((a+tem)*(qap+tem))
		app = ap+d*az
		bpp = bp+d*bz
		aold = az
		am = ap /bpp
		bm = bp /bpp
		az = app /bpp
		bz = 1.d0
		if ( dabs(az-aold).lt.eps* dabs(az)) goto 1
      enddo
      stop 'a or b too big or itmax too small in betacf'
1     betacf = az
      return
      end
	 
	  
c######################################################################################
      double precision function beta (a,b)
c######################################################################################
      implicit double precision (a-h,o-z)
      beta = dexp(gammln(a+b)-gammln(a)-gammln(b))
      return
      end
 
c######################################################################################    
      double precision function gammln (xx)
c######################################################################################
      implicit double precision (a-h,o-z)
      dimension cof(6)
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     > -1.231739516d0,.120858003e-2,-.536382e-5,2.50662827465d0/
      data half,r_one,fpf/0.5d0,1.0d0,5.5d0/
      x = xx - r_one
      tmp = x +fpf
      tmp = (x+half)*dlog(tmp)-tmp
      ser = r_one
      do j = 1,6
		x = x+r_one
		ser = ser +cof(j)/x
      enddo
      gammln = tmp + dlog(stp*ser)
      return
      end