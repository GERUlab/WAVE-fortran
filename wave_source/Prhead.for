c	version 3.0
c	author m. vanclooster 15/02/2002

c###################################################################################
      subroutine prhead
c     in   : conin, dmcap, dt, dx, dxinter, itertot, ncs, phsurf, rtex,pond_max
c     out  : flxs, iter, itertot, ph, ph1, wc, wc1,pond
c     calls: balance_error, calc_bboco, calc_uboco, calc_wc, check_bboco,
c            check_uboco, dmccon, fix_uboco, halve, rer
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      dimension r1(kt_comps),r2(kt_comps)
      logical bot_ok,top_ok, flxs_cond_bot,flxs_cond_top,
     $        free_drainage,balance_error, no_conv

      no_conv = .false.
      call dmccon
      call rer
      do i=1,ncs
		ph1(i)=ph(i)
		wc1(i)=wc(i)
      enddo
20    call calc_bboco(ncbot,phbot,flxsbot,flxs_cond_bot,free_drainage)
      call calc_uboco(phsurf,flxsbot,flxar,flxs_cond_top)
c     start of iteration loop (repeat if boundary condition is not correct) 
      itertop = 0
30    itertop= itertop + 1
c     tridiagonal algorithm : forward substitution 
      iter = 0
c     start of iteration loop
40    iter=iter +1
      itertot=itertot+1
      if(.not.flxs_cond_top)then
c		calculation of first node, pressure head condition
		h0=dt/dx
		h1=h0/dxinter(1)
		h2=h0/dxinter(2)
		a=h2*conin(2)
		c=h1*conin(1)
		b=dmcap(1)+a+c
		e=c*phsurf+dmcap(1)*ph(1)-wc(1)+wc1(1)
     $	+h0*(conin(1)-conin(2))-dt*rtex(1)
		r1(1)=a/b
		r2(1)=-e/a
      else
c		calculation of first node, flux condition
		h0=dt/dx
		h1=h0/dxinter(1)
		h2=h0/dxinter(2)
		a=h2*conin(2)
		b=dmcap(1)+a
		e=dmcap(1)*ph(1)-wc(1)+wc1(1)
     $		-h0*(flxar+conin(2))-dt*rtex(1)
		r1(1)=a/b
		r2(1)=-e/a
      endif
      do  i=2,ncbot-1
c		calculation of coefficients for 1<i<ncbot  
		h0=dt/dx
		h1=h0/dxinter(i)
		h2=h0/dxinter(i+1)
		a=h2*conin(i+1)
		c=h1*conin(i)
		b=dmcap(i)+a+c
		e=dmcap(i)*ph(i)-wc(i)+wc1(i)
     $	+h0*(conin(i)-conin(i+1))-dt*rtex(i)
c		forward substitution
		r1(i)=a/(b-c*r1(i-1))
		r2(i)=(c*r1(i-1)*r2(i-1)-e)/a
	enddo 
      if(.not.flxs_cond_bot)then
c		calculation of coefficients for i=ncbot,pressure head condition
		h0=dt/dx
		h1=h0/dxinter(ncbot)
		h2=h0/dxinter(ncbot+1)
		a=h2*conin(ncbot+1)
		c=h1*conin(ncbot)
		b=dmcap(ncbot)+a+c
		e=dmcap(ncbot)*ph(ncbot)-wc(ncbot)+wc1(ncbot)
     $	+h0*(conin(ncbot)-conin(ncbot+1))-dt*rtex(ncbot)
c		forward substitution
		r1(ncbot)=a/(b-c*r1(ncbot-1))
		r2(ncbot)=(c*r1(ncbot-1)*r2(ncbot-1)-e)/a
c		backward substitution for node ncbot  
		ph(ncbot)=r1(ncbot)*(phbot-r2(ncbot))
      else
c		calculation of coefficients for i=ncbot, flux condition
		h0=dt/dx
		h1=h0/dxinter(ncbot)
		c=h1*conin(ncbot)
		b=dmcap(ncbot)+c
		e=dmcap(ncbot)*ph(ncbot)-wc(ncbot)+wc1(ncbot)
     $	+h0*(conin(ncbot)+flxsbot)-dt*rtex(ncbot)
		ph(ncbot)=(e-c*r1(ncbot-1)*r2(ncbot-1))/(b-c*r1(ncbot-1))
		phbot=ph(ncbot)+dxinter(ncbot+1)*(flxsbot/conin(ncbot+1)+1)
      endif
c     tridiagonal algorithm : backward substitution   
      do j=ncbot,2,-1
		ph(j-1)=r1(j-1)*(ph(j)-r2(j-1))
	enddo
c     calculation of diff. moist. capacities and hydr. conductivities
      call dmccon
c	calculation of the sink/source	
      call rer
c     calculation of water contents and fluxes
      if(flxs_cond_top)then
		flxs(1)=flxar
      else
		flxs(1)=-conin(1)*( (phsurf-ph(1))/dxinter(1) +1.)
      endif
      do i=1,ncbot -1
		wc(i)=calc_wc(ph(i),i)
		flxs(i+1)=-conin(i+1)*( (ph(i)-ph(i+1))/dxinter(i+1) +1.)
	enddo
      wc(ncbot)=calc_wc(ph(ncbot),ncbot)
      flxs(ncbot+1)=-conin(ncbot+1)*((ph(ncbot)-phbot)/dxinter(ncbot+1)+1.)
c     if there was no convergence check balance error 
      if(.not.no_conv.and.balance_error(ncbot)) then
c		check if max. number of iterations is exceeded,
		if(iter.ge.maxiter)then
c			restore situation and  restart loop
			call halve(no_conv)
			do i=1,ncs
				ph(i)=ph1(i)
				wc(i)=wc1(i)
			enddo 
			call dmccon
			call rer
			goto 20
		else
c			calculation of flxsbot in case of free drainage
			if(free_drainage)flxsbot=-conin(ncbot+1)
			goto 40
c			end of iteration loop
		endif
      endif
c     exit newton-raphson iteration loop if internal balance errors show 
c     to be small enough  (balance_error = false)
      call check_bboco(phbot,flxsbot,flxs_cond_bot ,bot_ok)
      call check_uboco(phsurf,flxs_cond_top,top_ok)
      if(.not.(top_ok.and.bot_ok).and.itertop.lt.3) then
c		in case of incorrect upper boundary condition 
c		less than three iterations for top boundary condition
		if(itertop.eq.2) call fix_uboco (phsurf)
c		start again from previous pressure head profile      
		do i=1,ncs
			ph(i)=ph1(i)
			wc(i)=wc1(i)      
		enddo 
		call dmccon
		call rer
c		recalculate flxsbot in case of free drainage
		if(free_drainage)flxsbot=-conin(ncbot+1)
		goto 30
      else
c		in case of correct upper boundary condition 
c         three iterations for top boundary condition
		do k=ncbot+ 1,ncs
			wc(k)=calc_wc(0.d0,k)
			ph(k)=ph(k-1)+(flxs(k)/conin(k)+1.)*dxinter(k)
			flxs(k+1)=flxs(k)+((wc(k)-wc1(k))/dt+rtex(k))*dx
		enddo 
		pond =  dmin1(pond_max, dmax1(0.d0, (flxs(1) - flxa)*dt))
		return
      endif
c     end of iteration loop for top boundary condition
      end

c###################################################################################
      logical function balance_error(ncbot)
c     in   : devstop, dt, dx, flxa, flxs, histor_file, iter, ncbot, pr_wat_histor, 
c            rtex, t, wc, wc1
c     out  : balance_error
c     calls: -
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      errmax=0.d0
      ncerrmax=0
      do i=1,ncbot
		err= dabs((wc(i)-wc1(i))/dt-(-flxs(i)+flxs(i+1))/dx
     $	+rtex(i))
		if(err.gt.errmax)then
			ncerrmax=i
			errmax=err
		endif
	enddo 
      balance_error = errmax.gt.devstop
c     printing of iteration history
      if(pr_wat_histor)write(histor_file,20)iter,t,flxa, flxs(1),errmax,ncerrmax
20    format(5x,i3,'  t=',f9.4,'   flxa=',f12.7,' flxs(1)='
     $,f10.5,'  errmax=',e10.3,'  ncerrmax=',i3)
      return
      end

c###################################################################################
      subroutine dmccon
c     in   : ncs, ph, phsa, phsurf
c     out  : conduc, conin, dmcap
c     calls: calc_con, calc_dmc
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'

      do i=1,ncs
		dmcap(i)=calc_dmc (ph(i),i)
		conduc(i)= calc_con (ph(i),i)
	enddo 
      conin(1)=sqrt(calc_con( dmax1(phsa,phsurf),1)*conduc(1))
      do i=2,ncs
		conin(i)=sqrt(conduc(i)*conduc(i-1))
	enddo
      conin(ncs+1)=conduc(ncs)
      return
      end