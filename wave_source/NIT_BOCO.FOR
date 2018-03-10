c	version 3.0
c	author: m. vanclooster 24/02/2002
      subroutine nit_boco
c###################################################################################
c     in   : carborg, cchumo, cclito, ccmano, cnhumo, cnlito, cnmano,
c            fag, fhp, flr, harvest_date, idman, idplow, ncs, ncsplow,
c            nof, nplo, rdens_left, rfract_tot, rnitorg, rnmaxp, ro,
c            simplant, t, tflcorg, tflnorg
c     out  : carblit, carbman, cchumo, cclito, ccmano, cnhumo, cnlito,
c            cnmano, tflcorg, tflnorg
c     calls: -
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'gen.com'
      include  'sol.com'
      include  'nit.com'
      data itelof,itelplow/1,1/

      carblit=0.0d0
      carbman=0.0d0
      rnitman=0.0d0
      rnitlit=0.0d0

c     manuring
      if(nof.gt.0) then
		if(idint(t).gt. idman(itelof)) then
			itelof=itelof+1
			if(itelof.gt.nof) itelof=nof
		endif
		if(idint(t).eq. idman(itelof)) then
			carbman = carborg(itelof)*0.3d0
			rnitman = carbman*(1/ro)
			carblit = carborg(itelof)-carbman
			rnitlit = rnitorg(itelof)-rnitman
			cnmano(1) = cnmano(1)+rnitman
			ccmano(1) = ccmano(1)+carbman
			cnlito(1) = cnlito(1)+rnitlit
			cclito(1) = cclito(1)+carblit
		endif
      endif

c     harvest
      if (simplant .and. idint(t).eq.harvest_date) then
		do i=1,ncs
			rnitlit = (rdens_left(i)/rfract_tot)*(1.d0-fag-flr-fhp)*rnmaxp
			carblit = rnitlit*10.0d0
			cnlito(i) = cnlito(i)+rnitlit
			cclito(i) = cclito(i)+carblit
		enddo 
      endif

c     plowing
      if(nplo.gt.0) then
		if(idint(t).gt.idplow(itelplow)) then
			itelplow = itelplow +1
			if (itelplow.gt.nplo) itelplow = nplo
		endif
		if(idint(t) .eq. idplow(itelplow)) then
			plcnlit =  0.d0
			plcclit =  0.d0
			plcnman =  0.d0
			plccman =  0.d0
			plcnhum =  0.d0
			plcchum =  0.d0
			do i=1,ncsplow(itelplow)
				plcnlit = plcnlit+cnlito(i)
				plcclit = plcclit+cclito(i)
				plcnman = plcnman+cnmano(i)
				plccman = plccman+ccmano(i)
				plcnhum = plcnhum+cnhumo(i)
				plcchum = plcchum+cchumo(i)
			enddo
			if(simplant.and.itelplow.eq.1) then
				rnitlit  = flr*rnmaxp+fag*rnmaxp
				carblit  = flr*rnmaxp*10.0d0+fag*rnmaxp*50.0d0
			else
				rnitlit = 0.0d0
				carblit = 0.0d0
			endif
			plcnlit = plcnlit+rnitlit
			plcclit = plcclit+carblit
			do i=1,ncsplow(itelplow)
				cnlito(i) = plcnlit/ncsplow(itelplow)
				cclito(i) = plcclit/ncsplow(itelplow)
				cnmano(i) = plcnman/ncsplow(itelplow)
				ccmano(i) = plccman/ncsplow(itelplow)
				cnhumo(i) = plcnhum/ncsplow(itelplow)
				cchumo(i) = plcchum/ncsplow(itelplow)
			enddo
		endif
      endif
c     inflow
      tflnorg = tflnorg+rnitlit+rnitman
      tflcorg = tflcorg+carblit+carbman
      return
      end