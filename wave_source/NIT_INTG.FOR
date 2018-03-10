c	version 3.0
c	author: m. vanclooster 24/02/2002
      subroutine nit_intgr
c###################################################################################
c     in   : cchumo, cclito, ccmano, cnhumo, cnlito, cnmano, deccorg,
c            decnorg, dt, dx, isd, ncs, rdeniti, rdenitm, rhureai, rhuream,
c            rmini, rminm, rnitrii, rnitrim, rvoli, rvolm, suni, sunm,
c            tcorgs, tcorgs_isd, tflcorg, tflnorg, tnorgs, tnorgs_isd
c     out  : cberr_corg, cberr_corg_isd, cberr_norg, cberr_norg_isd,
c            cden, chyd, cmin, cnit, cupt, cvol, dcorg, dcorg_isd, dnorg,
c            dnorg_isd, ptcarborg, ptcarborg_isd, ptchum, ptchum_isd,
c            ptclit, ptclit_isd, ptcman, ptcman_isd, ptden, ptden_isd,
c            pthyd, pthyd_isd, ptmin, ptmin_isd, ptnhum, ptnhum_isd,
c            ptnit, ptnitorg, ptnitorg_isd, ptnit_isd, ptnlit, ptnlit_isd,
c            ptnman, ptnman_isd, ptscorg, ptscorg_isd, ptsnorg, ptsnorg_isd,
c            ptup, ptup_isd, ptvol, ptvol_isd
c     calls: -
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      include   'nit.com'

      volume=dx
c	profile totals 
      ptclit=0.d0
      ptnlit=0.d0
      ptcman=0.d0
      ptnman=0.d0
      ptchum=0.d0
      ptnhum=0.d0
      do i=1,ncs
		ptclit=ptclit+cclito(i)
		ptnlit=ptnlit+cnlito(i)
		ptcman=ptcman+ccmano(i)
		ptnman=ptnman+cnmano(i)
		ptchum=ptchum+cchumo(i)
		ptnhum=ptnhum+cnhumo(i)
		ptsnorg=ptsnorg+decnorg(i)*volume*dt
		ptscorg=ptscorg+deccorg(i)*volume*dt
		if (i.eq.isd) then
			ptclit_isd=ptclit
			ptnlit_isd=ptnlit
			ptcman_isd=ptcman
			ptnman_isd=ptnman
			ptchum_isd=ptchum
			ptnhum_isd=ptnhum
		endif
      enddo
      ptnitorg=ptnlit+ptnman+ptnhum
      ptcarborg=ptclit+ptcman+ptchum
c	isd totals      
	do i = 1, isd
		ptsnorg_isd=ptsnorg_isd+decnorg(i)*volume*dt
		ptscorg_isd=ptscorg_isd+deccorg(i)*volume*dt
      enddo
      ptnitorg_isd=ptnlit_isd+ptnman_isd+ptnhum_isd
      ptcarborg_isd=ptclit_isd+ptcman_isd+ptchum_isd

c     profile sink  
      do i=1,ncs
		ptup(2)=ptup(2)+(sunm(i,2)+suni(i,2))*volume*dt
		cupt(i,2)=cupt(i,2)+(sunm(i,2)+suni(i,2))*volume*dt
		ptmin(2)=ptmin(2)+(rminm(i,2)+rmini(i,2))*volume*dt
		cmin(i,2)=cmin(i,2)+(rminm(i,2)+rmini(i,2))*volume*dt
      enddo
      do i=1,ncs
		ptup(3)=ptup(3)+(sunm(i,3)+suni(i,3))*volume*dt
		cupt(i,3)=cupt(i,3)+(sunm(i,3)+suni(i,3))*volume*dt
		cunc_no3(i) = cunc_no3(i)+unc(i,3)
		cund_no3(i) = cund_no3(i)+und(i,3)
		ptmin(3)=ptmin(3)+(rminm(i,3)+rmini(i,3))*volume*dt
		cmin(i,3)=cmin(i,3)+(rminm(i,3)+rmini(i,3))*volume*dt

      enddo
      do i=1,ncs
		pthyd=pthyd+(rhuream(i)+rhureai(i))*volume*dt
		chyd(i)=chyd(i)+(rhuream(i)+rhureai(i))*volume*dt
		ptnit=ptnit+(rnitrim(i)+rnitrii(i))*volume*dt
		cnit(i)=cnit(i)+(rnitrim(i)+rnitrii(i))*volume*dt
		ptvol=ptvol+(rvolm(i)+rvoli(i))*volume*dt
		cvol(i)=cvol(i)+(rvolm(i)+rvoli(i))*volume*dt
		ptden=ptden+(rdenitm(i)+rdeniti(i))*volume*dt
		cden(i)=cden(i)+(rdenitm(i)+rdeniti(i))*volume*dt
      enddo
c	isd sink
      do i = 1, isd
		ptup_isd(2)=ptup_isd(2)+(sunm(i,2)+suni(i,2))*volume*dt
		ptmin_isd(2)=ptmin_isd(2)+(rminm(i,2)+rmini(i,2))*volume*dt
      enddo   
      do i = 1, isd
		ptup_isd(3)=ptup_isd(3)+(sunm(i,3)+suni(i,3))*volume*dt
		ptmin_isd(3)=ptmin_isd(3)+(rminm(i,3)+rmini(i,3))*volume*dt
      enddo        
      do i = 1, isd
		pthyd_isd=pthyd_isd+(rhuream(i)+rhureai(i))*volume*dt
		ptnit_isd=ptnit_isd+(rnitrim(i)+rnitrii(i))*volume*dt
		ptvol_isd=ptvol_isd+(rvolm(i)+rvoli(i))*volume*dt
		ptden_isd=ptden_isd+(rdenitm(i)+rdeniti(i))*volume*dt
      enddo
c	organic balance terms
      dnorg=ptnitorg - tnorgs
      dcorg=ptcarborg - tcorgs
      dnorg_isd=ptnitorg_isd - tnorgs_isd
      dcorg_isd=ptcarborg_isd - tcorgs_isd
c     organic mass balance error 
      cberr_norg= dnorg - (tflnorg + ptsnorg)
      cberr_corg= dcorg - (tflcorg + ptscorg)
      cberr_norg_isd= dnorg_isd - (tflnorg + ptsnorg_isd)
      cberr_corg_isd= dcorg_isd - (tflcorg + ptscorg_isd)
      return
      end