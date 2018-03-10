c	version 3.0
c	author: m. vanclooster 15/02/2002
c###################################################################################
      subroutine sol_intgr
c     in   : acsolio, acsolmo, acsolo, csolio, csolmo, csolo, decsoli,
c            decsolm, dt, dx, imobsw, isd, ncs, nla, nr_of_sol, parasol,
c            sflxsa, t, tcsolio, tcsolmo, tcsolo, tflsol, tsoli, 
c            tsoli_isd, wcio, wcmo, wco
c     out  : cberr_sol, cberr_sol_isd, dleasa, dleasa_isd, dsol, dsol_isd,
c            pleasa, pleasa_isd, rleasa, rleasa_isd, solinfl, tcsink, 
c            tcsink_isd, tcsolio, tcsolmo, tcsolo, tflsol, tsinki, tsinkm
c     calls: -
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      integer sp

      do sp = 1,nr_of_sol
c		calculate new total values of solute on the nodes
		volume=dx
		do i=1,ncs
			if(imobsw)then
				tcsolmo(i,sp)=(csolmo(i,sp)*wcmo(i)+(acsolmo(i,sp)*
     $			parasol(1,nla(i),sp)))*volume
				tcsolio(i,sp)=(csolio(i,sp)*wcio(i)+(acsolio(i,sp)*
     $			parasol(1,nla(i),sp)))*volume
				tcsolo(i,sp)=(tcsolmo(i,sp)+tcsolio(i,sp))
			else
				tcsolo(i,sp)=(csolo(i,sp)*wco(i)+(acsolo(i,sp)*
     $			parasol(1,nla(i),sp)))*volume
				tcsolmo(i,sp)=tcsolo(i,sp)
				tcsolio(i,sp)=0.0d0
			endif
		enddo 
c		calculate profile totals for this time 
		ptsol=0.0
		ptsoli=0.0
		ptsolm=0.0
		ptsinki=0.0
		ptsinkm=0.0
		do i=1,ncs
			ptsoli=ptsoli+tcsolio(i,sp)
			ptsolm=ptsolm+tcsolmo(i,sp)
			ptsol=ptsol+tcsolo(i,sp)
			if (i.eq.isd) ptsol_isd = ptsol
			tsinki(i,sp) = tsinki(i,sp) +(sinki(i,sp))*volume*dt
			tsinkm(i,sp) = tsinkm(i,sp) +(sinkm(i,sp)+
     $			parasol(8,nla(i),sp)
     $			*(csolmo(i,sp)-csolio(i,sp)))*volume*dt
			ptsinki = ptsinki +(sinki(i,sp))*volume*dt
			ptsinkm = ptsinkm +(sinkm(i,sp)+parasol(8,nla(i),sp)
     $			*(csolmo(i,sp)-csolio(i,sp)))*volume*dt
			if (i.eq.isd) ptsink_isd = ptsinki + ptsinkm 
		enddo 
		ptsink=ptsinki+ptsinkm

c		calculate cumulative values for the whole profile
c		netto flux
		rleasa(sp)=rleasa(sp)+sflxsa(ncs+1,sp)
c		cumulative downward flux (negativ)
		if(sflxsa(ncs+1,sp).le.0.0) dleasa(sp)=dleasa(sp)+sflxsa(ncs+1,sp)
c		cumulativ upward flux (positiv)
		if(sflxsa(ncs+1,sp).gt.0.0) pleasa(sp)=pleasa(sp)+sflxsa(ncs+1,sp)
c		inflow
		tflsol(sp) = tflsol(sp)+sflxsa(1,sp)
c		inflow during one day
		if(t- dint(t).eq.0.) solinfl(sp)=0.0
		solinfl(sp)=solinfl(sp)+sflxsa(1,sp)
c		sink first order
		tcsink(sp) = tcsink(sp) + ptsink
c		change in the system from the start till the present time step
		dsol(sp)=ptsol-tsoli(sp)
c		mass balance error 
		cberr_sol(sp)=dsol(sp) - (rleasa(sp)+tflsol(sp)+plo_in(sp)+tcsink(sp))

c		calculate cumulative values for the isd profile
c		netto flux
		rleasa_isd(sp)=rleasa_isd(sp)+sflxsa(isd+1,sp)
c		cumulative downward flux (negativ)
		if(sflxsa(isd+1,sp).le.0.0) 
     $	dleasa_isd(sp)=dleasa_isd(sp)+sflxsa(isd+1,sp)
c		cumulativ upward flux (positiv)
		if(sflxsa(isd+1,sp).gt.0.0)
     $	pleasa_isd(sp)=pleasa_isd(sp)+sflxsa(isd+1,sp)
c		sink first order
		tcsink_isd(sp) = tcsink_isd(sp) + ptsink_isd
c		change in the system from the start till the present time step
		dsol_isd(sp)=ptsol_isd-tsoli_isd(sp)
c		mass balance error 
		cberr_sol_isd(sp)=dsol_isd(sp) - (rleasa_isd(sp)+tflsol(sp)+
     $	tcsink_isd(sp)+plo_in(sp))
	enddo 
      return
      end
