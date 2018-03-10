c	Version 3.0
c	author : M. vanclooster 18/02/02

c###################################################################################
      subroutine soln
c     in   : acsolio, acsolmb, acsolmo, csolb, csolio, csolmb, csolmo,
c            csolo, diffus, dt, dx, flxsa, flxsah, imobsw, ncs, nla,
c            nr_of_sol, parasol, pvela, pvelah, pveloh, sinki, sinkm,
c            solsur, wcio, wciob, wcma, wcmah, wcmo, wcmob, wcob
c     out  : acsolib, acsolio, acsolmb, acsolmo, acsolo, csolb, csolib,
c            csolio, csolmb, csolmo, csolo, diffus, sflxsa
c     calls: -
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      integer sp
      dimension diffush(kt_comps+1)
      dimension ssb(kt_sol),ss(kt_sol)
      dimension betha(kt_comps),y(kt_comps)
      dimension diffusb(kt_comps +1 , kt_sol)


      do sp = 1, nr_of_sol
c		set actual values to passed values
		do i=1,ncs
			csolb(i,sp)=csolo(i,sp)
			csolmb(i,sp)=csolmo(i,sp)
			csolib(i,sp)=csolio(i,sp)
			acsolib(i,sp)=acsolio(i,sp)
			acsolmb(i,sp)=acsolmo(i,sp)
		enddo
c		the solution assumes use of mobile/immobile
		if(.not.imobsw) then
			do i=1,ncs
				csolmb(i,sp)=csolb(i,sp)
				csolib(i,sp)=0.0d0
			enddo
		endif  
c		set the lower boundary condition zero gradient 
		csolmb(ncs+1,sp)=csolmb(ncs,sp)
c		calculate -diffusion coefficients
		do i=1,ncs+1
			dp=parasol(3,nla(i),sp)*
     $		parasol(4,nla(i),sp)* dexp(parasol(5,nla(i),sp)*wcmo(i))
			dm= dmax1(0.d0,parasol(6,nla(i),sp)-(0.139d0*dx))*
     $			dabs(pvela(i))
			diffusb(i,sp) = diffus(i,sp)
			diffus(i,sp)=(dp/wcma(i))+(dm/wcma(i))
		enddo
		do i=1,ncs+1
			diffush(i)=(diffusb(i,sp)+diffus(i,sp))/2
		enddo
		diffush(1)=0.0d0
		ssb(sp)=ss(sp)
		if(flxsah(1).lt.0.d0) then
			ss(sp)=solsur(sp)/ dabs(flxsah(1))
		else
			ss(sp)= 0.0d0
		endif
c		starts of the loop for the solution of the transport equation 
		do i = 1, ncs
c			calculation of the betha's according to the flux direction 
			if(flxsa(i).lt.0.0d0) then
				betha1=1.0d0
				betha4=0.0d0
			else                  
				betha1=0.0d0
				betha4=1.0d0
			endif
			if(flxsa(i+1).lt.0.0d0)then
				betha2=1.0d0
				betha3=0.0d0
			else                  
				betha2=0.0d0
				betha3=1.0d0
			endif
c			coefficients for tridiagonal matrix for i=1
			if (i.eq.1) then
				ab=(dt*pveloh(1)*(wcmo(1)-wcmob(1))*flxsah(1)/
     $				(16.d0*dx*dx))
     $				+(diffush(1)*wcmah(1)/(2.d0*dx*dx))
				bb=(diffush(2)*wcmah(2)/(2.d0*dx*dx))
     $				+(dt*pveloh(1)*(wcmo(1)-wcmob(1))*pvelah(2)/
     $				(16.d0*dx*dx))
				cb=flxsah(1)/(2.d0*dx)
				db=flxsah(2)/(2.d0*dx)
				bdl=(wcmo(1)/dt + ab + bb + 
     $				(parasol(1,1,sp)*parasol(2,1,sp)*parasol(9,1,sp)/dt)
     $				-betha2*db +betha4*cb)
				cdl=(-bb-betha3*db)
				ddl=csolmb(1,sp)*(wcmob(1)/dt-ab-bb+
     $				(parasol(1,1,sp)*parasol(2,1,sp)*parasol(9,1,sp)/dt)
     $				+betha2*db-betha4*cb) + csolmb(2,sp)*(bb
     $				+betha3*db)+ ssb(sp)*(ab-betha1*cb)
     $				+ss(sp)*(ab-betha1*cb)+sinkm(1,sp)
c			coefficients for tridiagonal matrix for i=i
			else if (i.lt.ncs) then
				ab=(diffush(i)*wcmah(i)/(2.d0*dx*dx))
     $				+(dt*pveloh(i)*(wcmo(i)-wcmob(i))*pvelah(i)/
     $				(16.d0*dx*dx))
				bb=(diffush(i+1)*wcmah(i+1)/(2.d0*dx*dx))
     $			+(dt*pveloh(i)*(wcmo(i)-wcmob(i))*pvelah(i+1)/(16.d0*dx*dx))
				cb=flxsah(i )/(2.d0*dx)
				db=flxsah(i+1)/(2.d0*dx)
				adl=(-ab+betha1*cb)
				bdl=(wcmo(i)/dt+ab+bb+
     $				(parasol(1,nla(i),sp)*parasol(2,nla(i),sp)
     $				*parasol(9,nla(i),sp)/dt)
     $				-betha2*db +betha4*cb)
				cdl=(-bb-betha3*db)
				ddl=csolmb(i-1,sp)*(ab-betha1*cb)
     $				+csolmb(i,sp)*(wcmob(i)/dt-ab-bb+
     $				(parasol(1,nla(i),sp)*parasol(2,nla(i),sp)*
     $				parasol(9,nla(i),sp)/dt)
     $				+betha2*db-betha4*cb)
     $				+csolmb(i+1,sp)*(bb+betha3*db)+sinkm(i,sp)
			else
c			coefficients for tridiagonal matrix for i=ncs
				ab=(diffush(ncs)*wcmah(ncs)/(2.d0*dx*dx))
     $			+(dt*pveloh(ncs)*(wcmo(ncs)
     $			-wcmob(ncs))*pvelah(ncs)/(16.d0*dx*dx))
				bb=0.0d0
				cb=flxsah(ncs)/(2.d0*dx)
				db=flxsah(ncs+1)/(2.d0*dx)
				adl=(-ab+betha1*cb)
				bdl=(wcmo(ncs)/dt + ab + 
     $				parasol(1,nla(ncs),sp)*parasol(2,nla(ncs),sp)*
     $				parasol(9,nla(ncs),sp)/dt)
     $				-betha2*db + bb
     $				+betha4*cb
				ddl=csolmb(ncs-1,sp)*(ab-betha1*cb)
     $				+csolmb(ncs,sp)*(wcmob(ncs)/dt-ab+
     $				(parasol(1,nla(ncs),sp)*parasol(2,nla(ncs),sp)*
     $				parasol(9,nla(ncs),sp)/dt)
     $				+betha2*db
     $				-bb-betha4*cb)
     $				+csolmb(ncs+1,sp)*(bb + betha3*db)
			endif
c			solution of the tridiagionale matrix
			if (i.eq.1 )then
				alpha=bdl
			else
				alpha=bdl-adl*betha(i-1)
			endif
			betha(i)=cdl/alpha
			if (i.eq.1) then
				y(1)=ddl/bdl
			else
			y(i)=(ddl-adl*y(i-1))/(bdl-adl*betha(i-1))
			endif
		enddo
c		corrections
		if( dabs(y(ncs)-csolmb(ncs,sp)).gt.0.001d0.and.
     $		flxsah(ncs+1).gt.0.00000000001d0)then
			csolmo(ncs,sp)=((csolmb(ncs+1,sp)*wcmo(ncs+1))+(csolmb(ncs,sp)
     $		*wcmo(ncs)))/(wcmo(ncs+1)+wcmo(ncs))
		elseif(flxsah(ncs+1).gt.0.00000000001d0) then
			csolmo(ncs,sp)=csolmb(ncs+1,sp)
		else
			csolmo(ncs,sp)=y(ncs)
		endif     
		do i=ncs-1,1,-1
			csolmo(i,sp)=y(i)-betha(i)*csolmo(i+1,sp)
		enddo
		do i=1,ncs
			if(csolmo(i,sp).lt.0.0000000000000001d0) then
				csolmo(i,sp)=0.0d0
			endif
		enddo         

c		adsorbed mass 
		do i=1,ncs
			acsolmo(i,sp)=parasol(2,nla(i),sp)*
     $		parasol(9,nla(i),sp)*csolmo(i,sp)
		enddo

c		redistribution mobile/immobile
		do i=1,ncs
			if(imobsw)then
				ag=(wciob(i)+(1-parasol(9,nla(i),sp))*parasol(1,nla(i),sp)
     $				*parasol(2,nla(i),sp))/dt
     $				-parasol(8,nla(i),sp)/2
				bg=(wcio(i)+(1-parasol(9,nla(i),sp))*parasol(1,nla(i),sp)
     $				*parasol(2,nla(i),sp))/dt
     $				+parasol(8,nla(i),sp)/2
				csolio(i,sp)=(ag*csolib(i,sp)+(parasol(8,nla(i),sp)*
     $			csolmo(i,sp))/2.d0 +
     $				(parasol(8,nla(i),sp)*csolmb(i,sp))/2)/bg+sinki(i,sp)/bg
					acsolio(i,sp)=csolio(i,sp)*parasol(2,nla(i),sp)*
     $				(1-parasol(9,nla(i),sp))
				csolo(i,sp)=(csolio(i,sp)*wciob(i)+csolmo(i,sp)*
     $				wcmob(i))/wcob(i)
					acsolo(i,sp)=acsolio(i,sp)+acsolmo(i,sp)
			else
				csolio(i,sp)=0.0d0
				acsolio(i,sp)=0.0d0
				csolo(i,sp)=csolmo(i,sp)
				acsolo(i,sp)=acsolmo(i,sp)
			endif
		enddo

c		solute fluxes
		do i=2,ncs
			if(flxsah(i).lt.0)then
				sflxsa(i,sp)=(flxsah(i)*(csolmo(i-1,sp)+
     $				csolmb(i-1,sp))/2.d0)*dt
			else
				sflxsa(i,sp)=(flxsah(i)*(csolmo(i,sp)+
     $				csolmb(i,sp))/2.d0)*dt
			endif
		enddo
		if (flxsah(1).lt.0.0d0)then
			sflxsa(1,sp) = solsur(sp)*dt
		else
			sflxsa(1,sp)=0.0d0
		endif
		volume= dx
		if( dabs(csolmo(ncs,sp)-csolmb(ncs,sp)).gt.0.001d0.
     $	and.flxsah(ncs+1).gt.0.00000000001d0)then
			sflxsa(ncs+1,sp)=((-csolmb(ncs,sp)+csolmo(ncs,sp))*
     $		wcmo(ncs)*dx)-(acsolmb(ncs,sp)-acsolmo(ncs,sp))*
     $		parasol(1,nla(i),sp)*volume
		else
			if (flxsah(ncs+1).lt.0.0d0)then
				sflxsa(ncs+1,sp)=(flxsah(ncs+1)*
     $			((csolmo(ncs,sp)+ csolmb(ncs,sp))/2))*dt
			else
				sflxsa(ncs+1,sp)=(flxsah(ncs+1)*csolmb(ncs+1,sp))*dt
			endif
		endif
      enddo
      return
      end