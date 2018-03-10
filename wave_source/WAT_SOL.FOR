c	version 3.0
c	author:m vanclooster 19/02/02
c###################################################################################
      subroutine calc_wat_sol
c     in   : flxs flxsa, imobsw, ncs, nla, parasol, pvela, wc, wcio, wcma,
c            wcmo, wco
c     out  : flxsa, flxsah, pvela, pvelah, pveloh, wcio, wciob, wcma, wcmah,
c            wcmob, wcmo, wco, wcob
c     calls: -
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      include   'sol.com'
      dimension flxsab (kt_comps +1 ) ,pvelab(kt_comps +1 )
      dimension pvelo (kt_comps +1 ) ,pvelob(kt_comps +1 )
      dimension wcmab (kt_comps +1 )
      integer sp
      parameter (sp =1)
      logical first_time
      data first_time/.true./

      if (first_time) then
		first_time = .false.
c		calculation of the fluxes and water contents across the nodes
		do i=1,ncs+1
			flxsa(i)=flxs(i)
		enddo
		if(imobsw) then
			do i=1,ncs
				wcmo(i)=parasol(7,nla(i),sp)*wc(i)
				wcio(i)=wc(i)-wcmo(i)
			enddo
		else
			do i=1,ncs
				wcmo(i)=wc(i)
				wcio(i)=0.0d0
			enddo
		endif
		if(flxsa(1).le.0)then
			if(imobsw)then
				wcma(1)=0.999d0
			else
				wcma(1)=1.0d0
			endif
		else 
			if(imobsw)then
				wcma(1)=parasol(7,1,sp)*wc(1)
			else
				wcma(1)=wc(1)
			endif
		endif
		if(imobsw) then
			do i=2,ncs
				wcma(i)=(wcmo(i)+wcmo(i-1))/2
			enddo
			wcma(ncs+1)=wcma(ncs)
		else
			do i=2,ncs
				wcma(i)=wc(i)
			enddo
			wcma(ncs+1)=wcma(ncs)
		endif
		do i=1,ncs
			wco(i)=wc(i)
		enddo
c		calculation of the pore velocity across the nodes
		do i=1,ncs +1
			pvela(i)=flxsa(i)/wcma(i)
		enddo
		do i=1,ncs
			flxso=(flxsa(i)+flxsa(i+1))/2
			pvelo(i)=flxso/wcmo(i)
		enddo
      endif

c     calculation of the fluxes and water contents across the nodes
      do i = 1,ncs+1
		flxsab(i) = flxsa(i)
		flxsa(i)=flxs(i)
      enddo
      do i = 1, ncs
		wcmob(i) = wcmo(i)
		wciob(i) = wcio(i)
      enddo
      if(imobsw) then
		do i=1,ncs
			wcmo(i)=parasol(7,nla(i),sp)*wc(i)
			wcio(i)=wc(i)-wcmo(i)
		enddo
      else
		do i=1,ncs
			wcmo(i)=wc(i)
			wcio(i)=0.0d0
		enddo
      endif
      do i = 1,ncs +1
		wcmab(i) = wcma(i)
      enddo
      if(flxsa(1).le.0)then
		if(imobsw)then
			wcma(1)=0.999
		else
			wcma(1)=1.0d0
		endif
      else 
		if(imobsw)then
			wcma(1)=parasol(7,1,sp)*wc(1)
		else
			wcma(1)=wc(1)
		endif
      endif
      if(imobsw) then
		do i=2,ncs
			wcma(i)=(wcmo(i)+wcmo(i-1))/2
		enddo
		wcma(ncs+1)=wcma(ncs)
      else
		do i=2,ncs
			wcma(i)=wc(i)
		enddo
		wcma(ncs+1)=wcma(ncs)
      endif
      do i=1,ncs
		wcob(i) = wco(i)
		wco(i)=wc(i)
      enddo

c     calculation of the pore velocity across the nodes
      do i=1,ncs +1
		pvelab(i) = pvela(i)
		pvela(i)=flxsa(i)/wcma(i)
      enddo
      do i=1,ncs
		pvelob(i) = pvelo(i)
		flxso=(flxsa(i)+flxsa(i+1))/2
		pvelo(i)=flxso/wcmo(i)
      enddo
      do i=1,ncs+1
		flxsah(i)=(flxsab(i)+flxsa(i))/2
		wcmah(i)=(wcmab(i)+wcma(i))/2
		pvelah(i)=(pvelab(i)+pvela(i))/2
      enddo
      do i=1,ncs
		pveloh(i)=(pvelob(i)+pvelo(i))/2
      enddo

c     set the upper boundary condition
      if (flxsa(1).gt.0.0) then
		flxsa(1)=0.0
		flxsab(1)=0.0
		flxsah(1)=0.0
      endif
      end