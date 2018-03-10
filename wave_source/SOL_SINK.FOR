c	version 3.0
c	author: m. vanclooster 15/02/02

c###################################################################################
      subroutine sol_sink
c     in   : csolio, csolmo, dt, imobsw, ncs, nla, nr_of_sol, parasol,
c            rates, wcio, wcmo
c     out  : decsoli, decsolm, sinki, sinkm
c     calls: -
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      integer sp

      do sp = 1, nr_of_sol
		do i = 1, ncs -1
			decsolm(i,sp) = -rates(i,sp)*csolmo(i,sp)*wcmo(i)
			decsoli(i,sp) = -rates(i,sp)*csolio(i,sp)*wcio(i)
		enddo
      enddo
	do sp = 1, nr_of_sol
		do i=1,ncs-1
c			the mobile sinkterm
			sinkm(i,sp)=decsolm(i,sp)
     $		-parasol(8,nla(i),sp)*(csolmo(i,sp)-csolio(i,sp))
c			reduction of the mobile sink term
			defsolm=csolmo(i,sp)*(wcmo(i)+parasol(2,nla(i),sp)*
     $		parasol(1,nla(i),sp)*parasol(9,nla(i),sp))/dt+sinkm(i,sp)
			totsolm=-sinkm(i,sp)
			if(defsolm.lt.0.0d0.and.totsolm.gt.0.0d0)then
				facsolm= dmax1(0.d0, 1.d0+defsolm/totsolm)
				sinkm(i,sp)= sinkm(i,sp)*facsolm
				decsolm(i,sp)=decsolm(i,sp)*facsolm
			endif
c			the immobile sinkterm
			if(imobsw) then
				sinki(i,sp)=decsoli(i,sp)
c				reduction of the immobile sink term
				defsoli=csolio(i,sp)*(wcio(i)+parasol(2,nla(i),sp)*
     $			parasol(1,nla(i),sp)* (1-parasol(9,nla(i),sp)))/dt + 
     $			sinki(i,sp)
				totsoli=-sinki(i,sp)
				if(defsoli.lt.0.0d0.and.totsoli.gt.0.0d0)then
					facsoli= dmax1(0.d0, 1.d0+defsoli/totsoli)
					sinki(i,sp)=sinki(i,sp)*facsoli
					decsoli(i,sp)=decsoli(i,sp)*facsoli
				endif
			endif
		enddo
      enddo
      return
      end