c	Version 3.0
c	author: m. vanclooster 15/02/2002
c###################################################################################
      subroutine soiltemp
c     in   : dt, dx, ncs, nday, nla, t, tmax, tmin, wcto, wctob
c     out  : temp
c     calls: report_err, stop_simulation
c###################################################################################
      implicit double precision  (a-h,o-z)
      include  'constant'
      include   'gen.com'
      include   'clim.com'
      include   'temp.com'
      parameter ( pi=3.141592654)
      logical first_time
      dimension f1(kt_comps),g1(kt_comps),beta(kt_comps),aka(4),allam(kt_comps)
      dimension part(5),alamo(5),geoma(3),geomq(3),geomos(3),geomom(3)
c     1=quartz; 2=organic matter; 3=other solids; 4=air; 5=water
      data part/0.54, 0.045,0.015,0.0,0.0/
c     alamo mcal/cm/s/øc      
      data alamo/20.4, 0.60, 7.0, 0.0615, 0.176/
      data geomq/0.125, 0.125, 0.75/
      data geomos/0.125, 0.125, 0.75/
      data geomom/0.5, 0.5, 0./
      data first_time/.true./

      if (first_time) then
		first_time = .false.
		temp(1) = (tmin(nday)+tmax(nday))/2.d0
		temp(ncs) = 22.0d0
		if (ncs.eq.0) then
			call report_err ('number of compartments is zero')
			call stop_simulation ('program stopped: check err_file')
		endif
		grad = (temp(1) - temp(ncs))/ncs
		do  i = 1, ncs
			temp(i)= temp(1) - grad*i
		enddo
      endif
      do i=1,ncs
		if(i.eq.1) then
			tw=(wcto(i)+wctob(i))/2.
		elseif(i.gt.1 .and. i.le.(ncs)) then
			tw=(wcto(i)+wctob(i)+wcto(i-1)+wctob(i-1))/4.
		endif
		beta(i)=0.00084*bulk_dens(nla(i))+0.0041868*tw
		if (beta(i).eq.0.d0) then
			call report_err ('beta is zero')
			call stop_simulation ('program stopped: check err_file')
		endif
		if(wcto(i).ge..20)then
			alama=alamo(5)
		else
			alama=alamo(4)
		endif
		if (alama.eq.0.d0) then
			call report_err ('alama is zero')
			call stop_simulation ('program stopped: check err_file')
		endif
c	    determining k values used in calculating the thermal conductivity
c		of the soil
		aka(1)=0.d0
		do j=1,3
			aka(1)=(1.d0/(1.d0+((alamo(1)/alama)-1.d0)*geomq(j)))+aka(1)
		enddo
		aka(1)=aka(1)/3.d0
		aka(2)=0.d0
		do j=1,3
			aka(2)=(1.d0/(1.d0+((alamo(2)/alama)-1.d0)*geomom(j)))+aka(2)
		enddo 
		aka(2)=aka(2)/3.d0
		aka(3)=0.d0
		do j=1,3
			aka(3)=(1.d0/(1.d0+((alamo(3)/alama)-1.d0)*geomos(j)))+aka(3)
		enddo 
          aka(3)=aka(3)/3.d0
          if(wcto(i).ge. 0.20d0) then
			geoma(1)=(wcto(i)-.0845d0)/1.10d0
          else
			geoma(1)=(wcto(i) - .0333d0)/2.22d0
          endif
          geoma(2)=geoma(1)
          geoma(3)=1.d0 - geoma(1) - geoma(2)
          if (wcto(i).ge. 0.20d0) then
			alamv=.176d0
          else
			alamv=wcto(i)/.880d0
          endif
		aka(4)=0.d0
          do m=1,3
			aka(4)=(1.d0/(1.d0+(((alamo(4)+alamv)/alama)-1.d0)*geoma(m)))
     $		+aka(4)
		enddo 
          aka(4)=aka(4)/3.d0
          a12=0.d0
          a13=0.d0
          part(4)=1.d0-wcto(i)-part(1)-part(2)-part(3)
          do mm=1,4
			a12=aka(mm)*part(mm)*alamo(mm)+a12
			a13=aka(mm)*part(mm)+a13
		enddo
          if (wcto(i).ge. 0.20d0)then
			fx=tw
          else
			fx=1.d0-tw-part(1)-part(2)-part(3)
          endif
          a12=a12+(alama*fx)
          a13=a13+fx
c         calculating thermal conductivities at half depth, and half time steps
          if(wcto(i).le. 0.20d0)then
			allam(i)=(a12/a13)*1.65d0*3600.d0*24.d0*4.1868d0/
     $		(1000.d0*10.d0)
          else
			allam(i)=(a12/a13)*3600.d0*24.d0*4.1868d0/(1000.d0*10.d0)
          endif
	enddo
c     top boundary condition 
      temb=(tmax(nday)+tmin(nday))/2.d0+(tmax(nday)-tmin(nday))/2.d0
     $* dsin(2.d0*pi*(t- dint(t)-0.25d0))
c     calculating coefficients for the tridiagonal matrix.
      do i = 1, ncs-1
		aa=(dt*allam(i))/(beta(i)*2.d0*dx*dx)
          bb=(dt*allam(i+1))/(beta(i+1)*2.d0*dx*dx)
          a1=-aa
          b1=aa+bb+1.d0
          c1=-bb
          if (i.eq.1) then
			d1=aa*temb+(1.-aa-bb)*temp(1)+bb*temp(2)
			f1(1)=c1/b1
			g1(1)=(d1-a1*temb)/b1
          else
			d1=aa*temp(i-1)+(1.d0-aa-bb)*temp(i)+bb*temp(i+1)
			f1(i)=c1/(b1-f1(i-1)*a1)
			g1(i)=(d1-a1*g1(i-1))/(b1-a1*f1(i-1))
          endif
      enddo 
c     bottom boundary condition (zero gradient)
      temp(ncs)=temp(ncs-1)
c     backward substitution      
      do i = ncs -1, 1,-1
		temp(i)=g1(i)-f1(i)*temp(i+1)
      enddo
      return
      end

c###################################################################################
      subroutine what_is_soiltemp (temp_out)
c###################################################################################
      implicit double precision  (a-h,o-z)
      include  'constant'
      include   'gen.com'
      include   'clim.com'
      include   'temp.com'
      dimension temp_out (kt_comps)

      do i= 1, ncs
		temp_out(i) = temp(i)
	enddo 
      return
      end
