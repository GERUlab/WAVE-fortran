c	version 3.0
c	author: m.vanclooster
c###################################################################################
      subroutine wr_cropdata
c     in   : amx, asrqso, eai, ear, eff, harvest_date, ncrop, ncs, 
c            nr_of_tables, nsl, ntabel, plant_date, pr_crop, rgr, rkdf,
c            ssl, stor, table, tbase,
c     out  : -
c     calls: big_title, close_file, date, heading, open_file, small_title
c###################################################################################
      implicit double precision (a-h, o-z)
      include   'constant'
      include   'gen.com'
      include   'crop.com'
      integer res_file, pl_day, pl_month, pl_year, hv_day
      integer hv_year, hv_month
      character *(25) x_out, y_out
      logical table_nec(kt_table,kt_crop)
      data table_nec/
     $.true., .true.,.true. ,.true. ,.true. ,.true. ,.true. ,.false.,
     $.true.,.true. ,. true.,.true.,.true.,.true.,.true.,
     $.true., .true.,.true. ,.true. ,.true. ,.true. ,.true. ,.false.,
     $.false.,.true. ,.true.,.true.,.true.,.true.,.true.,
     $.true., .true.,.true.,.true. ,.true. ,.true. ,.true. ,.true.,
     $.false.,.true.,.false.,.true.,.true.,.true.,.true.,
     $.false.,.true.,.false.,.false.,.false.,.false.,.false.,.false.,
     $.false.,.true.,.false.,.true.,.true.,.true.,.true.,
     $.true., .true.,.false.,.false.,.true. ,.true. ,.true. ,.true. ,
     $.false.,.true.,.true.,.true.,.true.,.true.,.true./

      call open_file (res_file,'cropdata.res')
      call heading (res_file,'cropdata.res',
     $	'input for sucros model')
      call date_wave (plant_date, pl_year, pl_month, pl_day)
      call date_wave (harvest_date, hv_year, hv_month, hv_day)

c	crop phenological parameters 
      write(res_file,10) ncrop,nsl,pl_year, pl_month, pl_day, 
     $     hv_year, hv_month, hv_day, rgr,
     $     tbase,amx, eff,rkdf,scp,rmainso,asrqso,ear,rlat,rlaicr
      if (ncrop.ne.4) write(res_file,20) eai
      if (ncrop.eq.4) write(res_file,40) rmatr
      write(res_file,60) ssl

c     tables
      call big_title (res_file, 'tables')
      do j=1,nr_of_tables
		if (table_nec(j,ncrop)) then
			if (j.eq.1) then
				call small_title (res_file,
     $				'dvs dependent reduction factor of the light 
     $				assimilation rate')
				x_out = 'dvs'
				y_out = 'red_dvs_amax'
				ln_x = 3
				ln_y = 6
			else if (j.eq.2)then
				call small_title (res_file,
     $				'temperature dependent reduction factor of the light
     $				assimilation rate')
				x_out = 'effective temp'
				y_out = 'red_temp_amax'
				ln_x = 14
				ln_y = 6
			else if (j.eq.3)then
				call small_title (res_file,
     $				'development rate as a function of temperature (dvs < 1)')
				x_out = 'temperature'
	            y_out = 'dvr'
				ln_x = 10
				ln_y = 6
			else if (j.eq.4)then
				call small_title (res_file,
     $				'development rate as a function of temperature (dvs > 1)')
				x_out = 'temperature'
				y_out = 'dvr'
				ln_x = 10
				ln_y = 6
			else if (j.eq.5)then
				call small_title (res_file,
     $				'fraction of total dry matter allocated to the shoot')
				x_out = 'dvs'
				y_out = 'fsh'
				ln_x = 3
				ln_y = 8
			else if (j.eq.6)then
				call small_title (res_file,
     $				'fraction of the shoot dry matter allocated to leaves')
				x_out = 'dvs'
				y_out = 'flv'
				ln_x = 3
				ln_y = 8
			else if (j.eq.7)then
				if(ncrop.ne.4) then
					call small_title (res_file,
     $				'fraction of the shoot dry matter allocated to the stems')
					x_out = 'dvs'
					y_out = 'fst'
					ln_x = 3
					ln_y = 8
				endif
			else if (j.eq.8)then
				if (ncrop.eq.3) then
					call small_title (res_file,
     $				'fraction of the shoot dry matter allocated to the cob')
					x_out = 'dvs'
					y_out = 'fcob'
					ln_x = 3
					ln_y = 8
				else if (ncrop.eq.5) then
					call small_title (res_file,
     $				'fraction of total dry matter allocated to fibrous root')
					x_out = 'dvs'
					y_out = 'frt'
					ln_x = 3
					ln_y = 8
				endif
			else if (j.eq.9)then
				if (ncrop.eq.1) then
					call small_title (res_file,
     $				'reduction of the death rate of leaves ')
					x_out = 'dvs'
					y_out = 'weigh_kd_lai_senes'
					ln_x = 3
					ln_y = 8
				endif
			else if (j.eq.10)then
				call small_title (res_file,
     $			'the  specific leave area as a function of dvs')
				x_out = 'dvs'
				y_out = 'sla'
				ln_x = 3
				ln_y = 3
			else if (j.eq.11)then
				if (ncrop.le.2) then
					call small_title (res_file,
     $				'the death rate of the leaves as a function of eff temp')
					x_out = 'effective temp'
					y_out = 'kd_lai_senes'
					ln_x = 14
					ln_y = 12
				elseif (ncrop.eq.5) then
					call small_title (res_file,
     $				'the death rate per degree day as a function of the 
     $				development rate')
					x_out = 'dvs'
					y_out = 'kd_lai_senes_mod'
					ln_x = 3
					ln_y = 16
				endif
			else if (j.eq.12)then
				call small_title (res_file,
     $			'the nitrogen content in the leaves ')
				x_out = 'dvs'
				y_out = 'nitrogen content'
				ln_x = 3
				ln_y = 16
			else if (j.eq.13)then
				call small_title (res_file,
     $			'the nitrogen content in the stems ')
				x_out = 'dvs'
				y_out = 'nitrogen content'
				ln_x = 3
				ln_y = 16
			else if (j.eq.14)then
				call small_title (res_file,
     $			'the nitrogen content in the roots ')
				x_out = 'dvs'
				y_out = 'nitrogen content'
				ln_x = 3
				ln_y = 16
			else if (j.eq.15)then
				call small_title (res_file,
     $			'the nitrogen content in the storage organs ')
				x_out = 'dvs'
				y_out = 'nitrogen content'
				ln_x = 3
				ln_y = 16
			endif
			write(res_file,70) ntabel(j)
			write(res_file,80)x_out(1:ln_x),y_out(1:ln_y)
			do k = 1, ntabel(j)
				write(res_file,90)table(1,k,j),table(2,k,j)
			enddo
		endif
      enddo

c     parameters for root growth
      call big_title(res_file,'parameters for root growth')
      write(res_file,100) specweig
      write(res_file,110) rmrd

c     printing control variables
      call big_title(res_file,'output')
      if (pr_crop_gro)  then
		write (res_file,120)plus
      else
		write (res_file,120)neg
      endif
      if (pr_crop_dev)  then
		write (res_file,130)plus
      else
		write (res_file,130)neg
      endif
      call close_file(res_file)

10    format(//
     $,' crop number              :',i10,/
     $,' number of seedlings      :',i10,/
     $,' emergence date           :',i4,1x,i2,1x,i2,/
     $,' harvest date             :',i4,1x,i2,1x,i2,/
     $,' lai development rate during juvenile stage      :',f10.4,/
     $,' base temperature for crop growth                :',f10.4,/
     $,' co2 assimilation rate at light saturation       :',f10.4,/
     $,' initial light use efficiency                    :',f10.4,/
     $,' light extinction coefficient                    :',f10.4,/
     $,' scattering coefficient                          :',f10.4,/
     $,' maintenance factor for storage organs           :',f10.4,/
     $,' assim. requir. for dm to stor org.(kg/kg)       :',f10.4,/
     $,' ear area growth rate (ha ear/kg/day)            :',f10.7,/
     $,' latitude of the location                        :',f10.4,/
     $,' critical lai for leaf death due to self shading :',f10.4)
20    format(
     $' initial value for ear area index (two sided)    :',f10.4)
40    format(
     $' initial value of the maturity class             :',f10.4)
60    format(
     $' lai of seedling                                 :',f10.7)

70    format('number of elements: ', i10,/)
80    format(a16,1x,a16,/,'---------------- ----------------')
90    format(2x,f12.3,5x,f12.3)
100   format(' specific weight of roots in m/gr= ',f10.4)
110   format(' rmrd= ',f10.4)
120   format(1x,
     $ 'state variables for crop yield                    : ',a)
130   format(1x,
     $ 'state variables that characterise crop development: ',a)
      return
      end