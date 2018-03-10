c	Version 3.0
c	m. vanclooster 26/02/2002
c###################################################################################
      subroutine cropdata
c     in   : ncs, simnit
c     out  : amx, asrqso, eai, ear, eff, ncrop, nr_of_tables, nsl, ntabel,
c            scp, sla, specweig, ssl, table, tbase, tend
c     calls: check_et, end_file, hold_line, init_file, next_line, rd_int,
c            rd_logic, rd_real, report_err, sum
c###################################################################################
      implicit double precision (a-h, o-z)
      include   'constant'
      include   'gen.com'
      include   'crop.com'

      character *180 message
      integer rd_int
      logical rd_logic,table_nec(kt_table,kt_crop)
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

      call init_file('cropdata.in')

      ncrop = rd_int (1, kt_crop)
      nsl = rd_int (1, mx_int)
      rgr = rd_real(mn_real, mx_real)
      tbase= rd_real(-50.d0, 100.d0)
      amx= rd_real(mn_real, mx_real)
      eff= rd_real(-mx_real, mx_real)
      rkdf= rd_real(-mx_real, mx_real)
      scp= rd_real(-mx_real, mx_real)
      rmainso= rd_real(-mx_real, mx_real)
      asrqso= rd_real(-mx_real, mx_real)

      if(ncrop.le.2) ear= rd_real(-mx_real, mx_real)
      rlat= rd_real(-90.d0, 90.d0)
      rlaicr= rd_real (mn_real, mx_real)
      if(ncrop.le.2) eai= rd_real(0.d0, mx_real)
      if(ncrop.eq.4) rmatr=rd_real(0.d0, mx_real)
      ssl= rd_real(0.d0, mx_real)
      if (ncrop.eq.1.and.tend.lt.366.d0) then
		call report_err 
     $	('simulation of winter wheat should start before new year ')
      endif
      nr_of_err = 0
      if (simnit) then
		nr_of_tables = 15
      else
		nr_of_tables = 11
      endif
      do i=1,nr_of_tables
		if (table_nec(i,ncrop))then
			ntabel(i)= rd_int (0,kt_cropinfo)
			do j=1,ntabel(i)
				call hold_line
				do k = 1,2
					table(k,j,i) = rd_real(- mx_real, mx_real)
				enddo
				call next_line
			enddo
			call check_et
			do j=2,ntabel(i)
				if(table(1,j-1,i).ge.table(1,j,i)) nr_of_err=nr_of_err+1
			enddo
		endif
      enddo 
      if(nr_of_err.ne.0) then
		write (message,10)nr_of_err
		call report_err(message)
      end if
      specweig = rd_real (0.d0, mx_real)
      rmrd = rd_real(-mx_real, mx_real)
      pr_crop_gro = rd_logic()
      pr_crop_dev = rd_logic()
      call end_file
10	format('data for x values of crop growth tables are not', 
     $         ' in increasing order.number of errors :', i5 )
      return
      end