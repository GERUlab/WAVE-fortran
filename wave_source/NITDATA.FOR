c	version 3.0
c	author:m.vanclooster 26/2/2002
c###################################################################################
      subroutine nitdata
c     in   : isd, isucr, ncs, simplant, yearin
c     out  : alfa_rdens, carborg, cchumo, cclito, ccmano, cnhumo, cnlito,
c            cnmano, fag, fe, fh, fhp, flr, g, idman, nof, pr_nit_cum,
c            pr_nit_sum, pr_nit_tseries, pr_temp_for_nit, rdo, rkdenit,
c            rkhum, rkhyd, rklit, rkman, rknitri, rkvol, rnitorg, rnmaxp,
c            ro, rorad, tcorgs, tcorgs_isd, tnorgs, tnorgs_isd, w0_rdens
c     calls: check_et, end_file, hold_line, init_file, jdate, next_line, rd_date, 
c            rd_int, rd_logic, rd_real, report_err
c###################################################################################
      implicit double precision (a-h, o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      include   'nit.com'
      integer applicat
      integer day, month
      integer year_for
      integer rd_int
      logical rd_logic,rd_ext, warning_reported
      dimension test_min(kt_in_row-1),test_max(kt_in_row-1)
      dimension result (kt_in_row, 366)

      call init_file ('nitdata.in')
      rd_ext =rd_logic()
      if (rd_ext) call init_file('nitubc.nu')
c	crop residues 
      if (simplant) then  
		fag = rd_real(0.d0,1.d0)
		flr = rd_real(0.d0,1.d0)
		fhp = rd_real(0.d0,1.d0)
      endif
c     organic fertilizer
      nof = rd_int(0,kt_nof)
      do i = 1,nof
		call hold_line
		applicat = rd_int (1,kt_nof)
		if (applicat.ne.i) then
			call report_err
     $			("error in sequence of fertilizer applications")
		else
			call rd_date (month, day)
			idman(i) = jdate(year_for(month,day),month,day,yearin)
			carborg(i) = rd_real (0.d0, mx_real)
			rnitorg(i) = rd_real (0.d0, mx_real)
		endif
		call next_line
      enddo
      if (nof.gt. 0) call check_et
      if (rd_ext) call end_file
c     n sink term
      rd_ext =rd_logic()
      if (rd_ext) call init_file('nitsink.ns')
c     organic turn-over
      ro = rd_real(mn_real,mx_real)
      fe = rd_real(-mx_real,mx_real)
      fh = rd_real(-mx_real,mx_real)
c     uptake
      if (simplant) then
		rorad = rd_real(0.d0 ,mx_real)
		rdo = rd_real(mn_real, mx_real)
		if (.not.isucr) then
			g = rd_real(0.d0,1.d0)
			rnmaxp = rd_real(0.d0,mx_real)
			w0_rdens = rd_real(0.d0,mx_real)
		endif
      endif
      if (rd_ext) call end_file
c     n rate constants decomposition
      rd_ext =rd_logic()
      if (rd_ext) call init_file('nitdecom.nd')
	do i = 1, 7
		test_min(i) = 0.d0
		test_max(i) = mx_real
	enddo
      call interp_profile (result, ncs, 7, 
     $	'decomposition rates',test_min, test_max,9)
	do i = 1,ncs
		rknitri(i) = result(1,i)
		rkdenit(i) = result(2,i)
		rkhyd(i)   = result(3,i)
		rkvol(i) =  result(4,i)
		rklit(i) =  result(5,i)
		rkman(i) =  result(6,i)
		rkhum(i) =  result(7,i)
	enddo
	if (rd_ext) call end_file
c	checks
      warning_reported = .false.
      do sp = 1, nr_of_sol
		do i = 1, ncs
			if (rates(i,sp).ne.0.d0) then
				if (.not.warning_reported) then
					call warn_err 
     $					('nitrates: changed the rates from soldata.in to 0')
				endif
				rates(i,sp) = 0.d0
			endif
		enddo
      enddo
      if( rknitri(ncs).ne.0.d0 .or. rkdenit(ncs).ne.0.d0 .or.
     $	rkhyd(ncs).ne.0.d0 .or. rkvol(ncs).ne.0.d0 .or.
     $	 rklit(ncs).ne.0.d0.or. rkman(ncs).ne.0.d0 .or.
     $	rkhum(ncs).ne.0.d0)
     $		call warn_err
     $			('nitrates: transformation rates for bottom comp. not 0')
c     initial values
      rd_ext =rd_logic()
      if (rd_ext) call init_file('nitinit.ni')
	do  i = 1, 6
		test_min(i) = mn_real
		test_max(i) = mx_real
	enddo
      call interp_profile (result, ncs, 6, 
     $	'initial values organic fraction'
     $	,test_min, test_max,9)
	do i = 1,ncs
		ccmano(i)=  result(1,i)
		cnmano(i)=  result(2,i)
		cclito(i)=  result(3,i)
		cnlito(i)=  result(4,i)
		cchumo(i)=  result(5,i)
		cnhumo(i)=  result(6,i)
	enddo
	if (rd_ext) call end_file 
	tnorgs = 0.d0
	tcorgs = 0.d0
	do i = 1, ncs
		tcorgs = tcorgs + ccmano(i) + cclito(i) + cchumo(i)
		tnorgs = tnorgs + cnmano(i) + cnlito(i) + cnhumo(i)
		if (i.eq.isd) then
			tcorgs_isd = tcorgs
			tnorgs_isd = tnorgs
		endif
	enddo	 
c     print control
c     file with overview of main state variables (nit_sum.out)
      pr_nit_sum = rd_logic()
c     file with cumulative terms of the organic pool (nit_cum.out)
      pr_nit_cum = rd_logic()
c     file with n in organic litter  (orgnlit.out)
      pr_nit_tseries(1) = rd_logic()
c     file with n in organic manure  (orgnman.out)
      pr_nit_tseries(2) = rd_logic()
c     file with n in organic humus  (orgnhum.out)
      pr_nit_tseries(3) = rd_logic()
c     file with c in organic litter  (orgclit.out)
      pr_nit_tseries(4) = rd_logic()
c     file with c in organic manure  (orgcman.out)
      pr_nit_tseries(5) = rd_logic()
c     file with c in organic humus  (orgchum.out)
      pr_nit_tseries(6) = rd_logic()
c     file with nh4 uptake rates (nh4upt.out)
      pr_nit_tseries(7) = rd_logic()
c     file with no3 uptake rates (no3upt.out)
      pr_nit_tseries(8) = rd_logic()
c     file with denitrification rates (denitrif.out)
      pr_nit_tseries(9) = rd_logic()
c     file with nitrification rates (nitrific.out)
      pr_nit_tseries(10) = rd_logic()
c     file with hydrolysis rates (hydrolys.out)
      pr_nit_tseries(11) = rd_logic()
c     file with volatilisation rates (volatil.out)
      pr_nit_tseries(12) = rd_logic()
c     file with nh4 mineralisation rates (nh4min.out)
      pr_nit_tseries(13) = rd_logic()
c     file with no3 mineralisation rates (no3min.out)
      pr_nit_tseries(14) = rd_logic()
c     file with cumulative convective nitrate uptake rates (no3_cupt.out)
      pr_nit_tseries(15) = rd_logic()
c     file with cumulative diffusive nitrate uptake rates  (no3_dupt.out)
      pr_nit_tseries(16) = rd_logic()

      call end_file
      end