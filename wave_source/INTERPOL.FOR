c	version 3.0
c	author m. vanclooster 15/02/2002
c######################################################################################
      subroutine interp_profile (result, nr_rows, nr_cols, table_name
     $ ,test_min, test_max,nine)
c######################################################################################
	implicit double precision (a-h, o-z)
      include 'constant'
      dimension result (kt_in_row, kt_comps+1)
      dimension komp (kt_comps+1)
      dimension test_min(kt_in_row-1),test_max(kt_in_row-1)
      double precision obs(kt_in_row, kt_comps+1)
      character*(20) table_name 
      logical already_read(kt_comps+1) 
      logical eot
      logical error_occurred
      integer rd_int
      if (nr_rows.gt.kt_comps+1.or.nr_cols.gt.kt_in_row-1) then
		call report_err ('too many rows or columns in table ' //
     $	table_name)
		call stop_simulation('programme stopped: programming error')
      endif
      if (nine.ne.9) then
		call stop_simulation ('interp_profile: programming error')
      endif
      do i = 1, nr_rows
		already_read(i) = .false.
      enddo

      error_occurred = .false.
      nr_obs = 0
5     call hold_line

      if (.not.eot().and.nr_obs.le.nr_rows) then
		lyr = rd_int(1,nr_rows)
		if (already_read(lyr)) then
			call report_err ('there is a duplicate row in table ' //
     $		table_name)
			error_occurred = .true.
		else
			nr_obs = nr_obs + 1
			already_read(lyr) = .true.
		endif
		komp(nr_obs) = lyr
		if (nr_obs.gt.1) then 
			if (komp(nr_obs).lt.komp(nr_obs-1)) then
				call report_err ('the rows are not ascending in table ' //
     $			table_name)
				error_occurred = .true.
			endif
		endif
		do j = 1, nr_cols
			obs(j, nr_obs) = rd_real(test_min(j),test_max(j))
		enddo
		call next_line
		goto 5
      endif
      call next_line
      if (error_occurred) return
      if (nr_obs.lt.nr_rows) then
		call lin_interpol (result,komp,obs,nr_obs,nr_rows,nr_cols,table_name)
      else
		do i = 1, nr_rows
			do j = 1, nr_cols
				result(j,i) = obs(j,i)
			enddo
		enddo
      endif
      return
      end

c######################################################################################
      subroutine interp_tseries (result,nr_cols, table_name
     $ ,test_min, test_max,nine)
c######################################################################################
      implicit double precision (a-h, o-z)
      include 'constant'
      include 'gen.com'
      dimension result (kt_in_row, 366)
      integer d_nr (366), day_nr_no_yr , day
      dimension test_min(kt_in_row-1),test_max(kt_in_row-1)
      double precision obs(kt_in_row, 366)
      character*(20) table_name 
      logical already_read(366)
      logical eot
      logical error_occurred
      integer rd_int

      if (nr_cols.gt.kt_in_row-1) then
		call report_err ('too many columns in table ' //table_name)
		call stop_simulation('programme stopped: programming error')
      endif
      if (nine.ne.9) then
		call stop_simulation ('interp_profile: programming error')
      endif
      do i = 1, ndsim
		already_read(i) = .false.
      enddo
      error_occurred = .false.
      nr_obs = 0
5     call hold_line
      if (.not.eot().and.nr_obs.lt.ndsim) then
		month = rd_int (1,12)
		day = rd_int (1,31)
		nr_day =  day_nr_no_yr (month,day)
		if (already_read(nr_day)) then
			call report_err ('there is a duplicate row in table ' //
     $		table_name)
			error_occurred = .true.
		else
			nr_obs = nr_obs + 1
			already_read(nr_day) = .true.
		endif
		d_nr(nr_obs) = nr_day
		if (nr_obs.gt.1) then 
			if (d_nr(nr_obs).lt.d_nr(nr_obs-1)) then
				call report_err ('the rows are not ascending in table ' //
     $			table_name)
				call report_err ('the rows are not ascending in table ')
				error_occurred = .true.
			endif
		endif
		do j = 1, nr_cols
			obs(j, nr_obs) = rd_real(test_min(j),test_max(j))
		enddo
		call next_line
		goto 5
      endif
      call next_line
      if (error_occurred) return
      if (nr_obs.lt.ndsim) then
		call lin_interpol (result,d_nr,obs,nr_obs, ndsim, nr_cols,table_name)
      else
		do i = 1, ndsim
			do j = 1, nr_cols
				result(j,i) = obs(j,i)
			enddo
		enddo
      endif
      return
      end

c######################################################################################
      subroutine lin_interpol(result, identif, obs, 
     $ nr_obs, nr_rows, nr_cols,table_name)
c######################################################################################
      implicit double precision (a-h, o-z)
      include 'constant'
      dimension result (kt_in_row, 366)
      dimension identif (366)
      double precision obs(kt_in_row, 366)
      character*(20) table_name 
      integer index

      if (nr_obs.lt.1) then
		call report_err ('no observations in table ' //table_name)
		return
      endif
      if (nr_obs.gt.366) then
		call report_err
     $	('too many observations in table'//table_name//'(>366)')
		return
      endif
      do index = 1, identif(1)
		do j = 1, nr_cols
			result(j,index) = obs(j,1)
		enddo
      enddo
      do i = 2, nr_obs
10		if (identif(i).ge.index) then
			do j = 1, nr_cols
				result(j,index) = obs(j,i-1)+ dfloat(index - identif(i-1)) *
     $			(obs(j,i)-obs(j,i-1))/dfloat(identif(i)-identif(i-1))
			enddo
			index = index + 1
			goto 10
		endif   
      enddo
      do i = identif(nr_obs)+1, nr_rows
		do j = 1, nr_cols
			result(j,i) = obs(j,nr_obs)
		enddo
      enddo
      end