c	version 3.0
c	author m. vanclooster 15/02/2002
c	corrections: m. javaux 03/04/02 
c######################################################################################
      subroutine check_et
c######################################################################################
      include   'in_sys.com'
      character *(ms_length) message
      integer excess
      if (file_end) then
		end_of_table = .false.
		return
      endif
      excess = 0
      if (.not.end_of_table) call read_line
      if (.not.end_of_table) then
		call report_err ('missed end of table (ET)')
      end if
10    if (.not.end_of_table.and..not.file_end) then
		excess = excess + 1
		call read_line 
		goto 10
      endif
      if (excess.gt.0) then 
		write(message,20)  excess , 'excess lines'
		call report_err(message) 
20		format (i3,1x,a15)
      endif
      end_of_table = .false.
      end

c######################################################################################
      logical function eot ()
c######################################################################################
      include  'in_sys.com'
      eot = end_of_table
      if (eot) call check_et
      return
      end

      logical function input_error ()
      include  'in_sys.com'
      input_error = error_occured
      return
      end

c######################################################################################
      subroutine end_file
c######################################################################################
      include   'in_sys.com'
      character *(ms_length) message
      close (unit=input_unit(nr_of_files), err= 10)
      call free_unit(input_unit(nr_of_files))
      nr_of_files = nr_of_files - 1
      file_end = nr_of_files.lt.1
      return
10    write (message,20) fname(nr_of_files)
20	format ('unable to close file: ', a)
      call report_err (message)
      call free_unit(input_unit(nr_of_files))
      nr_of_files = nr_of_files - 1
      file_end = nr_of_files.lt.1
      return
      end

c######################################################################################
      subroutine hold_line
c######################################################################################
      include   'in_sys.com'
      call read_line
      hold = .true.
      end

c######################################################################################
      subroutine init_err
c######################################################################################
      include   'in_sys.com'
      logical first_time
      integer err_unit
      common /err_sys/ first_time,err_unit
      call get_unit(err_unit)
      open (unit=err_unit, file= 'err_file.err', status = 'new', err= 10)
      close (err_unit)
      return
10    open (unit=err_unit, file= 'err_file.err', status = 'old', err= 20)
      close (err_unit, status = 'delete')
      open (unit=err_unit, file= 'err_file.err', status = 'new', err= 20)
      close (err_unit)
      return 
20    write (*,*) 'unable to open the errorfile'
      call stop_simulation ('programme stopped: error during input')
      return
      end

c######################################################################################
      subroutine init_file (filename)
c######################################################################################
      include   'in_sys.com'
      character* (*) filename
      character *(ms_length) message
      logical first_time
      integer err_unit
      common /err_sys/ first_time,err_unit
      if (first_time) then
		first_time =.false.
		call init_err
      endif
      nr_of_files = nr_of_files + 1
      fname(nr_of_files) = filename
      nr_of_lines(nr_of_files) = 0
      hold = .false. 
      end_of_table = .false. 
      report_err_for_line = .true.
      call get_unit (input_unit(nr_of_files))
      open (unit=input_unit(nr_of_files),status ='old',file= filename, err= 10)
      file_end = .false.
      return
10    write (message,20) fname(nr_of_files)
20    format ('file does not exist: ',a)
      call report_err (message)
      file_end =  .true.
      end

c######################################################################################
      subroutine next_line
c######################################################################################
      include   'in_sys.com'
      hold = .false.
      report_err_for_line = .true.
      end

c######################################################################################
      subroutine rd_date (month,day)
c######################################################################################
      include   'in_sys.com'
      integer month,day, rd_int
      integer ndm(12)
      data ndm/31,29,31,30,31,30,31,31,30,31,30,31/
      month = rd_int (1, 12)
      day = rd_int (1,ndm(month))
      return
      end

c######################################################################################
      integer function rd_int (min,max)
c######################################################################################
      include   'in_sys.com'
      character *(ms_length) message
      character *(ln_length) line_rest
      integer min, max,res
      if (file_end)then
		rd_int = min
		return
      endif
      if(end_of_table) then 
		if (report_err_for_line) call report_err ('not enough lines in input')
		if (hold) report_err_for_line = .false.
		rd_int = min
		return
      end if
      if (.not. hold) then
		call read_line
		if (line.eq.empty_line) then
			call report_err('not enough lines in input')
			rd_int = min
          return
		endif
      else
		call extract_one(line,line_rest)
		if (line.eq.empty_line) then
			call report_err('not enough values on line')
			rd_int = min
          return
		endif
      endif
      read (line, fmt=20,err=10) res
20    format (i10)
      if (res.lt.min.or.res.gt.max) then
		if (check_range) then
			write(message,30) res, min,max
			call report_err (message)
			rd_int = min
		else
			write(message,30) res, min,max
30			format (i10, ' is not between ',i10,' and ', i10)
			call warn_err (message)
			rd_int = res
		endif
      else
		rd_int = res
      endif
      if (hold) line = line_rest
      return
10    call report_err ('error in this line (integer expected): '//line)   
      rd_int = min 
      return
      end

c######################################################################################
      logical function rd_logic ()
c######################################################################################
      include   'in_sys.com'
      character *(ln_length) line_rest
      integer i
      if (file_end)then
		rd_logic = .false.
		return
      endif
      if(end_of_table) then 
		if (report_err_for_line) call report_err ('not enough lines in input')
		if (hold) report_err_for_line = .false.
		rd_logic = .false.
		return
      endif
      if (.not. hold) then
		call read_line
		if (line.eq.empty_line) then
			if (report_err_for_line) call report_err
     $		('not enough lines in input')
			if (hold) report_err_for_line = .false.
			rd_logic =  .false.
			return
		endif          
      else
		call extract_one (line,line_rest)
		if (line.eq.empty_line) then
			call report_err ('not enough values on line')
			rd_logic = .false.
			return
		endif
      endif
      i=1
10    if (line(i:i).eq. ' ') then
		i=i+1
		goto 10
      endif
      if (line(i:i).eq.'y'.or.line(i:i).eq.'Y') then
		rd_logic = .true.
      elseif (line(i:i).eq.'n'.or.line(i:i).eq.'N') then
		rd_logic = .false.
      else
		call report_err ('y,Y,n,N expected on this line: '//line)
		rd_logic = .false. 
      endif
      if (hold) line = line_rest
      return
      end 

c######################################################################################
      double precision function rd_real (min,max)
c######################################################################################
      include   'in_sys.com'
      character *(ln_length) line_rest
      character *(ms_length) message
      double precision min, max,res
      if (file_end)then
		rd_real = min
		return
      endif
      if(end_of_table) then 
		if (report_err_for_line) call report_err ('not enough lines in input')
		if (hold) report_err_for_line = .false.
		rd_real = min
		return
      endif
      if (.not. hold) then
		call read_line
		if (line.eq.empty_line) then
			call report_err ('not enough lines in input')
			rd_real = min
			return
		endif
      else
		call extract_one(line,line_rest)
		if (line.eq.empty_line) then
			call report_err ('not enough values on line')
			rd_real = min
			return
		endif
      endif
      read (line, fmt=20,err=10) res
20    format(f12.0)
      if (res.lt.min.or.res.gt.max) then
		if (check_range) then
			write(message,30) res, min,max
			call report_err (message)
			rd_real = min
		else
			write(message,30) res, min,max
30			format (e12.3, ' is not between ',e12.3,' and ', e12.3)
			call warn_err (message)
			rd_real = res
		endif
      else
		rd_real = res
      endif  
      if (hold) line = line_rest
      return
10    call report_err ('error in this line (real expected): '//line)   
      rd_real = min 
      return
      end

c######################################################################################
      subroutine read_line
c######################################################################################
      include   'in_sys.com'
      character *(ln_length)  ln_rd
      character *(ms_length) message
      if (file_end) 	return
      if (end_of_table) then
		line = empty_line
		return
      endif
10    read (input_unit(nr_of_files), fmt = 30 ,end= 25, err = 20) ln_rd
      nr_of_lines(nr_of_files) = nr_of_lines(nr_of_files) + 1
      if (ln_rd .eq.empty_line ) goto 10
      call rm_comment(ln_rd)
      if (ln_rd(1:1).eq.'C'.or.ln_rd(1:1).eq.'c'.or.ln_rd.eq.empty_line) goto 10
      if (ln_rd(1:1).eq.'E' .or. ln_rd(1:1).eq.'e' )  then
		end_of_table = .true.
		ln_rd = empty_line
      endif
	line = ln_rd
      return
20    line = empty_line 
      call report_err ('error reading line from external device')
      return
25    line = empty_line 
      write (message,40)fname(nr_of_files)
40    format('reached end of file: ',a)
      call report_err (message)
      file_end = .true.
      return
30    format (a120)
	write (*,*) 'read_line336'
      end

c######################################################################################
      subroutine report_err (message)
c######################################################################################
      include   'in_sys.com'
      character *(*) message
      logical first_time
      integer err_unit
      common /err_sys/ first_time,err_unit
      if (first_time) then
		first_time =.false.
		call init_err
      endif
      open (unit=err_unit, file= 'err_file.err', status = 'old',
     $access = 'append', err= 10)
      if (nr_of_files.gt.0) then
		error_occured = .true.
		if (file_end) then 
			write(err_unit,*,err=10)  message
			write(err_unit,*)
		else 
			write(err_unit,*,err=10) 'problems with file: ',fname(nr_of_files) 
			if (nr_of_lines(nr_of_files).ne. 0)
     $		write(err_unit,fmt =40, err=10) nr_of_lines(nr_of_files)
			write(err_unit,*,err=10)  message
			write(err_unit,*)
		endif
	endif
	close(err_unit)
	return
10	close (err_unit, err=20)
20	write (*,*) ' cannot write to the errorfile'
40	format (' error on line: ', i4)
      end

c######################################################################################
      subroutine warn_err (message)
c######################################################################################
      include   'in_sys.com'
      character *(*) message
      logical first_time
      integer err_unit
      common /err_sys/ first_time,err_unit
      if (first_time) then
		first_time =.false.
		call init_err
      endif
      open (unit=err_unit, file= 'err_file.err', status = 'old',
     $access = 'append', err= 10)
      if (nr_of_files.gt.0) then
		if (.not.file_end) then
			write(err_unit,*,err=10) 'warning for file: ',fname(nr_of_files) 
			if (nr_of_lines(nr_of_files).ne. 0)
     $		write(err_unit,fmt =40, err=10) nr_of_lines(nr_of_files)
			write(err_unit,*,err=10)  message
			write(err_unit,*)
		endif
      else
		write(err_unit,*,err=10)  message
		write(err_unit,*)
      endif
      close(err_unit)
      return
10    close (err_unit, err=20)
20    write (*,*) ' cannot write to the errorfile'
40    format (' possible error on line: ', i4)
      end

      subroutine rm_comment (line)
      integer first, ln_line
      logical space
      character *(*) line
      character empty_line
      parameter (empty_line=' ')
      ln_line = len(line)
      first = 1
10    if (line(first:first).ne.':') then
		first= first+1
		if (first.le.ln_line)goto 10
      endif
      if (first.gt.ln_line) return
20    if (first.lt.ln_line) then
		first = first +1
		if (space(line(first:first))) goto 20
      endif
      if(space(line(first:first))) then
		line = empty_line
      else
		line (1: ln_line) = line(first:ln_line)
      endif
      return
      end

c######################################################################################
      subroutine set_range_checking
c######################################################################################
      include  'in_sys.com'
      check_range = .true.
      return
      end

c######################################################################################
      logical function space (rd_char)
c######################################################################################
      character rd_char
      space =ichar(rd_char).lt.39
      return
      end

c######################################################################################
      subroutine extract_one(head,tail)
c######################################################################################
      include   'in_sys.com'
      character*(ln_length)head,tail
      integer first,last
      logical space
      first = 1
10    if (space(head(first:first)))then
		first = first +1
		if (first.le.ln_length) goto 10
      endif
      if (first.gt.ln_length) then
		head(1:ln_length) = empty_line
		tail(1:ln_length) = empty_line
		return
      endif
      last = first
20    if (.not.space(head(last:last)))then
		last = last +1
		if (last.le.ln_length) goto 20
      endif
      if (last.gt.ln_length) then
		head(1:ln_length) = head(first:ln_length)
		tail(1:ln_length) = empty_line
		return
      endif
      tail (1:ln_length)=  head(last:ln_length)
      head(1:ln_length) = head(first:last -1)
      return
      end

c######################################################################################
      block data insys
c######################################################################################
      include  'in_sys.com'
      logical first_time
      integer err_file
      common /err_sys/ first_time,err_file
      data first_time/.true./
      data nr_of_files /0/ 
      data check_range /.false./
      data error_occured /.false./
      end