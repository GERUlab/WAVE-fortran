c	version 3.0
c	author m. vanclooster 15/02/2002

c######################################################################################
      subroutine wr_table (f_nr,x,time,nr_of_points,nr_of_digits, nr_of_dec)
c######################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'out_sys.com'
      integer	f_nr,nr_of_points,res,zero,time,l_t
      integer day, year, month
      parameter (l_t = 15)
      integer nodig,nodec,sp_col,space_in_column
      character *(3) blanks_l,blanks_r
      dimension x(*) 
      character *(60) basic_form
      character *(30) basic_desc
      character *(60) spec_form
      character *(30) spec_desc
      character *(3) ch_len_nodig, ch_len_nodec,ch_len_space
      character *(max_col_width*kt_comps)line

c     determine basic format 
      zero = ichar('0')
      column_size = max0(nr_of_digits +3,9)
      nodig= nr_of_digits +1
c     construct the format description
c     type conversion:idint-> char (ex: 7 -> '7')
      l_ndi = 0
10    if (nodig.ne.0) then
		l_ndi = l_ndi + 1
		res = mod(nodig,10)
		ch_len_nodig = char(res + zero) //ch_len_nodig
		nodig = nodig/10
		goto 10
      endif
      l_ndc = 0
      nodec= nr_of_dec
      if (nodec.eq.0) then
		l_ndc = l_ndc + 1
		ch_len_nodec = '0'
      else
20		res = mod(nodec,10)
		l_ndc = l_ndc + 1
		ch_len_nodec = char(res + zero) //ch_len_nodec
		nodec = nodec/10
		if (nodec.ne.0) goto 20
      endif
      basic_desc(1:1) = 'f'
      basic_desc(2:2+l_ndi-1) = ch_len_nodig
      basic_desc(2+l_ndi:2+l_ndi) = '.'
      basic_desc(3+l_ndi:3+l_ndi+l_ndc-1) = ch_len_nodec
c     determine number of blanks
c     calculate number of blanks on left and rigth of basic-form
      if  (nr_of_digits.ge.6) then
		nbl = 1
		nbr = 1
      else
		if (mod(nr_of_digits,2).eq.0) then    
			nbl = max0( min0(9,(8-nr_of_digits)/2),0)
			nbr = nbl
		else
			nbr = max0( min0(9,(7-nr_of_digits)/2),0)
			nbl = nbr +1
		endif
      endif
      zero = ichar('0')
      if (nbl.gt.0) then
		blanks_l = char(nbl+zero)//'x,'
      else
		blanks_l = ' '
      endif
      if (nbr.gt.0) then
		blanks_r = ','//char(nbr+zero)//'x'
      else
		blanks_r = ' '
      endif
      basic_form = '('//blanks_l//basic_desc(1:3+l_ndi+l_ndc-1)//
     $blanks_r//')'
      space_in_column = max0(nr_of_digits +1,7)
c     determine format for numbers equal or larger than 1.e(nr_of_digits - nr_ of_dec)
c     or numbers  smaller than  1.e(- nr_of_dec)
c     with as many decimals as possible (space_in_column - 6)
      l_spc = 0
      sp_col = space_in_column
30    if (sp_col.ne.0) then
		l_spc = l_spc + 1
		res = mod(sp_col,10)
		ch_len_space = char(res + zero) //ch_len_space
		sp_col = sp_col/10
		goto 30
      endif
      l_ndc = 0
      nodec= space_in_column - 6
40    res = mod(nodec,10)
      l_ndc = l_ndc + 1
      ch_len_nodec = char(res + zero) //ch_len_nodec
      nodec = nodec/10
      if (nodec.ne.0) goto 40
      spec_desc(1:1) = 'e'
      spec_desc(2:2+l_spc-1) = ch_len_space
      spec_desc(2+l_spc:2+l_spc) = '.'
      spec_desc(3+l_spc:3+l_spc+l_ndc-1) = ch_len_nodec
      spec_form = '(1x,'//spec_desc(1:3+l_spc+l_ndc-1)//',1x)'
      call date_wave (time, year, month,day) 
      write(line(1:l_t),50) year, month, day
50    format (1x,i4,2x,i2,3x,i2,1x)
      smallest_normal = 1.d0/(10.d0**(nr_of_dec))
      too_large = 10.d0**(nr_of_digits -nr_of_dec-1)
      do i = 1,nr_of_points
		if (
     $      ( dabs(x(i)).ge.smallest_normal.or..not.change_fm_small).and.
     $      (( x(i).lt.too_large.and.x(i).gt.(too_large/(-10.d0))).or.
     $	  x(i).eq.0.d0)) then
			write(line((i-1)*column_size+1+l_t:i*column_size+l_t),
     $		basic_form) x(i)
		else
			write(line((i-1)*column_size+1+l_t:i*column_size+l_t),
     $		spec_form) x(i)
		endif
	enddo 
      write(f_nr,60) line (1:column_size*nr_of_points+l_t)
60    format (a)
      return
      end

c######################################################################################
      subroutine avg  (x, out, factor,first,from,to,nr_of_points)
c######################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      dimension x(*), out(*)
      integer from, to,first
      dimension from(*), to(*)

      do i = 1, nr_of_points
		out(i)= 0.d0
		do  j = from(i), to(i)
			out(i) = x(j+first-1) +out(i)
		enddo
		out(i) = factor*out(i)/(to(i)-from(i)+1)
	enddo
      end

c######################################################################################
      subroutine sum  (x, out, factor,from,to,nr_of_points)
c######################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      dimension x(*), out(*)
      integer from, to
      dimension from(*), to(*)

      do i = 1, nr_of_points
		out(i)= 0.d0
		do j = from(i), to(i)
			out(i) = x(j) +out(i)
		enddo
		out(i) = factor*out(i)
	enddo
      end

c######################################################################################
      subroutine heading_table (f_nr,from,to,nr_of_points,nr_of_digits)
c######################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'out_sys.com'
      integer f_nr, from,to,column_size,zero,nr_of_digits
      dimension from(*),to(*)
      character *(8) basic_form1
      character *(10) basic_form2
      character *(3) blanks_l1,blanks_l2
      character *(22) form_desc1,form_desc2
      character *(max_col_width*kt_comps)heading
      integer l_t
      parameter (l_t = 15)
      data basic_form2 /'i3,'' /'',i3'/
      data basic_form1 /'2x,i3,2x'/

      column_size = max0(nr_of_digits + 3,9)
c     calculate number of blanks on left of basic-form1 (7 wide), nbl1
      nbl1 = column_size-8
      nbl2 = column_size-9
      zero = ichar('0')
      if (nbl2.eq.0) then
		form_desc1 = '(1x,'// basic_form1//',1x)'
		form_desc2 = '('//basic_form2//',1x)'
      else
		blanks_l1(1:1) = char(nbl1+zero)
		blanks_l1(2:3) = 'x,'
		blanks_l2(1:1) = char(nbl2+zero)
		blanks_l2(2:3) = 'x,'
		form_desc1 = '('//blanks_l1// basic_form1//',1x)'
		form_desc2 = '('//blanks_l2//basic_form2//',1x)'
      endif
      write (heading (1:l_t),10)'year month day'
10    format (a14,1x)
      do i = 1,nr_of_points
		if (from(i).eq.to(i)) then
			write(heading((i-1)*column_size+1+l_t:i*column_size+l_t),
     $		form_desc1) from(i)
		else
			write(heading((i-1)*column_size+1+l_t:i*column_size+l_t),
     $		form_desc2) from(i),to(i)
		endif
	enddo 
      write (f_nr,20) 'compartments'
20	format (6x,a20)
      write(f_nr,30)  heading (1:column_size*nr_of_points+l_t)
30    format (a)
      end

c######################################################################################
      subroutine heading(file,filename,what)
c######################################################################################
      include  'constant'
      include  'out_sys.com'
      integer file
      character*(*) filename
      character*(*) what
      character*(80) bar
      character bar_cont(80)
      equivalence (bar,bar_cont)
      integer ln
      data bar_cont/80*'*'/

      ln = max0(len(what ),len(filename))
      write (file,10) bar(1:ln),filename, what,bar(1:ln)
10	format(/,a,/, a ,':',/, a,/, a,/)
      return
      end 

c######################################################################################
      subroutine big_title(file,title)
c######################################################################################
      include  'constant'
      include  'out_sys.com'
      integer file
      character*(*) title
      character*(80) bar
      character bar_cont(80)
      equivalence (bar,bar_cont)
      integer ln
      data bar_cont/80*'*'/

      ln = len(title )
      write (file,10) title, bar(1:ln)
      prev_big = .true.
10    format(//, a,/, a,/)
      return
      end

c######################################################################################
      subroutine small_title(file,title)
c######################################################################################
      include  'constant'
      include  'out_sys.com'
      integer file
      character*(*) title
      character*(80) bar
      character bar_cont(80)
      equivalence (bar,bar_cont)
      integer ln
      data bar_cont/80*'-'/

      ln = len(title )
      if (prev_big) then
		prev_big = .false.
		write (file,10) title, bar(1:ln)
      else
		write (file,20) title, bar(1:ln)
      endif
10    format( a,/, a,/)
20    format( /,a,/, a,/)
      return
      end

c######################################################################################
      subroutine date_title(file,time)
c######################################################################################
      implicit double precision (a-h, o-z)
      include  'constant'
      include  'out_sys.com'
      include  'gen.com'
      integer file, day, month, year, time

      call date_wave (time, year, month, day)
      write (file,10) year, month, day
10    format (//,72('*'),/
     $           ,'date (y/m/d): ',  i4,1x,i2,1x,i2, /,
     $           '-------------')
      return
      end


c######################################################################################
      subroutine open_file (file,file_name)
c######################################################################################
      include  'constant'
      include  'out_sys.com'
      integer file
      character*(*) file_name
      character*70  message

      call get_unit(file) 
      fname(file) = file_name      
      open(unit = file, file = file_name, err= 10)
      return
10    write(message,20) file_name 
      call report_err(message)
20    format (' unable to open output file: ', a32)        
      call stop_simulation ('programme stopped: error during output')
      end


c######################################################################################
      subroutine append_file (file,file_name)
c######################################################################################
      include  'constant'
      include  'out_sys.com'
      integer file
      character*(*) file_name
      character*70  message

      call get_unit(file) 
      fname(file) = file_name      
      open(unit = file, file = file_name, err= 10, access='append')
      return
10    write(message,20) file_name 
      call report_err(message)
20    format (' unable to open output file: ', a32)        
      call stop_simulation('programme stopped: error during output')
      end


c######################################################################################
      subroutine close_file (file)
c######################################################################################
      include  'constant'
      include  'out_sys.com'
      integer file
      character*70  message

      close(unit = file, err= 10)
      call free_unit(file) 
      return
10    write(message,20) fname (file)
      call report_err(message)
20    format (' unable to close output file: ', a32)        
      call stop_simulation('programme stopped: error during output')
      end


c######################################################################################
      subroutine get_unit (file)
c######################################################################################
      include  'constant'
      include  'out_sys.com'
      integer file

      i= 1
10    if (free(i).and.i.ne.5.and.i.ne.6) then
		file = i 
      else if (i .lt.  max_nr_of_files) then
		i = i + 1
		goto 10
      else
		call report_err (' too many files')
		call stop_simulation ('programme stopped: too many files')
      endif
      free(i) = .false.
      end


c######################################################################################
      subroutine free_unit (file)
c######################################################################################
      include  'constant'
      include  'out_sys.com'
      integer file

      free(file) = .true.
      end


c######################################################################################
      subroutine stop_simulation (message)
c######################################################################################
      include  'constant'
      include  'out_sys.com'
      character*(*)  message

c     close all output
      i = 1
10    if (.not.free(i)) close(i, err= 20)
20    if (i .lt.  max_nr_of_files) then
		i = i + 1
		goto 10
      endif
      write (*,*) message
      stop
      end


c######################################################################################
      block data outsys
c######################################################################################
      include  'constant'
      include  'out_sys.com'
      data free /max_nr_of_files*.true./
      end