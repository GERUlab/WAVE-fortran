c     Version 3.0
c	author:m. vanclooster (18/02/02)
c###################################################################################
      subroutine temp_print
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'temp.com'

      if (pr_temp)call temp_tseries
      end



c###################################################################################
      subroutine temp_tseries
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'temp.com'
      integer temp_file
      integer time
      dimension out(kt_comps)
      logical first_time
      data first_time /.true./ 

      time = nday +idint(tinit)-1
      if (first_time) then
		first_time = .false.
		call open_file(temp_file,'temp.out')
		call heading(temp_file,'temp.out',
     $		'average soil temperature profile(degrees c)')
		call heading_table(temp_file,from,to,nr_of_points,7)
      endif
      call avg (temp,out,1.0d0,1,from,to,nr_of_points)
      call wr_table (temp_file,out,time,nr_of_points, 7,1)
      return
      end