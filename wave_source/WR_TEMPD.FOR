c	version 3.0
c	author: m vanclooster 22/02/2002
c###################################################################################
      subroutine wr_tempdata
c     in   : npl, paratemp, pr_temp
c     out  : -
c     calls: big_title, close_file, heading, open_file
c###################################################################################
      implicit double precision (a-h, o-z)
      include   'constant'
      include   'gen.com'
      include   'temp.com'
      integer res_file

      call open_file (res_file, 'tempdata.res')
      call heading (res_file, 'tempdata.res',
     $ 'input for mopeling heat transport')
      call big_title (res_file,'parameters for each soil layer')
      call big_title (res_file,'output')
      if ( pr_temp) then
		write (res_file,10)plus
      else
		write (res_file,10)neg
      endif
10    format ('the temperature profile (temp.out): ',a)
      call close_file(res_file)
      end