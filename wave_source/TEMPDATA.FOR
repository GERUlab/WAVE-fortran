c	version 3.0
c	m. vanclooster 18/02/02
c###################################################################################
      subroutine tempdata
c     in   : npl, simnit, simsol
c     out  : pr_temp
c     calls: bulk_density, check_et, end_file, hold_line, init_file, next_line,
c            rd_int, rd_logic, rd_real, report_err, temp_pr_for_nit
c  correction 10.3.94: reading value of pr_temp
c###################################################################################
      implicit double precision (a-h, o-z)
      include  'constant'
      include  'gen.com'
      include  'temp.com'
      logical rd_logic

      call init_file ('tempdata.in')
      pr_temp = rd_logic()
      call end_file
      end
	