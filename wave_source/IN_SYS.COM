c	  version 3.0
c	  author: m.vanclooster 22/02/2002

      integer nr_of_lines
      integer input_unit 
      integer ln_length , ms_length

      logical file_end, hold, end_of_table ,check_range
     $	,error_occured, report_err_for_line, ok1

      parameter (ln_length = 200 , ms_length = 80 ) 

      character *32 fname
      character *(ln_length) line
      character *(ln_length) empty_line

      parameter (empty_line = ' ',ok1 = .true.)

      common/bas_ios/ file_end,hold,end_of_table , 
     $		nr_of_lines(10),input_unit(10), fname(10),line, nr_of_files,
     $		check_range,error_occured, report_err_for_line
