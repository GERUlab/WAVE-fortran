c	  version 3.0
c	  author m. vanclooster 22/02/2002

      integer max_nr_of_files, max_col_width
      logical free,prev_big, change_fm_small
      parameter (max_nr_of_files  = 30, max_col_width=20)
      parameter (change_fm_small=.true.)
      character *32 fname

      common/out_sys_com/ free(max_nr_of_files), 
     $		fname(max_nr_of_files),prev_big
