c	  version 3.0
c	  author m. vanclooster 22/02/2002

      logical pr_temp

      common/temp_com/wcto(kt_comps),wctob(kt_comps), 
     $		paratemp(kt_paratemp,kt_layrs), temp(kt_comps),
     $		pr_temp

