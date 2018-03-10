c	  version 3.0
c	  author: m. vanclooster 22/02/2002

      integer yearin,from,to,tprint,plant_date,harvest_date
      logical isucr,simnit, simsol,simtemp,simplant
      logical iprct,printing_time,range_err ,mo_1_year

      common  

c profile dimensions 
     $dsp, dx, dxinter(kt_comps +1), x(kt_comps), bulk_dens(kt_layrs),

c time  variables
     $dt, t, tend, tinit,  tmd, 

c numerical aspects water flow
     $devstop, dthm, dtmax, dtmin, 

c profile dimensions
     $ncs, npl, nla(kt_comps+1), ncsl(kt_layrs),
 
c time variables
     $yearin, ndsim, nday, ndts, nts,
     $plant_date,harvest_date,mo_1_year,
     $init_month,init_year,init_day,

c printing control variables
     $from(kt_comps),to(kt_comps),tprint,
     $nr_of_points, npr, isd,
     $printing_time(366), iprct ,range_err,

c main switches
     $isucr, simnit, simsol,simtemp,simplant

