c	  version 3.0
c	  author: m.vanclooster 22/02/2002

      logical pr_crop_gro, pr_crop_dev 

      common/ sucr_com/amx, asrqso, dvs, eai, ear, eff,rgr,
     $					rkdf, rlaicr, rlat, rldf(kt_comps),
     $					rmainso, rmatr, rmrd,
     $					rtdep, scp, plai,sla, slai, slaig,
     $					specweig,sucr_red_wat, sucr_red_nit,
     $					ssl, table(2,kt_cropinfo,kt_table),
     $					tsum, tbase, wcrn, wkob, wlv, wlvg, wrt, 
     $					wso, wst, ncrop, nr_of_tables, nsl, ntabel(kt_table),
     $					pr_crop_gro, pr_crop_dev