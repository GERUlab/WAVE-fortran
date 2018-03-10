c	  version 3.0
c	  author: m. vanclooster 22/02/2002

      logical multimodal  
      logical ponded
      logical idvs, ilinr, ismxnl
      logical eq_grwt, igiveph
      logical pr_wat_histor,pr_wat_et, pr_wat_tseries, pr_wat_cum,pr_wat_sum
      integer histor_file

      double precision kc

      common/pfk_com/ sucreturn(kt_comps), wcreturn(kt_comps),
     $		parakon(kt_layrs,kt_parakon), parapf(kt_layrs,kt_parapf),
     $		hystlim, ihist(kt_comps), khys(kt_comps), modelk, modelpf, 
     $		nparkon, nparpf,multimodal
      common/wat_boco_com/ arel, brel, dxn, phsa, flxa, phsurf, 
     $		pond_max, pond, pond_from (366), pond_to (366), pond_maximum(366),
     $		gprh(366), gwl(366), qdeep(366), dstor, stor, rinf, ibboco, 
     $		ncomp, i_pond, ponded
      common/wat_upt_com/  akc, arer, brer, drz(366),dvs_kc(366),kc(366),
     $		epa, esa, epa_min_intc, ev_intc, p0, p1, p2h, p2l, p3, 
     $		rlai(366), rnam, rt_distr(kt_comps), rtex(kt_comps), tb, te, 
     $		transp_act, idvs, nr_of_kc_values, iday_kc(366),
     $		ilinr, ismxnl
      common/prhead_com/ dmcap(kt_comps), conduc(kt_comps), 
     $		conin(kt_comps+1), flxs(kt_comps+1), ph(kt_comps), ph1(kt_comps),
     $		wc(kt_comps), wc1(kt_comps), iter, itertot, nrepeat, 
     $		nr_no_conv
      common/wat_intgr_com/  cev_intc, cflbu, cflbup, 
     $		cflsd,  cflsdp, cinf, cirr, cprec, cpsev, cptra, 
     $		csev, crtex(kt_comps),  ctra, ctra1, cum_intc, 
     $		visd, voli, vol, 
     $		cwatupt_isd, cberr_wat, cberr_wat_isd,visdi
      common/wat_init_com/ eq_grwt, igiveph
      common/wat_out_com/ pr_wat_histor,pr_wat_et,pr_wat_tseries(5), 
     $		pr_wat_cum,  pr_wat_sum, histor_file
