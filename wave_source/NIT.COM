c	Version 3.0
c	author: m. vanclooster 26/02/2002
      common /nit_upt_com/  g, rdo, rfract_tot,g1a,
     $	rnmaxp, rorad, rdens_left(kt_comps), 
     $	suni(kt_comps,kt_sol),sunm(kt_comps,kt_sol),w0_rdens,
     $	unc(kt_comps,kt_sol),und(kt_comps,kt_sol),
     $	fr_n_st_le, fr_n_rt_st, fr_n_mn_le
      common /nit_boco_com/carblit, carbman, carborg(kt_nof), fag,
     $	fhp, flr, tflnorg, tflcorg, rnitorg(kt_nof), idman(kt_nof),
     $	nof
      common /nit_gen_com/cchumo(kt_comps), cclito(kt_comps), 
     $	ccmano(kt_comps),cnhumo(kt_comps), cnlito(kt_comps),
     $	cnmano(kt_comps), ro
      common /nit_sink_com/deccorg(kt_comps), decnorg(kt_comps),
     $	fe, fh, rcarbmin(kt_comps), rdeniti(kt_comps), rdenitm(kt_comps), 
     $	rhureai(kt_comps), rhuream(kt_comps), rkdenit(kt_comps), 
     $	rkhum(kt_comps), rkhyd(kt_comps), rklit(kt_comps), 
     $	rkman(kt_comps), rknitri(kt_comps), rkvol(kt_comps), 
     $	rmini(kt_comps,kt_sol), rminm(kt_comps,kt_sol), 
     $	rnitrii(kt_comps), rnitrim(kt_comps), rvoli(kt_comps), 
     $	rvolm(kt_comps)
      common /nit_intgr_com/ cden(kt_comps),chyd(kt_comps),
     $	cmin(kt_comps,kt_sol),cnit(kt_comps),cvol(kt_comps),
     $	cupt(kt_comps,kt_sol), ptden, pthyd, ptmin(kt_sol), ptnit,
     $	ptup(kt_sol), ptvol, tcorgs, tnorgs,ptclit,ptcman,ptchum,
     $	ptcarborg, ptnlit,ptnman,ptnhum,ptnitorg,dcorg,dnorg,
     $	cberr_corg,cberr_norg,ptsnorg,ptscorg,
     $	ptden_isd, pthyd_isd, ptmin_isd(kt_sol), ptnit_isd,
     $	ptup_isd(kt_sol), ptvol_isd, tcorgs_isd, tnorgs_isd,
     $	ptclit_isd,ptcman_isd,ptchum_isd,
     $	ptcarborg_isd, ptnlit_isd,ptnman_isd,ptnhum_isd,
     $	ptnitorg_isd,dcorg_isd,dnorg_isd,
     $	cberr_corg_isd,cberr_norg_isd,ptsnorg_isd,ptscorg_isd,
     $	cunc_no3(kt_comps),cund_no3(kt_comps)
      logical pr_nit_sum, pr_nit_cum, pr_nit_tseries
      common /nit_print_com/ pr_nit_sum, pr_nit_cum, 
     $pr_nit_tseries(16)