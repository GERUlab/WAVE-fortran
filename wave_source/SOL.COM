c	version 3.0
c	author: m. vanlooster 22/02/2002

      logical imobsw
      logical wet_and_dry
      logical pr_sol_sum, pr_sol_cum, pr_sol_tseries,out_sol

      common /sol_gen_com/ parasol(kt_parasol,kt_layrs,kt_sol),
     $		imobsw, nr_of_sol
      common /sol_wat_com/ flxsa (kt_comps+1), flxsah (kt_comps+1),
     $		pvela (kt_comps+1), pvelah (kt_comps+1), pveloh (kt_comps),
     $		wcio (kt_comps), wciob (kt_comps), wcma (kt_comps+1), 
     $		wcmah (kt_comps+1), wcmo (kt_comps), 
     $		wcmob (kt_comps), wco (kt_comps), wcob (kt_comps)
      common /soln_com/ 
     $		acsolib (kt_comps, kt_sol), acsolio(kt_comps, kt_sol),
     $		acsolmb(kt_comps, kt_sol), acsolmo(kt_comps, kt_sol),
     $		acsolo(kt_comps, kt_sol), csolb(kt_comps, kt_sol),
     $		csolib(kt_comps, kt_sol), csolio(kt_comps, kt_sol),
     $		csolmb(kt_comps, kt_sol), csolmo(kt_comps, kt_sol),
     $		csolo(kt_comps, kt_sol), diffus(kt_comps+1, kt_sol),
     $		sflxsa(kt_comps+1, kt_sol)
      common /sol_intgr_com/ cberr_sol(kt_sol), dleasa(kt_sol), 
     $		dsol(kt_sol), pleasa(kt_sol), rleasa(kt_sol),
     $		tcsolio(kt_comps, kt_sol), tcsolmo(kt_comps, kt_sol),
     $		tcsolo(kt_comps, kt_sol), tcsink(kt_sol), tflsol(kt_sol),
     $		tsoli(kt_sol), tsinki(kt_comps,kt_sol),
     $		tsinkm(kt_comps,kt_sol),
     $		rleasa_isd(kt_sol), dleasa_isd(kt_sol), pleasa_isd(kt_sol), 
     $		tcsink_isd(kt_sol), dsol_isd(kt_sol), cberr_sol_isd(kt_sol),
     $		tsoli_isd(kt_sol)
      common /sol_boco_com/ conirsol(kt_nirr,kt_sol), ddepsol(kt_sol),
     $		fsol(kt_ninof,kt_sol), wdepsol(kt_sol), ressol(kt_sol),
     $		solinfl(kt_sol), solsur(kt_sol), plo_in(kt_sol),
     $		idinof(kt_ninof),
     $		idirr(kt_nirr), idplow(kt_plo), ncsplow(kt_plo), ninof,
     $		nirr, nplo,wet_and_dry
      common / sol_sink_com/ decsoli(kt_comps,kt_sol),
     $		decsolm(kt_comps,kt_sol), rates(kt_comps,kt_sol),
     $		sinki(kt_comps,kt_sol), sinkm(kt_comps,kt_sol)
      common / sol_print_com/ pr_sol_sum, pr_sol_cum,
     $		pr_sol_tseries(3),out_sol(kt_sol)
