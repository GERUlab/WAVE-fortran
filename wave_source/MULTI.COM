c     version 3.0
c	  author: m. vanclooster 22/02/2002

      parameter (nr_por_class = 5)

      common/n_poros/ alpha(nr_por_class,kt_layrs),
     $	rn(nr_por_class,kt_layrs), rm(nr_por_class,kt_layrs),
     $	p(nr_por_class,kt_layrs), q(nr_por_class,kt_layrs),
     $	w(nr_por_class,kt_layrs), wcr(kt_layrs), wcs(kt_layrs),
     $	xexp(kt_layrs), csat(kt_layrs),total(kt_layrs),
     $	nr_poros(kt_layrs)

