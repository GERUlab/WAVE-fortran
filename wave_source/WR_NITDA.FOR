c	version 3.0
c	author: m.vanclooster 26/02/2002
c###################################################################################
      subroutine wr_nitdata
c     in   : alfa_rdens, carborg, cchumo, cclito, ccmano, cnhumo, cnlito,
c            cnmano, fag, fe, fh, fhp, flr, g, idman, isucr, ncs, nof,
c            pr_nit_cum, pr_nit_sum, pr_nit_tseries, pr_temp_for_nit, rdo,
c            rkdenit, rkhum, rkhyd, rklit, rkman, rknitri, rkvol, rnitorg,
c            rnmaxp, ro, rorad, simplant, w0_rdens
c     out  : -
c     calls: big_title, close_file, date, heading, open_file, small_title
c###################################################################################
      implicit double precision (a-h, o-z)
      include   'constant'
      include   'gen.com'
      include   'nit.com'
      integer appl
      integer res_file
      integer day, month, year

      call open_file (res_file, 'nitdata.res')
      call heading (res_file, 'nitdata.res','input for nitrogen modeling')
c     nitrate : upper boundary conditions
      call big_title(res_file,'upper boundary conditions')
      if (simplant) then
		call small_title(res_file,'plant contribution to the organic pool')
		write (res_file,3300) fag
		write (res_file,3302) flr
		write (res_file,3304) fhp
      endif
c     organic amendments 
      call small_title(res_file,'manure apllication')
      write (res_file,3062) nof
      if (nof.gt.1) write (res_file,3072)
      do  appl = 1,nof
		call date_wave (idman(appl),year, month,day) 
		write (res_file,3082)year, month, day, carborg(appl),rnitorg(appl)
      enddo
c	sink terms
      call big_title(res_file,'sink term')
c	organic matter tun over
      call small_title(res_file,'organic matter transformation')
      write (res_file,3200) ro
      write (res_file,3202) fe
      write (res_file,3204) fh
c     nitrogen uptake
      if (simplant) then
		call small_title(res_file,'nitrogen uptake')
		write (res_file,3210) rorad
		write (res_file,3212) rdo
		if (.not.isucr) then
			write (res_file,3206) g
			write (res_file,3208) rnmaxp
			write (res_file,3250) w0_rdens
		endif
      endif
c     n transformation rates
      call big_title(res_file,'transformation rates for organic matter')
      write (res_file,3150)
      do i = 1, ncs
		write (res_file,3160)i, rknitri(i), rkdenit(i), rkhyd(i), 
     $	rkvol(i), rklit(i), rkman(i), rkhum(i)
      enddo
      call big_title(res_file,'initial values')
      write (res_file,3190)
      do i = 1, ncs
        write (res_file,3195)i, ccmano(i), cnmano(i), cclito(i),
     $  cnlito(i),cchumo(i), cnhumo(i)
      enddo
c     print control
      call big_title(res_file, 'output')
      if (pr_nit_sum) then
		write (res_file,4000)plus
      else
		write (res_file,4000)neg
      endif
      if (pr_nit_cum) then
		write (res_file,4020)plus
      else
		write (res_file,4020)neg
      endif
      if (pr_nit_tseries(1))then
		write (res_file,4040)plus
      else
		write (res_file,4040)neg
      endif
      if (pr_nit_tseries(2))then
		write (res_file,4060)plus
      else
		write (res_file,4060)neg
      endif
      if (pr_nit_tseries(3))then
		write (res_file,4080)plus
      else
		write (res_file,4080)neg
      endif
      if (pr_nit_tseries(4))then
		write (res_file,5000)plus
      else
		write (res_file,5000)neg
      endif
      if (pr_nit_tseries(5))then
		write (res_file,5020)plus
      else
		write (res_file,5020)neg
      endif
      if (pr_nit_tseries(6))then
		write (res_file,5040)plus
      else
		write (res_file,5040)neg
      endif
      if (pr_nit_tseries(7))then
		write (res_file,5060)plus
      else
		write (res_file,5060)neg
      endif
      if (pr_nit_tseries(8))then
		write (res_file,5080)plus
      else
		write (res_file,5080)neg
      endif
      if (pr_nit_tseries(9))then
		write (res_file,6000)plus
      else
		write (res_file,6000)neg
      endif
      if (pr_nit_tseries(10))then
		write (res_file,6020)plus
      else
		write (res_file,6020)neg
      endif
      if (pr_nit_tseries(11))then
		write (res_file,6040)plus
      else
		write (res_file,6040)neg
      endif
      if (pr_nit_tseries(12))then
		write (res_file,6060)plus
      else
		write (res_file,6060)neg
      endif
      if (pr_nit_tseries(13))then
		write (res_file,6080)plus
      else
		write (res_file,6080)neg
      endif
      if (pr_nit_tseries(14))then
		write (res_file,7000)plus
      else
		write (res_file,7000)neg
      endif
      if (pr_nit_tseries(15))then
		write (res_file,7020)plus
      else
		write (res_file,7020)neg
      endif
      if (pr_nit_tseries(16))then
		write (res_file,7040)plus
      else
		write (res_file,7040)neg
      endif

      call close_file(res_file)

3300  format('fraction n in crop residue at harvest: ', f10.2)
3302  format('fraction n in living roots at harvest: ', f10.2)
3304  format('fraction n removed at harvest        : ', f10.2)
3062    format('the number of applications of manure : ',
     $   1x,i2)
3072    format ('year  mo day    carborg    rnitorg ',/,
     $          '               (mg/m**2)  (mg/m**2)')
3082  format (i4,2x,i2,2x,i2,1x,2(f10.1,1x))

3200  format('ro       : ', f10.2)
3202  format('fe       : ', f10.2)
3204  format('fh       : ', f10.2)
3206  format
     $('fraction of growing season before maximum uptake : ', f10.2)
3208  format
     $('maximum amount of n accumulated at harvest       : ', f10.2)
3210  format
     $('root radius                                      : ', f10.2,
     $ ' (mm)')
3212  format
     $('distance between root surface and free solution  : ', f10.2,
     $ ' (mm)')
3250  format 
     $('root densty in the first compartment             : ',f10.2)
3150  format (
     $'comp rknitri  rkdenit    rkhyd    rkvol   ',
     $' rklit    rkman    rkhum',/,
     $'     (day-1)  (day-1)  (day-1)  (day-1)   ',
     $'(day-1) (day-1)  (day-1)')
3160    format (i3,1x,8(e8.1,1x))
3190  format (
     $'comp  ccmano   cnmano   cclito   cnlito   cchum',
     $'o   cnhumo',/,
     $'     (mg/m**2)(mg/m**2)(mg/m**2)(mg/m**2)(mg/m**2',
     $')(mg/m**2) ')
3195    format (i3,1x,6(e8.1,1x))
4000  format 
     $('overview of main state variables (nit_sum.out)    : ',a)
4020  format  
     $('cumulative terms of the organic pool (nit_cum.out): ',a)
4040  format  
     $('organic litter n concentration (orgnlit.out)      : ',a)
4060  format  
     $('organic manure n concentration (orgnman.out)      : ',a)
4080  format  
     $('organic humus n  concentration (orgnhum.out)      : ',a)
5000  format  
     $('organic litter c concentration (orgclit.out)      : ',a)
5020  format  
     $('organic manure c concentration (orgcman.out)      : ',a)
5040  format  
     $('organic humus c concentration (orgchum.out)       : ',a)
5060  format  
     $('cumulative nh4 uptake (nh4upt.out)                : ',a)
5080  format  
     $('cumulative no3 uptake (no3upt.out)                : ',a)
6000  format  
     $('cumulative denitrification (denitrif.out)         : ',a)
6020  format  
     $('cumulative nitrification (nitrific.out)           : ',a)
6040  format  
     $('cumulative hydrolysis (hydrolys.out)              : ',a)
6060  format  
     $('cumulative volatilisation (volatil.out)           : ',a)
6080  format  
     $('cumulative nh4 mineralisation (nh4min.out)        : ',a)
7000  format  
     $('cumulative no3 mineralisation (no3min.out)        : ',a)
7020  format
     $('cumulative no3 convective uptake (no3_cupt.out)   : ',a)
7040  format
     $('cumulative no3 diffusive uptake (no3_dupt.out)    : ',a)

      end