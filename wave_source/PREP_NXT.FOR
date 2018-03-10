      subroutine prep_nxt_yr
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'


      call prep_nxt_yr_wat
      if (simsol) call prep_nxt_yr_sol
      if (simnit) call prep_nxt_yr_nit
      end 


      subroutine prep_nxt_yr_wat
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      integer file


      call open_file(file, 'watinit.wi')
      write (file,10)
      do i = 1, ncs
      write (file,20) i, ph(i)/10.d0
      enddo
      write (file,30)
10    format ('n',/,'y',/)
20    format (i3,1x,f24.8)
30    format ('et')
      call close_file (file)
      end

 
      subroutine prep_nxt_yr_sol
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'sol.com'
      integer file,sp


      call open_file(file, 'solinit.si')
      do i = 1, ncs
      if (imobsw)then
        write (file,10) i, 
     >     (tcsolmo(i,sp)+tcsolio(i,sp), sp = 1, nr_of_sol)
      else
        write (file,10) i, (tcsolo(i,sp), sp = 1, nr_of_sol)
      endif
      enddo
      write (file,20)
10    format (i3, 5(1x,f15.2))
20    format ('et')
      call close_file (file)
      end


 
      subroutine prep_nxt_yr_nit
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'nit.com'
      integer file
      call open_file(file, 'nitinit.ni')
      do i = 1, ncs
      write (file,10) i, dmax1(ccmano(i),0.00001d0),
     >   dmax1(cnmano(i),0.00001d0),
     >   dmax1(cclito(i),0.00001d0),
     >   dmax1(cnlito(i),0.00001d0),
     >   dmax1(cchumo(i),0.00001d0),
     >   dmax1(cnhumo(i),0.00001d0)
      enddo
      write (file,20)
10    format (i3, 6(1x,f16.5))
20    format ('et')
      call close_file (file)
      end
