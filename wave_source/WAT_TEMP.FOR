c	version 3.0
c	author:m vanclooster 21/02/02
c###################################################################################
      subroutine calc_wat_temp
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'wat.com'
      include   'temp.com'

      do i=1,ncs
		wcto(i)=wc(i)
		wctob(i) = wc1(i)
      enddo
      end
