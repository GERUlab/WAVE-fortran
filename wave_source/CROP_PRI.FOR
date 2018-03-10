c	version 3.0
c	author: m. vanclooster 22/02/2002
c	corrections: m javaux 03/04/2002 
c###################################################################################
      subroutine crop_print
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'crop.com'

      if(pr_crop_gro)call crop_out_gro
      if(pr_crop_dev)call crop_out_dev
      end

c###################################################################################
      subroutine crop_out_gro
c     in   : dvs, ncrop, nday, rtdep, slai, slaig, tinit, wcob, wcrn, wkob, 
c            wlv, wrt, wso, wst
c     out  :
c     calls: date, heading, open_file
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'crop.com'
      integer crop_file
      integer year, month, day
      logical first_time
      character *(15) crop_name(5)
      data crop_name /'winter wheat','summer wheat','maize',
     $'potatoes','sugar beet'/
      data first_time /.true./
      
	if (first_time) then
		first_time = .false.
		call open_file(crop_file, 'crop_gro.out')
		call heading(crop_file, 'crop_gro.out',
     $		'state variables for crop yield')
		write (crop_file,10) crop_name(ncrop)
		if (ncrop.eq.3) then
			write (crop_file,20)
		else if (ncrop.eq.5) then
			write (crop_file,30)
		else if (ncrop.eq.4) then
			write (crop_file,45)
		else
			write (crop_file,40)
		endif
      endif
      call date_wave(nday +idint(tinit) - 1, year, month,day)
c     mais > cobs
      if (ncrop.eq.3) then
		write (crop_file,50) year, month, day, slaig,slai - slaig,
     $	dabs(rtdep), wlv, wst, wrt, wso, wkob,dvs
c	sugarbeets > crown
      else if (ncrop.eq.5) then
		write (crop_file,50)year, month, day, slaig,slai - slaig,
     $	dabs(rtdep), wlv, wst, wrt, wso, wcrn,dvs
c     all other crops
      else
		write (crop_file,60)year, month, day, slaig,slai - slaig,
     $	dabs(rtdep), wlv, wst, wrt, wso,dvs
      endif

10    format('crop : ',a)
20    format
     $    ('year month day    laig        laid     rtdep      wlv',
     $     '       wst       wrt       wso      wcob        dvs',/,
     $     '              (m**2/m**2) (m**2/m**2) (cm)      (kg/ha)',
     $     '   (kg/ha)   (kg/ha)   (kg/ha)   (kg/ha)',/)
30    format
     $    ('year month day    laig        laid     rtdep      wlv',
     $     '       wst       wrt       wso      wcrn       dvs',/,
     $     '              (m**2/m**2) (m**2/m**2)  (cm)     (kg/ha)',
     $     '   (kg/ha)   (kg/ha)   (kg/ha)   (kg/ha)     (°c)',/)
40    format
     $    ('year month day    laig        laid     rtdep      wlv',
     $     '       wst       wrt       wso       dvs',/,
     $     '              (m**2/m**2) (m**2/m**2)  (cm)     (kg/ha)',
     $     '   (kg/ha)   (kg/ha)   (kg/ha)',/)
45    format
     $    ('year month day    laig        laid     rtdep      wlv',
     $     '       wst       wrt       wso         dvs',/,
     $     '              (m**2/m**2) (m**2/m**2)  (cm)     (kg/ha)',
     $     '   (kg/ha)   (kg/ha)   (kg/ha)     (øc)',/)


50    format (i4,2x,i2,3x,i2,1x,f8.2,3x,f9.2,1x,f8.2,1x,f9.2,1x,f9.2,
     $  1x,f9.2,1x,f9.2,1x,f9.2,4x,f10.2)
60    format (i4,2x,i2,3x,i2,1x,f8.2,3x,f9.2,1x,f8.2,1x,f9.2,1x,f9.2,
     $  1x,f9.2,1x,f9.2,4x,f10.2)
      return
      end

c###################################################################################
      subroutine crop_out_dev
c     in   : dvs, ncrop, nday, rtdep, slai, slaig, tinit, wcob, wcrn, wkob, 
c            wlv, wrt, wso, wst
c     out  :
c     calls: date, heading, open_file
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'crop.com'
      integer crop_file
      integer year, month, day
      logical first_time
      character *(15) crop_name(5)
      data crop_name /'winter wheat','summer wheat','maize',
     $'potatoes','sugar beet'/
      data first_time /.true./

      if (first_time) then
		first_time = .false.
		call open_file(crop_file, 'crop_dev.out')
		call heading(crop_file, 'crop_dev.out',
     $	'state variables that characterise crop development')
		write (crop_file,10) crop_name(ncrop)
		if (simnit) then
			if (ncrop.lt.4) then
				write (crop_file,20)
			else
				write (crop_file,30)
			endif
		else
			if (ncrop.lt.4) then
				write (crop_file,40)
			else
				write (crop_file,50)
			endif
		endif
	endif 

      call date_wave(nday +idint(tinit) - 1, year, month,day)
      if (simnit) then
		write (crop_file,60) year, month, day,
     $		slaig,slai - slaig, dabs(rtdep),sucr_red_wat,
     $		sucr_red_nit,dvs
      else
		write (crop_file,60) year, month, day,
     $		slaig,slai - slaig, dabs(rtdep),sucr_red_wat,
     $		dvs
      endif

10    format('crop : ',a)
20    format
     $  ('year month day   laig         laid    rtdep   red_wat',
     $   ' red_nit     dvs',/,
     $   '              (m**2/m**2) (m**2/m**2) (cm)      (-)  ',
     $   '  (-)         ')
30    format
     $  ('year month day   laig         laid    rtdep   red_wat',
     $   ' red_nit     dvs',/,
     $   '              (m**2/m**2) (m**2/m**2) (cm)      (-)  ',
     $   '  (-)        (øc)')
40    format
     $  ('year month day   laig         laid    rtdep   red_wat',
     $   '    dvs',/,
     $   '              (m**2/m**2) (m**2/m**2) (cm)      (-)  ' )
50    format
     $  ('year month day   laig         laid    rtdep   red_wat',
     $   '    dvs',/,
     $   '              (m**2/m**2) (m**2/m**2) (cm)      (-)  ',
     $   '    (øc)')


60    format (i4,2x,i2,3x,i2,1x, f8.2,3x,f9.2,1x,f8.2,1x,
     $     f8.4,1x,f8.4,1x,f10.2)
      return
      end