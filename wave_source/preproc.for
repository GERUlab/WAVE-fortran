c	version 3.0
c	author: m. vanclooster 21/02/2002

c###################################################################################
	subroutine program_version
c###################################################################################
c	prompt interface
	write (*,*) 'wave version 3.0'
      write (*,*) 'oct 2004 _ Version Marijn'
	return 
	end 


c###################################################################################
	subroutine reading
	implicit double precision (a-h,o-z)
      include 'constant'
      include 'gen.com'

c	read input
      call gendata
      call wr_gendata
      call watdata
      call wr_watdata
      if(simsol) then
		call soldata
		call wr_soldata
		if(simnit) then
			call nitdata
			call wr_nitdata
		endif
      endif
      if (isucr) then
		call cropdata
		call wr_cropdata
      endif
      if (simtemp) then
		call tempdata
		call wr_tempdata
      endif
      call climdata
      call wr_climdata
      if (input_error()) then
		call stop_simulation ('programme stopped: error during input')
      endif
	return
	end 
