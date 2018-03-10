c    Version 3.0
c    Author: m vanclooster 22/02/2002

      common /clim_com/ et0(366), prec(366), rirr(366),
     $					fin(366),tmin(366),tmax(366),hsh(366)
