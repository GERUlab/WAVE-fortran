
**********************************
watdata.res:
input for modeling water transport
**********************************



parameters for each soil layer
******************************

model type for the mrc                :    1
model type for the conductivity curve :    1
parameters for mrc model
------------------------

 layer    1          2          3          4          5
   1  0.471E-02  0.492E+00  0.182E-02  0.333E+00  0.100E+01
   2  0.886E-01  0.469E+00  0.216E-02  0.499E+00  0.100E+01
   3  0.859E-01  0.377E+00  0.108E-02  0.908E+00  0.100E+01
   4  0.859E-01  0.377E+00  0.108E-02  0.908E+00  0.100E+01
   5  0.718E-01  0.397E+00  0.420E-03  0.571E+00  0.100E+01

parameters for conductivity model
---------------------------------

 layer    1          2          3
   1  0.220E+03  0.143E+01  0.149E+01
   2  0.248E+03  0.177E+01  0.145E+01
   3  0.169E+03  0.137E+01  0.150E+01
   4  0.140E+03  0.121E+01  0.154E+01
   5  0.955E+02  0.133E+01  0.151E+01


upper boundary condition
************************


the minimum allowed pressure head at the soil surface    : -0.10E+08 (cm)
the maximum ponding depth   0.00 (mm)


sink (evaporation and plant water uptake)
*****************************************

kc values
---------

kc as a function of time
1990  1  1   1.000


bottom boundary condition
*************************

the groundwater level: 
year mo day   level
               (mm)
1990  1  1    -2000.0
1990  1  2    -2000.0
1990  1  3    -2000.0
1990  1  4    -2000.0
1990  1  5    -2000.0
1990  1  6    -2000.0
1990  1  7    -2000.0
1990  1  8    -2000.0
1990  1  9    -2000.0
1990  1 10    -2000.0


initial values
**************

the pressure head profile is calculated in equilibrium with the groundwater level
initial wc / ph profile
-----------------------

     depth         ph      wc
       (mm)       (cm) (m**3/m**3)
      -50.0      -195.0   0.290
     -150.0      -185.0   0.292
     -250.0      -175.0   0.294
     -350.0      -165.0   0.327
     -450.0      -155.0   0.329
     -550.0      -145.0   0.331
     -650.0      -135.0   0.334
     -750.0      -125.0   0.336
     -850.0      -115.0   0.339
     -950.0      -105.0   0.350
    -1050.0       -95.0   0.352
    -1150.0       -85.0   0.355
    -1250.0       -75.0   0.357
    -1350.0       -65.0   0.360
    -1450.0       -55.0   0.363
    -1550.0       -45.0   0.366
    -1650.0       -35.0   0.370
    -1750.0       -25.0   0.375
    -1850.0       -15.0   0.380
    -1950.0        -5.0   0.388
    -2050.0         0.0   0.397
    -2150.0         0.0   0.397


output
******

the iteration history (wat_hist.out)                 : -
overview of main state variables (wat_sum.out)       : +
cumulative terms of the water balance (wat_cum.out)  : +
evapotransp. and gr. water table level               : +
the pressure head (ph.out)                           : +
the water content (wc.out)                           : +
the cumulative root extraction (crtex.out)           : -
