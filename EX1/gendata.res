
***************************
gendata.res:
general input for the model
***************************



model switches
**************

plants     : -
temperature: -
solutes    : -
nitrogen   : -


profile development
*******************

compartment size            =   100.00 mm
number of soil layers       =     5
layer  1 has:   3 compartment(s)
layer  2 has:   2 compartment(s)
layer  3 has:   3 compartment(s)
layer  4 has:   1 compartment(s)
layer  5 has:  13 compartment(s)
bulk density (kg/l)
-------------------

layer    1
 1   0.0000
 2   0.0000
 3   0.0000
 4   0.0000
 5   0.0000


simulation time variables
*************************

start of calculations: 1990  1  1
end of calculations  : 1990  1 10


parameters concerning the numerical solution
********************************************

maximum time step                   :      1.00 days
minimum time step                   :    0.0010 days
maximum change of moisture content  :    0.0020 m**3/m**3
maximum balance error               :   0.01000 m**3/m**3/day


parameters concerning output
****************************

input testing
-------------

input is tested for range errors

output dates
------------

the time increment between printing of 
the summary table of the state variables is constant:   1 day(s)
the bottom compartment of the top layer, isd:   10

compartments for which there is output
--------------------------------------

  from   to
    1     3
    4     6
    7     9
   10    12
