# solpos.cs
c# port of solpos.c

The original is available here https://www.nrel.gov/grid/solar-resource/solpos.html The original structure is still intact. Only the minimum was done to make it compile. It’s still very C like. But its working, and Test function is also included. Precision is much higher with C# as you will see in the test function. More could be done to cSharpify but for now I’m happy to have the original structure remain intact. 

Original Disclaimer

By accessing these files, you agree to abide by the NREL data disclaimer. This software has been tested on a variety of platforms, but it is not guaranteed to work on yours. It is provided here as a convenience 
https://www.nrel.gov/disclaimer.html


Test Output:
```
***** TEST S_solpos: *****
Note that your final decimal place values may vary
based on your computer's floating-point storage and your
compiler's mathematical algorithms.  If you agree with
NREL's values for at least 5 significant digits, assume it works.

Note that S_solpos has returned the day and month for the
input daynum.  When configured to do so, S_solpos will reverse
this input/output relationship, accepting month and day as
input and returning the day-of-year in the daynum variable.

NREL    -> 1999.07.22, daynum 203, retval 0, amass 1.335752,         ampress 1.326522
SOLTEST -> 1999.07.22, daynum 203, retval 0, amass 1.33575576483888, ampress 1.3265254683395
NREL    -> azim 97.032875,        cosinc 0.912569,          elevref 48.409931
SOLTEST -> azim 97.0333143862794, cosinc 0.912569749169862, elevref 48.4097498637026
NREL    -> etr 989.668518,       etrn 1323.239868,      etrtilt 1207.547363
SOLTEST -> etr 989.665707776715, etrn 1323.23983749449, etrtilt 1207.54864659392
NREL    -> prime 1.037040,         sbcf 1.201910,         sunrise 347.173431
SOLTEST -> prime 1.03704002282382, sbcf 1.20191088527296, sunrise 347.174605315733
NREL    -> sunset 1181.111206,      unprime 0.964283,          zenref 41.590069
SOLTEST -> sunset 1181.11019305054, unprime 0.964282937969009, zenref 41.5902501362974

Raw airmass loop:
NREL    -> 37.92  5.59  2.90  1.99  1.55  1.30  1.15  1.06  1.02  1.00
SOLTEST -> 37.92  5.59  2.90  1.99  1.55  1.30  1.15  1.06  1.02  1.00 
```
