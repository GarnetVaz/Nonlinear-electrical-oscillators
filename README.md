Nonlinear-electrical-oscillators
================================

Code to obtain the results in the paper.
The code has only been tested on a linux platform with Ubuntu 12.04

**Requirements:**

** Python 2.7.3 along with the following modules:
  1. Networkx - 1.7
  2. Scipy >= 0.12.0
  3. Numpy >= 1.6.1
  4. Matplotlib - 1.1.1rc

** R (used just for creating plots)
  1. reshape - 0.8.4   
  2. ggplot2 - 0.9.3.1 
  3. plyr_1.8

The code "linearcapcode.py" is the main code holding the class 
definition along with all methods.  To obtain a quick summary of a 
single simulation, use:
$ run linearcapcode.py

To obtain all the results described in the "Results: Comparison of 
Steady-State Algorithms" subsection of the paper, use the code 
"numerical_comparison.py".

The R codes are called internally from Python and need not be run 
independently.

To run simulations and create the plots described in the "Results: 
Gap Tuning" subsection of the paper, one also needs the Python module 
multiprocessing >= 0.70a1.

The code "graphmulti.py" will run the simulations corresponding to 
the barabasi-albert graphs.  Graph type can be changed as documented 
in the code to erdos-renyi or watts-strogatz graphs.  Plots will be 
generated automatically using R.

The code "graphmulti.py" is set to use 10 processors. This can be 
reduced if needed in line 231 with the variable 'numparallel'. 
Running the complete set of simulations for one type of graph will 
take 6+ hours.

