Nonlinear-electrical-oscillators
================================

Code to obtain the results in the paper.
The code has only been tested on a linux platform with Ubuntu 12.04

*Requirements:

** Python 2.7.3 along with the following modules:
  1. Networkx - 1.7
  2. Scipy >= 0.12.0
  3. Numpy >= 1.6.1
  4. Matplotlib - 1.1.1rc

** R (used just for creating plots)
  1. reshape - 0.8.4   
  2. ggplot2 - 0.9.3.1 
  3. plyr_1.8

The code _linearcapcode.py_ is the main code holding the class definition along with all methods.
To obtain a quick summary of a single simulation you can use:
$ run linearcapcode.py

To obtain all the results described in Section 5 of the paper use the code *numerical_comparison.py*
