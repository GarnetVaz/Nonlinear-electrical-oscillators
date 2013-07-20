#!/usr/bin/python
"""Code to run the numerical solver for the grid graph 20x20.

This code will compute and solve all that is mentioned in the
numerical comparison section of the paper (Section 5).

Specifically :
1) Fixed point error comparison for all three methods. (fixedpoint.txt)
2) Energy conservation for all methods. (Table)
3) Comparison with numerical integrator. (Table)
4) Decay of the fourier coefficients. (fdecay.txt)
"""

import linearcapcode as cc
import numpy as np
import sys

def energy_comparison():
  # Table 4
  print "Energy balance for numerical solution    : {0:g}".format(data.NEDiff.real)
  print "Energy balance for perturbative solution : {0:g}".format(data.pertEDiff.real)
  print "Energy balance for iterative solution    : {0:g}".format(data.iterEDiff.real)

def fixedpointerror():
  # Table 3
  numerr = data.fixed_point_error(data.Nfk)
  perterr = data.fixed_point_error(data.FSol.T)
  itererr = data.fixed_point_error(data.alphamat)
  print "Fixed point error for numerical method    : {0:g}".format(numerr)
  print "Fixed point error for perturbative method : {0:g}".format(perterr)
  print "Fixed point error for iterative method    : {0:g}".format(itererr)

# Create any random graph. Graph structure will be overwritten with a grid later.
latsize = 20
gtype = 'erdos_renyi_graph'
param = dict(n = latsize**2, p = 0.7)

# Builds graph.
data = cc.Mygraph(gtype, param)

# Create all parameters for input.
Cap0 = np.ones(data.N)
G = np.ones(data.N) * 0.01
G[-latsize] = 1.0

Input_nodes = np.zeros(data.N)
Input_nodes[:latsize] = 1
Input_nodes[np.arange(latsize, size + 1, latsize) - 1] = 1 # +1 is to ensure that the bot-left node has two inputs.
Eps = 0.5
L = 1 * (np.ones(2*latsize*(latsize - 1) + sum(Input_nodes)))

# The forcing amplitude here is the numerator. 0.03 * Sin(Omega t)
forc_amp = 0.03/2
forc_coef = np.array((0, -1))
Ord_req = 21

# Initialize the graph.
tpoints = 64
data.init_graph(Cap0, Eps, L, G, Input_nodes, forc_amp, forc_coef, Ord_req, tpoints, Omega = -1.0, gtype = 'grid')

# Set the value of Omega.
data.Omega = 1.0
data.solve()
data.fouriersol()
# Solve the system using the numerical integrator.
data.numeric_solver(doplot = False, _cycles = 1500)

# Obtain a reference solution from the iterative method.
try :
  data.iterative_solver(_tol = 3.0e-16)
except ValueError:
  print "Value error caught while solving for the iterative solution."
  simfail = True
except FloatingPointError:
  print "Run time error from overflow."
  simfail = True

refSol = data.alphamat.copy()
numerr = data.fixed_point_error(data.Nfk) # Plot for line.

# Solve the perturbative system saving the fixed point error for the
# solution upto order i.
pertsoldiff = []

for i in xrange(1, Ord_req) :
  data.Coefcell = []
  data.Symcell = []
  data.Ord_req = i
  data.solve()
  data.fouriersol()
  talpha = np.zeros((latsize**2,Ord_req-1)).astype(complex)
  talpha[:,:i] = data.FSol.T
  newdiff = data.fixed_point_error(talpha)
  pertsoldiff.append(newdiff)

# Solve the system using the iterative method with maxiter = i and save the error.
itersoldiff = []

for i in xrange(1, 21) :
  data.iterative_solver(_maxiter = i,_tol = 3.0e-16)
  itersoldiff.append(data.fixed_point_error(data.alphamat))
  if itersoldiff[i-1] < 1e-16:
    break

fixedpointerror = np.zeros((len(itersoldiff),2))
fixedpointerror[:,0] = pertsoldiff[:len(itersoldiff)]
fixedpointerror[:,1] = itersoldiff

# Write out data for plot of fixed point error
# The first line containds the error in numerical method.
f = open('fixedpoint.txt','w')
f.write('{0:g}\t{1:g}\n'.format(0.0,numerr))
for i in xrange(len(itersoldiff)):
  f.write('{0:g}\t{1:g}\n'.format(fixedpointerror[i,0],fixedpointerror[i,1]))
f.close()

# Data for creating the decay of the Fourier coefficients.
fdecay = np.zeros((20, 3))
for i in xrange(20) :
    fdecay[i, 0] = np.log(np.linalg.norm(data.alphamat[:, i]))
    fdecay[i, 1] = np.log(np.linalg.norm(data.FSol[i, :]))
    fdecay[i, 2] = np.log(np.linalg.norm(data.Nfk[:, i]))

f = open('fdecay.txt','w')
np.savetxt(f,fdecay,delimiter='\t')
f.close()

# Create the plots shown in figures.
