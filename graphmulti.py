#!/usr/bin/env python

from __future__ import division
import numpy as np
import itertools as it
import linearcapcode as cc
import multiprocessing as mltp
import cPickle as pickle
import time as time
import os

np.seterr(all='raise')

def driverfunc(idnum, (size, num_wanteigs, delta, forc_amp, omega_ratio),
               graphtype, R_premaxvolt, R_postmaxvolt,
               R_prepertEDiff, R_postpertEDiff,
               R_preiterEDiff, R_postiterEDiff,
               R_failurenum,
               R_minL, R_maxL,
               R_oldEnergy, R_newEnergy,
               R_newtoniter,
               R_totnorm,
               R_oldconc, R_newconc) :

    # Saved variables.
    minL = 0.0; maxL = 0.0;
    newEnergy = 0.0; oldEnergy = 0.0;
    prepertEDiff = 0.0; postpertEDiff = 0.0;
    preiterEDiff = 0.0; postiterEDiff = 0.0;
    premaxvolt = 0.0; postmaxvolt = 0.0;
    failures = 0; totnorm = 0.0
    newtoniter = 0;
    oldconc = 0.0; newconc = 0.0;

    numrepeats = 100
    maxrepeats = 2*numrepeats
    repeats = 0
    numsuccess = 0

    while (((repeats < maxrepeats)) and (not(numsuccess == numrepeats))) :
        if repeats != 0 :
            del data
        repeats += 1
        simfail = False
        forc_coef = np.array((0,1))

        if graphtype == 'erdos':
            gtype = 'erdos_renyi_graph'
            param = dict(n = size, p = 0.25)
        elif graphtype == 'barabasi':
            gtype = 'barabasi_albert_graph'
            param = dict(n = size, m = 3)
        else:
            gtype == 'connected_watts_strogatz_graph'
            param = dict(n = size, k = 5, p = 0.3, tries = 200)

        data = cc.Mygraph(gtype, param)
        Cap0 = np.ones(data.N)
        G = 0.15 * np.ones(data.N)

        if graphtype == 'erdos':
            G = 0.5 * np.ones(data.N)

        # Number of input nodes is size/10. The last nodes are always chosen as input.
        # This does not make any difference since the nodes can be renumbered.
        numinputs = int(np.floor(data.N / 10))
        Input_nodes = np.zeros(data.N)
        Input_nodes[-numinputs:] = 1

        Eps = 0.5
        L = np.ones(data.rawBmat.shape[1] + sum(Input_nodes))
        Ord_req = 10

        tpoints = 501
        data.init_graph(Cap0, Eps, L, G, Input_nodes, forc_amp, forc_coef, Ord_req, tpoints)

        # iterative solution
        data.iterative_solver()
        simfail = data.simfail

        if simfail == False :                     # If there were no errors continue with simulation.
            # Save all 'pre' data.
            data.solve()
            data.fouriersol()
            oldpertEDiff = data.pertEDiff.real
            olditerEDiff = data.iterEDiff.real
            oldmaxvolt = np.max(data.itersolmat.real)
            nummodes = data.alphamat.shape[1]
            oldconcvec = np.sum(np.abs(data.alphamat[:,1:nummodes]),1) / np.sum(np.abs(data.alphamat),1)
            oldenergy = np.zeros(nummodes)

            for j in range(nummodes):
                oldenergy[j] = np.linalg.norm(data.alphamat[:,j])

            wanteigs = data.Omega**2 * np.ones(num_wanteigs)
            multfac = np.ones(len(wanteigs))
            multfac[0] = delta
            wanteigs = multfac * wanteigs

            # Optimization.
            data.optimalL(wanteigs)

            if data.newtoniter == 2000 :
                simfail = True

            if simfail == False:
                data.init_graph(Cap0, Eps, data.L, G, Input_nodes, forc_amp, forc_coef, Ord_req, tpoints, Omega = data.Omega * omega_ratio)
                # iterative solution
                data.iterative_solver()
                simfail = data.simfail
                data.solve()
                data.fouriersol()

            else :
                failures += 1

            # Setup basic check so that the statistics do not get skewed.
            TOL = 1e-6
            if simfail == False :
                if ((oldpertEDiff > TOL) or \
                        (olditerEDiff > TOL) or \
                        (np.max(data.L) > 1e+5) or \
                        (oldmaxvolt > 2) or \
                        (np.max(data.itersolmat.real) > 2) or \
                        (data.pertEDiff > TOL) or \
                        (data.iterEDiff > TOL) or \
                        (data.totnorm > TOL)) : simfail = True

            if simfail == False :
                newconcvec = np.sum(np.abs(data.alphamat[:,1:nummodes]),1) / np.sum(np.abs(data.alphamat),1)
                newenergy = np.zeros(nummodes)
                for j in range(nummodes):
                    newenergy[j] = np.linalg.norm(data.alphamat[:,j])

                # If the run reaches here, then all is well.
                postpertEDiff += data.pertEDiff.real
                postiterEDiff += data.iterEDiff.real
                postmaxvolt += np.max(data.itersolmat.real)
                premaxvolt += oldmaxvolt
                prepertEDiff += oldpertEDiff
                preiterEDiff += olditerEDiff
                oldEnergy += np.mean(oldconcvec)
                newEnergy += np.mean(newconcvec)
                maxL += np.max(data.L)
                minL += np.min(data.L)
                totnorm += data.totnorm
                newtoniter += data.newtoniter
                oldconc += sum(oldenergy[1:nummodes]) / sum(oldenergy)
                newconc += sum(newenergy[1:nummodes]) / sum(newenergy)

                # Successful trial
                numsuccess += 1

            else :
                failures += 1

        else :
            failures += 1

    if numsuccess == numrepeats :
        R_premaxvolt[idnum] = premaxvolt / numrepeats
        R_postmaxvolt[idnum] = postmaxvolt / numrepeats
        R_prepertEDiff[idnum] = prepertEDiff / numrepeats
        R_postpertEDiff[idnum] = postpertEDiff / numrepeats
        R_preiterEDiff[idnum] = preiterEDiff / numrepeats
        R_postiterEDiff[idnum] = postiterEDiff / numrepeats
        R_failurenum[idnum] = failures
        R_maxL[idnum] = maxL / numrepeats
        R_minL[idnum] = minL / numrepeats
        R_oldEnergy[idnum] = oldEnergy / numrepeats
        R_newEnergy[idnum] = newEnergy / numrepeats
        R_newtoniter[idnum] = int(newtoniter / numrepeats)
        R_totnorm[idnum] = totnorm / numrepeats
        R_oldconc[idnum] = oldconc / numrepeats
        R_newconc[idnum] = newconc / numrepeats
        return

    else :
        R_premaxvolt[idnum] = 0
        R_postmaxvolt[idnum] = 0
        R_prepertEDiff[idnum] = 0
        R_postpertEDiff[idnum] = 0
        R_preiterEDiff[idnum] = 0
        R_postiterEDiff[idnum] = 0
        R_failurenum[idnum] = failures
        R_maxL[idnum] = 0
        R_minL[idnum] = 0
        R_oldEnergy[idnum] = 0
        R_newEnergy[idnum] = 0
        R_newtoniter[idnum] = 0
        R_totnorm[idnum] = 0

        return

def main(graphtype = 'barabasi') :

    # Main parameters over which the loop will run.
    FORC_AMP = np.array((0.001, 0.005, 0.01, ))
    SIZE = np.arange(25, 201, 50)
    DELTA = np.linspace(0.25, 0.75, 10)
    OMEGA_RATIO = np.array((1.0,)) # Not used anymore
    NUM_WANTEIGS = np.array((2,))

    numFA = len(FORC_AMP); numS = len(SIZE); numD = len(DELTA);
    numOR = len(OMEGA_RATIO); numNW = len(NUM_WANTEIGS)

    # Total number of unique simulations. Each one will be run 'numrepeats' times and averaged over.
    numsims = numFA * numS * numD * numOR * numNW

    # Iterator to loop over all values.
    ITR = it.product(SIZE, NUM_WANTEIGS, DELTA, FORC_AMP, OMEGA_RATIO)

    # Results that will be saved.
    R_premaxvolt = mltp.Array('d', np.zeros(numsims))
    R_postmaxvolt = mltp.Array('d', np.zeros(numsims))
    R_prepertEDiff = mltp.Array('d', np.zeros(numsims))
    R_preiterEDiff = mltp.Array('d', np.zeros(numsims))
    R_postpertEDiff = mltp.Array('d', np.zeros(numsims))
    R_postiterEDiff = mltp.Array('d', np.zeros(numsims))
    R_failurenum = mltp.Array('i', np.zeros(numsims, dtype = 'i'))
    R_maxL = mltp.Array('d', np.zeros(numsims))
    R_minL = mltp.Array('d', np.zeros(numsims))
    R_oldEnergy = mltp.Array('d', np.zeros(numsims))
    R_newEnergy = mltp.Array('d', np.zeros(numsims))
    R_newtoniter = mltp.Array('i', np.zeros(numsims, dtype = 'i'))
    R_totnorm = mltp.Array('d', np.zeros(numsims))
    R_oldconc = mltp.Array('d', np.zeros(numsims))
    R_newconc = mltp.Array('d', np.zeros(numsims))

    siminfo = []
    numparallel = 10           # Each time we run 'numparallel' different runs.
    thrdlist = []
    numdone = 0

    for numdone in xrange(numparallel) :
        simiter = ITR.next()
        siminfo.append((numdone, simiter))
        plist = mltp.Process(target = driverfunc,
                         args = (numdone, simiter,graphtype,
                                 R_premaxvolt, R_postmaxvolt,
                                 R_prepertEDiff, R_postpertEDiff,
                                 R_preiterEDiff, R_postiterEDiff,
                                 R_failurenum,
                                 R_minL, R_maxL,
                                 R_oldEnergy, R_newEnergy,
                                 R_newtoniter,
                                 R_totnorm,
                                 R_oldconc, R_newconc))
        thrdlist.append(plist)
        plist.start()

    # This loop makes sure there are 10 jobs running at all times unlike
    # before.
    while True :
        if numdone == (numsims - 1) :
            break
        # Check every second.
        time.sleep(1)
        for i in xrange(numparallel) :
            if numdone == (numsims - 1) :
                break
            if not thrdlist[i].is_alive() :
                numdone += 1
                simiter = ITR.next()
                siminfo.append((numdone, simiter))
                plist = mltp.Process(target = driverfunc,
                                     args = (numdone, simiter,graphtype,
                                             R_premaxvolt, R_postmaxvolt,
                                             R_prepertEDiff, R_postpertEDiff,
                                             R_preiterEDiff, R_postiterEDiff,
                                             R_failurenum,
                                             R_minL, R_maxL,
                                             R_oldEnergy, R_newEnergy,
                                             R_newtoniter,
                                             R_totnorm,
                                             R_oldconc, R_newconc))

                thrdlist[i] = plist
                thrdlist[i].start()
    for thrd in thrdlist :
        thrd.join()

    # Save all data into a numpy array.
    # 19 columns of saved data for numsims rows.
    # Columns are :
    # Simulation ID : IDNUM
    # Inputs : 5 different INPUTS
    # Maximum voltages before optimization : R_premaxvolt, R_postmaxvolt
    # Energy balance before optimization (perturbative) : R_prepertEDiff, R_postpertEDiff
    # Energy balance before optimization (iterative) : R_preiterEDiff, R_postiterEDiff
    # Number of failures for each set of inputs : R_failurenum
    # Minimum / Maximum inductance after optimization : R_minL, R_maxL
    # Mean energy before and after optimization : R_oldEnergy, R_newEnergy
    # Number of iterations for the Newton method : R_newtoniter
    # The residual in the Newton method : R_totnorm,
    # Percentage energy in the higher harmonics : R_oldconc, R_newconc
    results = np.zeros((numsims, 21))

    results[:, 20] = np.array(R_newconc)
    results[:, 19] = np.array(R_oldconc)
    results[:, 18] = np.array(R_totnorm)
    results[:, 17] = np.array(R_newtoniter)
    results[:, 16] = np.array(R_newEnergy)
    results[:, 15] = np.array(R_oldEnergy)
    results[:, 14] = np.array(R_maxL)
    results[:, 13] = np.array(R_minL)
    results[:, 12] = np.array(R_failurenum)
    results[:, 11] = np.array(R_postiterEDiff)
    results[:, 10] = np.array(R_preiterEDiff)
    results[:,  9] = np.array(R_postpertEDiff)
    results[:,  8] = np.array(R_prepertEDiff)
    results[:,  7] = np.array(R_postmaxvolt)
    results[:,  6] = np.array(R_premaxvolt)

    for i in xrange(numsims) :
        results[i , 0] = i
        results[i, 1:6] = siminfo[i][1]

    filename = graphtype + '.txt'
    np.savetxt(filename,results,delimiter=',')

    return results
if __name__ == '__main__' :

    # Run simulations.
    # The choice for different graph types are
    # 1. 'barabasi' : Default
    # 2. 'watts_strogatz'
    # 3. 'erdos'

    # The function main will run ALL the simulations for the given graph type.
    graph = 'barabasi'
    main(graph)

    # Analyse results using R.
    # The code calls the multiplot.R script which plots the results similar to Section 7.
    resultFileName = graph + '.txt'
    command = 'R CMD BATCH --no-save --no-restore \'--args {0:s}\' analyse.R'.format(resultFileName)
    retcode = os.system(command)
