#!/usr/bin/env python

from __future__ import division
from scipy import linalg as linalg
from scipy.integrate import ode
import numpy as np
import networkx as nx
from inspect import getargspec
import matplotlib.pyplot as plt
import scipy.sparse as sps
from scipy.sparse import linalg as sla
import matplotlib.colors as colors
import matplotlib.cm as cmx
import scipy.fftpack as fftp

np.seterr(all='warn', over='raise')

class Mygraph :
    def __init__(self, graphtype, graphvars) :
        """ Initialize graph.

        graphtype : Based on networkx graph type.
            Ex : 'barabasi_albert'
        graphvars : Dictionary of all parameters for the graph as
        given in the documentation.
            Ex : For 'barabasi_albert' : dict(n = 10, m = 2, seed = 1)
        """

        try :
            gname = eval('nx.' + graphtype)
            allargs = getargspec(gname)
            impargs = allargs.args
            varlist = '('

            for i in impargs :
                if graphvars.has_key(i) :
                    varlist += i + '=' + str(graphvars.get(i)) + ','

            varlist = varlist[:-1]
            varlist += ')'

            while True :
                self.Gp = eval('nx.' + graphtype + varlist)
                if nx.is_connected(self.Gp) :
                    break


        except AttributeError :
            print "Unknown graph. Check name."
            raise
        except TypeError :
            print "Arguments in graphvars do not meet requirements."
            raise

        self.N = self.Gp.number_of_nodes()
        self.graphvars = graphvars
        self.graphtype = graphtype
        self.graphvars = graphvars
        self.rawBmat = nx.incidence_matrix(self.Gp, oriented = True)

    def init_graph(self, Cap0, Eps, L, G, Input_nodes, forc_amp, forc_coef , Ord_req, tpoints = 51, Omega = -1.0, gtype = 'random') :
        """ Set system parameters.

        Cap0 : array of C0 values for the capacitor (length = N)
        Eps : epsilon value
        L : array of values for inductors. (length = Bmat.shape[1])
        G : array of values for resistors. (length = N)
        Input_nodes : Binary array to indicate which nodes are input. (length = N)
        forc_amp : Corresponds to 'A' in the forcing term.
        forc_coef : Corresponds to x,y in 'x + jy'.
        Ord_req : Solution will be found upto Ord_req orders.
        Omega : if > 0, this value is used.
                else smallest eigenvalue of graphlap will be used.
        Can be overwritten as well.

        The Bmat matrix has been converted to a sparse format.
        Currently two eigenvalues are calculated.

        To call the code for a random graph use gtype = 'random'.
        This will do all the required work and create all the variables.

        For the square grid graph with forcing on the left and bottom boundaries
        use gtype = 'grid'. Remember to override the Omega variable.
        If you do override the Omega variable remember to run the following
        code afterwards:
        _times = np.linspace(0.0, 2.0*np.pi / self.Omega, tpoints)
        _times = _times[:-1]
        self.times = _times


        """
        if gtype is 'random' :
            tempeye = np.eye(self.N)

            # Need this if there are more than 1 inputs to a node.
            def retvec(x, ind) :
                return np.ones(x) * ind
            ind_inp = np.array([])
            for i in xrange(len(Input_nodes)) :
                ind_inp = np.hstack([ind_inp, retvec(Input_nodes[i], i)])
            ind_inp = ind_inp.astype(np.int32)
            self.Bmat = np.hstack((-self.rawBmat, tempeye[:, ind_inp]))
            self.Bmat = sps.lil_matrix(self.Bmat).tocsr()

            try :

                self.Cap0 = Cap0
                self.L = L
                self.G = G
                self.Input_nodes = Input_nodes
                self.Eps = Eps
                self.Ord_req = Ord_req
                self.Coefcell = []
                self.Symcell = []
                self.forc_amp = forc_amp
                self.forc_coef = forc_coef
                self.Graphlap = self.Bmat * sps.diags(1.0/self.L, 0, shape = ((len(self.L), len(self.L)))) *  self.Bmat.T
            except :
                print "Array size for parameters is not consistent in init_graph."

            if Omega > 0 :
                self.Omega = Omega

            else :
                eigvals, evecs = sla.eigsh(self.Graphlap, k = 2, which = 'LM', sigma = 0.0, maxiter = 2000, mode = 'normal')
                eigvals = np.real(eigvals)
                eigvals.sort()
                self.Omega = np.sqrt(eigvals[1])                     # OMEGA

            # move this here so it is independent of the solution method
            # also, this is consistent with the idea of initializing "self"
            # in this function
            _times = np.linspace(0.0, 2.0*np.pi / self.Omega, tpoints)
            _times = _times[:-1]
            self.times = _times

            Pmat = sps.lil_matrix((self.Bmat.shape[0], self.Bmat.shape[1]))
            sparseeye = sps.eye(self.Bmat.shape[0], self.Bmat.shape[0])
            sparseeye = sparseeye.tolil()
            Pmat[:, -int(sum(Input_nodes)) :] = sparseeye[:, ind_inp]
            _temp = np.zeros(self.Bmat.shape[0])
            _temp[Input_nodes > 0] = 1.0
            self.W = np.zeros(self.Bmat.shape[1]).astype(complex)
            self.W = Pmat.T * _temp.reshape(-1, 1) * forc_amp * (forc_coef[0] + 1.0j * forc_coef[1])

            self.Vin = self.Bmat * sps.diags(1.0/self.L, 0, shape = ((len(self.L), len(self.L)))) * self.W.reshape(-1, 1)
            return

        elif gtype is 'grid' :
            tempeye = np.eye(self.N)
            latsize = int(np.sqrt(len(Cap0)))
            # Need this is there are more than 1 inputs to a node.
            def retvec(x, ind) :
                return np.ones(x) * ind
            ind_inp = np.array([])
            for i in xrange(len(Input_nodes)) :
                ind_inp = np.hstack([ind_inp, retvec(Input_nodes[i], i)])
            ind_inp = ind_inp.astype(np.int32)

            # Create the Bmat matrix.
            bigmat = sps.lil_matrix((latsize**2, 2*(latsize * (latsize - 1))))

            b1 = sps.spdiags((np.ones(latsize), -np.ones(latsize)), (0, -1), 2*latsize, latsize - 1)
            b2 = sps.spdiags((-np.ones(latsize), np.ones(latsize)), (0, -latsize), 2*latsize, latsize)

            b = sps.hstack((b1, b2))

            for i in xrange(latsize - 1) :
                bigmat[i*latsize:i*latsize + 2*latsize, i*(2*latsize - 1):(i+1)*(2*latsize - 1)] = b

            bigmat[-latsize:, -latsize+1:] = sps.spdiags((np.ones(latsize), -np.ones(latsize)), (0, -1), latsize, latsize-1)

            # Fit the input nodes to the bottom and left boundaries.
            leftbound = sps.spdiags((np.ones(latsize)), 0, latsize**2, latsize)
            botbound = sps.identity(latsize**2).tocsr()

            botbound = botbound[:, np.arange(2 * latsize, latsize**2 + 1, latsize) - 1]

            self.Bmat = sps.hstack((bigmat, leftbound, botbound)).tocsr()

            try :
                self.M = latsize            # The size of the lattice.
                self.Cap0 = Cap0
                self.L = L
                self.G = G
                self.Input_nodes = Input_nodes
                self.Eps = Eps
                self.Ord_req = Ord_req
                self.Coefcell = []
                self.Symcell = []
                self.forc_amp = forc_amp
                self.forc_coef = forc_coef
                self.Graphlap = self.Bmat * sps.diags(1.0/self.L, 0, shape = ((len(self.L), len(self.L)))) *  self.Bmat.T

            except :
                print "Array size for parameters is not consistent in init_graph."

            if Omega < 0 :
                eigvals, evecs = sla.eigsh(self.Graphlap, k = 2, which = 'LM', sigma = 0.0, maxiter = 2000, mode = 'normal')

                eigvals = np.real(eigvals)
                eigvals.sort()
                self.Omega = np.sqrt(eigvals[1])                     # OMEGA
            else :
                self.Omega = Omega

            _times = np.linspace(0.0, 2.0*np.pi / self.Omega, tpoints)
            _times = _times[:-1]
            self.times = _times

            Pmat = sps.lil_matrix((self.Bmat.shape[0], self.Bmat.shape[1]))
            sparseeye = sps.eye(self.Bmat.shape[0], self.Bmat.shape[0])
            sparseeye = sparseeye.tolil()
            Pmat[:, -int(sum(Input_nodes)) :] = sparseeye[:, ind_inp]
            _temp = np.zeros(self.Bmat.shape[0])
            _temp[Input_nodes > 0] = 1.0
            self.W = np.zeros(self.Bmat.shape[1]).astype(complex)
            self.W = Pmat.T * _temp.reshape(-1, 1) * forc_amp * (forc_coef[0] + 1.0j * forc_coef[1])

            self.Vin = self.Bmat * sps.diags(1.0/self.L, 0, shape = ((len(self.L), len(self.L)))) * self.W.reshape(-1, 1)
            return

    def Malpha(self, RHS, alpha) :
        """ Forms and solves the LU system.


        Returns the new coeffcient. Does not automatically update any value.
        Internally called by method - solve.
        """
        Mmat = -alpha**2 * np.diag(self.Cap0) + (1.0j * alpha * np.diag(self.G)) + self.Graphlap
        coef = linalg.solve(Mmat, RHS)
        return coef

    def Higher_order(self, itr, _ord = 0) :
        """ Form higher order solutions.

        Will not return anything. Will update Coefcell
        and Symcell before exit. Internally called by method - solve.
        """

        _matcols = np.concatenate((np.array([2]), np.arange(5, (5 + ((self.Ord_req + _ord) - 1) * 2), 2)))
        _coef = np.zeros((self.N, _matcols[itr])).astype(complex)
        _sym = np.arange(-np.floor(_matcols[itr]/2), np.ceil(_matcols[itr]/2))

        # print _matcols

        for m in xrange(itr) :
            _mcols = _matcols[m]
            _kcols = _matcols[itr - m - 1]

            for ii in xrange(_mcols) :
                for jj in xrange(_kcols) :

                    _index = self.Symcell[m][ii] + self.Symcell[itr - m - 1][jj] + np.floor(_matcols[itr]/2)
                    _index = int(_index)
                    _coef[:, _index] += self.Cap0 \
                                                * ((1.0j * self.Omega * self.Symcell[m][ii]) * self.Coefcell[m][:, ii] \
                                                * (1.0j * self.Omega * self.Symcell[itr - m - 1][jj]) \
                                                * self.Coefcell[itr - m - 1][:, jj] \
                                                + self.Coefcell[m][:, ii] \
                                                * ((1.0j * self.Omega * self.Symcell[itr - m - 1][jj] ) ** 2 \
                                                   * self.Coefcell[itr - m - 1][:, jj]))

        lensym = len(_sym)
        lenlow = int(np.floor(lensym/2))
        for kk in xrange(lenlow, lensym) :
            _coef[:, kk] = self.Malpha(_coef[:, kk], _sym[kk] * self.Omega)

            if kk != lenlow :   # This corresponds to everything except 0.
                _coef[:, lensym - kk - 1] = _coef[:, kk].conjugate()

        self.Coefcell.append(_coef)
        self.Symcell.append(_sym)


    def solve(self) :
        """ Main solve function.

        Solution will be found up to order self.Ord_req. Internally calls
        both Malpha and higher_order.
        """

        _temp = np.array(self.Malpha(self.Vin.T.ravel(), self.Omega))   # Order 0 solution.
        _temp = np.vstack((_temp.conjugate(), _temp)).T
        self.Coefcell.append(_temp)
        self.Symcell.append(np.array([-1, 1]))

        for itr in xrange(1, self.Ord_req) :
            self.Higher_order(itr, _ord = 0)


    def extend(self, _ord = 1) :
        """ Finds the solution one more order (default), or upto '_ord' more orders.

        It automatically updates both the time_sol and the fouriersol.
        """
        for itr in xrange(self.Ord_req, self.Ord_req + _ord) :
            self.Higher_order(itr, _ord)

        self.timesol(more_ord = _ord)

    def timesol(self, more_ord = 0) :
        """ Finds the solution in the time domain.

        tpoints : Number of time points to evaluate the solution at.
        more_ord : Used to extend the order by using all the lower order calculations.
        """
        _epsvec = [1]

        for itr in xrange(1, self.Ord_req + more_ord) :
            _epsvec.append(_epsvec[itr - 1] * self.Eps)

        if more_ord == 0 :
            self.Sol = np.zeros((len(self.times), self.N), ).astype(complex)

            for i in xrange(len(self.times)) :
                _tsol = np.zeros(self.N).astype(complex)

                for j in xrange(len(_epsvec)) :
                    _tsol += _epsvec[j] * (np.sum(self.Coefcell[j] \
                                                      * (np.exp(1.0j * self.Omega * self.Symcell[j] * self.times[i])), axis = 1))

                    self.Sol[i,:] = _tsol
            # Make the solution real.
            self.Sol = self.Sol.real

        elif more_ord > 0 :
            for i in xrange(len(self.times)) :
                _tsol = np.zeros(self.N).astype(complex)

                for j in xrange(len(_epsvec) - self.Ord_req) :
                    _tsol += _epsvec[j] * (np.sum(self.Coefcell[j] \
                                                      * (np.exp(1.0j * self.Omega * self.Symcell[j] * self.times[i])), axis = 1))

                self.Sol[i,:] = _tsol

            # print "Order of solution increased from {} to {}".format(self.Ord_req, self.Ord_req + more_ord)
            self.Ord_req += more_ord # Updating the Ord_req to get to current order.

            # Update Fourier solution as well.
            self.fouriersol(more_ord = more_ord)

            # Make the solution real.
            self.Sol = self.Sol.real

        else :
            print "Order to increase cannot be negative."

        # maxvoltage = 1.0/self.Eps
        # if (np.max(np.abs(self.Sol.real)) > maxvoltage) :
        #     print "Solution is not right because capacitance function has crossed 0."

    def optimalL(self, deseigs = np.zeros(0), maxiters = 2000) :
        """Finds the optimal values of L using Newton's method.

        Requires an input of Omega values which are smaller than the number of nodes.
        The actual vector can be smaller than self.N.
        """

        if deseigs.shape[0] == 0 :
            deseigs = np.array([0.75*data.Omega**2, data.Omega**2]).flatten()

        _ind = np.argsort(deseigs)
        deseigs = deseigs[_ind[range(0, len(deseigs), 1)]] # To sort in ascending order.
        # deseigs = deseigs[_ind[range(len(deseigs) - 1, -1, -1)]] # To sort in descending order.

        # Newton loop
        Lvec = self.L
        ndes = len(deseigs)
        Bmat = self.Bmat
        numedges = Bmat.shape[1]
        normvec = np.zeros(maxiters)
        mycounter = 0

        ### TEMP
        # normvec = []

        for i in xrange(maxiters):

            # calculate current graph Laplacian
            # use actualLvec to encode the positivity constraint on Lvec
            # it is obvious that all entries of actualLvec will be >= 0
            actualLvec = Lvec**2
            graphlap = Bmat * sps.diags(actualLvec, 0, shape = ((numedges, numedges))) * Bmat.T

            # now figure out eigenpairs
            # d, v = linalg.eig(graphlap)

            # Just to debug the code.
            # graphlap = sps.lil_matrix(graphlap)
            # graphlap = sps.csc_matrix(graphlap)
            try :
                d, v = sla.eigsh(graphlap, k = ndes, which = 'LM', maxiter = 2000, sigma = 0.0, mode = 'normal')

            except :
                self.newtoniter = 2000
                return

            # we can do this because we *know* that graphlap is symmetric
            d = d.real
            numeigs = len(d)
            # Sorted version.
            sortid = np.argsort(d)
            v = v[:, sortid[:ndes]] # Pick the evecs and evals corresponding to the smallest evals.
            d = d[sortid[:ndes]]

            # Using the faster version. Can be verified that they both produce same results.
            _jac = np.zeros((ndes, numedges)).astype(complex)

            for j in xrange(ndes) : # column increase
                Btj = np.asarray(Bmat.T * np.asmatrix(v[:, j]).T).flatten()
                for k in xrange(numedges) : # row increase
                    _jac[j, k] = Btj[k] * (2.0*Lvec[k]) * Btj[k]

            # handle inequality constraints
            # user-defined lower and upper bounds
            _lb = 1.0e-3
            _ub = 1.0e3

            # this bit of code exists merely to record which entries of Lvec
            # violate the lower/upper bound constraints
            badsmall = []
            badbig = []
            for j in xrange(numedges):
                if (Lvec[j] < _lb):
                    badsmall.append(j)
                elif (Lvec[j] > _ub):
                    badbig.append(j)

            # it is useful to know how many constraints are violated, so that
            # we can correctly allocate space for constraints and Jconstraints
            constraints = np.zeros(len(badsmall) + len(badbig))
            Jconstraints = np.zeros((len(constraints), numedges))

            # now we fill in both the constraint vector and its Jacobian
            ctr = 0
            # lb constraints go first
            for j in badsmall:
                constraints[ctr] = 0.5*(Lvec[j]**2 - _lb)**2
                Jconstraints[ctr,j] = (Lvec[j]**2 - _lb)*2.0*Lvec[j]
                ctr += 1

            # ub constraints go second
            for j in badbig:
                constraints[ctr] = 0.5*(Lvec[j]**2 - _ub)**2
                Jconstraints[ctr,j] = (Lvec[j]**2 - _ub)*2.0*Lvec[j]
                ctr += 1

            # stack Jacobians
            _jac = np.vstack([_jac, Jconstraints])

            # compute F
            # mismatch between obtained and desired eigs
            func = d - deseigs

            # stack the mismatch on top of the constraints,
            # matching the stacking of the Jacobian above
            funcnorm = linalg.norm(func)
            totnorm = funcnorm
            if len(constraints) > 0:
                constrnorm = linalg.norm(constraints)
                totnorm += constrnorm

            self.totnorm = totnorm
            # normvec.append(totnorm)
            func = np.hstack([func, constraints])

            # termination criterion
            if (totnorm < 1.0e-12):
                break

            # take a step
            # note that there is no need to compute the full pinv,
            # if we are just going to multiply it by a vector anyway
            step, _resid, _rank, _sigma = linalg.lstsq(_jac, func)
            step = step.flatten(1)
            if linalg.norm(step.imag) < 1e-14 :
                step = np.real(step)
                Lvec -= 0.15 * step
                mycounter += 1

            else :
                mycounter = 2000
                break

        # self.nvec = normvec
        self.newtoniter = int(mycounter)
        self.L = 1.0 / actualLvec
        self.Graphlap = Bmat * sps.diags(actualLvec, 0, shape = ((numedges, numedges))) * Bmat.T

    def fouriersol(self, more_ord = 0) :
        """Finds the solution in Fourier space. Also computes the L2 norm.

        Sums up the Coefcell by considering the repective powers of \epsilon.
        NOTE: Computes only the +k\omega coeffcients.
        L2 norm is found for each node. Will help in finding
        which nodes to end to
        transfer energy into the higher harmonics better.
        NOTE: np.diff(self.l2norm, axis = 1) will be helpful here.
        """
        _epsvec = [1.0]

        # When we reach here, if called from extend, we have already updated the Ord_req.
        for itr in xrange(1, self.Ord_req) :
            _epsvec.append(_epsvec[itr - 1] * self.Eps)

        _maxcols = len(self.Symcell[self.Ord_req - 1][self.Symcell[self.Ord_req - 1] > 0])

        if more_ord == 0 :

            self.FSol = np.zeros((_maxcols, self.N)).astype(complex)
            self.l2norm = np.zeros((self.N, len(_epsvec)))
            self.ordnorm = np.zeros((self.N, len(_epsvec)))

            for i in xrange(self.Ord_req) :
                _tsol = np.zeros((self.N, _maxcols)).astype(complex)
                _poscols = self.Symcell[i] > 0
                _tsol = np.hstack((self.Coefcell[i][:, _poscols], np.zeros((self.N, _maxcols - sum(_poscols)))))
                self.FSol += _epsvec[i] * _tsol.T
                self.ordnorm[:, i] = np.real(np.sum(self.Coefcell[i]*self.Coefcell[i].conjugate(),1))
                self.l2norm[:, i] = np.real(np.sum(self.FSol.T * self.FSol.conjugate().T, 1))

        else :
            # Using pre-computed solutions.
            _nFSol = np.zeros((_maxcols, self.N)).astype(complex)
            _oldmaxcols = len(self.Symcell[self.Ord_req - 1 - more_ord][self.Symcell[self.Ord_req - 1 - more_ord] > 0])
            _nFSol[:_oldmaxcols, :] = self.FSol
            _nl2norm = np.zeros((self.N, len(_epsvec)))
            _nl2norm[:, :self.Ord_req - more_ord] = self.l2norm
            _nordnorm = np.zeros((self.N, len(_epsvec)))
            _nordnorm[:, :self.Ord_req - more_ord] = self.ordnorm

            # Updating part.
            for i in xrange(self.Ord_req - more_ord, self.Ord_req) :
                _tsol = np.zeros((self.N, _maxcols)).astype(complex)
                _poscols = self.Symcell[i] > 0
                _tsol = np.hstack((self.Coefcell[i][:, _poscols], np.zeros((self.N, _maxcols - sum(_poscols)))))
                _nFSol += _epsvec[i] * _tsol.T
                _nordnorm[:, i] = np.real(np.sum(self.Coefcell[i]*self.Coefcell[i].conjugate(),1))
                _nl2norm[:, i] = np.real(np.sum(self.FSol.T * self.FSol.conjugate().T, 1))

            self.FSol = _nFSol
            self.l2norm = _nl2norm
            self.ordnorm = _nordnorm

        self.maxmode = np.zeros(self.N)
        for j in range(self.N):
            self.maxmode[j] = np.argmax(np.abs(self.FSol[:, j]))

        self.energy_cons(method = 'perturbative')

    def energy_cons(self, method = 'perturbative') :
        """ Compute the energy put in to the system and the energy flowing out.
        UPDATE: Does the energy balance for both the perturbative solution and
        the iterative solution. Mention the type in the input.
        method = 'perturbative' or 'iterative'

        Computes the current coefficient of +\omega only at the Input_nodes,
        since this is all that is required to compute energyin and energyout.
        """

        if method == 'perturbative' :
            _cur = np.zeros((sum(self.Input_nodes > 0), 2)).astype(complex)
            _vtotal = np.zeros((sum(self.Input_nodes > 0))).astype(complex)

            _epsvec = [1.0]

            for itr in xrange(1, self.Ord_req) :
                _epsvec.append(_epsvec[itr - 1] * self.Eps)

            for i in xrange(self.Ord_req) :
                _vtotal += _epsvec[i] * self.Coefcell[i][self.Input_nodes > 0, self.Symcell[i] == 1]

            _vin = self.Vin[self.Input_nodes > 0]
            _vin = _vin.ravel()

            _cur[:, 1] = (_vin - _vtotal) / (1.0j * self.Omega)
            _cur[:, 0] = _cur[:, 1].conjugate()
            self.Cur = _cur

            _energyin = np.sum(_cur * np.vstack((_vin, _vin.conjugate())).T)

            _tfourier = self.FSol[:, self.G > 0]
            _energyout = 2 * np.sum(self.G * (_tfourier * _tfourier.conjugate())) # 2 because FSol only contains +ve frequencies.

            self.pertEDiff = (_energyin - _energyout).real

        elif method == 'iterative' :
            _cur = np.zeros((sum(self.Input_nodes > 0), 2)).astype(complex)
            _vtotal = self.alphamat[self.Input_nodes > 0, 0]

            _vin = self.Vin[self.Input_nodes > 0]
            _vin = _vin.ravel()

            _cur[:, 1] = (_vin - _vtotal) / (1.0j * self.Omega)
            _cur[:, 0] = _cur[:, 1].conjugate()
            self.Cur = _cur

            _energyin = np.sum(_cur * np.vstack((_vin, _vin.conjugate())).T)

            _tfourier = self.alphamat[self.G > 0, :].T
            _energyout = 2 * np.sum(self.G * (_tfourier * _tfourier.conjugate()))

            self.iterEDiff = _energyin - _energyout

        elif method == 'numerical' :
            _cur = np.zeros((sum(self.Input_nodes > 0), 2)).astype(complex)
            _vtotal = self.Nfk[self.Input_nodes > 0, 0]

            _vin = self.Vin[self.Input_nodes > 0]
            _vin = _vin.ravel()

            _cur[:, 1] = (_vin - _vtotal) / (1.0j * self.Omega)
            _cur[:, 0] = _cur[:, 1].conjugate()
            self.Cur = _cur

            _energyin = np.sum(_cur * np.vstack((_vin, _vin.conjugate())).T)

            _tfourier = self.Nfk[self.G > 0, :].T
            _energyout = 2 * np.sum(self.G * (_tfourier * _tfourier.conjugate()))

            self.NEDiff = _energyin - _energyout

    def iterative_solver(self, _maxmode = 20, _maxiter = 200):

        # internal tolerance
        _mytol = 3.0e-16

        _Vin = self.Vin.ravel()

        # set up two matrices
        alphamat = 1.0j*np.zeros((self.N, _maxmode))

        # precompute LU factorizations
        lumat = []
        pivmat = []
        for k in range(1,_maxmode+1):
            alpha = k*self.Omega
            _Mmat = -alpha**2 * np.diag(self.Cap0) + (1.0j * alpha * np.diag(self.G)) + self.Graphlap
            [_lu, _piv] = linalg.lu_factor(_Mmat, overwrite_a=True)
            lumat.append(_lu)
            pivmat.append(_piv)

        self.simfail = False
        # To count the number of iterations.
        # iterdrop = []
        # errornorm = []

        try :
            for i in range(_maxiter):

                alphamatnew = 1.0j*np.zeros((self.N, _maxmode))
                # k is an actual Fourier mode number, no off-by-1 correction
                for k in range(1,_maxmode+1):
                    rhs = 1.0j*np.zeros(self.N)
                    if (k==1):
                        rhs += _Vin

                    convolve = 1.0j*np.zeros(self.N)
                    # l is an actual Fourier mode number, no off-by-1 correction
                    for l in range(-_maxmode,_maxmode+1):
                        if (l > 0):
                            # correct off-by-1
                            alphal = alphamat[:,l-1]
                        elif (l < 0):
                            # correct off-by-1
                            alphal = (alphamat[:,-l-1]).conj()
                        else:
                            alphal = np.zeros(self.N)

                        if ( ((k-l) > 0) and ((k-l) <= _maxmode) ):
                            # correct off-by-1
                            alphakml = alphamat[:,k-l-1]
                        elif ( ((k-l) < 0) and ((k-l) >= -_maxmode) ):
                            # correct off-by-1
                            alphakml = (alphamat[:,l-k-1]).conj()
                        else:
                            alphakml = np.zeros(self.N)

                        convolve += alphal*alphakml

                    rhs -= (0.5*self.Eps*(self.Omega**2)*(k**2))*self.Cap0*convolve

                    # correct off-by-1
                    alphamatnew[:,k-1] = linalg.lu_solve((lumat[k-1],pivmat[k-1]), rhs, overwrite_b = False)

                # Check the differences in the alphamat's for each iteration.
                _normcheck = linalg.norm(alphamat - alphamatnew, np.inf)
                alphamat = alphamatnew.copy()

################################################################
                # # # This part of code strictly for calculating the erro between numerical solution and the iterative solution at each iteration.
                # iterdrop.append(_normcheck)
                # [maxmodemat, tptsmat] = np.meshgrid( np.arange(1, _maxmode+1), self.times)
                # mttmat = (np.exp(1.0j*self.Omega*maxmodemat*tptsmat)).T
                # tempsolmat = 2.0*np.real(np.dot(alphamat, mttmat))
                # errornorm.append(linalg.norm(tempsolmat - self.NSol))
################################################################

                if (_normcheck < _mytol):
                    break

        except FloatingPointError :
            print "Iterative method failure."
            self.simfail = True
            return

        # store the final solution
        self.alphamat = alphamat.copy()

        # store the last value of the 2-norm difference
        self.iternormcheck = _normcheck

        # convert to time-domain solution
        [maxmodemat, tptsmat] = np.meshgrid( np.arange(1, _maxmode+1), self.times)
        mttmat = (np.exp(1.0j*self.Omega*maxmodemat*tptsmat)).T
        self.itersolmat = 2.0*np.real(np.dot(self.alphamat, mttmat))
        self.energy_cons(method = 'iterative')
        # self.itererr = errornorm
        # self.iterdrop = iterdrop

        return

    def fixed_point_error(self, alphamat):

        """ Code to check what $|| \alpha^{j} - F(\alpha^{j}) ||$.
        alphamat is the \alphamat from the iterative solution or
        self.FSol from the perturbative solution.

        Make sure that the matrix is $N \times M$ where $M$ is the number of modes.
        Note: This method returns the value and does not store it.
        """
        _Vin = self.Vin.ravel()

        # set up two matrices
        alphamatnew = np.zeros_like(alphamat).astype(complex)
        _maxmode = alphamat.shape[1]

        # precompute LU factorizations
        lumat = []
        pivmat = []
        for k in range(1,_maxmode+1):
            alpha = k*self.Omega
            _Mmat = -alpha**2 * np.diag(self.Cap0) + (1.0j * alpha * np.diag(self.G)) + self.Graphlap
            [_lu, _piv] = linalg.lu_factor(_Mmat, overwrite_a=True)
            lumat.append(_lu)
            pivmat.append(_piv)

        # k is an actual Fourier mode number, no off-by-1 correction
        for k in range(1,_maxmode+1):
            rhs = 1.0j*np.zeros(self.N)
            if (k==1):
                rhs += _Vin

            convolve = 1.0j*np.zeros(self.N)
            # l is an actual Fourier mode number, no off-by-1 correction
            for l in range(-_maxmode,_maxmode+1):
                if (l > 0):
                    # correct off-by-1
                    alphal = alphamat[:,l-1]
                elif (l < 0):
                    # correct off-by-1
                    alphal = (alphamat[:,-l-1]).conj()
                else:
                    alphal = np.zeros(self.N)

                if ( ((k-l) > 0) and ((k-l) <= _maxmode) ):
                    # correct off-by-1
                    alphakml = alphamat[:,k-l-1]
                elif ( ((k-l) < 0) and ((k-l) >= -_maxmode) ):
                    # correct off-by-1
                    alphakml = (alphamat[:,l-k-1]).conj()
                else:
                    alphakml = np.zeros(self.N)

                convolve += alphal*alphakml

            rhs -= (0.5*self.Eps*(self.Omega**2)*(k**2))*self.Cap0*convolve

            # correct off-by-1
            alphamatnew[:,k-1] = linalg.lu_solve((lumat[k-1],pivmat[k-1]), rhs, overwrite_b = False)

        # Check the differences in the alphamat's for each iteration.
        _normcheck = linalg.norm(alphamat - alphamatnew, np.inf)

        return _normcheck

    def firstordersolve(self, maxmode = 10) :
        """Code to solve the same system using the first order equations.
        Solution

        Solution will be found up to Ord_req.
        """
        invlmat = np.diag(1.0/self.L).astype(complex)
        tr = sps.lil_matrix(invlmat * self.Bmat.T)

        ll = sps.diags(1.0/self.Cap0, 0) * sps.lil_matrix(self.Bmat)
        lr = sps.diags(1.0/self.Cap0 * self.G, 0)

        bigmat = sps.bmat([[None, -tr], [ll, -lr]])
        powers = np.arange(maxmode + 1)
        lumat = []

        edgeplusnode = self.Bmat.shape[0] + self.Bmat.shape[1]
        for k in xrange(1, maxmode + 1) :
            submat = sps.diags(np.ones(edgeplusnode), 0) * (1.0j * powers[k] * self.Omega)
            finalmat = sps.csc_matrix(-bigmat + submat)
            decomp = sla.splu(finalmat)
            lumat.append(decomp)

        self.firstordercoef = []
        self.firstordersym = []

        tempeye = np.eye(self.N)
        Pmat = np.hstack((np.zeros((self.N, self.Bmat.shape[1] - sum(self.Input_nodes))), \
                           tempeye[:, self.Input_nodes > 0]))

        forcing = self.forc_amp * ((self.forc_coef[0] * np.ones(self.N) * self.Input_nodes) \
                              + (self.forc_coef[1] * 1.0j * np.ones(self.N) * self.Input_nodes))
        W = np.dot(Pmat.T, forcing)

        orderzero = np.zeros(bigmat.shape[0]).astype(complex)
        orderzero[:len(W)] = 1.0/self.L * W

        # Solve for the order 0 case.
        orderzero = lumat[0].solve(orderzero)
        orderzero = np.vstack((orderzero.conjugate(), orderzero)).T
        self.firstordercoef.append(orderzero)

        self.firstordersym.append([-1, 1])

        def higherorder_firstorder(Cap0, coefcell, itr, symcell, edgeplusnode, _ord = 0) :

            _matcols = np.concatenate((np.array([2]), np.arange(5, (5 + ((self.Ord_req + _ord) - 1) * 2), 2)))
            _coef = np.zeros((edgeplusnode, _matcols[itr])).astype(complex)
            multvec = np.hstack((np.zeros(edgeplusnode - len(Cap0)), 1.0/Cap0 * np.ones(len(Cap0))))

            for m in xrange(itr) :
                _mcols = _matcols[m]
                _ncols = _matcols[itr - 1 - m]
                for ii in xrange(_mcols) :
                    for jj in xrange(_ncols) :

                        _index = symcell[m][ii] + symcell[itr - 1 - m][jj] + np.floor(_matcols[itr]/2)
                        _index = int(_index)
                        _coef[:, _index] +=  multvec * \
                            (coefcell[m][:, ii] * \
                        (1.0j * symcell[itr - 1 - m][jj] * \
                         coefcell[itr - 1 - m][:, jj]))

            lensym = len(symcell[itr])
            lenlow = int(np.ceil(lensym/2))
            for kk in xrange(lenlow, lensym) :
                _coef[:, kk] = lumat[kk - lenlow].solve(_coef[:, kk])
                _coef[:, lensym - kk - 1] = _coef[:, kk].conjugate()

            _coef[:, lenlow - 1] = 0.0
            return _coef

        for itr in xrange(1, self.Ord_req) :
            _coef = higherorder_firstorder(self.Cap0, self.firstordercoef, itr, self.Symcell, edgeplusnode)
            self.firstordercoef.append(_coef)

    def numeric_solver(self, ind = [], _integrator = 'dopri5', _cycles = 200, doplot = False) :
        """ Numerical solution to the problem.
        Code is based on the first order equations.
        The variable self.Nfk is the numerical solution in Fourier space.
        Makes plots for _indices if doplot = True.
        """

        def f_general(t, y):
            v = y[0:self.N]
            i = y[self.N:]
            zdot = np.zeros(len(y))
            c = self.Cap0 * (1.0 - self.Eps * v)
            zdot[:len(v)] = (self.Bmat * i - self.G * v) / c
            vin = 2 * self.W.real *  np.cos(self.Omega * t) - 2 * self.W.imag * np.sin(self.Omega * t)
            zdot[len(v):] = -(self.Bmat.T * v) + vin.ravel()
            return zdot

        r = ode(f_general).set_integrator(_integrator, nsteps = 20000, rtol = 1e-10, atol = 1e-12)

        Tvec = np.linspace(0, _cycles * 2 * np.pi / self.Omega, 40 * _cycles)
        LastT = np.linspace(Tvec[-1], (_cycles + 1) * 2 * np.pi / self.Omega, 65) # Sampling will be done for 64 points
        LastT = LastT[:-1]
        t0 = 0.0;
        y0 = np.zeros(self.N + self.Bmat.shape[1])
        r.set_initial_value(y0, t0)

        # Solutions saved for the last cycle.
        _numsol = np.zeros((self.N, len(LastT)))
        _time = np.zeros(len(LastT))

        # integrate forward in time by numcycles,
        # spinning up the solution to steady-state
        ind = 0
        while r.successful() and r.t < Tvec[len(Tvec) - 1]:
            ind += 1
            r.integrate(Tvec[ind])

        # now store the LastT solution values into numsol.
        for j in xrange(len(LastT)):
            _numsol[:, j] = np.real(r.y[: self.N])
            _time[j] = LastT[j]
            if j < (len(LastT) - 1) :
                r.integrate(LastT[j + 1])

        self.NSol = _numsol
        # Overriding the self.times variables here. The only reason
        # for finding the numerical solution is for comparison. In this
        # case we should first obtain the numerical solution and then
        # compute the iterative/perturbative solutions for these times.
        self.times = _time

        # Obtain numerical solution in Fourier basis.
        fk = fftp.fft(self.NSol)/self.NSol.shape[1]        # Obtains solution for all nodes. The solution is truncated to 20 modes.

        # Pick out only the first 20 +\omega harmonics so that the solution is exactly like the alphamat/FSol shape.
        self.Nfk = fk[:, 1:21]

        self.energy_cons(method = 'numerical')
        if doplot :
            if len(ind) != 0:
                for myvert in ind:
                    plt.figure(4)
                    plt.plot(self.times, np.real(self.Sol[:, myvert]),'r')
                    plt.plot(self.times + dt, _numsol[myvert, :], '*b')
                    plt.legend(('Perturbative', 'Numerical'), shadow = True)
                    plt.show()


    def plot(self, whichplots = {}, whichnodes = [], save2file = True) :
        """ Plot both solution and underlying graph.

        Plots the solution in the time domain along with the underlying graph.
        whichplots whould be a dictionary. If empty then no plot are made.
        The plots made for the keys are :
        1) Key : 'timesol'
        2) Key : 'Graph'
        3) Key : 'modeplot'
        4) Key : 'numericsol'
        Example usage : whichplots = {'timesol' = True, 'Graph' = True}
        will plot just the timesol and the underlying graph.
        'whichnodes' plots only a subset of nodes so it is easier to compare between all three solutions.
        If save2file = False, then plots are generated automatically.
        """


        if whichplots.has_key('timesol') :
            if whichplots['timesol'] == True :
                plt.figure(1)
                for i in xrange(self.N) :
                    if i in set(inp_nodes) :
                        plt.plot(self.times, np.real(self.Sol[:, i]), '-b')
                    else :
                        plt.plot(self.times, np.real(self.Sol[:, i]), '-r')

                        plt.xlim(min(self.times), max(self.times))
                        plt.xlabel('Time')
                        plt.ylabel('Solution')

                        # finalize pdf file
                        if save2file :
                            plt.savefig("timesol.pdf", format = "pdf")
                            plt.close(1)

        if whichplots.has_key('itersol') :
            if whichplots['itersol'] == True :
                inp_nodes = list(np.nonzero(self.Input_nodes)[0])
                oth_nodes = set(self.Gp.nodes()) - set(inp_nodes)
                plt.figure(10)
                for i in xrange(self.N) :
                    if i in set(inp_nodes) :
                        plt.plot(self.times, np.real(self.itersolmat[i,:]), '-b')
                    else :
                        plt.plot(self.times, np.real(self.itersolmat[i,:]), '-r')

                        plt.xlim(min(self.times), max(self.times))
                        plt.xlabel('Time')
                        plt.ylabel('Solution')

                        # finalize pdf file
                        if save2file :
                            plt.savefig("itersol.pdf", format = "pdf")
                            plt.close(10)

        if whichplots.has_key('Graph') :
            if whichplots['Graph'] == True :
                inp_nodes = list(np.nonzero(self.Input_nodes)[0])
                oth_nodes = set(self.Gp.nodes()) - set(inp_nodes)
                self.Gp.add_node(-1)
                inp_nodes.append(-1)
                plt.figure(2)
                shelllist = list(([-1], list(inp_nodes), list(oth_nodes)))
                pos = nx.shell_layout(self.Gp, nlist = shelllist)
                pos[-1] = np.array([0.0, 0.0])

                # pos = nx.spring_layout(self.Gp, iterations = 100)
                # pos[-1] = np.array([0.5, 0.5])

                maxvoltages = np.max(self.itersolmat, axis = 1)
                newmax = np.zeros(len(maxvoltages) + 1)
                newmax[:len(maxvoltages)] = maxvoltages
                newmax[len(newmax) - 1] = np.abs(self.forc_amp)

                dpos = dict()
                for k in inp_nodes :
                    # # Add edges between input and source.
                    if k != -1 :
                        self.Gp.add_edge(-1, k)

                nx.draw_networkx_nodes(self.Gp, pos = pos, nodelist = self.Gp.nodes(),
                                       node_color = newmax, node_cmap = cmx.Blues,
                                       vmin = np.min(newmax), vmax = np.max(newmax),
                                       alpha = 0.7)

                plt.colorbar(orientation='vertical', shrink = 0.75, anchor = (0, 0.85))

                nx.draw_networkx_edges(self.Gp, pos = pos, edge_color = np.log(self.L), edge_cmap = cmx.Oranges, width = 1, alpha = 0.8)
                plt.colorbar(orientation='vertical', shrink = 0.75, anchor = (0, 1.0))
                plt.axis('off')
                plt.suptitle('Color of edges is logarithm of inductance.', fontsize = 16)
                self.pos = pos

                # finalize pdf file
                if save2file :
                    plt.savefig("graph.pdf", format = "pdf")
                    plt.close(2)

        if whichplots.has_key('modeplot') :
            if whichplots['modeplot'] == true :
                plt.figure(3)
                mydeg = np.array(nx.degree(self.Gp).values())
                mydeg = 100 * np.sqrt(mydeg / 2)

                nx.draw_networkx_nodes(self.Gp, pos = pos, node_size = mydeg,
                                       node_color = self.maxmode, alpha = 0.9)

                nx.draw_networkx_edges(self.Gp, pos = pos, width = 0.3,
                                       alpha = 0.5)
                plt.axis('off')

                # finalize pdf file
                if save2file :
                    plt.savefig("modeplot.pdf", format = "pdf")
                    plt.close(3)

        if whichplots.has_key('fouriersol') :
            if whichplots['fouriersol'] == True :
                plt.figure(4)
                for i in xrange(self.N) :
                    plot(range(self.FSol.shape[0]) + 1, self.l2norm[i, :], '-b')
                    plt.xlabel('Mode')
                    plt.ylabel('magnitude')
                    plt.title('Magnitude for each node in the graph.')

                    # finalize pdf file
                    if save2file :
                        plt.savefig("fouriersol.pdf", format = "pdf")
                        plt.close(4)

        if whichplots.has_key('numericsol') :
            if whichplots['numericsol'] == True :
                plt.figure(5)
                for i in xrange(self.N) :
                    if i in set(inp_nodes) :
                        plt.plot(self.times, np.real(self.NSol[i, :]), '-b')
                    else :
                        plt.plot(self.times, np.real(self.NSol[i, :]), '-r')

                        plt.xlim(min(self.times), max(self.times))
                        plt.xlabel('Time')
                        plt.ylabel('Time Stepper Solution')

                        # finalize pdf file
                        if save2file :
                            plt.savefig("Numericsol.pdf", format = "pdf")
                            plt.close(5)

        if whichplots.has_key('comparesol') :
            if whichplots['comparesol'] == True :
                plt.figure(6)
                for i in whichnodes :
                    plt.plot(self.times, np.real(self.NSol[i, :]), '-b')
                    plt.plot(self.times, np.real(self.Sol[1:, i]), '-r', \
                             self.times, np.real(self.Sol[1:, i]), 'r*')
                    plt.plot(self.times, np.real(self.itersolmat[i, :]), '-k')

                    plt.xlim(min(self.times), max(self.times))
                    plt.xlabel('Time')
                    plt.ylabel('All Solutions')
                    plt.legend(('Numerical', 'Perturbative', 'Iterative'), shadow = True)

                # finalize pdf file
                if save2file :
                    plt.savefig("Numericsol.pdf", format = "pdf")
                    plt.close(5)

        plt.show()

    def summary(self) :
        """Function to print some summary of the results.
        """

        print "Difference in energy in the perturbative solution is :\t{}".format(np.abs(self.pertEDiff))
        print "Difference in energy in the iterative solution is    :\t{}".format(np.abs(self.pertEDiff))


        _max = np.max(np.abs(self.Sol.real))
        print "Maximum amplitude in computed solution over all nodes in perturbative solution :\t{}".format(_max)
        _max = np.max(np.abs(self.itersolmat))
        print "Maximum amplitude in computed solution over all nodes in iterative solution    :\t{}".format(_max)

        # print "Solution found upto order {}".format(self.Ord_req)
        # print "Maxmode summary :\t", Counter(self.maxmode + 1)
        # _adjmat = nx.adj_matrix(self.Gp)
        # _adjmat = np.asarray(np.sum(_adjmat, 0))[0]
        # _ind = np.nonzero(self.maxmode)
        # _deg = _adjmat[_ind]
        # print "Degree's for nodes with maxmode > 1 are :\n", _deg
        # _ideg = _adjmat[np.nonzero(self.Input_nodes)[0]]
        # print "Degree's for input nodes are :\n", _ideg

if __name__ == '__main__' :

    # Pick Networkx type graph.
    gtype = 'barabasi_albert_graph'
    # Assign atleast the required parameters for the graph.
    param = dict(n = 25, m = 3)

    # Builds graph.
    data = Mygraph(gtype, param)

    # Create all parameters for input.
    G = 0.05 * np.ones(data.N)                      # Need to change G.

    mydegrees = nx.degree(data.Gp).values()
    maxdegree = np.max(mydegrees)

    for j in range(data.N):
        if (mydegrees[j]==maxdegree):
            G[j] += 1

    Cap0 = G
    Eps = 0.5

    # Make sure there is atleast 1 input node.
    while True:
        Input_nodes = np.random.rand(data.N) > 0.9
        if (sum(Input_nodes) > 0) :
            break

    Input_nodes = [int(x) for x in Input_nodes]
    Input_nodes = np.array(Input_nodes)

    L = 1 * np.ones(data.rawBmat.shape[1] + sum(Input_nodes))
    forc_amp = 1e+0 * 0.2
    forc_coef = np.array((1.0, 1.0))
    Ord_req = 10
    tpoints = 201

    # Initialize the graph.
    data.init_graph(Cap0, Eps, L, G, Input_nodes, forc_amp, forc_coef, Ord_req, tpoints)

    # Solves for the parameters upto order Ord_req.
    print "Solving original problem"
    data.solve()

    # Time domain solution along with the solution in the Fourier domain
    data.timesol()
    data.fouriersol()
    data.iterative_solver()
    data.summary()
    data.firstordersolve()

    # Resolve by finding optimal Lvec.
    omega = data.Omega
    wanteigs = np.array([4 * data.Omega**2, data.Omega**2])

    try :
        # data.optimalL()
        data.optimalL(wanteigs)
        print "Optimal Lvec found."
    except sla.ArpackNoConvergence :
        print "Arpack convergence problem."
        pass

    # Reset the graph with the optimal values. This is very important.
    newdata = Mygraph(gtype, param)
    newdata.Gp = data.Gp
    newdata.N = data.N
    newdata.graphtype = data.graphtype
    newdata.graphvars = data.graphvars
    newdata.rawBmat = nx.incidence_matrix(data.Gp, oriented = True)

    newdata.init_graph(Cap0, Eps, data.L, G, Input_nodes, forc_amp, forc_coef, Ord_req, tpoints, Omega = omega)
    # del data
    newdata.times = np.linspace(0.0, 2.0*np.pi / omega, tpoints)

    # Solves for the parameters upto order Ord_req.
    newdata.solve()
    newdata.timesol()
    newdata.fouriersol()
    #### Numerical solution.
    # newdata.numeric_solver(_cycles = 1000, doplot = False)
    newdata.iterative_solver()
    newdata.summary()

    # Some plots. Note : # The second 'comparesol' requires a calculation of all three solutions (numerical, iterative and perturbative)
    # newdata.plot(whichplots = {'timesol' : True, 'numericsol' : True}, save2file = False)
    # newdata.plot(whichplots = {'comparesol' : True}, whichnodes = [13], save2file = False)

    # Solve the first order system of equations.
    # newdata.firstordersolve()
