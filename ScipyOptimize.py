#!/usr/bin/env python

import numpy as np
from scipy.optimize import linprog
from numpy.linalg import svd

accuracy = 10e-9

def scipy_optimize_linprog(m, DmAk):
    #print('Python function scipy_optimize_linprog() called')

    DmAk = DmAk.reshape((m, m))

    ALP = np.c_[ -np.eye(m), np.ones((m, 1)) ]
    B = np.zeros((m, 1))
    Aeq = np.r_[ DmAk, np.ones((1, m)) ]
    Aeq = np.c_[ Aeq, np.zeros((m+1, 1)) ]
    beq = np.r_[ B, np.matrix([1]) ]
    ff = np.c_[ np.zeros((1, m)), np.matrix([-1]) ]
    ff = np.array(ff.tolist()[0])

    res = linprog(c=ff, A_ub=ALP, b_ub=B, A_eq=Aeq, b_eq=beq)
    
    if res.success is False:
        return (-1, -1)

    xstar = res.x[:-1]
    g = res.fun

    return (xstar.tolist(), g)

def scipy_optimize_null(m, A, atol=accuracy, rtol=0):
    A = A.reshape((m, m))

    A = np.atleast_2d(A)
    u, s, vh = svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T

    # Matrix NULL as list
    #print(ns.flatten().tolist())

    nn = ns.shape[1]
    b = ns[:,0]

    return (b.tolist(), nn)

def scipy_optimize_chol(m, Hk):
    Hk = Hk.reshape((m, m))

    try:
        np.linalg.cholesky(Hk + np.eye(m) * accuracy)
        return 0
    except np.linalg.linalg.LinAlgError:
        return 1
