import numpy as np
import scipy.optimize as opt
import csv
import matplotlib.pyplot as plt
import sys
from pathlib import Path

# FITPARAMS
# A1 s1 A2 s2 A3 s3 ...

# Dodgy L2 norm of model vs. data
#
FLAG_PLOT = False
FLAG_DEBUG = False

def model(fitparams, local_x):
    Y = np.zeros_like(local_x)
    for i in range(len(fitparams)//2):
        Y += fitparams[i*2]*np.exp(-(local_x/fitparams[i*2+1])**2)
    return Y

def ReadAndFit(filename, legend):
    fdists = np.genfromtxt(filename)
    # These correspond to the meaning of the FormFactor.txt entries themselves
    KMIN = 0
    KMAX = 2
    dim = len(fdists.shape)

    kgrid = np.linspace(KMIN,KMAX,fdists.shape[0 if dim == 1 else 1])

    num_gaussians = 6

    retval = {}
    if dim == 1:
        params = fit_gauss(fdists, kgrid, num_gaussians)
        retval[legend[0]]=params
        if FLAG_PLOT:
            smoothK = np.linspace(KMIN,KMAX,100)
            # make the fit look prettier
            plt.clf()
            plt.title('Configuration 1')
            plt.plot(kgrid, fdists[i], label="data")
            plt.plot(smoothK, model(params, smoothK), label="fit, %d Gaussians" % num_gaussians)
            plt.legend()
            plt.show()


    else:
        for i in range(len(fdists)):
            params = fit_gauss(fdists[i], kgrid, num_gaussians)
            retval[legend[i]]=params
            if FLAG_PLOT:
                smoothK = np.linspace(KMIN,KMAX,100)
                # make the fit look prettier
                plt.clf()
                plt.title('Configuration ' + legend[i])
                plt.plot(kgrid, fdists[i], label="data")
                plt.plot(smoothK, model(params, smoothK), label="fit, %d Gaussians" % num_gaussians)
                plt.legend()
                plt.show()

    return retval

def trapezoid_rint(func, R):
    X = np.append(R, 2*R[-1]-R[-2])
    return np.sum(func*R**2*(X[1:]-X[:-1]))

def fit_gauss(F, K, NGAUSS):
    # minimum spread
    h = 0.0001
    # maximum spread
    H = 100

    # Gaussian weight max
    A = 10

    guess = []
    smallest = []
    largest = []

    def modelint(fitparams):
        I = 0
        for i in range(NGAUSS):
            I += fitparams[i*2]*fitparams[i*2+1]
        return I

    def objective(fitparams, *args):
        return np.sum((model(fitparams, K) - F)**2)

    # Each tuple specifies amplitude, width
    for n in range(1,NGAUSS+1):
        guess.append(((-1)**n * 10/n, 1))
        smallest.append((-A*n, h))
        largest.append((A*n, H))

    if FLAG_DEBUG: print("Guess: ", guess)

    guess_params = np.array(guess).flatten()
    min_params =   np.array(smallest).flatten()
    max_params =   np.array(largest).flatten()

    INTEGRAL =trapezoid_rint(F, K)
    # Changes by factor of sqrt(pi) to compensate
    integralnum = 2*INTEGRAL/1.7724538509

    # NL constraint on the Gaussians
    normalisation = opt.NonlinearConstraint(modelint, integralnum, integralnum)
    # constraints on Gaussian parameter space
    hypercube = opt.Bounds(min_params, max_params)

    res = opt.minimize(objective,guess_params,constraints=[normalisation], bounds=hypercube)


    if not res.success:
        print("Optimisation failure.")
    if FLAG_DEBUG:
        print(res.message)
        print("Termination after ", res.nit, " iterations")
        print("Result: ",res.x)
    return res.x



p = Path('./output/')
paramList = []
for x in p.iterdir():
    if x.is_dir():
        f = x / 'Rates/Form_Factor.txt'
        if f.is_file():
            print("===================================================")
            print("Fitting data for dir " + str(x))
            legend = []
            with open(x / 'index.txt', 'r') as fp:
                for line in fp:
                    legend.append(line)
            params = ReadAndFit(f, legend)
            paramList.append({'parameters':params, 'name':str(x)})
        else:
            print('Missing form factors for ', x)
print("===================================================\n\n")
print(paramList)
