# random_variate_headers.py
# 
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Function declarations should be here


import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scp

def compute_params(F,foverfp,pm):
    pcand = pm
    i = 1.0

    while(F(pcand) > 0):
        i *= 2.0
        pcand = i*pm

    pplus = scp.bisect(F,pm,pcand)

    pcand = pm
    i = -1.0

    while(F(pcand) > 0):
        i *= 2.0
        pcand = i*pm

    pminus = scp.bisect(F,pm,pcand)

    #P = np.linspace(pminus,pplus,101)
    #F = f(P)
    #F /= sum(F)*(pplus-pminus)/101.

    lambdap = -foverfp(pplus)
    lambdam = foverfp(pminus)

    qmi = lambdam/(pplus-pminus)
    qp = lambdap/(pplus-pminus)
    qm = 1 - (qp+qmi)

    return np.array([pplus,pminus,lambdap,lambdam,qmi,qp,qm,pm])

def generate_variate(f,params):
    pplus = params[0]
    pminus = params[1]
    lambdap = params[2]
    lambdam = params[3]
    qmi = params[4]
    qp = params[5]
    qm = params[6]
    pm = params[7]
    flag = 0
    while(flag == 0):
	# Create random variates U and V
        U = np.random.uniform(low=0.0,high=1.0)
        V = np.random.uniform(low=0.0,high=1.0)
        if(U <= qm):
            Y = U/qm
            X = (1-Y)*(pminus+lambdam)+Y*(pplus-lambdap)
            if(V <= (f(X)/f(pm))):
                flag = 1
        elif(U <= qm+qp):
            E = -np.log(abs((U-qm)/qp))
            X = pplus - lambdap*(1-E)
            if(V<=np.exp(E)*f(X)/f(pm)):
                flag = 1
        else:
            E = -np.log(abs((U-(qm+qp))/qmi))
            X = pminus + lambdam*(1-E)
            if(V <= np.exp(E)*f(X)/f(pm)):
                flag = 1
    return X
