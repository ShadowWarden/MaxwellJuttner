# random_variate.py 
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Algorithm 1 in Swisdak for Equation 8
#

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scp

A = 1.0
u = 0.1
gamma_u = 1./np.sqrt(1-u**2)
p_u = u*gamma_u
pm =  p_u/A*(1+np.sqrt(u**2+A**2))
p = np.linspace(0,12,1001)

def f(p):
    return (1+A*gamma_u*np.sqrt(1+p**2))*np.exp(-A*(p-p_u)**2/abs(np.sqrt(1+p**2)*gamma_u+p*p_u+1))

def f_prime(p):
    gamma_p = np.sqrt(1+p**2)
    return f(p)*(A*gamma_u/(2*gamma_p*(1+A*gamma_u*gamma_p))-A*(0.5*gamma_u/gamma_p-gamma_u*u))

def f2m(pp):
    gamma_p = np.sqrt(1+pp**2)
    return np.sqrt(0.5/(A*gamma_u*gamma_p)**2*(1+np.sqrt(1+(2*A*gamma_u*gamma_p)**2)))

def Gamma(p):
    return np.sqrt(1+p**2)
#pminus = -p[np.where(f(p) == max(f(p)))[0][0]]

def F(p):
    return f(p)/f(pm)-(1/np.e)

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

lambdap = -f(pplus)/f_prime(pplus)
lambdam = f(pminus)/f_prime(pminus)

qmi = lambdam/(pplus-pminus)
qp = lambdap/(pplus-pminus)
qm = 1 - (qp+qmi)

print(qp,qm,qmi)
print(f(pminus)-f(pm)/np.e)
print(lambdap,lambdam)

X = 0.0
E = 0.0
N = 1000
Flag = 0
x = np.zeros(N)
xs = np.zeros(N)
print("N = ",N)
pplus_perp_start = 1.0
Reject = 0

i = 0
while(i<N):
    flag = 0
    print("Iteration:",i,"Rejects:",Reject)
    while(flag == 0):
	# Create random variates U and V
        U = np.random.uniform(low=0.0,high=1.0)
        V = np.random.uniform(low=0.0,high=1.0)
        if(U <= qm):
            Y = U/qm
            X = (1-Y)*(pminus+lambdam)+Y*(pplus-lambdap)
            if(V <= f(X)/f(pm)):
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
#	print("Iteration ",i)
#	print(X)
    x[i] = X

    pp = x[i]
    p2m = f2m(pp)
    # Compute for pperp
    def f2(ps):
        gamma_p = np.sqrt(1+pp**2)
        gamma_s = np.sqrt(1+ps**2)
        return ps*np.exp(-A*((pp-p_u)**2+gamma_p**2*gamma_u**2*ps**2)/(gamma_p*gamma_u*gamma_s+pp*p_u+1))
    
    def f2_prime(ps):
        gamma_p = np.sqrt(1+pp**2)
        gamma_s = np.sqrt(1+ps**2)
        V = gamma_p*gamma_u*gamma_s+pp*p_u+1
        U = gamma_p**2*gamma_u**2*ps**2+(pp-p_u)**2
        fpr = 2*ps**2*(V*gamma_p**2*gamma_u**2-0.5*U/gamma_s)/V**2
        return (1-fpr)*np.exp(-A*((pp-p_u)**2+gamma_p**2*gamma_u**2*ps**2)/(gamma_p*gamma_u*gamma_s+pp*p_u+1))

    def Ff(ps):
        return f2(ps)/f2(p2m) - (1/np.e)

    pcand = p2m
    j = 1.0

    while(Ff(pcand) > 0):
        j *= 2.0
        pcand = j*p2m

    pplus_perp = scp.bisect(Ff,p2m,pcand)

    pcand = p2m
    j = -1.0

    while(Ff(pcand) > 0):
        j *= 2.0
        pcand = j*p2m

    pminus_perp = scp.bisect(Ff,p2m,pcand)

    lambdap_perp = -f2(pplus_perp)/f2_prime(pplus_perp)
    lambdam_perp = f2(pminus_perp)/f2_prime(pminus_perp)
    qmi_perp = lambdam_perp/(pplus_perp-pminus_perp)
    qp_perp = lambdap_perp/(pplus_perp-pminus_perp)
    qm_perp = 1 - (qp_perp+qmi_perp)

    flag = 0
    # Create random variates U and V
    U = np.random.uniform(low=0.0,high=1.0)
    V = np.random.uniform(low=0.0,high=1.0)
    if(U <= qm_perp):
        Y = U/qm_perp
        X = (1-Y)*(pminus_perp+lambdam_perp)+Y*(pplus_perp-lambdap_perp)
        if(V <= f2(X)/f2(p2m)):
            flag = 1
    elif(U <= qm_perp+qp_perp):
        E = -np.log(abs((U-qm_perp)/qp_perp))
        X = pplus_perp - lambdap_perp*(1-E)
        if(V<=np.exp(E)*f2(X)/f2(p2m)):
            flag = 1
    else:
        E = -np.log(abs((U-(qm_perp+qp_perp))/qmi_perp))
        X = pminus_perp + lambdam_perp*(1-E)
        if(V <= np.exp(E)*f2(X)/f2(p2m)):
            flag = 1
#	print("Iteration ",i)
#	print(X)
    if(flag == 0):
        continue
    xs[i] = X
    uperp = xs[i]*np.sqrt(1+x[i]**2)
    upll = x[i]
#    if(np.sqrt(((upll/Gamma(upll))-u)**2 + (uperp/Gamma(uperp))**2) >= 1.):
        # Speed of light violated. Reject
#        print(np.sqrt((upll/Gamma(upll))**2 + (uperp/Gamma(uperp))**2))
#        Reject += 1
#        continue
    i += 1

xperp = xs*np.sqrt(1+x**2)     
