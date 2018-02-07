# random_variate.py 
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Algorithm 1 in Swisdak for Equation 8
#

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scp
import random_variate_headers as randvar
import analyse

A = 1.0
u = 0.7
gamma_u = 1./np.sqrt(1-u**2)
p_u = u*gamma_u
pm =  p_u/A*(1+np.sqrt(u**2+A**2))
p = np.linspace(0,12,1001)

def f(p):
    return (1+A*gamma_u*np.sqrt(1+p**2))*np.exp(-A*(p-p_u)**2/(np.sqrt(1+p**2)*gamma_u+p*p_u+1))

def f2m(pp):
    gamma_p = np.sqrt(1+pp**2)
    return np.sqrt(0.5/(A*gamma_u*gamma_p)**2*(1+np.sqrt(1+(2*A*gamma_u*gamma_p)**2)))

def Gamma(p):
    return np.sqrt(1+p**2)
#pminus = -p[np.where(f(p) == max(f(p)))[0][0]]

def foverfp(p):
    gamma_p = Gamma(p)
    f = p/(1+A*gamma_u*gamma_p)/(2*gamma_p)-(2*(gamma_p*gamma_u+p*p_u+1)*A*(p-p_u)-A*(p-p_u)**2*(gamma_u*p/gamma_p+p_u))/(gamma_p*gamma_u+p*p_u+1)**2
    return 1./f

def F(p):
    return f(p)/f(pm)-(1/np.e)**2

params = randvar.compute_params(F,foverfp,pm)

print(params)

X = 0.0
E = 0.0
N = 10000
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

    x[i] = randvar.generate_variate(f,params)

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
        return f2(ps)/f2(p2m) - (1/np.e)**2

    params_perp = randvar.compute_params(Ff,foverfp,p2m)

    params_perp[2] = -f2(params_perp[0])/f2_prime(params_perp[0])
    params_perp[3] = f2(params_perp[1])/f2_prime(params_perp[1])
    params_perp[4] = params_perp[3]/(params_perp[0]-params_perp[1])
    params_perp[5] = params_perp[2]/(params_perp[0]-params_perp[1]) 
    params_perp[6] = 1 - (params_perp[5]+params_perp[4])



    xs[i] = randvar.generate_variate(f2,params_perp)
    uperp = xs[i]*np.sqrt(1+x[i]**2)
    upll = x[i]
#    if(np.sqrt(((upll/Gamma(upll))-u)**2 + (uperp/Gamma(uperp))**2) >= 1.):
        # Speed of light violated. Reject
#        print(np.sqrt((upll/Gamma(upll))**2 + (uperp/Gamma(uperp))**2))
#        Reject += 1
#        continue
    i += 1

xperp = xs*np.sqrt(1+x**2)

# Plot generated contour
fig = plt.figure(1)
analyse.plot_contour(x,xperp,200,A,u)

# Plot generated histogram
fig2 = plt.figure(2)
N = analyse.plot_histogram(x,xperp,100)
analyse.plot_analytical_mag(N[1][0],N[1][-1],A,u,100,N)
plt.show()
