# analyse.py
#
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Plotting and the like

import matplotlib.pyplot as plt
import numpy as np

def plot_analytical_mag(plower,pupper,A,u,Npts,N):
    gamma_u = 1/np.sqrt(1-u**2)
    p_u = u*gamma_u
    P = np.linspace(plower,pupper,Npts)
    F = P*np.sinh(A*p_u*P)/(A*p_u)*np.exp(-A*(gamma_u*np.sqrt(1+P**2)-1))
    F *= np.sum(N[0])/np.sum(F)
    plt.plot(P,F,'b-.',label="Analytical Prediction")
    plt.xlabel("$|p|$")
    plt.legend()
    plt.draw()

def plot_histogram(x,xperp,Nbins):
    X = (x**2+xperp**2)**0.5
    N = plt.hist(X,Nbins,normed=True,color='red',label="Generated Histogram")
    plt.draw()
    return N

def plot_contour(x,xperp,Nbins,A,u):
    xedges = np.linspace(-max(x),max(x),Nbins)
    yedges = np.linspace(-max(xperp),max(xperp),Nbins)
    H = np.histogram2d(xperp,x,bins=(xedges,yedges),normed=True)
    X,Y = np.meshgrid(xedges,yedges)
    plt.pcolormesh(X,Y,np.log(H[0]),cmap='RdBu')
    plt.title(r"Plot of $\log(f(x_{\parallel},x_{\perp}))$ for $A= %.2f$ and $u = %.2f$" % (A,u))
    plt.xlabel(r"$x_{\parallel}$")
    plt.ylabel(r"$x_{\perp}$")
    plt.grid()
    plt.colorbar()
    plt.draw()
