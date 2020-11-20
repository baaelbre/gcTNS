# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 19:43:44 2020

@author: Gebruiker
"""
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def E_integrated(mu,nu):
    def energy_density_p(p): return 1/(2*np.pi)*(excitation(p, mu, nu) - (p**2 + mu))
    return quad(energy_density_p,0,np.inf,epsabs=10**-25,limit=10000)[0]
def C_00_exact(p,mu,nu):
    return 0.5*(-1 + 1/np.sqrt(1-(4 *nu**2)/(p**2 + mu)**2))
def T_exact(p,mu,nu):
    return 0.5*(-p**2 + p**2/np.sqrt(1-(4 *nu**2)/(p**2 + mu)**2))
def C_01_exact(p,mu,nu):
    return nu/((p**2+mu)*np.sqrt(1-(4 *nu**2)/(p**2 + mu)**2))
def C_10_exact(p,mu,nu):
    #zelfde uitdrukking
    return C_01_exact(p,mu,nu)
def excitation(p, mu, nu):
    return np.sqrt((p**2+mu)**2-4*nu**2)

if __name__ == '__main__':
    mu, nu = 1, 0.4
    print(E_integrated(mu,nu))
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlabel(r'$\mu$', fontsize=15)
    ax.set_ylabel(r'$\nu$', fontsize=15)
    ax.set_zlabel(r'$\epsilon(0, \mu, \nu)$', fontsize=15)
    mu_, nu_ = np.linspace(0,2,100), np.linspace(-1,1,100)
    mu_, nu_ = np.meshgrid(mu_, nu_)
    ax.plot_surface(mu_, nu_, excitation(0, mu_, nu_))
    fig_exc = plt.figure(figsize=(8,6))
    ax_exc = fig_exc.add_subplot(111)
    ax_exc.grid()
    ax_exc.set_xlabel(r'$p$', fontsize=15)
    ax_exc.set_ylabel(r'$\epsilon(p)$', fontsize=15)
    p = np.linspace(-2,2,100)
    ax_exc.plot(p, excitation(p,mu,nu), label = r'1p')
    #ax_exc.plot(p, 2*excitation(p/2,mu,nu), label = r'2p')
    #ax_exc.plot(p, 3*excitation(p/3,mu,nu), label = r'3p')
    ax_exc.legend(fontsize=15)
