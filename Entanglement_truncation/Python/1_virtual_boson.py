import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigs, LinearOperator, gmres
from scipy.optimize import minimize

# cMPS code for the truncated gaussian cTNS
# obtained by a harmonic virtual hilbert space

def kron_delta(m,n):
    if m==n:
        return 1
    else:
        return 0

##########################
# CONSTRUCT CMPS MATRICES
##########################
e_exact = -0.0865028
alpha = 0.62219402
# V should be bigger than alpha**2
V = 0.80029003
Omega = np.sqrt(V**2-alpha**4)
energies = []
chis = range(2,10,1)
# bond dimension
for chi in chis: 
    
    x_exp, xsq_exp, psq_exp = np.empty((chi,chi)), np.empty((chi,chi)), np.empty((chi,chi))
    for m in range(chi):
        for n in range(chi):
            x_exp[m][n] = np.sqrt(n)*kron_delta(n,m+1)+np.sqrt(n+1)*kron_delta(n,m-1)
            xsq_exp[m][n] = (2*n+1)*kron_delta(n,m) + np.sqrt((n+2)*(n+1))*kron_delta(n, m-2) + np.sqrt(n*(n-1))*kron_delta(n,m+2)
            psq_exp[m][n] = (2*n+1)*kron_delta(n,m) - np.sqrt((n+2)*(n+1))*kron_delta(n,m-2) - np.sqrt(n*(n-1))*kron_delta(n,m+2)
    x_exp /= np.sqrt(2*Omega)
    xsq_exp /= (2*Omega)
    psq_exp *= Omega/2
    R = alpha*x_exp
    Q = - (psq_exp/2 + V*xsq_exp )
    
    ########################
    # CALCULATE FIXED POINTS
    ########################
    def FindSSL(Q, R):
    
        chi = Q.shape[0]
        
        # construct transfer matrix handle and cast to LinearOperator
        def transferLeftHandle(rho):
            rho = rho.reshape(chi,chi)
            return np.reshape(rho@Q+np.conj(Q).T@rho+np.conj(R).T@rho@R, chi**2)
        transferLeft = LinearOperator((chi ** 2, chi ** 2), matvec=transferLeftHandle)
    
        # calculate fixed point
        lam, l = eigs(transferLeft, k=1, which='LM')
    
        return l.reshape(chi, chi)
    
    def FindSSR(Q, R):
    
        chi = Q.shape[0]
    
        # construct transfer matrix handle and cast to LinearOperator
        def transferRightHandle(rho):
            rho = rho.reshape(chi,chi)
            return np.reshape(rho@np.conj(Q).T+Q@rho+R@rho@np.conj(R).T, chi**2)
        transferRight = LinearOperator((chi ** 2, chi ** 2), matvec=transferRightHandle)
    
        # calculate fixed point
        lam, r = eigs(transferRight, k=1, which='LM')
    
        return r.reshape(chi, chi)
    #################
    # OPTIMIZATION
    #################
    mu, nu = 1, 0.4
    c = 0
    
    rhoL = FindSSL(Q, R)
    rhoR = FindSSL(Q, R)
    # make sure density matrices are hermitian
    rhoL = (rhoL+np.conj(rhoL).T)/2
    rhoR = (rhoR+np.conj(rhoR).T)/2
    #normalize
    rhoR = rhoR/np.trace(rhoL@rhoR)
    QR = Q@R - R@Q
    ekin = lambda Q, R: np.trace(rhoL@QR@rhoR@np.conj(QR).T)
    density = lambda Q, R: mu * np.trace(rhoL@R@rhoR@np.conj(R).T)
    pair = lambda Q, R: - nu * np.trace(rhoL@R@rhoR@R + rhoL@np.conj(R).T@rhoR@np.conj(R).T)
    e = lambda Q,R: ekin(Q,R) + density(Q,R) + pair(Q,R)
    #eint = c * np.trace(rhoL@R**2@rhoR@np.conj(R).T**2)
    energies.append(e(Q,R)-e_exact)
    print(e(Q,R))
plt.plot(chis, energies)
plt.xlabel('bond dimension')
plt.ylabel('energy error')
print(e(Q,R))
