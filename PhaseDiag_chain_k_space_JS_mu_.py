# -*- coding: utf-8 -*-
"""
Created on Fri May  8 12:03:32 2020

@author: nodar
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 11})
#from numpy import linalg as LA
#from numpy.linalg import inv
#diag = LA.eig


'''Parameters'''
import param_phase_diag as pa
a = pa.a
N = pa.N
k = np.linspace(-np.pi, np.pi, N)/a#k vector

t = pa.t
Delta = pa.Delta#
K = pa.K
alpha = pa.alpha
theta = pa.theta
phi = pa.phi

#####2D parameter space
N_1 = pa.N_1
N_2 = pa.N_2

mu_1 = pa.mu_1
mu_2 = pa.mu_2
mu = np.linspace(mu_1, mu_2, N_1)

j1= pa.j1
j2 = pa.j2
J = np.linspace(j1, j2, N_2)


U = 1/(np.sqrt(2))*np.array([[1,0,1,0],[0,1,0,1],[1j,0,-1j,0],[0,1j,0,-1j]])
U_inv = 1/(np.sqrt(2))*np.array([[1,0,-1j,0],[0,1,0,-1j],[1,0,1j,0],[0,1,0,1j]])



import topo_inv as topo     

Q_2D = np.zeros([N_1, N_2], dtype=float)
W = np.zeros([N_1, N_2], dtype=float)
min_gap = np.zeros([N_1, N_2], dtype=float)

for N_i in range(N_1):
    for N_j in range(N_2):
        
        Q_2D[N_i, N_j] = topo.topological_inv(U, U_inv, mu[N_i], J[N_j], t,Delta, alpha, a, K, theta, phi)       
        W[N_i, N_j] = topo.winding_number(k, U, U_inv, mu[N_i], J[N_j], t,Delta, alpha, a, K, theta, phi)
        min_gap[N_i, N_j] = topo.diagonalization(k,mu[N_i], J[N_j], t,Delta,alpha,K, theta, phi, a)        
        
        #print(N_i, N_j)
        
#%%

'''Plot 2D phase diagram''' 
m = Q_2D.T 
array_plot = np.flip(m, axis=0)
np.savetxt('2D_topo.txt', m)
      
plt.figure(1)
plt.imshow(m, origin='lower', aspect='auto', extent=[mu_1, mu_2, j1, j2])
plt.xlabel(r'$\mu$')
plt.ylabel('J')
plt.title(r'Topological invariant $Q$, $\Delta=0.5$, t=1.0, $\alpha=3.0$, U=0.0')
plt.savefig('Q_2D.png', dpi = 260, bbox_inches='tight')

gap = min_gap.T 
np.savetxt('min_gap.txt', gap)
      
plt.figure(2)
plt.imshow(gap, origin='lower', aspect='auto', extent=[mu_1, mu_2, j1, j2], cmap = 'bone')
plt.xlabel(r'$\mu$')
plt.ylabel('J')
plt.title(r'Energy gap, $\Delta=0.5$, t=1.0, $\alpha=3.0$, U=0.0')
plt.savefig('gap_min.png', dpi = 260, bbox_inches='tight')

mm = W.T 
np.savetxt('2D_W.txt', mm)
      
plt.figure(3)
plt.imshow(mm, vmin=-1, vmax=1, origin='lower', aspect='auto', extent=[mu_1, mu_2, j1, j2], cmap = 'PiYG')
plt.xlabel(r'$\mu$')
plt.ylabel('J')
plt.title(r'Winding number $W$, $\Delta=0.5$, t=1.0, $\alpha=3.0$, U=0.0')
plt.colorbar(orientation='vertical', fraction=.1)
plt.savefig('W_2D.png', dpi = 260, bbox_inches='tight')

#save parameters
np.savetxt('mu_vector.txt', mu)
np.savetxt('J_vector.txt', J)

