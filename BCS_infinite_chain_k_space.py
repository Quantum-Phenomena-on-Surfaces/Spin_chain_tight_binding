# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 12:18:44 2020

@author: nodar
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 11})
from numpy import linalg as LA
diag = LA.eig
from numpy.linalg import inv 

'''Parameters'''
import param_inf as pa
N = pa.N
a = pa.a
k = np.linspace(-np.pi, np.pi, N)/a#k
mu = pa.mu
t = pa.t
Delta = pa.Delta
B = pa.JS
alpha = pa.alpha

'''Matrix'''
def matrix(k,mu,B,t,Delta,alpha):
    H = np.zeros([4, 4], dtype=complex)
    H[0, 0] = - mu + B - t*np.exp(1j*k*a) - t*np.exp(-1j*k*a)
    H[1, 1] = - mu - B - t*np.exp(1j*k*a) - t*np.exp(-1j*k*a)
    H[2, 2] = + mu - B + t*np.exp(1j*k*a) + t*np.exp(-1j*k*a)
    H[3, 3] = + mu + B + t*np.exp(1j*k*a) + t*np.exp(-1j*k*a)
        
    H[0, 3] = -Delta
    H[1, 2] = Delta
    H[2, 1] = Delta
    H[3, 0] = -Delta         
        
    H[0, 1] =  alpha*np.exp(1j*k*a) - alpha*np.exp(-1j*k*a)
    H[1, 0] = -alpha*np.exp(1j*k*a) + alpha*np.exp(-1j*k*a)
    H[2, 3] = -alpha*np.exp(1j*k*a) + alpha*np.exp(-1j*k*a)
    H[3, 2] =  alpha*np.exp(1j*k*a) - alpha*np.exp(-1j*k*a)   
    
    return(H)
    
def matrix_diff(k,mu,B,t,Delta,alpha):
    H = np.zeros([4, 4], dtype=complex)
    H[0, 0] = -1j*a*t*np.exp(1j*k*a) + 1j*a*t*np.exp(-1j*k*a)
    H[1, 1] = -1j*a*t*np.exp(1j*k*a) + 1j*a*t*np.exp(-1j*k*a)
    H[2, 2] = +1j*a*t*np.exp(1j*k*a) - 1j*a*t*np.exp(-1j*k*a)
    H[3, 3] = +1j*a*t*np.exp(1j*k*a) - 1j*a*t*np.exp(-1j*k*a)
        
    H[0, 3] = 0.0
    H[1, 2] = 0.0
    H[2, 1] = 0.0
    H[3, 0] = 0.0        
        
    H[0, 1] =  alpha*1j*a*np.exp(1j*k*a) + alpha*1j*a*np.exp(1j*k*a)
    H[1, 0] = -alpha*1j*a*np.exp(1j*k*a) - alpha*1j*a*np.exp(1j*k*a)
    H[2, 3] = -alpha*1j*a*np.exp(1j*k*a) - alpha*1j*a*np.exp(1j*k*a)
    H[3, 2] =  alpha*1j*a*np.exp(1j*k*a) + alpha*1j*a*np.exp(1j*k*a)
    
    return(H)
    

E_total = np.zeros([4, N], dtype=complex)
E_total_sorted = np.zeros([4, N], dtype=complex)
psi_total = np.zeros([4,4,N], dtype=complex)

band1 = np.zeros(N, dtype = complex)
band2 = np.zeros(N, dtype = complex)
band3 = np.zeros(N, dtype = complex)
band4 = np.zeros(N, dtype = complex)

'''Diagonalize for each k'''

Hi_k = np.zeros([4, 4, N], dtype=complex)
Hi_k_diff = np.zeros([4, 4, N], dtype=complex)

for N_i in range(N):                   
    
    H_i = matrix(k[N_i],mu,B,t,Delta,alpha)
    Hi_k[:,:,N_i] =  H_i
    Hi_k_diff[:,:,N_i] = matrix_diff(k[N_i],mu,B,t,Delta,alpha)
                            
    #Diagonalize                       
    H_m = np.matrix(H_i)    
    (E, psi) = diag(H_m)
    E_total[:, N_i] = E
    
    E_order = np.sort(E)##ordenar
    E_order = E##no ordenar

   
    E_total_sorted[:, N_i] = E_order
    psi_total[:,:,N_i] = psi

'''Topological invariant calculation'''
from pfapack import pfaffian as pf

#unitary transformation to Majorana basis
U = 1/(np.sqrt(2))*np.array([[1,0,1,0],[0,1,0,1],[1j,0,-1j,0],[0,1j,0,-1j]])
U_inv = 1/(np.sqrt(2))*np.array([[1,0,-1j,0],[0,1,0,-1j],[1,0,1j,0],[0,1,0,1j]])


#unitary rotation
r = np.sin(np.pi/4)
U_pi = np.array([[r,0,r,0],[0,r,0,r],[-r,0,r,0],[0,-r,0,r]])
U_pi_inv = inv(U_pi)

H_00 = matrix(0, mu, B, t, Delta, alpha)
H_0 = 1j*U@H_00@U_inv
pfa1 = pf.pfaffian(H_0)
    
H_pii = matrix(np.pi/a, mu, B, t, Delta, alpha)
H_pi = 1j*U@H_pii@U_inv
pfa2 = pf.pfaffian(H_pi)

#topological invariant    
Q = np.sign(pfa1*pfa2)
print('Topological invariant Q=', int(np.real(Q)))

'''Plot all En versus k'''
plt.figure(1)
for n in range(4):    
    plt.plot(k, np.real(E_total_sorted[n,:]), '.C0', ms = 1.0)

plt.xlabel(r'$k$')
plt.ylabel(r'$E$')
#plt.title(r'$\mu=2.0$, $t=1.0$, $\Delta = 1.0$, $JS = 0.0$, $\alpha = 1.0$')
plt.xlim([-np.pi/a, np.pi/a])
plt.show()
#plt.axvline(np.pi/3*a, color='k', linestyle='--', linewidth= 0.7)
#plt.axvline(-np.pi/3*a, color='k', linestyle='--', linewidth= 0.7)
plt.savefig('bands.png')

#####save data
np.savetxt('bands.txt', np.real(E_total_sorted))
np.savetxt('k_vector.txt', k)

'''Winding number'''
#tau_z = np.array([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]])

H_ii = np.zeros(N - 1, dtype=complex)

for N_i in range(1, N):         
        
    H_i = Hi_k[:,:,N_i]
    H_i_1 = Hi_k[:,:,N_i-1]
    
    A_i =U_pi@H_i@U_pi_inv
    A_i_1 =U_pi@H_i_1@U_pi_inv
    
    dA = (A_i - A_i_1)/(k[N_i] - k[N_i-1])    
    
    t = dA@inv(A_i)
    H_ii[N_i - 1] = t[0,0] + t[1,1]
        
w = np.trapz(H_ii, k[1 : N])
W =1.0/(2*np.pi*1j)*w
#print(w)
print('Winding', W)
#print('Winding module', abs(W))


































    
