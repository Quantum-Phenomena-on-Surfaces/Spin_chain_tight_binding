# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 15:23:22 2020

@author: nodar
"""

import numpy as np
#import matplotlib as mpl
#mpl.use('Agg') ####to run in Oberon/other clusters
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.figsize'] = [8, 6]
from numpy import linalg as LA
diag = LA.eig
from tqdm import trange

'''Choose the basis to solve'''
basis = 'eh'#Maj or eh

'''Parameters'''
import param as pa
t = pa.t
mu = pa.mu
Delta = pa.Delta
J = pa.J
S = pa.S
B = J*S#Zeeman term
alpha = pa.alpha
U = pa.U
state = pa.state

'''size and k-vector'''
N_atoms = pa.N_atoms
a = pa.a
a = N_atoms*a#cell size
N_imp = pa.N_imp
step = pa.step
Nk = pa.Nk
k = np.linspace(-np.pi, np.pi, Nk)/a#k vector

'''impurities'''
if (state == 'FM'):
        thetaS = np.zeros(N_imp, dtype = 'float')
elif (state == 'AF'):
    thetaS = np.zeros(N_imp, dtype = 'float')
    for i in range(int(len(thetaS)/2)):
        thetaS[2*i+1] = np.pi
        if N_imp % 2 == 0:
            thetaS[-1] = np.pi
elif (state == 'inplane'):
    thetaS = np.full(N_imp, - np.pi/2.0, dtype = 'float')
        
phi = np.zeros(N_imp)

E_total = np.zeros([4*N_atoms, Nk], dtype = complex)
psi_total = np.zeros([4*N_atoms,4*N_atoms,Nk], dtype = complex)

'''Majorana basis change'''
U1 = 1/(np.sqrt(2))*np.array([[1,0,1,0],[0,1,0,1],[1j,0,-1j,0],[0,1j,0,-1j]])
U_inv1 = 1/(np.sqrt(2))*np.array([[1,0,-1j,0],[0,1,0,-1j],[1,0,1j,0],[0,1,0,1j]])
unity = np.eye(N_atoms)

U_u = np.kron(unity, U1)
U_inv = np.kron(unity, U_inv1)
    

'''Diagonalize H'''
import BdG_general_matrix as mt  
for N_i in trange(Nk):
    
    H_i = mt.matrix(k[N_i],mu,J,S,t,Delta,alpha,N_atoms,N_imp,a,step,U,thetaS,phi)
    
    if (basis == 'eh'):
        H_i2 = H_i
        
    elif (basis == 'Maj'):  
        H_i2 = U_u@H_i@U_inv
        
    else:
        print('basis= Maj or eh')
                            
    #Diagonalize                       
    H_m = np.matrix(H_i2)    
    (E, psi) = diag(H_m)
    
    #E_order = np.sort(E)
    E_total[:, N_i] = E
    psi_total[:,:,N_i] = psi
    
#%%
   
'''u, v components E=0'''
import BdG_wavefun as wf
(u_up, v_up, u_down, v_down, u_up2, v_up2, u_down2, v_down2,\
           min_E, min_E2) = wf.wave_fun(E_total, psi_total, N_atoms, Nk, N_imp)   

print('min value 1', np.real(min_E))
print('min value 2', np.real(min_E2))


'''Plot all E_n versus k'''
plt.figure(1)
for n in range(4*N_atoms):    
    plt.plot(k, np.real(E_total[n,:]), 'r.', ms = 1.0)
    
plt.xlabel(r'$k$')
plt.ylabel(r'$E$')
plt.xlim([-np.pi/a,np.pi/a])        
plt.ylim([-1.0, 1.0])  
plt.savefig('bands.png', dpi=300, bbox_inches='tight')
#np.savetxt('Energies.txt', E_total)
     
        
'''LDOS calculation in first and last atom'''
import LDOS_calc_new as cal
omega = np.linspace(-1.5*Delta, 1.5*Delta, 1001)#omega vector  r  
D = pa.Dynes


(LDOS_M_u_up, LDOS_M_u_down, LDOS_M_v_up, LDOS_M_v_down, LDOS_N_u_up, LDOS_N_u_down, LDOS_N_v_up, LDOS_N_v_down, M, N,\
 LDOS_tot_uup, LDOS_tot_udn, LDOS_tot_vup, LDOS_tot_vdn) = cal.LDOS_calc(E_total, psi_total, D, omega, N_atoms, N_imp, Nk)

#%%

'''plot u-v'''
import plot_uv_chain as uv
uv.plot_uv(u_up, v_up, u_down, v_down, u_up2, v_up2, u_down2, v_down2,\
           min_E, min_E2, M, N, N_atoms)


'''plot LDOSat both ends'''
import plot_LDOS as LD
LD.LDOS(LDOS_M_u_up, LDOS_M_u_down, LDOS_M_v_up, LDOS_M_v_down, LDOS_N_u_up, LDOS_N_u_down, LDOS_N_v_up, LDOS_N_v_down, M, N, omega)






























     