#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 09:48:27 2020

@author: cristina
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 11})

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return (array[idx], idx)

def wave_fun(E_total, psi_total, N, Nk, Nimp):
    
    nk = int(Nk/2.0)#check in the middle of k-vector
    #nk = 0# check im k=-pi/a
    min_E = min(abs(E_total[:,nk]))#min energy 
    
    E_s = sorted(abs(E_total[:,nk]))
    min_E_exc = E_s[3]
    print('min', min_E, 'exited', min_E_exc)
    
    #select positive and negative min energy
    nE = np.where(abs(E_total[:,nk]) == min_E)
    (min_E2, nE2) = find_nearest(E_total[:,nk], -min_E)
    psi_k0 = psi_total[:,:,nk]
    psi_E = psi_k0[:,nE[0]]#eigen vector corresponding to min energy
    psi_E2 = psi_k0[:,nE2]#eigen vector corresponding to min energy negative
    print('min_E index:', nE)
    print('min_E2_index:', nE2)
    
    #convert to Majorana states
    '''Majorana basis change'''
    U1 = 1/(np.sqrt(2))*np.array([[1,0,1,0],[0,1,0,1],[1j,0,-1j,0],[0,1j,0,-1j]])
    U_inv1 = 1/(np.sqrt(2))*np.array([[1,0,-1j,0],[0,1,0,-1j],[1,0,1j,0],[0,1,0,1j]])
    unity = np.eye(N)

    U_u = np.kron(unity, U1)
    #U_inv = np.kron(unity, U_inv1)
    
    psi_EM = np.dot(U_u,psi_E)
    psi_E2M = np.dot(U_u,psi_E2)
    
    #select u, v
    u_up = np.zeros(N, dtype = complex)
    v_up = np.zeros(N, dtype = complex)
    u_down = np.zeros(N, dtype = complex)
    v_down = np.zeros(N, dtype = complex)
    
    u_up2 = np.zeros(N, dtype = complex)
    v_up2 = np.zeros(N, dtype = complex)
    u_down2 = np.zeros(N, dtype = complex)
    v_down2 = np.zeros(N, dtype = complex)
    

    for i in range(N):

        u_up[i] = psi_EM[0::4][i]  
        v_up[i] = psi_EM[2::4][i]        
        u_down[i] = psi_EM[1::4][i]
        v_down[i] = psi_EM[3::4][i] 
        
        u_up2[i] = psi_E2M[0::4][i]  
        v_up2[i] = psi_E2M[2::4][i] 
        u_down2[i] = psi_E2M[1::4][i]  
        v_down2[i] = psi_E2M[3::4][i]
        
        
        
    #save data
    #np.savetxt('u_up.txt', u_up)
    #np.savetxt('u_down.txt', u_down)
    #np.savetxt('v_up.txt', v_up)
    #np.savetxt('v_dn.txt', v_down)
    
    #np.savetxt('u2_up.txt', u_up2)
    #np.savetxt('u2_down.txt', u_down2)
    #np.savetxt('v2_up.txt', v_up2)
    #np.savetxt('v2_dn.txt', v_down2)        
        

    return(u_up, v_up, u_down, v_down, u_up2, v_up2, u_down2, v_down2, min_E, min_E2)




