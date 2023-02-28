#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 10:18:37 2020

@author: cristina
"""


import numpy as np

def matrix(k,mu,J,S,t,Delta,alpha,N,N_imp,a,step,U,thetaS,phi):    
    H = np.zeros([N, N, 4, 4], dtype=complex)
    Hi = np.zeros([4*N, 4*N], dtype=complex)
    
    #diagonal terms
    for i_atom in range(N):        
        H[i_atom, i_atom, 0, 0] = - mu 
        H[i_atom, i_atom, 1, 1] = - mu 
        H[i_atom, i_atom, 2, 2] = + mu 
        H[i_atom, i_atom, 3, 3] = + mu 
        
        H[i_atom, i_atom, 0, 3] = -Delta
        H[i_atom, i_atom, 1, 2] = Delta
        H[i_atom, i_atom, 2, 1] = Delta
        H[i_atom, i_atom, 3, 0] = -Delta

    
    #off-diagonal terms
    for i_atom in range(N-1):
        #hopping
        H[i_atom, i_atom + 1, 0, 0] = -t
        H[i_atom, i_atom + 1, 1, 1] = -t
        H[i_atom, i_atom + 1, 2, 2] = +t
        H[i_atom, i_atom + 1, 3, 3] = +t
        
        H[i_atom + 1, i_atom, 0, 0] = -t
        H[i_atom + 1, i_atom, 1, 1] = -t
        H[i_atom + 1, i_atom, 2, 2] = +t
        H[i_atom + 1, i_atom, 3, 3] = +t
        
        #SOC
        H[i_atom, i_atom + 1, 0, 1] = +alpha
        H[i_atom, i_atom + 1, 1, 0] = -alpha
        H[i_atom, i_atom + 1, 2, 3] = -alpha
        H[i_atom, i_atom + 1, 3, 2] = +alpha
        
        H[i_atom + 1, i_atom, 0, 1] = -alpha
        H[i_atom + 1, i_atom, 1, 0] = +alpha
        H[i_atom + 1, i_atom, 2, 3] = +alpha
        H[i_atom + 1, i_atom, 3, 2] = -alpha
        
    
        
    ####cell edges        
    #hopping    
    H[0, N-1, 0, 0] = -t*np.exp(-1j*k*a)
    H[0, N-1, 1, 1] = -t*np.exp(-1j*k*a)
    H[0, N-1, 2, 2] = +t*np.exp(-1j*k*a)
    H[0, N-1, 3, 3] = +t*np.exp(-1j*k*a)
    
    H[N-1, 0, 0, 0] = -t*np.exp(1j*k*a)
    H[N-1, 0, 1, 1] = -t*np.exp(1j*k*a)
    H[N-1, 0, 2, 2] = +t*np.exp(1j*k*a)
    H[N-1, 0, 3, 3] = +t*np.exp(1j*k*a)
        
    #SOC
    H[0, N-1, 0, 1] = - alpha*np.exp(-1j*k*a)
    H[0, N-1, 1, 0] = + alpha*np.exp(-1j*k*a)
    H[0, N-1, 2, 3] = + alpha*np.exp(-1j*k*a)
    H[0, N-1, 3, 2] = - alpha*np.exp(-1j*k*a)
    
    H[N-1, 0, 0, 1] = + alpha*np.exp(1j*k*a)
    H[N-1, 0, 1, 0] = - alpha*np.exp(1j*k*a)
    H[N-1, 0, 2, 3] = - alpha*np.exp(1j*k*a)
    H[N-1, 0, 3, 2] = + alpha*np.exp(1j*k*a)
    
    ###add magnetic impurities
    for i in range(0, N_imp, step):
        i_imp = int((N - N_imp)/2) + i
        
        theta_i = thetaS[i]
        phi_i = phi[i]
        
        H[i_imp, i_imp, 0, 0] = H[i_imp, i_imp, 0, 0] + J*S*np.cos(theta_i) - U
        H[i_imp, i_imp, 1, 1] = H[i_imp, i_imp, 1, 1] - J*S*np.cos(theta_i) - U
        H[i_imp, i_imp, 2, 2] = H[i_imp, i_imp, 2, 2] - J*S*np.cos(theta_i) + U
        H[i_imp, i_imp, 3, 3] = H[i_imp, i_imp, 3, 3] + J*S*np.cos(theta_i) + U
        
        H[i_imp, i_imp, 0, 1] = H[i_imp, i_imp, 0, 1] + J*S*np.sin(theta_i)*np.exp(-1j*phi_i)  
        H[i_imp, i_imp, 1, 0] = H[i_imp, i_imp, 1, 0] + J*S*np.sin(theta_i)*np.exp(1j*phi_i)
        H[i_imp, i_imp, 2, 3] = H[i_imp, i_imp, 2, 3] - J*S*np.sin(theta_i)*np.exp(1j*phi_i)
        H[i_imp, i_imp, 3, 2] = H[i_imp, i_imp, 3, 2] - J*S*np.sin(theta_i)*np.exp(-1j*phi_i)
        
       
    #reshape matrix
    for i in range(N):
        for j in range(N):
            for t_i in range(4):
                for t_j in range(4):
                    Hi[(i) * 4 + t_i, (j) * 4 + t_j] = H[i, j, t_i, t_j]
                    
                    
    return(Hi)    