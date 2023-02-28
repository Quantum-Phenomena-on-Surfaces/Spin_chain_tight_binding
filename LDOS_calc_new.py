#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 11:25:58 2021

@author: cristina
"""

import numpy as np
#import LDOS_func as f
pi = np.pi

def LDOS_calc(E_total, psi_total, D, omega, N_atoms, N_imp, Nk):
    
    #select k=0
    nk = int(Nk/2.0)#check in the middle of k-vector
    E = E_total[:,nk]#min energy    
    index_plus = np.where(E>0)#positive energy values
    #E_plus = E[index_plus]    
    
    M = int((N_atoms - N_imp)/2)#first atom index
    N = M+int(N_imp/2.0)#middle of the chain
    psi = psi_total[:,:,nk]    
    
    #first atom
    u_Mup = np.zeros(int(len(E)), dtype=complex)
    u_Mdown = np.zeros(int(len(E)), dtype=complex)
    v_Mup = np.zeros(int(len(E)), dtype=complex)
    v_Mdown = np.zeros(int(len(E)),dtype=complex)
    #last atom
    u_Nup = np.zeros(int(len(E)), dtype=complex)
    u_Ndown = np.zeros(int(len(E)), dtype=complex)
    v_Nup = np.zeros(int(len(E)), dtype=complex)
    v_Ndown = np.zeros(int(len(E)),dtype=complex)
    #first atom
    u_Mup[::] = psi[4*M, :]
    u_Mdown[::] = psi[4*M+1, :]
    v_Mup[::] = psi[4*M+2, :]
    v_Mdown[::] = psi[4*M+3, :]
    #last atom
    u_Nup[::] = psi[4*N, :]
    u_Ndown[::] = psi[4*N+1, :]
    v_Nup[::] = psi[4*N+2, :]
    v_Ndown[::] = psi[4*N+3, :]

    #first LDOS
    LDOS_M_u_up = np.zeros(len(omega))
    LDOS_M_u_down = np.zeros(len(omega))
    LDOS_M_v_up = np.zeros(len(omega))
    LDOS_M_v_down = np.zeros(len(omega))
    
    #last LDOS
    LDOS_N_u_up = np.zeros(len(omega))
    LDOS_N_u_down = np.zeros(len(omega))
    LDOS_N_v_up = np.zeros(len(omega))
    LDOS_N_v_down = np.zeros(len(omega))
    
    def dirac(x, xi):
        t = xi/(pi*(x**2 + xi**2))
        return(t)
    
    
    def LDOS_u(omega, E, u, xi):
        t = sum(u*np.conjugate(u)*dirac(omega - E, xi))
        return(t)
    
    
    def LDOS_v(omega, E, v, xi):
        t = sum(v*np.conjugate(v)*dirac(omega + E, xi))
        return(t)    
    
    for i in range(len(omega)):

        LDOS_M_u_up[i] = LDOS_u(omega[i], E,u_Mup, D)    
        LDOS_M_u_down[i] = LDOS_u(omega[i], E, u_Mdown, D)
        LDOS_M_v_up[i] = LDOS_v(omega[i], E,v_Mup, D)    
        LDOS_M_v_down[i] = LDOS_v(omega[i], E, v_Mdown, D)
        
        LDOS_N_u_up[i] = LDOS_u(omega[i], E, np.real(u_Nup), D)    
        LDOS_N_u_down[i] = LDOS_u(omega[i], E, np.real(u_Ndown), D)
        LDOS_N_v_up[i] = LDOS_v(omega[i], E, np.real(v_Nup), D)    
        LDOS_N_v_down[i] = LDOS_v(omega[i], E, np.real(v_Ndown), D)
        
        
        
    ######LDOS everywhere in the chain
    LDOS_tot_uup = np.zeros([N_atoms, len(omega)], dtype = complex)
    LDOS_tot_udn = np.zeros([N_atoms, len(omega)], dtype = complex)
    LDOS_tot_vup = np.zeros([N_atoms, len(omega)], dtype = complex)
    LDOS_tot_vdn = np.zeros([N_atoms, len(omega)], dtype = complex)
    
    uup_tot = psi[::4, :]
    udn_tot = psi[1::4, :]
    vup_tot = psi[2::4, :]
    vdn_tot = psi[3::4, :]
    
    for Ni in range(N_atoms):
        for n_omega in range(len(omega)):
            
            LDOS_tot_uup[Ni, n_omega] = LDOS_u(omega[n_omega], + E, uup_tot[Ni], D)
            LDOS_tot_udn[Ni, n_omega] = LDOS_u(omega[n_omega], + E, udn_tot[Ni], D)
            LDOS_tot_vup[Ni, n_omega] = LDOS_v(omega[n_omega], + E, vup_tot[Ni], D)
            LDOS_tot_vdn[Ni, n_omega] = LDOS_v(omega[n_omega], + E, vdn_tot[Ni], D)
    
    
    

    
    return(LDOS_M_u_up, LDOS_M_u_down, LDOS_M_v_up, LDOS_M_v_down, LDOS_N_u_up, LDOS_N_u_down, LDOS_N_v_up, LDOS_N_v_down, M, N,\
           LDOS_tot_uup, LDOS_tot_udn, LDOS_tot_vup, LDOS_tot_vdn)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    