# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 18:04:09 2020

@author: nodar
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 11})
plt.rcParams['figure.figsize'] = [12.8, 9.6]

def LDOS(LDOS_M_u_up, LDOS_M_u_down, LDOS_M_v_up, LDOS_M_v_down, LDOS_N_u_up, LDOS_N_u_down, LDOS_N_v_up, LDOS_N_v_down, M, N, omega):
    
    #omega = omega*27211.6
    
    fig = plt.figure(4)
    fig.suptitle('First impurity of the chain')
    #ax = fig.add_subplot(111)# The big subplot
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.plot(omega, LDOS_M_u_up, label = 'u_up')
    ax1.plot(omega, LDOS_M_u_down, label = 'u_down')
    ax1.plot(omega, LDOS_M_v_up, label = 'v_up')
    ax1.plot(omega, LDOS_M_v_down, label = 'v_down')
    ax1.legend()
    ax1.set_ylabel('LDOS')
    ax1.set_xlim(min(omega),max(omega))
    
    ax2.plot(omega, LDOS_M_u_up+LDOS_M_v_up, label = 'u+v_up')
    ax2.plot(omega, LDOS_M_u_down+LDOS_M_v_down, label = 'u+v_down')
    ax2.legend()
    ax2.set_xlim(min(omega),max(omega))
    
    ax3.plot(omega, LDOS_M_u_up+LDOS_M_v_down, label = 'u_up+v_down')
    ax3.plot(omega, LDOS_M_u_down+LDOS_M_v_up, label = 'u_down+v_up')
    ax3.legend()
    ax3.set_xlabel('Energy')
    ax3.set_ylabel('LDOS')
    ax3.set_xlim(min(omega),max(omega))
    
    ax4.plot(omega, LDOS_M_u_up+LDOS_M_u_down+LDOS_M_v_up+LDOS_M_v_down, label = 'total')
    ax4.legend()
    ax4.set_xlabel('Energy')
    ax4.set_xlim(min(omega),max(omega))
    
    #np.savetxt('total_LDOS.txt', LDOS_M_u_up+LDOS_M_u_down+LDOS_M_v_up+LDOS_M_v_down)
    #np.savetxt('omega.txt', omega)
    
    #save fig
    #plt.savefig('Figures/LDOS_first_imp.png', dpi=300, bbox_inches='tight')
    plt.savefig('LDOS_first_imp.png', dpi=300, bbox_inches='tight')
    
    fig = plt.figure(5)
    fig.suptitle('Middle of the chain')
    #ax = fig.add_subplot(111)# The big subplot
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.plot(omega, LDOS_N_u_up, label = 'u_up')
    ax1.plot(omega, LDOS_N_u_down, label = 'u_down')
    ax1.plot(omega, LDOS_N_v_up, label = 'v_up')
    ax1.plot(omega, LDOS_N_v_down, label = 'v_down')
    ax1.legend()
    ax1.set_ylabel('LDOS')
    ax1.set_xlim(min(omega),max(omega))
    
    ax2.plot(omega, LDOS_N_u_up+LDOS_N_v_up, label = 'u+v_up')
    ax2.plot(omega, LDOS_N_u_down+LDOS_N_v_down, label = 'u+v_down')
    ax2.legend()
    ax2.set_xlim(min(omega),max(omega))
    
    ax3.plot(omega, LDOS_N_u_up+LDOS_N_v_down, label = 'u_up+v_down')
    ax3.plot(omega, LDOS_N_u_down+LDOS_N_v_up, label = 'u_down+v_up')
    ax3.legend()
    ax3.set_xlabel('Energy')
    ax3.set_ylabel('LDOS')
    ax3.set_xlim(min(omega),max(omega))
    
    ax4.plot(omega, LDOS_N_u_up+LDOS_N_u_down+LDOS_N_v_up+LDOS_N_v_down, label = 'total')
    ax4.legend()
    ax4.set_xlabel('Energy')
    ax4.set_xlim(min(omega),max(omega))
    
    #save fig
    #plt.savefig('Figures/LDOS_last_imp.png', dpi=300, bbox_inches='tight')
    plt.savefig('LDOS_last_imp.png', dpi=300, bbox_inches='tight')
    
    
    ###save data
    #np.savetxt('LDOS_uup_M.txt', LDOS_M_u_up)
    #np.savetxt('LDOS_udn_M.txt', LDOS_M_u_down)
    #np.savetxt('LDOS_vup_M.txt', LDOS_M_v_up)
    #np.savetxt('LDOS_vdn_M.txt', LDOS_M_v_down)
    
    #np.savetxt('LDOS_uup_N.txt', LDOS_N_u_up)
    #np.savetxt('LDOS_udn_N.txt', LDOS_N_u_down)
    #np.savetxt('LDOS_vup_N.txt', LDOS_N_v_up)
    #np.savetxt('LDOS_vdn_N.txt', LDOS_N_v_down)
    
    #np.savetxt('omega.txt', omega)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    