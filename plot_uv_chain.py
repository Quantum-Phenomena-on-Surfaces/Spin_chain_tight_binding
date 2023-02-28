# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 16:46:37 2020

@author: nodar
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 11})
plt.rcParams['figure.figsize'] = [12.8, 9.6]

def plot_uv(u_up, v_up, u_down, v_down, u_up2, v_up2, u_down2, v_down2,\
            min_E, min_E2, M, N, N_atoms):
    
    min_E = np.real(min_E)
    min_E2 = np.real(min_E2)
    
    fig = plt.figure(2)
    #ax = fig.add_subplot(111)# The big subplot
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    
    ax1.plot(np.real(u_up), '.-', label = 'u up')
    ax1.plot(np.real(u_down), '.-', label = 'u down')
    ax1.plot(np.real(v_up), '.-', label = 'v up')
    ax1.plot(np.real(v_down), '.-', label = 'v down')
    #ax1.title('Real part 1')
    ax1.legend()
    ax1.title.set_text('Eigenenergy=%f'%min_E)
    ax1.set_ylabel('Real part')
    ax1.axvline(x=M, linestyle = '--', linewidth=1.0, color='k')
    ax1.axvline(x=N, linestyle = '--', linewidth=1.0, color='k')
    ax1.set_xlim(-3, N_atoms+2)
    
    ax3.plot(np.imag(u_up), '.-', label = 'u up')
    ax3.plot(np.imag(u_down), '.-', label = 'u down')
    ax3.plot(np.imag(v_up), '.-', label = 'v up')
    ax3.plot(np.imag(v_down), '.-', label = 'v down')
    ax3.set_ylabel('Imag part')
    ax3.legend()
    ax3.axvline(x=M, linestyle = '--', linewidth=1.0, color='k')
    ax3.axvline(x=N, linestyle = '--', linewidth=1.0, color='k')
    ax3.set_xlim(-3, N_atoms+2)
    
    ax2.plot(np.real(u_up2), '.-', label = 'u up')
    ax2.plot(np.real(u_down2), '.-', label = 'u down')
    ax2.plot(np.real(v_up2), '.-', label = 'v up')
    ax2.plot(np.real(v_down2), '.-', label = 'v down')
    #ax3.title('Real part 2')
    ax2.legend()
    ax2.title.set_text('Eigenenergy=%f'%min_E2)
    ax2.axvline(x=M, linestyle = '--', linewidth=1.0, color='k')
    ax2.axvline(x=N, linestyle = '--', linewidth=1.0, color='k')
    ax2.set_xlim(-3, N_atoms+2)
    
    ax4.plot(np.imag(u_up2), '.-', label = 'u up')
    ax4.plot(np.imag(u_down2), '.-', label = 'u down')
    ax4.plot(np.imag(v_up2), '.-', label = 'v up')
    ax4.plot(np.imag(v_down2), '.-', label = 'v down')
    #ax4.title('Imag part 2')
    ax4.legend()
    ax4.set_xlabel('site')
    ax4.axvline(x=M, linestyle = '--', linewidth=1.0, color='k')
    ax4.axvline(x=N, linestyle = '--', linewidth=1.0, color='k')
    ax4.set_xlim(-3, N_atoms+2)
    
    #save fig
    #plt.savefig('Figures/u_v_real_img.png', dpi=300, bbox_inches='tight')
    plt.savefig('u_v_real_img.png', dpi=300, bbox_inches='tight')
    
    fig = plt.figure(3)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)    
    
    ax1.plot(np.abs(u_up), '.-', label = 'u up')
    ax1.plot(np.abs(u_down), '.-', label = 'u down')
    ax1.plot(np.abs(v_up), '.-', label = 'v up')
    ax1.plot(np.abs(v_down), '.-', label = 'v down')
    ax1.legend()
    ax1.set_title('Eigenenergy=%f'%min_E)
    ax1.set_ylabel('Modulus')
    ax1.axvline(x=M, linestyle = '--', linewidth=1.0, color='k')
    ax1.axvline(x=N, linestyle = '--', linewidth=1.0,color='k')
    ax1.set_xlim(-3, N_atoms+2)    
    
    ax2.plot(np.abs(u_up2), '.-', label = 'u up')
    ax2.plot(np.abs(u_down2), '.-', label = 'u down')
    ax2.plot(np.abs(v_up2), '.-', label = 'v up')
    ax2.plot(np.abs(v_down2), '.-', label = 'v down')
    #ax2.title('Modulus 2')
    ax2.legend()
    ax2.set_title('Eigenenergy=%f'%min_E2)
    ax2.axvline(x=M, linestyle = '--', linewidth=1.0, color='k')
    ax2.axvline(x=N, linestyle = '--', linewidth=1.0, color='k')
    ax2.set_xlim(-3, N_atoms+2)    
    
    ax3.plot(np.abs(u_up)+np.abs(v_down), '.-', label = 'u up+v_down')
    ax3.plot(np.abs(u_down)+np.abs(v_up), '.-', label = 'u down+v_up')
    #ax3.title('Modulus 1')
    ax3.legend()
    ax3.set_xlabel('site')
    ax3.set_ylabel('Modulus')
    ax3.axvline(x=M, linestyle = '--', linewidth=1.0, color='k')
    ax3.axvline(x=N, linestyle = '--', linewidth=1.0, color='k')
    ax3.set_xlim(-3, N_atoms+2)    
    
    ax4.plot(np.abs(u_up2)+np.abs(v_down2), '.-', label = 'u up+v_down')
    ax4.plot(np.abs(u_down2)+np.abs(v_up2), '.-', label = 'u down+v_up')
    #ax3.title('Modulus 1')
    ax4.legend()
    ax4.set_xlabel('site')
    ax4.axvline(x=M, linestyle = '--', linewidth=1.0, color='k')
    ax4.axvline(x=N, linestyle = '--', linewidth=1.0, color='k')
    ax4.set_xlim(-3, N_atoms+2)
    
    #save fig
    #plt.savefig('Figures/u_v_mod.png', dpi=300, bbox_inches='tight')
    plt.savefig('u_v_mod.png', dpi=300, bbox_inches='tight')    
    
    ###savedata
    #u_up_mod =  np.abs(u_up)
    #u_dn_mod =  np.abs(u_down)
    #v_up_mod =  np.abs(v_up)
    #v_dn_mod =  np.abs(v_down)
    
    #np.savetxt('u_up_mod.txt', u_up_mod)
    #np.savetxt('u_down_mod.txt', u_dn_mod)
    #np.savetxt('v_up_mod.txt', v_up_mod)
    #np.savetxt('v_dn_mod.txt', v_dn_mod)
    
    
 
    
    
    
    
    
    
    
    
    
    
    
    
    
    