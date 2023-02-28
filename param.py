#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 11:06:53 2023

@author: cristina
"""

t = 1.0#hopping term
mu = 2.0#chemical potential
Delta = 0.5#SC gap
J = 2.0#magnetic coupling
S = 1.0#spin
alpha = 1.0#SOC
U = 0.0#potential scattering
state = 'FM'#chain spin

'''size and k-vector'''
N_atoms = 51#number of atoms
a = 1.0#lattice param

N_imp = 21#number of impurities
step = 1#space between impurities
Nk = 301#number of k-points

Dynes = 0.01