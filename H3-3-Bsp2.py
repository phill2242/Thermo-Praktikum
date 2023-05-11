# -*- coding: utf-8 -*-
"""
Created on Wed May 10 20:44:32 2023

@author: Philip
"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as cp
import ht as ht

mA = 3000
d = 7e-3
A = np.pi/4 * d**2 
p = 10.05e5
T = cp.PropsSI('T', 'P', p, 'Q', 0, 'water')


m = mA*A
x = 0.5
D = d
Te = 200+273.15 - T
rhol = cp.PropsSI('D', 'P', p, 'Q', 0, 'water')
rhog = cp.PropsSI('D', 'P', p, 'Q', 1, 'water')
mul = cp.PropsSI('V', 'P', p, 'Q', 0, 'water')
mug = cp.PropsSI('V', 'P', p, 'Q', 1, 'water')
kl = cp.PropsSI('CONDUCTIVITY', 'P', p, 'Q', 0, 'water')
kg = cp.PropsSI('CONDUCTIVITY', 'P', p, 'Q', 1, 'water')
Cpl = cp.PropsSI('CPMASS', 'P', p, 'T', T-1, 'water')
Cpg = cp.PropsSI('CP0MASS', 'P', p, 'Q', 1, 'water')
Hvap = cp.PropsSI('H', 'P', p, 'Q', 1, 'water') - cp.PropsSI('H', 'P', p, 'Q', 0, 'water')
sigma = cp.PropsSI('SURFACE_TENSION', 'P', p, 'Q', 0, 'water')
dPsat = cp.PropsSI('P', 'T', Te + T, 'Q', 0.5, 'water') - cp.PropsSI('P', 'T', 180+273.15, 'Q', 0.5, 'water')
Psat = p

Pc = cp.PropsSI('PCRIT', 'water')
MW = cp.PropsSI('MOLARMASS', 'water') / 1000

alpha = ht.boiling_flow.Chen_Bennett(m, x, D, rhol, rhog, mul, mug, kl, Cpl, Hvap, sigma, dPsat, Te)
print('Chen_Bennett')
print(alpha)
faktor = alpha/88000
print(f'faktor: {faktor}')
alpha = ht.boiling_flow.Chen_Edelstein(m, x, D, rhol, rhog, mul, mug, kl, Cpl, Hvap, sigma, dPsat, Te)
print('Chen_Edelstein')
print(alpha)
faktor = alpha/88000
print(f'faktor: {faktor}')
alpha =  ht.boiling_flow.Li_Wu(m, x, D, rhol, rhog, mul, kl, Hvap, sigma, q=None, Te = Te)
print('Li_Wu')
print(alpha)
faktor = alpha/88000
print(f'faktor: {faktor}')
alpha = ht.boiling_flow.Liu_Winterton(m, x, D, rhol, rhog, mul, kl, Cpl, MW, p, Pc, Te)
print('Liu_Winterton')
print(alpha)
faktor = alpha/88000
print(f'faktor: {faktor}')
alpha =  ht.boiling_flow.Sun_Mishima(m, D, rhol, rhog, mul, kl, Hvap, sigma, q=None, Te=Te)
print('Sun_Mishima')
print(alpha)
faktor = alpha/88000
print(f'faktor: {faktor}')
#alpha = ht.boiling_flow.Thome(m, x, D, rhol, rhog, mul, mug, kl, kg, Cpl, Cpg, Hvap, sigma, Psat, Pc, q=None, Te=Te)
#print('Thome')
#print(alpha)
alpha = ht.boiling_flow.Yun_Heo_Kim(m, x, D, rhol, mul, Hvap, sigma, Te = Te)
print('Yun_Heo_Kim')
print(alpha)
faktor = alpha/88000
print(f'faktor: {faktor}')