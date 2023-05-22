# -*- coding: utf-8 -*-
"""
Beispiel 1 - H3.4 - VDI Wärmeatlass S. 931

"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as cp
import ht as ht

# Angaben Beispiel
da = 30e-3 #[m]
s = 1.5e-3 #[m]
di = da - 2*s #[m]
A = np.pi/4 * di**2 

p = 4.82e5 #[Pa]

q = 60000 #[W/m²]

mA = 250 #[kg/m²s]
m = mA * A

D = di
rhol = 640 # [kg/m³]
rhog = 12.5 # [kg/m³]
mul = 2.308e-4 #[Pas]
mug = 1.14e-5 #[Pas]
kl = 1.126e-1 #[W/mK]
kg = 2.82e-2 #[W/mK]
Cpl = 4420 #[J/kgK]
Cpg = 2140 #[J/kgK]
Hvap = 509700 # [J/kg]
sigma = 1.23e-2 #[kg/s²]
Psat = p
# Annahmewert für Te !!! Chen_Bennett Chen_Edelstein Liu_Winterton !!!
Te = 10
#Te = 0

# Ethanol als grobe Annäherung !!! Chen_Bennett Chen_Edelstein !!!
T = cp.PropsSI('T', 'P', p, 'Q', 1, 'Ethanol')
dPsat = cp.PropsSI('P', 'T', T + Te, 'Q', 1, 'Ethanol') - cp.PropsSI('P', 'T', T, 'Q', 1, 'Ethanol')
#dPsat =0

# Weitere Recherche
MW = 120.91 #[g/mol]
Pc = 41.4e5 #[Pa]

# Vergleichspunkte
alpha_lsg=[4899, 5042, 5215, 5409, 5616, 5833, 6055]
x_lsg = np.linspace(0, 0.3, 7)

# Berechnungspunkte
x_calc = np.linspace(0, 0.3, 100)

alpha_out = []
# , 
Methoden = ['Chen_Bennett', 'Chen_Edelstein', 'Li_Wu', 'Liu_Winterton', 'Sun_Mishima', 'Yun_Heo_Kim']

for i,methode in enumerate(Methoden):
    
    alpha_print = []
    
    for x in x_calc:
        if methode == 'Chen_Bennett':
            alpha = ht.boiling_flow.Chen_Bennett(m, x, D, rhol, rhog, mul, mug, kl, Cpl, Hvap, sigma, dPsat, Te)
        elif methode == 'Chen_Edelstein':
            alpha = ht.boiling_flow.Chen_Edelstein(m, x, D, rhol, rhog, mul, mug, kl, Cpl, Hvap, sigma, dPsat, Te)
        elif methode == 'Li_Wu':
            alpha =  ht.boiling_flow.Li_Wu(m, x, D, rhol, rhog, mul, kl, Hvap, sigma, q=q, Te=None)
        elif methode == 'Liu_Winterton':
            alpha = ht.boiling_flow.Liu_Winterton(m, x, D, rhol, rhog, mul, kl, Cpl, MW, p, Pc, Te)
        elif methode == 'Sun_Mishima':
            alpha =  ht.boiling_flow.Sun_Mishima(m, D, rhol, rhog, mul, kl, Hvap, sigma, q=q, Te=None)
        elif methode == 'Yun_Heo_Kim':
            alpha = ht.boiling_flow.Yun_Heo_Kim(m, x, D, rhol, mul, Hvap, sigma, q=q,Te=None)
        elif methode == 'Thome':
            alpha = ht.boiling_flow.Thome(m, x, D, rhol, rhog, mul, mug, kl, kg, Cpl, Cpg, Hvap, sigma, Psat, Pc, q=q, Te=None)
        else:
            alpha = 1000000
            print('Fehler')
        
        alpha_print.append(alpha)
    
    alpha_out.append(alpha_print)
    
    plt.plot(x_calc, alpha_print, label = methode)

plt.plot(x_lsg, alpha_lsg, 'ok')    
plt.title(f'Wärmeübergangskoeffizienten Te={Te} \n- berechnet, o gegeben')
plt.ylabel('alpha(x) [W/m²K]')
plt.xlabel('Dampfgehalt')
plt.legend()
    
            
          