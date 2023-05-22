# -*- coding: utf-8 -*-
"""
Praktikum Thermodynamik - Philip Besant

Ortaufgelöste Modellierung eines Verdampfers
"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as cp
import ht as ht
import pandas as pd
import heat_transfer_evaporator as hte
from scipy.integrate import solve_ivp
from scipy import constants

"""Einlesen Messdaten"""
#excel = pd.read_excel("C:\\Users\\phili\\OneDrive - Universitaet Duisburg-Essen\\Thermo Praktikum\\Daten Wärmepumpe\\Alex Daten Wärmepumpe\\Daten Maurits.xlsx", skiprows=0)
excel = pd.read_excel("E:\\OneDrive - Universitaet Duisburg-Essen\\Thermo Praktikum\\Daten Wärmepumpe\\Alex Daten Wärmepumpe\\Daten Maurits.xlsx", skiprows=0)

# Auswahl Messreihe und Berechnungsmethode
MR = 9 # 9 bis 57

""" 
0 - VDI-WA nach Alexandra Welp
 1 - Chen_Bennett
 2 - Chen_Edelstein
 3 - Li_Wu
 4 - Liu_Winterton
 5 - Sun_Mishima
 6 - Thome
 7 - Yun_Heo_Kim
 """
BR = 1

p_mdKMein = excel.loc[MR]['Druck Verdampfer ein'] * 1e5
m_mdKM = excel.loc[MR]['Massenstrom KM']
m_mdW = excel.loc[MR]['Massenstrom Verdampfer SF']

x_md = [0, 1.4, 2.8, 4.2, 5.6, 7, 8.4, 9.8, 11.2, 12.6, 14]

# Messdaten KM
T9 = excel.loc[MR]['Temperatur Verdampfer ein']
T21 = excel.loc[MR]['Temperatur 21']
T22 = excel.loc[MR]['Temperatur 22']
T23 = excel.loc[MR]['Temperatur 23']
T24 = excel.loc[MR]['Temperatur 24']
T25 = excel.loc[MR]['Temperatur 25']
T13 = excel.loc[MR]['Temperatur 13 Wasserseite']
T14 = excel.loc[MR]['Temperatur 14 Wasserseite']
T15 = excel.loc[MR]['Temperatur 15 Wasserseite']
T16 = excel.loc[MR]['Temperatur 16 Wasserseite']
T7 = excel.loc[MR]['Temperatur Verdampfer aus']
T_mdKM = [T9, T21, T13, T22, T14, T23, T15, T24, T16, T25, T7]

# Messdaten Wasser
T1 = excel.loc[MR]['Temperatur Verdampfer Vorlauf']
T30 = excel.loc[MR]['Temperatur 30']
T29 = excel.loc[MR]['Temperatur 29']
T28 = excel.loc[MR]['Temperatur 28']
T27 = excel.loc[MR]['Temperatur 27']
T26 = excel.loc[MR]['Temperatur 26']
T20 = excel.loc[MR]['Temperatur 20 Gasseite']
T19 = excel.loc[MR]['Temperatur 19 Gasseite']
T18 = excel.loc[MR]['Temperatur 18 Gasseite']
T17 = excel.loc[MR]['Temperatur 17 Gasseite']
T2 = excel.loc[MR]['Temperatur Verdampfer Rücklauf']
T_mdW = [T2, T26, T17, T27, T18, T28, T19, T29, T20, T30, T1]

xbutan = excel.loc[MR]['Molenbruch Isobutan']
xpropan = 1 - xbutan

"""Eingabe"""
### Abmessungen ###
# Rohre Dis Venzik S.70
di = 10e-3 #[m]
da = 12e-3 #[m]
Di = 16e-3 #[m]  
Da = 18e-3 #[m]
L_verdampfer = 14 #[m]

### Randbedingungen ###
Medien = ['Water','n-Propane','IsoButane'] 

# Wasser
#Vdot_w = 6 / 60 / 1000 #[m³/s]
#Tw_ein = 60 + 273.15 #[K]
#Tw_ein = 25.17 + 273.15 #[K]
#Tw_aus = 21.65 + 273.15 #[K]
p_w = 3e5 #[Pa]

Tw_aus = T_mdW[0] + 273.15 #[K]
Tw_ein = T_mdW[-1] + 273.15 #[K]

# Kältemittel 
#p_KMein = 2.75e5 #[Pa]
#mdot_KM = 0.0063 #[kg/s]
#T_KMein = 15 + 273.15 #[K]

p_KMein = p_mdKMein
mdot_KM = m_mdKM
T_KMein = T_mdKM[0] + 273.15 #[K]
#T_KMein = -5 + 273.15 #[K]

### Stoffdaten ###
h_w = cp.PropsSI('H', 'P', p_w, 'T', Tw_ein, 'water')

rho_w = cp.PropsSI('D', 'P', p_w, 'H', h_w, 'water') #[kg/m³] 
cp_w = cp.PropsSI('CPMASS', 'P', p_w, 'H', h_w, 'water') #[J/(kg*K)]
lam_w = cp.PropsSI('CONDUCTIVITY', 'P', p_w, 'H', h_w, 'water') #[W/(m*K)]
eta_w = cp.PropsSI('VISCOSITY', 'P', p_w, 'H', h_w, 'water') # [Pa*s]
Pr_w =  cp.PropsSI('PRANDTL', 'P', p_w, 'H', h_w, 'water')
ny_w = eta_w/rho_w # [m²/s]


Vdot_w = m_mdW / rho_w


# Stahl 1.4571 / Werkstoff37
rho_rohr = 7.98e3 #[kg/m³] VDI-WA D6.1-Tab2
c_rohr = 470 #[J/kg] VDI-WA D6.1-Tab4
lam_rohr = 15.4 #[W/(m*K)]


"""Funktionen"""

# Druckverlust
"""
def DruckverlustND(rho, rho_L, rho_G, m, x, l, tetha = 0):
    
    eps = x/rho_G * ( (1 + 0.12*(1-x)) * (x/rho_G + (1-x)/rho_L) + \
           (1.18*(1-x)*(constants.g*constants.Boltzmann*(rho_L-rho_G))**0.25) \
            / (m*rho_L**0.5))**-1
    
    dp_statisch = l * (rho_L*(1-eps) + rho_G*eps)*constants.g*np.sin(tetha) # H3.2-(3)
    
    dp_besch = m**2 * ((x[1]**2 / (eps*rho_G) + (1-x[1])/((1-eps)*rho*l)) \
                    - (x[0]**2 / (eps*rho_G) + (1-x[0])/((1-eps)*rho*l)))
    
    dp_ges = dp_statisch + dp_besch
    
    return dp_ges
    
"""

def DruckverlustFluid(Da, u, rho, ny, L, Di = 0):
    """
    Berechnung des Druckverlusts in glatten Durchströmten Rohren mit 
    Kreisquerschnitt
    
    Gleichungen nach VDI-Wärmeatlas, Kap.L1.2
    """
    
    D = Da-Di # L1.2-(13) Vereinfacht für Rohr/Ringspalt
    Re = u*D/ny # B1-(10)
    
    if Re > 2320:
        if Re>100000:
            print('Reynoldszahl zu hoch')
        zeta = 0.3164/(Re**(1/4)) # L1.2-(5)
        
    else:
        zeta = 64/Re # L1.2-(4)
    
    return zeta * L/D * rho*u**2 /2 # [Pa] L1.2-(1)

# Wärmeübertragungskoeffizienten
def alpha_i(d, l, Vdot, ny, lam, Pr, di=0):
    """
    Berechnung des Wärmeübergangskoeffizienten auf der Innenseite 
    eines durchströmten Rohres
    
    Gleichungen nach VDI-Wärmeatlas, Kap.G1, 11.Aufl., Springer Vieweg, 2013
    ISBN:  978-3-642-19980-6 / 978-3-642-19981-3
    """

    u = Vdot / (((d/2)**2 - (di/2)**2 )* np.pi) # [m/s]
    
    dh = d - di
    
    Re = u * dh / ny # B1-(10)

    if Re > 10000:

        Xi = (1.8 * np.log10(Re) - 1.5)**-2 # G1-(27)

        Nu = ((Xi/8) * Re * Pr)/(1 + 12.7 * (Xi/8)**0.5 \
                * (Pr**(2/3) -1))\
                * (1+(dh/l)**(2/3)) # G1-(26)

    elif Re > 2300:
        gamma = (Re - 2300)/(1e4 - 2300) # G1-(30)

        Nu_wmq2 = 1.615 * (2300 * Pr * dh/l)**(1/3) # G1-(32)
        Nu_wmq3 = (2/(1+22 * Pr))**(1/6) * (2300 * dh/l)**0.5 # G1-(33)
        Nu_wL = (49.371 + (Nu_wmq2-0.7)**3 + Nu_wmq3**3)**(1/3) # G1-(31)

        Nu_wT = ((0.0308/8) * 1e4 * Pr)/ \
                (1 + 12.7 * (0.0308/8)**0.5 * (Pr**(2/3) -1)) # G1-(37)

        Nu = (1-gamma)*Nu_wL + gamma * Nu_wT # G1-(29)
        
    else:
        Nu = 3.66 # G1-(1)

    alpha_i = Nu * lam / dh # [W/m²K] B1-(9)
       
    return alpha_i

def alpha_km(KM, m, T, P, P_w, l, Methode, xlocal, D = di):
    """
    Berechnung des Wärmeübergangskoeffizienten auf der KM-Seite
    Methode bei ND:
        0 - VDI-WA nach Alexandra Welp
        1 - Chen_Bennett
        2 - Chen_Edelstein
        3 - Li_Wu
        4 - Liu_Winterton
        5 - Sun_Mishima
        6 - Thome
        7 - Yun_Heo_Kim
    """
    
    #print(f'Enthalpie: {h}')
    T_KM = T[0]
    T_W = T[1]
    
    # Ethalpie der Flüssigkeit und des Gasen im Nassdampfgebiet in [J/kg]
    H_KM_G = cp.PropsSI('H', 'P', P, 'Q', 1, KM)
    H_KM_L = cp.PropsSI('H', 'P', P, 'Q', 0, KM)
    
    # Stoffeigenschaften in abhängigkeit der Enthalpie
    rho_KM = cp.PropsSI('D', 'T', T_KM, 'P', P, KM)
    eta_KM = cp.PropsSI('V', 'P', P, 'T', T_KM, KM)
    lam_KM = cp.PropsSI('CONDUCTIVITY', 'P', P, 'T', T_KM, KM)
    Pr_KM = cp.PropsSI('PRANDTL', 'P', P, 'T', T_KM, KM)
    
    Tsiede0_KM = cp.PropsSI('T', 'P', p_KMein, 'Q', 0, KM)
    Tsiede1_KM = cp.PropsSI('T', 'P', p_KMein, 'Q', 1, KM)
    
    if T_KM > Tsiede0_KM and T_KM < Tsiede1_KM:
        
        #print('alpha in ND')
        
        x = cp.PropsSI('Q', 'P', P, 'T', T_KM, KM) # Dampfgehalt x
        
        # Dichte der Flüssigkeit und des Gasen im Nassdampfgebiet in [kg/m³]
        rhol = cp.PropsSI('D', 'P', P, 'Q', 0, KM)
        rhog = cp.PropsSI('D', 'P', P, 'Q', 1, KM)
        
        # Dyn. Viskosität der Flüssigkeit und des Gasen im Nassdampfgebiet in [Pa*s]
        mul = cp.PropsSI('V', 'P', P, 'Q', 0, KM)
        mug = cp.PropsSI('V', 'P', P, 'Q', 1, KM)
        
        # Wärmeleitfähigkeit der Flüssigkeit und des Gasen im Nassdampfgebiet in [W/(m*K)]
        kl = cp.PropsSI('CONDUCTIVITY', 'P', P, 'Q', 0, KM)
        kg = cp.PropsSI('CONDUCTIVITY', 'P', P, 'Q', 1, KM)
        
        # Wärmekapazitäten (isobar) der Flüssigkeit und des Gasen im Nassdampfgebiet in [J/(kg*K)]
        Cpl = cp.PropsSI('CPMASS', 'P', P, 'T', Tsiede0_KM-0.1, KM)
        Cpg = cp.PropsSI('CP0MASS', 'P', P, 'Q', 1, KM)
        
        # Verdampfungsenthalpie in [J/kg]
        Hvap = H_KM_G - H_KM_L
        
        # Oberflächenspannung der Flüssigkeit in [N/m]
        #sigma = cp.PropsSI('SURFACE_TENSION', 'P', P, 'Q', 0, KM)
        sigma = xbutan * cp.PropsSI('SURFACE_TENSION', 'P', P, 'Q', 0, 'IsoButane') \
                + xpropan * cp.PropsSI('SURFACE_TENSION', 'P', P, 'Q', 0, 'n-Propane')
        
        # Siedetemperatur des KM
        Tsiede_KM = T_KM
        
        # Überhitzung der Wand in [K]
        Te = T_W - T_KM
        
        # Differenz im Sättigungsdruck in [Pa]
        dPsat = cp.PropsSI('P', 'T', T_W, 'Q', 1, KM) - cp.PropsSI('P', 'T', Tsiede_KM, 'Q', 1, KM)
        Psat = P
        
        # Kritischer Druck des Kältemittels
        Pc = cp.PropsSI('PCRIT', KM)
        
        # Molare Masse des Fluids in [g/mol]
        MW = cp.PropsSI('MOLARMASS', KM) / 1000
        
        # Prandtlzahl der Flüssigkeit und des Gasen im Nassdampfgebiet
        Pr_KM_G = cp.PropsSI('PRANDTL', 'P', P, 'Q', 1, KM)
        Pr_KM_L = cp.PropsSI('PRANDTL', 'P', P, 'Q', 0, KM)
        
        if Methode == 0:
            
            p1 = hte.PointND(T_KM, x, KM, 'Testpunkt Isobutan')
            alpha_KM = p1.alpha_sieden(m, D, xlocal)
        
        elif Methode == 1:
            
            alpha_KM = ht.boiling_flow.Chen_Bennett(m, x, D, rhol, rhog, mul, mug, kl, Cpl, Hvap, sigma, dPsat, Te)
            
        elif Methode == 2:
            
            alpha_KM = ht.boiling_flow.Chen_Edelstein(m, x, D, rhol, rhog, mul, mug, kl, Cpl, Hvap, sigma, dPsat, Te)
        
        elif Methode == 3:
            
            alpha_KM =  ht.boiling_flow.Li_Wu(m, x, D, rhol, rhog, mul, kl, Hvap, sigma, q=None, Te = Te)
            
        elif Methode == 4:
            
            alpha_KM = ht.boiling_flow.Liu_Winterton(m, x, D, rhol, rhog, mul, kl, Cpl, MW, P, Pc, Te)
            
        elif Methode == 5:
            
            alpha_KM =  ht.boiling_flow.Sun_Mishima(m, D, rhol, rhog, mul, kl, Hvap, sigma, q=None, Te=Te)
            
        elif Methode == 6:
            
            alpha_KM = ht.boiling_flow.Thome(m, x, D, rhol, rhog, mul, mug, kl, kg, Cpl, Cpg, Hvap, sigma, Psat, Pc, Te = Te)
            
        elif Methode == 7:
            
            alpha_KM = ht.boiling_flow.Yun_Heo_Kim(m, x, D, rhol, mul, Hvap, sigma, Te = Te)
        
        else:
            print("Auswahl Berechnungsmethode für Nassdampfgebiet Fehlerhaft")
        
    else:
        alpha_KM = alpha_i(D, l, m/rho_KM, eta_KM/rho_KM, lam_KM, Pr_KM)
        #print('alpha nach alpha_i')
    # print(alpha_KM)
    return alpha_KM

# Wärmeübertrager
def Wärmeübertrager(x, T, KM, p_km, p_w, mdot_KM, mdot_W, L, Methode, di=di, da=da):
    
    T_KM = T[0]
    T_W = T[1]
    delta_T = T_W - T_KM
    if delta_T.any() < 0:
        raise ValueError("no evaporator")
    
    alphaKM = alpha_km(KM, mdot_KM, T, p_km, p_w, L, Methode, x)
    #print(alphaKM)
    
    Rl_W = 1/(alpha_w*np.pi*da)
    Rl_Rohr = np.log(da/di) / (2*np.pi*lam_rohr)
    Rl_KM = 1/(alphaKM*np.pi*di)

    Rl_ges = Rl_W + Rl_Rohr + Rl_KM
    
    dT_KM = 1/Rl_ges * delta_T / mdot_KM / cp.PropsSI('CPMASS', 'P', p_km, 'T', T_KM, KM) 
    dT_W = 1/Rl_ges * delta_T / mdot_W / cp_w
    
    return np.array([dT_KM, dT_W])

      
"""Berechnung"""    

alpha_w = alpha_i(Di, L_verdampfer, Vdot_w, ny_w, lam_w, Pr_w, da)
print(f'Wärmeübergangskoeffizient Wasser: {alpha_w}')

Auflösung = 100
x_calc = np.linspace(0, L_verdampfer, Auflösung)

KM = f'IsoButane[{xbutan}]&n-Propane[{xpropan}]'

Tsiede_KM = cp.PropsSI('T', 'P', p_KMein, 'Q', 0, KM) - 273.15
print(f'Siedetemperatur start KM: {Tsiede_KM}')
Ttau_KM = cp.PropsSI('T', 'P', p_KMein, 'Q', 1, KM) - 273.15
print(f'Siedetemperatur ende KM: {Ttau_KM}')

hKM_ein = cp.PropsSI('H', 'P', p_KMein, 'T', T_KMein, KM)
hW_aus = cp.PropsSI('H', 'P', p_w, 'T', Tw_aus, Medien[0])

T0 = np.array([T_KMein, Tw_aus])


# def Wärmeübertrager(x, T, KM, p_km, p_w, mdot_KM, mdot_W, L, Methode, di=di, da=da):
# solve_ivp(Kollektor,(t_calc[0],t_calc[-1]),T0,t_eval=t_calc,args=(),atol=1e-9,rtol=1e-6)
res = solve_ivp(lambda x, T: Wärmeübertrager(x, T, KM, p_KMein, p_w, mdot_KM, m_mdW, L_verdampfer, BR),(x_calc[0],x_calc[-1]),T0,t_eval=x_calc,args=(),atol=1e-9,rtol=1e-6)

x_res = res.t
TKM_res = res.y[0]
TW_res = res.y[1]




print('')
T_KMausber = TKM_res[-1]
T_Weinber = TW_res[-1]
print(f'berechnete Wassereintrittstemperatur: {T_Weinber} K')
h_KMausber = cp.PropsSI('H', 'P', p_KMein, 'T', T_KMausber, KM)
h_Weinber = cp.PropsSI('H', 'P', p_w, 'T', T_Weinber, 'water')
Q_Wber = (h_Weinber - hW_aus) * m_mdW
Q_KMber = (h_KMausber - hKM_ein) * m_mdKM
print(f'berechneter Wärmestrom Wasser: H_w = {Q_Wber} W')
print(f'berechneter Wärmestrom Kältemittel: H_km = {Q_KMber} W')

hKM_ausmd = cp.PropsSI('H', 'P', p_KMein, 'T', T_mdKM[-1] + 273.15, KM)
Q_Wmd = (h_w - hW_aus) * m_mdW
Q_KMmd = (hKM_ausmd - hKM_ein) * m_mdKM
print(f'gemessener Wärmestrom Wasser: H_w = {Q_Wmd} W')
print(f'gemessener Wärmestrom Kältemittel: H_km = {Q_KMmd} W')



    
plt.figure(0)
plt.plot(x_res, TKM_res - 273.15, 'b', label = 'Kältemittel')
plt.plot(x_res, TW_res - 273.15, 'r', label = 'Wasser')
plt.plot(x_md, T_mdKM, 'ob', label = 'MD-Kältemittel')
plt.plot(x_md, T_mdW, 'or', label = 'MD-Wasser',)
plt.hlines(Tsiede_KM, 0, L_verdampfer, 'k')
plt.hlines(Ttau_KM, 0, L_verdampfer, 'k')

plt.title(f'Temperaturverläufe Im Verdampfer - BR: {BR}')
plt.ylabel('Temperatur [°C]')
plt.xlabel('Länge [m]')
plt.legend()
