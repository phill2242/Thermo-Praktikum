"""
Praktikum Thermodynamik - Philip Besant

Ortaufgelöste Modellierung eines Verdampfers
"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as cp
import ht as ht
#import heat_transfer_evaporator as hte
from scipy.integrate import solve_bvp
from scipy import constants
import pandas as pd

#### CoolProp
# Eigenschaften
# http://www.coolprop.org/coolprop/HighLevelAPI.html#table-of-string-inputs-to-propssi-function
# Fluide
# http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids

"""Einlesen Messdaten"""
#excel = pd.read_excel("C:\\Users\\phili\\OneDrive - Universitaet Duisburg-Essen\\Thermo Praktikum\\Daten Wärmepumpe\\Alex Daten Wärmepumpe\\Daten Maurits.xlsx", skiprows=0)
excel = pd.read_excel("E:\\OneDrive - Universitaet Duisburg-Essen\\Thermo Praktikum\\Daten Wärmepumpe\\Alex Daten Wärmepumpe\\Daten Maurits.xlsx", skiprows=0)

# Auswahl Messreihe und Berechnungsmethode
MR = 0 # 0 bis 6 ohne 3
BR = 7

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

### Stoffdaten ###
h_w = cp.PropsSI('H', 'P', p_w, 'T', Tw_aus, 'water')

print(f'Eintrittsenthalpie W: {h_w}')

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

# Zurücksetzen Zähler
k = 0


"""Funktionen"""

# Druckverlust

def DruckverlustND(rho, rho_L, rho_G, m, x, l, tetha = 0):
    
    eps = x/rho_G * ( (1 + 0.12*(1-x)) * (x/rho_G + (1-x)/rho_L) + \
           (1.18*(1-x)*(constants.g*constants.Boltzmann*(rho_L-rho_G))**0.25) \
            / (m*rho_L**0.5))**-1
    
    dp_statisch = l * (rho_L*(1-eps) + rho_G*eps)*constants.g*np.sin(tetha) # H3.2-(3)
    
    dp_besch = m**2 * ((x[1]**2 / (eps*rho_G) + (1-x[1])/((1-eps)*rho*l)) \
                    - (x[0]**2 / (eps*rho_G) + (1-x[0])/((1-eps)*rho*l)))
    
    dp_ges = dp_statisch + dp_besch
    
    return dp_ges
    

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

def alpha_i(d, l, Vdot, ny, lam, Pr, di=0):
    """
    Berechnung des Wärmeübergangskoeffizienten auf der Innenseite 
    eines durchströmten Rohres
    
    Gleichungen nach VDI-Wärmeatlas, Kap.G6, 11.Aufl., Springer Vieweg, 2013
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

def alpha_km(KM, m, h, P, P_w, l, Methode, D = di):
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
    
    global k
    
    #print(f'Enthalpie: {h}')
    H_KM = h[0]
    H_W = h[1]
    
    # Ethalpie der Flüssigkeit und des Gasen im Nassdampfgebiet in [J/kg]
    H_KM_G = cp.PropsSI('H', 'P', P, 'Q', 1, KM)
    H_KM_L = cp.PropsSI('H', 'P', P, 'Q', 0, KM)
    
    # Stoffeigenschaften in abhängigkeit der Enthalpie
    rho_KM = cp.PropsSI('D', 'H', H_KM, 'P', P, KM)
    eta_KM = cp.PropsSI('V', 'P', P, 'H', H_KM, KM)
    lam_KM = cp.PropsSI('CONDUCTIVITY', 'P', P, 'H', H_KM, KM)
    Pr_KM = cp.PropsSI('PRANDTL', 'P', P, 'H', H_KM, KM)
    
    
    if H_KM > H_KM_L and H_KM < H_KM_G:
        
        x = cp.PropsSI('Q', 'P', P, 'H', H_KM, KM) # Dampfgehalt x
        
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
        Cpl = cp.PropsSI('CPMASS', 'P', P, 'H', H_KM_L-1, KM)
        Cpg = cp.PropsSI('CP0MASS', 'P', P, 'Q', 1, KM)
        
        # Verdampfungsenthalpie in [J/kg]
        Hvap = H_KM_G - H_KM_L
        
        # Oberflächenspannung der Flüssigkeit in [N/m]
        sigma = cp.PropsSI('SURFACE_TENSION', 'P', P, 'Q', 0, KM)
        
        # Siedetemperatur des KM
        Tsiede_KM = cp.PropsSI('T', 'P', P, 'H', H_KM, KM)
        
        # Wassertemperatur 
        T_W = cp.PropsSI('T', 'H', H_W, 'P', P_w, 'Water')
        
        # Überhitzung der Wand in [K]
        Te = T_W - Tsiede_KM
        
        # Differenz im Sättigungsdruck in [Pa]
        dPsat = cp.PropsSI('P', 'T', T_W, 'Q', 0.5, KM) - cp.PropsSI('P', 'T', Tsiede_KM, 'Q', 0.5, KM)
        Psat = P
        
        # Kritischer Druck des Kältemittels
        Pc = cp.PropsSI('PCRIT', KM)
        
        # Molare Masse des Fluids in [g/mol]
        MW = cp.PropsSI('MOLARMASS', KM) / 1000
        
        # Prandtlzahl der Flüssigkeit und des Gasen im Nassdampfgebiet
        Pr_KM_G = cp.PropsSI('PRANDTL', 'P', P, 'Q', 1, KM)
        Pr_KM_L = cp.PropsSI('PRANDTL', 'P', P, 'Q', 0, KM)
        
        if Methode == 0:
            
            #alpha_KM = hte.alpha_sieden()
            alpha_KM = 1000
        
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
        k+=1
        #print('alpha nach alpha_i')
    #print(alpha_KM)
    return alpha_KM

def Wärmeübertrager(h, KM, p_km, p_w, mdot_KM, mdot_W, L, Aufl, Methode, di=di, da=da):
    
    T_KM = cp.PropsSI('T', 'P', p_km, 'H', h[0], KM)
    T_W = cp.PropsSI('T', 'P', p_w, 'H', h[1], 'Water')
    delta_T = T_W - T_KM
    
    alphaKM = alpha_km(KM, mdot_KM, h, p_km, p_w, L, Methode)
    #print(f'alphaKM: {alphaKM}')
    
    Rl_W = 1/(alpha_w*np.pi*da)
    Rl_Rohr = np.log(da/di) / (2*np.pi*lam_rohr)
    Rl_KM = 1/(alphaKM*np.pi*di)

    Rl_ges = Rl_W + Rl_Rohr + Rl_KM    
    
    dh_KM = 1/Rl_ges * L/Aufl * delta_T / mdot_KM
    dh_W = 1/Rl_ges * L/Aufl * delta_T / mdot_W
    
    return np.array([dh_KM, dh_W, alphaKM])


      
"""Berechnung"""    


Auflösung = 2000
x_calc = np.linspace(0, L_verdampfer, Auflösung)
L_stück = L_verdampfer / Auflösung

A_Stück = (Di + da)/2 * np.pi * L_stück 

KM = 'n-Propane'

alpha_w = alpha_i(Di, L_verdampfer, Vdot_w, ny_w, lam_w, Pr_w, di)
print(f'Wärmeübergangskoeffizient Wasser: {alpha_w}')

Tsiede_KM = cp.PropsSI('T', 'P', p_KMein, 'Q', 0, KM) - 273.15
print(f'Siedetemperatur start KM: {Tsiede_KM}')
Tsiede_KM = cp.PropsSI('T', 'P', p_KMein, 'Q', 1, KM) - 273.15
print(f'Siedetemperatur ende KM: {Tsiede_KM}')

hKM_ein = cp.PropsSI('H', 'P', p_KMein, 'T', T_KMein, KM)
print(f'Eintrittsenthalpie KM: {hKM_ein}')
hW_ein = cp.PropsSI('H', 'P', p_w, 'T', Tw_ein, Medien[0])
hW_aus = cp.PropsSI('H', 'P', p_w, 'T', Tw_aus, Medien[0])
print(f'Austrittsenthalpie W: {hW_aus}')

mdot_W = Vdot_w * rho_w

h = np.zeros((2, x_calc.size))
h[0, :] = hKM_ein
h[1, :] = hW_aus
alpha_calc = []

for n in range(0, x_calc.size - 1):
    
    h_calc = np.array([h[0][n], h[1][n]])
    
    WÜ = Wärmeübertrager(h_calc, KM, p_KMein, p_w, mdot_KM, mdot_W, L_verdampfer, Auflösung, BR)
    
    h[0][n+1] = h[0][n] + WÜ[0]
    h[1][n+1] = h[1][n] + WÜ[1]
    alpha_calc.append(WÜ[2])
    
# Auslesen Ergebnissse
hKM_res = h[0]
hW_res = h[1]

# Berechnung Temperatur
TKM_res = cp.PropsSI('T', 'P', p_KMein, 'H', hKM_res, KM) - 273.15
TW_res = cp.PropsSI('T', 'P', p_w, 'H', hW_res, 'Water') - 273.15

# Darstellung Temperaturverlauf
plt.figure(0)
plt.plot(x_calc, TKM_res, 'b', label = 'Kältemittel')
plt.plot(x_calc, TW_res, 'r',label = 'Wasser')

plt.plot(x_md, T_mdKM, 'ob', label = 'MD-Kältemittel')
plt.plot(x_md, T_mdW, 'or', label = 'MD-Wasser',)

plt.title(f'Temperaturverläufe Im Verdampfer MR{MR} BR{BR}')
plt.ylabel('Temperatur [°C]')
plt.xlabel('Länge [m]')
plt.legend()

# Kontrolle zu Nassdampfanteil
ND_Anteil = round((1 - k/Auflösung) *100)
print(f'Anteil ND: {ND_Anteil}%')