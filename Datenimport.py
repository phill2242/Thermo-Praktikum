# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 09:35:20 2023

@author: Philip
"""

import pandas as pd
import matplotlib.pyplot as plt
import CoolProp.CoolProp as cp

#csvdata= pd.read_csv("E:\\OneDrive - Universitaet Duisburg-Essen\\Thermo Praktikum\\Daten Wärmepumpe\\Alex Daten Wärmepumpe\\Simon Schmeing\\Messdaten\\Messung 20\\3 - Kopie.csv", header = None)

#names = ['Zeit', 'Temp 9', 'Druck 4', 'Temp 7', 'Druck 1', 'Temp 1', 'Temp 2', 'Temp 17', 'Temp 18', 'Temp 19', 'Temp 20', 'Temp 13', 'Temp 14', 'Temp 15', 'Temp 16', 'Temp 6', 'Druck 3', 'Temp 10', 'Druck 2', 'Temp 3', 'Temp 5', 'Temp 11', 'Temp 8', 'Temp 12', 'Temp 4', 'Leistung', 'Druck 6', 'Druck 5', 'Verd.aus Temp', 'Temp 33', 'Temp 21', 'Temp 26', 'Temp 28', 'Temp 27', 'Temp 29', 'Temp 25', 'Temp 22', 'Temp 23', 'Temp 24', 'Temp 30', 'Temp 31', 'Temp 32', 'Temp 35', 'Temp 34']

#excel = pd.read_excel("E:\\OneDrive - Universitaet Duisburg-Essen\\Thermo Praktikum\\Daten Wärmepumpe\\Alex Daten Wärmepumpe\\Daten Maurits.xlsx", skiprows=0)
excel = pd.read_excel("C:\\Users\\phili\\OneDrive - Universitaet Duisburg-Essen\\Thermo Praktikum\\Daten Wärmepumpe\\Alex Daten Wärmepumpe\\Daten Maurits.xlsx", skiprows=0)

MR = 9

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

plt.figure(0)
plt.plot(x_md, T_mdKM, 'ob', label = 'MD-Kältemittel')
plt.plot(x_md, T_mdW, 'or', label = 'MD-Wasser',)

test = excel.loc[3]['Molenbruch Isobutan']

for MR in range(75):
    xbutan = excel.loc[MR]['Molenbruch Isobutan']
    
    if xbutan >= 0:
        
        xpropan = 1 - xbutan
        p = excel.loc[MR]['Druck Verdampfer ein'] * 1e5
        T_kmein = excel.loc[MR]['Temperatur Verdampfer ein']
        
        
        KM = f'IsoButane[{xbutan}]&n-Propane[{xpropan}]'
        Tsiede_KM = cp.PropsSI('T', 'P', p, 'Q', 0, KM) - 273.15
        Tkondens_KM = cp.PropsSI('T', 'P', p, 'Q', 1, KM) -273.15
        
        print(f'Zeile {MR +2}:')
        print(f'x-Butan: {xbutan}')
        print(f'Siedetemperatur: {Tsiede_KM} C')
        print(f'Kondensationstemperatur: {Tkondens_KM} C')
        print(f'Temp. KM-ein: {T_kmein} C')
        
        if T_kmein <  Tkondens_KM and T_kmein > Tsiede_KM:
            print('Eintritt bereits im ND-Gebiet')            
        elif T_kmein < Tsiede_KM:
            print('Eintritt als unterkühlte Flüssigkeit')
        elif T_kmein > Tkondens_KM:
            print('!!!KEINE VERDAMPFUNG!!!')
        else:
            print('Fehler in Fallunterscheidung')
            
        print('')
        
    else:
        print('')