# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 09:35:20 2023

@author: Philip
"""

import pandas as pd

csvdata= pd.read_csv("E:\\OneDrive - Universitaet Duisburg-Essen\\Thermo Praktikum\\Daten Wärmepumpe\\Alex Daten Wärmepumpe\\Simon Schmeing\\Messdaten\\Messung 20\\3 - Kopie.csv", header = None)

names = ['Zeit', 'Temp 9', 'Druck 4', 'Temp 7', 'Druck 1', 'Temp 1', 'Temp 2', 'Temp 17', 'Temp 18', 'Temp 19', 'Temp 20', 'Temp 13', 'Temp 14', 'Temp 15', 'Temp 16', 'Temp 6', 'Druck 3', 'Temp 10', 'Druck 2', 'Temp 3', 'Temp 5', 'Temp 11', 'Temp 8', 'Temp 12', 'Temp 4', 'Leistung', 'Druck 6', 'Druck 5', 'Verd.aus Temp', 'Temp 33', 'Temp 21', 'Temp 26', 'Temp 28', 'Temp 27', 'Temp 29', 'Temp 25', 'Temp 22', 'Temp 23', 'Temp 24', 'Temp 30', 'Temp 31', 'Temp 32', 'Temp 35', 'Temp 34']

excel = pd.read_excel("E:\\OneDrive - Universitaet Duisburg-Essen\\Thermo Praktikum\\Daten Wärmepumpe\\Alex Daten Wärmepumpe\\Messdaten Propan und Isobutan.xlsx", skiprows=2)

sheetname="Tabelle1"