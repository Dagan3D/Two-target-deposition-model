# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 10:25:43 2023

@author: butma
"""

from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pylab as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.size'] = 10
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
plt.style.use('seaborn-whitegrid')
mpl.rc('font', size=13)
mpl.rcParams['font.family'] = "serif"
mpl.rcParams['axes.linewidth'] = 2

input_df =  pd.read_excel(".\Data\Spray_rate.xls")
input_df["P"] = input_df["P"] * 133.322368
#input_df["P"] = pd.to_numeric(input_df['P'], errors='coerce')
df = pd.DataFrame()
df["P_O2"] = np.arange(start = 0.0, stop = 0.003 * 133.322368, step = 0.01)
Si_interpolation = interpolate.PchipInterpolator(input_df["P"], input_df["Si"])
Mo_interpolation = interpolate.interp1d(input_df["P"], input_df["Mo"], kind="quadratic")

df["k_S_Si"] = Si_interpolation(df["P_O2"]) 
df["k_S_Mo"] = Mo_interpolation(df["P_O2"])
 
plt.plot(input_df["P"], input_df["Si"], label = "Experiment")
plt.plot(df["P_O2"], df["k_S_Si"], label = "Interpolate")
plt.legend()
plt.title("Изменение интенсивности распыения Si")
plt.show()
plt.plot(input_df["P"], input_df["Mo"], label = "Experiment")
plt.plot(df["P_O2"], df["k_S_Mo"], label = "Interpolate")
plt.legend()
plt.title("Изменение интенсивности распыения Mo")
plt.show()