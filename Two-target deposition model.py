# -*- coding: utf-8 -*-
#%% Настройка питона
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.constants as cnst
import matplotlib.pylab as mpl
import target_setup as ts
import Model

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.size'] = 10
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
plt.style.use('seaborn-whitegrid')
mpl.rc('font', size=13)
mpl.rcParams['font.family'] = "serif"
mpl.rcParams['axes.linewidth'] = 2

#%% Переменные
# A_Ti = 0.031        #Площадь титановой мишени (м2)
# A_Cr = 0.031        #Площадь хромовой мишени (м2)
# A_chamber_Ti = 0.2  #Подверженность поверхности камеры воздействию потока частиц Ti и TiO2 (м2)
# A_chamber_Cr = 0.2  #Подверженность поверхности камеры воздействию потока частиц Cr и Cr2O3 (м2)
# S_Ti = 0.4          #Коэффициент распыления Ti для E(Ar+) = 400 эВ
# S_Ti02 = 0.015      #Коэффициент распыления Ti02 для E(Ar+) = 400 эВ
# S_Cr = 1.0          #Коэффициент распыления Cr для E(Ar+) = 400 эВ
# S_Cr2O3 = 0.25      #Коэффициент распыления Cr2O3 для E(Ar+) = 400 эВ
# S_O2 = 50/1000      #Скорость перекачки (м3/с)
# alpha0_Ti = 1.0     #Коэффициент задерживания молекулы кислорода на прореагировавшей поверхности Ti
# alpha0_TiO2 = 0.01  #Коэффициент задерживания молекулы кислорода на прореагировавшей поверхности TiO2
# alpha0_Cr = 1.0     #Коэффициент задерживания молекулы кислорода на прореагировавшей поверхности Cr
# alpha0_Cr2O3 = 0.01 #Коэффициент задерживания молекулы кислорода на прореагировавшей поверхности Cr2O3
# J_Ti = 60           #Плотность ионов аргона, вызывающих распыление с поверхности мишени Ti
# J_Cr = 40           #Плотность ионов аргона, вызывающихщих распыление с поверхности мишени Cr
#T = 300

#%% Расчёт

df = pd.DataFrame([])

df["P_O2"] = np.arange(start = 0, stop = 0.32, step = 0.00001)

Ti = ts.target(k = 2, n = 2, A = 0.031, A_chamber = 0.2,
            S = 0.4, S_compound = 0.015,
            alpha0 = 1.0, alpha0_compound = 0.01, J = 60)

Cr = ts.target(k = 4, n = 3, A = 0.031, A_chamber = 0.2,
            S = 1.0, S_compound = 0.25,
            alpha0 = 1.0, alpha0_compound = 0.01, J = 40)

Ti_Cr = Model.model(target_1 = Ti, target_2 = Cr, T = 300, S = 50/1000, moleclar_mass = 32/1000/cnst.Avogadro)

fl_rate = pd.DataFrame()
fl_rate["P_O2"] = df.P_O2

Sl = pd.DataFrame()
Sl["P_O2"] = df.P_O2

temp_range = [10, 100, 200, 273, 300, 1000, 2000, 3000]

for t in temp_range:
    Ti_Cr.T = t
    df["q_O2"] = Ti_Cr.q_of_F_t_c(df.P_O2)
    df["flow_rate"] = df.q_O2 * 1E3
    df["dq/dP"] = Ti_Cr.dq_dp(df.P_O2)
    df["S_l"] = (Ti_Cr.S - df["dq/dP"])*1000
    fl_rate[str(t) + " °K"] = df.flow_rate
    Sl[str(t) + " °K"] = df.S_l


#%% Отображение
df[["P_O2", "S_l"]].plot(x = "P_O2", xlim = [0.0001, 0.05],  ylim = [0, 950], legend = None, ylabel = r"Characteristic function $S_L$ $(L*s^{-1})$", xlabel = "Oxygen partial pressure (Pa)")
df.plot(x = "flow_rate", y = "P_O2", legend = None, ylabel = r"$P_{O_2}$ (Pa)", xlabel = "Oxygen flow rate (sccm)")


for t in temp_range:
    plt.plot(fl_rate[str(t) + " °K"], fl_rate["P_O2"], label = str(t) + " °K")
plt.legend()
plt.ylabel = r"Characteristic function $S_L$ $(L*s^{-1})$"
plt.xlabel = "Oxygen partial pressure (Pa)"
plt.show()
    

for t in temp_range:
    plt.plot(Sl["P_O2"], Sl[str(t) + " °K"], label = str(t) + " °K")
plt.legend()
plt.xlim(-0.00001, 0.05)
plt.ylim(0, 1200)
plt.ylabel = r"Characteristic function $S_L$ $(L*s^{-1})$"
plt.xlabel = "Oxygen partial pressure (Pa)"
plt.show()