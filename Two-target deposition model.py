# -*- coding: utf-8 -*-
#%%
import Sup_functions as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.constants as cnst

A_Ti = 0.031        #Площадь титановой мишени (м2)
A_Cr = 0.031        #Площадь хромовой мишени (м2)
A_chamber_Ti = 0.2  #Подверженность поверхности камеры воздействию потока частиц Ti и TiO2 (м2)
A_chamber_Cr = 0.2  #Подверженность поверхности камеры воздействию потока частиц Cr и Cr2O3 (м2)
S_Ti = 0.4          #Коэффициент распыления Ti для E(Ar+) = 400 эВ
S_Ti02 = 0.015      #Коэффициент распыления Ti02 для E(Ar+) = 400 эВ
S_Cr = 1.0          #Коэффициент распыления Cr для E(Ar+) = 400 эВ
S_Cr2O3 = 0.25      #Коэффициент распыления Cr2O3 для E(Ar+) = 400 эВ
S_O2 = 50/1000           #Скорость перекачки (м3/с)
alpha0_Ti = 1.0     #Коэффициент задерживания молекулы кислорода на прореагировавшей поверхности Ti
alpha0_TiO2 = 0.01  #Коэффициент задерживания молекулы кислорода на прореагировавшей поверхности TiO2
alpha0_Cr = 1.0     #Коэффициент задерживания молекулы кислорода на прореагировавшей поверхности Cr
alpha0_Cr2O3 = 0.01 #Коэффициент задерживания молекулы кислорода на прореагировавшей поверхности Cr2O3
J_Ti = 60           #Плотность ионов аргона, вызывающих распыление с поверхности мишени Ti
J_Cr = 40           #Плотность ионов аргона, вызывающих распыление с поверхности мишени Cr

df = pd.DataFrame([])

moleclar_mass_O2 = 32/cnst.Avogadro

df["P_O2"] = np.arange(start = 0.3, stop = 0, step = -0.001)
df["F"] = sp.F_of_P(df.P_O2, 273, moleclar_mass_O2)

df["tetha_t_Ti"] = sp.theta_t_of_alpha(k = 2, n = 2, alpha0_target = alpha0_Ti,
                                       alpha0_compound = alpha0_TiO2, F = df.F,
                                       J = J_Ti, S_compound = S_Ti02)

df["tetha_t_Cr"] = sp.theta_t_of_alpha(k = 2, n = 1.5, alpha0_target = alpha0_Cr,
                                       alpha0_compound = alpha0_Cr2O3, F = df.F,
                                       J = J_Cr, S_compound = S_Cr2O3)

df["tetha_c_Ti"] = sp.theta_c_of_alpha(k = 2, n = 2, tetha_t = df.tetha_t_Ti,
                                       alpha0_target = alpha0_Ti, alpha0_compound = alpha0_TiO2,
                                       F = df.F, J = J_Ti, S_target = S_Ti, S_compound = S_Ti02, A_t = A_Ti, A_c = A_chamber_Ti)

df["tetha_c_Cr"] = sp.theta_c_of_alpha(k = 2, n = 1.5, tetha_t = df.tetha_t_Cr,
                                       alpha0_target = alpha0_Cr, alpha0_compound = alpha0_Cr2O3,
                                       F = df.F, J = J_Cr, S_target = S_Cr, S_compound = S_Cr2O3, A_t = A_Cr, A_c = A_chamber_Cr)
#Поток O2 на Ti мишень (Pa*m^3/s)
df["q_t_Ti"] = sp.q_of_tetha(alpha0_target = alpha0_Ti, alpha0_compound = alpha0_TiO2, F = df.F, theta = df.tetha_t_Ti, A = A_Ti)

#Поток O2 на Cr мишень (Pa*m^3/s)
df["q_t_Cr"] = sp.q_of_tetha(alpha0_target = alpha0_Cr, alpha0_compound = alpha0_Cr2O3, F = df.F, theta = df.tetha_t_Cr, A = A_Cr)

#Поток O2 на подложкодержатель и стенки, обращенные к Ti мишени (Pa*m^3/s)
df["q_c_Ti"] = sp.q_of_tetha(alpha0_target = alpha0_Ti, alpha0_compound = alpha0_TiO2, F = df.F, theta = df.tetha_c_Ti, A = A_Ti)

#Поток O2 на подложкодержатель и стенки, обращенные к Cr мишени (Pa*m^3/s)
df["q_c_Cr"] = sp.q_of_tetha(alpha0_target = alpha0_Cr, alpha0_compound = alpha0_Cr2O3, F = df.F, theta = df.tetha_c_Cr, A = A_Cr)

#Поток O2 в насос
df["q_P"] = sp.q_of_P(P = df.P_O2, S = S_O2)

#Суммарный поток O2
df["q_O2"] = df.q_t_Ti + df.q_t_Cr + df.q_c_Ti + df.q_c_Cr + df.q_P

df["flow_rate"] = df.q_O2 * 1E3

df["R_Ti"] = sp.R_of_tetha(J = J_Ti, S_compound = S_Ti02, S_target = S_Ti, tetha_t = df.tetha_t_Ti)

df["R_Cr"] = sp.R_of_tetha(J = J_Cr, S_compound = S_Cr2O3, S_target = S_Cr, tetha_t = df.tetha_t_Cr)

df.plot(x = "flow_rate", y = "P_O2")
df.plot(x = "R_Ti", y = "P_O2")
df.plot(x = "R_Cr", y = "P_O2")

#%%
