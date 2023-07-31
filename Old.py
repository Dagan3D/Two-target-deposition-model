# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 17:56:19 2023

@author: butma
"""


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

df["q_t_c"] = df.q_t_Ti + df.q_t_Cr + df.q_c_Ti + df.q_c_Cr

#Поток O2 в насос
df["q_P"] = sp.q_of_P(P = df.P_O2, S = S_O2)

#Суммарный поток O2
df["q_O2"] = df.q_t_Ti + df.q_t_Cr + df.q_c_Ti + df.q_c_Cr + df.q_P

df["flow_rate"] = df.q_O2 / df.P_O2 * 1E3

df["R_Ti"] = sp.R_of_tetha(J = J_Ti, S_compound = S_Ti02, S_target = S_Ti, tetha_t = df.tetha_t_Ti)

df["R_Cr"] = sp.R_of_tetha(J = J_Cr, S_compound = S_Cr2O3, S_target = S_Cr, tetha_t = df.tetha_t_Cr)

df.plot(x = "flow_rate", y = "P_O2")
# df.plot(x = "R_Ti")
# df.plot(x = "R_Cr")
df[["P_O2", "R_Ti", "R_Cr"]].plot(x = "P_O2")

#%%



df["t_Ti"] = sp.t_of_F(k = 2, n = 2, J = J_Ti, alpha0_target = alpha0_Ti, F = df.F, S_compound = S_Ti02)

df["t_Cr"] = sp.t_of_F(k = 4, n = 3, J = J_Cr, alpha0_target = alpha0_Cr, F = df.F, S_compound = S_Cr2O3)

df["c_Ti"] = sp.c_of_F(k = 2, n = 2, alpha0_target = alpha0_Ti,
                       alpha0_compound = alpha0_TiO2, F = df.F, J = J_Ti, S_target = S_Ti,
                       S_compound = S_Ti02, A_t = A_Ti, A_c = A_chamber_Ti)

df["c_Cr"] = sp.c_of_F(k = 4, n = 3, alpha0_target = alpha0_Cr,
                       alpha0_compound = alpha0_Cr2O3, F = df.F, J = J_Cr, S_target = S_Cr,
                       S_compound = S_Cr2O3, A_t = A_Cr, A_c = A_chamber_Cr)

K = sp.K_calc(T, M0 = moleclar_mass_O2)

df["q_O2"] = sp.q_of_F_t_c(F = df.F, t_1 = df.t_Ti, t_2 = df.t_Cr, c_1 = df.c_Ti, c_2 = df.c_Cr,
                           A_t1 = A_Ti, A_t2 = A_Cr, A_c1 = A_chamber_Ti, A_c2 = A_chamber_Cr,
                           alpha0_target1 = alpha0_Ti, alpha0_target2 = alpha0_Cr, S_a = S_O2*K)

df[["P_O2", "t_Ti", "t_Cr", "c_Ti", "c_Cr"]].plot(x = "P_O2")
df["flow_rate"] = df.q_O2 * 1E3
df.plot(x = "flow_rate", y = "P_O2")


df["dq/dP"] = sp.dq_dp(S = S_O2, K = K,
                        t_1 = df.t_Ti, t_2 = df.t_Cr,
                        c_1 = df.c_Ti, c_2 = df.c_Cr, A_t1 = A_Ti, A_t2 = A_Cr,
                        A_c1 = A_chamber_Ti, A_c2 = A_chamber_Cr,
                        S_target_1 = S_Ti, S_compound_1 = S_Ti02,
                        S_target_2 = S_Cr, S_compound_2 = S_Cr2O3,
                        alpha0_target1 = alpha0_Ti, alpha0_target2 = alpha0_Cr)

df["S_l"] = (S_O2 - df["dq/dP"])*1000
df["q_O2"] = df["dq/dP"].cumsum() * 1e-5