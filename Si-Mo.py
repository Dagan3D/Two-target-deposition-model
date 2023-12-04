# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 15:ЭЛГЕАН27:28 2023

@author: butmahh``
"""

# -*- coding: utf-8 -*-
# %% Настройка питона
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.constants as cnst
import matplotlib.pylab as mpl
from scipy import interpolate
import target_setup as ts
import Model
import matplotlib.colors as colors

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.size'] = 10
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
plt.style.use('seaborn-whitegrid')
mpl.rc('font', size=13)
mpl.rcParams['font.family'] = "serif"
mpl.rcParams['axes.linewidth'] = 2

# %% Переменные
A_Si = 0.004  # Площадь кремниевой мишени (м2)
A_Mo = 0.004  # Площадь молибденовой мишени (м2)
# Подверженная воздействию потока частиц площадь поверхности камеры  Si и Si3N4 (м2)
A_chamber_Si = 0.3
# Подверженная воздействию потока частиц площадь поверхности камеры Mo и MoN (м2)
A_chamber_Mo = 0.3
S_Si = 0.9    # Коэффициент распыления Si для E(Ar+) = 400 эВ
S_Si3N4 = 0.3  # Коэффициент распыления SiN для E(Ar+) = 400 эВ
S_Mo = 0.9  # Коэффициент распыления Mo для E(Ar+) = 400 эВ
S_MoN = 0.3 # Коэффициент распыления MoN для E(Ar+) = 400 эВ
S_O2 = 1200/1000  # Скорость перекачки (м3/с)
alpha0_Si = 1  # Коэффициент задерживания молекулы кислорода на непрореагировавшей поверхности Si
# Коэффициент задерживания молекулы кислорода на прореагировавшей поверхности Si3N4
alpha0_Si3N4 = 0.01
# Коэффициент задерживания молекулы кислорода на прореагировавшей поверхности Mo
alpha0_Mo = 0.1
# Коэффициент задерживания молекулы кислорода на прореагировавшей поверхности MoN
alpha0_MoN = 0.01
# Плотность потока ионов аргона, вызывающих распыление с поверхности мишени Si (А*м^-2)
J_Si = 165
# Плотность потока ионов аргона, вызывающихщих распыление с поверхности мишени Mo (А*м^-2)
J_Mo = 85
T = 273  # Температура в К

# %% Определение пределов моделирования
df = pd.DataFrame([])
df["P_O2"] = np.arange(start=0.0, stop=0.003 * 133.322368, step=0.00001)

# %% Интерполяция силы тока из экспериментаьных данных по Si
x_Si = np.array([0, 0.0001, 0.0003, 0.0006, 0.001, 0.0015,
                0.002, 0.0025, 0.003]) * 133.322368
y_Si = np.array([0.66, 0.69, 0.72, 0.79, 0.95, 1.0, 1.05, 1.1, 1.14]) / A_Si
spl = interpolate.UnivariateSpline(x_Si, y_Si, k=2)
df["J_Si"] = spl(df.P_O2)

# %% Интерполяция силы тока из экспериментаьных данных по Mo
x_Mo = np.array([0.0001, 0.0003, 0.0006, 0.001, 0.003]) * 133.322368
y_Mo = np.array([0.27, 0.27, 0.27, 0.26, 0.25]) / A_Mo

linear_model = np.polyfit(x_Mo, y_Mo, 1)
linear_model_fn = np.poly1d(linear_model)

df["J_Mo"] = linear_model_fn(df.P_O2)

# %% Интерполяция изменения эффективности распыленя мишеней
input_df = pd.read_excel(".\Data\Spray_rate.xls")
input_df["P_O2"] = input_df["P"] * 133.322368
input_df["P"] *= 1000

Si_interpolation = interpolate.PchipInterpolator(
    input_df["P_O2"], input_df["Si"])
Mo_interpolation = interpolate.interp1d(
    input_df["P_O2"], input_df["Mo"], kind="quadratic")
df["k_S_Si"] = Si_interpolation(df["P_O2"])
df["k_S_Mo"] = Mo_interpolation(df["P_O2"])

df["k_S_Si_orig"] = Si_interpolation(df["P_O2"])
df["k_S_Mo_orig"] = Mo_interpolation(df["P_O2"])

input_df["Si"] /= input_df["Si"][0]
input_df["Mo"] /= input_df["Mo"][0]
df["k_S_Mo"] /= df["k_S_Mo"][0]
df["k_S_Si"] /= df["k_S_Si"][0]

# %% Инициализация мишеней
Si = ts.target(k=2, n=4/3, t = 7, A=A_Si, A_chamber=A_chamber_Si,
               S=S_Si * df.k_S_Si, S_compound= S_Si3N4 * df.k_S_Si,
               alpha0=0.5, alpha0_compound=0.00002,         #Азот к Si на мишени
               alpha0_c=0.001, alpha0_compound_c=0.000001,   #Азот к Si на поверхности
               alpha0_O2 = 0.12, alpha0_O2_compound=0.16,     #Кислород к Si
               J=J_Si, Ro=3210, M = 140.2833/1000)

Mo = ts.target(k=2, n=1, t = 4, A=A_Mo, A_chamber=A_chamber_Mo,
               S=S_Mo * df.k_S_Si, S_compound=S_MoN * df.k_S_Si,
               alpha0=0.1, alpha0_compound=0.0001,       #Азот к молибдену на мишени
               alpha0_c=0.017, alpha0_compound_c=0.000001, #Азот к молибдену на поверхности
               alpha0_O2 = 0.01, alpha0_O2_compound=0.17,    #Кислород к молибдену
               J=J_Mo, Ro=9300, M = 109.95/1000)

Si_Mo = Model.model(target_1=Si, target_2=Mo, T=273, S = 1200 * 1E-3,
                    moleclar_mass=28/1000/cnst.Avogadro)

# %% Расчёт характеристической функции и потока
Sl = pd.DataFrame()
Sl["P_O2"] = df.P_O2
df["dq/dP"] = Si_Mo.dq_dp(df.P_O2)
df["S_l"] = (Si_Mo.S - df["dq/dP"])*1000

df["q_O2"] = Si_Mo.q_of_F_t_c(df.P_O2)
df["flow_rate"] = df.q_O2 * 1E3

# %% Расчёт коэфициента распыления
df["F"] = Si_Mo.F_of_P(df.P_O2)
df["R_Si"] = Si.R_of_F(Si.N_target_of_F(df.F), in_nm_s=True)
df["R_Mo"] = Mo.R_of_F(Mo.N_target_of_F(df.F), in_nm_s=True)

# %% Вычисление количества веществ в образце
df["N_Si"]= Si.N_target_of_F(df.F)

df["N_Mo"] = Mo.N_target_of_F(df.F)
df["N_O2"] = Si_Mo.N_O2_of_P(df.P_O2, P_residual=7e-5*133)
df["N_N"] = Si_Mo.N_gase_of_P(df.P_O2)[-1]
#n_N = pd.DataFrame()
#n_N["a_1"], n_N["a_2"], n_N["b_1"], n_N["b_2"], n_N["c_1"], n_N["c_2"], df["N_N"] = Si_Mo.N_gase_of_P(df.P_O2)
#n_N["sum_1"] = n_N["a_1"] + n_N["b_1"] + n_N["c_1"]
#n_N["sum_2"] = n_N["a_2"] + n_N["b_2"] + n_N["c_2"]

# %% Расчёт соотношения элементов
df["N"] = df.N_Si + df.N_Mo + df.N_O2 + df.N_N

df["C_Si"] = df.N_Si / df.N * 100
df["C_Mo"] = df.N_Mo / df.N * 100
df["C_N"] = df.N_N / df.N * 100
df["C_O2"] = df.N_O2 / df.N * 100

df["P_O2_torr"] = df.P_O2 * 0.0075 * 1000
df["P_O2_from_q"] = df.q_O2/Si_Mo.S * 0.0075 * 1000

# %% Зарузка эксперементаьных значений концентрации
exp_C = pd.read_excel("./Data/Atomic_content.xlsx")
exp_R = pd.read_excel("./Data/Spray_rate_original.xlsx")

# %% Настройи отображения
colors_list = list(colors.TABLEAU_COLORS.values())

#%% Вывод тока на кремниевой мишени
plt.plot(x_Si, y_Si)
plt.plot(df.P_O2, df.J_Si)
plt.legend(["Experiment Si", "Model"])
plt.show()

#%% Вывод тока на кремниевой мишени 
plt.plot(x_Mo, y_Mo)
plt.plot(df.P_O2, df.J_Mo)
plt.legend(["Experiment Mo", "Model"])
plt.show()

#%% Вывод скорости осаждения кремния и молибдена
plt.scatter(exp_R["P"], y=exp_R.Si, c=colors_list[0])
plt.scatter(exp_R["P"], y=exp_R.Mo, c=colors_list[1])
df.R_Si /= 600
df.R_Mo /= 42
plt.plot(df["P_O2_torr"], df["R_Si"], label="Si", c=colors_list[0])
plt.plot(df["P_O2_torr"], df["R_Mo"], label="Mo", c=colors_list[1])

plt.title("Скорость осаждения Si и Mo")
plt.xlabel("Давление, мТорр")
plt.ylabel("Скорость распыления, нм/c")
plt.legend()
plt.show()

#%% ЧБ Вывод скорости осаждения кремния и молибдена
plt.scatter(exp_R["P"], y=exp_R.Si, c="black", linestyle = "-")
plt.scatter(exp_R["P"], y=exp_R.Mo, c="black", linestyle = "--")
plt.plot(df["P_O2_torr"], df["R_Si"], label="Si", c=colors_list[0])
plt.plot(df["P_O2_torr"], df["R_Mo"], label="Mo", c=colors_list[1])

plt.title("Скорость осаждения Si и Mo")
plt.xlabel("Давление, мТорр")
plt.ylabel("Скорость распыления, нм/c")
plt.legend()
plt.show()

#%% Вывод коэфициентов для коэфициента распыления 
plt.scatter(exp_R["P"], y=exp_R.Si, c=colors_list[0])
plt.scatter(exp_R["P"], y=exp_R.Mo, c=colors_list[1])
plt.plot(df["P_O2_torr"], df.k_S_Si_orig, label="Si", c=colors_list[0])
plt.plot(df["P_O2_torr"], df.k_S_Mo_orig, label="Mo", c=colors_list[1])

plt.title("Скорость осаждения Si и Mo")
plt.xlabel("Давление, мТорр")
plt.ylabel("Скорость распыления, нм/c")
plt.legend()
plt.show()

#%% Вывод эксериментальных значений скорости распыления
plt.scatter(exp_R["P"], y=exp_R.Si, c="black", marker = '*')
plt.scatter(exp_R["P"], y=exp_R.Mo, c="black", marker = 'o')
plt.plot(df["P_O2_torr"], df.k_S_Si_orig, label="Si", c="black",  linestyle = "-",)
plt.plot(df["P_O2_torr"], df.k_S_Mo_orig, label="Mo", c="black",  linestyle = "--",)

plt.title("Скорость осаждения Si и Mo")
plt.xlabel("Давление, мТорр")
plt.ylabel("Скорость распыления, нм/c")
plt.legend()
plt.savefig("./image/ЧБ Скорость осаждения Si и Mo.tiff", format = "tiff", dpi = 300)
plt.show()

#%% Вывод расчитаных зменение интенсивности распыления Si
# plt.scatter(input_df["P"], input_df["Si"],
#             label="Experiment", c=colors_list[5])
plt.plot(df["P_O2_torr"], df["k_S_Si"], label=r"$S_{Si}$", c=colors_list[0])
plt.plot(df["P_O2_torr"], df["k_S_Mo"], label=r"$S_{Mo}$", c=colors_list[1])
plt.legend()
plt.title("Множитель при коэффициенте распыления")
plt.xlabel("Давление, мТорр")
plt.ylabel("Множитель")
plt.show()

#%% Вывод расчитаных зменение интенсивности распыления Mo

plt.scatter(input_df["P"], input_df["Mo"],
            label="Experiment", c=colors_list[5])
plt.plot(df["P_O2_torr"], df["k_S_Mo"], label="Interpolate", c=colors_list[1])
plt.legend()
plt.title("Изменение интенсивности распыления Mo")
plt.xlabel("Давление, мТорр")
plt.ylabel("Эффективность распыления Mo")
plt.show()

#%% Вывод харатеристичской функции и потока азота в зависимости от давения азота
df[["P_O2", "S_l"]].plot(x = "P_O2",  legend = None, ylabel = r"$S_L$, $Л*с^{-1}$", xlabel = "Парциальное давление азота, Па", )
df.plot(x = "flow_rate", y = "P_O2", legend = None, ylabel = r"$P_{N_2}$, Па", xlabel = "Расход азота, sccm")

#%% Вывод количество осаденных атосмов разного сорта
df[["P_O2_torr", "N_Si", "N_N", "N_O2", "N_Mo"]].plot(x="P_O2_torr", xlabel = "Парциальное давление азота, Па", ylabel = "Cкорость осаждения, $см^{-2}с^{-1}$") 
plt.show()

df[["P_O2_torr", "N_Si", "N_Mo"]].plot(x="P_O2_torr", xlabel = "Парциальное давление азота, Па", ylabel = "Cкорость осаждения, $см^{-2}с^{-1}$") 
df[["P_O2_torr", "N_N"]].plot(x="P_O2_torr", xlabel = "Парциальное давление азота, Па", ylabel = "Cкорость осаждения, $см^{-2}с^{-1}$") 
df[["P_O2_torr", "N_O2"]].plot(x="P_O2_torr", xlabel = "Парциальное давление азота, Па", ylabel = "Cкорость осаждения, $см^{-2}с^{-1}$") 
plt.show()
#%% Цветная версия вывода расёта содержания элементов

plt.plot(df["P_O2_torr"], df.C_Si, label="Si", c=colors_list[0])
plt.plot(df["P_O2_torr"], df.C_Mo, label="Mo", c=colors_list[1])
plt.plot(df["P_O2_torr"], df.C_N,  label="N", c=colors_list[2])
plt.plot(df["P_O2_torr"], df.C_O2, label="O2", c=colors_list[3])

plt.scatter(exp_C["мТорр"], y=exp_C.Si, c=colors_list[0])
plt.scatter(exp_C["мТорр"], y=exp_C.Mo, c=colors_list[1])
plt.scatter(exp_C["мТорр"], y=exp_C.N,  c=colors_list[2])
plt.scatter(exp_C["мТорр"], y=exp_C.O,  c=colors_list[3])

plt.xlabel("Париальное давление азота, мТорр")
plt.ylabel("Содержание элемента, Ат %")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncols=4, mode="expand", borderaxespad=0.)

plt.show()

#%% ЧБ версия вывода расёта содержания элементов

plt.plot(df["P_O2_torr"], df.C_Si, label="Si", c="black", linestyle = "-")
plt.scatter(exp_C["мТорр"], y=exp_C.Si, label="Si", c="black", marker = '*')

plt.plot(df["P_O2_torr"], df.C_Mo, label="Mo", c="black", linestyle = "--")
plt.scatter(exp_C["мТорр"], y=exp_C.Mo, label="Mo", c="black", marker = 'o')

plt.plot(df["P_O2_torr"], df.C_N,  label="N", c="black", linestyle = "-.")
plt.scatter(exp_C["мТорр"], y=exp_C.N,  label="N", c="black", marker = 'v')

plt.plot(df["P_O2_torr"], df.C_O2, label="O2", c="black", linestyle = ":")
plt.scatter(exp_C["мТорр"], y=exp_C.O,  label="O2", c="black", marker = 'h')

plt.xlabel("Давление, мТорр")
plt.ylabel("Содержание элемента, Ат %")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncols=4, mode="expand", borderaxespad=0.)
plt.savefig("./image/ЧБ Расчёт содержания элекментов.tiff", format = "tiff", dpi = 300)
plt.show()

#%% Экспериментальное содержание элекментов ЧБ
plt.plot(exp_C["мТорр"], exp_C.Si, label="Si", c="black", marker = '*', linestyle = "-")
plt.plot(exp_C["мТорр"], exp_C.Mo, label="Mo", c="black", marker = 'o', linestyle = "--")
plt.plot(exp_C["мТорр"], exp_C.N,  label="N", c="black", marker = 'v', linestyle = "-.")
plt.plot(exp_C["мТорр"], exp_C.O,  label="O2", c="black", marker = 'h', linestyle = ":")

plt.xlabel("Давление, мТорр")
plt.ylabel("Содержание элемента, Ат %")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncols=4, mode="expand", borderaxespad=0.)
plt.savefig("./image/ЧБ Экспериментальное содержание элекментов.tiff", format = "tiff", dpi = 300)
plt.show()

#%% Перерасчёт с разными значениями остаточного авления кислорода
s = ["10 ", "30 ", "70 "]
for i, P_O2 in enumerate([1e-5*133, 3e-5*133, 7e-5 * 133]):
        
    df["N_Si"]= Si.N_target_of_F(df.F)
    
    df["N_Mo"] = Mo.N_target_of_F(df.F)
    df["N_O2"] = Si_Mo.N_O2_of_P(df.P_O2, P_residual=P_O2)
    n_N = pd.DataFrame()
    n_N["a_1"], n_N["a_2"], n_N["b_1"], n_N["b_2"], n_N["c_1"], n_N["c_2"], df["N_N"] = Si_Mo.N_gase_of_P(
        df.P_O2)
    n_N["sum_1"] = n_N["a_1"] + n_N["b_1"] + n_N["c_1"]
    n_N["sum_2"] = n_N["a_2"] + n_N["b_2"] + n_N["c_2"]
    
    df["N"] = df.N_Si + df.N_Mo + df.N_O2 + df.N_N
    
    df["C_Si"] = df.N_Si / df.N * 100
    df["C_Mo"] = df.N_Mo / df.N * 100
    df["C_N"] = df.N_N / df.N * 100
    df["C_O2"] = df.N_O2 / df.N * 100
    
    df["P_O2_torr"] = df.P_O2 * 0.0075 * 1000
    df["P_O2_from_q"] = df.q_O2/Si_Mo.S * 0.0075 * 1000
    
    if i == 0:   
        plt.plot(df["P_O2_torr"], df.C_Si, label="Si", alpha = 0.3, c="black", linestyle = "-")    
        plt.plot(df["P_O2_torr"], df.C_Mo, label="Mo", alpha = 0.3, c="black", linestyle = "--")
        plt.plot(df["P_O2_torr"], df.C_N,  label="N", alpha = 0.3, c="black", linestyle = "-.")
        plt.plot(df["P_O2_torr"], df.C_O2, label="O2", c="black", linestyle = ":")
    else:
        plt.plot(df["P_O2_torr"], df.C_Si, alpha = 0.3, c="black", linestyle = "-")    
        plt.plot(df["P_O2_torr"], df.C_Mo, alpha = 0.3, c="black", linestyle = "--")
        plt.plot(df["P_O2_torr"], df.C_N,  alpha = 0.3, c="black", linestyle = "-.")
        plt.plot(df["P_O2_torr"], df.C_O2, c="black", linestyle = ":")
    
    plt.text(3.03, df.C_Si[39996] - 0.7, s[i], size = 7)
    plt.text(3.03, df.C_Mo[39996] - 0.7, s[i], size = 7) 
    plt.text(3.03, df.C_N[39996] - 0.7, s[i], size = 7)     

plt.xlabel("Давление, мТорр")
plt.ylabel("Содержание элемента, Ат %")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncols=4, mode="expand", borderaxespad=0.)
plt.savefig("./image/ЧБ Экспериментальное содержание элекментов", format = "tiff", dpi = 300)

plt.text(2, 3, s[0] +" мкТорр")
plt.text(0.7, 19, s[1] + " мкТорр")
plt.text(1, 35, s[2] + " мкТорр")

plt.xlim(-0.01, 3.3)
plt.ylim(-0.01, 60)
plt.savefig("./image/ЧБ Разные остаточные давления кислорода.tiff", format = "tiff", dpi = 300)
plt.show()    
    
#%% Перерасчёт с разными значениями тока кремниевой мишени
for i, J_Si in enumerate([50, 120, 200]):
    Si.J = J_Si
    
    df["N_Si"]= Si.N_target_of_F(df.F)
    
    df["N_Mo"] = Mo.N_target_of_F(df.F)
    df["N_O2"] = Si_Mo.N_O2_of_P(df.P_O2, P_residual=0.933e-2)
    df["N_N"] = Si_Mo.N_gase_of_P(df.P_O2)[-1]    
    df["N"] = df.N_Si + df.N_Mo + df.N_O2 + df.N_N
    
    df["C_Si"] = df.N_Si / df.N * 100
    df["C_Mo"] = df.N_Mo / df.N * 100
    df["C_N"] = df.N_N / df.N * 100
    df["C_O2"] = df.N_O2 / df.N * 100
    
    df["P_O2_torr"] = df.P_O2 * 0.0075 * 1000
    
    if i == 0:    
        plt.plot(df["P_O2_torr"], df.C_Si, label="Si",  c="black", linestyle = "-")    
        plt.plot(df["P_O2_torr"], df.C_Mo, label="Mo", alpha = 0.3, c="black", linestyle = "--")
        plt.plot(df["P_O2_torr"], df.C_N,  label="N", alpha = 0.3, c="black", linestyle = "-.")
        plt.plot(df["P_O2_torr"], df.C_O2, label="O2", alpha = 0.3, c="black", linestyle = ":")
    else:
        plt.plot(df["P_O2_torr"], df.C_Si, c="black", linestyle = "-")    
        plt.plot(df["P_O2_torr"], df.C_Mo, alpha = 0.3, c="black", linestyle = "--")
        plt.plot(df["P_O2_torr"], df.C_N,  alpha = 0.3, c="black", linestyle = "-.")
        plt.plot(df["P_O2_torr"], df.C_O2, alpha = 0.3, c="black", linestyle = ":")
    
    plt.text(3.03, df.C_Mo[39996] - 0.7, str(J_Si), size = 7)
    plt.text(3.03, df.C_O2[39996] - 0.7, str(J_Si), size = 7) 
    plt.text(3.03, df.C_N[39996] - 0.7, str(J_Si), size = 7)   
    
    
plt.xlabel("Давление, мТорр")
plt.ylabel("Содержание элемента, Ат %")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncols=4, mode="expand", borderaxespad=0.)
plt.savefig("./image/ЧБ Экспериментальное содержание элекментов", format = "tiff", dpi = 300)
plt.text(2, 6, r"$J_{Si}$ = 50 $A/м^{2}$")
plt.text(2.1, 29, r"$J_{Si}$ = 120 $A/м^{2}$")
plt.text(2.1, 38, r"$J_{Si}$ = 200 $A/м^{2}$")
plt.xlim(-0.01, 3.3)
plt.ylim(-0.01, 60)
plt.savefig("./image/ЧБ Разные токи Si.tiff", format = "tiff", dpi = 300)
plt.show()   

#%% Перерасчёт с разными значениями тока молибденой мишени
Si.J = 165
for i, J_Mo in enumerate([50, 120, 200]):
    Mo.J = J_Mo
    
    df["N_Si"]= Si.N_target_of_F(df.F)
    
    df["N_Mo"] = Mo.N_target_of_F(df.F)
    df["N_O2"] = Si_Mo.N_O2_of_P(df.P_O2, P_residual=0.933e-2)
    df["N_N"] = Si_Mo.N_gase_of_P(df.P_O2)[-1]    
    df["N"] = df.N_Si + df.N_Mo + df.N_O2 + df.N_N
    
    df["C_Si"] = df.N_Si / df.N * 100
    df["C_Mo"] = df.N_Mo / df.N * 100
    df["C_N"] = df.N_N / df.N * 100
    df["C_O2"] = df.N_O2 / df.N * 100
    
    df["P_O2_torr"] = df.P_O2 * 0.0075 * 1000
    
    if i == 0:    
        plt.plot(df["P_O2_torr"], df.C_Si, label="Si", alpha = 0.3, c="black", linestyle = "-")    
        plt.plot(df["P_O2_torr"], df.C_Mo, label="Mo", c="black", linestyle = "--")
        plt.plot(df["P_O2_torr"], df.C_N,  label="N", alpha = 0.3, c="black", linestyle = "-.")
        plt.plot(df["P_O2_torr"], df.C_O2, label="O2", alpha = 0.3, c="black", linestyle = ":")
    else:
        plt.plot(df["P_O2_torr"], df.C_Si, alpha = 0.3, c="black", linestyle = "-")    
        plt.plot(df["P_O2_torr"], df.C_Mo, c="black", linestyle = "--")
        plt.plot(df["P_O2_torr"], df.C_N, alpha = 0.3, c="black", linestyle = "-.")
        plt.plot(df["P_O2_torr"], df.C_O2, alpha = 0.3, c="black", linestyle = ":")
    
    plt.text(3.03, df.C_Si[39996] - 0.7, str(J_Mo), size = 7)
    plt.text(3.03, df.C_O2[39996] - 0.7, str(J_Mo), size = 7) 
    plt.text(3.03, df.C_N[39996] - 0.7, str(J_Mo), size = 7)
    
plt.xlabel("Давление, мТорр")
plt.ylabel("Содержание элемента, Ат %")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncols=4, mode="expand", borderaxespad=0.)
plt.savefig("./image/ЧБ Экспериментальное содержание элекментов", format = "tiff", dpi = 300)
plt.text(0.5, 3, r"$J_{Mo}$ = 50 $A/м^{2}$")
plt.text(0.6, 9.5, r"$J_{Mo}$ = 120 $A/м^{2}$")
plt.text(0.8, 16, r"$J_{Mo}$ = 200 $A/м^{2}$")
plt.xlim(-0.01, 3.3)
plt.ylim(-0.01, 60)
plt.savefig("./image/ЧБ Разные токи Mo.tiff", format = "tiff", dpi = 300)
plt.show()     

