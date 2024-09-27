import streamlit as st

import plotly.graph_objects as go
import numpy as np
from Model import *
import Target_setup as ts
import Deposition_model
import scipy.constants as cnst


@st.cache_data
def calc_error(_model_Si_Mo, P_N, P_O, J_Si_min, J_Si_max, J_Mo_min, J_Mo_max):
    P_N = P_N
    P_O = P_O
    J_Si = np.linspace(J_Si_min, J_Si_max, 25)
    J_Mo = np.linspace(J_Mo_min, J_Mo_max, 25)

    P_N_, J_Si_, J_Mo_ = np.meshgrid(P_N, J_Si, J_Mo)

    predict = _model_Si_Mo.predict(
        P_N2=P_N, P_O2=P_O, J_Si=J_Si_, J_Mo=J_Mo_)

    errors = point.loss(predict, purpose)
    return (errors, J_Si_, J_Mo_)


@st.cache_data
def calc_error_array(_model_Si_Mo, P_min, P_max, P_O, J_Si_min, J_Si_max, J_Mo_min, J_Mo_max):
    X = np.linspace(P_min, P_max, 200)
    min_errors = np.array([])
    for P in X:
        errors, _, _ = calc_error(_model_Si_Mo, P, P_O, J_Si_min,
                                  J_Si_max, J_Mo_min, J_Mo_max)
        min_errors = np.append(min_errors, np.min(errors))
    return min_errors, X


C_Si_btn = st.number_input("Si, %", value=40.0)
C_Mo_btn = st.number_input("Mo, %", value=6.0)
C_N_btn = st.number_input("N, %", value=53.0)
C_O_btn = st.number_input("O, %", value=1.0)

J_Si_min, J_Si_max = st.slider(
    "Пределы моделирования по току Si мишени, А", 0.0001, 150.0, (1.0, 10.0))
J_Mo_min, J_Mo_max = st.slider(
    "Пределы моделирования по току Mo мишени, А", 0.0001, 150.0, (1.0, 10.0))

P_N_slider = st.slider("Давление азота, Па", 0.00001, 20.0, 1.5)
P_O_input = st.number_input("Остаточное давление в камере, мкПа")
P_O_input *= 1e-6

Si = ts.target(k=2, n=4/3, t=7, A=0.004, A_chamber=0.3,
               S=0.9, S_compound=0.3,
               alpha0=0.5, alpha0_compound=0.00002,  # Азот к Si на мишени
               alpha0_c=0.001, alpha0_compound_c=0.000001,  # Азот к Si на поверхности
               alpha0_O2=0.12, alpha0_O2_compound=0.16,  # Кислород к Si
               J=165, Ro=3210, M=140.2833/1000)

Mo = ts.target(k=2, n=1, t=4, A=0.004, A_chamber=0.3,
               S=0.9, S_compound=0.3,
               alpha0=0.1, alpha0_compound=0.0001,  # Азот к молибдену на мишени
               alpha0_c=0.017, alpha0_compound_c=0.000001,  # Азот к молибдену на поверхности
               alpha0_O2=0.01, alpha0_O2_compound=0.17,  # Кислород к молибдену
               J=85, Ro=9300, M=109.95/1000)

Si_Mo = Deposition_model.deposition_model(
    target_1=Si, target_2=Mo, T=273, S=1200 * 1E-3,
    moleclar_mass=28/1000/cnst.Avogadro)

model_Si_Mo = model(Si_Mo)

purpose = point(C_Si=C_Si_btn, C_Mo=C_Mo_btn, C_N=C_N_btn, C_O2=C_O_btn)

errors, J_Si_, J_Mo_ = calc_error(model_Si_Mo, P_N_slider, P_O_input, J_Si_min,
                                  J_Si_max, J_Mo_min, J_Mo_max)
st.text("Минимальная ошибка при данном давлении " + str(np.min(errors)))

sizes = errors.flatten()

fig = go.Figure(data=[go.Scatter3d(
    x=J_Si_.flatten(),
    y=J_Mo_.flatten(),
    z=sizes,
    mode='markers',
    marker=dict(
        size=2,
        colorscale='Viridis_r',
        color=sizes,  # Цвета также могут зависеть от размера
        line=dict(width=0.5)
    )
)])

# Настраиваем оси
fig.update_layout(scene=dict(
    xaxis_title='Ток Si мишени',
    yaxis_title='Ток Mo мишени',
    zaxis_title='Ошибка'),)

# Показываем график
st.plotly_chart(fig, use_container_width=True)

if st.checkbox('Вычислить наименьшие ошибка для всех давлений'):
    P_min, P_max = st.slider(
        "Пределы моделирования по давлений азота", 0.0001, 20.0, (0.0001, 1.9))
    min_errors, X = calc_error_array(
        model_Si_Mo, P_min, P_max, P_O_input, J_Si_min, J_Si_max, J_Mo_min, J_Mo_max)
    df = pd.DataFrame([])
    df["Давление азота, Па"] = X
    print(min_errors)
    df["Минимальная ошибка, %"] = min_errors
    st.line_chart(df, x="Давление азота, Па")
