import plotly.graph_objects as go
import numpy as np
from Model import *
import Target_setup as ts
import Deposition_model
import scipy.constants as cnst

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

purpose = point(C_Si=40, C_Mo=6, C_N=53, C_O2=1)

P_N = 1.5
P_O = 0.000001
J_Si = np.arange(0.001, 2, 0.01)
J_Mo = np.arange(0.001, 2, 0.01)

P_N_, J_Si_, J_Mo_ = np.meshgrid(P_N, J_Si, J_Mo)

errors = point.loss(model_Si_Mo.predict(
    P_N2=P_N, P_O2=P_O, J_Si=J_Si_, J_Mo=J_Mo_), purpose)


sizes = errors.flatten() / np.max(errors) * 10

fig = go.Figure(data=[go.Scatter3d(
    x=J_Si_.flatten(),
    y=J_Mo_.flatten(),
    z=sizes,
    mode='markers',
    marker=dict(
        size=2,
        color=sizes,  # Цвета также могут зависеть от размера
        line=dict(width=0.5)
    )
)])

# Настраиваем оси
fig.update_layout(scene=dict(
    xaxis_title='X Axis',
    yaxis_title='Y Axis',
    zaxis_title='Z Axis'),)

# Показываем график
fig.show()
