# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 13:43:39 2023

@author: butma
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import UnivariateSpline



x = np.array([0, 0.0001, 0.0003, 0.0006, 0.001, 0.0015, 0.002, 0.0025, 0.003])
x = x * 133.322368
y = np.array([0.66, 0.69, 0.72, 0.79, 0.95, 1.0 ,1.05, 1.1, 1.14])
plt.plot(x, y)
spl = UnivariateSpline(x, y, k = 3)


x_pred = np.arange(start = 0.0, stop = 0.003 * 133.322368, step = 0.01)
y_perd = spl(x_pred)
y_lin = np.interp(x_pred, [0.001 * 133.322368, 0.003 * 133.322368] , [0.95, 1.14])

df = pd.DataFrame()
df["x_pred"] = x_pred
df["y_pred"] = y_perd
df["y_lin"] = y_lin
a = df.loc[df["x_pred"] <= 0.13332236800000002]["y_pred"]
b = df.loc[df["x_pred"] > 0.13332236800000002]["y_lin"]
c = pd.concat([a,b])
c = c.reset_index(drop=True)
df["pred"] = c

plt.plot(x, y)
plt.plot(df["x_pred"], df["pred"])
plt.show()

 