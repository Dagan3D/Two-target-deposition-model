import numpy as np
import pandas as pd
import Deposition_model
import Target_setup as ts
from Point import *
import scipy.constants as cnst
from scipy import interpolate


class model:

    """
    Класс в котором описаны методы вычисления
    параметров получаемого слоя исходя из параметров проведения процесса
    """

    def __init__(self, model):
        self.model = model

    def K_match(self, N_pressure, path=".\Data\Spray_rate.xls"):
        """Возвращает коэффициенты эффективности распыления мишеней

        Args:
            N_pressure (float): Давление азота в Па. Значение не должно входить за границ в файле .
            path (str, optional): Путь до файла со скоростью распыления. Defaults to ".\Data\Spray_rate.xls".

        Returns:
            dict: Словарь с коэффициентами распыления "Si": K_S_Si, "Mo": K_S_Mo.
        """

        input_df = pd.read_excel(path)
        input_df["P_N2"] = input_df["P"] * 133.322368
        input_df["P"] *= 1000

        Si_interpolation = interpolate.PchipInterpolator(
            input_df["P_N2"], input_df["Si"])
        Mo_interpolation = interpolate.PchipInterpolator(
            input_df["P_N2"], input_df["Mo"])

        Si_K = Si_interpolation(N_pressure)
        Si_K /= Si_interpolation(input_df["P_N2"][0])
        Mo_K = Mo_interpolation(N_pressure)
        Mo_K /= Mo_interpolation(input_df["P_N2"][0])

        k_S_Si = Si_K
        k_S_Mo = Mo_K

        return {"Si": float(k_S_Si), "Mo": float(k_S_Mo)}

    def predict(self, P_N2, P_O2, J_Si, J_Mo):
        K = self.K_match(P_N2)
        K_S_Si, K_S_Mo = K["Si"], K["Mo"]
        t1 = self.model.target_1
        t2 = self.model.target_2
        S_t1 = t1.S
        Sc_t1 = t1.S_compound
        t1.S *= K_S_Si
        t1.S_compound *= K_S_Si
        t1.J = J_Si

        S_t2 = t2.S
        Sc_t2 = t2.S_compound
        t2.S *= K_S_Mo
        t2.S_compound *= K_S_Mo
        t2.J = J_Mo

        F = self.model.F_of_P(P_N2)
        N_Si = t1.N_target_of_F(F)
        N_Mo = t2.N_target_of_F(F)
        N_O2 = self.model.N_O2_of_P(P_N2, P_residual=P_O2)
        N_N = self.model.N_gase_of_P(P_N2)[-1]

        N = N_Si + N_Mo + N_O2 + N_N
        C_Si = N_Si / N * 100
        C_Mo = N_Mo / N * 100
        C_N = N_N / N * 100
        C_O2 = N_O2 / N * 100

        predict = point(C_Si=C_Si,
                        C_Mo=C_Mo,
                        C_N=C_N,
                        C_O2=C_O2
                        )

        t1.S = S_t1
        t1.S_compound = Sc_t1
        t2.S = S_t2
        t2.S_compound = Sc_t2

        return predict


def fit(self):
    pass


if __name__ == "__main__":

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
    prediction = model_Si_Mo.predict(
        P_N2=0.5, P_O2=0.000001, J_Si=7.09, J_Mo=2.97)
    print(F"{prediction}\n\nОшибка - {point.loss(prediction, purpose):.3f}")
