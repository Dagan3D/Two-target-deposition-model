# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 10:02:02 2023.

@author: Danil
"""


import math
import scipy.constants as cnst

def P_to_q (P, S):
    """
    Возвращает скорость потка газа, зная давление и скорость накачки.
    
    Parameters
    ----------
    P : TYPE
        Давление подаваемого газа (Па).
    S : TYPE
        Скорость накачки подаваемого газа (м3/c).

    Returns
    -------
    q : TYPE
        Скорость потока газа (Па*м3/с).
    """    
    q = P*S
    return q

def F_of_P (P, T, M):
    """
    Возвращает молекулярнй поток реагирующего газа.

    Parameters
    ----------
    P : TYPE
        Давление подаваемого газа (Па).
    T : TYPE
        Абсолютная температура реакционного газа (K).
    M : TYPE
        Молекулярная масса реактивного газа (кг).

    Returns
    -------
    F : TYPE
        Молекулярный поток реактивного газа (молекул/(м2*с))
    """
    F = P/math.sqrt(2 * cnst.pi * cnst.Boltzmann * T * M)
    return F


def theta_of_alpha(k, n, alpha0_metal, alpha0_compound, F, J, S_compound):
    """
    Возвращает долю прореагировавшей поверхности из коэффициент задерживания молекулы окислителя на поверхности мишени.
    
    Parameters
    ----------
    k : TYPE
        Количество атомов мишени в соединение (в Cr2O3 - 2).
    n : TYPE
        Количество атомов окислителя в соединение на один атом мишени (в Cr2O3 - 1.5) 
    alpha0_metal : TYPE
        Коэффициент задерживания молекулы окислителя на непрореагировавшей поверхности мишени.
    alpha0_compound : TYPE
        Коэффициент задерживания молекулы окислителя на прореагировавшей поверхности мишени.
    F : TYPE
        Молекулярный поток реактивного газа (молекул/(м2*с)) .
    J : TYPE
        Плотность потока ионов аргона распыляющих мешень (молекул/(м2*с)).
    S_compound : TYPE
        Прореагировавшая плотность мишени (м2).

    Returns
    -------
    tetha : TYPE
        Доля прореагировавшей поверхности.

    """
    kn = k/n    
    tetha = (kn * alpha0_metal * F)/(F * (alpha0_metal - alpha0_compound) + (J/cnst.elementary_charge)*S_compound)
    return tetha




if __name__ == "__main__":

    q = P_to_q(10, 20)
    print(q)

    F = F_of_P(10, 273, 16)
    print(F)
    
    t = theta_of_alpha(1, 1, 10, 20, 30, 1, 20)
    print(t)
    
    
    