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



def theta_t_of_alpha(k, n, alpha0_metal, alpha0_compound, F, J, S_compound):
    """
    Возвращает долю прореагировавшей поверхности мишени из коэффициент задерживания молекулы окислителя на поверхности мишени.
    
    Parameters
    ----------
    k : TYPE
        Количество атомов мишени в соединение (в Cr2O3 - 2).
    n : TYPE
        Количество атомов окислителя в соединение на один атом мишени (в Cr2O3 - 1.5). 
    alpha0_metal : TYPE
        Коэффициент задерживания молекулы окислителя на непрореагировавшей поверхности мишени.
    alpha0_compound : TYPE
        Коэффициент задерживания молекулы окислителя на прореагировавшей поверхности мишени.
    F : TYPE
        Молекулярный поток реактивного газа (молекул/(м2*с)) .
    J : TYPE
        Плотность потока ионов аргона распыляющих мешень (молекул/(м2*с)).
    S_compound : TYPE
        Коэффициент распыления прореагировавшего материала.

    Returns
    -------
    tetha : TYPE
        Доля прореагировавшей поверхности.

    """
    kn = k/n    
    a = kn * alpha0_metal * F
    b = (kn * F * (alpha0_metal - alpha0_compound)) + (J / cnst.elementary_charge * S_compound)
    tetha_t =  a / b
    
    return tetha_t



def theta_c_of_alpha(k, n, tetha_t, alpha0_metal, alpha0_compound, F, J, S_target, S_compound, A_t, A_c):
    """
    Возвращает долю прореагировавшей поверхности стенок

    Parameters
    ----------
    k : TYPE
        Количество атомов мишени в соединение (в Cr2O3 - 2).
    n : TYPE
        Количество атомов окислителя в соединение на один атом мишени (в Cr2O3 - 1.5). 
    tetha_t : TYPE
        Доля прореагировавшей поверхности мишени.
    alpha0_metal : TYPE
        DESCRIPTION.
    alpha0_compound : TYPE
        DESCRIPTION.
    alpha0_metal : TYPE
        Коэффициент задерживания молекулы окислителя на непрореагировавшей поверхности мишени.
    alpha0_compound : TYPE
        Коэффициент задерживания молекулы окислителя на прореагировавшей поверхности мишени.
    F : TYPE
        Молекулярный поток реактивного газа (молекул/(м2*с)) .
    J : TYPE
        Плотность потока ионов аргона распыляющих мешень (молекул/(м2*с)).
    S_target : TYPE
        Коэффициент распыления материала мишени.
    S_compound : TYPE
        Коэффициент распыления прореагировавшего материала.
    A_t : TYPE
        Пплощади поверхности мишеней.
    A_c : TYPE
        Площади поверхности стенок камеры.

    Returns
    -------
    tetha_c : TYPE
        Непрореагировашая площади поверхности стенок камеры.

    """
    
    kn = k/n    
    Je = J / cnst.elementary_charge
    A_tc = A_t / A_c
    
    a = (kn * alpha0_metal * F) + (Je * S_compound * tetha_t * A_tc) 
    b = (kn * alpha0_metal * F) + (Je * S_compound * tetha_t * A_tc) - (kn * alpha0_compound * F) + (Je * S_target * (1 - tetha_t) * A_tc)
    tetha_c =  a / b
    
    return tetha_c
    

if __name__ == "__main__":

    q = P_to_q(10, 20)
    print(q)

    F = F_of_P(10, 273, 16)
    print(F)
    
    theta_t = theta_t_of_alpha(1, 1, 10, 20, 30, 1, 20)
    print(theta_t)
    
    theta_с = theta_c_of_alpha(k = 4, n = 3, tetha_t = 10, alpha0_metal = 20,
                               alpha0_compound = 40, F = 1, J = 1,
                               S_target = 1, S_compound = 1, A_t = 2, A_c = 1)
    print(theta_с)
    
    
    