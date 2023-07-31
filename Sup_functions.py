# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 10:02:02 2023.

@author: Danil
"""

import math
import scipy.constants as cnst

def q_of_P (P, S):
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


def theta_t_of_alpha(k, n, alpha0_target, alpha0_compound, F, J, S_compound):
    """
    Возвращает долю прореагировавшей поверхности мишени из коэффициент задерживания молекулы окислителя на поверхности мишени.
    
    Parameters
    ----------
    k : TYPE
        Количество атомов мишени в соединение (в Cr2O3 - 2).
    n : TYPE
        Количество атомов окислителя в соединение на один атом мишени (в Cr2O3 - 1.5). 
    alpha0_target : TYPE
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
        Доля прореагировавшей поверхности мишени.

    """
    kn = k/n    
    a = kn * alpha0_target * F
    b = (kn * F * (alpha0_target - alpha0_compound)) + (J / cnst.elementary_charge * S_compound)
    tetha_t =  a / b
    
    return tetha_t



def theta_c_of_alpha(k, n, tetha_t, alpha0_target, alpha0_compound, F, J, S_target, S_compound, A_t, A_c):
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
    alpha0_target : TYPE
        Коэффициент задерживания молекулы окислителя на непрореагировавшей поверхности мишени.
    alpha0_compound : TYPE
        Коэффициент задерживания молекулы окислителя на прореагировавшей поверхности мишени.
    F : TYPE
        Молекулярный поток реактивного газа (молекул/(м2*с)).
    J : TYPE
        Плотность потока ионов аргона распыляющих мешень (молекул/(м2*с)).
    S_target : TYPE
        Коэффициент распыления материала мишени.
    S_compound : TYPE
        Коэффициент распыления прореагировавшего материала.
    A_t : TYPE
        Площади поверхности мишеней
    A_c : TYPE
        Площадь поверхности держателя подложки + стенки камеры.

    Returns
    -------
    tetha_c : TYPE
        Доля прореагировавшей поверхности стенок камеры.

    """
    
    kn = k/n    
    Je = J / cnst.elementary_charge
    A_tc = A_t / A_c
    
    a = (kn * alpha0_target * F) + (Je * S_compound * tetha_t * A_tc) 
    b = (kn * alpha0_target * F) + (Je * S_compound * tetha_t * A_tc) - (kn * alpha0_compound * F) + (Je * S_target * (1 - tetha_t) * A_tc)
    tetha_c =  a / b
    
    return tetha_c
    
    
def q_of_tetha (alpha0_target, alpha0_compound, F, theta, A, K1 = 3.7e-21):
    """
    Возврщает поток потребляемого реакционного газа.

    Parameters
    ----------
    alpha0_target : TYPE
        Коэффициент задерживания молекулы окислителя на непрореагировавшей поверхности мишени.
    alpha0_compound : TYPE
        Коэффициент задерживания молекулы окислителя на прореагировавшей поверхности мишени.
    F : TYPE
        Молекулярный поток реактивного газа (молекул/(м2*с)).
    theta : TYPE
        Доля прореагировавшей поверхности.
    A : TYPE
        Площадь поверхности.
    K1 : TYPE, optional
        Коэффициент пересчёта. The default is 3.7e-21.

    Returns
    -------
    Поток потребляемого реакционного газа (Pa*m^3/s).

    """
    q = K1 * ((alpha0_target * F * (1 - theta)) + (alpha0_compound * F * theta)) * A
    return q


def R_of_tetha (J, S_compound, S_target, tetha_t):
    """
    Возвращает скорость распления мишени

    Parameters
    ----------
    J : TYPE
        Плотность потока ионов аргона распыляющих мешень (молекул/(м2*с)).
    S_compound : TYPE
        Коэффициент распыления прореагировавшего материала.
    S_target : TYPE
        Коэффициент распыления материала мишени.
    tetha_t : TYPE
        Доля прореагировавшей поверхности мишени.

    Returns
    -------
    R : TYPE
        Скорость распления мишени (молекул/(м2*с)).

    """
        
    Je = J / cnst.elementary_charge
    R = Je * (S_compound * tetha_t  +  S_target * (1 - tetha_t))
    return R

def t_of_F (k, n, J, alpha0_target, F, S_compound):
    """
    Возвращает t = (1 - theta_t)

    Parameters
    ----------
    k : TYPE
        Количество атомов мишени в соединение (в Cr2O3 - 2).
    n : TYPE
        Количество атомов окислителя в соединение на один атом мишени (в Cr2O3 - 1.5).
    J : TYPE
        Плотность потока ионов аргона распыляющих мешень (молекул/(м2*с)).
    alpha0_target : TYPE
        Коэффициент задерживания молекулы окислителя на непрореагировавшей поверхности мишени.
    F : TYPE
        Молекулярный поток реактивного газа (молекул/(м2*с)).
    S_compound : TYPE
        Коэффициент распыления прореагировавшего материала.

    Returns
    -------
    t : TYPE
        t = (1 - theta_t)

    """
    kn = k / n
    Je = J / cnst.elementary_charge
    
    a = kn * alpha0_target * F
    b = Je * S_compound
    
    t = (a/b + 1)**-1
    return t

def c_of_F (k, n,  alpha0_target, alpha0_compound, F, J, S_target, S_compound, A_t, A_c):
    """
    Возвращает c = (1 - theta_c).

    Parameters
    ----------
    k : TYPE
        Количество атомов мишени в соединение (в Cr2O3 - 2).
    n : TYPE
        Количество атомов окислителя в соединение на один атом мишени (в Cr2O3 - 1.5). 
    alpha0_target : TYPE
        Коэффициент задерживания молекулы окислителя на непрореагировавшей поверхности мишени.
    alpha0_compound : TYPE
        Коэффициент задерживания молекулы окислителя на прореагировавшей поверхности мишени.
    F : TYPE
        Молекулярный поток реактивного газа (молекул/(м2*с)).
    J : TYPE
        Плотность потока ионов аргона распыляющих мешень (молекул/(м2*с)).
    S_target : TYPE
        Коэффициент распыления материала мишени.
    S_compound : TYPE
        Коэффициент распыления прореагировавшего материала.
    A_t : TYPE
        Площади поверхности мишеней
    A_c : TYPE
        Площадь поверхности держателя подложки + стенки камеры.

    Returns
    -------
    c : TYPE
        c = (1 - theta_c).

    """
    kn = k / n
    eF = cnst.elementary_charge * F  
    A_t_c = A_t / A_c
    
    a = (S_target/S_compound) * (A_t_c)
    b1 = (kn * alpha0_target * eF / (J * S_compound))**2
    b2 = eF / (J * S_compound) * (kn * alpha0_target + A_t_c * kn * alpha0_target)
    b3 = S_target/S_compound * A_t_c
    
    c = a / (b1 + b2 + b3)
    return c

def K_calc (T, M0, K1 = 3.7e-21):
    """
    Возвращает К, нужный для характеристической функции

    Parameters
    ----------
    T : TYPE
        Абсолютная температура реакционного газа (K).
    M : TYPE
        Молекулярная масса реактивного газа (кг).
    K1 : TYPE, optional
        Коэфициент пересчета. The default is 3.7e-21.

    Returns
    -------
    K : TYPE
        DESCRIPTION.

    """
    K = math.sqrt(2 * cnst.pi * cnst.Boltzmann * T * M0)/K1
    return K

def q_of_F_t_c (F, t_1, t_2, c_1, c_2, 
                A_t1, A_t2, A_c1, A_c2,
                alpha0_target1, alpha0_target2, S_a, K1 = 3.7e-21):
    """
    Вычисляет поток кислорада (Pa*m^3/s)

    Parameters
    ----------
    F : TYPE
        Молекулярный поток реактивного газа (молекул/(м2*с)).
    t_1 : TYPE
        DESCRIPTION.
    t_2 : TYPE
        DESCRIPTION.
    c_1 : TYPE
        DESCRIPTION.
    c_2 : TYPE
        DESCRIPTION.
    A_t1 : TYPE
        Площади поверхности мишени 1.
    A_t2 : TYPE
        Площади поверхности мишени 2.
    A_c1 : TYPE
        Площадь поверхности держателя подложки + стенки камеры у мишени 1.
    A_c2 : TYPE
        Площадь поверхности держателя подложки + стенки камеры у мишени 2.
    alpha0_target1 : TYPE
        Коэффициент задерживания молекулы окислителя на непрореагировавшей поверхности мишени 1.
    alpha0_target2 : TYPE
        Коэффициент задерживания молекулы окислителя на непрореагировавшей поверхности мишени 2.
    S_a : TYPE
        Пересчтаный коэфициент откачки S_a = KS.
    K1 : TYPE, optional
        Коэфициент пересчета. The default is 3.7e-21.

    Returns
    -------
    Поток кислорада (Pa*m^3/s).

    """
    
    a1 = t_1 * alpha0_target1 * A_t1 
    a2 = t_2 * alpha0_target2 * A_t2
    a3 = c_1 * alpha0_target1 * A_c1
    a4 = c_2 * alpha0_target2 * A_c2
       
    q_O2 = K1 * F * (a1 + a2 + a3 + a4 + S_a)
    return q_O2



def dq_dp (S, K, t_1, t_2, c_1, c_2,
               A_t1, A_t2, A_c1, A_c2,
               S_target_1, S_compound_1,
               S_target_2, S_compound_2,
               alpha0_target1, alpha0_target2):
    """
    Вычисляет дифирениал потока кислорода от изменения его давления (Pa*m^3/s).

    Parameters
    ----------
    S : TYPE
        Скорость перекачки (м3/с)
    K : TYPR
        Коэфицент пересчёта с чётом тмпературы и массы молекулы газа см. K_calc
    t_1 : TYPE
        DESCRIPTION.
    t_2 : TYPE
        DESCRIPTION.
    c_1 : TYPE
        DESCRIPTION.
    c_2 : TYPE
        DESCRIPTION.
    A_t1 : TYPE
        Площади поверхности мишени 1.
    A_t2 : TYPE
        Площади поверхности мишени 2.
    A_c1 : TYPE
        Площадь поверхности держателя подложки + стенки камеры у мишени 1.
    A_c2 : TYPE
        Площадь поверхности держателя подложки + стенки камеры у мишени 2.
    S_target : TYPE
        Коэффициент распыления материала мишени 1.
    S_compound : TYPE
        Коэффициент распыления прореагировавшего материала 1.
    S_target : TYPE
        Коэффициент распыления материала мишени 2.
    S_compound : TYPE
        Коэффициент распыления прореагировавшего материала 2.
    alpha0_target1 : TYPE
        Коэффициент задерживания молекулы окислителя на непрореагировавшей поверхности мишени 1.
    alpha0_target2 : TYPE
        Коэффициент задерживания молекулы окислителя на непрореагировавшей поверхности мишени 2.
    S_a : TYPE
        Пересчтаный коэфициент откачки S_a = KS.
    K1 : TYPE, optional
        Коэфициент пересчета. The default is 3.7e-21.

    Returns
    -------
    dq/dp.

    """

    
    a1 = alpha0_target1 * A_c1 * c_1**2 * ( (S_compound_1/S_target_1) * (A_c1/A_t1) * ((1 - t_1)/t_1)**2 - 1)
    a2 = alpha0_target2 * A_c2 * c_2**2 * ( (S_compound_2/S_target_2) * (A_c2/A_t2) * ((1 - t_2)/t_2)**2 - 1)
    a3 = alpha0_target1 * A_t1 * t_1**2
    a4 = alpha0_target2 * A_t2 * t_2**2
    
    dq_dp = S - (1/K) * (a1 + a2 - a3 - a4)
    return dq_dp


if __name__ == "__main__":

    q = q_of_P(10, 20)
    print(q)

    F = F_of_P(10, 273, 16)
    print(F)

    theta_t = theta_t_of_alpha(1, 1, 10, 20, 30, 1, 20)
    print(theta_t)

    theta_с = theta_c_of_alpha(k=4, n=3, tetha_t=10, alpha0_target=20,
                               alpha0_compound=40, F=1, J=1,
                               S_target=1, S_compound=1, A_t=2, A_c=1)
    print(theta_с)
