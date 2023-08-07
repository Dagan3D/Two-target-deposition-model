# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 16:51:26 2023

@author: butma
"""


import math
import scipy.constants as cnst
import target_setup as ts

class model:
      
    def __init__(self, target_1: ts.target , target_2: ts.target, T: float,
                 moleclar_mass: float, S: float, K1: float = 3.7e-21):
        """
        

        Parameters
        ----------
        target_1 : ts.target
            Первая мишень.
        target_2 : ts.target
            Вторая мишень.
        T : float
            Абсолютная температура реакционного газа (K).
        moleclar_mass : float
            Молекулярная масса реактивного газа (кг).
        S : float
            Скорость накачки подаваемого газа (м3/c).
        K1 : float, optional
            Коэффициент пересчёта. The default is 3.7e-21.

        Returns
        -------
        None.

        """
        self.target_1 = target_1
        self.target_2 = target_2
        self.T = T
        self.moleclar_mass = moleclar_mass
        self.S = S
        self.K1 = K1
        
            
    
    def q_of_P (P: float, S: float):
        """
        Возвращает скорость потка газа, зная давление и скорость накачки.
        
        Parameters
        ----------
        P : float
            Давление подаваемого газа (Па).
        S : float
            Скорость накачки подаваемого газа (м3/c).

        Returns
        -------
        q : float
            Скорость потока газа (Па*м3/с).
        """    
        q = P*S
        return q
    
    
    def F_of_P (self, P: float):
        """
        Возвращает молекулярнй поток реагирующего газа.

        Parameters
        ----------
        P : float
            Давление подаваемого газа (Па).
        T : float
            Абсолютная температура реакционного газа (K).
        M : float
            Молекулярная масса реактивного газа (кг).

        Returns
        -------
        F : float
            Молекулярный поток реактивного газа (молекул/(м2*с))
        """
        F = P/math.sqrt(2 * cnst.pi * cnst.Boltzmann * self.T * self.moleclar_mass)
        return F
    
    def K_calc (self):
        """
        Возвращает К, нужный для характеристической функции

        Parameters
        ----------
        T : float
            Абсолютная температура реакционного газа (K).
        M : float
            Молекулярная масса реактивного газа (кг).
        K1 : float, optional
            Коэфициент пересчета. The default is 3.7e-21.

        Returns
        -------
        K : float
            DESCRIPTION.

        """
        K = math.sqrt(2 * cnst.pi * cnst.Boltzmann * self.T * self.moleclar_mass)/self.K1
        return K
    
    def q_of_F_t_c (self, P):
        """
        Вычисляет поток кислорада (Pa*m^3/s)

        Parameters
        ----------
        P : TYPE
            Давление подаваемого газа (Па).
    
        Returns
        -------
        Поток кислорада (Pa*m^3/s).

        """
        F = self.F_of_P(P)
        
        S_a = self.S * self.K_calc()
        
        a1 = self.target_1.t_of_F(F) * self.target_1.alpha0 * self.target_1.A
        a2 = self.target_2.t_of_F(F) * self.target_2.alpha0 * self.target_2.A
        a3 = self.target_1.c_of_F(F) * self.target_1.alpha0 * self.target_1.A_chamber
        a4 = self.target_2.c_of_F(F) * self.target_2.alpha0 * self.target_2.A_chamber
           
        q_O2 = self.K1 * F * (a1 + a2 + a3 + a4 + S_a)
        return q_O2
    
    def dq_dp (self, P):
        
        F = self.F_of_P(P)
        
        t_1 = self.target_1.t_of_F(F)
        t_2 = self.target_2.t_of_F(F)
        c_1 = self.target_1.c_of_F(F)
        c_2 = self.target_2.c_of_F(F)
        
        a1 = self.target_1.alpha0 * self.target_1.A_chamber * c_1**2 * ( (self.target_1.S_compound/self.target_1.S) * (self.target_1.A_chamber/self.target_1.A) * ((1 - t_1)/t_1)**2 - 1)
        a2 = self.target_2.alpha0 * self.target_2.A_chamber * c_2**2 * ( (self.target_2.S_compound/self.target_2.S) * (self.target_2.A_chamber/self.target_2.A) * ((1 - t_2)/t_2)**2 - 1)
        a3 = self.target_1.alpha0 * self.target_1.A * t_1**2
        a4 = self.target_2.alpha0 * self.target_1.A * t_2**2
        
        dq_dp = self.S - (1/self.K_calc()) * (a1 + a2 - a3 - a4)
        return dq_dp
    
    def qt_of_P(self, P, target: ts.target):
        F = self.F_of_P(P)
        t = target.t_of_F(F)
        qt = F * t * target.alpha0 * target.A
        return qt
        
    def qc_of_P(self, P, target: ts.target):
        F = self.F_of_P(P)
        c = target.c_of_F(F)
        qc = F * c * target.alpha0 * target.A_chamber
        return qc
    
    def qS_of_P(self, P):
        S_a = self.S * self.K_calc()
        qS = self.F_of_P(P)*S_a
        return qS
        

    
    
   