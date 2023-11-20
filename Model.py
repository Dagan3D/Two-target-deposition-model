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
            Скорость откачки подаваемого газа (м3/c).
        alpha0_O2 : float
            Коэффициент задерживания молекулы кислорода.
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
        """
        Возвращает производную потока от давления.

        Parameters
        ----------
        P : TYPE
            Давление (Па).

        Returns
        -------
        dq_dp : TYPE
            Производную потока от давления.

        """
        
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
        """
        Часть потока идущая на мишень.

        Parameters
        ----------
        P : TYPE
            Давление (Па).
        target : ts.target
            Мишень.

        Returns
        -------
        qt : TYPE
            Поток на мишень.

        """
        F = self.F_of_P(P)
        t = target.t_of_F(F)
        qt = F * t * target.alpha0 * target.A
        return qt
        
    def qc_of_P(self, P, target: ts.target):
        """
        Часть потока идущая на камеру.

        Parameters
        ----------
        P : TYPE
            Давление (Па).
        target : ts.target
            Мишень.

        Returns
        -------
        qt : TYPE
            Поток на камеру.

        """
        F = self.F_of_P(P)
        c = target.c_of_F(F)
        qc = F * c * target.alpha0 * target.A_chamber
        return qc
    
    def qS_of_P(self, P):
        """
        Часть потока идущая в систему откачки.

        Parameters
        ----------
        P : TYPE
            Давление (Па).

        Returns
        -------
        qt : TYPE
            Поток на камеру.

        """
        S_a = self.S * self.K_calc()
        qS = self.F_of_P(P)*S_a
        return qS
    
    def N_gase_of_P(self, P):
        """
        Возвращает колличество реакционного газа в образце.

        Parameters
        ----------
        P : TYPE
            Давление.

        Returns
        -------
        Колличество реакционного газа в образце.

        """
        
        t_1 = self.target_1
        t_2 = self.target_2
        
        je_1 = t_1.J/cnst.elementary_charge
        je_2 = t_2.J/cnst.elementary_charge
        
        S_1 = t_1.S_compound
        S_2 = t_2.S_compound
        
        F = self.F_of_P(P)
        
        tetha_t_1 = 1 - t_1.t_of_F(F)
        tetha_t_2 = 1 - t_2.t_of_F(F)
        tetha_c_1 = 1 - t_1.c_of_F(F)
        tetha_c_2 = 1 - t_2.c_of_F(F)
        
        At_Ac_1 = t_1.A/t_1.A_chamber
        At_Ac_2 = t_2.A/t_2.A_chamber
       
        a_1 = je_1 * S_1 * tetha_t_1 * At_Ac_1                  #Азот из Si3N4 с мишени
        a_2 = je_2 * S_2 * tetha_t_2 * At_Ac_2                  #Азот из MoN с мишени
        
        b_1 = t_1.alpha0_c * F * (1 - tetha_c_1)                  #Азот захвченыый Si в полёте
        b_2 = t_2.alpha0_c * F * (1 - tetha_c_2)                  #Азот захвченыый Mo в полёте
        
        c_1 = t_1.alpha0_compound_c * F * tetha_c_1               #Азот захвченыый Si3N4 в полёте  
        c_2 = t_2.alpha0_compound_c * F * tetha_c_2               #Азот захвченыый MoN в полёте
                
        N = a_1 + a_2 + b_1 + b_2 + c_1 + c_2
        #N = N/2
        return (a_1, a_2, b_1, b_2, c_1, c_2, N)
        
    
    def N_O2_of_P (self, P, P_residual = 0.733e-3, D_0 = 20.95):
        """
        Возвращает колличество остаточноего кисорода в образце.

        Parameters
        ----------
        P : TYPE
            Давление.

        Returns
        -------
        Колличество реакционного газа в образце.

        """
        t_1 = self.target_1
        t_2 = self.target_2
        F = self.F_of_P(P)
        tetha_t_1 = 1 - t_1.t_of_F(F)
        tetha_t_2 = 1 - t_2.t_of_F(F)
        tetha_c_1 = 1 - t_1.c_of_F(F)
        tetha_c_2 = 1 - t_2.c_of_F(F)
        
        F = D_0/100 * (P_residual) / math.sqrt(2 * cnst.pi * cnst.Boltzmann * self.T * 32/1000/cnst.Avogadro)        
        
        F  = F * (1 - t_1.alpha0_O2 * (1 - tetha_t_1))
        
        O_Si_1 = F * ( t_1.alpha0_O2 * (1 - tetha_c_1))
        O_Si_2 = F * t_1.alpha0_O2_compound * tetha_c_1
        
        O_Mo_1 = F * (t_2.alpha0_O2 * (1 - tetha_c_2))
        O_Mo_2 = t_2.alpha0_O2_compound * F * tetha_c_2
                
        N = (O_Si_1 + O_Si_2 + O_Mo_1 + O_Mo_2)
        return N

    
    
   