# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 16:22:23 2023

@author: butma
"""
import math
import scipy.constants as cnst

class target:
    
    def __init__(self, k, n, A, A_chamber, S, S_compound, alpha0, alpha0_compound, J ):
        """
        Инициалиирует мешени, камеру вокруг нее и молекулярный поток направленый в нее

        Parameters
        ----------
        k : TYPE
            Количество атомов мишени в соединение (в Cr2O3 - 2).
        n : TYPE
            Количество атомов окислителя в соединение на один атом мишени (в Cr2O3 - 1.5).
        A : TYPE
            Площадь мишени (м2).
        A_chamber : TYPE
            Подверженность поверхности камеры воздействию потока частиц (м2).
        S : TYPE
            Коэффициент распыления материала мишени.
        S_compoud : TYPE
            Коэффициент распыления прореагировавшего материала.
        alpha0 : TYPE
            Коэффициент задерживания молекулы окислителя на непрореагировавшей поверхности мишени.
        alpha0_compound : TYPE
            Коэффициент задерживания молекулы окислителя на прореагировавшей поверхности мишени.
        J : TYPE
            Плотность потока ионов аргона распыляющих мешень (молекул/(м2*с)).

        Returns
        -------
        None.

        """
        self.k = k
        self.n = n
        self.A = A
        self.A_chamber = A_chamber
        self.S = S
        self.S_compound = S_compound
        self.alpha0 = alpha0
        self.alpha0_compound = alpha0_compound
        self.J = J
    
    
    def t_of_F (self, F):
        """
        Возвращает t = (1 - theta_t)

        Parameters
        ----------
        F : TYPE
            Молекулярный поток реактивного газа (молекул/(м2*с)).

        Returns
        -------
        t : TYPE
            t = (1 - theta_t)

        """
        kn = self.k / self.n
        Je = self.J / cnst.elementary_charge
        
        a = kn * self.alpha0 * F
        b = Je * self.S_compound
        
        t = (a/b + 1)**-1
        return t
    
    def c_of_F (self, F):
        """
        Возвращает c = (1 - theta_c).

        Parameters
        ----------

        F : TYPE
            Молекулярный поток реактивного газа (молекул/(м2*с)).

        Returns
        -------
        c : TYPE
            c = (1 - theta_c).

        """
        kn = self.k / self.n
        eF = cnst.elementary_charge * F  
        A_t_c = self.A / self.A_chamber
        
        a = (self.S/self.S_compound) * (A_t_c)
        b1 = (kn * self.alpha0 * (eF / (self.J * self.S_compound)))**2
        b2 = eF / (self.J * self.S_compound) * (kn * self.alpha0 + A_t_c * kn * self.alpha0)
        b3 = self.S/self.S_compound * A_t_c
        
        c = a / (b1 + b2 + b3)
        return c
    
    
if __name__ == "__main__":
   Ti = target(k = 2, n = 2, A = 0.031, A_chamber = 0.2,
               S = 0.4, S_compound = 0.015,
               alpha0 = 1.0, alpha0_compound = 0.01, J = 60)
   
   print(Ti.t_of_F(10))
   print(Ti.c_of_F(10))
   
   