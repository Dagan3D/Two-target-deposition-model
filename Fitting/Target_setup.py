# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 16:22:23 2023

@author: butma
"""
import scipy.constants as cnst

class target:
    
    def __init__(self, k, n, t, A, A_chamber, S, S_compound, alpha0,
                 alpha0_compound, alpha0_c, alpha0_compound_c, alpha0_O2, alpha0_O2_compound, J, Ro, M):
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
        alpha0_c: TYPE
            alpha0 для поверхности осаждения.
        alpha0_compound : TYPE
            Коэффициент задерживания молекулы окислителя на прореагировавшей поверхности мишени.
        alpha0_compound_с : TYPE
            alpha0_compound для поверхности осаждения.
        J : TYPE
            Плотность потока ионов аргона распыляющих мешень (молекул/(м2*с)).
        Ro : Type
            Плотность вещества мишени (кг/м^3).
        M : Type
            Молярная масса вещества мишени (кг/м).

        Returns
        -------
        None.

        """
        self.k = k
        self.n = n
        self.t = t
        self.A = A
        self.A_chamber = A_chamber
        self.S = S
        self.S_compound = S_compound
        self.alpha0 = alpha0
        self.alpha0_c = alpha0_c
        self.alpha0_compound = alpha0_compound
        self.alpha0_compound_c = alpha0_compound_c
        self.alpha0_O2 = alpha0_O2
        self.alpha0_O2_compound = alpha0_O2_compound 
        self.J = J
        self.Ro = Ro
        self.M = M
    
    
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
        
        alpha0 = self.alpha0_c
        
        a = (self.S/self.S_compound) * (A_t_c)
        b1 = (kn * alpha0 * (eF / (self.J * self.S_compound)))**2
        b2 = eF / (self.J * self.S_compound) * (kn * alpha0 + A_t_c * kn * alpha0)
        b3 = self.S/self.S_compound * A_t_c
        
        c = a / (b1 + b2 + b3)
        return c
    
    def theta_t_of_F(self, F):
        """
        Возвращает долю прореагировавшей поверхности мишени из коэффициент задерживания молекулы окислителя на поверхности мишени.
        
        Parameters
        ----------
        F : TYPE
            Молекулярный поток реактивного газа (молекул/(м2*с)) .

        Returns
        -------
        tetha : TYPE
            Доля прореагировавшей поверхности мишени.

        """
        
        k = self.k
        n = self.n
        alpha0_target = self.alpha0 
        alpha0_compound = self.alpha0_compound
        J = self.J
        S_compound = self.S_compound
        
        kn = k/n    
        a = kn * alpha0_target * F
        b = (kn * F * (alpha0_target - alpha0_compound)) + (J / cnst.elementary_charge * S_compound)
        tetha_t =  a / b

        
        return tetha_t
    
    def theta_c_of_F(self, F):
        kn = self.k / self.n    
        Je = self.J / cnst.elementary_charge
        A_tc = self.A/ self.A_chamber
        
        a = (kn * self.alpha0 * F) + (Je * self.S_compound * self.theta_t_of_F(F) * A_tc) 
        b = (kn * self.alpha0 * F) + (Je * self.S_compound * self.theta_t_of_F(F) * A_tc)
        c = -1 * (kn * self.alpha0_compound * F) + (Je * self.S * (1 - self.theta_t_of_F(F)) * A_tc)
        tetha_c =  a / (b + c)
        return tetha_c
    
    
    def R_of_F (self, F, in_nm_s = False):
        """
        Возвращает распылённый поток из мишени
        
        Parameters
        ----------
        F : TYPE
            Молекулярный поток реактивного газа (молекул/(м2*с)) .

        Returns
        -------
        R : TYPE
            Распылённый из мишени поток (Атомы) / (ионы врезавшиеся в мишень * секунды * метр^2).
        """
        
        
        S_compoud = self.S_compound
        tetha_t = self.theta_t_of_F(F)
        S = self.S
        t = self.t_of_F(F)
        
        Je = self.J / cnst.elementary_charge
        
        R = Je * (S_compoud * tetha_t + S * t)
        if in_nm_s:
            R = (R / cnst.Avogadro * self.M) / self.Ro * 1e9
            return R
        return R
    
    
    def  N_target_of_F (self, F):
        """
        Возвращает колличество осевщего вещества мишении.
        
        Parameters
        ----------
        F : TYPE
            Молекулярный поток реактивного газа (молекул/(м2*с)) .

        Returns
        -------
        N : TYPE
            Колличество осевщего вещества мишении.
        """
        
        Je = self.J / cnst.elementary_charge
        S_compoud = self.S_compound
        tetha_t = 1 - self.t_of_F(F)
        S = self.S
        A_t = self.A
        A_c = self.A_chamber
        Je = self.J / cnst.elementary_charge
        
        N = (Je * S_compoud * tetha_t * self.t + Je * S * (1-tetha_t) ) * A_t/A_c
        return N
        
        
    
if __name__ == "__main__":
   Ti = target(k = 2, n = 2, A = 0.031, A_chamber = 0.2,
               S = 0.4, S_compound = 0.015,
               alpha0 = 1.0, alpha0_compound = 0.01, J = 60)
   
   print(Ti.t_of_F(10))
   print(Ti.c_of_F(10))
   
   