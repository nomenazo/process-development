#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import pandas as pd
from pandas import Series,DataFrame
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.special import logsumexp
from math import exp,expm1


# In[49]:


class compostage:
    def __init__(self, echelle, procede, emissions, Qair, µ_hmax=5, S_r=10, K_s=1, flim_O2h=20, Nav=5, Klim = 5,
                pO2 = 2, K_O2 = 3, k_hr = 1, NXrb = 2, k_hs = 3, NX_sb= 4, tNXh = 3, µ_a=2, NX_a=3, Y_NO3=1,
                Kh=5, pIntLG=1, GP=2, Ts=1, rho_airsec=1, flim_hum=1,Vair=2, b_h=2, b_a = 3, f_Iaero = 0.2, t_NXi = 1, 
                 pmax_denit=2, f_limNO3 = 1, f_limTdenit = 1, WFPS = 1, pWFPS_denit = 0.5, pN2O_nit = 1, pN2O_denit= 0.3, 
                 Si = 1, Y_CH4=2, temps = 10, v_max = 1, Km = 0.5, CH4_out = 1 , beta = 0.5, f_OUR = 0.5, Y_CO2 =0.5,
                 Y_h = 0.5):
        self.echelle = echelle
        self.procede = procede
        self.emissions = emissions
        self.Qair= Qair
        self.µ_hmax = µ_hmax
        self.S_r=S_r
        self.K_s=K_s
        self.Nav=Nav
        self.Klim=Klim
        self.pO2=pO2
        self.K_O2=K_O2
        self.k_hr = k_hr #cte d'hydrolyse de la fraction rapidement biodégradable
        #pour la production de méthane voir s'il s'agit de la même cte d'hydrolyse
        self.NXrb =NXrb
        self.k_hs = k_hs #cte d'hydrolyse de la fraction lentement biodégradable
        self.NX_sb=NX_sb
        self.tNXh = tNXh #teneur en azote de la biomasse hétérotrophe
        self.µ_a=µ_a
        self.NX_a=NX_a
        self.Y_NO3=Y_NO3 #rendement spécifique de production de nitrate
        self.Kh=Kh #constante de Henry, qui va encore dépendre de Ts donc à résoudre
        self.pIntLG=pIntLG #paramètre d'échange
        self.GP=GP #constante des GP
        self.Ts=Ts #température sortie
        self.rho_airsec=rho_airsec
        self.flim_hum=flim_hum #fonction limite humidité, à résoudre
        self.Vair=Vair#Volume de l'air dans l'andain, à résoudre
        self.b_h = b_h
        self.b_a =b_a #taux spécifique de mortalité
        self.f_Iaero = f_Iaero #Proportion de matière inerte issue du décès de la biomasse
        self.t_NXi = t_NXi #teneur en N de la M.O inerte
        self.pmax_denit = pmax_denit #denitrification maximale, émission de N2 et N20
        self.f_limNO3 = f_limNO3 #fonction de limitation par le stock de nitrates
        self.f_limTdenit = f_limTdenit #fonction de limitation de la dénitrification par la température
        self.WFPS = WFPS
        self.pWFPS_denit = pWFPS_denit #seuil pour une dénitrification
        self.pN2O_nit = pN2O_nit #part d'émission de N2O sur l'ammonium nitrifié
        self.pN2O_denit = pN2O_denit #part maximale d'émission de N2O sur l'émission N2+N2O
        self.Si = Si #substrat insoluble initial????
        self.Y_CH4 = Y_CH4 #coefficient de rendement de CH4
        self.temps = temps
        self.v_max = v_max #maximum methane oxidation rate 
        self.CH4_out = CH4_out #methane dans la partie aerobic, à étudier : diffusion de CH4 de la partie anaérobique vers aérobique
        self.Km = Km #half saturation constant
        self.beta = beta #mass conversion factor of VS to VS_aero
        self.f_OUR = f_OUR #fonction de limitattion de l'oxydation du méthane par l'O2
        self.Y_CO2 = Y_CO2 #rendement de production en CO2
        self.Y_h = Y_h #rendement spécifique de la biomasse sur le substrat
        
        print('essai codage', self.echelle)
        
    #-----------------------------------------------
    "Durée du compostage"
    
    def durée(self):
        self.temps = np.arange(0,90,1)
        return [self.temps[0], self.temps[89]]
    
    
    "Emission de NH3"
    def resolution(self): #azote disponible pour les microorganismes
        def system(temps, Y): 
            self.µ_h=Y[0]
            self.flimNav=Y[1]
            self.azote=Y[2]
            self.Xh = Y[3]
            
            self.µ_h = self.µ_hmax*(self.S_r/(self.S_r+self.K_s))*self.fonction_limO2h()*self.flimNav
            
            self.flimNav = (self.azote)/(logsumexp(self.azote)+self.Klim)
            
            self.dazote_dtemps = ((self.k_hr*self.NXrb)+(self.k_hs*self.NX_sb))
            -((self.µ_h*self.Xh*self.tNXh)+((self.µ_a*self.Nmicro_auto())*(1+self.Y_NO3)))
            
            self.dXh_dtemps = np.exp(logsumexp(self.µ_h*self.Xh - self.b_h*self.Xh)) 
                               
            return [self.µ_h, self.flimNav, self.dazote_dtemps, self.dXh_dtemps]
        
        self.solution = solve_ivp(system, self.durée() , [0, 1, 3, 2], method = 'RK45')
        self.µ_h = self.solution.y[0] 
        self.flimNav=self.solution.y[1]
        self.azote=self.solution.y[2]
        self.Xh = self.solution.y[3]
        
        #solutions_system = [self.µ_h, self.flimNav, self.azote, self.Xh]
        return np.array([[self.µ_h],[self.flimNav],[self.azote],[self.Xh]])       
    
    def fonction_limO2h(self): #à faire : calculer pO2
        self.flim_O2h = self.pO2/(self.pO2+self.K_O2)
        return self.flim_O2h
    
    def Nmicro_auto(self):#azote dans les micro autotrophes
        #à faire : calculer µ_a
        def equationNmicro_auto(t, NXa):
            return self.µ_a*NXa - self.b_a*NXa
        
        self.solutionNmicro_auto = solve_ivp(equationNmicro_auto, [0,100], [1], method = 'RK45', t_eval = np.arange(0,99,1), max_step = 100)
        return self.solutionNmicro_auto.y[0][98]
    
    def Nmicro_auto_death(self): #cinétique de décès des micro autotrophes
        return - self.b_a * self.Nmicro_auto()
    
    def Nmicro_hetero_death(self): #cinétique de décès des micro hétérotrophes contenant les N
        return - self.b_h * self.resolution()[3][0][55] * self.tNXh
    
    def Ninerte(self): #azote dans la M.O inerte
        def equationNinerte(t,N_Xi):
            return (self.b_h*self.resolution()[3][0][55]*self.f_Iaero*self.t_NXi)*(1 + (self.Nmicro_auto_death()/self.Nmicro_hetero_death()))
            
        self.solutionNinerte = solve_ivp(equationNinerte, [0,100], [10], method = 'RK45', t_eval = np.arange(0,99,1))
        return self.solutionNinerte.y[0][98]
    
    #def azoteRB(self):#dynamique de l'azote dans la fraction RB des M.O
        #def equationN_rb(t, NX_rb):
            #return -self.k_hr*self.NX_rb - (self.b_h*self.resolution()[3][0][55] * self.tNXh) - (self.b_a * self.Nmicro_auto()) - (self.Ninerte())
        #self.solutionazoteRB = solve_ivp(equationN_rb, [0,100], [10], method = 'RK45', t_eval = np.arange(0,99,1))
        #return self.solutionazoteRB.y[0][98]
   
    def ammoniaque(self):#émissions ammoniacales
        def equationammoniac(t, NH3):
            return (self.Qair*0.9*self.resolution()[2][0][55]*self.Kh*self.pIntLG*self.flim_hum)/(self.rho_airsec*self.Vair*self.GP*self.Ts)
        
        self.solutionNH3 = solve_ivp(equationammoniac, [0,100],[0], method='RK45', t_eval=np.arange(0,99,1))
        return self.solutionNH3.y[0][98]
    
     
    #dynamique de l'azote dans la biomasse microbienne hétérotrophe NXh : proportionnelle à la teneur N
    #def Nmicro_hetero(self):
        #def equationNmicro_hetero ()
    
    #------------------------
    
    "Emission de N2 et N2O"
    
    #Production de NO3 associée à la croissance de la biomsse autotrophe
    
    def prodNO3(self):
        def equationNO3(t, NO3):
            return self.Y_NO3*self.µ_a*self.Nmicro_auto()
        
        self.solutionNO3 = solve_ivp(equationNO3, [0,100],[0], method='RK45', t_eval=np.arange(0,99,1))
        return self.solutionNO3.y[0][98]
    
    #Consommation de NO3 par dénitrification
    
    def consNO3(self):
        def equationNO3cons(t, NO3):
            return self.pmax_denit*NO3*self.f_limNO3* self.f_limTdenit
        self.solutionNO3cons = solve_ivp(equationNO3cons, [0,100],[2], method='RK45', t_eval=np.arange(0,99,1))
        return self.solutionNO3cons.y[0][98]
    
    #Emissions de N2O
    def emmissionN2O(self):
        def equationN2O(t, N2O):
            if self.WFPS<self.pWFPS_denit:
                return (self.pN2O_nit*self.prodNO3())+(self.pN2O_denit*self.consNO3())
            else:
                return self.pN2O_denit*(self.pN2O_nit*self.prodNO3())+(self.pN2O_denit*self.consNO3())
        self.solutionN2O = solve_ivp(equationN2O, [0,100],[2], method='RK45', t_eval=np.arange(0,99,1))
        return self.solutionN2O.y[0][98]
    
    #Emissions de N2
    def emissionN2(self):
        def equationN2(t, N2):
            if self.WFPS<self.pWFPS_denit:
                return (1-self.pN2O_denit)*self.prodNO3()
            else:
                return (1-self.pN2O_denit*(self.pN2O_nit*self.prodNO3()))+(1-self.pN2O_denit*self.consNO3())
        self.solutionN2 = solve_ivp(equationN2, [0,10],[2], method='RK45', t_eval=np.arange(0,9,1))
        return self.solutionN2.y[0][8]
 #---------------------------------------------------------   
    
    "Emission de CH4, Jinyi Ge et al. 2016"
    "L'émission de CH4 durant le compostage est la différence entre la génération de CH4 dans la partie anaérobique et l'oxydation de CH4 produit en CO2"
#Production de CH4
    def hydr_rate(self):
        self.Rsi = self.k_hr * self.Si* exp(-self.k_hr*self.temps) #constante d'hydrolyse à étudier
        return self.Rsi

    def prodCH4(self):
        self.CH4gen = self.Y_CH4*self.hydr_rate()
        return self.CH4gen
#oxydation de CH4
    def oxyCH4(self):
        self.v_oxi = (self.v_max*self.prodCH4()*self.f_OUR*self.beta)/(self.Km+self.prodCH4())
        #la fonction de limitation par la température devrait être introduite, mais nous n'avons pas les données
        return self.v_oxi
    
    def CH4emi(self):
        return self.prodCH4() - self.oxyCH4()
#-------------------------------------------------------------

    "Emission de CO2, Oudart"
    def CO2emis(self): #il s'agit du CO2 rpoduit par la croissance microbienne sans l'oxydation du méthane
        self.CO2 = (self.Y_CO2/self.Y_h) *self.resolution()[0][0][55]*self.resolution()[3][0][55]
        return self.CO2         


# In[50]:


compostage1=compostage('grande', 'andain retourné', 'NH3', 100)


# In[51]:


compostage1.durée()


# In[11]:


compostage1.prodNO3()


# In[12]:


compostage1.prodCH4()


# In[13]:


compostage1.consNO3()


# In[14]:


compostage1.oxyCH4()


# In[8]:


compostage1.CH4emi()


# In[15]:


compostage1.CO2emis()


# In[58]:


compostage1.emmissionN2O()


# In[59]:


compostage1.emissionN2()


# In[4]:


compostage1.Nmicro_auto_death()


# In[5]:


compostage1.Nmicro_hetero_death()


# In[14]:


compostage1.Ninerte()


# In[9]:


compostage1.resolution()


# In[53]:


compostage1.resolution()[1]


# In[6]:


compostage1.resolution().shape


# In[61]:


compostage1.resolution()[2][0][55]


# In[40]:


compostage1.Nmicro_auto()


# In[41]:


compostage1.fonction_limO2h()


# In[67]:


compostage1.ammoniaque()


# In[2]:


system


# In[167]:


compo=compostage(20, 'normale', 'NO3')


# In[16]:


class homecomposting(compostage): #classe fille
    def __init__(self, echelle, procede, emissions, types):
        compostage.__init__(self, echelle, procede, emissions) #hérite des pptés de la classe compostage
        self.types = types


# In[17]:


HC1=homecomposting(100, 'par retournement', 'NH3', 'petit')


# In[18]:


print ('Le premier procédé est le homecomposting avec un volume de', (HC1.echelle), 'avec un procede', HC1.procede, 'de type',
      HC1.types, ',avec des émissions de', HC1.emissions)


# In[19]:


isinstance(HC1, homecomposting)


# In[20]:


isinstance(HC1, compostage)


# In[21]:


issubclass(homecomposting, compostage)


# In[22]:


issubclass(compostage, homecomposting)


# In[23]:


#fonction limite par l'oxygène
pO2_biofilm = 0.5 #pression partielle du biofilm
K02_h =1 #cte de demi-saturation de l'oxygène pour la biomase hétérotrophe
flim_O2h = pO2_biofilm/(pO2_biofilm+K02_h)
flim_O2h


# In[47]:


KNh = 2
flimNav = Nav/(Nav+KNh)
plt.plot (flimNav, Nav)


# In[48]:


plt.plot (flimNav, temps)


# In[28]:


#ammoniac gazeux dans le tas
Kh = 2
Pechange = 2
Rgp = 12
Ts = 13
NH3gaz = 0.9 *Nav*Kh*Pechange/(Rgp*Ts)
plt.plot (NH3gaz, Nav)


# In[51]:


plt.plot (temps, NH3gaz)


# In[54]:


#volatilization de NH3
def emissionNH3(NH3, temps):
    Qair = 1 #débit d'air
    flim_hum = 1 #fontion limite par l'humidité
    rho_airsec = 12
    Vair = 10 #volume d'air dans le tas
    return self(Qair*flim_hum*NH3gaz)/(rho_airsec*Vair)

