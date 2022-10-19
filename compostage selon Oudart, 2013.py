#!/usr/bin/env python
# coding: utf-8

# In[4]:


import numpy as np
import pandas as pd
from pandas import Series,DataFrame
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.special import logsumexp


# In[5]:


class compostage:
    def __init__(self, echelle, procede, emissions, Qair, µ_hmax=5, S_r=10, K_s=1, flim_O2h=20, Nav=5, Klim = 5,
                pO2 = 2, K_O2 = 3, k_hr = 1, NXrb = 2, k_hs = 3, NX_sb= 4, tNXh = 3, µ_a=2, NX_a=3, Y_NO3=1,
                Kh=5, pIntLG=1, GP=2, Ts=1, rho_airsec=1, flim_hum=1,Vair=2, b_h=2, b_a = 3, f_Iaero = 0.2, t_NXi = 1):
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
        self.k_hr = k_hr
        self.NXrb =NXrb
        self.k_hs = k_hs
        self.NX_sb=NX_sb
        self.tNXh = tNXh #teneur en azote de la biomasse hétérotrophe
        self.µ_a=µ_a
        self.NX_a=NX_a
        self.Y_NO3=Y_NO3
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
        print('essai codage', self.echelle)
    
    "Emission de NH3"
    #---------------------------------
    def resolution(self): #azote disponible pour les microorganismes
        def system(temps, Y): 
            self.µ_h=Y[0]
            self.flimNav=Y[1]
            self.azote=Y[2]
            self.Xh = Y[3]
            self.NXrb = Y[4]
            
            self.µ_h = self.µ_hmax*(self.S_r/(self.S_r+self.K_s))*self.fonction_limO2h()*self.flimNav
            
            self.flimNav = (self.azote)/(logsumexp(self.azote)+self.Klim)
            
            self.dazote_dtemps = ((self.k_hr*self.azoteRB())+(self.k_hs*self.NX_sb))
            -((self.µ_h*self.Xh*self.tNXh)+((self.µ_a*self.Nmicro_auto())*(1+self.Y_NO3)))
            
            self.dXh_dtemps = np.exp(logsumexp(self.µ_h*self.Xh - self.b_h*self.Xh)) 
                               
            return [self.µ_h, self.flimNav, self.dazote_dtemps, self.dXh_dtemps, dNXrb_temps]
        
        self.solution = solve_ivp(system, [0, 100], [0, 1, 3, 2,0], method = 'RK45')
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
    
    def azoteRB(self):#dynamique de l'azote dans la fraction RB des M.O
        def equationN_rb(t, NXrb):
            return -self.k_hr*self.NXrb - (self.b_h*self.resolution()[3][0][55] * self.tNXh) - (self.b_a * self.Nmicro_auto()) - (self.Ninerte())
        self.solutionazoteRB = solve_ivp(equationN_rb, [0,100], [10], method = 'RK45', t_eval = np.arange(0,99,1))
        return self.solutionazoteRB.y[0][98]
   
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
    
    "Emission de CH4"
    
    "Emission de CO2"
    
    
    
            


# In[6]:


compostage1=compostage('grande', 'andain retourné', 'NH3', 100)


# In[7]:


compostage1.Nmicro_auto_death()


# In[18]:


compostage1.Nmicro_hetero_death()


# In[19]:


compostage1.Ninerte()


# In[ ]:


compostage1.resolution()


# In[ ]:


compostage1.resolution()[1]


# In[39]:


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


# In[55]:


NH3 = solve_ivp(emissionNH3, [0,100], [0])
plt.plot (NH3.temps, NH3.y[0])


# In[40]:



a = -1
def test(x,temps):
    return a
x=odeint(test, [20], temps)


# In[41]:


plt.plot(temps,x)


# In[ ]:




