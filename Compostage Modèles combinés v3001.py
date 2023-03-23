#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
from pandas import Series,DataFrame
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.special import logsumexp
from math import exp,expm1
from math import sqrt
from scipy import interpolate
from scipy.interpolate import interp1d


# In[2]:


pip install xlrd


# In[3]:


data ='Data.xlsx'


# In[4]:


operation = pd.read_excel(data, 'tech', index_col=0)
operation


# In[5]:


X= operation["u"]
X.loc['hc']


# In[6]:


operation.loc['CC','u']


# In[7]:


operation["u"]


# In[8]:


operation.iloc[0,:]


# In[9]:


operation.loc['hc',:]


# In[10]:


import xlrd


# In[11]:


#temperature profile for different types of aeration, interpolation from Rasapoor and al.
#1-Natural ventilation
NVS = pd.read_excel("NVS.xlsx")
x=NVS['time']
y=NVS['T']
courbe_inter = interp1d(x,y) 
NVS_inter = courbe_inter(np.arange(1,120,1))
plt.plot(NVS_inter)


# In[52]:


class Compostage:
    def __init__(self, data_path, aer, scale, region, QO2 = 3):
        
        self.scale = scale
        self.region = region
        self.QO2 = QO2 #Débit d'air entrant, modèles différents pour chaque type d'aération, à résoudre
        self.aer = aer
        
        with pd.ExcelFile(data_path) as f:
            
            #variables d'état 
            composition = pd.read_excel(data_path, 'Variables', index_col=0).fillna(0.0)
            mass = [composition[c] for c in composition.columns]
            self.mass = mass
            
            #Paramètres cinétiques
            parameters = pd.read_excel(data_path, 'kinetics', index_col=0).fillna(0.0)
            Valeurs = [parameters[c] for c in parameters.columns]
            self.Valeurs = Valeurs
            
            #paramètres stochiométriques
            stoech = pd.read_excel(data_path, 'stoe', index_col=0).fillna(0.0)
            val = [stoech[c] for c in stoech.columns]
            self.val=val
            
            #paramètres d'opérations
            oper = pd.read_excel(data_path, 'tech', index_col=0).fillna(0.0)
            self.cond = oper['u']
            self.area = oper['A']
            self.air = oper['Qair']
            
            #paramètres régionaux
            reg = pd.read_excel(data_path, 'reg', index_col=0).fillna(0.0)
            self.Temp = reg['Ta'] #température ambiante
            
            
        with pd.ExcelFile(aer) as f:
            temp = pd.read_excel(aer, 'temp')
            x=temp['time']
            y=temp['T']
            self.courbe_inter = interp1d(x,y)
            self.temp_inter = self.courbe_inter(np.arange(1,120,1))
            
    
    def v(self): #variables d'état
        self.massetotale = self.mass[0]
        return self.massetotale
    
    def ki(self): #paramètres cinétiques
        self.vitesse=self.Valeurs[0]
        return self.vitesse
    
    def sto(self): #coeff stoechiométriques
        self.stoe=self.val[0]
        return self.stoe
    
    def U(self): #heat transfer coeff
        if self.scale == 'home composting':
            return self.cond['hc']
        elif self.scale == 'community composting':
            return self.cond['CC']
        elif self.scale =='industrial composting':
            return self.cond['IC']
        else:
            print ("Unknown technology")
    
    def A(self): #surface area
        if self.scale == 'home composting':
            return self.area['hc']
        elif self.scale == 'community composting':
            return self.area['CC']
        elif self.scale =='industrial composting':
            return self.area['IC']
        else:
            print ("Unknown technology")
    
    
    def F(self): #airflow
        if self.scale == 'home composting':
            return self.air['hc']
        elif self.scale == 'community composting':
            return self.air['CC']
        elif self.scale =='industrial composting':
            return self.air['IC']
        else:
            print ("Unknown technology")
    
    def Ta(self): #ambiant temperature
        if self.region == 'France':
            return self.Temp['FR']
        else:
            return self.Temp['GLO']
        
    def Tint(self): #temperature profile
        return self.temp_inter
        
        #if self.aer == 'HT': #heap turning
            #return self.HT_int
        #elif self.aer == 'FA': #forced aeration
            #return self.FA_int
        #elif self.aer == 'NV': #natural ventilation
            #return self.NV_int
        #else:
            #print ("Unknown technology")
        
    
    def resolution(self):
        def system(t, Y):
            
            self.C=Y[0]
            self.P=Y[1]
            self.L=Y[2]
            self.H=Y[3]
            self.CE=Y[4]
            self.LG=Y[5]
            self.Xi=Y[6]
            self.Sc=Y[7]
            self.Sp=Y[8]
            self.Sl=Y[9]
            self.Sh=Y[10]
            self.Slg=Y[11]
            self.Xmb=Y[12]
            self.Xtb=Y[13]
            self.Xma=Y[14]
            self.Xta=Y[15]
            self.Xmf=Y[16]
            self.Xtf=Y[17]
            self.Xdb=Y[18]
            self.sO2=Y[19]
            #fonction de limitation
            self.fO2=Y[20]
            #paramètres cinétiques affectés
            self.µmb=Y[21]
            self.µtb=Y[22]
            self.µma=Y[23]
            self.µta=Y[24]
            self.µmf=Y[25]
            self.µtf=Y[26]
            #constantes d'hydrolyses affectés
            self.kh1C=Y[27]
            self.kh2P=Y[28]
            self.kh3L=Y[29]
            self.kh4C=Y[30]
            self.kh5P=Y[31]
            self.kh6L=Y[32]
            self.kh7H=Y[33]
            self.kh8CE=Y[34]
            self.kh9LG=Y[35]
            self.kh10H=Y[36]
            self.kh11CE=Y[37]
            self.kh12LG=Y[38]
            self.faB_C=Y[39]
            self.faB_L=Y[40]
            self.faA_C=Y[41]
            self.faA_L=Y[42]
            self.faA_H=Y[43]
            self.faF_C=Y[44]
            self.faF_L=Y[45]
            self.faF_H=Y[46]
            self.faF_LG=Y[47]
            #emissions
            self.CO2=Y[48]
            self.flT1=Y[49]
            self.flT2=Y[50]
            self.Xa = Y[51]
            self.Sn=Y[52]#dissolved ammonium
            self.NO3= Y[53]
            self.N2O = Y[54]
            self.N2 = Y[55]
            self.T = Y[56]
            #Evolution de la température
            #self.Qcond=Y[49]
            #self.Qbio=Y[50]
            #self.Qconv=Y[51]
            #self.T=Y[52]
            
          
            
            #Hydrolysis of substrate by microorganisms
            self.kh1C = self.ki()[0]*(self.Xmb/((self.ki()[29]*self.Xmb)+self.C))
            self.kh2P = self.ki()[1]*(self.Xmb/((self.ki()[29]*self.Xmb)+self.P))
            self.kh3L = self.ki()[2]*(self.Xmb/((self.ki()[29]*self.Xmb)+self.L))
            self.kh4C = self.ki()[3]*(self.Xtb/((self.ki()[29]*self.Xtb)+self.C))
            self.kh5P = self.ki()[4]*(self.Xtb/((self.ki()[29]*self.Xtb)+self.P))
            self.kh6L = self.ki()[5]*(self.Xtb/((self.ki()[29]*self.Xtb)+self.L))
            self.kh7H = self.ki()[6]*(self.Xta/((self.ki()[29]*self.Xta)+self.H))
            self.kh8CE = self.ki()[7]*(self.Xtf/((self.ki()[29]*self.Xtf)+self.CE))
            self.kh9LG = self.ki()[8]*(self.Xtf/((self.ki()[29]*self.Xtf)+self.LG))
            self.kh10H = self.ki()[9]*(self.Xma/((self.ki()[29]*self.Xma)+self.H))
            self.kh11CE = self.ki()[10]*(self.Xmf/((self.ki()[29]*self.Xmf)+self.CE))
            self.kh12LG = self.ki()[11]*(self.Xmf/((self.ki()[29]*self.Xmf)+self.LG))
            
            #Growth of microorganisms
            self.µmb = self.ki()[12]*self.fO2*self.flT1
            self.µtb = self.ki()[13]*self.fO2*self.flT2
            self.µma = self.ki()[14]*self.fO2*self.flT1
            self.µta = self.ki()[15]*self.fO2*self.flT2
            self.µmf = self.ki()[16]*self.fO2*self.flT1
            self.µtf = self.ki()[17]*self.fO2*self.flT2
            
            #Growth limitations
            #(1)Oxygen
            self.fO2 = self.sO2 /(self.ki()['kO2']+self.sO2)
            
            #(2)Substrate availability
            ##for bacteries:
            self.faB_C = self.Sc/np.exp(logsumexp(self.Sc+self.Sl))  ###availability of Sc
            self.faB_L = self.Sl/np.exp(logsumexp(self.Sc+self.Sl))  ###availability of Sc
            ##for actinomycetes
            self.faA_C = self.Sc/np.exp(logsumexp(self.Sc+self.Sl+self.Sh))
            self.faA_L= self.Sl/np.exp(logsumexp(self.Sc+self.Sl+self.Sh))
            self.faA_H= self.Sh/np.exp(logsumexp(self.Sc+self.Sl+self.Sh))
            ##for fungi
            self.faF_C = self.Sc/np.exp(logsumexp(self.Sc+self.Sl+self.Sh+self.Slg))
            self.faF_L= self.Sl/np.exp(logsumexp(self.Sc+self.Sl+self.Sh+self.Slg))
            self.faF_H= self.Sh/np.exp(logsumexp(self.Sc+self.Sl+self.Sh+self.Slg))
            self.faF_LG= self.Slg/np.exp(logsumexp(self.Sc+self.Sl+self.Sh+self.Slg))
    
            
            #(3)Temperature mesophilic (1) and thermophilic (2)
            self.flT1 = np.exp(logsumexp(((self.T - self.ki()['Tmax,1']) *((self.T - self.ki()['Tmin,1'])*(self.T - self.ki()['Tmin,1'])))
            /((self.ki()['Topt,1'] - self.ki()['Tmin,1'])*(((self.ki()['Topt,1'] - self.ki()['Tmin,1'])\
            *(self.T - self.ki()['Topt,1']))-((self.ki()['Topt,1'] - self.ki()['Tmax,1'])*\
                                              (self.ki()['Topt,1']+self.ki()['Tmin,1']-(2*self.T)))))))
            
            
            self.flT2 = np.exp(logsumexp(((self.T - (self.ki()['Tmax,2']))*((self.T - self.ki()['Tmin,2'])**2))            /((self.ki()['Topt,2'] - self.ki()['Tmin,2'])*(((self.ki()['Topt,2'] - self.ki()['Tmin,2'])            *(self.T - self.ki()['Topt,2']))-((self.ki()['Topt,2'] - self.ki()['Tmax,2'])*                                              (self.ki()['Topt,2']+self.ki()['Tmin,2']-(2*self.T)))))))
            
            
            #equations différentielles
            #aerobic degradation
            self.dC_dt= self.C * (-self.kh1C - self.kh4C)
            
            self.dP_dt = - (self.kh2P)*self.P - (self.kh5P)*self.P + (1 - self.ki()['fi'])*(self.ki()[24])*self.Xdb
            
            self.dL_dt = self.L * (-self.kh3L - self.kh6L)
            
            self.dH_dt = self.H * (-self.kh7H - self.kh10H)
            
            self.dCE_dt= self.CE * (-self.kh8CE - self.kh11CE)
            
            self.dLG_dt= self.LG * (-self.kh9LG - self.kh12LG)
            
            self.dXi_dt = self.ki()['fi'] * self.ki()[24]*self.Xdb
            
            self.dSc_dt = self.C * (self.kh1C + self.kh4C) +self.CE * (self.kh8CE + self.kh11CE)            -(self.µmb*self.faB_C*self.Xmb)-(self.µtb*self.faB_C*self.Xtb)-(self.µma*self.faA_C*self.Xma)-(self.µta*self.faA_C*self.Xta)            - (self.µmf*self.faF_C*self.Xmf) - (self.µtf*self.faF_C*self.Xtf)
            
            self.dSp_dt = self.P *(self.kh2P+self.kh5P) - (self.µmb*self.Xmb)-(self.µtb*self.Xtb)-            (self.µma*self.Xma)-(self.µta*self.Xta)-(self.µmf*self.Xmf) - (self.µtf*self.Xtf)
            
            
            self.dSl_dt = self.L *(self.kh3L + self.kh6L)-(self.µmb*self.faB_L*self.Xmb)-(self.µtb*self.faB_L*self.Xtb)-            (self.µma*self.faA_L*self.Xma)-(self.µta*self.faA_L*self.Xta)-(self.µmf*self.faF_L*self.Xmf) - (self.µtf*self.faF_L*self.Xtf)
            
            self.dSh_dt = self.H*(self.kh7H + self.kh10H)-(self.µma*self.faA_H*self.Xma)-(self.µta*self.faA_H*self.Xta)            -(self.µmf*self.faF_H*self.Xmf) - (self.µtf*self.faF_H*self.Xtf)
            
            self.dSlg_dt = self.LG * (self.kh9LG + self.kh12LG)-(self.µmf*self.faF_LG*self.Xmf) - (self.µtf*self.faF_LG*self.Xtf)
            
            self.dXmb_dt = (self.Xmb * self.µmb)*((self.faB_C*self.sto()['a'])+self.sto()['b']+(self.faB_L*self.sto()['c']))-(self.ki()[18]*self.Xmb)
            
            self.dXtb_dt = (self.Xtb * self.µtb)*((self.faB_C*self.sto()['a'])+self.sto()['b']+(self.faB_L*self.sto()['c']))-(self.ki()[19]*self.Xtb)
            
            self.dXma_dt = (self.Xma * self.µma)*((self.faA_C*self.sto()['a'])+self.sto()['b']+(self.faA_L*self.sto()['c'])+(self.faA_H*self.sto()['d']))            -(self.ki()[20]*self.Xma)
            
            self.dXta_dt = (self.Xta * self.µta)*((self.faA_C*self.sto()['a'])+self.sto()['b']+(self.faA_L*self.sto()['c'])+(self.faA_H*self.sto()['d']))            -(self.ki()[21]*self.Xta)
            
            self.dXmf_dt = (self.Xmf*self.µmf)*((self.faF_C*self.sto()['e'])+(self.sto()['f'])+(self.faF_L*self.sto()['g'])+(self.faF_H*self.sto()['h'])+(self.faF_LG*self.sto()['i']))            -(self.sto()['z']*self.ki()[22]*self.Xmf)
            
            self.dXtf_dt= (self.Xtf *self.µtf)*((self.faF_C*self.sto()['e'])+self.sto()['f']+(self.faF_L*self.sto()['g'])+            (self.faF_H*self.sto()['h'])+(self.faF_LG*self.sto()['i']))-(self.sto()['z']*self.ki()[23]*self.Xtf)
            
            self.dXdb_dt =(self.ki()[18]*self.Xmb)+(self.ki()[19]*self.Xtb)+(self.ki()[20]*self.Xma)+(self.ki()[21]*self.Xta)+(self.ki()[22]*self.Xmf)+(self.ki()[23]*self.Xtf)
            
            self.dsO2_dt = self.QO2 - (((self.µmb*self.Xmb)+(self.µtb*self.Xtb))*((self.faB_C*(6-(5*self.sto()['a'])))+((33/2)-(5*self.sto()['b']))+(self.faB_L*((134/4)-(5*(self.sto()['c']))))))-            (((self.µma*self.Xma)+(self.µta*self.Xta))*((self.faA_C*(6-(5*self.sto()['a'])))+((33/2)-(5*self.sto()['b']))+(self.faA_L*((134/4)-(5*self.sto()['c'])))+                                                      (self.faA_H*(10-(5*self.sto()['d'])))))-(((self.µmf*self.Xmf)+(self.µtf*self.Xtf))            *((self.faF_C*(6-((21/2)*self.sto()['e'])))+((33/2)-((21/2)*self.sto()['f']))+(self.faF_L*((139/4)-((21/2)*self.sto()['g'])))+(self.faF_H*(10-((21/2)*self.sto()['h'])))+             (self.faF_LG*((49/2)-((21/2)*self.sto()['i'])))))
            
            self.dCO2_dt = ((self.µmb*self.Xmb)+(self.µtb*self.Xtb))*((self.faB_C*(6-(5*self.sto()['a'])))+(16-(5*self.sto()['b']))+(self.faB_L*(25-(5*(self.sto()['c'])))))+            ((self.µma*self.Xma)+(self.µta*self.Xta))*((self.faA_C*(6-(5*self.sto()['a'])))+(16-(5*self.sto()['b']))+(self.faA_L*(25-(5*self.sto()['c'])))+                                                      (self.faA_H*(10-(5*self.sto()['d']))))+((self.µmf*self.Xmf)+(self.µtf*self.Xtf))            *((self.faF_C*(6-(10*self.sto()['e'])))+(16-(10*self.sto()['f']))+(self.faF_L*(25-(10*self.sto()['g'])))+(self.faF_H*(10-(10*self.sto()['h'])))+             (self.faF_LG*(20-(10*self.sto()['i'])))) #Sole-Mauri considers a transfer of CO2 dissolved in liquid phase to gas phase
            
            ##Nitrification - denitrification
            self.dXa_dt = self.Xa *((self.sto()['l']*self.ki()['µa'])-self.ki()['ba']) #flim of growth of autotroph micro
            
            self.dSn_dt = (((self.µmb*self.Xmb)+(self.µtb*self.Xtb))*((self.faB_C*(-self.sto()['a']))+(4-self.sto()['b'])-(self.faB_L*self.sto()['c'])))+            (((self.µma*self.Xma)+(self.µta*self.Xta))*((-self.faA_C*self.sto()['a'])+(4-self.sto()['b'])-(self.faA_L*self.sto()['c'])-            (self.faA_H*self.sto()['d'])))+(((self.µmf*self.Xmf)+(self.µtf*self.Xtf))            *(-(self.faF_C*self.sto()['e'])+(4-self.sto()['f'])-(self.faF_L*self.sto()['g'])-(self.faF_H*self.sto()['h'])-             (self.faF_LG*self.sto()['i']))) - self.ki()['µa'] #-NH3 émis
            
            self.dNO3_dt = (self.sto()['n']*self.ki()['µa']*self.Xa)-(self.ki()['pmaxdenit']*self.NO3)#fonction de limitation NO3 et T à insérer
            
            self.dN2O_dt = self.ki()['pmaxdenit'] * self.ki()['pdenit']*self.NO3 #flimNO3 et T
            
            self.dN2_dt = self.ki()['pmaxdenit'] * (1-self.ki()['pdenit'])*self.NO3
            
            self.dT_dt = self.Tint() #à résoudre
            
            
            ##Methane production
            #a-Methane generation
            #self.dCH4gen_dt = self.Ych4*(-(self.C*self.khC)-(self.P*self.khp)-(self.L*self.khl)-(self.H*self.khh)\
                                         #-(self.CE*self.khCE))
            #b-Methane oxidation
            #self.dCH4oxi_dt = self.vmax*(self.CH4gen/(self.Km+self.CH4gen))
            
            #c-emission of methane
            #self.CH4emis = self.CH4gen - self.CH4oxi
            
            
            
            ##Energy balance
            #a-Heat transfer by conduction
            #self.dQcond_dt = self.U() * self.A() *(self.T - self.Ta())
            #b-biological heat 
            #self.dQbio_dt = self.ki()['HcO2'] * self.CO2 
            #c-convection by air flow (sensible heat change)
            #self.dQconv_dt = self.F()*(self.ki()['Ca']*(self.Ta()-self.T))
            #d-latent heat change : evaporation
            #self.dQvap_dt = self.F()*(self.ki()['hv']*(self.ki()['hout'] -self.ki()['hin'])
            #c-temperature of compost
            #self.dT_dt = (self.Qbio - self.Qcond - self.Qconv)\
            #/((self.ki()['Cp_C']*(self.C))+(self.ki()['Cp_L']*(self.L))\
              #+(self.ki()['Cp_P']*(self.P))+(self.ki()['Cp_LG']*(self.LG))+(self.ki()['Cp_HE']*(self.H)))
            
            ##Liquide-gas transfer  
            
            
            return [self.dC_dt, self.dP_dt, self.dL_dt, self.dH_dt, self.dCE_dt,self.dLG_dt,self.dXi_dt, self.dSc_dt,                     self.dSp_dt, self.dSl_dt, self.dSh_dt, self.dSlg_dt,self.dXmb_dt, self.dXtb_dt, self.dXma_dt,                    self.dXta_dt,self.dXmf_dt,self.dXtf_dt,self.dXdb_dt, self.dsO2_dt, self.fO2, self.µmb, self.µtb, self.µma,                   self.µta, self.µmf, self.µtf, self.kh1C, self.kh2P, self.kh3L, self.kh4C, self.kh5P, self.kh6L, self.kh7H,                   self.kh8CE, self.kh9LG, self.kh10H, self.kh11CE, self.kh12LG,self.faB_C,self.faB_L,self.faA_C,self.faA_L,self.faA_H,self.faF_C,self.faF_L,                   self.faF_H,self.faF_LG, self.dCO2_dt, self.flT1,self.flT2, self.dXa_dt, self.Sn, self.dNO3_dt, self.dN2O_dt, self.dN2_dt,self.T]
            
        self.solution = solve_ivp(system, [0,119], ([self.v()['G'],self.v()['P'],self.v()['L'],self.v()['HE'],self.v()['CE'],self.v()['LG'], self.v()['Xi'],                                                    self.v()['Sc'],self.v()['Sp'],self.v()['Sl'],self.v()['Sh'],self.v()['Slg'],                                             self.v()['Xmb'],self.v()['Xtb'],self.v()['Xma'],self.v()[15],self.v()['Xmf'],                                                    self.v()['Xtf'],self.v()['Xdb'],self.v()['S(O2)'], 0, self.ki()[12],                                                  self.ki()[13],self.ki()[14],self.ki()[15],self.ki()[16],self.ki()[17],self.ki()[0],self.ki()[1],self.ki()[2],                                                  self.ki()[3],self.ki()[4],self.ki()[5],self.ki()[6],self.ki()[7],self.ki()[8],self.ki()[9],self.ki()[10],                                                  self.ki()[11],0,0,0,0,0,0,0,0,0,0,1,1,self.v()['Xa'], self.v()['S(NH3)'], 0, 0                                                  ,0,self.Tint()[0]]), t_eval = np.arange(0,110) ,method = 'RK45')
        return self.solution.y[56]
        #return plt.plot(self.solution.t, self.solution.y[56])
                           
                           
            
      


# In[53]:


t_eval("=np.arange(3,110,1)")


# In[54]:


Data1 = 'Data.xlsx'


# In[55]:


essai = Compostage(Data1, 'NVS.xlsx', 'home composting', 'France')


# In[56]:


essai.ki()[0]


# In[57]:


essai.U()


# In[58]:


essai.Tint()


# In[59]:


essai.Tint().shape


# In[60]:


plt.plot(essai.Tint())


# In[61]:


t_eval = np.arange(1,121,1)
t_eval


# In[62]:


essai.Tint()[0]


# In[63]:


essai.Tint().shape


# In[64]:


essai.A()


# In[65]:


essai.resolution()


# In[51]:


essai.resolution().shape


# In[216]:


essai.A()


# In[259]:


essai.ki()['Tmax,1']


# In[99]:


i=555/12
i


# In[ ]:




