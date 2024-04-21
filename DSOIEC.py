#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 19:10:21 2024

@author: pm.deoliveiradejes
"""

import math as math
import cmath
from pyomo.environ import *
from pyomo.opt import SolverFactory
import numpy as np
# from gurobipy import *
#import gurobipy as gp

import matplotlib.pyplot as plt
#from pyomo.util.infeasible import log_infeasible_constraints
# LOAD CONFIG MODEL
Sbase= 10 # MVA
Vbase= 69 # kV
Ebase= 10 # MWh
PQ = False
# PQ = True # Congestion Habilitada

# CONG_LIMIT = True
CONG_LIMIT = False # Congestion Deshabilitada
    
T = 24

# ---------------- HIDROGENO ----------------------



# cons_energ = 50 #KWh/Kg
# cons_agua = 0.015 #m3/Kg
# precio_agua = 5 #$/m3

# PH2O = 1 #MW
# PH2O_standby = 0.05 #5% Warm Start Seconds
# PH2O_idle = 0.005 # 0.5% Cold Start Seconds


# --------------------------------------------------
PG2_max = np.array([0,0,0,0,0,0,0.811697096,1.82361281,2.884230349,3.804153725,4.480513858,5.043290512,5.254331757,4.848483209,4.588740138,3.809565039,2.873407721,1.829024124,0.768406584,0,0,0,0,0]);
lambda1 = np.array([1042.90,957.10,957.10,952.81,940.59,947.19,1000.00,1196.04,1319.14,1324.09,1273.60,1179.87,1221.12,1217.16,1162.38,1133.00,1101.32,1244.88,1333.66,1353.47,1359.74,1358.42,1314.52,1148.51]); #Wholesale spot market


l1_max = np.max(lambda1)
l1_min = np.min(lambda1)
lambda2_max  = 10000;#$/pu.h
lambda3_max  = 10000;#$/pu.h
epsilon2=(lambda2_max-1000)/0.5;#$/pu^2.h
epsilon3=(lambda3_max-1000)/3.5;#$/pu^2.h

PD2_max=lambda2_max/epsilon2;
PD3_max=lambda3_max/epsilon3;
PD2_min=(lambda2_max-1000)/epsilon2
PD3_min=(lambda3_max-1000)/epsilon3
# PD2_min=0
# PD3_min=0
# elasticity2_min=-1/(epsilon2*PD2_min/1000);#inelastic
# elasticity3_min=-1/(epsilon3*PD3_min/1000);#inelastica
# PD2_0min=(lambda2_max-min(lambda1))/epsilon2; 
# PD3_0min=(lambda3_max-min(lambda1))/epsilon3; 
# elasticity2_min=-inv(epsilon2*PD2_0min/min(lambda1));%inelastic
# elasticity3_min=-inv(epsilon3*PD3_0min/min(lambda1));%inelastica
# PD2_0max=(lambda2_max-max(lambda1))/epsilon2;
# PD3_0max=(lambda3_max-max(lambda1))/epsilon3;
# elasticity2_max=-inv(epsilon2*PD2_0max/max(lambda1));
# elasticity3_max=-inv(epsilon3*PD3_0max/max(lambda1));

PF = 0.8

beta=0
alpha=0 # OPEX
gamma=0 #CAPEX

# BATTERY MODEL
C_3 = 10 #pu
rho=0.2 # RHO pct charge2
PG3_max=3 #pu

## 69 kV Network model
# Phase: 266.8 MCM 65/7 ACSR 460A
# Shield wire/neutral: Copperweld 3/8"
# ATP Analisis: 69kV_7.acp
#
#   G *<--3m--> --
#     |
#     |________ 6m
#     |     A ! --
#     |________ 3m
#     |     B ! --
#     |________ 3m
#     |     C ! --
#     |
#     |         15.2m
#     |
#     |___________
#      \\\\\\\\\\\
#
# Sbase=10MVA, Vbase=69kV
zline=0.240544+0.481652j; #ohms/km
Zbase=(Vbase**2)/Sbase;#ohms
zpu=zline/Zbase;#pu/km
L12=30;#km
L13=10;#km
L23=15;#km
y12=1/(L12*zpu);#siemens
y13=1/(L13*zpu);#siemens
y23=1/(L23*zpu);#siemens
Ybus = np.array([[y12+y13,-y12,-y13],[-y12,y12+y23,-y23], [-y13,-y23,y13+y23]])#siemens
G=Ybus.real
B=Ybus.imag

if CONG_LIMIT:
    Smax = 2;
else:
    Smax = PF*math.sqrt(3)*.460*Vbase/Sbase #pu

# -------------------  Modelo  -----------------------
model = ConcreteModel()

# -------------------  Variables -----------------------
model.T = RangeSet(1,T)
model.T1 = RangeSet(1,T+1)

range_PG2 = dict.fromkeys(range(1,T+1),(0,1))
for i in range(0,T):
    range_PG2[i+1] = (0,PG2_max[i])
# print(range_PG2)
model.PG2 = Var(model.T, initialize = 0, bounds = range_PG2)

model.PD2 = Var(model.T, bounds = (0,PD2_max))
# model.QD2 = Var(model.T)

model.PG3 = Var(model.T, bounds = (-PG3_max,PG3_max),initialize = 0)

model.PD3 = Var(model.T, bounds = (0, PD3_max))
# model.QD3 = Var(model.T)

if PQ:
    model.QG3 = Var(model.T, initialize = 0)
    model.QG2 = Var(model.T, initialize = 0)

model.E3 = Var(model.T1, bounds = (C_3*rho,C_3),initialize = C_3*rho)

model.L2 = Var(model.T, bounds = (0,lambda2_max))
model.L3 = Var(model.T, bounds = (0,lambda3_max))

range_L1 = dict.fromkeys(range(1,T+1),1)
for i in range(0,T):
    range_L1[i+1] = lambda1[i]
model.L1 = Param(model.T, initialize = range_L1)

model.k = RangeSet(1,3)

lb = { 1: 1, 2: 0.8, 3:0.8}
ub = { 1: 1, 2: 1.2, 3:1.2}

def fb(model, k, t):
    return (lb[k], ub[k])

model.V = Var(model.k, model.T, bounds=fb, initialize = 1)

lb = { 1: 0, 2: None, 3:None}
ub = { 1: 0, 2: None, 3:None}

def fb_th(model, k, t):
    return (lb[k], ub[k])

model.Th = Var(model.k, model.T, bounds=fb_th, initialize = 0)

# -------------------  Objetivo  -----------------------

def PG1_func(i):
    return model.V[1,i]*(sum((model.V[k,i]*(G[0,k-1]*cos((model.Th[1,i] - model.Th[k,i])) + B[0,k-1]*sin(model.Th[1,i] - model.Th[k,i]))) for k in model.k))
#     return 1
## FUncion Objetuvo SW Social Welfare
csp_exp = sum((((lambda2_max*model.PD2[i]-(epsilon2*(model.PD2[i]**2)/2))-(model.L2[i]*model.PD2[i]) ) + ((lambda3_max*model.PD3[i]-(epsilon3*(model.PD3[i]**2)/2))-(model.L3[i]*model.PD3[i]))) for i in model.T) # Consumer surplus
psp_exp = sum(((model.L2[i]*model.PG2[i] - (((1/2)*beta*(model.PG2[i]**2)+alpha*model.PG2[i]+gamma)))) for i in model.T) # Producer surplus
ssp_exp = sum((model.L3[i]*model.PG3[i]) for i in model.T) # Storage surplus
nsp_exp = sum((model.L1[i]*(-PG1_func(i)) + model.L2[i]*(model.PD2[i]-model.PG2[i]) + model.L3[i]*(model.PD3[i]-model.PG3[i])) for i in model.T) #Network surplus

obj_exp = psp_exp + csp_exp + ssp_exp + nsp_exp
model.obj = Objective(expr = obj_exp, sense=maximize)

# -------------------  Restricciones  -----------------------

def Sij_cons(model, t, i, j):
    Pij = model.V[i, t] * model.V[j, t] * (G[i-1,j-1] * cos( model.Th[i, t]-model.Th[j, t] ) + B[i-1,j-1] * sin( model.Th[i, t]-model.Th[j, t] )) - G[i-1,j-1] * model.V[i, t]**2
    Pji = model.V[i, t] * model.V[j, t] * (G[i-1,j-1] * cos( model.Th[j, t]-model.Th[i, t] ) + B[i-1,j-1] * sin( model.Th[j, t]-model.Th[i, t] )) - G[i-1,j-1] * model.V[j, t]**2
    Qij = model.V[i, t] * model.V[j, t] * (G[i-1,j-1] * sin( model.Th[i, t]-model.Th[j, t] ) - B[i-1,j-1] * cos( model.Th[i, t]-model.Th[j, t] )) + B[i-1,j-1] * model.V[i, t]**2;
    Qji = model.V[i, t] * model.V[j, t] * (G[i-1,j-1] * sin( model.Th[j, t]-model.Th[i, t] ) - B[i-1,j-1] * cos( model.Th[j, t]-model.Th[i, t] )) + B[i-1,j-1] * model.V[j, t]**2;
    pLij = abs(Pij-Pji)/2;
    qLij = abs(Qij-Qji)/2;
    return (pLij**2 + qLij**2) <= Smax**2;


def S12_cons(model, t):
    return Sij_cons(model, t, 1, 2)
model.S12_cons = Constraint(model.T, rule = S12_cons)

def S13_cons(model, t):
    return Sij_cons(model, t, 1, 3)
model.S13_cons = Constraint(model.T, rule = S13_cons)

def S23_cons(model, t):
    return Sij_cons(model, t, 2, 3)
model.S23_cons = Constraint(model.T, rule = S23_cons)

def PGj_cons(model,i,j):
    return ((model.V[j,i]*(sum((model.V[k,i]*(G[j-1,k-1]*cos((model.Th[j,i] - model.Th[k,i])) + B[j-1,k-1]*sin(model.Th[j,i] - model.Th[k,i]))) for k in model.k))))
def QGj_cons(model,i,j):
    return ((model.V[j,i]*(sum((model.V[k,i]*(G[j-1,k-1]*sin((model.Th[j,i] - model.Th[k,i])) - B[j-1,k-1]*cos(model.Th[j,i] - model.Th[k,i]))) for k in model.k))))

def PG2_cons(model,i):
    return PGj_cons(model,i,2) + model.PD2[i] == model.PG2[i]
model.PG2_cons = Constraint(model.T, rule = PG2_cons)

def PG3_cons(model,i):
    return PGj_cons(model,i,3) + model.PD3[i] == model.PG3[i]
model.PG3_cons = Constraint(model.T, rule = PG3_cons)

def QG2_cons(model,i):
    return (QGj_cons(model,i,2) + (model.PD2[i]*math.tan(math.acos(PF))) == model.QG2[i])
if PQ:
    model.QG2_cons = Constraint(model.T, rule = QG2_cons)

def QG3_cons(model,i):
    return (QGj_cons(model,i,3) + (model.PD3[i]*math.tan(math.acos(PF))) == model.QG3[i])
if PQ:
    model.QG3_cons = Constraint(model.T, rule = QG3_cons)

def L2_cons(model,i):
    return (lambda2_max - epsilon2*model.PD2[i]) == model.L2[i]
model.L2_cons = Constraint(model.T, rule = L2_cons)

def L3_cons(model,i):
    return (lambda3_max - epsilon3*model.PD3[i]) == model.L3[i]
model.L3_cons = Constraint(model.T, rule = L3_cons)

model.cons_e3 = ConstraintList()
model.cons_e3.add(model.E3[1] == model.E3[25])
model.cons_e3.add(model.E3[1] == C_3*rho)
for i in range(1,25):
    model.cons_e3.add(model.E3[i+1] == model.E3[i] - model.PG3[i])
    

# SolverFactory('mindtpy').solve(model, mip_solver='glpk', nlp_solver='ipopt')

SolverFactory('mindtpy').solve(model,
#                                    strategy='FP',
#                                    init_strategy='FP',
                                 #  mip_solver='glpk',
                                  # nlp_solver='ipopt',
#                                    add_regularization='level_L1',
#                                    solution_pool=True,
#                                    num_solution_iteration=10, # default=5
                                   tee=True
                                   )

model.obj.display()

# ----------------- Graficas -----------------

def PG1_func(i):
    return value(model.V[1,i])*(sum((value(model.V[k,i])*(G[0,k-1]*cos((value(model.Th[1,i]) - value(model.Th[k,i]))) + B[0,k-1]*sin(value(model.Th[1,i]) - value(model.Th[k,i])))) for k in model.k))

PG1_values = []
for i in range(1,25):
    PG1_values.append(PG1_func(i))
#-----------------------------------------------------
def Csp2_exp(i):
    U2 = (lambda2_max*value(model.PD2[i])-(epsilon2*(value(model.PD2[i])**2)/2))
    CD2 = (value(model.L2[i])*value(model.PD2[i]))
    return (U2 - CD2) # Consumer surplus
def Csp3_exp(i):
    U3 = (lambda3_max*value(model.PD3[i])-(epsilon3*(value(model.PD3[i])**2)/2))
    CD3 = (value(model.L3[i])*value(model.PD3[i]))
    return (U3 - CD3) # Consumer surplus
def Csp_exp(i):
    return Csp2_exp(i) + Csp3_exp(i) # Consumer surplus
#------------------------------------------------------
def Psp2_exp(i):
#     return value(model.L2[i])*value(model.PG2[i]) - (((1/2)*beta*(value(model.PG2[i])**2)+alpha*value(model.PG2[i])+gamma))
    return value(model.L2[i])*value(model.PG2[i])
# def Psp3_exp(i):
#     return (value(model.L3[i])*value(model.PG3[i])) # Producer surplus
def Psp_exp(i):
    return  Psp2_exp(i)  # Producer surplus
# ---------------------------------------------------
def Ssp_exp(i):
    return value(model.L3[i])*value(model.PG3[i]) # Storage surplus

def Nsp_exp(i):
    nsp_1 = value(model.L1[i])*(-PG1_values[i-1])
    nsp_2 = value(model.L2[i])*(value(model.PD2[i])-value(model.PG2[i]))
    nsp_3 = value(model.L3[i])*(value(model.PD3[i])-value(model.PG3[i]))
#     print(nsp_1)
#     print(nsp_2)
#     print(nsp_3)
    return (nsp_1 + nsp_2 + nsp_3) #Network surplus

def Obj_exp(i):
    return Csp_exp(i) + Psp_exp(i) + Ssp_exp(i) + Nsp_exp(i)

Obj_values = []
prev = 0
for i in range(1,25):
    Obj_values.append(Obj_exp(i)+prev)
    prev += Obj_exp(i)
    
Csp_values = []
prev = 0
for i in range(1,25):
    Csp_values.append(Csp_exp(i)+prev)
    prev += Csp_exp(i)
#----------------
Psp_values = []
prev = 0
for i in range(1,25):
    Psp_values.append(Psp_exp(i)+prev)
    prev += Psp_exp(i)
Psp2_values = []
prev = 0
for i in range(1,25):
    Psp2_values.append(Psp2_exp(i))
# Psp3_values = []
# prev = 0
# for i in range(1,25):
#     Psp3_values.append(Psp3_exp(i))
#-----------------------
Ssp_values = []
prev = 0
for i in range(1,25):
    Ssp_values.append(Ssp_exp(i)+prev)
    prev += Ssp_exp(i)

Nsp_values = []
prev = 0
for i in range(1,25):
    Nsp_values.append(Nsp_exp(i)+prev)
    prev += Nsp_exp(i)

# print(list(range(1,24)))
# print(Obj_exp(24))
# print(Ssp_exp(24))
# print(Psp_exp(24))
# for i in range(1,25):
#     print(Csp_exp(i))
# print(Nsp_exp(24))

E3ypoints = list(model.E3.extract_values().values())
E3xpoints = list(range(1,26))
plt.plot(E3xpoints, E3ypoints, label = "E3")

PG3ypoints = list(model.PG3.extract_values().values())
PG3xpoints = list(range(1,25))
plt.plot(PG3xpoints, PG3ypoints, label = "PG3")

PG2ypoints = list(model.PG2.extract_values().values())
PG2xpoints = list(range(1,25))
plt.plot(PG2xpoints, PG2ypoints, label = "PG2")

PG1ypoints = PG1_values
PG1xpoints = list(range(1,25))
plt.plot(PG1xpoints, PG1ypoints, label = "PG1")

PD3ypoints = list(model.PD3.extract_values().values())
PD3xpoints = list(range(1,25))
plt.plot(PD3xpoints, PD3ypoints, label = "PD3")

PD2ypoints = list(model.PD2.extract_values().values())
PD2xpoints = list(range(1,25))
plt.plot(PD2xpoints, PD2ypoints, label = "PD2")

plt.ylabel("Potencia PU")
plt.xlabel("Horas")
plt.legend(ncol=6)
plt.title("Potencia Activa")
plt.grid(which = "both")
plt.minorticks_on()
plt.tick_params(which = "minor", bottom = False, left = False)
plt.show()

if PQ:
    plt.figure()
    QG3ypoints = list(model.QG3.extract_values().values())
    QG3xpoints = list(range(1,25))
    plt.plot(QG3xpoints, QG3ypoints, label = "QD3")

    QG2ypoints = list(model.QG2.extract_values().values())
    QG2xpoints = list(range(1,25))
    plt.plot(QG2xpoints, QG2ypoints, label = "QD2")

    plt.ylabel("Potencia PU")
    plt.xlabel("Horas")
    plt.legend(ncol=6)
    plt.title("Potencia Reactiva")
    plt.grid(which = "both")
    plt.minorticks_on()
    plt.tick_params(which = "minor", bottom = False, left = False)
    plt.show()

plt.figure()
L3ypoints = list(model.L3.extract_values().values())
L3xpoints = list(range(1,25))
plt.plot(L3xpoints, L3ypoints, label = "L3")

L2ypoints = list(model.L2.extract_values().values())
L2xpoints = list(range(1,25))
plt.plot(L2xpoints, L2ypoints, label = "L2")

L1ypoints = list(model.L1.extract_values().values())
L1xpoints = list(range(1,25))
plt.plot(L1xpoints, L1ypoints, label = "L1")

plt.ylabel("Precio   $ PU")
plt.xlabel("Horas")
plt.legend(ncol=3)
plt.title("Lambda")
plt.grid(which = "both")
plt.minorticks_on()
plt.tick_params(which = "minor", bottom = False, left = False)
plt.show()

plt.figure()

dict1 = list(model.V.extract_values().values()) 
print(len(dict1))


V1 = dict1[0:24]
V2 = dict1[24:48]
V3 = dict1[48:72]

L3xpoints = list(range(1,25))
plt.plot(L3xpoints, V3, label = "V3")

L2xpoints = list(range(1,25))
plt.plot(L2xpoints, V2, label = "V2")

plt.ylabel("PU")
plt.xlabel("Horas")
plt.legend(ncol=3)
plt.title("Voltaje")
plt.grid(which = "both")
plt.minorticks_on()
plt.tick_params(which = "minor", bottom = False, left = False)
plt.show()

plt.figure()

dict1 = list(model.Th.extract_values().values()) 
print(len(dict1))


Th1 = dict1[0:24]
Th2 = dict1[24:48]
Th3 = dict1[48:72]

L3xpoints = list(range(1,25))
plt.plot(L3xpoints, Th3, label = "Th3")

L2xpoints = list(range(1,25))
plt.plot(L2xpoints, Th2, label = "Th2")

plt.ylabel("PU")
plt.xlabel("Horas")
plt.legend(ncol=3)
plt.title("Angulo")
plt.grid(which = "both")
plt.minorticks_on()
plt.tick_params(which = "minor", bottom = False, left = False)
plt.show()


plt.figure()

plt.plot(range(1,25), Obj_values, label = "Social Welfare")
plt.plot(range(1,25), Csp_values, label = "ConsumerSP")
plt.plot(range(1,25), Psp_values, label = "ProducerSP")
plt.plot(range(1,25), Ssp_values, label = "StorageSP")
plt.plot(range(1,25), Nsp_values, label = "NetworkSP")

plt.ylabel("Beneficio $")
plt.xlabel("Horas")
plt.legend(ncol=5)
plt.title("Beneficios")
plt.grid(which = "both")
plt.minorticks_on()
plt.tick_params(which = "minor", bottom = False, left = False)
plt.show()
print(value(model.obj))
# print(Obj_values[23])
# print(Csp_values[23])
# print(Psp_values[23])
# print(Ssp_values[23])
# print(Nsp_values[23])

# EP2 = sum(model.PD3.extract_values().values())
# EP3 = sum(model.PD3.extract_values().values())

EG1 = sum(model.PD3.extract_values().values())
EG2 = sum(model.PG2.extract_values().values())
EG3 = sum(model.PG3.extract_values().values())

#Probar init E3 diferente de 0
# ----------------- Reporte -----------------

SolarSold = round(sum(L2ypoints[i] * PG2ypoints[i] for i in range(0,24)),2)

PG3_bought = 0;
PG3_sold = 0;
EG3_bought = 0;
EG3_sold = 0;
for i in range(1,25):
    val = value(model.PG3[i])
    if(val>0):
        PG3_sold += val*value(model.L3[i])
        EG3_sold += val
    else:
        PG3_bought += val*value(model.L3[i])
        EG3_bought += val

PG1_bought = 0;
PG1_sold = 0;
EG1_bought = 0;
EG1_sold = 0;
for i in range(1,25):
    val = PG1_values[i-1]
    if(val>0):
        PG1_sold += val*value(model.L1[i])
        EG1_sold += val
    else:
        PG1_bought += val*value(model.L1[i])
        EG1_bought += val

DemandUtility = round(sum(((lambda2_max*PD2ypoints[i])-(epsilon2*(PD2ypoints[i]**2)/2)) + (lambda3_max*PD3ypoints[i]-(epsilon3*(PD3ypoints[i]**2)/2)) for i in range(0,24)),2)
DemandCost = round(sum(L2ypoints[i] * PD2ypoints[i] + L3ypoints[i] * PD3ypoints[i] for i in range(0,24)),2)
ECost = round(sum(PD3ypoints[i] + PD2ypoints[i] for i in range(0,24))*Sbase,2)


#Error signo B en PQ
def perdidas(t):
    P12 = perdidas_ij(1,2,t)
    P13 = perdidas_ij(1,3,t)
    P23 = perdidas_ij(2,3,t)
  
    return P12 + P13 + P23

def perdidas_ij(i,j,t):#1,2
    Pij = value(model.V[i, t]) * value(model.V[j, t]) * (G[i-1,j-1] * cos( value(model.Th[i, t])-value(model.Th[j, t]) ) + B[i-1,j-1] * sin( value(model.Th[i, t])-value(model.Th[j, t]) )) - G[i-1,j-1] * value(model.V[i, t])**2
    Pji = value(model.V[i, t]) * value(model.V[j, t]) * (G[i-1,j-1] * cos( value(model.Th[j, t])-value(model.Th[i, t]) ) + B[i-1,j-1] * sin( value(model.Th[j, t])-value(model.Th[i, t]) )) - G[i-1,j-1] * value(model.V[j, t])**2
    return Pij + Pji

Perdi_value = []
for i in range(0,24):
    Perdi_value.append(perdidas(i+1))
# print(Perdi_value)
Perdi_tot = sum(Perdi_value)
# print(Perdi_tot)

    
print('Optimization results:\n') 
print('Robust Energy Community Social Welfare ', round(value(model.obj),2) ,' Eur/day\n')
print('-------------------------------------------------------\n')

print('Solar PV Energy Sold                Eur/day', SolarSold ,', MWh/day', round(EG2*Sbase,2),'\n')
print('Solar PV Cost                       Eur/day 0 \n')
print('Solar PV Surplus                    Eur/day ',SolarSold,' \n')
print('-------------------------------------------------------\n')

print('Energy sold by the storage          Eur/day ', round(PG3_sold,2),', MWh/day', round(EG3_sold*Sbase,2),'  \n')
print('Energy bought by the storage        Eur/day ', round(PG3_bought,2),', MWh/day', round(EG3_bought*Sbase,2),'  \n')
print('-------------------------------------------------------\n')

# print('Storage Surplus                     Eur/day ', round(Ssp_values[23],2),' \n')
# print('Producer Surplus                    Eur/day ', round(Psp_values[23],2),' \n')
print('Demand Benefit                      Eur/day ',DemandUtility,' \n')
print('Demand Energy bougth                Eur/day ',DemandCost,', MWh/day',ECost,'  \n')
# print('Consumer Surplus                    Eur/day ',round(Csp_values[23],2),' \n')
# print('Network Surplus                     Eur/day ',round(Nsp_values[23],2),'\n')
print('-------------------------------------------------------\n')
print('Energy sold by the market          Eur/day ', round(PG1_sold,2),', MWh/day', round(EG1_sold*Sbase,2),'  \n')
print('Energy bought by the market        Eur/day ', round(PG1_bought,2),', MWh/day', round(EG1_bought*Sbase,2),'  \n')

print('-------------------------------------------------------\n')
print('------------------     Surplus    ---------------------\n')

print('Storage Surplus                     Eur/day ', round(Ssp_values[23],2),' \n')
print('Producer Surplus                    Eur/day ', round(Psp_values[23],2),' \n')
print('Consumer Surplus                    Eur/day ', round(Csp_values[23],2),' \n')
print('Network Surplus                     Eur/day ', round(Nsp_values[23],2),'\n')
print('Community Social Welfare            Eur/day ', round(Ssp_values[23] + Nsp_values[23] + Psp_values[23] + Csp_values[23],2),'\n')
print('  \n')

print('-------------------------------------------------------\n')
print('------------------ Energy Balance ---------------------\n')
print('Energy Sold by Solar PV             MWh/day ', round(EG2*Sbase,2),'\n')
print('Energy sold by the storage          MWh/day ', round(EG3_sold*Sbase,2),'  \n')
print('Energy sold by the market           MWh/day ', round(EG1_sold*Sbase,2),'  \n')
print('Total Vendido                       MWh/day ', round((EG1_sold+EG3_sold+EG2)*Sbase,2),'  \n')
print('  \n')

print('Demand Energy bougth                MWh/day ', ECost,'  \n')
print('Energy bought by the storage        MWh/day ', round(EG3_bought*-1*Sbase,2),'  \n')
print('Energy bought by the market         MWh/day ', round(EG1_bought*-1*Sbase,2),'  \n')
print('Perdidas Red                        MWh/day ', round(Perdi_tot*Sbase,2),'  \n')
print('Total Consumido                     MWh/day ', round((ECost/10 + Perdi_tot -EG1_bought-EG3_bought)*Sbase,2),'  \n')
print('  \n')

# print('COMMUNITY SW                         Eur/day %6.2f \n',SWelfare)
print('*******************************************************')