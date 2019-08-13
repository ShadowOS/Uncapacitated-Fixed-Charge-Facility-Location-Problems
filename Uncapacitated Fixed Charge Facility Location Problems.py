# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 02:53:53 2019

@author: Karnika
"""


from gurobipy import*
import os
import xlrd
import xlsxwriter
import numpy as np
from random import randint

p=50                        # Number of Chromosome
S=4                         # Number of Generation
Digit=[]
chro=[]
chromosome=[]
Facility=[]
Demand={} 
distance = {}
fixed_cost = {}
chromosome_genes={}


i=0
for i in  range(p):
    if i<=p:
        Digit.append(i)
        i = i + 1
chro.append("chromosome")

for i in range(p):
    if i<=p:
        a=str(Digit[i])
        ab= chro[0] +" "+ a
        chromosome.append(ab)
        i = i+1

book = xlrd.open_workbook(os.path.join("Fixed_charge_1.xlsx"))

sh = book.sheet_by_name("cost")

i = 1
while True:
    try:
        sp = sh.cell_value(i,0)
        Facility.append(sp)
        fixed_cost[sp] = sh.cell_value(i,1)
        Demand[sp]=sh.cell_value(i,2)
        i = i + 1   
    except IndexError:
        break
sh = book.sheet_by_name("distance")

i = 1
for P in Facility:
    j = 1
    for Q in Facility:
        distance[P,Q] = sh.cell_value(i,j)
        j += 1
    i += 1

m = Model("GA_FC")

m.modelSense=GRB.MINIMIZE

X_j = m.addVars(Facility,vtype=GRB.BINARY,name='X_j')
Y_ij = m.addVars(Facility,Facility,vtype=GRB.CONTINUOUS,name='Y_ij')

m.setObjective(sum(fixed_cost[j]*X_j[j] for j in Facility) + sum((Demand[i]*distance[i,j]*Y_ij[i,j]) for i in Facility for j in Facility))
    
for i in Facility:
    m.addConstr(sum(Y_ij[i,j] for j in Facility) == 1)

for i in Facility:
    for j in Facility:
        m.addConstr(Y_ij[i,j] <= X_j[j])

m.optimize()

for v in m.getVars():
    if v.x > 0.01:
        print(v.varName, v.x)
print('Objective:',m.objVal)