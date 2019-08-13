# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 13:57:04 2019

@author: Karnika
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 20:46:02 2019

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
    
## GA module 
# GA paramater 
# this is cromosome  
Con1 = {}
N=0
###fitness of the cromosomes###
fit_crom_1 = 0
X_values = {}
Chromo_X=[]

def Fitness_Value_of_chromosome(): 
    k = 0
    for j in Facility:
        Con1[j] = m.addConstr(X_j[j] == Chromo_X[k])
        k = k + 1 
    m.update()    
    m.optimize()
    fit_crom_1 = m.objVal   
    for j in Facility:
        m.remove(Con1[j])  
    return(fit_crom_1)
    
def Chromosome():
    Chromo_X.clear()
    for j in Facility:
        a = np.random.random(1)[0]
        if a >= 0.5:
            b = 1
        else:
            b = 0
        Chromo_X.append(b)
    N=Fitness_Value_of_chromosome()
    Chromo_X.append(N)
    print('Chromo_X is ',Chromo_X)
    return(Chromo_X)

k=0 
   
c=Chromosome()  
l=c[-1]
ga=0
i=0
while i <= p:
    try:
        c=Chromosome()
        k=c[-1]
        if k < l:    
            x=chromosome[ga]
            g=0
            for j in Facility:
                d=c[g]
                chromosome_genes[x,j]=d
                g=g+1
            d=c[g]    
            chromosome_genes[x,'Fitness']=d
            ga=ga+1       
        else:
            i = i - 1
     
    except IndexError:
        break      
#### crossover ###

Chromo_1=[]
Chromo_2=[]
Chromo_3=[]
Chromo_4=[]
Chromo_5=[]
fa=0
fb=0

def Fitness_Value_of_chrom_3(): 
    k = 0
    for j in Facility:
        Con1[j] = m.addConstr(X_j[j] == Chromo_3[k])
        k = k + 1 
    m.update()    
    m.optimize()
    fit_crom_1 = m.objVal   
    for j in Facility:
        m.remove(Con1[j])  
    return(fit_crom_1)
def Fitness_Value_of_chrom_4(): 
    k = 0
    for j in Facility:
        Con1[j] = m.addConstr(X_j[j] == Chromo_4[k])
        k = k + 1 
    m.update()    
    m.optimize()
    fit_crom_1 = m.objVal   
    for j in Facility:
        m.remove(Con1[j])  
    return(fit_crom_1)

def crossover():         
    
    x = randint(0,50)   ##  x defines a crossover point 
    y= randint(0,p-1)   ##  y defines randomly selected chromosome 1
    z=randint(0,p-1)    ##  z defines randomly selected chromosome 2
    d=chromosome[y]
    e=chromosome[z]
    print(x,y,z)
    Chromo_1.clear()
    Chromo_2.clear()
    Chromo_3.clear()
    Chromo_4.clear()

    for i in Facility:
        
        Chromo_1.append(chromosome_genes[d,i])
        Chromo_2.append(chromosome_genes[e,i])
        
    for i in range(0, x):
        b = Chromo_1[i]
        Chromo_3.append(b)
        a = Chromo_2[i]
        Chromo_4.append(a)
    for i in range(x,50):
        b = Chromo_2[i]
        Chromo_3.append(b)
        a = Chromo_1[i]
        Chromo_4.append(a)
    
    print("the formulated Chromo_3 is")
    print(Chromo_3)   
    fa=Fitness_Value_of_chrom_3()
    print("fa value is ",fa)
    print("the formulated Chromo_4 is")
    print(Chromo_4)
    fb=Fitness_Value_of_chrom_4()
    print("fb value is ",fb)
    if fa <= fb:
        Chromo_5=Chromo_3
        Chromo_5.append(fa)
    else:
       Chromo_5=Chromo_4
       Chromo_5.append(fb)
    return(Chromo_5)

for x in range(S):
    k=0 
    c=crossover()
    l=c[-1]
    ga=0
    i=0
    while i <= p:
        try:
            c=crossover()
            k=c[-1]
            if k < l:    
                x=chromosome[ga]
                g=0
                for j in Facility:
                    d=c[g]
                    chromosome_genes[x,j]=d
                    g=g+1
                d=c[g]    
                chromosome_genes[x,'Fitness']=d
                ga=ga+1       
            else:
                i = i - 1
         
        except IndexError:
            break 
workbook=xlsxwriter.Workbook('Result.xlsx')
worksheet=workbook.add_worksheet('chromosome_genes')            
i=1
k=0
for x in chromosome:
    j=1
    l=0
    for y in Facility:
        e=chromosome_genes[x,y]
        worksheet.write(i,j,e)
        f=chromosome[k]
        worksheet.write(i,0,f)
        g=Facility[l]
        worksheet.write(0,j,g)
        l+=1
        j+=1
    e=chromosome_genes[x,'Fitness']
    worksheet.write(i,j,e)
    g='Fitness'
    worksheet.write(0,j,g)    
    i+=1
    k+=1
workbook.close() 
