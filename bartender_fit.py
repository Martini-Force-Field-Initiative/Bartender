#!/usr/bin/env python

import numpy as np
import scipy.optimize



def hooke(x,eq,k):
    return 0.5 * k * (x-eq)**2

def simple_periodic(x,eq,k,n):
    return k*(1+np.cos(n*x - eq))

def ryck_belle(x,p0,p1,p2,p3,p4,p5):
    cos=np.cos
    psi=x-np.pi
    return p0+p1*cos(psi) + p1*(cos(psi)**2) + p3*(cos(psi)**3) + p4*(cos(psi)**4) + p4*(cos(psi)**5)

def reb(x,eq,k):		
    return 0.5 * k * np.power((np.cos(x)-eq), 2.0) * (1 / np.power(np.sin(x), 2))

def cosangle(x,eq,k):		
    return 0.5 * k * np.power((np.cos(x) - np.cos(eq)), 2)


fname="task.fit"


fin=open(fname,"r")

fittype=fin.readline().rstrip("\n")


td=[]

td.append(fin.readline().replace("]","").replace("[","")) # x

td.append(fin.readline().replace("]","").replace("[","")) # y

td.append(fin.readline()) #initial guess

fin.close()

x=[]
y=[]
param=[]


print("Y QUE WEAAAAAAAA")

xtext=td[0].split()

for i in xtext:
    x.append(float(i))

ytext=td[1].split()

for i in ytext:
    y.append(float(i))

paramtext=td[2].split()

for i in paramtext:
    param.append(float(i))

x = np.array(x)
y = np.array(y)

f = ""
print(fittype)
if fittype=="hooke":
    f=hooke
if fittype=="simple_periodic":
    f=simple_periodic
if fittype=="ryck_belle":
    f=ryck_belle
if fittype=="reb":
    f=reb
if fittype=="cosa":
    f=cosangle

popt,pcov=scipy.optimize.curve_fit(f,x,y,param)

res = (y - f(x, *popt))**2

rmsd=np.sqrt(res.mean())



foutname=fittype+".out"

fout=open(foutname,"w")
fout.write(str(popt).replace("\n","").replace("[","").replace("]","")+"\n")
fout.write(str(rmsd)+"\n")
fout.close()
