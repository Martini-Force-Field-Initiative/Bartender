#!/usr/bin/env python

#This one is for testing/internal use, not meant to be in the
#final product

import sys, os
import matplotlib.pyplot as plt
import numpy as np

itp=open("gmx_out.itp","r")

reading=False
for i in itp:
    if i.startswith(";") or len(i)<3:
        continue
    if "[dihedrals]" in i:
        reading=True
        continue
    if not reading:
        continue
    f=i.split()
    if f[4]!="1":
        continue
    beads=[int(f[0]),int(f[1]),int(f[2]),int(f[3])]
    ph=float(f[5])
    kd=float(f[6])
    n=round(float(f[7]))
    xvgname="Simple_periodic_%d-%d-%d-%d.xvg"%(beads[0],beads[1],beads[2],beads[3])
    xvg=open("xvg/"+xvgname,"r")
    x=[]
    ydist=[]
    ybt=[]
    ypy=[]
    ypydeg=[]
    for i in xvg:
        ir=i.lstrip()
        if ir.startswith("#") or ir.startswith("@"):
            continue
        f=i.split()
        if len(f)<3:
            print("Error!", i)
            sys.exit(1)
        x.append(float(f[0]))
        ydist.append(float(f[1]))
        ybt.append(float(f[2]))
        angle=np.radians(x[-1])
        eq=np.radians(ph)
        a=kd*(1+np.cos(n*angle-eq))
        ypy.append(a)
        ypydeg.append(kd*(1+np.cos(n*x[-1]-eq)))
        print("rad,deg",eq,ph)#########
    xvg.close()
    plt.plot(x,ydist,label="distribution")
    plt.plot(x,ybt,label="BT fit")
    plt.plot(x,ypy,"m",label="itp file")
    plt.plot(x,ypy,"c",label="itp file deg")
    plt.legend()
    plt.show()

