#!/usr/bin/env python

import os, sys
import matplotlib.pyplot as plt
import numpy as np

    
#Feel free to change this name to a less stupid one :-D
#this one does the procedure we decided on (I hope!) The "flooding" itself isdone by the xflood function
#this one just runs that for all datasets and then "trims" so we don't have a lot of empty bins 
#on the edges.
def histoflood(xs,ys,inc, ran):
    nx=np.arange(ran[0],ran[1],inc)
    nys=[]
    for i in range(len(xs)):
        nys.append(np.zeros(len(nx)))
    for xid,x in enumerate(xs):
        floodx(x,nx,ys[xid],nys[xid]) #this one does the actual flooding
    #we now need to trim the data, as our initial range is probably a lot larger than
    #the actual data range
    begin=-1
    for i,v in enumerate(nys[0]):
        for j in nys:
            if j[i]>0:
                begin=i
                break
        if begin!=-1:
            break
    if begin==-1:
        begin=0
    end=-2
    for i,v in enumerate(reversed(nys[0])):
        for j in nys:
            actualindex=len(ys[0])-i-1
            if j[actualindex]>0:
                end=actualindex
                break
        if end!=-2:
            break
    print("range", begin,end)
    nnys=[]
    if end==-2:
        for i,v in enumerate(nys):
            nnys.append(v[begin:])
        return nx[begin:],nnys
    for i,v in enumerate(nys):
        nnys.append(v[begin:end])
    return nx[begin:end],nnys


#the "flooding" function. It just assigns the value for a bin in the original histogram
#to the closest bin in the new histogram.
def floodx(x,nx,y,ny):
    for i,v in enumerate(x):
        prevnx=0
        currnx=0
        nxindex=0
        Found=False
        for nxid,nxv in enumerate(nx):
            if v<nxv:
                currnx=nxv
                nxindex=nxid
                Found=True
                break
          #  prevnxv=nxv
            prevnx=nxv
        if not Found:
            continue
        if np.abs(v-prevnx)<=np.abs(v-currnx):
            ny[nxindex-1]+=y[i]
        else:
            ny[nxindex]+=y[i]


##For the actual program!
#We first read the input file. That file contains N (>=2) lines. Each line contains 2 to 3 fields:
#a file name F, a column number C and a label for that dataset. We will process the files F and output the Cth column of 
# the corresponding file. The legend in the first line will be ignored, and that line will be labeled "Reference"
#The file/column in the first line will be assumed to be the reference
#We will do something pretty simple. 
# All the columns will be assumed to be Y-axes sharing the same x-axis (column 0 of the first file)



# 1. All  the columns will be plotted together

# 2. The absolute differences between the last N-1 columns, and the first one, will be plotted together

# 3. The integrated absolute differences between the last N-1 columns, and the first one, 
# will be printed to a file.

###
####CAREFUL HERE! This dictionary contains the bin size and range for each property.
######
### You might want to play with the values. Do consider that I set up the angles range intentionally
##broad because, depending on the fit, bartender plots some angles centered on 0 (so there are 
# negative angles) and some centered on 180 (so there are angles > 180). The Flooding function
#will later trim the histograms to the actual data range, so there should be no problems.
binsdata={"angle":[1,[-180,360]],"bond":[0.001,[0,2.5]],"dihe":[10,[-180,360]]}




show=True
if "--noshow" in sys.argv:
    show=False

glyphs=["ro","b-","g-","m-","k-","c-","k--","b^-","ro-","g.-","c:"]

if "--glyphs" in sys.argv:
    glyphs=sys.argv[sys.argv.index("--glyphs")+1].split()


finp=open(sys.argv[1],'r')



xs=[]
ys=[]
labels=[]
reads=0
#I'm a bit oldschool with my Python.
for i in finp:
    fields=i.split()
    fdata=open(fields[0],'r')
    c=int(fields[1])
    xs.append([])
    ys.append([])
    for i in fdata:
        if i.startswith("#") or i.startswith("@"):
            continue #ignore comments and xvg data
        fields2=i.split()
        if len(fields2)<=c:
            raise Exception("File "+fields[0]+" does not have enough columns")
        xs[-1].append( float(fields2[0]))
        ys[-1].append(float(fields2[c]))
    if len(fields)>2:
        labels.append(fields[2])  
    else:
        labels.append("data "+str(reads))   
    fdata.close()
finp.close()
xs=np.array(xs)
ys=np.array(ys)



fig, ax = plt.subplots()

inc=binsdata[sys.argv[2]][0]
ran=binsdata[sys.argv[2]][1]
x,ys=histoflood(xs,ys,inc, ran)

print("toplot", x, ys,labels)

#now we have the x axis in x,and all the y data in ys

#1. All  the columns will be plotted together

for i,v in enumerate(ys):
    ax.plot(x,v,glyphs[i],label=labels[i])
ax.legend(loc='lower right')
plt.savefig("Fits.png",dpi=600)

if show:
    plt.show()

# We now prepare the absolute differences

diffs=[]

for i in ys[1:]:
    diffs.append(np.abs(i-ys[0]))

fig, ax = plt.subplots()
#We plot the differences

for i,v in enumerate(diffs):
    ax.plot(x,v,glyphs[i+1],label=labels[i+1])
ax.legend(loc='lower right')
plt.savefig("Diffs.png",dpi=600)

if show:
    plt.show()

#We now put the integrated absolute differences in a file (they are actually sums)

fout=open("InterangedDiffs.txt",'w')

sums =  []
for i,v in enumerate(diffs):
    sums.append(np.sum(v))
    outstr=labels[i+1]+" "+str(sums[-1])+"\n"
    fout.write(outstr)
    print(outstr.replace("\n",""))

fout.close()



#We now plot the integrated absolute differences
x=np.arange(len(sums))
plt.bar(x,sums) #one can change the plot type here.
plt.savefig("InterangedDiffs.png",dpi=600)

if show:
    plt.show()







