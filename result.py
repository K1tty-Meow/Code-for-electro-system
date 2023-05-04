import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd


tau=np.loadtxt("eqstep.txt")
Rtau=tau[np.argsort(-tau[:,0])]
xi=np.loadtxt("sequ1.txt")
#x=np.arange(8.0,0.8,-0.2)
x=xi[np.argsort(-xi)]
#x=Rtau[:,0]/10
print(x)
nlen=len(x)
meanf=[]
sigmaf=[]
for i in range(nlen):
	#print(Rtau[index,1])
	index=x[i]*10
	data=np.loadtxt("npmf"+str(int(index))+".dat")
	nsamp=len(data)
	f=0
	sigma=0
	for j in range(nsamp):
		f=f+data[j]
	meanf.append(f/nsamp)
	for k in range(nsamp):
		#index=int((80-i)/2)
		sigma=sigma+(data[k]-meanf[i])**2
	sigmaf.append(np.sqrt(sigma/nsamp/(nsamp-1)))


point=[]
pmf=[]
sigmapmf=[]
for i in range(nlen):
	point.append(meanf[i]-meanf[0])
for i in range(nlen):
	sumall=0
	sigpmf=0
	for k in range(i+1):
		sumall=sumall+point[k]*abs(x[k]-x[k-1])
		sigpmf=sigpmf+sigmaf[k]*abs(x[k]-x[k-1])
	pmf.append(sumall)
	sigmapmf.append(sigpmf)

x1=x[0:31]
pmf1=pmf[0:31]
sigmaf1=sigmaf[0:31]
#x=np.arange(1.0,8.2,0.2)
plt.errorbar(x,pmf,yerr=sigmaf,fmt='o',ecolor='r',color='b',elinewidth=2,capsize=4,ms=2)
plt.axhline(y=0.0,c="k",ls="--",lw=0.5)
plt.plot(x,pmf,label=r"$c_{1:1}=20mM, c_{3:3}=5mM$")
plt.xlabel(r"R/$\sigma$")
plt.ylabel(r"$\Delta U/k_{B}T$")
plt.xlim(1,8)
plt.title(r"Fixed $\gamma$=$0\degree$")
plt.legend()
plt.savefig("fmy0801.png")
print(meanf)
np.savetxt("pmf0801.txt",pmf)
np.savetxt("error0801.txt",sigmaf)
np.savetxt("meanf0801.txt",meanf)
#print(point)
#print(x)
#print(pmf)
#print(meanf)
