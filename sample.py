import numpy as np
import matplotlib.pyplot as plt
import emcee
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
import sys

order=int(float(sys.argv[1])*10)
samp=np.loadtxt("num"+str(order)+".dat")
#plt.plot(samp)
#plot_acf(u,lags=40000)
#plt.show()
try:
	out=emcee.autocorr.integrated_time(samp, c=5, tol=50, quiet=False)
	with open("eqstep.txt","a+") as f:
		np.savetxt(f,np.array([[order,out]]))
	print(0)
	print(out)
	print(order)
except emcee.autocorr.AutocorrError:
	print(1)
