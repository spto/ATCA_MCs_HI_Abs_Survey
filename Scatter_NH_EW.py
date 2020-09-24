import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import csv
csvpath="/Users/boyangliu/Dropbox/ATCA_SMC/190624/fitting_190624/"
with open(csvpath+'N-EW.csv', newline='') as csvfile:
    data = np.array(list(csv.reader(csvfile)))
line_count=len(data[:,0])-1

Ref=np.zeros(line_count)
Nuncor=np.zeros(line_count)
Nuncor_err=np.zeros(line_count)
EW=np.zeros(line_count)
EW_err=np.zeros(line_count)

for n in np.arange(line_count):
    Ref[n]=data[n+1,0]
    Nuncor[n]=data[n+1,2]
    Nuncor_err[n]=data[n+1,3]
    EW[n]=data[n+1,4]
    EW_err[n]=data[n+1,5]

plt.figure(figsize=(7, 5))
slct=np.where(Ref==2020)
plt.errorbar(Nuncor[slct], EW[slct], xerr=Nuncor_err[slct], yerr=EW_err[slct], fmt=',',alpha=0.3, color='black')
slct=np.where(Ref==1994)
plt.errorbar(Nuncor[slct], EW[slct], xerr=Nuncor_err[slct], yerr=EW_err[slct], fmt='.',alpha=0.3, color='red')
slct=np.where(Ref==2000)
plt.errorbar(Nuncor[slct], EW[slct], xerr=Nuncor_err[slct], yerr=EW_err[slct], fmt='.',alpha=0.3, color='blue')

linex=np.arange(71)-5
liney=linex/(0.01832*10)
plt.plot(linex, liney, color='gray', alpha=0.3)
plt.text(s='10 K', x=0.5, y=25)
linex=np.arange(71)-5
liney=linex/(0.01832*100)
plt.plot(linex, liney, color='gray', alpha=0.3)
plt.text(s='100 K', x=38, y=24)
linex=np.arange(71)-5
liney=linex/(0.01832*1000)
plt.plot(linex, liney, color='gray', alpha=0.3)
plt.text(s='$<T_{s}>=$1000 K', x=50, y=3.5)

plt.plot(linex, 0*linex, color='gray', alpha=0.1)
plt.text(s='EW=0', x=58, y=-1.5, alpha=0.2)

plt.text(s='D1994', x=55, y=-5, color='red')
plt.text(s='M2000', x=55, y=-7, color='blue')
plt.text(s='L2020', x=55, y=-9, color='black')

plt.xlim(0, 65)
plt.ylim(-10, 30)
plt.xlabel(r'$N_{H,uncor}\ (10^{20} cm^{-2})$', size='large')
plt.ylabel('EW $(km\ s^{-1})$', size='large')

# plt.show()
outpath='/Users/boyangliu/Dropbox/ATCA_SMC/190624/hist_scatter/'
plt.savefig(outpath+'NH_EW.pdf',dpi=300)
