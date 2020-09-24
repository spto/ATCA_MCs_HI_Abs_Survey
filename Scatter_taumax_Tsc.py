import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import csv
csvpath="/Users/boyangliu/Dropbox/ATCA_SMC/190624/fitting_190624/"

with open(csvpath+'Table1_190624.csv', newline='') as csvfile:
    data = np.array(list(csv.reader(csvfile)))
line_count=len(data[:,0])-1

IDtable1=[]
tau_peak=np.zeros(line_count)

for n in np.arange(line_count):
    IDtable1.append(data[n+1,1])
    tau_peak[n]=data[n+1,5]

with open(csvpath+'Table2_190624.csv', newline='') as csvfile:
    data = np.array(list(csv.reader(csvfile)))
line_count=len(data[:,0])-1

IDtable2=[]
for n in np.arange(line_count):
    m=float(data[n+1,3])
    m_err=float(data[n+1,4])
    if (m>0) & (m>3*m_err):
        IDtable2.append(np.array(data[n+1,0]))

print('length', len(IDtable2))
line_count=len(IDtable2)

Tsc0=np.zeros(line_count)
Tsc0_err=np.zeros(line_count)
for n in np.arange(line_count):
    Tsc0[n]=data[n+1,7]
    Tsc0_err[n]=data[n+1,8]

IDtable1=np.array(IDtable1)
IDtable2=np.array(IDtable2)

result=np.zeros((line_count,2))
result[:,0]=Tsc0
for n in np.arange(line_count):
    result[n,1]=tau_peak[np.where(IDtable1==IDtable2[n])[0]]

plt.figure(figsize=(7, 5))
# plt.scatter(result[:,1], result[:,0], s=50, alpha=0.5)
plt.errorbar(result[:,1], result[:,0], yerr=Tsc0_err, fmt='o',alpha=0.7)
plt.ylim(0, 50)
plt.xlabel(r'$\tau_{max}$', size='large')
plt.ylabel('$T_{c}(K)$', size='large')

# plt.show()
outpath='/Users/boyangliu/Dropbox/ATCA_SMC/190624/hist_scatter/'
plt.savefig(outpath+'taumax_Tsc.pdf',dpi=300)
