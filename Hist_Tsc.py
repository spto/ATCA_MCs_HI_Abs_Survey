import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from uncertainties import unumpy

import csv
csvpath="/Users/boyangliu/Dropbox/ATCA_SMC/190624/fitting_190624/"
with open(csvpath+'Table2_190624.csv', newline='') as csvfile:
    data = np.array(list(csv.reader(csvfile)))
line_count=len(data[:,0])-1

m=np.zeros(line_count)
m_err=np.zeros(line_count)
Tsc0=np.zeros(line_count)
Tsc1=np.zeros(line_count)
Tsc2=np.zeros(line_count)
Tsc0_err=np.zeros(line_count)
Tsc1_err=np.zeros(line_count)
Tsc2_err=np.zeros(line_count)

for n in np.arange(line_count):
    m[n]=data[n+1,3]
    m_err[n]=data[n+1,4]
    Tsc0[n]=data[n+1,9]
    Tsc0_err[n]=data[n+1,10]
    Tsc1[n]=data[n+1,7]
    Tsc1_err[n]=data[n+1,8]
    Tsc2[n]=data[n+1,11]
    Tsc2_err[n]=data[n+1,12]

Tsc0_unp=unumpy.uarray(Tsc0, Tsc0_err)
Tsc1_unp=unumpy.uarray(Tsc1, Tsc1_err)
Tsc2_unp=unumpy.uarray(Tsc2, Tsc2_err)

Tsc0_select=Tsc0[np.where((Tsc0>0) & (m>0) & (m>3*m_err))]
Tsc1_select=Tsc1[np.where((Tsc0>0) & (m>0) & (m>3*m_err))]
Tsc2_select=Tsc2[np.where((Tsc0>0) & (m>0) & (m>3*m_err))]
mean1=np.mean(Tsc0_select)
mean2=np.mean(Tsc1_select)
mean3=np.mean(Tsc2_select)
Ntotal=str(len(Tsc0_select))
print(Ntotal)
print(mean1,mean2,mean3)

plt.figure(figsize=(3.5, 5))

left, width = 0.15, 0.8
bottom, height = 0.65, 0.25
spacing = 0.0

rect_1 = [left, bottom, width, height]
rect_2 = [left, bottom-spacing-height, width, height]
rect_3 = [left, bottom-(spacing+height)*2, width, height]

ax_1 = plt.axes(rect_1)
plt.xticks([])
plt.title('$N_{total}=42$')
plt.ylabel('N')
plt.ylim(0,13)
ax_2 = plt.axes(rect_2)
plt.xticks([])
plt.ylabel('N')
plt.ylim(0,13)
ax_3 = plt.axes(rect_3)
plt.ylabel('N')
plt.xlabel('$T_{c} (K)$')
plt.ylim(0,13)

bins = np.arange(0, 100, 5)
ax_1.hist(Tsc0_select, bins=bins, orientation='vertical', edgecolor='black', facecolor='tab:blue')
ax_2.hist(Tsc1_select, bins=bins, orientation='vertical', edgecolor='black', facecolor='lightcoral')
ax_3.hist(Tsc2_select, bins=bins, orientation='vertical', edgecolor='black', facecolor='tab:grey')

ax_1.axvline(x=mean1,linewidth=2, color='blue',linestyle='dashed')
ax_1.axvline(x=mean2,linewidth=2, color='red',linestyle='dashed')
ax_1.axvline(x=mean3,linewidth=2, color='black',linestyle='dashed')

ax_2.axvline(x=mean1,linewidth=2, color='blue',linestyle='dashed')
ax_2.axvline(x=mean2,linewidth=2, color='red',linestyle='dashed')
ax_2.axvline(x=mean3,linewidth=2, color='black',linestyle='dashed')

ax_3.axvline(x=mean1,linewidth=2, color='blue',linestyle='dashed')
ax_3.axvline(x=mean2,linewidth=2, color='red',linestyle='dashed')
ax_3.axvline(x=mean3,linewidth=2, color='black',linestyle='dashed')

ax_1.text(s='q=0.25', x=80, y=11)
ax_2.text(s='q=0.50', x=80, y=11)
ax_3.text(s='q=0.75', x=80, y=11)

outpath='/Users/boyangliu/Dropbox/ATCA_SMC/190624/hist_scatter/'
plt.savefig(outpath+'Tsc.pdf',dpi=300)

plt.show()
