import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import csv
csvpath="/Users/boyangliu/Dropbox/ATCA_SMC/190624/fitting_190624/"
with open(csvpath+'Table2_190624.csv', newline='') as csvfile:
    data = np.array(list(csv.reader(csvfile)))
line_count=len(data[:,0])-1

m=np.zeros(line_count)
m_err=np.zeros(line_count)
Tew=np.zeros(line_count)
Tew_err=np.zeros(line_count)

for n in np.arange(line_count):
    m[n]=data[n+1,3]
    m_err[n]=data[n+1,4]
    Tew[n]=data[n+1,5]
    Tew_err[n]=data[n+1,6]

filter=np.where(np.absolute(m) > 3*m_err)[0]
print(len(m))
m=m[filter]
print(len(m))
Tew=Tew[filter]

left, width = 0.1, 0.35
bottom, height = 0.15, 0.7
spacing = 0.1

rect_1 = [left, bottom, width, height]
rect_2 = [left + width + spacing, bottom, 0.35, height]

plt.figure(figsize=(10, 3.5))

ax_1 = plt.axes(rect_1)
ax_2 = plt.axes(rect_2)

bins = np.arange(-30, 30, 5)
ax_1.hist(m, bins=bins, orientation='vertical', edgecolor='black')
ax_1.set_xlabel('slope (K)',  size='large')
ax_1.set_ylabel('N',  size='large')
ax_1.set_xticks(np.arange(7)*10-30)

bins = np.arange(5, 90, 5)
ax_2.hist(Tew, bins=bins, orientation='vertical', edgecolor='black')
ax_2.set_xlabel('$T_{ew}$ (K)',  size='large')
ax_2.set_ylabel('N',  size='large')
ax_2.set_xticks(np.arange(10)*10)
ax_2.set_yticks(np.arange(11))

# plt.show()
outpath='/Users/boyangliu/Dropbox/ATCA_SMC/190624/hist_scatter/'
plt.savefig(outpath+'m_Tew.pdf',dpi=300)
