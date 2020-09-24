import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import csv
csvpath="/Users/boyangliu/Dropbox/ATCA_SMC/190624/fitting_190624/"
with open(csvpath+'Table1_190624.csv', newline='') as csvfile:
    data = np.array(list(csv.reader(csvfile)))
line_count=len(data[:,0])-1

Scont=np.zeros(line_count)
tau_peak=np.zeros(line_count)
sigma_tau_peak=np.zeros(line_count)
SNratio=np.zeros(line_count)
N_H_corr_iso=np.zeros(line_count)
Ts=np.zeros(line_count)
Ts_aver_err=np.zeros(line_count)
Ts_lowerlimit=np.zeros(line_count)

for n in np.arange(line_count):
    Scont[n]=data[n+1,4]
    # tau_peak[n]=data[n+1,5]
    # sigma_tau_peak[n]=data[n+1,6]
    # SNratio[n]=tau_peak[n]/sigma_tau_peak[n]
    N_H_corr_iso[n]=np.float(data[n+1,8])*np.float(data[n+1,10])
    Ts[n]=data[n+1,11]
    Ts_aver_err[n]=data[n+1,12]
    Ts_lowerlimit[n]=data[n+1,14]
    SNratio[n]=Ts[n]/Ts_aver_err[n]

sort_index=np.flip(np.argsort(Scont),0)

Scont=Scont[sort_index]
tau_peak=tau_peak[sort_index]
sigma_tau_peak=sigma_tau_peak[sort_index]
SNratio=SNratio[sort_index]
N_H_corr_iso=N_H_corr_iso[sort_index]
Ts=Ts[sort_index]
Ts_aver_err=Ts_aver_err[sort_index]
Ts_lowerlimit=Ts_lowerlimit[sort_index]

# SNr_gt3=np.where(SNratio > 3)[0]
# SNr_lt3=np.where(SNratio < 3)[0]
SNr_gt3=np.where(Ts_lowerlimit == 0)[0]
SNr_lt3=np.where(Ts_lowerlimit !=0)[0]
sigma_tau_peak_bad=np.where(sigma_tau_peak > 0.1)[0]

# definitions for the axes
left, width = 0.1, 0.6
bottom, height = 0.15, 0.7
spacing = 0.005

rect_scatter = [left, bottom, width, height]
rect_histy = [left + width + spacing, bottom, 0.2, height]

# start with a rectangular Figure
plt.figure(figsize=(8, 6))

ax_scatter = plt.axes(rect_scatter)
ax_scatter.tick_params(direction='in', top=True, right=True)
ax_histy = plt.axes(rect_histy)
ax_histy.tick_params(direction='in', labelleft=False)

# ax_scatter.errorbar(N_H_corr_iso[sigma_tau_peak_bad], Ts[sigma_tau_peak_bad], yerr=sigma_tau_peak[sigma_tau_peak_bad], fmt='none', color='black', alpha=0.3)
ax_scatter.scatter(N_H_corr_iso[SNr_gt3], Ts[SNr_gt3], s=50, alpha=0.5)
ax_scatter.scatter(N_H_corr_iso[SNr_lt3], Ts_lowerlimit[SNr_lt3], edgecolor='black', s=50, alpha=0.3, marker='^')
ax_scatter.errorbar(N_H_corr_iso[SNr_gt3], Ts[SNr_gt3], yerr=Ts_aver_err[SNr_gt3], fmt='none', color='black', alpha=0.3)

ax_scatter.set_xlim((0, 9.99E21))
ax_scatter.set_ylim((0, 1500))
ax_scatter.set_xlabel(r'$N_{H,cor,iso}(cm^{-2})$', size='large')
ax_scatter.set_ylabel(r'$<T_{s}>$(K)', size='large')

ax_scatter.axhline(y=299,linewidth=1, color='r',linestyle='dashed')
ax_scatter.axhline(y=197,linewidth=1, color='b',linestyle='dashed')
# ax_scatter.axhline(y=236,linewidth=1, color='r')
# ax_scatter.axhline(y=233.28,linewidth=1, color='b')

ax_scatter.axhline(y=1400,xmin=0.25,xmax=0.35,linewidth=1, color='r',linestyle='dashed')
ax_scatter.axhline(y=1350,xmin=0.25,xmax=0.35,linewidth=1, color='b',linestyle='dashed')
# ax_scatter.axhline(y=1600,xmin=0.25,xmax=0.35,linewidth=1, color='r')
# ax_scatter.axhline(y=1500,xmin=0.3,xmax=0.4,linewidth=1, color='b')

ax_scatter.text(s='Mean $<T_{s}>$ = 299 $\pm$ 8 K', x=3.7E21, y=1380)
ax_scatter.text(s='N weighted Mean $<T_{s}>$ = 197 $\pm$ 5 K', x=3.7E21, y=1330)
# ax_scatter.text(s='Median $<T_{s}>$ = 236 $\pm$ 10 K', x=3.7E21, y=1580)
# ax_scatter.text(s='N weighted  Median $<T_{s}>$ = 233 K', x=4.1E21, y=1480)

# the histogram
bins = np.arange(0, 2000, 100)
ax_histy.hist(Ts, bins=bins, orientation='horizontal', edgecolor='black')
ax_histy.set_ylim(ax_scatter.get_ylim())
ax_histy.set_xlabel('counts',  size='large')

# plt.show()
outpath='/Users/boyangliu/Dropbox/ATCA_SMC/190624/hist_scatter/'
plt.savefig(outpath+'NH_Ts_1500.pdf',dpi=200)
