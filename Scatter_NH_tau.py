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

for n in np.arange(line_count):
    Scont[n]=data[n+1,4]
    tau_peak[n]=data[n+1,5]
    sigma_tau_peak[n]=data[n+1,6]
    SNratio[n]=np.float(data[n+1,7])/np.float(data[n+1,13])
    N_H_corr_iso[n]=np.float(data[n+1,8])*np.float(data[n+1,10])

sort_index=np.flip(np.argsort(Scont),0)

Scont=Scont[sort_index]
tau_peak=tau_peak[sort_index]
sigma_tau_peak=sigma_tau_peak[sort_index]
SNratio=SNratio[sort_index]
N_H_corr_iso=N_H_corr_iso[sort_index]

SNr_gt3=np.where(SNratio > 3)[0]
SNr_lt3=np.where(SNratio < 3)[0]
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

# the scatter plot:
z = Scont
# scaled_z = (z - z.min()) / z.ptp()
scaled_z = (z - 0.0) / 0.16
colors = plt.cm.hsv(scaled_z)

ax_scatter.errorbar(N_H_corr_iso[sigma_tau_peak_bad], tau_peak[sigma_tau_peak_bad], yerr=sigma_tau_peak[sigma_tau_peak_bad], fmt='none', color='black', alpha=0.3)
ax_scatter.scatter(N_H_corr_iso[SNr_gt3], tau_peak[SNr_gt3], facecolor=colors[SNr_gt3], s=50, alpha=0.5)
ax_scatter.scatter(N_H_corr_iso[SNr_lt3], tau_peak[SNr_lt3], facecolor=colors[SNr_lt3], edgecolor='black', s=50, alpha=0.3, marker='v')

ax_scatter.set_xlim((0, 9.99E21))
ax_scatter.set_ylim((0, 8))
ax_scatter.set_xlabel(r'$N_{H,cor,iso}(cm^{-2})$', size='large')
ax_scatter.set_ylabel(r'$\tau_{peak}$', size='large')

axins1 = inset_axes(ax_scatter, width="30%", height="3%", loc='upper center')
norm = mpl.colors.Normalize(vmin=0,vmax=0.16)
plt.colorbar(cm.ScalarMappable(norm=norm, cmap='hsv'), cax=axins1, orientation="horizontal", ticks=[0,0.08,0.16], alpha=0.5)
axins1.set_xlabel('$S_{cont}\ $(Jy)')

# the histogram
bins = np.arange(0, 10, 0.25)
ax_histy.hist(tau_peak, bins=bins, orientation='horizontal', edgecolor='black')
ax_histy.set_ylim(ax_scatter.get_ylim())
ax_histy.set_xlabel('counts',  size='large')

# plt.show()
outpath='/Users/boyangliu/Dropbox/ATCA_SMC/190624/hist_scatter/'
plt.savefig(outpath+'NH_tau.pdf',dpi=200)
