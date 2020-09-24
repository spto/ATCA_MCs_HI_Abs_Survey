import aplpy
import numpy as np
import csv
csvpath="/Users/boyangliu/Dropbox/ATCA_SMC/190624/fitting_190624/"
with open(csvpath+'Table1_190624.csv', newline='') as csvfile:
    data = np.array(list(csv.reader(csvfile)))
line_count=len(data[:,0])-1
print(line_count)

RA=np.zeros(line_count)
Dec=np.zeros(line_count)
Scont=np.zeros(line_count)
tau_peak=np.zeros(line_count)
sigma_tau_peak=np.zeros(line_count)
EW=np.zeros(line_count)
EW_err=np.zeros(line_count)

for n in np.arange(line_count):
    RAhms=np.array(data[n+1,2].split(':'))
    RAdeg=np.float(RAhms[0])*15+np.float(RAhms[1])*15/60+np.float(RAhms[2])*15/3600
    RA[n]=RAdeg

    Decdms=np.array(data[n+1,3].split(':'))
    Decdeg=np.float(Decdms[0])-np.float(Decdms[1])/60-np.float(Decdms[2])/3600
    Dec[n]=Decdeg

    Scont[n]=data[n+1,4]
    tau_peak[n]=data[n+1,5]
    sigma_tau_peak[n]=data[n+1,6]
    EW[n]=data[n+1,7]
    EW_err[n]=data[n+1,13]

sort_index=np.flip(np.argsort(Scont),0)

RA=RA[sort_index]
Dec=Dec[sort_index]
Scont=Scont[sort_index]
tau_peak=tau_peak[sort_index]
sigma_tau_peak=sigma_tau_peak[sort_index]
EW=EW[sort_index]
EW_err=EW_err[sort_index]

path='/Users/boyangliu/Dropbox/MAGMA/'
fig = aplpy.FITSFigure(path+'HIcol.fits')

fig.show_grayscale(invert=1)                #white background
fig.recenter(79.75,-69.5,radius=5)          #radius in degree
fig.add_grid()
fig.grid.set_linewidth(0.5)

fig.tick_labels.set_xformat('hh:mm')
fig.tick_labels.set_yformat('dd')

fig.add_colorbar()
fig.colorbar.show()
fig.colorbar.set_axis_label_text('$N_{H,uncor}(cm^{-2})$')

fig.show_contour(path+'HIcol.fits', levels=[1E21], colors='r', linewidths=0.1,alpha=0.3)

for i in np.arange(line_count):
    xw=RA[i]
    yw=Dec[i]
    radius=np.sqrt(Scont[i])/2
    print(Scont[i],radius)
    taupeak=tau_peak[i]
    sigmataupeak=sigma_tau_peak[i]
    EWvalue=EW[i]
    EWerr=EW_err[i]
    if EWvalue > EWerr*3:
        fig.show_circles(xw, yw, radius,linewidth=0.1,edgecolor='b',facecolor=[0.8, (np.log10(taupeak)+1.35)/2.3, 0])
    else:
        fig.show_circles(xw, yw, radius,linewidth=0.5,edgecolor='w',facecolor=[0.8, (np.log10(taupeak)+1.35)/2.3, 0])

fig.show_circles(93, -73, np.sqrt(1)/2, edgecolor='black')
fig.show_circles(93, -73, np.sqrt(0.316)/2, edgecolor='black')
fig.show_circles(93, -73, np.sqrt(0.1)/2, edgecolor='black')
fig.show_circles(93, -73, np.sqrt(0.0316)/2, edgecolor='black')
fig.show_circles(93, -73, np.sqrt(0.01)/2, edgecolor='black')
fig.add_label(0.12, 0.04, '0.01,0.032,0.1,0.32,1 Jy', relative=True, weight='light')
fig.add_label(0.12, 0.02, 'Source flux', relative=True, weight='light')

for i in np.arange(110):
    print(i,(np.log10(0.05*np.power(10,i/50.))+1.35)/2.3)
    fig.show_rectangles(1525+i*3, 350, 3, 50, coords_frame='pixel', facecolor=[0.8, (np.log10(0.05*np.power(10,i/50.))+1.35)/2.3, 0])
fig.add_label(0.8, 0.09, '0.05', relative=True, weight='light')
# fig.add_label(0.875, 0.09, '0.5', relative=True, weight='light')
fig.add_label(0.875+0.022577, 0.09, '1', relative=True, weight='light')
fig.add_label(0.95, 0.09, '5', relative=True, weight='light')
fig.add_label(0.97, 0.09, '8', relative=True, weight='light')
fig.add_label(0.885, 0.06, r'$\tau_{peak}$', relative=True, weight='light')

fig.add_label(086.08712, -69.48975, '30 Dor Region', color='white')
fig.add_label(072.66154, -69.34411, 'West Arm Region', color='white')
fig.add_label(079.49894, -66.79701, 'Northern Shells Region', color='white')
fig.add_label(082.84706, -72.78386, 'Magellanic Bridge', color='white')

outpath='/Users/boyangliu/Dropbox/ATCA_SMC/190624/maps/'
fig.save(outpath+'Map_sources.pdf',dpi=300)
