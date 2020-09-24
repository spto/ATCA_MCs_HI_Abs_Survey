import aplpy
import numpy as np

import csv
csvpath="/Users/boyangliu/Dropbox/ATCA_SMC/190624/fitting_190624/"
with open(csvpath+'Table1_190624.csv', newline='') as csvfile:
    data = np.array(list(csv.reader(csvfile)))
line_count=len(data[:,0])-1

RA=np.zeros(line_count)
Dec=np.zeros(line_count)
Ts_aver=np.zeros(line_count)
SNratio=np.zeros(line_count)
Ts_lowerlimit=np.zeros(line_count)

for n in np.arange(line_count):
    RAhms=np.array(data[n+1,2].split(':'))
    RAdeg=np.float(RAhms[0])*15+np.float(RAhms[1])*15/60+np.float(RAhms[2])*15/3600
    RA[n]=RAdeg

    Decdms=np.array(data[n+1,3].split(':'))
    Decdeg=np.float(Decdms[0])-np.float(Decdms[1])/60-np.float(Decdms[2])/3600
    Dec[n]=Decdeg

    Ts_aver[n]=data[n+1,11]
    SNratio[n]=np.float(data[n+1,7])/np.float(data[n+1,13]) #EW SN
    Ts_lowerlimit[n]=data[n+1,14]

sort_index=np.flip(np.argsort(SNratio),0)

RA=RA[sort_index]
Dec=Dec[sort_index]
Ts_aver=Ts_aver[sort_index]
SNratio=SNratio[sort_index]
Ts_lowerlimit=Ts_lowerlimit[sort_index]

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

fig.show_contour(path+'HIcol.fits', levels=[1E21], colors='r', linewidths=0.1, alpha=0.3)

for i in np.arange(line_count):
    xw=RA[i]
    yw=Dec[i]
    SNr=SNratio[i]
    Tsaver=Ts_aver[i]
    Tslowerlimit=Ts_lowerlimit[i]
    if SNr > 3:
        radius=(np.log10(SNr/3)+0.5)/7.5
        lgTs=np.log10(Tsaver)
        if Tsaver < 100:
            fig.show_circles(xw, yw, radius,linewidth=0.1,edgecolor='b',facecolor=[1.2-0.2*lgTs, 0, 2-lgTs])
        if (Tsaver > 100) & (Tsaver < 1000):
            fig.show_circles(xw, yw, radius,linewidth=0.1,edgecolor='b',facecolor=[0.8, lgTs-2, 0])
        if Tsaver > 1000:
            fig.show_circles(xw, yw, radius,linewidth=0.1,edgecolor='b',facecolor=[0.8, 1, lgTs-3])
    else:
        lgTslowerlimit=np.log10(Tslowerlimit)
        # print(i,Tsaver, Tslowerlimit)
        if Tslowerlimit < 100:
            fig.show_circles(xw, yw, 1/15.,linewidth=1,edgecolor='w',facecolor=[1.2-0.2*lgTslowerlimit, 0, 2-lgTslowerlimit])
        if (Tslowerlimit > 100) & (Tslowerlimit < 1000):
            fig.show_circles(xw, yw, 1/15.,linewidth=1,edgecolor='w',facecolor=[0.8, lgTslowerlimit-2, 0])
        if Tslowerlimit > 1000:
            fig.show_circles(xw, yw, 1/15.,linewidth=1,edgecolor='w',facecolor=[0.8, 1, lgTslowerlimit-3])

fig.show_circles(93, -73, (np.log10(3/3)+0.5)/7.5, edgecolor='black')
fig.show_circles(93, -73,(np.log10(30/3)+0.5)/7.5, edgecolor='black')
fig.show_circles(93, -73, (np.log10(300/3)+0.5)/7.5, edgecolor='black')
fig.add_label(0.11, 0.05, '$SNR_{EW}$: 3,30,300', relative=True, weight='light')

for i in np.arange(100):
    fig.show_rectangles(1525+i, 350, 3, 50, coords_frame='pixel', facecolor=[1-i/100.*0.2, 0, 1-i/100])
fig.add_label(0.795, 0.09, r'$10^{1}$', relative=True, weight='light')
fig.add_label(0.850, 0.09, r'$10^{2}$', relative=True, weight='light')
for i in np.arange(100):
    fig.show_rectangles(1625+i, 350, 3, 50, coords_frame='pixel', facecolor=[0.8, i/100., 0])
fig.add_label(0.905, 0.09, r'$10^{3}$', relative=True, weight='light')
for i in np.arange(100):
    fig.show_rectangles(1725+i, 350, 3, 50, coords_frame='pixel', facecolor=[0.8, 1, i/100.])
fig.add_label(0.96, 0.09, r'$10^{4}$', relative=True, weight='light')

fig.add_label(0.88, 0.06, r'$<T_{s}>$(K)', relative=True, weight='light')

outpath='/Users/boyangliu/Dropbox/ATCA_SMC/190624/maps/'
fig.save(outpath+'Map_Ts.pdf',dpi=300) #it takes 1129.708s to produce an eps! use png for quick view!
