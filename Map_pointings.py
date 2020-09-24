import aplpy
import numpy as np
import csv
csvpath="/Users/boyangliu/Dropbox/ATCA_SMC/190624/fitting_190624/"
with open(csvpath+'Map_pointings.csv', newline='') as csvfile:
    data = np.array(list(csv.reader(csvfile)))
line_count=len(data[:,0])-1

Name=np.zeros(line_count)
RAfield=np.zeros(line_count)
Decfield=np.zeros(line_count)
RAsource=np.zeros(line_count)
Decsource=np.zeros(line_count)
Scont=np.zeros(line_count)

for n in np.arange(line_count):
    Name=data[n+1,3]
    Scont[n]=data[n+1,8]

    RAfieldhms=np.array(data[n+1,4].split(':'))
    RAfielddeg=np.float(RAfieldhms[0])*15+np.float(RAfieldhms[1])*15/60+np.float(RAfieldhms[2])*15/3600
    RAfield[n]=RAfielddeg

    Decfielddms=np.array(data[n+1,5].split(':'))
    Decfielddeg=np.float(Decfielddms[0])-np.float(Decfielddms[1])/60-np.float(Decfielddms[2])/3600
    Decfield[n]=Decfielddeg

    RAsourcehms=np.array(data[n+1,6].split(':'))
    RAsourcedeg=np.float(RAsourcehms[0])*15+np.float(RAsourcehms[1])*15/60+np.float(RAsourcehms[2])*15/3600
    RAsource[n]=RAsourcedeg

    Decsourcedms=np.array(data[n+1,7].split(':'))
    Decsourcedeg=np.float(Decsourcedms[0])-np.float(Decsourcedms[1])/60-np.float(Decsourcedms[2])/3600
    Decsource[n]=Decsourcedeg

# sort_index=np.flip(np.argsort(Scont),0)

# RA=RA[sort_index]
# Dec=Dec[sort_index]
# Scont=Scont[sort_index]
# tau_peak=tau_peak[sort_index]
# sigma_tau_peak=sigma_tau_peak[sort_index]
# EW=EW[sort_index]
# EW_err=EW_err[sort_index]

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

fig.show_contour(path+'HIcol.fits', levels=[1E21], colors='r', linewidths=0.1,alpha=0.2)

for i in np.arange(line_count):
    print(i)

    xf=RAfield[i]
    yf=Decfield[i]
    fieldradius=17./60. #34arcmin beam
    figurewidth=2048*2./3600. #2048pixels, 2''/pixel
    fig.show_circles(xf, yf, fieldradius,linewidth=1,edgecolor='w')
    fig.show_circles(xf, yf, figurewidth/2,linewidth=0.1,edgecolor='w', alpha=0.5)
    # fig.show_rectangles(xf, yf, 4096./3600., 4096./3600., coords_frame='world', edgecolor='w', linewidth=0.1)

    xs=RAsource[i]
    ys=Decsource[i]
    sourceradius=(np.log10(Scont[i])-0.5)/30
    fig.show_circles(xs, ys, sourceradius, linewidth=0, facecolor='white')

# fig.show_circles(80.92054167, -70.83855556, 0.01, linewidth=0, facecolor='red')
# fig.show_circles(80.921, -70.85644444, 0.01, linewidth=0, facecolor='red')

# legend
fig.show_circles(95, -73, (np.log10(50)-0.5)/30, linewidth=0, facecolor='white')
fig.show_circles(94, -73.06, (np.log10(150)-0.5)/30, linewidth=0, facecolor='white')
fig.show_circles(93, -73.12, (np.log10(500)-0.5)/30, linewidth=0, facecolor='white')
fig.show_circles(92, -73.18, (np.log10(1500)-0.5)/30, linewidth=0, facecolor='white')
fig.add_label(0.12, 0.06, '50,150,500,1500 mJy', relative=True, weight='ultralight', size=9)
fig.add_label(0.12, 0.04, 'Catalogued 20cm flux', relative=True, weight='ultralight', size=9)

outpath='/Users/boyangliu/Dropbox/ATCA_SMC/190624/maps/'
fig.save(outpath+'Map_pointings.pdf',dpi=300)
