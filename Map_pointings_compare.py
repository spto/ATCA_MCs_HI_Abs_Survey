import aplpy
import numpy as np
import csv
csvpath="/Users/boyangliu/Dropbox/ATCA_SMC/190624/fitting_190624/"
with open(csvpath+'Dic94MZ2000Pointings.csv', newline='') as csvfile:
    data = np.array(list(csv.reader(csvfile)))
line_count=len(data[:,0])-1

Ref=np.zeros(line_count)
RAsource=np.zeros(line_count)
Decsource=np.zeros(line_count)
Scont=np.zeros(line_count)

for n in np.arange(line_count):
    Ref[n]=data[n+1,0]
    RAsource[n]=data[n+1,2]
    Decsource[n]=data[n+1,3]
    Scont[n]=data[n+1,4]

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
    print(i,Ref[i])
    xs=RAsource[i]
    ys=Decsource[i]
    sourceradius=(np.log10(Scont[i])-0.5)/30
    if Ref[i]==1994:
        fig.show_circles(xs, ys, sourceradius, linewidth=0, facecolor='red', alpha=1)
    elif Ref[i]==2000:
        fig.show_circles(xs, ys, sourceradius, linewidth=0, facecolor='blue', alpha=1)
    elif Ref[i]==2020:
        fig.show_circles(xs, ys, sourceradius, linewidth=0, facecolor='white', alpha=0.7)

# legend
fig.show_circles(95, -73, (np.log10(15)-0.5)/30, linewidth=0, facecolor='black')
fig.show_circles(94, -73.06, (np.log10(50)-0.5)/30, linewidth=0, facecolor='black')
fig.show_circles(93, -73.12, (np.log10(150)-0.5)/30, linewidth=0, facecolor='black')
fig.show_circles(92, -73.18, (np.log10(500)-0.5)/30, linewidth=0, facecolor='black')
fig.add_label(0.12, 0.06, '15,50,150,500 mJy', relative=True, weight='ultralight', size=9)
fig.add_label(0.12, 0.04, 'Observed continuum flux', relative=True, weight='ultralight', size=9)

fig.add_label(0.92, 0.08, 'D1994', relative=True, weight='ultralight', size=9, color='red')
fig.add_label(0.92, 0.06, 'M2000', relative=True, weight='ultralight', size=9, color='blue')
fig.add_label(0.92, 0.04, 'L2020', relative=True, weight='ultralight', size=9, color='white')

outpath='/Users/boyangliu/Dropbox/ATCA_SMC/190624/maps/'
fig.save(outpath+'Map_pointings_compare.pdf',dpi=300)
