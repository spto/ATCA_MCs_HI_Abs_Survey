import numpy as np
import matplotlib.pyplot as plt
import csv
import argparse

from astropy.io import fits
from astropy.io import votable
from astropy.coordinates import SkyCoord, Angle
from astropy.wcs import WCS
from astropy.table import Table, Column
from astropy.io.votable.tree import Param,Info
from astropy.io.votable import from_table, writeto
from astropy.io.votable import parse_single_table
from astropy import units as u

def main():

    path = "/Users/boyangliu/Dropbox/ATCA_SMC/190624/"

    with open(path+'fitting_190624/Table2_190624.csv', mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        for row in csv_reader: #must align with csv_reader
            # if line_count == 0:
            #     print(f'Column names are {", ".join(row)}')
            #     line_count += 1
            # else:
            sourcename=row["sourcename"]
            peakveloinkms=float(row["peakveloinkms"])
            pointpairs=int(row["pointpairs-1"])
            print('\t',sourcename,' at ',peakveloinkms,' km/s has ',pointpairs,' pairs.')

            votable=sourcename+"_opacity.votable.xml"
            votablefile=path+votable

            table = parse_single_table(votablefile)
            velo=table.array['velocity']/1000.
            abs=1-table.array['smooth_opacity']
            abserr=table.array['sigma_tau_smooth']    #Noise in the absorption profile, after smoothing
            em=table.array['em_mean']
            emerr=table.array['em_std']

            # Get current plot area size
            fig_size = plt.rcParams["figure.figsize"]
            print("Default size:", fig_size)

            # Set plot area width and height
            fig_size[0] = 18
            fig_size[1] = 5
            plt.rcParams["figure.figsize"] = fig_size
            print("Changed size:", fig_size)

            plt.figure()

            # subplot 1
            plt.subplot(1, 3, 1)
            plt.plot(velo, abs, 'o-', color='black', linewidth=1, markersize=2)
            plt.xlabel(r'V$_{LSRK}$ (km s$^{-1})$')
            plt.ylabel(r'$1-e^{-\tau}$')
            tau_max_smooth = + abserr # 3sigma!
            tau_min_smooth = - abserr
            plt.fill_between(velo, tau_min_smooth, tau_max_smooth, facecolor='cornflowerblue', color='cornflowerblue', zorder=3, alpha=0.5)
            tau_max_smooth = + 3*abserr # 3sigma!
            tau_min_smooth = - 3*abserr
            plt.fill_between(velo, tau_min_smooth, tau_max_smooth, facecolor='cornflowerblue', color='cornflowerblue', zorder=3, alpha=0.3)
            tau_max_smooth = + 5*abserr # 3sigma!
            tau_min_smooth = - 5*abserr
            plt.fill_between(velo, tau_min_smooth, tau_max_smooth, facecolor='cornflowerblue', color='cornflowerblue', zorder=3, alpha=0.1)
            plt.axvline(x=peakveloinkms,linewidth=1, color='r',linestyle='dashed')
            channelwidth=0.10305609955
            print(peakveloinkms-pointpairs*channelwidth,peakveloinkms+pointpairs*channelwidth)
            plt.axvspan(peakveloinkms-pointpairs*channelwidth,peakveloinkms+pointpairs*channelwidth, color='r',alpha=0.5)

            # subplot 2
            plt.subplot(1, 3, 2)
            plt.plot(velo, em, '-', color='black')
            plt.xlabel(r'V$_{LSRK}$ (km s$^{-1})$')
            plt.ylabel(r'$T_{B}$ (K)')
            em_max = em + emerr
            em_min = em - emerr
            plt.fill_between(velo, em_min, em_max, facecolor='lightgray', color='lightgray')
            plt.axvline(x=peakveloinkms,linewidth=1, color='r',linestyle='dashed')
            channelwidth=0.10305609955
            print(peakveloinkms-pointpairs*channelwidth,peakveloinkms+pointpairs*channelwidth)
            plt.axvspan(peakveloinkms-pointpairs*channelwidth,peakveloinkms+pointpairs*channelwidth, color='r',alpha=0.5)

            # subplot 3
            plt.subplot(1, 3, 3)
            # fig, ax = plt.subplots()
            peakveloinms=peakveloinkms*1000
            peakat = np.where(np.absolute(velo*1000-peakveloinms)<52)[0][0]
            pointpairs=pointpairs+1
            meanpoints=np.zeros((pointpairs,4))
            for n in range(pointpairs):
                meanabs=(abs[peakat+n]+abs[peakat-n])/2
                meanem=(em[peakat+n]+em[peakat-n])/2
                meanabserr=np.sqrt(np.power(abserr[peakat+n],2)+np.power(abserr[peakat-n],2))/2 #Not sure how to: the errors are not indepentdent with each other.
                meanemerr=np.sqrt(np.power(emerr[peakat+n],2)+np.power(emerr[peakat-n],2))/2
                meanpoints[n,0]=meanabs
                meanpoints[n,1]=meanem
                meanpoints[n,2]=meanabserr
                meanpoints[n,3]=meanemerr
            z = np.polyfit(meanpoints[:,0], meanpoints[:,1], 1)
            Tew=z[1]
            m=z[0]
            q=0.5
            Tsc=(1-q)*Tew+m
            xmax=np.max(abs[(peakat-pointpairs+1):(peakat+pointpairs)])
            xmin=np.min(abs[(peakat-pointpairs+1):(peakat+pointpairs)])
            ymax=np.max(em[(peakat-pointpairs+1):(peakat+pointpairs)])
            ymin=np.min(em[(peakat-pointpairs+1):(peakat+pointpairs)])
            plt.scatter(abs[(peakat-pointpairs+1):(peakat+pointpairs)], em[(peakat-pointpairs+1):(peakat+pointpairs)],color='blue')
            plt.plot([xmin,xmax], [z[0]*xmin+z[1],z[0]*xmax+z[1]])
            plt.errorbar(meanpoints[:,0], meanpoints[:,1], xerr=meanpoints[:,2], yerr=meanpoints[:,3], fmt='o',color='black',alpha=0.5)
            # plt.text(0.80,0.95, 'm='+str(np.around(z[0], decimals=1)),transform=ax.transAxes)
            # plt.text(0.80,0.90, 'Tew='+str(np.around(z[1], decimals=1))+'K',transform=ax.transAxes)
            # plt.text(0.80,0.85, 'Tsc='+str(np.around(Tsc, decimals=1))+'K',transform=ax.transAxes)
            plt.xlabel(r'Smoothed Opacity')
            plt.ylabel(r'Emission (K)')

            plt.suptitle(sourcename)
            filename=path+'components/'+sourcename+'_'+str(peakveloinkms)+'_'+str(pointpairs)+'.pdf'
            plt.savefig(filename)

            line_count += 1
        print(f'Processed {line_count} lines.')

if __name__ == "__main__":
    exit(main())
