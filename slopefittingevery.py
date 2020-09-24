#This fitting uses all the channels of the smoothed absorption spectra.

from __future__ import print_function, division

import argparse
import os
import sys
import time
import csv
import subprocess
from astropy.io import fits
from astropy.io import votable
from astropy.coordinates import SkyCoord, Angle
from astropy.wcs import WCS
from astropy.table import Table, Column
from astropy.io.votable.tree import Param,Info
from astropy.io.votable import from_table, writeto
from astropy.io.votable import parse_single_table
from astropy import units as u
import matplotlib.pyplot as plt
import math
import numpy as np
import numpy.core.records as rec
from scipy.interpolate import interp1d
from string import Template
import aplpy
from astropy.convolution import convolve, Gaussian1DKernel, CustomKernel
from scipy.odr import *
from uncertainties import ufloat

def main():
    path = "/Users/boyangliu/Dropbox/ATCA_SMC/190624/"
    #input, as Table2_input.csv

    source_name=['0439-6722_src_043856.17-672153.65','0439-6722_src_043856.17-672153.65','0445-6747_src_044451.16-674636.42','0450-6940_src_044949.02-694014.48','0452-6823_src_045221.02-682306.25','0452-6823_src_045221.02-682306.25','0452-6823_src_045221.02-682306.25','0454-7040_src_045332.96-704027.83','0456-6850_src_045618.66-684948.82','0456-6850_src_045618.66-684948.82','0456-6906_src_045426.16-691101.65','0456-6906_src_045541.96-690629.59','0456-6906_src_045551.56-690210.16','0456-6906_src_045551.56-690210.16','0506-6555_src_050531.12-655515.98','0507-6657_src_050501.47-664516.79','0507-6657_src_050529.68-665556.44','0507-6657_src_050529.68-665556.44','0507-6657_src_050831.19-670621.86','0510-6803_src_051129.60-680617.78','0514-6707_src_051340.07-670640.96','0514-6707_src_051340.07-670640.96','0514-6707_src_051537.51-672128.43','0514-6707_src_051537.51-672128.43','0514-6707_src_051537.51-672128.43','0514-6707_src_051537.51-672128.43','0514-6707_src_051537.51-672128.43','0517-7149_src_051641.61-714905.31','0517-7149_src_051641.61-714905.31','0517-7149_src_051641.61-714905.31','0517-7149_src_051641.61-714905.31','0521-6960_src_051832.77-693521.22','0521-6960_src_052105.17-695940.06','0521-6960_src_052105.17-695940.06','0522-7038_src_052229.47-703755.45','0522-7038_src_052340.83-705018.61','0527-6549_src_052626.21-655623.79','0527-6549_src_052626.21-655623.79','0527-6549_src_052631.85-654906.93','0527-6549_src_052631.85-654906.93','0527-6549_src_052631.85-654906.93','0528-6759_src_052635.36-674908.36','0528-6759_src_052635.36-674908.36','0528-6759_src_052635.36-674908.36','0528-6759_src_052635.36-674908.36','0528-6759_src_052745.92-675925.20','0533-7232_src_053344.75-721624.58','0536-6855_src_053527.99-691611.45','0536-6855_src_053527.99-691611.45','0536-6855_src_053537.53-685507.14','0536-6855_src_053537.53-685507.14','0536-6855_src_053537.53-685507.14','0536-6855_src_053604.76-691845.42','0538-6954_src_053804.67-695336.78','0538-6954_src_053804.67-695336.78','0538-6954_src_053804.67-695336.78','0538-6954_src_053937.53-694526.69','0538-6954_src_053946.00-693839.19','0538-6954_src_053946.00-693839.19','0538-6954_src_053946.00-693839.19','0538-6954_src_053946.00-693839.19','0538-6954_src_053946.00-693839.19','0538-6954_src_054004.51-694438.42','0538-6954_src_054004.51-694438.42','0538-6954_src_054004.51-694438.42','0538-6954_src_054004.51-694438.42','0538-6954_src_054004.51-694438.42','0540-6718_src_054010.40-671814.37','0540-6906_src_053657.09-691328.10','0540-6906_src_053657.09-691328.10','0540-6906_src_053657.09-691328.10','0540-6906_src_053657.09-691328.10','0543-6845_src_054314.55-684435.27','0543-6845_src_054314.55-684435.27','0551-6905_src_055130.16-691631.49','0552-6814_src_055205.83-681439.44','0552-6814_src_055205.83-681439.44']
    peak_velo_kms=[257.1,262.4,271.7,257.6,259.4,277.5,281.5,272.2,270,275.4,255.1,253.6,266.5,271.5,298.2,299.4,292.2,297.1,303.3,315.9,302.2,310.6,285.6,296.6,309,312.3,318.6,256.3,246.1,237.5,240.6,232.9,258.6,277,249.9,237.3,308.2,318.9,316.8,308.5,302,250.7,286.4,289.9,308.3,295.1,241,278.7,283.2,262.9,267.9,277.1,245.6,257.5,261.4,282.2,257,292.5,249.7,252.6,276.7,245.5,255.3,249.7,237.3,232.4,239.9,301,260.7,268.1,272.3,264.2,272.9,279.3,315.2,298.5,303.2]
    point_pairs=[28,11,10,18,18,29,10,22,15,17,39,28,13,15,26,23,32,17,12,11,32,19,30,22,15,15,15,15,32,15,14,17,39,25,16,13,17,16,36,35,27,17,10,23,37,22,19,31,9,31,14,25,12,25,10,9,69,11,10,19,21,27,20,20,9,9,17,22,17,18,23,17,29,23,11,20,19]

    sort_index=np.argsort(source_name)
    # print(sort_index)

    for n in np.arange(len(point_pairs)):
        fullname=source_name[sort_index[n]]
        sourcename=source_name[sort_index[n]][14:22]+source_name[sort_index[n]][23:30]
        peakveloinkms=peak_velo_kms[sort_index[n]]
        pairsnumber=point_pairs[sort_index[n]]              #excluding the peak point
        pointpairs=point_pairs[sort_index[n]]+1             #including the peak point
        vstart=peakveloinkms-0.103*pairsnumber
        vend=peakveloinkms+0.103*pairsnumber

        votable=fullname+"_opacity.votable.xml"
        votablefile=path+votable
        #print(votablefile)
        table = parse_single_table(votablefile)
        velo=table.array['velocity']
        abs=1-table.array['smooth_opacity']
        abserr=table.array['sigma_tau_smooth']    #Noise in the absorption profile, after smoothing
        em=table.array['em_mean']
        emerr=table.array['em_std']

        #print("Find peak location")
        peakveloinms=peakveloinkms*1000
        peakat = np.where(np.absolute(velo-peakveloinms)<52)[0][0]

        meanabs=np.zeros(pointpairs)
        meanem=np.zeros(pointpairs)
        meanabserr=np.zeros(pointpairs)
        meanemerr=np.zeros(pointpairs)
        for n in range(pointpairs):
            meanabs[n]=(abs[peakat+n]+abs[peakat-n])/2
            meanem[n]=(em[peakat+n]+em[peakat-n])/2
            meanabserr[n]=np.sqrt(np.power(abserr[peakat+n],2)+np.power(abserr[peakat-n],2))/2 #Not sure how to: the errors are not indepentdent with each other.
            meanemerr[n]=np.sqrt(np.power(emerr[peakat+n],2)+np.power(emerr[peakat-n],2))/2

        q0=0.5
        q1=0.25
        q2=0.75
        def lin_func(p, x):
             m, Tew = p
             return m*x + Tew
        lin_model = Model(lin_func)
        data = RealData(meanabs, meanem, sx=meanabserr, sy=meanemerr)
        odr = ODR(data, lin_model, beta0=[0., 1.])
        out = odr.run()
        m=out.beta[0]
        Tew=out.beta[1]
        m_err=out.sd_beta[0]
        Tew_err=out.sd_beta[1]
        Tsc0=(1-q0)*Tew+m
        Tsc0_err=np.sqrt(np.power((1-q0)*Tew_err,2)+np.power(m_err,2))
        Tsc1=(1-q1)*Tew+m
        Tsc1_err=np.sqrt(np.power((1-q1)*Tew_err,2)+np.power(m_err,2))
        Tsc2=(1-q2)*Tew+m
        Tsc2_err=np.sqrt(np.power((1-q2)*Tew_err,2)+np.power(m_err,2))

        m_unp=ufloat(m, m_err)
        tew_unp=ufloat(Tew, Tew_err)
        tsc_unp=ufloat(Tsc0, Tsc0_err)
        peakabs=ufloat(abs[peakat], abserr[peakat])
        # print(sourcename)

        with open("/Users/boyangliu/Dropbox/ATCA_SMC/190624/fitting_190624/Summary20190624_v2.txt", 'a') as f:
            print(sourcename,peakveloinkms,pointpairs-1,m,m_err,Tew,Tew_err,Tsc0,Tsc0_err,Tsc1,Tsc1_err,Tsc2,Tsc2_err,abs[peakat],abserr[peakat], file=f)  # Python 3.x

        # make Table2 for LaTex
        if (m+m_err) > 0:
            print(sourcename,'&', format(peakveloinkms, '.1f'), '&', pointpairs, '&', format(vstart, '.1f'), '-', format(vend, '.1f'), '&', str(format(peakabs, '.2f')).replace("+/-", "$\pm$"), '&', str(format(m_unp, '.1f')).replace("+/-", "$\pm$"), '&', str(format(tew_unp, '.1f')).replace("+/-", "$\pm$"), '&', str(format(tsc_unp, '.1f')).replace("+/-", "$\pm$"), '\\\\')
        else:
            print(sourcename,'&', format(peakveloinkms, '.1f'), '&', pointpairs, '&', format(vstart, '.1f'), '-', format(vend, '.1f'), '&', str(format(peakabs, '.2f')).replace("+/-", "$\pm$"), '&', str(format(m_unp, '.1f')).replace("+/-", "$\pm$"), '&', str(format(tew_unp, '.1f')).replace("+/-", "$\pm$"), '&', '\\nodata', '\\\\')


        # xmax=np.max(abs[(peakat-pointpairs+1):(peakat+pointpairs)])
        # xmin=np.min(abs[(peakat-pointpairs+1):(peakat+pointpairs)])
        # ymax=np.max(em[(peakat-pointpairs+1):(peakat+pointpairs)])
        # ymin=np.min(em[(peakat-pointpairs+1):(peakat+pointpairs)])

        # fig, ax = plt.subplots()
        # plt.scatter(meanabs, meanem)
        # plt.plot([xmin,xmax], [m*xmin+Tew, m*xmax+Tew])
        # plt.errorbar(meanabs, meanem, xerr=meanabserr, yerr=meanemerr, fmt='ro')
        # plt.text(0.80,0.95, 'm='+str(np.around(m, decimals=1)),transform=ax.transAxes)
        # plt.text(0.80,0.90, 'Tew='+str(np.around(Tew, decimals=1))+'K',transform=ax.transAxes)
        # plt.text(0.80,0.85, 'Tsc='+str(np.around(Tsc0, decimals=1))+'K',transform=ax.transAxes)
        # # plt.text(xmin+(xmax-xmin)*0.75, ymin+(ymax-ymin)*0.9, 'Tsc_all='+str(np.around(Tsc1, decimals=1))+'K')
        # plt.xlabel(r'Smoothed Opacity')
        # plt.ylabel(r'Emission (K)')
        # plt.savefig(path+'fitting_190624/'+votable[0:32]+'_'+str(peakveloinkms)+'_'+str(pointpairs-1)+'allchannels.png',format='png')
        # #plt.show()

    print("Finished")

if __name__ == "__main__":
    exit(main())
