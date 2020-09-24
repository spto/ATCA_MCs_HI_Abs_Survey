#Make Table1
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
from astropy.io.votable import parse
from astropy import units as u
import matplotlib.pyplot as plt
import math
import numpy as np
import numpy.core.records as rec

from scipy.interpolate import interp1d

from string import Template

import aplpy

from astropy.convolution import convolve, Gaussian1DKernel, CustomKernel

from uncertainties import unumpy
from uncertainties import ufloat

def main():
    path = "/Users/boyangliu/Dropbox/ATCA_SMC/190624/"
    plotpath=path#+"plot_files/"

    filelist=[]

    for root,dirs,files in os.walk(plotpath):
        for file in files:
            # change the extension from '.mp3' to
            # the one of your choice.
            if file.endswith('_region.png'):
                #print(str(file))
                filelist.append(str(file))

    print('Table 1_no_warning')
    print("Field", "Source_Name", "Scont", "tau_peak", "sigma_tau_peak", "EW", "EW_err", "N_H_uncor", "N_H_uncor_error", "f_H_cor_iso", "Ts_aver", "Ts_aver_err")
    for n in range(len(filelist)):
    #for n in [2]:
        Field=filelist[n][0:9]
        # Source_Name=filelist[n][14:33]
        Source_Name=filelist[n][14:22]+filelist[n][23:30]
        RAHH=filelist[n][14:16]
        RAMM=filelist[n][16:18]
        RASS=filelist[n][18:23]
        DecDD=filelist[n][23:26]
        DecMM=filelist[n][26:28]
        DecSS=filelist[n][28:33]
        RAHMS=RAHH+':'+RAMM+':'+RASS
        DecDMS=DecDD+':'+DecMM+':'+DecSS

        votablefile=path+filelist[n][0:33]+"_opacity.votable.xml"
        table = parse_single_table(votablefile)
        velo=table.array['velocity']
        smooth_opacity=table.array['smooth_opacity']    #exp(-tau)
        sigma_tau_smooth=table.array['sigma_tau_smooth']
        emission=table.array['em_mean']
        emission_error=table.array['em_std']

        votablefileoj=parse(votablefile)
        Scont=float(votablefileoj.infos[3].value)

        # if len(np.where(smooth_opacity < 0)[0]) >0:
            # print("For "+Source_Name+" the smoothed opacity is sometimes < 0!")
        tau=-1*np.log(smooth_opacity)
        tau_peak=max(tau)
        tau_peak_channel=np.where(tau == tau_peak)
        sigma_at_tau_peak=sigma_tau_smooth[tau_peak_channel[0][0]]
        sigma_tau_cont=sigma_tau_smooth[0]

        channel_width=(velo[1]-velo[0])/1000
        abs=1-smooth_opacity  #1-exp(-tau)
        EW=np.sum(abs)*channel_width
        smooth_opacity_err=sigma_tau_smooth
        abserr=smooth_opacity_err
        EWerr=np.sqrt(np.sum(abserr*abserr))*channel_width

        if math.isnan(emission[500]):
            print("No emission spectra found!")

        N_H_uncor=1.823*np.power(10, 18)*np.sum(emission)*channel_width
        N_H_uncor_error=1.823*np.power(10, 18)*np.sum(emission_error)*channel_width

        N_H_cor_iso=1.823*np.power(10, 18)*np.sum(emission*tau/abs)*channel_width
        f_H_cor_iso=N_H_cor_iso/N_H_uncor

        Ts_aver=np.sum(emission)*channel_width/EW

        Emierr=np.sqrt(np.sum(emission_error*emission_error))
        tauerr=np.sqrt(np.sum(sigma_tau_smooth*sigma_tau_smooth))
        sumemi=np.sum(emission)
        sumabs=np.sum(abs)
        Ts_aver_err=sumemi/sumabs*(Emierr/sumemi+tauerr/sumabs)

        emiwitherr=unumpy.uarray(emission, emission_error)
        abswitherr=unumpy.uarray(abs, abserr)
        Ts_aver_unp=np.sum(emiwitherr)/np.sum(abswitherr)
        # print(str(format(Ts_aver_unp, '.f')).replace("+/-", "$\pm$"))

        EW_unp=np.sum(abswitherr)*channel_width
        # print(Source_Name, str(format(EW_unp, '.f')).replace("+/-", "$\pm$"))

        N_H_uncor_unp=1.823*np.sum(emiwitherr)*channel_width/1000.                  #devided by 1E21
        # print(Source_Name, str(format(N_H_uncor_unp, '.f')).replace("+/-", "$\pm$"))

        tau_peak_unp=ufloat(tau_peak, sigma_at_tau_peak)
        # print(Source_Name, str(format(tau_peak_unp, '.f')).replace("+/-", "$\pm$"))

        # print(Field, format(Scont*1000, '.1f'))#, Source_Name, format(tau_peak_unp, '.f'), format(EW_unp, '.f'), format(N_H_uncor_unp, '.f'), format(Ts_aver_unp, '.1f'))
        # print(Source_Name, sigma_tau_cont)
        EW_unp_err=float(format(EW_unp, '.3f').split("+/-")[1])

        # Table 1 for LaTex:
        if EW_unp < (3*EW_unp_err) :
            Tslowerlimit=sumemi*channel_width/(3*EW_unp_err)
            # print(Field, '&', Source_Name, '&', format(Scont*1000, '.1f'), '&', str(format(tau_peak_unp, '.f')).replace("+/-", "$\pm$"), '&', str(format(EW_unp, '.f')).replace("+/-", "$\pm$"), '&', str(format(N_H_uncor_unp, '.f')).replace("+/-", "$\pm$"), '&', format(f_H_cor_iso, '.2f'), '&', '$> $'+str(format(Tslowerlimit,'.0f')), '\\\\')
            print(Field,  Source_Name, RAHMS, DecDMS, Scont, tau_peak, sigma_at_tau_peak, EW, N_H_uncor, N_H_uncor_error, f_H_cor_iso, Ts_aver, Ts_aver_err, EWerr, Tslowerlimit)

        else:
            # print(Field, '&', Source_Name, '&', format(Scont*1000, '.1f'), '&', str(format(tau_peak_unp, '.f')).replace("+/-", "$\pm$"), '&', str(format(EW_unp, '.f')).replace("+/-", "$\pm$"), '&', str(format(N_H_uncor_unp, '.f')).replace("+/-", "$\pm$"), '&', format(f_H_cor_iso, '.2f'), '&', str(format(Ts_aver_unp, '.f')).replace("+/-", "$\pm$"), '\\\\')
            print(Field,  Source_Name, RAHMS, DecDMS, Scont, tau_peak, sigma_at_tau_peak, EW, N_H_uncor, N_H_uncor_error, f_H_cor_iso, Ts_aver, Ts_aver_err, EWerr, 0)

        #print(Source_Name, smooth_opacity_err[0])

        # print(Source_Name, float(str(N_H_uncor_unp).split("+/-")[0])*10, float(str(N_H_uncor_unp).split("+/-")[1])*10, str(EW_unp).split("+/-")[0], str(EW_unp).split("+/-")[1])

        #fig = plt.figure()
        #plt.plot(velo,tau)
        #plt.show()
        #print(Field, '&', Source_Name, '&', RAHMS, '&', DecDMS, '&', format(Scont*1000, '.1f'), '&', format(tau_peak, '.1f'), '&', format(sigma_at_tau_peak, '.1f'), '&', format(EW, '.2f')+'$\pm$'+format(EWerr, '.2f'), '&', format(N_H_uncor/1E21, '.2f')+'$\pm$'+format(N_H_uncor_error/1E21, '.2f'), '&', format(f_H_cor_iso, '.2f'), '&', format(Ts_aver, '.1f')+'$\pm$'+format(Ts_aver_err, '.1f'), '\\\\')
        # print(Field, '&', Source_Name, '&', format(Scont*1000, '.1f'), '&', format(tau_peak, '.2f'), '&', format(sigma_at_tau_peak, '.2f'), '&', format(EW, '.2f')+'$\pm$'+format(EWerr, '.2f'), '&', format(N_H_uncor/1E21, '.2f')+'$\pm$'+format(N_H_uncor_error/1E21, '.2f'), '&', format(f_H_cor_iso, '.2f'), '&', format(Ts_aver, '.1f')+'$\pm$'+format(Ts_aver_err, '.1f'), '\\\\')
        # print(Field,  Source_Name,  format(Scont*1000, '.1f'),  format(tau_peak, '.2f'),  format(sigma_at_tau_peak, '.2f'),  format(EW, '.2f'),format(EWerr, '.2f'),  format(N_H_uncor/1E21, '.2f'),format(N_H_uncor_error/1E21, '.2f'),  format(f_H_cor_iso, '.2f'),  format(Ts_aver, '.1f'),format(Ts_aver_err, '.1f'))
        # print(Field,  Source_Name, RAHMS, DecDMS, Scont, tau_peak, sigma_at_tau_peak, EW, N_H_uncor, N_H_uncor_error, f_H_cor_iso, Ts_aver, Ts_aver_err, EWerr, Tslowerlimit)
        #print(Field, '&', Source_Name, '&', format(Scont*1000, '.1f'), '&', str(format(tau_peak_unp, '.f')).replace("+/-", "$\pm$"), '&', str(format(EW_unp, '.f')).replace("+/-", "$\pm$"), '&', str(format(N_H_uncor_unp, '.f')).replace("+/-", "$\pm$"), '&', format(f_H_cor_iso, '.2f'), '&', str(format(Ts_aver_unp, '.f')).replace("+/-", "$\pm$"), '\\\\')

    print("Finished")

if __name__ == "__main__":
    exit(main())
