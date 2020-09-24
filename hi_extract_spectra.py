#!/usr/bin/env python -u

# Find sources in the ATCA continuum data and produce absorption and emission spectra
# for each suitable continuum source.

# Author James Dempsey
# Date 28 Aug 2016

# modified by Katie Jameson
# 1 Aug 2017

# mod by Boyang Liu 20190617
# For 1000channels fits in 20190601 (aftslfcl)
#Changed the weighting method

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
from astropy import units as u
import matplotlib.pyplot as plt
import math
import numpy as np
import numpy.core.records as rec

from scipy.interpolate import interp1d

from string import Template

import aplpy

from astropy.convolution import convolve, Gaussian1DKernel, CustomKernel

def main():
    """
    Main script for analyse_data
    :return: The exit code
    """
    dir_name = '/mnt/science1/bliu/ATCA/results20190601/'
    dir_name_mfs = '/mnt/science1/bliu/ATCA/results20190601/'

    # Read day parameter
    args = parseargs()
    field_name = args.field_name
    start = time.time()

    print ("#### Started source finding in field %s at %s ####" % \
          (field_name, time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))))
    error_list = []

    # Find the sources for the field
    if not args.extract_only:
        error_list.extend(find_sources(dir_name, field_name))

    # For each file, extract spectra
    produce_spectra(dir_name,field_name)

    # Report
    end = time.time()
    print ('#### Processing completed at %s ####' \
          % (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))))

    if len(error_list) == 0:
        print ("Hooray! No errors found.")
    else:
        print ("%d errors were encountered:" % (len(error_list)))
        for err in error_list:
            print (err)
    return 0

class IslandRange(object):
    def __init__(self, isle_id):
        self.isle_id = isle_id


def parseargs():
    """
    Parse the command line arguments
    :return: An args map with the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Find sources in the data for a field and produce spectra for each suitable source.")

    parser.add_argument("field_name", help="The field number to be analysed.")
    parser.add_argument("--extract_only", help="Use the previous source finding results to extract spectra", default=False,
                        action='store_true')

    args = parser.parse_args()
    return args

class CommandFailedError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

def run_os_cmd(cmd, failOnErr=True):
    """
    Run an operating system command ensuring that it finishes successfully.
    If the comand fails, the program will exit.
    :param cmd: The command to be run
    :return: None
    """
    #print >> cmd
    sys.stdout.flush()
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode != 0:
            message = "Command '"+cmd+"' failed with code " + str(retcode)
            #print >> sys.stderr, message
            if failOnErr:
                raise CommandFailedError(message)
    except OSError as e:
        message = "Command '" + cmd + "' failed " + e
        #print >> sys.stderr, message
        if failOnErr:
            raise CommandFailedError(message)
    return None


def find_sources(dir_name, field_name):
    """
    Search a continuum file for sources using the Aegean source finder. A
    VOTable file containing the list of discovered sources will be written out
    for the field. This function will use the Aegean source finder
    ( https://github.com/PaulHancock/Aegean ) to identify the sources.

    :param day_dir_name: The name of the day's directory.
    :param field_name:  The name fo the field to be searched for sources.
    :return: A list of error messages, if any
    """
    error_list = []

    dir_name_mfs = '/mnt/science1/bliu/ATCA/results20190601/'
    cont_file = dir_name_mfs + field_name + '/' + field_name + ".1419.5.2.mfs.restor_aftslfcl.fits"
    table_file = dir_name + 'extracted_spectra/' + field_name + "_src.vot"

    try:
        print ("##--## Searching continuum image " + cont_file + " ##--##")
        run_os_cmd('BANE ' + cont_file)
        aegean_cmd = 'aegean ' + cont_file + ' --autoload --telescope=ATCA ' \
                     '--cores=1 --island --table=' + table_file
        run_os_cmd(aegean_cmd)
    except CommandFailedError as e:
        error_list.append(str(e))
    return error_list


def read_sources(filename, dir_name, field_name):
    print ("Extracting sources from " + filename)
    sources = []
    dir_name_mfs = '/mnt/science1/bliu/ATCA/results20190601/'
    cont_file = dir_name_mfs + field_name + '/' + field_name + ".1419.5.2.mfs.restor_aftslfcl.fits"   #compare with line 142?
    cont_im = fits.open(cont_file)
    cont_im_hdr = cont_im[0].header
    cont_im_sz1 = np.int(cont_im_hdr['NAXIS1'])
    cont_im_sz2 = np.int(cont_im_hdr['NAXIS2'])

    if not os.path.exists(filename):
        print ("Warning: File %s does not exist, skipping source read." % \
               filename)
        return sources

    src_votable = votable.parse(filename, pedantic=False)
    results = src_votable.get_first_table().array
    for row in results:
        id = str(row['island']) + "-" + str(row['source'])
        ra = row['ra']
        dec = row['dec']
        ra_str = row['ra_str']
        ra_str_split = ra_str.split(":")
        dec_str = row['dec_str']
        dec_str_split = dec_str.split(":")
        id_name = ra_str_split[0]+ra_str_split[1]+ra_str_split[2]+dec_str_split[0]+dec_str_split[1]+dec_str_split[2]
        rms = row['local_rms']
        flux = row['peak_flux']
        sn = flux / rms
        print ("Found source %s at %.4f, %.4f with flux %.4f and rms of %.4f "
               "giving S/N of %.4f" % (id, ra, dec, flux, rms, sn))
        #determine if the source is near the edge of the image, if yes, exclude
        w = WCS(cont_im_hdr)
        pix = w.wcs_world2pix(ra, dec, 0, 0, 1)
        x_coord = int(np.round(pix[0])) - 1
        y_coord = int(np.round(pix[1])) - 1
        if sn > 10 and flux > 0.02 and 10 < x_coord < cont_im_sz1-10 and 10 < y_coord < cont_im_sz2-10:
            src = dict(zip(results.dtype.names,row))
            src['id'] = id_name
            src['sn'] = sn

            #sources.append([ra, dec, id, flux, row['island']])
            sources.append(src)

            #make image of the source
            name_prefix = field_name + '_src_' + src['id']
            src_img_name = dir_name + 'extracted_spectra/' + name_prefix + "_region.png"
            cont_im = fits.open(cont_file)
            cont_im_hdr = cont_im[0].header
            cdelt2 = cont_im_hdr['CDELT2']
            width = row['a']*cdelt2
            height = row['b']*cdelt2
            pa = row['pa']

            #load image, recenter to the source and zoom in
            src_im = aplpy.FITSFigure(cont_file)
            src_im.show_grayscale(vmin=-0.01,vmax=0.1)
            src_im.recenter(ra, dec, radius=0.025)
            print ("Plotting source ellipse of source %s with width = %0.4f deg and height = %0.4f deg" % (id,width,height))
            src_im.show_ellipses(ra, dec, width, height, angle=-pa ,edgecolor='red')
            src_im.add_colorbar()
            src_im.colorbar.show()
            src_im.colorbar.set_axis_label_text('MFS Flux (Jy/beam)')
            src_im.add_beam()
            src_im.show_beam()
            src_im.axis_labels.show()
            src_im.save(src_img_name)

        else:
            print ("Ignoring source at %.4f, %.4f due to low S/N of %.4f or "
                   "flux of %.4f" % (ra, dec, sn, flux))

    #create image of source

    return sources


def read_islands(filename):
    print ("Extracting islands from " + filename)
    islands = {}

    if not os.path.exists(filename):
        print ("Warning: File %s does not exist, skipping island read." % \
               filename)
        return {}

    isle_votable = votable.parse(filename, pedantic=False)
    results = isle_votable.get_first_table().array
    for row in results:
        islands[row['island']] = row
    return islands


def calc_island_ranges(islands, pixel_size):
    island_ranges = []
    for island in islands.values():
        ir = IslandRange(island['island'])
        ra = island['ra']
        dec = island['dec']
        ra_width = abs(island['x_width'] * pixel_size[0])
        dec_width = abs(island['y_width'] * pixel_size[1])
        ir.min_ra = ra - (ra_width/2)
        ir.max_ra = ra + (ra_width/2)
        ir.min_dec = dec - (dec_width/2)
        ir.max_dec = dec + (dec_width/2)
        print("Island %d goes from %f to %f (%d*%f)/ %f to %f (%d*%f)" % (
            island['island'], ir.min_ra, ir.max_ra, island['x_width'], pixel_size[0], ir.min_dec, ir.max_dec,
            island['y_width'], pixel_size[1]))
        island_ranges.append(ir)
    return island_ranges


def find_edges(fluxes, num_edge_chan):
    """
    Seek from the edges to find where the data starts for this set of fluxes.
    This accounts for an optional number of channels in the data which have no
    data recorded.
    :param fluxes: The array of fluxes to be checked.
    :param num_edge_chan: The number of edge channels with data to be skipped
    :return: The index of the first and last cell to have data.
    """

    l_edge = 0
    r_edge = len(fluxes)-1

    while fluxes[l_edge] == 0 and l_edge < len(fluxes)-1:
        l_edge += 1

    while fluxes[r_edge] == 0 and r_edge > 0:
        r_edge -= 1

    return l_edge + num_edge_chan, r_edge - num_edge_chan


def extract_spectra(dir_name, field_name):

    print ("Extracting spectra...")

    num_edge_chan = 10
    fits_filename = dir_name + field_name + "/" + field_name + ".1419.5.2.restor.fits"
    src_filename = dir_name + 'extracted_spectra/' + field_name + "_src_comp.vot"
    isle_filename = dir_name + 'extracted_spectra/' + field_name + "_src_isle.vot"

    spectra = dict()
    source_ids = dict()
    if not os.path.exists(fits_filename):
        print ("Warning: File %s does not exist, skipping extraction." % \
              fits_filename)
        return spectra, source_ids, []

    sources = read_sources(src_filename, dir_name, field_name)
    islands = read_islands(isle_filename)
    hdulist = fits.open(fits_filename)
    image = hdulist[0].data
    header = hdulist[0].header
    CDELT = np.abs(header['CDELT2'])
    w = WCS(header)
    index = np.arange(header['NAXIS3'])
    beam_maj = header['BMAJ'] * 60 * 60
    beam_min = header['BMIN'] * 60 * 60
    beam_size = beam_maj/CDELT
    beam_area = math.radians(header['BMAJ']) * math.radians(header['BMIN'])
    print ("Beam was %f x %f arcsec giving area of %f radians^2." % (beam_maj, beam_min, beam_area))
    ranges = calc_island_ranges(islands, (header['CDELT1'], header['CDELT2']))
    velocities = w.wcs_pix2world(10,10,index[:],0,0)[2]
    for src in sources:
        c = SkyCoord(src['ra'], src['dec'], frame='icrs', unit="deg")
        rms = get_channel_rms(image, w, src, beam_size)
        img_slice = get_integrated_spectrum(image, w, src)
        l_edge, r_edge = find_edges(img_slice, num_edge_chan)
        print("Using data range %d - %d out of %d channels." % (
            l_edge, r_edge, len(img_slice)))

        # plotSpectrum(np.arange(slice.size), slice)
        spectrum_array = rec.fromarrays(
            [np.arange(img_slice.size)[l_edge:r_edge],
             velocities[l_edge:r_edge],
             img_slice[l_edge:r_edge],
             rms[l_edge:r_edge]],
             names='plane,velocity,flux,rms')
        spectra[c.ra] = spectrum_array

        # isle = islands.get(src['island'], None)
        src_map = {'id': src['id'], 'flux': src['peak_flux'], 'pos': c, 'beam_area': beam_area, 'ra_str': src['ra_str'], 'dec_str': src['dec_str']}
        src_map['a'] = src['a']
        src_map['b'] = src['b']
        src_map['pa'] = src['pa']
        print (src_map)
        source_ids[c.ra] = src_map

        CDELT3=header['CDELT3']/1000.
        CRVAL3=header['CRVAL3']/1000.

    del image
    del header
    hdulist.close()

    return spectra, source_ids, ranges, CDELT3, CRVAL3


def get_weighting_array(data):
    """
    Calculate the mean of the continuum values. This is based on precalculated regions where there is no gas expected.
    :param data: A cubelet to be analysed, should be a 3D of flux values.
    :param planes: A umpy array of plane, and velocity values.
    :return: A 2D array of weighting values for the
    """

    channelno=np.arange(len(data[:,0,0]))
    continuum_range = np.ravel(np.where((channelno<260)|(channelno>740)))
    continuum_sample = data[continuum_range,:,:]
    mean_cont = np.nanmean(continuum_sample, axis=0)    #baseline flux for each point

    #REMOVE NANs that are due to negative mean continuum values
    #set to 0
    nans = np.where(np.isnan(mean_cont))
    mean_cont[nans] = 0.0
    sum_mean_cont_sq = np.sum(mean_cont ** 2)
    weighting = mean_cont / sum_mean_cont_sq * np.nanmax(mean_cont)

    print ("Got weighting of {} from {} and {}".format(weighting, mean_cont, sum_mean_cont_sq))
    return weighting


def get_integrated_spectrum(image, w, src):
    """
    Calculate the integrated spectrum of the component.
    :param image: The image's data array
    :param w: The image's world coordinate system definition
    :param src: The details of the component being processed
    :return: An array of average flux/pixel across the component at each velocity step
    """
    pix = w.wcs_world2pix(src['ra'], src['dec'], 0, 0, 1)
    x_coord = int(np.round(pix[0])) - 1  # 266
    y_coord = int(np.round(pix[1])) - 1  # 197
    print("Translated %.4f, %.4f to %d, %d" % (
        src['ra'], src['dec'], x_coord, y_coord))
    radius = 4
    y_min = y_coord - radius
    y_max = y_coord + radius
    x_min = x_coord - radius
    x_max = x_coord + radius
    data = np.copy(image[0, :, y_min:y_max+1, x_min:x_max+1])

    origin = SkyCoord(src['ra'], src['dec'], frame='icrs', unit="deg")
    pa_rad = math.radians(src['pa'])
    total_pixels = (y_max-y_min +1) * (x_max-x_min +1)
    outside_pixels = 0
    for i in range(x_min, x_max+1):
        for j in range(y_min, y_max+1):
            eq_pos = w.wcs_pix2world(i+1, j+1, 0, 0, 1)
            point = SkyCoord(eq_pos[0], eq_pos[1], frame='icrs', unit="deg")
            if not point_in_ellipse(origin, point, src['a'], src['b'], pa_rad):
                data[:, i-x_min, j-y_min] = 0
                outside_pixels += 1
                #if ((i<data.shape[1]) & (j<data.shape[2])):
                #    data[:, i-x_min, j-y_min] = 0
                #    outside_pixels += 1

    print("Found {} pixels out of {} inside the component {} at {} {}".format(total_pixels - outside_pixels, total_pixels,
                                                                       src['id'],
                                                                       point.ra,
                                                                       point.dec))
    weighting = get_weighting_array(data)
    integrated = np.sum(data * weighting, axis=(1, 2))
    inside_pixels = total_pixels - outside_pixels
    if inside_pixels <= 0:
        print ("Error: No data for component!")
    # else:
    #     integrated /= inside_pixels

    channelno=np.arange(len(integrated))
    continuum_range = np.where((channelno<260)|(channelno>740))
    x=np.ravel(continuum_range)
    print("shape(x)",np.shape(x))
    continuum_sample = integrated[continuum_range]
    y=np.ravel(continuum_sample)
    print("shape(y)",np.shape(y))
    mean_cont = np.nanmean(continuum_sample)
    print("mean_cont",mean_cont)
    z = np.polyfit(x, y, 1)
    print("baseline fitting result",z)
    baseline=channelno*z[0]+z[1]
    integrated=integrated-baseline+mean_cont

    return integrated

def get_channel_rms(image, w, src, beam_size):
    """
    Calculate the RMS of each channel.
    :param image: The image's data array
    :param w: The image's world coordinate system definition
    :param src: The details of the component being processed
    :return: An array of RMS at each velocity step
    """

    #get the x and y coordinates of the center of the source
    pix = w.wcs_world2pix(src['ra'], src['dec'], 0, 0, 1)
    x_center = int(np.round(pix[0])) - 1
    y_center = int(np.round(pix[1])) - 1

    #go 5 beam lengths away and calculate the RMS from a 10 x 10 pixel area
    x_offset = int(np.round(x_center + 5*beam_size))
    y_offset = int(np.round(y_center + 5*beam_size))
    x_min = x_offset - 5
    x_max = x_offset + 5
    y_min = y_offset - 5
    y_max = y_offset + 5

    data = image[0,:,y_min:y_max+1,x_min:x_max+1]

    rms = np.std(data,axis=(1,2))

    return rms

def get_mean_continuum(spectrum):
    """
    Calculate the mean of the continuum values. This is based on precalculated regions where there is no gas expected.
    :param spectrum: The spectrum to be analysed, should be a numpy array of
        plane, velocity and flux values.
    :return: A single float which is the mean continuum flux.
    """

    # velocities = spectrum.velocity
    # velo1=-123.063627344+(start_channel+50)*0.103056351495
    # velo2=-123.063627344+(start_channel+300)*0.103056351495
    # velo3=-123.063627344+(start_channel+1700)*0.103056351495
    # velo4=-123.063627344+(start_channel+1950)*0.103056351495
    # continuum_range = np.where(((velocities>velo1*1000.)&(velocities<velo2*1000.))|((velocities>velo3*1000.)&(velocities<velo4*1000.)))

    #
    # # print ("...gave sample of", continuum_sample)
    # mean_cont = np.nanmean(continuum_sample)

    print("len(spectrum.flux[:])",len(spectrum.flux[:]))
    channelno=np.arange(len(spectrum.flux[:]))
    continuum_range = np.where((channelno<260)|(channelno>740))
    x=np.ravel(continuum_range)
    print("shape(x)",np.shape(x))
    continuum_sample = spectrum.flux[continuum_range]
    y=np.ravel(continuum_sample)
    print("shape(y)",np.shape(y))
    mean_cont = np.nanmean(continuum_sample)
    print("mean_cont",mean_cont)

    opacity_sample = continuum_sample/mean_cont
    sd_cont = np.std(opacity_sample)

    #same with line 950
    hann_window = np.hanning(31)
    hann_kernel = CustomKernel(hann_window)
    smooth_opacity_sample = convolve(opacity_sample, hann_kernel , boundary='extend')
    sd_cont_smooth = np.nanstd(smooth_opacity_sample)

    print ("Found standard deviation of the smoothed spectrum of {}".format(sd_cont_smooth))

    return mean_cont, sd_cont, sd_cont_smooth


def get_opacity(spectrum, mean):
    """
    Calculate the opacity profile for the spectrum. This simply divides the
    spectrum's flux by the mean.

    :param spectrum: The spectrum to be processed
    :param mean: The mean background flux, representing what the backlighting sources average flux.
    :return: The opacity (e^(-tau)) at each velocity step.
    """
    # print spectrum.flux
    # print spectrum.flux/mean
    return spectrum.flux/mean


def get_temp_bright(spectrum, beam_area, wavelen=0.210996048):
    """
    Calculate the brightness temperature (T_B) for the spectrum. This effectively converts the spectrum from units
    of Jy/beam to K.

    :param spectrum: The spectrum to be processed
    :param beam_area: The beam area in radian^2
    :return: The brightness temperature at each velocity step.
    """
    k = 1.3806503E-23  # boltzmann constant in J K^-1
    jy_to_si = 1E-26  # J s^-1 m^-2 Hz^-1

    factor = (wavelen**2 / (2*k)) * jy_to_si / (np.pi*beam_area/4)
    print (factor)
    return factor * spectrum.flux

def plot_spectrum(x, y, filename, title, sigma_tau, sigma_tau_smooth):
    """
    Output a plot of opacity vs LSR velocity to a specified file.

    :param x: The velocity data
    :param y: The opacity values for each velocity step
    :param rms: array of channel RMS values
    :param filename: The file the plot should be written to. Should be
         an .eps or .pdf file.
    :param title: The title for the plot
    :param con_start_vel: The minimum velocity that the continuum was measured at.
    :param con_end_vel: The maximum velocity that the continuum was measured at.
    """
    fig = plt.figure()
    plt.plot(x/1000, y, color='gray', zorder=2)

    #overplot a smoothed spectrum
    #smooth the spectrum with hanning of 31 channels (1.6 km/s)
    hann_window = np.hanning(31)
    hann_kernel = CustomKernel(hann_window)
    smooth_y= convolve(y, hann_kernel, boundary='extend')
    plt.plot(x/1000, smooth_y, color='blue', zorder=4)

    if len(sigma_tau) > 0:
        tau_max = 1 + sigma_tau
        tau_min = 1 - sigma_tau
        plt.fill_between(x/1000, tau_min, tau_max, facecolor='silver', color='silver', alpha=0.75, zorder=1)
        tau_max_smooth = 1 + sigma_tau_smooth
        tau_min_smooth = 1 - sigma_tau_smooth
        plt.fill_between(x/1000, tau_min_smooth, tau_max_smooth, facecolor='cornflowerblue', color='cornflowerblue', zorder=3, alpha=0.5)

    plt.axhline(1, color='black')

    plt.xlabel(r'Velocity (km/s)')
    plt.ylabel(r'$e^{(-\tau)}$')
    plt.title(title)
    plt.grid(True)
    plt.savefig(filename)
    #plt.show()
    plt.close()
    return

def plot_original_spectrum(x, y, filename, title):
    """
    Output a plot of flux vs LSR velocity to a specified file.

    :param x: The velocity data
    :param y: The opacity values for each velocity step
    :param rms: array of channel RMS values
    :param filename: The file the plot should be written to. Should be
         an .eps or .pdf file.
    :param title: The title for the plot
    :param con_start_vel: The minimum velocity that the continuum was measured at.
    :param con_end_vel: The maximum velocity that the continuum was measured at.
    """
    fig = plt.figure()
    plt.plot(x/1000, y, color='gray', zorder=2)

    #overplot a smoothed spectrum
    #smooth the spectrum with hanning of 31 channels (1.6 km/s)
    hann_window = np.hanning(31)
    hann_kernel = CustomKernel(hann_window)
    smooth_y= convolve(y, hann_kernel, boundary='extend')
    plt.plot(x/1000, smooth_y, color='blue', zorder=4)

    # if len(sigma_tau) > 0:
    #     tau_max = 1 + sigma_tau
    #     tau_min = 1 - sigma_tau
    #     plt.fill_between(x/1000, tau_min, tau_max, facecolor='silver', color='silver', alpha=0.75, zorder=1)
    #     tau_max_smooth = 1 + sigma_tau_smooth
    #     tau_min_smooth = 1 - sigma_tau_smooth
    #     plt.fill_between(x/1000, tau_min_smooth, tau_max_smooth, facecolor='cornflowerblue', color='cornflowerblue', zorder=3, alpha=0.5)

    plt.xlabel(r'Velocity (km/s)')
    plt.ylabel(r'Flux (Jy)')
    plt.title(title)
    plt.grid(True)
    plt.savefig(filename)
    #plt.show()
    plt.close()
    return

def plot_emission_spectrum(velocity, em_mean, em_std, filename, title):
    """
    Output a plot of emission vs LSR velocity to a specified file.

    :param velocity: The velocity data
    :param em_mean: The mean temperature values for each velocity step
    :param em_std: The standard deviation in temperature values for each velocity step
    :param filename: The file the plot should be written to. Should be
         an .eps or .pdf file.
    :param title: The title for the plot
    :param con_start_vel: The minimum velocity that the continuum was measured at.
    :param con_end_vel: The maximum velocity that the continuum was measured at.
    """

    if len(em_mean) == 0:
        if os.path.exists(filename):
            os.remove(filename)
        return

    fig = plt.figure()
    plt.plot(velocity/1000, em_mean)

    em_max = em_mean + em_std
    em_min = em_mean - em_std
    plt.fill_between(velocity/1000, em_min, em_max, facecolor='lightgray', color='lightgray')

    plt.xlabel(r'Velocity relative to LSR (km/s)')
    plt.ylabel(r'$T_B$ (K)')
    plt.title(title)
    plt.grid(True)
    plt.savefig(filename)
    #plt.show()
    plt.close()
    return


def output_spectra(spectrum, opacity, filename, longitude, latitude, em_mean, em_std, temp_bright, beam_area,
                   sigma_tau, sigma_tau_smooth, mean):
    """
    Write the spectrum (velocity, flux and opacity) to a votable format file.

    :param spectrum: The spectrum to be output.
    :param opacity:  The opacity to be output.
    :param filename:  The filename to be created
    :param longitude: The galactic longitude of the target object
    :param latitude: The galactic latitude of the target object
    """
    table = Table(meta={'name': filename, 'id': 'opacity'})
    table.add_column(Column(name='plane', data=spectrum.plane))
    table.add_column(Column(name='velocity', data=spectrum.velocity, unit='m/s',description='velocity relative to LSRK'))
    table.add_column(Column(name='opacity', data=opacity))

    #following smooth should be using same parameter with line 950
    hann_window = np.hanning(31)
    hann_kernel = CustomKernel(hann_window)
    smooth_y= convolve(opacity, hann_kernel, boundary='extend')
    table.add_column(Column(name='smooth_opacity', data=smooth_y, description='opacity smooth with hanning window 31 channels.'))

    table.add_column(Column(name='flux', data=spectrum.flux, unit='Jy', description='Flux per beam'))
    table.add_column(Column(name='temp_brightness', data=temp_bright, unit='K'))
    table.add_column(Column(name='sigma_tau', data=sigma_tau, description='Noise in the absorption profile, before smoothing'))
    table.add_column(Column(name='sigma_tau_smooth', data=sigma_tau_smooth, description='Noise in the absorption profile, after smoothing'))
    if len(em_mean) > 0:
        # The emission may not be available, so don't include it if not
        table.add_column(Column(name='em_mean', data=em_mean, unit='K'))
        table.add_column(Column(name='em_std', data=em_std, unit='K'))

    votable = from_table(table)
    votable.infos.append(Info('ra', 'longitude', longitude.value))
    votable.infos.append(Info('dec', 'latitude', latitude.value))
    votable.infos.append(Info('beam_area', 'beam_area', beam_area))
    votable.infos.append(Info('cont', 'continuum', mean))
    writeto(votable, filename)


def output_emission_spectra(filename, longitude, latitude, velocity, em_mean,
                            em_std, ems):
    """
    Write the emission spectrum (velocity, flux and opacity) to a votable format
    file.

    :param filename:  The filename to be created
    :param longitude: The galactic longitude of the target object
    :param latitude: The galactic latitude of the target object
    :param velocity:
    :param em_mean:
    :param em_std:
    :param ems:
    """
    table = Table(meta={'name': filename, 'id': 'emission'})
    table.add_column(Column(name='velocity', data=velocity, unit='m/s',description='velocity relative to LSRK'))
    table.add_column(Column(name='em_mean', data=em_mean, unit='K'))
    table.add_column(Column(name='em_std', data=em_std, unit='K'))
    for i in range(len(ems)):
        table.add_column(Column(name='em_'+str(i), data=ems[i].flux, unit='K'))

    votable = from_table(table)
    votable.infos.append(Info('ra', 'longitude', longitude.value))
    votable.infos.append(Info('dec', 'latitude', latitude.value))
    writeto(votable, filename)


def point_in_ellipse(origin, point, a, b, pa_rad):
    # Convert point to be in plane of the ellipse
    p_ra_dist = point.icrs.ra.degree - origin.icrs.ra.degree
    p_dec_dist = point.icrs.dec.degree - origin.icrs.dec.degree
    x = p_ra_dist * math.cos(pa_rad) + p_dec_dist * math.sin(pa_rad)
    y = - p_ra_dist * math.sin(pa_rad) + p_dec_dist * math.cos(pa_rad)

    a_deg = a / 3600
    b_deg = a / 3600

    # Calc distance from origin relative to a/b
    dist = math.sqrt((x / a_deg) ** 2 + (y / b_deg) ** 2)
    print("Point %s is %f from ellipse %f, %f, %f at %s." % (point, dist, a, b, math.degrees(pa_rad), origin))
    return dist <= 1.0


def point_in_island(point, islands):
    ra = point.icrs.ra.degree
    dec = point.icrs.dec.degree
    for island in islands:
        if island.min_ra <= ra <= island.max_ra and island.min_dec <= dec <= island.max_dec:
            print("Point %s in island %d at %f, %f" % (point, island.isle_id, island.min_ra, island.min_dec))
            return True
    print("Point %f, %f not in any of %d islands" % (ra, dec, len(islands)))
    return False


def calc_offset_points(longitude, latitude, beam_size, a, b, pa, islands, num_points=6, max_dist=0.04):
    spacing = 2.0 * math.pi / float(num_points)
    origin = SkyCoord(longitude, latitude, frame='icrs', unit="deg")
    pa_rad = math.radians(pa)
    points = []
    for i in range(0, num_points):
        angle = spacing * i
        mult = 0.5
        inside_component = True
        while inside_component:
            if mult*beam_size > max_dist:
                coord = None
                break;
            g_l = longitude + math.sin(angle)*beam_size*mult
            g_b = latitude + math.cos(angle)*beam_size*mult
            coord = SkyCoord(g_l, g_b, frame='icrs', unit="deg")

            inside_component = point_in_ellipse(origin, coord, a, b, pa_rad) or point_in_island(coord, islands)
            mult += 0.5
        if coord is None:
            print("Point could not be found for angle %f within max dist of %f (mult %f)" % (
                math.degrees(angle), max_dist, mult))
        else:
            print ("Point at angle %f is %s with mult %f" % (math.degrees(angle), str(coord), mult-0.5))
            points.append(coord)

    return points

def calc_offset_points_pb(longitude, latitude, beam_size, a, b, pa, islands, num_points=84):
    spacing = 2.0 * math.pi / float(num_points)
    origin = SkyCoord(longitude, latitude, frame='icrs', unit="deg")
    pa_rad = math.radians(pa)
    points = []
    for i in range(0, num_points):
        angle = spacing * i
        mult = 0.5
        inside_component = True
        while inside_component:
            g_l = longitude + math.sin(angle)*beam_size*mult
            g_b = latitude + math.cos(angle)*beam_size*mult
            coord = SkyCoord(g_l, g_b, frame='icrs', unit="deg")

            inside_component = point_in_ellipse(origin, coord, a, b, pa_rad) or point_in_island(coord, islands)
            mult += 0.5
        if coord is None:
            print("Point could not be found for angle %f within max dist of %f (mult %f)" % (
                math.degrees(angle), max_dist, mult))
        else:
            print ("Point at angle %f is %s with mult %f" % (math.degrees(angle), str(coord), mult-0.5))
            points.append(coord)

    return points

class Em_Spectrum(object):
    """
    A comtainer for a spectrum
    """

    def __init__(self, coord, velocity, flux):
        self.coord = coord
        self.velocity = velocity
        self.flux = flux

    def __str__(self):
        return self.coord

def get_emission_spectra(centre, abs_velocities, dir_name, field_name, a, b, pa, islands):
    """
    Extract emission spectra within a beam around a central point.
    Also extract the average emission spectrum in the 34' ATCA primary beam
    to use to scale the noise in the absorption spectra.

    :param centre: A SkyCoord containing the location of the central point
    :param velocities: The velocities list so that the emission data can be matched.
    :param file_list: A list of dictionaries describing the SGPS files.
    :paeram a: semi-major axis length in arcsec of the component ellipse
    :paeram b: semi-minor axis length in arcsec of the component ellipse
    :paeram pa: parallactic angle of the component ellipse
    :return: An array fo the mean and standard deviation of emission at each velocity.
    """

    print ("Extracting emission spectra...")

    filename = dir_name + 'extracted_spectra/' + field_name + '_emission.votable.xml'
    fits_filename = '/mnt/science1/bliu/HINSA/data/HI.fits'

    #get the beam size of the emission data
    hdulist = fits.open(fits_filename)
    image = hdulist[0].data
    header = hdulist[0].header
    em_beam_size = header['BMAJ']

    #donut
    coords = list(set(calc_offset_points(centre.ra.value, centre.dec.value, 2*em_beam_size, a, b, pa, islands))-set(calc_offset_points(centre.ra.value, centre.dec.value, 1*em_beam_size, a, b, pa, islands)))

    #primary beam size
    pb_size = 34./60.
    coords_PB = calc_offset_points_pb(centre.ra.value, centre.dec.value, pb_size, a, b, pa, islands)

    #ems = sgps.extract_spectra(coords, file_list)
    #took this from sgps.extract
    w = WCS(header)
    index = np.arange(header['NAXIS3'])
    velocities = w.wcs_pix2world(0, 0, index[:], 0)[2]    #deleted the last 0, back to older version with 4 args

    ems = []
    ems_pb = []


    # For each target coord in this file
    for i in range(0,len(coords)):

        coord = coords[i]
        # convert desired coord to pixel coords
        pix = w.wcs_world2pix(coord.ra, coord.dec, 0, 1)   #deleted the last 0, back to older version with 4 args
        x_coord = int(np.round(pix[0])) - 1
        y_coord = int(np.round(pix[1])) - 1
        # print("Translated %.4f, %.4f to %d, %d" % (
        #     coord.ra.value, coord.dec.value, x_coord, y_coord))

        if (y_coord<image.shape[1]) & (x_coord<image.shape[2]):    #back to older version
            # Extract slice
            slice = image[:, y_coord, x_coord]     #back to older version
            print('y_coord,x_coord',y_coord,x_coord)
            print('type(slice)',type(slice))
            print('len(slice)',len(slice))
            print('slice[42:52]', slice[42:52])

            #create an interpolation so the spectrum can be matched
            #to the absorption spectrum resolution
            print('type(velocities)',type(velocities))
            print('len(velocities)',len(velocities))
            print('velocities[42:52]', velocities[42:52])

            print('type(abs_velocities)',type(abs_velocities))
            print('len(abs_velocities)',len(abs_velocities))
            print('abs_velocities[1000:1010]', abs_velocities[1000:1010])

            f_interp = interp1d(velocities, slice, bounds_error=False)
            match_slice = f_interp(abs_velocities)
            print('type(match_slice)',type(match_slice))
            print('len(match_slice)',len(match_slice))
            print('match_slice[1000:1010]', match_slice[1000:1010])

            # Create spectra object
            spectrum = Em_Spectrum(coord, abs_velocities, match_slice)
            print('type(spectrum.flux)',type(spectrum.flux))
            print('len(spectrum.flux)',len(spectrum.flux))
            print('emission spectrum [1000:1010]', spectrum.flux[1000:1010])

            # Add to result list
            ems.append(spectrum)
            print('after ems.append(spectrum): len(ems)=',len(ems))


        else:
            print("WARNING: Unable to find emission data for " + str(centre))

    for i in range(0,len(coords_PB)):

        coord = coords_PB[i]

        # convert desired coord to pixel coords
        pix = w.wcs_world2pix(coord.ra, coord.dec, 0, 1)    #back to older version
        x_coord = int(np.round(pix[0])) - 1
        y_coord = int(np.round(pix[1])) - 1
        print("Translated %.4f, %.4f to %d, %d" % (
            coord.ra.value, coord.dec.value, x_coord, y_coord))

        if (y_coord<image.shape[1]) & (x_coord<image.shape[2]):     #back to older version
            # Extract slice
            slice = image[:, y_coord, x_coord]         #back to older version

            #create an interpolation so the spectrum can be matched
            #to the absorption spectrum resolution
            f_interp = interp1d(velocities, slice, bounds_error=False)
            match_slice = f_interp(abs_velocities)

            # Create spectra object
            spectrum = Em_Spectrum(coord, abs_velocities, match_slice)

            # Add to result list
            ems_pb.append(spectrum)
            print('after ems_pb.append(spectrum): len(ems_pb)=',len(ems_pb))

    del image
    del header
    hdulist.close()


    print("Found {} emission points from {} coords for point RA={}, Dec={}".format(len(ems), len(coords),
                                                                                centre.ra.value,
                                                                                centre.dec.value))
    if ems:
        print('len(ems)',len(ems))
        all_em = np.array([ems[i].flux for i in range(len(ems))])
        print('type(all_em)',type(all_em))
        print('len(all_em)',len(all_em))
        print('all_em', all_em)
        em_std = np.std(all_em, axis=0)
        em_mean = np.mean(all_em, axis=0)
        print('type(em_mean)',type(em_mean))
        print('len(em_mean)',len(em_mean))
        print('em_mean', em_mean)
        # f_em_std_interp=interp1d(ems[0].velocity, em_std, kind='quadratic')
        # em_std_interp = f_em_std_interp(abs_velocities)
        # f_em_mean_interp=interp1d(ems[0].velocity, em_mean, kind='quadratic')
        # em_mean_interp = f_em_mean_interp(abs_velocities)

        output_emission_spectra(filename,
                                centre.ra, centre.dec,
                                ems[0].velocity, em_mean, em_std, ems)


        all_em_pb = np.array([ems_pb[i].flux for i in range(len(ems_pb))])
        em_pb_std = np.std(all_em_pb, axis=0)
        em_pb_mean = np.mean(all_em_pb, axis=0)
        # f_em_pb_mean_interp=interp1d(ems[0].velocity, em_pb_mean, kind='quadratic')
        # em_pb_mean_interp = f_em_pb_mean_interp(abs_velocities)

        output_emission_spectra(filename,centre.ra, centre.dec, ems_pb[0].velocity, em_pb_mean, em_pb_std, ems_pb)   #filename_pb not defined

        # return em_mean_interp, em_std_interp, em_pb_mean_interp
        return em_mean, em_std, em_pb_mean

    print("WARNING: Unable to find emission data for " + str(centre))
    if os.path.exists(filename):
        os.remove(filename)
    return [], [], []


def calc_sigma_tau(opacity, cont_sd, cont_sd_smooth, pb_em_mean):
   """
   Calculate the noise in the absorption profile at each velocity step. Where emission data is available, this is
   based on the increased antenna temperature due to the received emission.

   :param opacity: The optical depth spectrum, used only for the shape of the data
   :param cont_sd: The standard deviation of the opacity spectrum away from the line
   :param cont_sd_smooth: The standard deviation of the smoothed opacity spectrum away from the line
   :param em_mean: The mean emission brightness temperature in K for the primary beam of ATCA
   :return: A numpy array containing the noise level in the opacity data at each velocity step.
   """

   # Tsys of receiver
   # website says 44.7 K
   tsys = 44.7
   # the brightness temp of the HI emission in PB needs to be
   # scaled by the antenna efficiency (found on ATCA sensitivity calculator website)
   ant_eff = 0.50

   #scale the off-line noise by the emission Tb (see Murray+ 2015)
   #compared to the Tsys of the ATCA observations
   if len(pb_em_mean) > 0:
       floor = np.zeros(pb_em_mean.shape)
       sigma_tau = cont_sd * ((tsys + np.fmax(floor, pb_em_mean)*ant_eff) / tsys)
       sigma_tau_smooth = cont_sd_smooth * ((tsys + np.fmax(floor, pb_em_mean)*ant_eff) / tsys)
   else:
       sigma_tau = np.full(opacity.shape, cont_sd)
       sigma_tau_smooth = np.full(opacity.shape, cont_sd_smooth)
   return sigma_tau, sigma_tau_smooth

def plot_em_abs(em_mean, em_std, opacity, sigma_tau_smooth, filename, CDELT3, CRVAL3):
    # read in the emission and absorption spectra and
    # plot an emission-absortion diagram

    if len(em_mean) == 0:
        if os.path.exists(filename):
            os.remove(filename)
        return

    #smooth the spectrum with hanning of 31 channels (1.6 km/s)
    hann_window = np.hanning(31)
    hann_kernel = CustomKernel(hann_window)
    smooth_opacity= convolve(opacity, hann_kernel, boundary='extend')
    smooth_abs = 1.0-smooth_opacity

    good = np.where(smooth_abs>(3.0*sigma_tau_smooth))

    super_good = np.where(smooth_abs>(5.0*sigma_tau_smooth))
    #print('sigma_tau_smooth=',sigma_tau_smooth)
    #print('super_good',super_good)
    sg=np.ravel(np.array(super_good))
    #print('sg=',sg)
    velo=CRVAL3+(sg+10)*CDELT3                              #revised 20181023
    velo=np.around(velo,decimals=1)
    #print('velo',velo)
    nsuper=len(sg)

    fig = plt.figure()
    plt.plot(smooth_abs[good], em_mean[good], color='lightsteelblue', markerfacecolor='blue',
        markeredgecolor='black',marker='.')

    plt.plot(smooth_abs[super_good], em_mean[super_good],markeredgecolor='deeppink',
        markerfacecolor='black', marker='.', linestyle='')

    ax = fig.add_subplot(111)
    for i in np.arange(nsuper):
        if (velo[i]%1) == 0:
            ax.annotate('%s' % velo[i], xy=(smooth_abs[sg[i]], em_mean[sg[i]]), textcoords='data')

    plt.xlabel(r'1-$e^{(-\tau)}$')
    plt.ylabel(r'$T_B$ (K)')
    plt.grid(True)
    plt.savefig(filename)
    #plt.show()
    plt.close()
    return

def produce_spectra(dir_name, field_name):

    with open(dir_name + 'extracted_spectra/' + field_name + '_spectra_web.html', 'w') as spectra_idx:
        t = Template(
            '<html>\n<head><title>D$field_name Spectra</title></head>\n'
            + '<body>\n<h1>Spectra previews for field $field_name</h1>\n<table>\n')
        spectra_idx.write(t.substitute(field_name=field_name))

        neg_mean = 0
        no_mean = 0
        all_cont_sd = []
        all_opacity = []

        spectra, source_ids, islands, CDELT3, CRVAL3 = extract_spectra(dir_name, field_name)
        t = Template('<tr><td colspan=5><b>Field: ${field_name}</b></td></tr>\n' +
            '<tr><td>Details</td>' +
            '<td>Source Image</td><td>Absorption</td><td>Emission</td><td>Emission-Absorption</td></tr>\n')
        spectra_idx.write(t.substitute(field_name=field_name))

        idx = 0
        for longitude in sorted(spectra.keys()):
            spectrum = spectra.get(longitude)
            src_data = source_ids.get(longitude)
            name_prefix = field_name + '_src_' + src_data['id']
            idx += 1
            mean, cont_sd, cont_sd_smooth = get_mean_continuum(spectrum)
            if mean is None:
                print("WARNING: Skipped spectrum %s with no continuum data" % (name_prefix, mean))
                no_mean += 1
                continue

            if mean < 0:
                print(("WARNING: Skipped spectrum %s with negative " +
                    "mean: %.5f") % (name_prefix, mean))
                neg_mean += 1
                continue

            loc = src_data['pos']
            spectrum_name = src_data['id']
            print('Continuum mean of %s (%s) is %.5f Jy, sd %.5f' % (spectrum_name, name_prefix, mean, cont_sd))
            all_cont_sd.append(cont_sd)
            opacity = get_opacity(spectrum, mean)
            peak_tau = np.nanmax(-np.log(opacity))
            temp_bright = get_temp_bright(spectrum, src_data['beam_area'])

            em_mean, em_std , em_pb_mean= get_emission_spectra(src_data['pos'],
                spectrum.velocity, dir_name, name_prefix,
                src_data['a'], src_data['b'], src_data['pa'], islands)      #change field_name to name_prefix 20190617
            # print opacity
            sigma_tau, sigma_tau_smooth = calc_sigma_tau(opacity, cont_sd, cont_sd_smooth, em_pb_mean)
            img_name = dir_name + 'extracted_spectra/' + name_prefix + "_plot.png"
            plot_spectrum(spectrum.velocity, opacity, img_name,"Spectra for source {}".format(spectrum_name), sigma_tau, sigma_tau_smooth)

            img_name = dir_name + 'extracted_spectra/' + name_prefix + "_flux.png"
            plot_original_spectrum(spectrum.velocity, spectrum.flux, img_name,"Spectra for source {}".format(spectrum_name))

            filename = dir_name + 'extracted_spectra/' + name_prefix + '_opacity.votable.xml'
            latitude = src_data['pos'].dec
            em_img_name = dir_name + 'extracted_spectra/' + name_prefix + "_emission.png"

            print('before plot_emission_spectrum')
            print('type(em_mean)',type(em_mean))
            print('len(em_mean)',len(em_mean))
            print('em_mean[1000:1010]', em_mean[1000:1010])
            print('type(spectrum.velocity)',type(spectrum.velocity))
            print('len(spectrum.velocity)',len(spectrum.velocity))
            print('spectrum.velocity[1000:1010]', spectrum.velocity[1000:1010])
            plot_emission_spectrum(spectrum.velocity, em_mean, em_std, #em_pb_mean,
                dir_name + 'extracted_spectra/' + name_prefix + "_emission.png",
                "Emission around {0}".format(spectrum_name))
            output_spectra(spectrum, opacity, filename, longitude, latitude,
                em_mean, em_std, temp_bright, src_data['beam_area'], sigma_tau, sigma_tau_smooth, mean)
            all_opacity.append(opacity)

            em_abs_name = dir_name + 'extracted_spectra/' + name_prefix + "_abs_em_plot.png"
            plot_em_abs(em_mean, em_std, opacity, sigma_tau_smooth, em_abs_name, CDELT3, CRVAL3)

            src_img_name = dir_name + 'extracted_spectra/' + name_prefix + "_region.png"

            # set up so that it writes the location on my Mac
            # mac_dir_name = '/Users/kjameson/Desktop/research/MC_HI_abs/spectra/SMC/'
            # src_img_name_mac = mac_dir_name + 'extracted_spectra/' + name_prefix + "_region.png"
            # em_img_name_mac = mac_dir_name + 'extracted_spectra/' + name_prefix + "_emission.png"
            # img_name_mac = mac_dir_name + 'extracted_spectra/' + name_prefix + "_plot.png"
            # em_abs_name_mac = mac_dir_name + 'extracted_spectra/' + name_prefix + "_abs_em_plot.png"

            src_img_name_web = './plot_files/' + name_prefix + "_region.png"
            em_img_name_web = './plot_files/' + name_prefix + "_emission.png"
            img_name_web = './plot_files/' + name_prefix + "_plot.png"
            em_abs_name_web = './plot_files/' + name_prefix + "_abs_em_plot.png"

            t = Template('<tr><td>${name}<br/>Peak&nbsp;Flux:&nbsp;${peak_flux}&nbsp;Jy<br/>' +
                             'Mean:&nbsp;${mean}&nbsp;Jy<br/>Peak&nbsp;tau:&nbsp;${peak_tau}<br/>'
                             'Cont&nbsp;SD:&nbsp;${cont_sd}</td><td><a href="${src_img}">' +
                             '<img src="${src_img}" width="500px"></a></td><td><a href="${img}">' +
                             '<img src="${img}" width="500px"></a></td><td><a href="${em_img}">' +
                             '<img src="${em_img}" width="500px"></a></td><td><a href="${em_abs_img}">' +
                             '<img src="${em_abs_img}" width="500px"></a></td></tr>\n')
            #spectra_idx.write(t.substitute(src_img=src_img_name_mac, img=img_name_mac, em_img=em_img_name_mac, peak_flux=src_data['flux'],
            #                                mean=mean, cont_sd=cont_sd, name=spectrum_name, peak_tau=peak_tau, em_abs_img=em_abs_name_mac))
            spectra_idx.write(t.substitute(src_img=src_img_name_web, img=img_name_web, em_img=em_img_name_web, peak_flux=src_data['flux'],
                                           mean=mean, cont_sd=cont_sd, name=spectrum_name, peak_tau=peak_tau, em_abs_img=em_abs_name_web))

        spectra_idx.write('</table></body></html>\n')

        if no_mean > 0:
            print("Skipped %d spectra with no continuum data." % no_mean)

        print("Skipped %d spectra with negative mean continuum." % neg_mean)
        print("Produced %d spectra with continuum sd of %.5f." % (
            len(all_cont_sd), np.mean(all_cont_sd)))
        return all_opacity


# Run the script if it is called from the command line
if __name__ == "__main__":
    exit(main())
