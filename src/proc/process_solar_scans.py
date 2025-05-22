'''
This script process the solar scans automatically. It also performs a 2D paraboloid
fit in order to get the location of the Sun in the image and the hpbw.

Author: Andrea Francesco Battaglia

History:
  - 2025/03/27 [Andrea F. Battaglia]: Created.
  - 2025/04/02 [Andrea F. Battaglia]: Script adapted for the 6 solar scans (one every hour).
'''

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
#import time
import configparser
import traceback  # for tracing errors
import argparse
import datetime
import glob

import warnings
warnings.filterwarnings('ignore')

from astropy.io import fits
from scipy.optimize import curve_fit


#**************************************************************************

def log_execution(message, plot_separation=True):
    """
    Description:
        This is to store in an external file a log, whether the script run
        successfully or an error occurred.
    """
    with open('C:\\xrt\\output\\solar_scan\\log_process_solar_scan.txt', 'a') as log_file:
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_file.write(f'\n{timestamp} - {message}')
        if plot_separation == True: 
            log_file.write(f'\n************************************************\n')

#**************************************************************************

def readfit(path):
    
    with fits.open(path) as hdu:
        # https://docs.astropy.org/en/stable/io/fits/
        dict_fitfile = {
            'data'       : hdu[0].data, # .astype(np.uint8)
            'freqs'      : hdu[1].data  ['Frequency'][0], # extract frequency axis
            'timeax'     : hdu[1].data  ['Time'][0],      # extract time axis
            'dT'         : hdu[0].header['CDELT1'],       # extract time resolution
            'date'       : hdu[0].header['DATE-OBS'],     # take first file
            'T0'         : hdu[0].header['TIME-OBS'],     # take first file
            'instrument' : hdu[0].header['INSTRUME'],
            'content'    : hdu[0].header['CONTENT']}
        
    return dict_fitfile

#**************************************************************************

def convert_to_hhmmss(time_df):
    '''
    This function converts the times of the log file of the motors in hh:mm:ss
    '''
    hhmmss = []
    for t in time_df:
        
        # Extract hours
        hh = int(t)
        # Extract minutes
        mm = int((t - hh) * 60)
        # Extract seconds
        ss = int(((t - hh) * 60 - mm) * 60)
        # Extract milliseconds
        mil = int(((t - hh) * 60 - mm - ss/60) * 60000)
        # Format time string
        time_str = f"{str(hh).zfill(2)}:{str(mm).zfill(2)}:{str(ss).zfill(2)}.{str(mil).zfill(3)}"
        hhmmss.append(time_str)
    
    return hhmmss

#**************************************************************************

def gaussian_2d(coords, amplitude, x0, y0, sigma_x, sigma_y):
    '''
    This function is used to fit a 2D Gaussian to the solar scan image.
    This function is now (2025-03-31) obsolete. We now use the 2D paraboloid.
    '''
    x, y = coords
    return amplitude * np.exp(-(((x-x0)**2)/(2*sigma_x**2) + ((y-y0)**2)/(2*sigma_y**2)))

#**************************************************************************

def paraboloid_2d(coords, amplitude, x0, y0, sigma_x, sigma_y):
    '''
    This function is used to fit a 2D paraboloid to the solar scan image.
    '''
    x, y = coords
    return amplitude - 4*np.log(2) * ((x-x0)**2/sigma_x**2 + (y-y0)**2/sigma_y**2)

#**************************************************************************

def main():
  
    ### Path to the folder where to get the FIT files
    path_fit_folder = 'C:\\xrt\\output\\data\\raw\\FITfiles'

    ### Path to the folder where to get the log of the motors
    path_log_motor = 'C:\\xrt\\src\\PythonScripts\\TrackingSun'

    ### Path where to store the diagnostics images
    folder_images = 'C:\\xrt\\output\\solar_scan\\maps_solar-scan'

    ### Directlry where to store the pointing offsets
    directory_csv_offsets = 'C:\\xrt\\output\\solar_scan'

    ### Path to the configsun file, where we get the aziref and eleref values
    path_configsun = 'C:\\xrt\\src\\PythonScripts\\TrackingSun\\configsun.ini'
    
    ### Frequency to consider for the analysis
    freq = 10637 # MHz

    ### If the distance of the center of the fitted paraboloid is larger than this value,
    # we consider the measure as "off-pointing"
    dist_quality_check = 0.7 # deg
    
    ### scan times
    ## "before" and "after" refers to the times we use to estimate the daily flux 3 times a day
    ref_time_scan_morning_before = ['08:59:00', '09:00:00']
    ref_time_scan_morning_after = ['09:59:00', '10:00:00']
    ref_time_scan_noon_before = ['10:59:00', '11:00:00']
    ref_time_scan_noon_after = ['11:59:00', '12:00:00']
    ref_time_scan_afternoon_before = ['12:59:00', '13:00:00']
    ref_time_scan_afternoon_after = ['13:59:00', '14:00:00']

    ###########################################################################
    
    ### The user can define a day in the past for re-analysis
    # By default, today is taken
    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter, 
    epilog='')
    
    parser.add_argument(
        '-day', type=str, metavar=str.__name__, 
        action='store', required=False, default=None,
        help='Choose day to analyse (default: today), format: YYYYMMDD')
    
    args = parser.parse_args()
    day = args.day
    
    ### If a day was not specified, use today's date
    if not day: day = datetime.datetime.now().strftime('%Y%m%d')
    filename_format = 'meteoswiss_%s_*_01.fit' % day #MyFile = 'SWISS-METEO_20231127_*_01.fit'
    filepath_format = os.path.join(path_fit_folder, filename_format)
    list_fit = sorted(glob.glob(filepath_format))
    list_fit.sort()

    ### scan time ranges
    scan_time_morning_before = [pd.to_datetime(day + ' ' + ref_time_scan_morning_before[0]), pd.to_datetime(day + ' ' + ref_time_scan_morning_before[-1])]
    scan_time_morning_after = [pd.to_datetime(day + ' ' + ref_time_scan_morning_after[0]), pd.to_datetime(day + ' ' + ref_time_scan_morning_after[-1])]
    scan_time_noon_before = [pd.to_datetime(day + ' ' + ref_time_scan_noon_before[0]), pd.to_datetime(day + ' ' + ref_time_scan_noon_before[-1])]
    scan_time_noon_after = [pd.to_datetime(day + ' ' + ref_time_scan_noon_after[0]), pd.to_datetime(day + ' ' + ref_time_scan_noon_after[-1])]
    scan_time_afternoon_before = [pd.to_datetime(day + ' ' + ref_time_scan_afternoon_before[0]), pd.to_datetime(day + ' ' + ref_time_scan_afternoon_before[-1])]
    scan_time_afternoon_after = [pd.to_datetime(day + ' ' + ref_time_scan_afternoon_after[0]), pd.to_datetime(day + ' ' + ref_time_scan_afternoon_after[-1])]

    ### Read configuration file and get aziref and eleref
    config = configparser.ConfigParser()
    config.read(path_configsun)
    aziref    = config.getfloat('Tracker','aziref')
    eleref    = config.getfloat('Tracker','eleref')
    
    try:
        
        ### read txt file
        date_motor = day[0:4]+'-'+day[4:6]+'-'+day[6:]
        filename_format_log = 'DISEQ-%s-Sun.txt' % date_motor
        path_log_motor = os.path.join(path_log_motor, filename_format_log)
        df_motor = pd.read_csv(path_log_motor, sep=',', header=3)
        time_head = df_motor.columns[0]
        azi_head = df_motor.columns[1]
        ele_head = df_motor.columns[2]
        msg_head = df_motor.columns[3]
        time_df = df_motor[time_head]

        ### Read the values from the log file
        df_motor_morning_before = df_motor[(df_motor[time_head] > 8.98) & (df_motor[time_head] < 9.0)]
        df_motor_morning_after = df_motor[(df_motor[time_head] > 9.98) & (df_motor[time_head] < 10.0)]
        df_motor_noon_before = df_motor[(df_motor[time_head] > 10.98) & (df_motor[time_head] < 11.0)]
        df_motor_noon_after = df_motor[(df_motor[time_head] > 11.98) & (df_motor[time_head] < 12.0)]
        df_motor_afternoon_before = df_motor[(df_motor[time_head] > 12.98) & (df_motor[time_head] < 13.0)]
        df_motor_afternoon_after = df_motor[(df_motor[time_head] > 13.98) & (df_motor[time_head] < 14.0)]
        msg_morning_before = df_motor_morning_before[msg_head].values
        msg_morning_after = df_motor_morning_after[msg_head].values
        msg_noon_before = df_motor_noon_before[msg_head].values
        msg_noon_after = df_motor_noon_after[msg_head].values
        msg_afternoon_before = df_motor_afternoon_before[msg_head].values
        msg_afternoon_after = df_motor_afternoon_after[msg_head].values
        time_df_morning_before = df_motor_morning_before[time_head].values
        time_df_morning_after = df_motor_morning_after[time_head].values
        time_df_noon_before = df_motor_noon_before[time_head].values
        time_df_noon_after = df_motor_noon_after[time_head].values
        time_df_afternoon_before = df_motor_afternoon_before[time_head].values
        time_df_afternoon_after = df_motor_afternoon_after[time_head].values
        azi_morning_before = df_motor_morning_before[azi_head].values
        azi_morning_after = df_motor_morning_after[azi_head].values
        azi_noon_before = df_motor_noon_before[azi_head].values
        azi_noon_after = df_motor_noon_after[azi_head].values
        azi_afternoon_before = df_motor_afternoon_before[azi_head].values
        azi_afternoon_after = df_motor_afternoon_after[azi_head].values
        ele_morning_before = df_motor_morning_before[ele_head].values
        ele_morning_after = df_motor_morning_after[ele_head].values
        ele_noon_before = df_motor_noon_before[ele_head].values
        ele_noon_after = df_motor_noon_after[ele_head].values
        ele_afternoon_before = df_motor_afternoon_before[ele_head].values
        ele_afternoon_after = df_motor_afternoon_after[ele_head].values

        ### Convert the times from hours to hhmmss
        hhmmss = convert_to_hhmmss(time_df)
        hhmmss_morning_before = convert_to_hhmmss(time_df_morning_before)
        hhmmss_morning_after = convert_to_hhmmss(time_df_morning_after)
        hhmmss_noon_before = convert_to_hhmmss(time_df_noon_before)
        hhmmss_noon_after = convert_to_hhmmss(time_df_noon_after)
        hhmmss_afternoon_before = convert_to_hhmmss(time_df_afternoon_before)
        hhmmss_afternoon_after = convert_to_hhmmss(time_df_afternoon_after)

        ### Create the time axes of the scan
        #filename_motor = path_log_motor.split('/')[-1]
        #yyyy_motor = filename_motor.split('-')[1]
        #mm_motor = filename_motor.split('-')[2]
        #dd_motor = filename_motor.split('-')[3]
        #date_motor = yyyy_motor + '-' + mm_motor + '-' + dd_motor
        time_axis_scan = [pd.to_datetime(date_motor + 'T' + hhmmss[ii]) for ii in range(len(time_df))]
        time_axis_scan_morning_before = [pd.to_datetime(date_motor + 'T' + hhmmss_morning_before[ii]) for ii in range(len(time_df_morning_before))]
        time_axis_scan_morning_after = [pd.to_datetime(date_motor + 'T' + hhmmss_morning_after[ii]) for ii in range(len(time_df_morning_after))]
        time_axis_scan_noon_before = [pd.to_datetime(date_motor + 'T' + hhmmss_noon_before[ii]) for ii in range(len(time_df_noon_before))]
        time_axis_scan_noon_after = [pd.to_datetime(date_motor + 'T' + hhmmss_noon_after[ii]) for ii in range(len(time_df_noon_after))]
        time_axis_scan_afternoon_before = [pd.to_datetime(date_motor + 'T' + hhmmss_afternoon_before[ii]) for ii in range(len(time_df_afternoon_before))]
        time_axis_scan_afternoon_after = [pd.to_datetime(date_motor + 'T' + hhmmss_afternoon_after[ii]) for ii in range(len(time_df_afternoon_after))]

        #### Find indices scan positions
        idx_start_scanning_morning_before = np.where(msg_morning_before == ' New scanning59 position')[0][0]
        idx_end_scanning_morning_before = np.where(msg_morning_before == ' New scanning59 position')[0][-1]
        idx_start_scanning_morning_after = np.where(msg_morning_after == ' New scanning59 position')[0][0]
        idx_end_scanning_morning_after = np.where(msg_morning_after == ' New scanning59 position')[0][-1]
        idx_start_scanning_noon_before = np.where(msg_noon_before == ' New scanning59 position')[0][0]
        idx_end_scanning_noon_before = np.where(msg_noon_before == ' New scanning59 position')[0][-1]
        idx_start_scanning_noon_after = np.where(msg_noon_after == ' New scanning59 position')[0][0]
        idx_end_scanning_noon_after = np.where(msg_noon_after == ' New scanning59 position')[0][-1]
        idx_start_scanning_afternoon_before = np.where(msg_afternoon_before == ' New scanning59 position')[0][0]
        idx_end_scanning_afternoon_before = np.where(msg_afternoon_before == ' New scanning59 position')[0][-1]
        idx_start_scanning_afternoon_after = np.where(msg_afternoon_after == ' New scanning59 position')[0][0]
        idx_end_scanning_afternoon_after = np.where(msg_afternoon_after == ' New scanning59 position')[0][-1]

        #### Just select the time of interests for the scan
        time_scan_morning_before = time_axis_scan_morning_before[idx_start_scanning_morning_before:idx_end_scanning_morning_before+1]
        azi_scan_morning_before = azi_morning_before[idx_start_scanning_morning_before:idx_end_scanning_morning_before+1]
        ele_scan_morning_before = ele_morning_before[idx_start_scanning_morning_before:idx_end_scanning_morning_before+1]
        time_scan_morning_after = time_axis_scan_morning_after[idx_start_scanning_morning_after:idx_end_scanning_morning_after+1]
        azi_scan_morning_after = azi_morning_after[idx_start_scanning_morning_after:idx_end_scanning_morning_after+1]
        ele_scan_morning_after = ele_morning_after[idx_start_scanning_morning_after:idx_end_scanning_morning_after+1]
        #msg_scan_morning = msg_morning[idx_start_scanning_morning:idx_end_scanning_morning+1]

        time_scan_noon_before = time_axis_scan_noon_before[idx_start_scanning_noon_before:idx_end_scanning_noon_before+1]
        azi_scan_noon_before = azi_noon_before[idx_start_scanning_noon_before:idx_end_scanning_noon_before+1]
        ele_scan_noon_before = ele_noon_before[idx_start_scanning_noon_before:idx_end_scanning_noon_before+1]
        time_scan_noon_after = time_axis_scan_noon_after[idx_start_scanning_noon_after:idx_end_scanning_noon_after+1]
        azi_scan_noon_after = azi_noon_after[idx_start_scanning_noon_after:idx_end_scanning_noon_after+1]
        ele_scan_noon_after = ele_noon_after[idx_start_scanning_noon_after:idx_end_scanning_noon_after+1]
        #msg_scan_noon = msg_noon[idx_start_scanning_noon:idx_end_scanning_noon+1]

        time_scan_afternoon_before = time_axis_scan_afternoon_before[idx_start_scanning_afternoon_before:idx_end_scanning_afternoon_before+1]
        azi_scan_afternoon_before = azi_afternoon_before[idx_start_scanning_afternoon_before:idx_end_scanning_afternoon_before+1]
        ele_scan_afternoon_before = ele_afternoon_before[idx_start_scanning_afternoon_before:idx_end_scanning_afternoon_before+1]
        time_scan_afternoon_after = time_axis_scan_afternoon_after[idx_start_scanning_afternoon_after:idx_end_scanning_afternoon_after+1]
        azi_scan_afternoon_after = azi_afternoon_after[idx_start_scanning_afternoon_after:idx_end_scanning_afternoon_after+1]
        ele_scan_afternoon_after = ele_afternoon_after[idx_start_scanning_afternoon_after:idx_end_scanning_afternoon_after+1]
        #msg_scan_afternoon = msg_afternoon[idx_start_scanning_afternoon:idx_end_scanning_afternoon+1]

        ### find in listfit the file containing 084500, 104500, 124500
        idx_file_morning_before = [ii for ii in range(len(list_fit)) if '0845' in list_fit[ii]][0]
        idx_file_morning_after = [ii for ii in range(len(list_fit)) if '0945' in list_fit[ii]][0]
        idx_file_noon_before = [ii for ii in range(len(list_fit)) if '1045' in list_fit[ii]][0]
        idx_file_noon_after = [ii for ii in range(len(list_fit)) if '1145' in list_fit[ii]][0]
        idx_file_afternoon_before = [ii for ii in range(len(list_fit)) if '1245' in list_fit[ii]][0]
        idx_file_afternoon_after = [ii for ii in range(len(list_fit)) if '1345' in list_fit[ii]][0]

        ### open the first FIT for some info
        dict_morning_before = readfit(list_fit[idx_file_morning_before])
        data_morning_before = dict_morning_before['data']
        timeax_morning_before = dict_morning_before['timeax']
        T0_morning_before = dict_morning_before['T0']
        freqs_morning_before = dict_morning_before['freqs']
        date_morning_before = dict_morning_before['date']
        time_axis_FIT_morning_before = pd.to_datetime(date_morning_before + ' ' + T0_morning_before) + pd.to_timedelta(timeax_morning_before, unit='s')

        dict_morning_after = readfit(list_fit[idx_file_morning_after])
        data_morning_after = dict_morning_after['data']
        timeax_morning_after = dict_morning_after['timeax']
        T0_morning_after = dict_morning_after['T0']
        freqs_morning_after = dict_morning_after['freqs']
        date_morning_after = dict_morning_after['date']
        time_axis_FIT_morning_after = pd.to_datetime(date_morning_after + ' ' + T0_morning_after) + pd.to_timedelta(timeax_morning_after, unit='s')

        dict_noon_before = readfit(list_fit[idx_file_noon_before])
        data_noon_before = dict_noon_before['data']
        timeax_noon_before = dict_noon_before['timeax']
        T0_noon_before = dict_noon_before['T0']
        freqs_noon_before = dict_noon_before['freqs']
        date_noon_before = dict_noon_before['date']
        time_axis_FIT_noon_before = pd.to_datetime(date_noon_before + ' ' + T0_noon_before) + pd.to_timedelta(timeax_noon_before, unit='s')

        dict_noon_after = readfit(list_fit[idx_file_noon_after])
        data_noon_after = dict_noon_after['data']
        timeax_noon_after = dict_noon_after['timeax']
        T0_noon_after = dict_noon_after['T0']
        freqs_noon_after = dict_noon_after['freqs']
        date_noon_after = dict_noon_after['date']
        time_axis_FIT_noon_after = pd.to_datetime(date_noon_after + ' ' + T0_noon_after) + pd.to_timedelta(timeax_noon_after, unit='s')

        dict_afternoon_before = readfit(list_fit[idx_file_afternoon_before])
        data_afternoon_before = dict_afternoon_before['data']
        timeax_afternoon_before = dict_afternoon_before['timeax']
        T0_afternoon_before = dict_afternoon_before['T0']
        freqs_afternoon_before = dict_afternoon_before['freqs']
        date_afternoon_before = dict_afternoon_before['date']
        time_axis_FIT_afternoon_before = pd.to_datetime(date_afternoon_before + ' ' + T0_afternoon_before) + pd.to_timedelta(timeax_afternoon_before, unit='s')

        dict_afternoon_after = readfit(list_fit[idx_file_afternoon_after])
        data_afternoon_after = dict_afternoon_after['data']
        timeax_afternoon_after = dict_afternoon_after['timeax']
        T0_afternoon_after = dict_afternoon_after['T0']
        freqs_afternoon_after = dict_afternoon_after['freqs']
        date_afternoon_after = dict_afternoon_after['date']
        time_axis_FIT_afternoon_after = pd.to_datetime(date_afternoon_after + ' ' + T0_afternoon_after) + pd.to_timedelta(timeax_afternoon_after, unit='s')

        ### Get the index of the frequency to look at
        idx_freq = np.argmin(np.abs(freqs_morning_before - freq))
        freq = freqs_morning_before[idx_freq]

        ##### only select times relevant for the scan in the FIT files
        morning_mask_before = (time_axis_FIT_morning_before >= scan_time_morning_before[0]) & (time_axis_FIT_morning_before <= scan_time_morning_before[1])
        morning_mask_after = (time_axis_FIT_morning_after >= scan_time_morning_after[0]) & (time_axis_FIT_morning_after <= scan_time_morning_after[1])
        noon_mask_before = (time_axis_FIT_noon_before >= scan_time_noon_before[0]) & (time_axis_FIT_noon_before <= scan_time_noon_before[1])
        noon_mask_after = (time_axis_FIT_noon_after >= scan_time_noon_after[0]) & (time_axis_FIT_noon_after <= scan_time_noon_after[1])
        afternoon_mask_before = (time_axis_FIT_afternoon_before >= scan_time_afternoon_before[0]) & (time_axis_FIT_afternoon_before <= scan_time_afternoon_before[1])
        afternoon_mask_after = (time_axis_FIT_afternoon_after >= scan_time_afternoon_after[0]) & (time_axis_FIT_afternoon_after <= scan_time_afternoon_after[1])
        data_FIT_morning_before = data_morning_before[idx_freq,morning_mask_before]
        data_FIT_morning_after = data_morning_after[idx_freq,morning_mask_after]
        data_FIT_noon_before = data_noon_before[idx_freq,noon_mask_before]
        data_FIT_noon_after = data_noon_after[idx_freq,noon_mask_after]
        data_FIT_afternoon_before = data_afternoon_before[idx_freq,afternoon_mask_before]
        data_FIT_afternoon_after = data_afternoon_after[idx_freq,afternoon_mask_after]
        time_axis_FIT_morning_before = time_axis_FIT_morning_before[morning_mask_before]
        time_axis_FIT_morning_after = time_axis_FIT_morning_after[morning_mask_after]
        time_axis_FIT_noon_before = time_axis_FIT_noon_before[noon_mask_before]
        time_axis_FIT_noon_after = time_axis_FIT_noon_after[noon_mask_after]
        time_axis_FIT_afternoon_before = time_axis_FIT_afternoon_before[afternoon_mask_before]
        time_axis_FIT_afternoon_after = time_axis_FIT_afternoon_after[afternoon_mask_after]

        ### we add +pd.Timedelta(1.1, unit='s') to the motor time to get the time at the middle of the scan
        time_scan_morning_before = [time_scan_morning_before[ii] + pd.Timedelta(1.1, unit='s') for ii in range(len(time_scan_morning_before))]
        time_scan_morning_after = [time_scan_morning_after[ii] + pd.Timedelta(1.1, unit='s') for ii in range(len(time_scan_morning_after))]
        time_scan_noon_before = [time_scan_noon_before[ii] + pd.Timedelta(1.1, unit='s') for ii in range(len(time_scan_noon_before))]
        time_scan_noon_after = [time_scan_noon_after[ii] + pd.Timedelta(1.1, unit='s') for ii in range(len(time_scan_noon_after))]
        time_scan_afternoon_before = [time_scan_afternoon_before[ii] + pd.Timedelta(1.1, unit='s') for ii in range(len(time_scan_afternoon_before))]
        time_scan_afternoon_after = [time_scan_afternoon_after[ii] + pd.Timedelta(1.1, unit='s') for ii in range(len(time_scan_afternoon_after))]
        
        ### average the data over 1.5 s around each center position of time_array_Dtime_morning, time_array_Dtime_noon, time_array_Dtime_afternoon
        data_array_Dtime_morning_before = []
        data_array_Dtime_morning_after = []
        data_array_Dtime_noon_before = []
        data_array_Dtime_noon_after = []
        data_array_Dtime_afternoon_before = []
        data_array_Dtime_afternoon_after = []
        for nn in range(len(azi_scan_noon_before)):
            idx_morning_before = np.argmin(np.abs(pd.to_datetime(time_axis_FIT_morning_before) - time_scan_morning_before[nn]))
            idx_noon_before = np.argmin(np.abs(pd.to_datetime(time_axis_FIT_noon_before) - time_scan_noon_before[nn]))
            idx_afternoon_before = np.argmin(np.abs(pd.to_datetime(time_axis_FIT_afternoon_before) - time_scan_afternoon_before[nn]))
            idx_morning_after = np.argmin(np.abs(pd.to_datetime(time_axis_FIT_morning_after) - time_scan_morning_after[nn]))
            idx_noon_after = np.argmin(np.abs(pd.to_datetime(time_axis_FIT_noon_after) - time_scan_noon_after[nn]))
            idx_afternoon_after = np.argmin(np.abs(pd.to_datetime(time_axis_FIT_afternoon_after) - time_scan_afternoon_after[nn]))
            data_array_Dtime_morning_before.append(np.mean(data_FIT_morning_before[idx_morning_before-3:idx_morning_before+3]))
            data_array_Dtime_noon_before.append(np.mean(data_FIT_noon_before[idx_noon_before-3:idx_noon_before+3]))
            data_array_Dtime_afternoon_before.append(np.mean(data_FIT_afternoon_before[idx_afternoon_before-3:idx_afternoon_before+3]))
            data_array_Dtime_morning_after.append(np.mean(data_FIT_morning_after[idx_morning_after-3:idx_morning_after+3]))
            data_array_Dtime_noon_after.append(np.mean(data_FIT_noon_after[idx_noon_after-3:idx_noon_after+3]))
            data_array_Dtime_afternoon_after.append(np.mean(data_FIT_afternoon_after[idx_afternoon_after-3:idx_afternoon_after+3]))

        ### Define the axes of the maps, being the delta of the azimuth and elevation
        delta_azi_scan_morning_before = azi_scan_morning_before - np.mean(azi_scan_morning_before)
        delta_ele_scan_morning_before = ele_scan_morning_before - np.mean(ele_scan_morning_before)
        delta_azi_scan_noon_before = azi_scan_noon_before - np.mean(azi_scan_noon_before)
        delta_ele_scan_noon_before = ele_scan_noon_before - np.mean(ele_scan_noon_before)
        delta_azi_scan_afternoon_before = azi_scan_afternoon_before - np.mean(azi_scan_afternoon_before)
        delta_ele_scan_afternoon_before = ele_scan_afternoon_before - np.mean(ele_scan_afternoon_before)

        delta_azi_scan_morning_after = azi_scan_morning_after - np.mean(azi_scan_morning_after)
        delta_ele_scan_morning_after = ele_scan_morning_after - np.mean(ele_scan_morning_after)
        delta_azi_scan_noon_after = azi_scan_noon_after - np.mean(azi_scan_noon_after)
        delta_ele_scan_noon_after = ele_scan_noon_after - np.mean(ele_scan_noon_after)
        delta_azi_scan_afternoon_after = azi_scan_afternoon_after - np.mean(azi_scan_afternoon_after)
        delta_ele_scan_afternoon_after = ele_scan_afternoon_after - np.mean(ele_scan_afternoon_after)

        ### Here we correct the delta azimuth for the solid angle effects due to the elevation
        delta_azi_scan_morning_before = [delta_azi_scan_morning_before[ii] * np.cos(ele_scan_morning_before[ii]*np.pi/180) for ii in range(len(delta_azi_scan_morning_before))]
        delta_azi_scan_noon_before = [delta_azi_scan_noon_before[ii] * np.cos(ele_scan_noon_before[ii]*np.pi/180) for ii in range(len(delta_azi_scan_noon_before))]
        delta_azi_scan_afternoon_before = [delta_azi_scan_afternoon_before[ii] * np.cos(ele_scan_afternoon_before[ii]*np.pi/180) for ii in range(len(delta_azi_scan_afternoon_before))]
        delta_azi_scan_morning_after = [delta_azi_scan_morning_after[ii] * np.cos(ele_scan_morning_after[ii]*np.pi/180) for ii in range(len(delta_azi_scan_morning_after))]
        delta_azi_scan_noon_after = [delta_azi_scan_noon_after[ii] * np.cos(ele_scan_noon_after[ii]*np.pi/180) for ii in range(len(delta_azi_scan_noon_after))]
        delta_azi_scan_afternoon_after = [delta_azi_scan_afternoon_after[ii] * np.cos(ele_scan_afternoon_after[ii]*np.pi/180) for ii in range(len(delta_azi_scan_afternoon_after))]
        
        ### This is just to define the ranges of the plot, plus some other stuff
        vmax = np.max([np.abs(delta_azi_scan_noon_after), np.abs(delta_ele_scan_noon_after)])
        vmax += 0.2
        vmin = -vmax
        idx_max_morning_before = np.argmax(data_array_Dtime_morning_before)
        idx_max_noon_before = np.argmax(data_array_Dtime_noon_before)
        idx_max_afternoon_before = np.argmax(data_array_Dtime_afternoon_before)
        idx_max_morning_after = np.argmax(data_array_Dtime_morning_after)
        idx_max_noon_after = np.argmax(data_array_Dtime_noon_after)
        idx_max_afternoon_after = np.argmax(data_array_Dtime_afternoon_after)

        ### Do the fit
        ### Create meshgrid for smooth plotting
        x = np.linspace(-vmax, vmax, 100)
        y = np.linspace(-vmax, vmax, 100)
        X, Y = np.meshgrid(x, y)

        plot_fit_morning_before = True
        quality_check_morning_before = 1
        try:
            #params_morning, _ = curve_fit(gaussian_2d, (delta_azi_scan_morning, delta_ele_scan_morning), data_array_Dtime_morning, 
            #                            p0=[np.max(data_array_Dtime_morning), 0, 0, 1.5, 1.5])
            #Z_morning = gaussian_2d((X, Y), *params_morning)
            params_morning_before, _ = curve_fit(paraboloid_2d, (delta_azi_scan_morning_before, delta_ele_scan_morning_before), data_array_Dtime_morning_before, 
                                        p0=[np.max(data_array_Dtime_morning_before), 0, 0, 1.5, 1.5],
                                        bounds=([0.01, -3, -3, 0.1, 0.1],
                                                [2*np.max(data_array_Dtime_morning_before), 3, 3, 6, 6]))
            Z_morning_before = paraboloid_2d((X, Y), *params_morning_before)
            Z_max_morning_before = np.max(Z_morning_before)
            labels_morning_before = {
                Z_max_morning_before*0.9999: 'max',
                Z_max_morning_before-3: 'hpbw'
            }
            dist_morning_before = np.sqrt(params_morning_before[1]**2+params_morning_before[2]**2)
            if dist_morning_before > dist_quality_check: quality_check_morning_before = 0
        except:
            plot_fit_morning_before = False
            params_morning_before = [np.nan, np.nan, np.nan, np.nan, np.nan]
            quality_check_morning_before = 0
            log_execution('The fit of the scan at 08:59 has failed!')
        
        plot_fit_morning_after = True
        quality_check_morning_after = 1
        try:
            #params_morning, _ = curve_fit(gaussian_2d, (delta_azi_scan_morning, delta_ele_scan_morning), data_array_Dtime_morning, 
            #                            p0=[np.max(data_array_Dtime_morning), 0, 0, 1.5, 1.5])
            #Z_morning = gaussian_2d((X, Y), *params_morning)
            params_morning_after, _ = curve_fit(paraboloid_2d, (delta_azi_scan_morning_after, delta_ele_scan_morning_after), data_array_Dtime_morning_after, 
                                        p0=[np.max(data_array_Dtime_morning_after), 0, 0, 1.5, 1.5],
                                        bounds=([0.01, -3, -3, 0.1, 0.1],
                                                [2*np.max(data_array_Dtime_morning_after), 3, 3, 6, 6]))
            Z_morning_after = paraboloid_2d((X, Y), *params_morning_after)
            Z_max_morning_after = np.max(Z_morning_after)
            labels_morning_after = {
                Z_max_morning_after*0.9999: 'max',
                Z_max_morning_after-3: 'hpbw'
            }
            dist_morning_after = np.sqrt(params_morning_after[1]**2+params_morning_after[2]**2)
            if dist_morning_after > dist_quality_check: quality_check_morning_after = 0
        except:
            plot_fit_morning_after = False
            params_morning_after = [np.nan, np.nan, np.nan, np.nan, np.nan]
            quality_check_morning_after = 0
            log_execution('The fit of the scan at 09:59 has failed!')

        plot_fit_noon_before = True
        quality_check_noon_before = 1
        try:   
            #params_noon, _ = curve_fit(gaussian_2d, (delta_azi_scan_noon, delta_ele_scan_noon), data_array_Dtime_noon,
            #                        p0=[np.max(data_array_Dtime_noon), 0, 0, 1.5, 1.5])
            #Z_noon = gaussian_2d((X, Y), *params_noon)
            params_noon_before, _ = curve_fit(paraboloid_2d, (delta_azi_scan_noon_before, delta_ele_scan_noon_before), data_array_Dtime_noon_before,
                                    p0=[np.max(data_array_Dtime_noon_before), 0, 0, 1.5, 1.5],
                                        bounds=([0.01, -3, -3, 0.1, 0.1],
                                                [2*np.max(data_array_Dtime_noon_before), 3, 3, 6, 6]))
            Z_noon_before = paraboloid_2d((X, Y), *params_noon_before)
            Z_max_noon_before = np.max(Z_noon_before)
            labels_noon_before = {
                Z_max_noon_before*0.9999: 'max',
                Z_max_noon_before-3: 'hpbw'
            }
            dist_noon_before = np.sqrt(params_noon_before[1]**2+params_noon_before[2]**2)
            if dist_noon_before > dist_quality_check: quality_check_noon_before = 0
        except:
            plot_fit_noon_before = False
            params_noon_before = [np.nan, np.nan, np.nan, np.nan, np.nan]
            quality_check_noon_before = 0
            log_execution('The fit of the scan at 10:59 has failed!')
        
        plot_fit_noon_after = True
        quality_check_noon_after = 1
        try:   
            #params_noon, _ = curve_fit(gaussian_2d, (delta_azi_scan_noon, delta_ele_scan_noon), data_array_Dtime_noon,
            #                        p0=[np.max(data_array_Dtime_noon), 0, 0, 1.5, 1.5])
            #Z_noon = gaussian_2d((X, Y), *params_noon)
            params_noon_after, _ = curve_fit(paraboloid_2d, (delta_azi_scan_noon_after, delta_ele_scan_noon_after), data_array_Dtime_noon_after,
                                    p0=[np.max(data_array_Dtime_noon_after), 0, 0, 1.5, 1.5],
                                        bounds=([0.01, -3, -3, 0.1, 0.1],
                                                [2*np.max(data_array_Dtime_noon_after), 3, 3, 6, 6]))
            Z_noon_after = paraboloid_2d((X, Y), *params_noon_after)
            Z_max_noon_after = np.max(Z_noon_after)
            labels_noon_after = {
                Z_max_noon_after*0.9999: 'max',
                Z_max_noon_after-3: 'hpbw'
            }
            dist_noon_after = np.sqrt(params_noon_after[1]**2+params_noon_after[2]**2)
            if dist_noon_after > dist_quality_check: quality_check_noon_after = 0
        except:
            plot_fit_noon_after = False
            params_noon_after = [np.nan, np.nan, np.nan, np.nan, np.nan]
            quality_check_noon_after = 0
            log_execution('The fit of the scan at 11:59 has failed!')
            
        plot_fit_afternoon_before = True
        quality_check_afternoon_before = 1
        try:    
            #params_afternoon, _ = curve_fit(gaussian_2d, (delta_azi_scan_morning, delta_ele_scan_afternoon), data_array_Dtime_afternoon,
            #                            p0=[np.max(data_array_Dtime_afternoon), 0, 0, 1.5, 1.5])
            #Z_afternoon = gaussian_2d((X, Y), *params_afternoon)
            params_afternoon_before, _ = curve_fit(paraboloid_2d, (delta_azi_scan_morning_before, delta_ele_scan_afternoon_before), data_array_Dtime_afternoon_before,
                                        p0=[np.max(data_array_Dtime_afternoon_before), 0, 0, 1.5, 1.5],
                                        bounds=([0.01, -3, -3, 0.1, 0.1],
                                                [2*np.max(data_array_Dtime_afternoon_before), 3, 3, 6, 6]))
            Z_afternoon_before = paraboloid_2d((X, Y), *params_afternoon_before)
            Z_max_afternoon_before = np.max(Z_afternoon_before)
            labels_afternoon_before = {
                Z_max_afternoon_before*0.9999: 'max',
                Z_max_afternoon_before-3: 'hpbw'
            }
            dist_afternoon_before = np.sqrt(params_afternoon_before[1]**2+params_afternoon_before[2]**2)
            if dist_afternoon_before > dist_quality_check: quality_check_afternoon_before = 0
        except:
            plot_fit_afternoon_before = False
            params_afternoon_before = [np.nan, np.nan, np.nan, np.nan, np.nan]
            quality_check_afternoon_before = 0
            log_execution('The fit of the scan at 12:59 has failed!')
            
        plot_fit_afternoon_after = True
        quality_check_afternoon_after = 1
        try:    
            #params_afternoon, _ = curve_fit(gaussian_2d, (delta_azi_scan_morning, delta_ele_scan_afternoon), data_array_Dtime_afternoon,
            #                            p0=[np.max(data_array_Dtime_afternoon), 0, 0, 1.5, 1.5])
            #Z_afternoon = gaussian_2d((X, Y), *params_afternoon)
            params_afternoon_after, _ = curve_fit(paraboloid_2d, (delta_azi_scan_morning_after, delta_ele_scan_afternoon_after), data_array_Dtime_afternoon_after,
                                        p0=[np.max(data_array_Dtime_afternoon_after), 0, 0, 1.5, 1.5],
                                        bounds=([0.01, -3, -3, 0.1, 0.1],
                                                [2*np.max(data_array_Dtime_afternoon_after), 3, 3, 6, 6]))
            Z_afternoon_after = paraboloid_2d((X, Y), *params_afternoon_after)
            Z_max_afternoon_after = np.max(Z_afternoon_after)
            labels_afternoon_after = {
                Z_max_afternoon_after*0.9999: 'max',
                Z_max_afternoon_after-3: 'hpbw'
            }
            dist_afternoon_after = np.sqrt(params_afternoon_after[1]**2+params_afternoon_after[2]**2)
            if dist_afternoon_after > dist_quality_check: quality_check_afternoon_after = 0
        except:
            plot_fit_afternoon_after = False
            params_afternoon_after = [np.nan, np.nan, np.nan, np.nan, np.nan]
            quality_check_afternoon_after = 0
            log_execution('The fit of the scan at 13:59 has failed!')


        ### Function to format the labels in the contour plot
        def custom_fmt_morning_before(x):
            return labels_morning_before.get(x, '')
        def custom_fmt_morning_after(x):
            return labels_morning_after.get(x, '')

        def custom_fmt_noon_before(x):
            return labels_noon_before.get(x, '')
        def custom_fmt_noon_after(x):
            return labels_noon_after.get(x, '')

        def custom_fmt_afternoon_before(x):
            return labels_afternoon_before.get(x, '')
        def custom_fmt_afternoon_after(x):
            return labels_afternoon_after.get(x, '')

        
        
        ###############################################################################
        ######### First diagnostics plot: time profiles with extracted values #########
        ###############################################################################

        dt_scan = 0.75
        
        plt.rcParams.update({'font.size': 18})
        fig = plt.figure(figsize=(25,14))

        plt.suptitle('Solar scan at '+str(int(freq))+' MHz on '+date_morning_before+' (aziref = '+str(aziref)+' deg; eleref = '+str(eleref)+' deg)', y=0.98)

        ax_morning_before = fig.add_subplot(321)
        ax_morning_before.plot(time_axis_FIT_morning_before, data_FIT_morning_before, label='Morning', lw=2)
        ax_morning_before.scatter(time_scan_morning_before, data_array_Dtime_morning_before, c='r', s=100, label='Morning (averaged)')
        ax_morning_before.set_ylabel('dBadu')
        ax_morning_before.set_xlabel('UT')
        ax_morning_before.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M:%S'))
        ax_morning_before.grid()
        ax_morning_before.tick_params(axis='x', direction='in', top=True, bottom=True, width=1., length=7)
        ax_morning_before.tick_params(axis='y', direction='in', left=True, right=True, width=1., length=7)
        ax_morning_before.legend()
        for i in range(len(time_scan_morning_before)):
            ax_morning_before.vlines(time_scan_morning_before[i], 0, 255, color='k', linestyle='dotted', alpha=0.5)
            ax_morning_before.axvspan(time_scan_morning_before[i]-pd.to_timedelta(dt_scan/2, unit='s'), time_scan_morning_before[i]+pd.to_timedelta(dt_scan/2, unit='s'), color='gray', alpha=0.4)
        ax_morning_before.set_ylim([np.min(data_FIT_morning_before)-3, np.max(data_FIT_morning_before)+3])
        ax_morning_before.xaxis.set_minor_locator(plt.matplotlib.dates.SecondLocator())
        ax_morning_before.tick_params(axis='x', which='minor', direction='in', top=True, bottom=True, width=1, length=4)

        ax_morning_after = fig.add_subplot(322)
        ax_morning_after.plot(time_axis_FIT_morning_after, data_FIT_morning_after, label='Morning', lw=2)
        ax_morning_after.scatter(time_scan_morning_after, data_array_Dtime_morning_after, c='r', s=100, label='Morning (averaged)')
        ax_morning_after.set_ylabel('dBadu')
        ax_morning_after.set_xlabel('UT')
        ax_morning_after.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M:%S'))
        ax_morning_after.grid()
        ax_morning_after.tick_params(axis='x', direction='in', top=True, bottom=True, width=1., length=7)
        ax_morning_after.tick_params(axis='y', direction='in', left=True, right=True, width=1., length=7)
        ax_morning_after.legend()
        for i in range(len(time_scan_morning_after)):
            ax_morning_after.vlines(time_scan_morning_after[i], 0, 255, color='k', linestyle='dotted', alpha=0.5)
            ax_morning_after.axvspan(time_scan_morning_after[i]-pd.to_timedelta(dt_scan/2, unit='s'), time_scan_morning_after[i]+pd.to_timedelta(dt_scan/2, unit='s'), color='gray', alpha=0.4)
        ax_morning_after.set_ylim([np.min(data_FIT_morning_after)-3, np.max(data_FIT_morning_after)+3])
        ax_morning_after.xaxis.set_minor_locator(plt.matplotlib.dates.SecondLocator())
        ax_morning_after.tick_params(axis='x', which='minor', direction='in', top=True, bottom=True, width=1, length=4)
        
        ax_noon_before = fig.add_subplot(323)
        ax_noon_before.plot(time_axis_FIT_noon_before, data_FIT_noon_before, label='Noon', lw=2)
        ax_noon_before.scatter(time_scan_noon_before, data_array_Dtime_noon_before, c='r', s=100, label='Morning (averaged)')
        ax_noon_before.set_ylabel('dBadu')
        ax_noon_before.set_xlabel('UT')
        ax_noon_before.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M:%S'))
        ax_noon_before.grid()
        ax_noon_before.tick_params(axis='x', direction='in', top=True, bottom=True, width=1., length=7)
        ax_noon_before.tick_params(axis='y', direction='in', left=True, right=True, width=1., length=7)
        ax_noon_before.legend()
        for i in range(len(time_scan_noon_before)):
            ax_noon_before.vlines(time_scan_noon_before[i], 0, 255, color='k', linestyle='dotted', alpha=0.5)
            ax_noon_before.axvspan(time_scan_noon_before[i]-pd.to_timedelta(dt_scan/2, unit='s'), time_scan_noon_before[i]+pd.to_timedelta(dt_scan/2, unit='s'), color='gray', alpha=0.4)  
        ax_noon_before.set_ylim([np.min(data_FIT_noon_before)-3, np.max(data_FIT_noon_before)+3])
        ax_noon_before.xaxis.set_minor_locator(plt.matplotlib.dates.SecondLocator())
        ax_noon_before.tick_params(axis='x', which='minor', direction='in', top=True, bottom=True, width=1, length=4)

        ax_noon_after = fig.add_subplot(324)
        ax_noon_after.plot(time_axis_FIT_noon_after, data_FIT_noon_after, label='Noon', lw=2)
        ax_noon_after.scatter(time_scan_noon_after, data_array_Dtime_noon_after, c='r', s=100, label='Morning (averaged)')
        ax_noon_after.set_ylabel('dBadu')
        ax_noon_after.set_xlabel('UT')
        ax_noon_after.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M:%S'))
        ax_noon_after.grid()
        ax_noon_after.tick_params(axis='x', direction='in', top=True, bottom=True, width=1., length=7)
        ax_noon_after.tick_params(axis='y', direction='in', left=True, right=True, width=1., length=7)
        ax_noon_after.legend()
        for i in range(len(time_scan_noon_after)):
            ax_noon_after.vlines(time_scan_noon_after[i], 0, 255, color='k', linestyle='dotted', alpha=0.5)
            ax_noon_after.axvspan(time_scan_noon_after[i]-pd.to_timedelta(dt_scan/2, unit='s'), time_scan_noon_after[i]+pd.to_timedelta(dt_scan/2, unit='s'), color='gray', alpha=0.4)  
        ax_noon_after.set_ylim([np.min(data_FIT_noon_after)-3, np.max(data_FIT_noon_after)+3])
        ax_noon_after.xaxis.set_minor_locator(plt.matplotlib.dates.SecondLocator())
        ax_noon_after.tick_params(axis='x', which='minor', direction='in', top=True, bottom=True, width=1, length=4)

        ax_afternoon_before = fig.add_subplot(325)
        ax_afternoon_before.plot(time_axis_FIT_afternoon_before, data_FIT_afternoon_before, label='Afternoon', lw=2)
        ax_afternoon_before.scatter(time_scan_afternoon_before, data_array_Dtime_afternoon_before, c='r', s=100, label='Afternoon (averaged)')
        ax_afternoon_before.set_ylabel('dBadu')
        ax_afternoon_before.set_xlabel('UT')
        ax_afternoon_before.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M:%S'))
        ax_afternoon_before.grid()
        ax_afternoon_before.tick_params(axis='x', direction='in', top=True, bottom=True, width=1., length=7)
        ax_afternoon_before.tick_params(axis='y', direction='in', left=True, right=True, width=1., length=7)
        ax_afternoon_before.legend()
        for i in range(len(time_scan_afternoon_before)):
            ax_afternoon_before.vlines(time_scan_afternoon_before[i], 0, 255, color='k', linestyle='dotted', alpha=0.5)
            ax_afternoon_before.axvspan(time_scan_afternoon_before[i]-pd.to_timedelta(dt_scan/2, unit='s'), time_scan_afternoon_before[i]+pd.to_timedelta(dt_scan/2, unit='s'), color='gray', alpha=0.4)
        ax_afternoon_before.set_ylim([np.min(data_FIT_afternoon_before)-3, np.max(data_FIT_afternoon_before)+3])
        ax_afternoon_before.xaxis.set_minor_locator(plt.matplotlib.dates.SecondLocator())
        ax_afternoon_before.tick_params(axis='x', which='minor', direction='in', top=True, bottom=True, width=1, length=4)

        ax_afternoon_after = fig.add_subplot(326)
        ax_afternoon_after.plot(time_axis_FIT_afternoon_after, data_FIT_afternoon_after, label='Afternoon', lw=2)
        ax_afternoon_after.scatter(time_scan_afternoon_after, data_array_Dtime_afternoon_after, c='r', s=100, label='Afternoon (averaged)')
        ax_afternoon_after.set_ylabel('dBadu')
        ax_afternoon_after.set_xlabel('UT')
        ax_afternoon_after.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M:%S'))
        ax_afternoon_after.grid()
        ax_afternoon_after.tick_params(axis='x', direction='in', top=True, bottom=True, width=1., length=7)
        ax_afternoon_after.tick_params(axis='y', direction='in', left=True, right=True, width=1., length=7)
        ax_afternoon_after.legend()
        for i in range(len(time_scan_afternoon_after)):
            ax_afternoon_after.vlines(time_scan_afternoon_after[i], 0, 255, color='k', linestyle='dotted', alpha=0.5)
            ax_afternoon_after.axvspan(time_scan_afternoon_after[i]-pd.to_timedelta(dt_scan/2, unit='s'), time_scan_afternoon_after[i]+pd.to_timedelta(dt_scan/2, unit='s'), color='gray', alpha=0.4)
        ax_afternoon_after.set_ylim([np.min(data_FIT_afternoon_after)-3, np.max(data_FIT_afternoon_after)+3])
        ax_afternoon_after.xaxis.set_minor_locator(plt.matplotlib.dates.SecondLocator())
        ax_afternoon_after.tick_params(axis='x', which='minor', direction='in', top=True, bottom=True, width=1, length=4)

        plt.tight_layout()
        plt.savefig(os.path.join(folder_images,'solar-scan59_time-profiles_'+date_motor+'.png'))
        plt.close()

        ###############################################################################

        ###############################################################################
        ################## Second diagnostics plot: map of the scans ##################
        ###############################################################################

        fig = plt.figure(figsize=(26,16))

        fig.suptitle('Solar scan: '+date_morning_before+'; aziref = '+str(aziref)+' deg; eleref = '+str(eleref)+' deg', y=0.98, fontsize=22, bbox=dict(facecolor='gray', alpha=0.2))

        ax_m_before = plt.subplot(231)
        ax_m_before.plot(delta_azi_scan_morning_before, delta_ele_scan_morning_before, 'k--', alpha=0.4, lw=2)
        if plot_fit_morning_before:
            contours_m_before = ax_m_before.contour(X, Y, Z_morning_before, 
                                    levels=[Z_max_morning_before-3, Z_max_morning_before*0.9999], 
                                    colors=['k', 'k'], alpha=0.5, linewidths=2)
            ax_m_before.clabel(contours_m_before, inline=True, fontsize=16, fmt=custom_fmt_morning_before)
            ax_m_before.text(params_morning_before[1]+0.1, params_morning_before[2]+0.1, '({:.2f}, {:.2f})'.format(params_morning_before[1], params_morning_before[2]), fontsize=18, color='gray')
        cf_m_before = ax_m_before.scatter(delta_azi_scan_morning_before, delta_ele_scan_morning_before, c=data_array_Dtime_morning_before, s=100, cmap='hot', 
                        edgecolor='k', linewidth=1)
        ax_m_before.set_xlabel(r'$\delta_{azimuth}$ [deg]')
        ax_m_before.set_ylabel(r'$\delta_{elevation}$ [deg]')
        ax_m_before.set_title('Morning: '+ref_time_scan_morning_before[0] + ' - ' + ref_time_scan_morning_before[1] + ' UT')
        ax_m_before.set_xlim([-vmax, vmax])
        ax_m_before.set_ylim([-vmax, vmax])
        ax_m_before.set_xticks(np.arange(-vmax, vmax+0.1, 0.4))
        ax_m_before.set_yticks(np.arange(-vmax, vmax+0.1, 0.4))
        ax_m_before.tick_params(axis='x', direction='in', top=True, bottom=True, width=1., length=7)
        ax_m_before.tick_params(axis='y', direction='in', left=True, right=True, width=1., length=7)
        ax_m_before.grid()
        cbar_m_before = plt.colorbar(cf_m_before)
        cbar_m_before.set_label('dBadu', labelpad=-59)
        ax_m_before.add_patch(plt.Circle((delta_azi_scan_morning_before[idx_max_morning_before], delta_ele_scan_morning_before[idx_max_morning_before]), 0.12, color='r', fill=False, lw=2))
        ax_m_before.text(delta_azi_scan_morning_before[idx_max_morning_before], delta_ele_scan_morning_before[idx_max_morning_before]+0.15, 'Max', color='r', fontsize=14, ha='center')
        if quality_check_morning_before == 0: ax_m_before.text(0, 0.3, 'off-pointing', fontsize=40, color='red', fontweight='bold', alpha=0.5, ha='center', va='center')
        
        ax_m_after = plt.subplot(234)
        ax_m_after.plot(delta_azi_scan_morning_after, delta_ele_scan_morning_after, 'k--', alpha=0.4, lw=2)
        if plot_fit_morning_after:
            contours_m_after = ax_m_after.contour(X, Y, Z_morning_after, 
                                    levels=[Z_max_morning_after-3, Z_max_morning_after*0.9999], 
                                    colors=['k', 'k'], alpha=0.5, linewidths=2)
            ax_m_after.clabel(contours_m_after, inline=True, fontsize=16, fmt=custom_fmt_morning_after)
            ax_m_after.text(params_morning_after[1]+0.1, params_morning_after[2]+0.1, '({:.2f}, {:.2f})'.format(params_morning_after[1], params_morning_after[2]), fontsize=18, color='gray')
        cf_m_after = ax_m_after.scatter(delta_azi_scan_morning_after, delta_ele_scan_morning_after, c=data_array_Dtime_morning_after, s=100, cmap='hot', 
                        edgecolor='k', linewidth=1)
        ax_m_after.set_xlabel(r'$\delta_{azimuth}$ [deg]')
        ax_m_after.set_ylabel(r'$\delta_{elevation}$ [deg]')
        ax_m_after.set_title('Morning: '+ref_time_scan_morning_after[0] + ' - ' + ref_time_scan_morning_after[1] + ' UT')
        ax_m_after.set_xlim([-vmax, vmax])
        ax_m_after.set_ylim([-vmax, vmax])
        ax_m_after.set_xticks(np.arange(-vmax, vmax+0.1, 0.4))
        ax_m_after.set_yticks(np.arange(-vmax, vmax+0.1, 0.4))
        ax_m_after.tick_params(axis='x', direction='in', top=True, bottom=True, width=1., length=7)
        ax_m_after.tick_params(axis='y', direction='in', left=True, right=True, width=1., length=7)
        ax_m_after.grid()
        cbar_m_after = plt.colorbar(cf_m_after)
        cbar_m_after.set_label('dBadu', labelpad=-59)
        ax_m_after.add_patch(plt.Circle((delta_azi_scan_morning_after[idx_max_morning_after], delta_ele_scan_morning_after[idx_max_morning_after]), 0.12, color='r', fill=False, lw=2))
        ax_m_after.text(delta_azi_scan_morning_after[idx_max_morning_after], delta_ele_scan_morning_after[idx_max_morning_after]+0.15, 'Max', color='r', fontsize=14, ha='center')
        if quality_check_morning_after == 0: ax_m_after.text(0, 0.3, 'off-pointing', fontsize=40, color='red', fontweight='bold', alpha=0.5, ha='center', va='center')
        
        ax_n_before = plt.subplot(232)
        ax_n_before.plot(delta_azi_scan_noon_before, delta_ele_scan_morning_before, 'k--', alpha=0.4)
        if plot_fit_noon_before:
            contours_n_before = ax_n_before.contour(X, Y, Z_noon_before, 
                                    levels=[Z_max_noon_before-3, Z_max_noon_before*0.9999], 
                                    colors=['k', 'k', 'k'], alpha=0.5, linewidths=2)
            ax_n_before.clabel(contours_n_before, inline=True, fontsize=16, fmt=labels_noon_before)
            ax_n_before.text(params_noon_before[1]+0.1, params_noon_before[2]+0.1, '({:.2f}, {:.2f})'.format(params_noon_before[1], params_noon_before[2]), fontsize=18, color='gray')
        cf_n_before = ax_n_before.scatter(delta_azi_scan_noon_before, delta_ele_scan_noon_before, c=data_array_Dtime_noon_before, s=100, cmap='hot', 
                        edgecolor='k', linewidth=1)
        ax_n_before.set_xlabel(r'$\delta_{azimuth}$ [deg]')
        ax_n_before.set_ylabel(r'$\delta_{elevation}$ [deg]')
        ax_n_before.set_title('Noon: ' + ref_time_scan_noon_before[0] + ' - ' + ref_time_scan_noon_before[1] + ' UT')
        ax_n_before.set_xlim([-vmax, vmax])
        ax_n_before.set_ylim([-vmax, vmax])
        ax_n_before.set_xticks(np.arange(-vmax, vmax+0.1, 0.4))
        ax_n_before.set_yticks(np.arange(-vmax, vmax+0.1, 0.4))
        ax_n_before.tick_params(axis='x', direction='in', top=True, bottom=True, width=1., length=7)
        ax_n_before.tick_params(axis='y', direction='in', left=True, right=True, width=1., length=7)
        ax_n_before.grid()
        cbar_n_before = plt.colorbar(cf_n_before)
        cbar_n_before.set_label('dBadu', labelpad=-59)
        ax_n_before.add_patch(plt.Circle((delta_azi_scan_noon_before[idx_max_noon_before], delta_ele_scan_noon_before[idx_max_noon_before]), 0.12, color='r', fill=False, lw=2))
        ax_n_before.text(delta_azi_scan_noon_before[idx_max_noon_before], delta_ele_scan_noon_before[idx_max_noon_before]+0.15, 'Max', color='r', fontsize=14, ha='center')
        if quality_check_noon_before == 0: ax_n_before.text(0, 0.3, 'off-pointing', fontsize=40, color='red', fontweight='bold', alpha=0.5, ha='center', va='center')
        
        ax_n_after = plt.subplot(235)
        ax_n_after.plot(delta_azi_scan_noon_after, delta_ele_scan_morning_after, 'k--', alpha=0.4)
        if plot_fit_noon_after:
            contours_n_after = ax_n_after.contour(X, Y, Z_noon_after, 
                                    levels=[Z_max_noon_after-3, Z_max_noon_after*0.9999], 
                                    colors=['k', 'k', 'k'], alpha=0.5, linewidths=2)
            ax_n_after.clabel(contours_n_after, inline=True, fontsize=16, fmt=labels_noon_after)
            ax_n_after.text(params_noon_after[1]+0.1, params_noon_after[2]+0.1, '({:.2f}, {:.2f})'.format(params_noon_after[1], params_noon_after[2]), fontsize=18, color='gray')
        cf_n_after = ax_n_after.scatter(delta_azi_scan_noon_after, delta_ele_scan_noon_after, c=data_array_Dtime_noon_after, s=100, cmap='hot', 
                        edgecolor='k', linewidth=1)
        ax_n_after.set_xlabel(r'$\delta_{azimuth}$ [deg]')
        ax_n_after.set_ylabel(r'$\delta_{elevation}$ [deg]')
        ax_n_after.set_title('Noon: ' + ref_time_scan_noon_after[0] + ' - ' + ref_time_scan_noon_after[1] + ' UT')
        ax_n_after.set_xlim([-vmax, vmax])
        ax_n_after.set_ylim([-vmax, vmax])
        ax_n_after.set_xticks(np.arange(-vmax, vmax+0.1, 0.4))
        ax_n_after.set_yticks(np.arange(-vmax, vmax+0.1, 0.4))
        ax_n_after.tick_params(axis='x', direction='in', top=True, bottom=True, width=1., length=7)
        ax_n_after.tick_params(axis='y', direction='in', left=True, right=True, width=1., length=7)
        ax_n_after.grid()
        cbar_n_after = plt.colorbar(cf_n_after)
        cbar_n_after.set_label('dBadu', labelpad=-59)
        ax_n_after.add_patch(plt.Circle((delta_azi_scan_noon_after[idx_max_noon_after], delta_ele_scan_noon_after[idx_max_noon_after]), 0.12, color='r', fill=False, lw=2))
        ax_n_after.text(delta_azi_scan_noon_after[idx_max_noon_after], delta_ele_scan_noon_after[idx_max_noon_after]+0.15, 'Max', color='r', fontsize=14, ha='center')
        if quality_check_noon_after == 0: ax_n_after.text(0, 0.3, 'off-pointing', fontsize=40, color='red', fontweight='bold', alpha=0.5, ha='center', va='center')
        
        ax_a_before = plt.subplot(233)
        ax_a_before.plot(delta_azi_scan_afternoon_before, delta_ele_scan_afternoon_before, 'k--', alpha=0.4)
        if plot_fit_afternoon_before:
            contours_a_before = ax_a_before.contour(X, Y, Z_afternoon_before, 
                                    levels=[Z_max_afternoon_before-3, Z_max_afternoon_before*0.9999], 
                                    colors=['k', 'k', 'k'], alpha=0.5, linewidths=2)
            ax_a_before.clabel(contours_a_before, inline=True, fontsize=15, fmt=labels_afternoon_before)
            ax_a_before.text(params_afternoon_before[1]+0.1, params_afternoon_before[2]+0.1, '({:.2f}, {:.2f})'.format(params_afternoon_before[1], params_afternoon_before[2]), fontsize=18, color='gray')
        cf_a_before = ax_a_before.scatter(delta_azi_scan_afternoon_before, delta_ele_scan_afternoon_before, c=data_array_Dtime_afternoon_before, s=100, cmap='hot', 
                        edgecolor='k', linewidth=1)
        ax_a_before.set_xlabel(r'$\delta_{azimuth}$ [deg]')
        ax_a_before.set_ylabel(r'$\delta_{elevation}$ [deg]')
        ax_a_before.set_title('Afternoon: ' + ref_time_scan_afternoon_before[0] + ' - ' + ref_time_scan_afternoon_before[1] + ' UT')
        ax_a_before.set_xlim([-vmax, vmax])
        ax_a_before.set_ylim([-vmax, vmax])
        ax_a_before.set_xticks(np.arange(-vmax, vmax+0.1, 0.4))
        ax_a_before.set_yticks(np.arange(-vmax, vmax+0.1, 0.4))
        ax_a_before.tick_params(axis='x', direction='in', top=True, bottom=True, width=1., length=7)
        ax_a_before.tick_params(axis='y', direction='in', left=True, right=True, width=1., length=7)
        ax_a_before.grid()
        cbar_a_before = plt.colorbar(cf_a_before)
        cbar_a_before.set_label('dBadu', labelpad=-59)
        ax_a_before.add_patch(plt.Circle((delta_azi_scan_afternoon_before[idx_max_afternoon_before], delta_ele_scan_afternoon_before[idx_max_afternoon_before]), 0.12, color='r', fill=False, lw=2))
        ax_a_before.text(delta_azi_scan_afternoon_before[idx_max_afternoon_before], delta_ele_scan_afternoon_before[idx_max_afternoon_before]+0.15, 'Max', color='r', fontsize=14, ha='center')
        if quality_check_afternoon_before == 0: ax_a_before.text(0, 0.3, 'off-pointing', fontsize=40, color='red', fontweight='bold', alpha=0.5, ha='center', va='center')
        
        ax_a_after = plt.subplot(236)
        ax_a_after.plot(delta_azi_scan_afternoon_after, delta_ele_scan_afternoon_after, 'k--', alpha=0.4)
        if plot_fit_afternoon_after:
            contours_a_after = ax_a_after.contour(X, Y, Z_afternoon_after, 
                                    levels=[Z_max_afternoon_after-3, Z_max_afternoon_after*0.9999], 
                                    colors=['k', 'k', 'k'], alpha=0.5, linewidths=2)
            ax_a_after.clabel(contours_a_after, inline=True, fontsize=15, fmt=labels_afternoon_after)
            ax_a_after.text(params_afternoon_after[1]+0.1, params_afternoon_after[2]+0.1, '({:.2f}, {:.2f})'.format(params_afternoon_after[1], params_afternoon_after[2]), fontsize=18, color='gray')
        cf_a_after = ax_a_after.scatter(delta_azi_scan_afternoon_after, delta_ele_scan_afternoon_after, c=data_array_Dtime_afternoon_after, s=100, cmap='hot', 
                        edgecolor='k', linewidth=1)
        ax_a_after.set_xlabel(r'$\delta_{azimuth}$ [deg]')
        ax_a_after.set_ylabel(r'$\delta_{elevation}$ [deg]')
        ax_a_after.set_title('Afternoon: ' + ref_time_scan_afternoon_after[0] + ' - ' + ref_time_scan_afternoon_after[1] + ' UT')
        ax_a_after.set_xlim([-vmax, vmax])
        ax_a_after.set_ylim([-vmax, vmax])
        ax_a_after.set_xticks(np.arange(-vmax, vmax+0.1, 0.4))
        ax_a_after.set_yticks(np.arange(-vmax, vmax+0.1, 0.4))
        ax_a_after.tick_params(axis='x', direction='in', top=True, bottom=True, width=1., length=7)
        ax_a_after.tick_params(axis='y', direction='in', left=True, right=True, width=1., length=7)
        ax_a_after.grid()
        cbar_a_after = plt.colorbar(cf_a_after)
        cbar_a_after.set_label('dBadu', labelpad=-59)
        ax_a_after.add_patch(plt.Circle((delta_azi_scan_afternoon_after[idx_max_afternoon_after], delta_ele_scan_afternoon_after[idx_max_afternoon_after]), 0.12, color='r', fill=False, lw=2))
        ax_a_after.text(delta_azi_scan_afternoon_after[idx_max_afternoon_after], delta_ele_scan_afternoon_after[idx_max_afternoon_after]+0.15, 'Max', color='r', fontsize=14, ha='center')
        if quality_check_afternoon_after == 0: ax_a_after.text(0, 0.3, 'off-pointing', fontsize=40, color='red', fontweight='bold', alpha=0.5, ha='center', va='center')
        
        plt.tight_layout()
        plt.savefig(os.path.join(folder_images,'solar-scan59_map_'+date_motor+'.png'), dpi=300)
        plt.close()
        
        #############################################################################################
        
        ### If everything runs successfully
        log_execution("Successfully run")
        
    except Exception as e:
        
        ### Save NaNs in the csv file if it fails
        params_morning_before = [np.nan, np.nan, np.nan, np.nan, np.nan]
        params_noon_before = [np.nan, np.nan, np.nan, np.nan, np.nan]
        params_afternoon_before = [np.nan, np.nan, np.nan, np.nan, np.nan]
        params_morning_after = [np.nan, np.nan, np.nan, np.nan, np.nan]
        params_noon_after = [np.nan, np.nan, np.nan, np.nan, np.nan]
        params_afternoon_after = [np.nan, np.nan, np.nan, np.nan, np.nan]
        quality_check_morning_before = 0
        quality_check_morning_after = 0
        quality_check_noon_before = 0
        quality_check_noon_after = 0
        quality_check_afternoon_before = 0
        quality_check_afternoon_after = 0
        
        ### If an error occurs, log the error message and stack trace
        error_message = f"Error: {str(e)}\n{traceback.format_exc()}"
        log_execution(error_message)

    
    ### Find the csv with the pointing offsets
    YYYY = str(scan_time_noon_before[0].date().year)
    filename_csv_offsets = 'pointing-offsets_%s.csv' % YYYY
    filepath_csv_offsets = os.path.join(directory_csv_offsets, filename_csv_offsets)
    # if it does not exist, create it
    if os.path.exists(filepath_csv_offsets):
        df_offsets = pd.read_csv(filepath_csv_offsets)
    else:
        df_offsets = pd.DataFrame(columns=['start_scan_UTC', 'end_scan_UTC', 'amplitude_fit', 'delta_azi_fit', 'delta_ele_fit', 'sigma_azi_fit', 'sigma_ele_fit', 'aziref', 'eleref', 'quality_check'])
        df_offsets.to_csv(filepath_csv_offsets, index=False)
    
    
    ### Store the results in the update the csv file in df_flux by adding the new values
    df_offsets = pd.read_csv(filepath_csv_offsets)
    new_row = pd.DataFrame({'start_scan_UTC': str(scan_time_morning_before[0]), 
                            'end_scan_UTC': str(scan_time_morning_before[1]), 
                            'amplitude_fit': params_morning_before[0], 
                            'delta_azi_fit': params_morning_before[1],
                            'delta_ele_fit': params_morning_before[2], 
                            'sigma_azi_fit': params_morning_before[3], 
                            'sigma_ele_fit': params_morning_before[4],
                            'aziref': aziref,
                            'eleref': eleref,
                            'quality_check': quality_check_morning_before}, index=[0])
    df_offsets = pd.concat([df_offsets,new_row], ignore_index=True)
    new_row = pd.DataFrame({'start_scan_UTC': str(scan_time_morning_after[0]), 
                            'end_scan_UTC': str(scan_time_morning_after[1]), 
                            'amplitude_fit': params_morning_after[0], 
                            'delta_azi_fit': params_morning_after[1],
                            'delta_ele_fit': params_morning_after[2], 
                            'sigma_azi_fit': params_morning_after[3], 
                            'sigma_ele_fit': params_morning_after[4],
                            'aziref': aziref,
                            'eleref': eleref,
                            'quality_check': quality_check_morning_after}, index=[0])
    df_offsets = pd.concat([df_offsets,new_row], ignore_index=True)
    new_row = pd.DataFrame({'start_scan_UTC': str(scan_time_noon_before[0]), 
                            'end_scan_UTC': str(scan_time_noon_before[1]), 
                            'amplitude_fit': params_noon_before[0], 
                            'delta_azi_fit': params_noon_before[1],
                            'delta_ele_fit': params_noon_before[2], 
                            'sigma_azi_fit': params_noon_before[3], 
                            'sigma_ele_fit': params_noon_before[4],
                            'aziref': aziref,
                            'eleref': eleref,
                            'quality_check': quality_check_noon_before}, index=[0])
    df_offsets = pd.concat([df_offsets,new_row], ignore_index=True)
    new_row = pd.DataFrame({'start_scan_UTC': str(scan_time_noon_after[0]), 
                            'end_scan_UTC': str(scan_time_noon_after[1]), 
                            'amplitude_fit': params_noon_after[0], 
                            'delta_azi_fit': params_noon_after[1],
                            'delta_ele_fit': params_noon_after[2], 
                            'sigma_azi_fit': params_noon_after[3], 
                            'sigma_ele_fit': params_noon_after[4],
                            'aziref': aziref,
                            'eleref': eleref,
                            'quality_check': quality_check_noon_after}, index=[0])
    df_offsets = pd.concat([df_offsets,new_row], ignore_index=True)
    new_row = pd.DataFrame({'start_scan_UTC': str(scan_time_afternoon_before[0]), 
                            'end_scan_UTC': str(scan_time_afternoon_before[1]), 
                            'amplitude_fit': params_afternoon_before[0], 
                            'delta_azi_fit': params_afternoon_before[1],
                            'delta_ele_fit': params_afternoon_before[2], 
                            'sigma_azi_fit': params_afternoon_before[3], 
                            'sigma_ele_fit': params_afternoon_before[4],
                            'aziref': aziref,
                            'eleref': eleref,
                            'quality_check': quality_check_afternoon_before}, index=[0])
    df_offsets = pd.concat([df_offsets,new_row], ignore_index=True)
    new_row = pd.DataFrame({'start_scan_UTC': str(scan_time_afternoon_after[0]), 
                            'end_scan_UTC': str(scan_time_afternoon_after[1]), 
                            'amplitude_fit': params_afternoon_after[0], 
                            'delta_azi_fit': params_afternoon_after[1],
                            'delta_ele_fit': params_afternoon_after[2], 
                            'sigma_azi_fit': params_afternoon_after[3], 
                            'sigma_ele_fit': params_afternoon_after[4],
                            'aziref': aziref,
                            'eleref': eleref,
                            'quality_check': quality_check_afternoon_after}, index=[0])
    df_offsets = pd.concat([df_offsets,new_row], ignore_index=True)
    df_offsets = df_offsets.sort_values(by='start_scan_UTC')
    df_offsets.to_csv(filepath_csv_offsets, index=False)
    
    
    
if __name__ == '__main__':
    main()