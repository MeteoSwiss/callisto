# -*- coding: utf-8 -*-
"""
This script tracks the Sun in azimuthal mode (azimuth and elevation)
Carefully check parameter in external file config.ini
It is sufficient to execute this script once every minute
Created on Thu Sep 14 20:13:50 2017
Updated: 09.05.2020, 05.05.2023, 13.05.2023, 04.08.2023, 23.08.2023
Updated; 23.12.2024 fixed positions for Tcold and Thot
Updated [Andrea]: 18.03.2025 added the scan of the sun for assessing the pointing offset
@author: Chr. Monstein
"""
#------------------------------------------------------------------------------

import ephem
import serial # was serial
import datetime
import math
import os.path
import numpy as np
import time
import configparser

#------------------------------------------------------------------------------
# Read configuration file
config = configparser.ConfigParser()
config.read('configsun.ini')
  
# Location parameter
MyLocation          = ephem.Observer()
MyLocation.lon      = config.get('Location','longitude')
MyLocation.lat      = config.get('Location','latitude')
MyLocation.elev     = config.getfloat('Location','elevation')
MyLocation.temp     = config.getfloat('Location','temperature')
MyLocation.pressure = config.getfloat('Location','pressure')

#------------------------------------------------------------------------------
# Tracker parameter
MaxRange  = config.getfloat('Tracker','MaxRange')
AziMin    = config.getfloat('Tracker','AziMin')
AziMax    = config.getfloat('Tracker','AziMax')
EleMin    = config.getfloat('Tracker','EleMin')
EleMax    = config.getfloat('Tracker','EleMax')
AziPark   = config.getfloat('Tracker','AziPark')
ElePark   = config.getfloat('Tracker','ElePark')
aziref    = config.getfloat('Tracker','aziref')
eleref    = config.getfloat('Tracker','eleref')
ele4aziref = config.getfloat('Tracker','ele4aziref')
azidir    = config.getfloat('Tracker','azidir')
eledir    = config.getfloat('Tracker','eledir')
planecorr = config.getfloat('Tracker','planecorr')
MyComport = config.get('Tracker','MyComport')

#------------------------------------------------------------------------------
# Calibration parameter
CaliProc = config.getboolean('Calibration','CaliProc')
AziCold  = config.getfloat('Calibration','AziCold')
EleCold  = config.getfloat('Calibration','EleCold')
AziHot   = config.getfloat('Calibration','AziHot')
EleHot   = config.getfloat('Calibration','EleHot')
TimeCali = config.getint('Calibration','TimeCali')

#------------------------------------------------------------------------------
# Scanning parameter
def getlist(option, sep=',', chars=None):
    """Return a list from a ConfigParser option. By default, 
       split on a comma and strip whitespaces."""
    return [ chunk.strip(chars) for chunk in option.split(sep) ]

scanning = config.getboolean('Scanning','scanning')
Dazi = np.array(getlist(config.get('Scanning','Dazi')), dtype=np.float32)
Dele = np.array(getlist(config.get('Scanning','Dele')), dtype=np.float32)

################# ADDED BY ANDREA #############################################
#------------------------------------------------------------------------------
# Scanning parameter at minute 59

scanning59 = config.getboolean('Scanning59','scanning59')
Dazi59 = np.array(getlist(config.get('Scanning59','Dazi59')), dtype=np.float32)
Dele59 = np.array(getlist(config.get('Scanning59','Dele59')), dtype=np.float32)
###############################################################################

#------------------------------------------------------------------------------

TelescopeParked = False # to prevent endless parking during the night

#------------------------------------------------------------------------------
# Park telescope after sun-set or after source is out of view

def Park():
    RotorAzi = azidir*(AziPark - aziref) # conversion 0° ... 360° -> +/- 180°
    RotorEle = eledir*(ElePark - eleref)
    print ("Parking at Azimuth  = %7.2f  Elevation   = %6.2f" % (AziPark,ElePark))
    
    try:
        DiSEqC = serial.Serial(
             port     = MyComport,
             baudrate = 9600,
             bytesize = serial.EIGHTBITS,
             parity   = serial.PARITY_NONE,
             timeout  = 2)
        
        print(DiSEqC)
             
        if (DiSEqC.isOpen()):
            print ("Successfully connected to antenna tracker at: "+DiSEqC.portstr)
            cmd = 'max{:6.2f}\r'.format(MaxRange) 
            DiSEqC.write(cmd.encode()) # set MaxRange
            time.sleep(0.02)

            cmd = 'azi{:8.4f}\r'.format(RotorAzi) 
            DiSEqC.write(cmd.encode()) # set azimuth or hour angle
            time.sleep(0.02)

            cmd = 'ele{:8.4f}\r'.format(RotorEle) 
            DiSEqC.write(cmd.encode()) # set elevation or declination
            time.sleep(0.02)
            
            DiSEqC.close()
        
    except Exception as ee: #IOError:
        print("Problem communication with tracker!")
        print(ee)
        raise
    
#------------------------------------------------------------------------------
# Logging events

def Log(x,y,msg):
    dat = datetime.datetime.now().strftime("%Y-%m-%d") # "%Y-%m-%d %H:%M:%S"
    filename ='DISEQ-' + dat + '-Sun.txt'
    if (os.path.exists(filename)):
        with open(filename, "a") as fp:   # Save x/y-data in file
            hms = datetime.datetime.now().strftime("%H:%M:%S.%f")[:-3] # "%Y-%m-%d %H:%M:%S"
            Tdec = float(hms[0:2]) + float(hms[3:5])/60 + float(hms[6:8])/3600 + float(hms[9:])/3600000
            st = '{:8.5f},'.format(Tdec) + '{:8.3f},'.format(x) + '{:8.3f}, '.format(y) + msg
            fp.write(st+'\n') 
    else:
        with open(filename, "a") as fp:   # Save x/y-data in file
            fp.write('DiSEqC radio telescope tracker\n') # write header information
            fp.write('Author: Christian Monstein, HB9SCT\n') # write header information
            fp.write('Version: 2023-04-27\n') # write header information
            fp.write('\n')
            hms = datetime.datetime.now().strftime("%H:%M:%S") # "%Y-%m-%d %H:%M:%S"
            Tdec = float(hms[0:2]) + float(hms[3:5])/60 + float(hms[6:8])/3600
            fp.write('TIME [UT], Azimuth[deg], Elevation[deg], Message\n') # write header information
            st = '{:8.5f},'.format(Tdec) + '{:8.3f},'.format(aziref) + '{:8.3f}, '.format(eleref) + 'Reference positions.\n'
            fp.write(st + '\n') 

#------------------------------------------------------------------------------
# Ephemerides evaluiation

dt = datetime.datetime.now() # PC must run on UT or GPS-time
MyLocation.date = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:02d}'.format(
   dt.year,dt.month,dt.day,dt.hour,dt.minute,dt.second)
#MyLocation.date = '2023/03/14 07:24:00' # example for testing at any date&time
print ('Selected date-time: ',MyLocation.date,' UT')

sun = ephem.Sun()
sun.compute(MyLocation)
azi = math.degrees(sun.az)
ele = math.degrees(sun.alt)
print("Sun data: Azimuth   = %7.2f  Elevation   = %6.2f" % (azi,ele))

lst = MyLocation.sidereal_time()
ha = (lst - sun.ra)/math.pi*180.0
dec = math.degrees(sun.dec)
print("Sun data: Hourangle = %7.2f  Declination = %6.2f" % (ha,dec))

CalibrationActive = False
if (CaliProc): # find out which position for calibration
    if ((dt.minute%TimeCali) == 0): # calibration cold sky up or down
        azi = AziCold
        ele = EleCold
        Log(azi,ele,'Cold sky reference position (Tref = Tcold)')
        CalibrationActive = True
        
    if (((dt.minute-1)%TimeCali) == 0): # calibration hot ground at low elevation
        azi = AziHot
        ele = EleHot
        Log(azi,ele,'Hot reference position e.g. on ground (To=Thot)')
        CalibrationActive = True
        
    if (((dt.minute-2)%TimeCali) == 0): # Calibration finished, back to the sun
        Log(azi,ele,'Back to the Sun') # keep original position of sun
        CalibrationActive = True

#------------------------------------------------------------------------------
# Correction of non-planarity

def ElevationRevision(myazi,myele):
    # constant value depends on local conditions and may change
    # constant can be found by measureing elevation at all azimuth by linear regression
    return (myele - planecorr*myazi)

#------------------------------------------------------------------------------
# track-/scan-/calibration-block

# I think it is more complicated than what I've tried here...
horizontal_squint = (aziref-180)*np.cos(ele4aziref*np.pi/180)
RotorAzi = azidir*(azi - (180 + horizontal_squint/np.cos(ele)))
# RotorAzi = azidir*(azi - aziref) # conversion 0° ... 360° -> +/- 180°
RotorEle = eledir*(ele - eleref) # conversion 0° ... 90° -> +/- 45°
RotorEle = ElevationRevision(RotorAzi,RotorEle) # +/-1.25°

if ((azi > AziMax) or (azi < AziMin) or 
    (ele < EleMin) or (ele > EleMax) or
    (np.abs(RotorEle) > MaxRange)    or
    (np.abs(RotorAzi) > MaxRange)):
    if (TelescopeParked == False):
        Park() # no more tracking activity useful
        print("Rotor position or sun out of range")
        Log(RotorAzi,RotorEle,"Rotor position possibly out of range")
        Log(azi,ele,"Sun position possibly out of range")
        TelescopeParked = True
else:
    print("Rotor and sun position in range: Azimuth   = %7.2f  Elevation   = %6.2f" % (RotorAzi,RotorEle))
    
    try:
        DiSEqC = serial.Serial(
             port     = MyComport,
             baudrate = 9600,
             bytesize = serial.EIGHTBITS,
             parity   = serial.PARITY_NONE,
             timeout  = 2)
        
        print(DiSEqC)
             
        if (DiSEqC.isOpen()):
            print ("Successfully connected to antenna tracker at: "+DiSEqC.portstr)
            cmd = 'max{:6.2f}\r'.format(MaxRange) 
            DiSEqC.write(cmd.encode()) # set MaxRange
            time.sleep(0.02)
            cmd = 'ele{:8.3f}\r'.format(RotorEle) 
            DiSEqC.write(cmd.encode()) # set elevation or declination
            time.sleep(0.02)

            cmd = 'azi{:8.3f}\r'.format(RotorAzi) 
            DiSEqC.write(cmd.encode()) # set azimuth or hour angle
            time.sleep(0.02)

            Log(azi,ele,'New tracking position')   

            if (CalibrationActive == False):
                if (scanning == True):
                    print ('Start scanning')
                    for i in range(0,len(Dazi)):
                        RotorAzi = azidir*(azi - aziref + Dazi[i]) # conversion 0° ... 360° -> +/- 180°
                        cmd = 'azi{:8.3f}\r'.format(RotorAzi) 
                        DiSEqC.write(cmd.encode()) # set azimuth or hour angle
                        time.sleep(1.1)
                        
                        RotorEle = eledir*(ele - eleref + Dele[i])
                        RotorEle = ElevationRevision(RotorAzi,RotorEle) # +/-1.25°
                        cmd = 'ele{:8.3f}\r'.format(RotorEle) 
                        DiSEqC.write(cmd.encode()) # set elevation or declination
                        time.sleep(2.2)
                        print("New scanning position Azimuth   = %7.2f  Elevation   = %6.2f" % (RotorAzi,RotorEle))                    
                        Log(azi + Dazi[i],ele + Dele[i],'New scanning position')   
            
            ########## ADDED BY ANDREA #########################################
            pos_azi = []
            pos_ele = []
            for i in range(len(Dazi59)):
                if i % 2 == 0:
                    pos_azi = np.concatenate((pos_azi, Dazi59))
                else:
                    pos_azi = np.concatenate((pos_azi, Dazi59[::-1]))
                pos_ele = np.concatenate((pos_ele, np.ones(len(Dele59)) * Dele59[i]))

            if ((scanning == False) and 
               (scanning59 == True) and 
               (dt.minute == 59) and
               ((dt.hour == 8) or (dt.hour == 9) or (dt.hour == 10) or (dt.hour == 11) or (dt.hour == 12) or (dt.hour == 13))):
                for i in range(0,len(pos_azi)):
                    RotorAzi = azidir*(azi - aziref + pos_azi[i]) # conversion 0° ... 360° -> +/- 180°
                    cmd = 'azi{:8.3f}\r'.format(RotorAzi)
                    DiSEqC.write(cmd.encode()) # set azimuth or hour angle
                    if i == 0:
                        time.sleep(0.4)
                    else:
                        time.sleep(0.03)

                    RotorEle = eledir*(ele - eleref + pos_ele[i])
                    RotorEle = ElevationRevision(RotorAzi,RotorEle) # +/-1.25°
                    cmd = 'ele{:8.3f}\r'.format(RotorEle)
                    DiSEqC.write(cmd.encode()) # set elevation or declination
                    time.sleep(2.2)
                    print("New scanning position Azimuth   = %7.2f  Elevation   = %6.2f" % (RotorAzi,RotorEle))                    
                    Log(azi + pos_azi[i],ele + pos_ele[i],'New scanning59 position')
                
                RotorAzi = azidir*(azi - aziref) # conversion 0° ... 360° -> +/- 180°
                cmd = 'azi{:8.3f}\r'.format(RotorAzi)
                DiSEqC.write(cmd.encode()) # set azimuth or hour angle
                time.sleep(0.3)

                RotorEle = eledir*(ele - eleref)
                RotorEle = ElevationRevision(RotorAzi,RotorEle) # +/-1.25°
                cmd = 'ele{:8.3f}\r'.format(RotorEle)
                print("End of the scan. Back to the Sun")
                Log(azi,ele,'End of the scan59. Back to the Sun')
            ####################################################################
            
            DiSEqC.close()
        
    except IOError:
        print ("Problem communication with elevation tracker!")
        raise

#------------------------------------------------------------------------------
