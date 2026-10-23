import pandas as pd
import numpy as np
import datetime
import glob

from PyBSI.tools import *

def convert_folder(path, nFBG : int=14) :
    files = glob.glob(path)
    files.sort()    
    for i, file in enumerate(files):
        if i == 0:
            absTime, refTime, OutWav = read(file, nFBG=nFBG)
            continue
        try :
            TabsTime, TrefTime, TOutWav = read(file, nFBG=nFBG)
            absTime, OutWav = np.append(absTime, TabsTime), np.append(OutWav, TOutWav, axis=1)
        except :
            print('Error file : %s' %file)
        print(absTime[0], absTime[-1])
        if verif_date(absTime) == True :
            IdxEndTime = np.where(absTime.astype('datetime64[D]')==absTime[-1].astype('datetime64[D]'))[0][0]
            OabsTime, OOutWav = absTime[:IdxEndTime], OutWav[:,:IdxEndTime]
            absTime, OutWav = absTime[IdxEndTime:], OutWav[:,IdxEndTime:]
            export_nc(OabsTime, OOutWav, refwav=None)
    try :
        export_nc(absTime, refTime, OutWav, refwav=None)
    except : 
        print('Error exporting netCDF file')
        pass

def read(filepath, refwav=None, refTime=None,  nFBG : int=14):
    #need to implement the refwav
    print('Reading file : %s' %filepath)
    info = filepath.split('/')
    sttime = info[-1][:16]

    date = datetime.datetime.strptime(sttime, '%Y%m%dT%H%M%SZ')
    timing = np.datetime64(date, 'ms')
    header = ['Time [sec]'] + ['FBG_%02d' %i for i in range(1,nFBG)]
    data = pd.read_csv(filepath, delimiter='\t', decimal=',', names=header, skiprows=1)
    data = data.to_numpy()
    rtime = data[:,0]
    data = data[:,1:].T

    urtime = (rtime*1000).astype(int)
    if refTime is None:
        refTime = timing
    print(refTime)
    absTime = np.array(refTime) + np.array(urtime) - urtime[0]


    return absTime, refTime, data

def export_nc(absTime, OutWav, refwav=None, path='./'):
    if refwav is None :
        refwav = OutWav[:,0]
 
    data = xr.Dataset(data_vars=dict(Wavelength_nm = (["FBG","Time"], OutWav),
                                    Calibration_Wavelength = (['FBG'], refwav)),
                      coords = dict(Time = absTime),
                      attrs=dict(description="FBG data -  Dimension of the data (FBG*Time)"))
    OutName = str(absTime[0]).replace('-','.').replace(':','.')
    data.to_netcdf(path+'/FBG_HSV_%s.nc' %(OutName))