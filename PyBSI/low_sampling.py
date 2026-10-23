import pandas as pd
import numpy as np
import warnings
import os
import csv 
import glob

from PyBSI.tools import *


def read(filepath, refwav=None, nChannel=8, maxFBGs=12, extended=False):
    """
    Read BSI FBG text files with multichannel data and return a 3D array of wavelengths over time.

    This function reads a file containing multichannel FBG (Fiber Bragg Grating) data, processes the 
    wavelengths for each channel, and returns the data in a structured format. It can optionally 
    return the data in a dictionary format with additional timing information.

    Parameters
    ----------
    filepath : str
        Path to the BSI text file containing the FBG data.
    
    refwav : np.ndarray, optional
        A 2D array used as the reference wavelength for the channels and FBGs. If not provided, the first
        entry in the data file will be used as the reference (default is None).
    
    nChannel : int, optional
        Number of channels to read from the file (default is 8).
    
    maxFBGs : int, optional
        Maximum number of FBGs per channel to read (default is 12).
    
    extended : bool, optional
        If True, the function returns dictionaries with channel-wise data and timestamps. If False, it returns
        simple arrays (default is False).

    Returns
    -------
    If extended is False:
        absTime : np.ndarray
            Array of absolute timestamps in `datetime64` format.
        
        refTime : np.ndarray
            Array of reference times in float format.
        
        OutWav : np.ndarray
            3D array of wavelengths with dimensions (Channel x FBG x Time).
        
        refwav : np.ndarray
            Reference wavelength used for comparison.
    
    If extended is True:
        dict_data : dict
            Dictionary containing channel-wise FBG data over time, indexed by channel.
        
        dict_time : dict
            Dictionary containing absolute timestamps, indexed by channel.
    
    Notes
    -----
    - If no reference wavelength (`refwav`) is provided, the first set of data is used as the reference.
    - A warning is raised if a shift greater than 3 nm is detected between the current data and the reference, 
      and the function attempts to auto-correct it using `fbg_verifier`.
    - The output 3D array `OutWav` has dimensions corresponding to (Channels, FBGs, Time).
    - The function assumes that each row in the input file contains tab-separated values, where the first column 
      is the timestamp and the following columns are wavelength values.

    Examples
    --------
    >>> filepath = 'data/fbg_measurements.txt'
    >>> absTime, refTime, OutWav, refwav = read_LS_FBG(filepath)
    
    >>> # With extended output (channel-wise dictionary)
    >>> dict_data, dict_time = read_LS_FBG(filepath, extended=True)
    
    """
    absTime = []
    refTime = []
    test = 0
    with open(filepath, mode='r', newline='') as file:
        reader = csv.reader(file)
        # Iterate over each row in the CSV file
        for i, row in tqdm(enumerate(reader)):
            if i == 0:
                continue
            lst = str(row[0]).split('\t')
            twav = np.array(np.hsplit(np.array(lst[2:]), nChannel))
            wav = np.full((nChannel, maxFBGs, 1), np.nan)
            nch, nspect = nChannel,maxFBGs
            if twav.shape[1] < maxFBGs :
                if test == 0 :
                    print('Fewer FBGs than the provided parameter maxFBGs %i' %maxFBGs)
                test +=1
                nch, nspect = nChannel,twav.shape[1]
            wav[:nch, :nspect, 0] = twav[:nChannel, :maxFBGs]
            wav = wav.astype(float)


            nbrFBG_step = np.sum(np.where(wav>0,1,0), axis=1)
            #print(nbrFBG)
            #wav = wav.reshape((wav.shape[0],wav.shape[1],1)).astype(float)
            
            if i == 1 :
                if refwav is None :
                    refwav = wav ## if no reference, designate the first value as reference. Used to check if no peak is missing
                nbrFBG = np.sum(np.where(refwav>0,1,0), axis=1)
                
            if (nbrFBG == nbrFBG_step).all() == False :
                try :
                    wav = verif_aligment(wav.reshape(nChannel, maxFBGs),refwav.reshape(nChannel, maxFBGs))
                    wav = wav.reshape((nChannel, maxFBGs, 1)).astype(float)
                except :
                    print('Error while performing the alignment of the FBG')
                    continue    
                
            wav = np.where(wav==0,np.nan,wav)
            if i == 1 : # initiate the output array
                print('Starting file : '+filepath)
                lst = str(row[0]).split('\t')
                absTime.append(lst[0]), refTime.append(float(lst[1]))
                OutWav = wav                
                continue
                
            try :
                OutWav = np.append(OutWav, wav, axis=2)
                absTime.append(lst[0]), refTime.append(float(lst[1]))
            except :
                print(lst[0], wav.shape)
                
        absTime = np.array(absTime, dtype='datetime64')
    if extended :
        dict_data = dict()
        dict_time = dict()
        for i in range(1,nChannel+1):
            dict_data['Channel_%i'%i] = (['Channel %i : FBGs'%i, 'Channel_%i_time'%i], OutWav[i-1])
            dict_time['Channel_%i_time'%i] = absTime
        return dict_data, dict_time
    return absTime, np.array(refTime), OutWav, refwav


def convertPerCH(path, nChannel=8, maxFBGs=12):
    
    if not os.path.isdir(path):
        raise ValueError(f"The provided path '{path}' does not exist or is not a directory.")
    files = glob.glob(path+'/*Z.txt')
    files.sort()
    ji = 0
    
    for file in files:
        try:
            # Read the current FBG file
            dict_data, dict_time = read(path + file, nChannel=nChannel, maxFBGs=maxFBGs, extended=True)
            mode = 'PerCH'
        except Exception as e:
            print(f"Error reading file {file}: {e}")
            continue  # Skip the file and move to the next

        if ji == 0:
            ALL_dict_data, ALL_dict_time = dict_data.copy(), dict_time.copy()
            ji += 1
            continue
        
        for i, name in enumerate(dict_data):
            ch = i + 1
            name_time = f'Channel_{ch}_time'
            
            # Merge the time arrays and the data arrays
            ALL_dict_time[name_time] = np.append(ALL_dict_time[name_time], dict_time[name_time])
            out_chann_arr = np.append(ALL_dict_data[name][1], dict_data[name][1], axis=1)
            
            ALL_dict_data[name] = (['Channel %i : FBGs' % ch, "Channel_%i_time" % ch], out_chann_arr)
            
            del out_chann_arr
        
        # If a day change is detected, process and export the data up to that point
        if verif_date(ALL_dict_time['Channel_1_time']):
            absTime = np.array(ALL_dict_time['Channel_1_time'], dtype='datetime64[ns]')
            IdxEndTime = np.where(absTime.astype('datetime64[D]') == absTime[-1].astype('datetime64[D]'))[0][0]
            
            OUT_dict_time, OUT_dict_data = dict(), dict()
            
            for ch in range(1, nChannel + 1):
                # Split the time and data arrays at the detected day change
                OUT_dict_time[f'Channel_{ch}_time'] = ALL_dict_time[f'Channel_{ch}_time'][:IdxEndTime].copy()
                OUT_dict_data[f'Channel_{ch}'] = (['Channel %i : FBGs' % ch, "Channel_%i_time" % ch], ALL_dict_data[f'Channel_{ch}'][1][:, :IdxEndTime])
                
                out_chann_arr = ALL_dict_data[f'Channel_{ch}'][1][:, IdxEndTime:].copy()
                
                ALL_dict_time[f'Channel_{ch}_time'] = ALL_dict_time[f'Channel_{ch}_time'][IdxEndTime:].copy()
                ALL_dict_data[f'Channel_{ch}'] = (['Channel %i : FBGs' % ch, "Channel_%i_time" % ch], out_chann_arr)
            
            try:
                export_nc(OUT_dict_data, OUT_dict_time, mode=mode)
                print(f"Exported data for {file} (up to day change)")
            except Exception as e:
                print(f"Error exporting data for {file}: {e}")
            
            del OUT_dict_data, OUT_dict_time, out_chann_arr
    
    # Try to export remaining data after processing all files
    try:
        export_nc(ALL_dict_data, ALL_dict_time, mode=mode)
        print(f"Final export completed for folder: {path}")
    except Exception as e:
        print(f"Error in final export: {e}")

def minimal(path, nChannel=8, maxFBGs=12):
    """
    Processes and exports data from text files in a specified directory.
    This function reads data from text files in the given directory, processes the data, 
    and exports it to NetCDF format. It handles multiple files, appends their data, 
    and detects date changes to split and export data accordingly.
    Args:
        path (str): The directory path containing the text files to process.
        nChannel (int, optional): The number of channels to process. Defaults to 8.
        maxFBGs (int, optional): The maximum number of Fiber Bragg Gratings (FBGs) to process. Defaults to 12.
    Notes:
        - The function assumes that the text files in the directory end with 'Z.txt'.
        - The function relies on external functions `read`, `verif_date`, and `export_nc` for reading, 
          verifying date changes, and exporting data, respectively.
        - If a date change is detected in the data, the function splits the data at the date boundary 
          and exports the earlier portion before continuing with the remaining data. Consequently,
          multiple NetCDF files may be created for a single directory and create 1 file per day of data. 
    Raises:
        Exception: If an error occurs while appending data from a file, the file name is printed.
    Outputs:
        - NetCDF files are exported for processed data segments, with the reference wavelength (`refwav`) included.
    """

    files = glob.glob(path+'/*Z.txt')
    files.sort()
    files = files[:-1]
    
    i=0
    for file in files:
        TabsTime, TrefTime, TOutWav, refwav = read(file, refwav=refwav, nChannel=nChannel, maxFBGs=maxFBGs, extended=False)
        if i == 0 :
            absTime, refTime, OutWav = TabsTime, TrefTime, TOutWav
            i+=1 
            continue
            
        print('Extracted data shapes :', TabsTime.shape, TrefTime.shape, TOutWav.shape)
        try :
            absTime, refTime, OutWav = np.append(absTime, TabsTime), np.append(refTime, TrefTime), np.append(OutWav, TOutWav, axis=2)
        except :
            print('Error file : %s' %file)
        #print(refwav)
        
        if verif_date(absTime) == True :
            print('Date change detected')
            IdxEndTime = np.where(absTime.astype('datetime64[D]')==absTime[-1].astype('datetime64[D]'))[0][0]
            OabsTime, OrefTime, OOutWav = absTime[:IdxEndTime], refTime[:IdxEndTime], OutWav[:,:,:IdxEndTime]
            absTime, refTime, OutWav = absTime[IdxEndTime:], refTime[IdxEndTime:], OutWav[:,:,IdxEndTime:]
            
            export_nc(OabsTime, OrefTime, OOutWav, refwav=refwav)
    export_nc(absTime, refTime, OutWav, refwav=refwav)



def convert2nc(path, nChannel=8, maxFBGs=12, extended=False):
    """
    Converts data from the specified directory into a NetCDF format.
    Parameters:
        path (str): The path to the directory containing the data to be converted.
                    Must be a valid and existing directory.
        nChannel (int, optional): The number of channels set in the BSI parameters. Defaults to 8.
        maxFBGs (int, optional): The maximum number of Fiber Bragg Gratings (FBGs) to process. Defaults to 12.
        extended (bool, optional): Determines the conversion method to use. If True, uses the extended 
                                   conversion method (`convertPerCH`). If False, uses the minimal 
                                   conversion method (`minimal`). Defaults to False.
    Raises:
        ValueError: If the provided `path` is not a valid directory or does not exist.
    Notes:
        - The `convertPerCH` function is used for extended conversion. It will produce 2D arrays for each channel.
        - The `minimal` function is used for minimal conversion. It will produce 3D arrays for all channels.
    """

    # Determine which conversion method to use based on the `extended` parameter.
    if not os.path.isdir(path):
        raise ValueError(f"The provided path '{path}' is not a valid directory.")
    if not os.path.isdir(path):
        raise ValueError(f"The provided path '{path}' does not exist or is not a directory.")
    
    if extended:
        convertPerCH(path, nChannel=nChannel, maxFBGs=maxFBGs) 
    else: 
        minimal(path, nChannel=nChannel, maxFBGs=maxFBGs)     
    