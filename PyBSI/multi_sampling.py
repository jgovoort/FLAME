import pandas as pd
import numpy as np
import warnings
import os

from PyBSI.tools import *

def bsi_hsv_dual_161_reader(filename, maxFBGs=8, nChannel=8):
    """
    Reads BSI files in dual mode (version 1.6.1) with high and variable sampling frequencies per channel.
    
    Parameters
    ----------
    filename : str
        Path to the BSI data file.
        
    maxFBGs : int, optional
        Maximum number of FBGs per channel (default is 8).
        
    nChannel : int, optional
        Maximum number of channels to process (default is 8).
    
    Returns
    -------
    dict_data : dict
        Dictionary containing the data for each channel.
        
    dict_time : dict
        Dictionary containing the corresponding time arrays for each channel.
    """
    try:
        # Read the file into a pandas DataFrame
        data_bsi = pd.read_csv(filename, header=None, skiprows=1, sep='\t', decimal=',', 
                               names=np.append(np.array(['Time', 'Channel']), np.arange(1, maxFBGs+1)))
    except Exception as e:
        raise ValueError(f"Error reading file {filename}: {e}")
    
    dict_data, dict_time = {}, {}
    
    # Extract reference time and absolute time
    refTime = data_bsi['Time']
    try:
        absTime_start = np.datetime64('%s-%s-%sT%s:%s:%s' % (
            filename[-24:-20], filename[-20:-18], filename[-18:-16],
            filename[-15:-13], filename[-13:-11], filename[-11:-5].replace(',', '.')), 'ms')
    except ValueError as e:
        raise ValueError(f"Error parsing timestamp from filename {filename}: {e}")
    
    if refTime.empty:
        warnings.warn(f"File {filename} is empty!")
        for ch in range(1, nChannel+1):
            dict_data[f'Channel_{ch}'] = (['Channel %i : FBGs' % ch, "Channel %i : Over time" % ch], np.empty((maxFBGs, 0)))
            dict_time[f'Channel_{ch}_time'] = np.empty((0,))
        return dict_data, dict_time

    # Calculate absolute time based on refTime
    dtime = (refTime - refTime.iloc[0]) * 1000  # Convert to milliseconds
    absTime = absTime_start + np.array(dtime, dtype='timedelta64[ms]')
    
    for ch in range(1, nChannel+1):
        cond = data_bsi['Channel'] == ch
        if cond.sum() == 0:
            warnings.warn(f"No data or no existing channel {ch}")
            continue
        
        channel_data = data_bsi[cond]
        if channel_data.empty:
            warnings.warn(f"No data for Channel {ch} in file {filename}")
            continue
        
        # Initialize arrays for time and data
        out_chann_time = np.array(absTime[cond])
        out_chann_arr = np.empty((maxFBGs, len(channel_data)))
        
        for fbg in range(maxFBGs):
            data_ch = np.array(channel_data.iloc[:, fbg+2])
            median_fbg = np.nanmedian(data_ch)
            
            # Clean data by removing outliers (values beyond +/- 1 nm from the median)
            clean_fbg = np.where((data_ch > median_fbg + 1) | (data_ch < median_fbg - 1), np.nan, data_ch)
            out_chann_arr[fbg] = clean_fbg
        
        dict_data[f'Channel_{ch}'] = (['Channel %i : FBGs' % ch, "Channel %i : Over time" % ch], out_chann_arr)
        dict_time[f'Channel_{ch}_time'] = out_chann_time

    return dict_data, dict_time

def folder_reader_DUAL_FBG(path, nChannel=8, maxFBGs=12):
    """
    Reads multiple BSI files in dual mode from a specified directory and compiles the data.
    
    Parameters
    ----------
    path : str
        Path to the directory containing the BSI data files.
        
    nChannel : int, optional
        Maximum number of channels to process (default is 8).
        
    maxFBGs : int, optional
        Maximum number of FBGs per channel (default is 12).
    
    Returns
    -------
    None
    """
    files = os.listdir(path)
    files.sort()
    
    ALL_dict_data, ALL_dict_time = {}, {}
    
    for file in files:
        if not file.endswith('Z.txt'):
            continue
            
        try:
            dict_data, dict_time = bsi_hsv_dual_161_reader(os.path.join(path, file), maxFBGs=maxFBGs, nChannel=nChannel)
        except Exception as e:
            print(f"Error reading file {file}: {e}")
            continue  # Skip this file and continue with the next
        
        mode = 'DUAL'
        
        # Initialize ALL_dict_data and ALL_dict_time on the first iteration
        if not ALL_dict_data:
            ALL_dict_data, ALL_dict_time = dict_data.copy(), dict_time.copy()
            continue
            
        # Append new data to existing data structures
        for i, name in enumerate(dict_data):
            ch = i + 1
            name_time = f'{name}_time'
            
            # Append time data
            ALL_dict_time[name_time] = np.append(ALL_dict_time[name_time], dict_time[name_time])
            
            # Append channel data
            ALL_dict_data[name] = (['Channel %i : FBGs' % ch, "Channel %i : Over time" % ch], 
                                    np.append(ALL_dict_data[name][1], dict_data[name][1], axis=1))
        
        # Check if the date is valid and handle data accordingly
        if verif_date(ALL_dict_time['Channel_1_time']):
            absTime = np.array(ALL_dict_time['Channel_1_time'], dtype='datetime64')
            IdxEndTime = np.where(absTime.astype('datetime64[D]') == absTime[-1].astype('datetime64[D]'))[0][0]

            # Create output dictionaries
            OUT_dict_time = {name: ALL_dict_time[name][:IdxEndTime].copy() for name in ALL_dict_time}
            OUT_dict_data = {name: (ALL_dict_data[name][0], ALL_dict_data[name][1][:, :IdxEndTime]) for name in ALL_dict_data}

            # Export the data
            export_nc(OUT_dict_data, OUT_dict_time, mode=mode)
            
            # Update ALL_dicts for remaining data
            for i, name in enumerate(dict_data):
                ch = i + 1
                if ch == 1:
                    ALL_dict_time[name_time] = ALL_dict_time[name_time][IdxEndTime:].copy()
                    ALL_dict_data[name] = (ALL_dict_data[name][0], ALL_dict_data[name][1][:, IdxEndTime:])
                else:
                    ALL_dict_time[name_time] = np.array([])
                    ALL_dict_data[name] = (['Channel %i : FBGs' % ch, "Channel %i : Over time" % ch], np.empty((maxFBGs, 0)))

    # Final export of all collected data
    try:
        export_nc(ALL_dict_data, ALL_dict_time, mode=mode)
    except Exception as e:
        print(f"Error exporting final data: {e}")