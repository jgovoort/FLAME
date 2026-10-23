import PyBSI.low_sampling as ls
import PyBSI.multi_sampling as ms

def convert(path, folder=True, refwav=None, version='1.3', nChannel=8, maxFBGs=12):
    """
    Convert FBG data from different versions of BSI software into NetCDF files.
    
    Parameters
    ----------
    path : str
        Path to the file or folder containing the FBG data.
    
    folder : bool, optional
        If True, processes a folder of files. If False, processes a single file. Default is True.
    
    refwav : np.ndarray, optional
        Reference wavelength for comparison. Default is None.
    
    version : str, optional
        Version of the BSI software. Options are:
            - '1.3' for low sampling frequency data (default).
            - '1.6.1DUAL' for data with dual interrogation times per channel.
    
    nChannel : int, optional
        Number of channels to process. Default is 8.
    
    maxFBGs : int, optional
        Maximum number of FBGs per channel. Default is 12.
    
    Returns
    -------
    None
    """
    
    try:
        # If a single file is to be processed
        if not folder:
            if version == '1.3':
                print(f"Processing single file with version {version} (Low Sampling)...")
                dict_data, dict_time = ls.read(path, nChannel=nChannel, maxFBGs=maxFBGs, extended=True)
                mode = 'PerCH'
            elif version == '1.6.1DUAL':
                print(f"Processing single file with version {version} (Dual Sampling)...")
                dict_data, dict_time = ms.bsi_hsv_dual_161_reader(path, maxFBGs=maxFBGs, nChannel=nChannel)
                mode = 'PerCH'
            else:
                raise ValueError(f"Unsupported version: {version}")
            
            # Export the data to NetCDF
            export_nc(dict_data, dict_time, mode=mode)
            print(f"Export completed for file: {path}")
        
        # If a folder of files is to be processed
        else:
            if version == '1.3':
                print(f"Processing folder with version {version} (Low Sampling)...")
                ls.folder_read(path, nChannel=nChannel, maxFBGs=maxFBGs)
            elif version == '1.6.1DUAL':
                print(f"Processing folder with version {version} (Dual Sampling)...")
                ms.folder_reader_DUAL_FBG(path, nChannel=nChannel, maxFBGs=maxFBGs)
            else:
                raise ValueError(f"Unsupported version: {version}")
        
        print(f"Processing completed for path: {path}")
    
    except Exception as e:
        print(f"An error occurred while processing {path}: {e}")