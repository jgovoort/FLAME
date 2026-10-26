# PyBSI

PyBSI contains Python utilities for reading Fiber Bragg Grating (FBG) files
exported by BSI interrogators and converting them into analysis-friendly
structures or NetCDF files.

The package currently targets two BSI data layouts:

- **Low sampling (`1.3`)**: multichannel text files containing wavelength
  measurements for several FBGs per channel.
- **Dual sampling (`1.6.1DUAL`)**: files with a channel identifier and
  channel-dependent sampling times.

The package version defined in [`__version__.py`](__version__.py) is `0.2.0`.

## Installation

PyBSI is currently used directly from this repository and does not yet include
packaging metadata such as `pyproject.toml` or `setup.py`. From the repository
root, install its runtime dependencies and keep the repository root on
`PYTHONPATH` by running scripts or notebooks from the root directory:

```bash
python -m pip install numpy pandas xarray
```

The low-sampling reader uses `tqdm` for progress reporting:

```bash
python -m pip install tqdm
```

## Low-sampling files

Use `PyBSI.low_sampling.read` to read one BSI text file. The reader returns
absolute timestamps, reference times, a wavelength array, and the reference
wavelength array:

```python
from PyBSI.low_sampling import read

absolute_time, reference_time, wavelengths, reference_wavelength = read(
    "path/to/measurement.txt",
    nChannel=8,
    maxFBGs=12,
)

print(wavelengths.shape)  # (channels, FBGs, time)
```

If `refwav` is omitted, the first measurement is used as the reference. To
retain channel-specific dictionaries and timestamps, set `extended=True`:

```python
from PyBSI.low_sampling import read

channel_data, channel_time = read(
    "path/to/measurement.txt",
    nChannel=8,
    maxFBGs=12,
    extended=True,
)
```

For a directory containing files ending in `Z.txt`, the module also provides
`convert2nc`. With `extended=False` it selects the compact conversion path;
with `extended=True` it selects the channel-oriented path:

```python
from PyBSI.low_sampling import convert2nc

convert2nc(
    "path/to/low_sampling_directory",
    nChannel=8,
    maxFBGs=12,
    extended=True,
)
```

These helpers are intended to split data at day boundaries and write NetCDF
files through `PyBSI.tools.export_nc`. Validate the input naming and output
paths against the target BSI export before processing a large directory.

## Dual-sampling files

For a single dual-sampling file, use
`PyBSI.multi_sampling.bsi_hsv_dual_161_reader`:

```python
from PyBSI.multi_sampling import bsi_hsv_dual_161_reader

channel_data, channel_time = bsi_hsv_dual_161_reader(
    "path/to/dual_measurement.txt",
    nChannel=8,
    maxFBGs=12,
)
```

To process a directory of dual-sampling files:

```python
from PyBSI.multi_sampling import folder_reader_DUAL_FBG

folder_reader_DUAL_FBG(
    "path/to/dual_sampling_directory",
    nChannel=8,
    maxFBGs=12,
)
```

Dual-sampling filenames must contain the timestamp pattern expected by the
reader. Files are read as tab-separated data with a time column, a channel
column, and one column per FBG.

## Utility functions

`PyBSI.tools` provides small helpers used by the readers:

- `find_nearest(array, value)` returns the index of the closest array value.
- `verif_date(absTime)` detects whether timestamps span more than one calendar
  day.
- `export_nc(...)` writes wavelength data and calibration values to NetCDF.

The resulting datasets use wavelength values in nanometres and include time
coordinates. Low-sampling exports use dimensions for channel, FBG, and time;
per-channel and dual-sampling exports retain channel-specific time
coordinates.

## Module overview

| Module | Main purpose |
| --- | --- |
| `low_sampling.py` | Read and convert low-sampling multichannel BSI files. |
| `multi_sampling.py` | Read and combine dual-sampling BSI files. |
| `per_channel.py` | Channel-oriented conversion helpers. |
| `spectrum.py` | Combine single-channel spectrum text files. |
| `tools.py` | Timestamp, nearest-value, and NetCDF export utilities. |
| `HSV.py` | Legacy reader for HSV-style FBG files. |

## Input and output notes

- `nChannel` and `maxFBGs` must match the acquisition configuration, or missing
  FBG values may be padded with `NaN`.
- Low-sampling input is expected to contain a timestamp, a reference-time
  value, and tab-separated wavelength values.
- Output filenames are generated from the first timestamp in the data.
- NetCDF output is written to the current directory unless a different output
  path is passed to the export helper.
- The readers are designed for the formats used by the FLAME project; validate
  timestamps and channel ordering when using data from another BSI export.

## License

PyBSI is part of FLAME and is distributed under the GNU General Public License
v3.0. See the repository [`LICENSE`](../LICENSE).