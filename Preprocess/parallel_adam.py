import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import splat
import io
import contextlib
import re

normal_list = ['M0.0','M1.0','M2.0','M3.0','M4.0','M5.0','M6.0', 'M7.0', 'M8.0', 'M9.0', 'L0.0', 'L1.0', 'L2.0', 'L3.0', 'L4.0',
                  'L5.0', 'L6.0', 'L7.0', 'L8.0', 'L9.0', 'T0.0', 'T1.0', 'T2.0', 'T3.0', 'T4.0', 'T5.0', 'T6.0', 'T7.0', 
               'T8.0', 'T9.0']
dsd_list = ['d/sdM4.0','d/sdM5.0','d/sdM6.0','d/sdM7.0','d/sdM8.0','d/sdM9.0','d/sdL0.0','d/sdL1.0','d/sdL7.0']
sd_list = ['sdM2.0','sdM4.0','sdM5.0','sdM5.5','sdM6.0','sdM7.0','sdM8.0','sdM9.5','sdL0.0','sdL3.5','sdL4.0']
esd_list = ['esdM0.0','esdM4.0','esdM5.0','esdM6.5','esdM7.5','esdM8.5']
vlg_list = ['M6.0gamma','M7.0gamma','M8.0gamma','M9.0gamma','L0.0gamma','L1.0gamma','L2.0gamma','L3.0gamma','L4.0gamma','L6.0gamma']
intg_list = ['M8.0beta','L0.0beta','L1.0beta','L2.0beta','L3.0beta']

spectral_types = sd_list+esd_list+vlg_list+intg_list+normal_list+dsd_list

def assign_dwarf_type(spectral_type):
    if spectral_type in normal_list:
        return 0
    elif spectral_type in dsd_list:
        return -1
    elif spectral_type in sd_list:
        return -2
    elif spectral_type in esd_list:
        return -3
    elif spectral_type in vlg_list:
        return 1
    elif spectral_type in intg_list:
        return 2
    else:
        return -99  # return some default value for spectral types not included in your lists

wav_002 = [(0.85,0.87),(0.87,0.89),(0.89, 0.91),(0.91, 0.93),(0.93, 0.95),(0.95, 0.97),(0.97, 0.99),(0.99, 1.01),
                       (1.01, 1.03),(1.03, 1.05),(1.05, 1.07),(1.07, 1.09),(1.09, 1.11), (1.11, 1.13),(1.13,1.15),(1.15, 1.17),
                       (1.17, 1.19),(1.19, 1.21),(1.21, 1.23),(1.23, 1.25),(1.25, 1.27),(1.27, 1.29),(1.29, 1.31),(1.31, 1.33),
                       (1.33, 1.35),(1.42,1.44),(1.44,1.46),(1.46,1.48),(1.48,1.50),(1.50,1.52), (1.52,1.54),(1.54,1.56),(1.56,1.58),
                       (1.58,1.60),(1.60,1.62),(1.62,1.64),(1.64,1.66),(1.66,1.68),(1.68,1.7),(1.7,1.72),(1.72,1.74),(1.74,1.76),
                       (1.76,1.78),(1.78,1.80),(1.95,1.97),(1.97,1.99),(1.99,2.01),(2.01,2.03),(2.03,2.05),(2.05,2.07),(2.07,2.09),
                       (2.09,2.11),(2.11,2.13),(2.13,2.15),(2.15,2.17),(2.17,2.19),(2.19,2.21),(2.21,2.23),(2.23,2.25),(2.25,2.27),
                       (2.27,2.29),(2.29,2.31),(2.31,2.33),(2.33,2.35),(2.35,2.37),(2.37,2.39),(2.39,2.41),(2.41,2.43),(2.43,2.45)]

# +
# import io
# import contextlib
# import re
# import numpy as np
# import splat

# # Initialize your spectrum and standards
# sp = splat.Spectrum(file='spex-prism_twa8b.fits', instrument='SPEX-PRISM', name='TWA 8B (M6.0gamma Std)')
# splat.initiateStandards(vlg=True)
# splat.STDS_VLG_SPEX['M6.0gamma'] = sp

# def process_spectrum(name, wav_ranges=None):
#     """
#     Process a spectrum given by 'name' and compute the mean flux in each wavelength bin.
    
#     Parameters:
#         name (str): The name identifier for the spectrum.
#         wav_ranges (list of tuple, optional): A list of wavelength intervals (start, end) over which to compute the flux.
#             If not provided, defaults to the wav_002 bins.
    
#     Returns:
#         tuple: Two lists:
#             - data: Contains tuples of (spectrum name, best-fit spectral type, chi-square, SNR).
#             - flux_means: Contains lists of [spectrum name, spectral type, flux_mean for each wavelength bin].
#     """
#     # Use default wavelength bins if none are provided
#     if wav_ranges is None:
#         wav_ranges = wav_002

#     data = []
#     flux_means = []

#     try:
#         spectra_list = splat.getSpectrum(name=name)
#         if not spectra_list:
#             return None, None  # No spectra found for this name

#         best_spectrum = None
#         highest_snr = -1
#         chi_squares = []

#         # Select the spectrum with the highest SNR
#         for sp in spectra_list:
#             snr = sp.snr  # Assuming SNR is an attribute; adjust if necessary
#             if snr > highest_snr:
#                 highest_snr = snr
#                 best_spectrum = sp

#         if best_spectrum is None:
#             return None, None

#         # Classify the best spectrum (redirecting output to capture classification info)
#         buffer = io.StringIO()
#         with contextlib.redirect_stdout(buffer):
#             classification = splat.classifyByStandard(best_spectrum, fit_ranges=[0.87, 2.39], verbose=True, all=True)
            
#         # Normalize the best spectrum
#         best_spectrum.normalize([1.27, 1.28])
        
#         output = buffer.getvalue()
#         matches = re.findall(r'Type (.*): statistic = (.*?),', output)
#         statistics = [(t, float(stat)) for t, stat in matches]
#         best_fit = min(statistics, key=lambda x: x[1])
#         chi_squares.append(best_fit[1])

#         # Ensure spectral_types is defined in your environment before using it
#         if best_fit[0] in spectral_types:
#             data.append((best_spectrum.name, best_fit[0], best_fit[1], highest_snr))
#             spectral_type = classification[0]
#             fluxes = []
#             for wl in wav_ranges:
#                 mask = (best_spectrum.wave.value >= wl[0]) & (best_spectrum.wave.value <= wl[1])
#                 if np.any(mask):
#                     flux_mean = np.mean(best_spectrum.flux.value[mask])
#                 else:
#                     flux_mean = np.nan  # or handle missing data as needed
#                 fluxes.append(flux_mean)
#             flux_means.append([best_spectrum.name, spectral_type] + fluxes)
    
#     except Exception as e:
#         print(f"An error occurred with {name}. Error: {e}. Skipping this spectrum.")

#     return data, flux_means
# -

wav_004 = [
    (0.85, 0.89),
    (0.89, 0.93),
    (0.93, 0.97),
    (0.97, 1.01),
    (1.01, 1.05),
    (1.05, 1.09),
    (1.09, 1.13),
    (1.13, 1.17),
    (1.17, 1.21),
    (1.21, 1.25),
    (1.25, 1.29),
    (1.29, 1.33),
    
    (1.42, 1.46),
    (1.46, 1.50),
    (1.50, 1.54),
    (1.54, 1.58),
    (1.58, 1.62),
    (1.62, 1.66),
    (1.66, 1.70),
    (1.70, 1.74),
    (1.74, 1.78),
    
    (1.95, 1.99),
    (1.99, 2.03),
    (2.03, 2.07),
    (2.07, 2.11),
    (2.11, 2.15),
    (2.15, 2.19),
    (2.19, 2.23),
    (2.23, 2.27),
    (2.27, 2.31),
    (2.31, 2.35),
    (2.35, 2.39),
    (2.39, 2.43)
]

wav_001 = [(0.85, 0.86),
 (0.86, 0.87),(0.87, 0.88),(0.88, 0.89),(0.89, 0.9), (0.9, 0.91),(0.91, 0.92), (0.92, 0.93), (0.93, 0.94), (0.94, 0.95),
 (0.95, 0.96), (0.96, 0.97), (0.97, 0.98), (0.98, 0.99), (0.99, 1.0), (1.0, 1.01), (1.01, 1.02),(1.02, 1.03), (1.03, 1.04),
(1.04, 1.05), (1.05, 1.06), (1.06, 1.07), (1.07, 1.08), (1.08, 1.09), (1.09, 1.1), (1.1, 1.11), (1.11, 1.12), (1.12, 1.13),
 (1.13, 1.14), (1.14, 1.15), (1.15, 1.16),(1.16, 1.17), (1.17, 1.18), (1.18, 1.19), (1.19, 1.2), (1.2, 1.21), (1.21, 1.22),
 (1.22, 1.23), (1.23, 1.24), (1.24, 1.25), (1.25, 1.26), (1.26, 1.27),(1.27, 1.28), (1.28, 1.29), (1.29, 1.3), (1.3, 1.31), (1.31, 1.32), (1.32, 1.33),
 (1.33, 1.34), (1.34, 1.35),(1.42, 1.43), (1.43, 1.44), (1.44, 1.45), (1.45, 1.46), (1.46, 1.47), (1.47, 1.48), (1.48, 1.49),
 (1.49, 1.5), (1.5, 1.51), (1.51, 1.52), (1.52, 1.53), (1.53, 1.54), (1.54, 1.55), (1.55, 1.56), (1.56, 1.57), (1.57, 1.58),
 (1.58, 1.59), (1.59, 1.6), (1.6, 1.61), (1.61, 1.62), (1.62, 1.63), (1.63, 1.64), (1.64, 1.65), (1.65, 1.66), (1.66, 1.67),
 (1.67, 1.68), (1.68, 1.69), (1.69, 1.7), (1.7, 1.71), (1.71, 1.72), (1.72, 1.73),(1.73, 1.74), (1.74, 1.75), (1.75, 1.76),
(1.76, 1.77), (1.77, 1.78), (1.78, 1.79), (1.79, 1.8), (1.95, 1.96), (1.96, 1.97), (1.97, 1.98), (1.98, 1.99), (1.99, 2.0),
 (2.0, 2.01), (2.01, 2.02), (2.02, 2.03), (2.03, 2.04), (2.04, 2.05), (2.05, 2.06), (2.06, 2.07), (2.07, 2.08),(2.08, 2.09),
 (2.09, 2.1), (2.1, 2.11), (2.11, 2.12), (2.12, 2.13),(2.13, 2.14), (2.14, 2.15), (2.15, 2.16), (2.16, 2.17), (2.17, 2.18),
 (2.18, 2.19), (2.19, 2.2), (2.2, 2.21), (2.21, 2.22), (2.22, 2.23), (2.23, 2.24), (2.24, 2.25), (2.25, 2.26), (2.26, 2.27),
 (2.27, 2.28), (2.28, 2.29), (2.29, 2.3), (2.3, 2.31),(2.31, 2.32), (2.32, 2.33), (2.33, 2.34), (2.34, 2.35), (2.35, 2.36),
 (2.36, 2.37), (2.37, 2.38), (2.38, 2.39), (2.39, 2.4), (2.4, 2.41), (2.41, 2.42), (2.42, 2.43), (2.43, 2.44), (2.44, 2.45)]

wav_010 = [
    # Segment 1
    (0.85, 0.95),
    (0.95, 1.05),
    (1.05, 1.15),
    (1.15, 1.25),
    (1.25, 1.35),
    
    # Segment 2
    (1.42, 1.52),
    (1.52, 1.62),
    (1.62, 1.72),
    (1.72, 1.82),
    
    # Segment 3
    (1.95, 2.05),
    (2.05, 2.15),
    (2.15, 2.25),
    (2.25, 2.35),
    (2.35, 2.45)
]


wav_003 = [
    # Range 1: 0.85 to 1.35
    (0.85, 0.88), (0.88, 0.91), (0.91, 0.94), (0.94, 0.97),
    (0.97, 1.0), (1.0, 1.03), (1.03, 1.06), (1.06, 1.09),
    (1.09, 1.12), (1.12, 1.15), (1.15, 1.18), (1.18, 1.21),
    (1.21, 1.24), (1.24, 1.27), (1.27, 1.3), (1.3, 1.33),
    (1.33, 1.35),
    
    # Range 2: 1.42 to 1.80 (skipping 1.35 to 1.42)
    (1.42, 1.45), (1.45, 1.48), (1.48, 1.51), (1.51, 1.54),
    (1.54, 1.57), (1.57, 1.6), (1.6, 1.63), (1.63, 1.66),
    (1.66, 1.69), (1.69, 1.72), (1.72, 1.75), (1.75, 1.78),
    (1.78, 1.8),
    
    # Range 3: 1.95 to 2.46 (skipping 1.80 to 1.95)
    (1.95, 1.98), (1.98, 2.01), (2.01, 2.04), (2.04, 2.07),
    (2.07, 2.1), (2.1, 2.13), (2.13, 2.16), (2.16, 2.19),
    (2.19, 2.22), (2.22, 2.25), (2.25, 2.28), (2.28, 2.31),
    (2.31, 2.34), (2.34, 2.37), (2.37, 2.4), (2.4, 2.43),
    (2.43, 2.46)
]

wav_005 = [
    # Segment 1
    (0.85, 0.90),
    (0.90, 0.95),
    (0.95, 1.00),
    (1.00, 1.05),
    (1.05, 1.10),
    (1.10, 1.15),
    (1.15, 1.20),
    (1.20, 1.25),
    (1.25, 1.30),
    (1.30, 1.35),
    
    # Segment 2
    (1.42, 1.47),
    (1.47, 1.52),
    (1.52, 1.57),
    (1.57, 1.62),
    (1.62, 1.67),
    (1.67, 1.72),
    (1.72, 1.77),
    
    # Segment 3
    (1.95, 2.00),
    (2.00, 2.05),
    (2.05, 2.10),
    (2.10, 2.15),
    (2.15, 2.20),
    (2.20, 2.25),
    (2.25, 2.30),
    (2.30, 2.35),
    (2.35, 2.40),
    (2.40, 2.45)
]


# +
import io
import contextlib
import re
import numpy as np
import splat

sp = splat.Spectrum(file='spex-prism_twa8b.fits',instrument='SPEX-PRISM',name='TWA 8B (M6.0gamma Std)')
splat.initiateStandards(vlg=True)
splat.STDS_VLG_SPEX['M6.0gamma'] = sp

def process_spectrum_001(name):
    # Initialize an empty list to hold the data
    data = []
    flux_means = []

    try:
        spectra_list = splat.getSpectrum(name=name)
        if not spectra_list:
            return None, None  # Return or handle empty spectra_list as needed

        best_spectrum = None
        highest_snr = -1
        chi_squares = []

        # Iterate through the spectra to find the one with the highest SNR
        for sp in spectra_list:
            snr = sp.snr  # Assuming SNR can be accessed directly; adjust if necessary
            if snr > highest_snr:
                highest_snr = snr
                best_spectrum = sp

        if best_spectrum is None:
            return None, None  # Handle case where no valid spectrum was found

        # Perform classification on the best spectrum
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            classification = splat.classifyByStandard(best_spectrum, fit_ranges=[0.87,2.39], verbose=True, all=True)
            
        # Normalize the best spectrum
        best_spectrum.normalize([1.27, 1.28])
        
        output = buffer.getvalue()
        matches = re.findall(r'Type (.*): statistic = (.*?),', output)
        statistics = [(type_, float(statistic)) for type_, statistic in matches]
        best_fit = min(statistics, key=lambda x: x[1])
        chi_squares.append(best_fit[1])

        if best_fit[0] in spectral_types:
            data.append((best_spectrum.name, best_fit[0], best_fit[1], highest_snr))
            spectral_type = classification[0]
            fluxes = [np.mean(best_spectrum.flux.value[np.where((best_spectrum.wave.value >= wl[0]) & (best_spectrum.wave.value <= wl[1]))[0]]) for wl in wav_001]
            flux_means.append([best_spectrum.name, spectral_type] + fluxes)
    
    except Exception as e:
        print(f"An error occurred with {name}. Error: {e}. Skipping this spectrum.")

    return data, flux_means


# +
import io
import contextlib
import re
import numpy as np
import splat

sp = splat.Spectrum(file='spex-prism_twa8b.fits',instrument='SPEX-PRISM',name='TWA 8B (M6.0gamma Std)')
splat.initiateStandards(vlg=True)
splat.STDS_VLG_SPEX['M6.0gamma'] = sp

def process_spectrum_010(name):
    # Initialize an empty list to hold the data
    data = []
    flux_means = []

    try:
        spectra_list = splat.getSpectrum(name=name)
        if not spectra_list:
            return None, None  # Return or handle empty spectra_list as needed

        best_spectrum = None
        highest_snr = -1
        chi_squares = []

        # Iterate through the spectra to find the one with the highest SNR
        for sp in spectra_list:
            snr = sp.snr  # Assuming SNR can be accessed directly; adjust if necessary
            if snr > highest_snr:
                highest_snr = snr
                best_spectrum = sp

        if best_spectrum is None:
            return None, None  # Handle case where no valid spectrum was found

        # Perform classification on the best spectrum
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            classification = splat.classifyByStandard(best_spectrum, fit_ranges=[0.87,2.39], verbose=True, all=True)
            
        # Normalize the best spectrum
        best_spectrum.normalize([1.27, 1.28])
        
        output = buffer.getvalue()
        matches = re.findall(r'Type (.*): statistic = (.*?),', output)
        statistics = [(type_, float(statistic)) for type_, statistic in matches]
        best_fit = min(statistics, key=lambda x: x[1])
        chi_squares.append(best_fit[1])

        if best_fit[0] in spectral_types:
            data.append((best_spectrum.name, best_fit[0], best_fit[1], highest_snr))
            spectral_type = classification[0]
            fluxes = [np.mean(best_spectrum.flux.value[np.where((best_spectrum.wave.value >= wl[0]) & (best_spectrum.wave.value <= wl[1]))[0]]) for wl in wav_010]
            flux_means.append([best_spectrum.name, spectral_type] + fluxes)
    
    except Exception as e:
        print(f"An error occurred with {name}. Error: {e}. Skipping this spectrum.")

    return data, flux_means


# +
import io
import contextlib
import re
import numpy as np
import splat

sp = splat.Spectrum(file='spex-prism_twa8b.fits',instrument='SPEX-PRISM',name='TWA 8B (M6.0gamma Std)')
splat.initiateStandards(vlg=True)
splat.STDS_VLG_SPEX['M6.0gamma'] = sp

def process_spectrum_005(name):
    # Initialize an empty list to hold the data
    data = []
    flux_means = []

    try:
        spectra_list = splat.getSpectrum(name=name)
        if not spectra_list:
            return None, None  # Return or handle empty spectra_list as needed

        best_spectrum = None
        highest_snr = -1
        chi_squares = []

        # Iterate through the spectra to find the one with the highest SNR
        for sp in spectra_list:
            snr = sp.snr  # Assuming SNR can be accessed directly; adjust if necessary
            if snr > highest_snr:
                highest_snr = snr
                best_spectrum = sp

        if best_spectrum is None:
            return None, None  # Handle case where no valid spectrum was found

        # Perform classification on the best spectrum
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            classification = splat.classifyByStandard(best_spectrum, fit_ranges=[0.87,2.39], verbose=True, all=True)
            
        # Normalize the best spectrum
        best_spectrum.normalize([1.27, 1.28])
        
        output = buffer.getvalue()
        matches = re.findall(r'Type (.*): statistic = (.*?),', output)
        statistics = [(type_, float(statistic)) for type_, statistic in matches]
        best_fit = min(statistics, key=lambda x: x[1])
        chi_squares.append(best_fit[1])

        if best_fit[0] in spectral_types:
            data.append((best_spectrum.name, best_fit[0], best_fit[1], highest_snr))
            spectral_type = classification[0]
            fluxes = [np.mean(best_spectrum.flux.value[np.where((best_spectrum.wave.value >= wl[0]) & (best_spectrum.wave.value <= wl[1]))[0]]) for wl in wav_005]
            flux_means.append([best_spectrum.name, spectral_type] + fluxes)
    
    except Exception as e:
        print(f"An error occurred with {name}. Error: {e}. Skipping this spectrum.")

    return data, flux_means


# +
import io
import contextlib
import re
import numpy as np
import splat

sp = splat.Spectrum(file='spex-prism_twa8b.fits',instrument='SPEX-PRISM',name='TWA 8B (M6.0gamma Std)')
splat.initiateStandards(vlg=True)
splat.STDS_VLG_SPEX['M6.0gamma'] = sp

def process_spectrum_004(name):
    # Initialize an empty list to hold the data
    data = []
    flux_means = []

    try:
        spectra_list = splat.getSpectrum(name=name)
        if not spectra_list:
            return None, None  # Return or handle empty spectra_list as needed

        best_spectrum = None
        highest_snr = -1
        chi_squares = []

        # Iterate through the spectra to find the one with the highest SNR
        for sp in spectra_list:
            snr = sp.snr  # Assuming SNR can be accessed directly; adjust if necessary
            if snr > highest_snr:
                highest_snr = snr
                best_spectrum = sp

        if best_spectrum is None:
            return None, None  # Handle case where no valid spectrum was found

        # Perform classification on the best spectrum
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            classification = splat.classifyByStandard(best_spectrum, fit_ranges=[0.87,2.39], verbose=True, all=True)
            
        # Normalize the best spectrum
        best_spectrum.normalize([1.27, 1.28])
        
        output = buffer.getvalue()
        matches = re.findall(r'Type (.*): statistic = (.*?),', output)
        statistics = [(type_, float(statistic)) for type_, statistic in matches]
        best_fit = min(statistics, key=lambda x: x[1])
        chi_squares.append(best_fit[1])

        if best_fit[0] in spectral_types:
            data.append((best_spectrum.name, best_fit[0], best_fit[1], highest_snr))
            spectral_type = classification[0]
            fluxes = [np.mean(best_spectrum.flux.value[np.where((best_spectrum.wave.value >= wl[0]) & (best_spectrum.wave.value <= wl[1]))[0]]) for wl in wav_004]
            flux_means.append([best_spectrum.name, spectral_type] + fluxes)
    
    except Exception as e:
        print(f"An error occurred with {name}. Error: {e}. Skipping this spectrum.")

    return data, flux_means


# +
import io
import contextlib
import re
import numpy as np
import splat

sp = splat.Spectrum(file='spex-prism_twa8b.fits',instrument='SPEX-PRISM',name='TWA 8B (M6.0gamma Std)')
splat.initiateStandards(vlg=True)
splat.STDS_VLG_SPEX['M6.0gamma'] = sp

def process_spectrum_003(name):
    # Initialize an empty list to hold the data
    data = []
    flux_means = []

    try:
        spectra_list = splat.getSpectrum(name=name)
        if not spectra_list:
            return None, None  # Return or handle empty spectra_list as needed

        best_spectrum = None
        highest_snr = -1
        chi_squares = []

        # Iterate through the spectra to find the one with the highest SNR
        for sp in spectra_list:
            snr = sp.snr  # Assuming SNR can be accessed directly; adjust if necessary
            if snr > highest_snr:
                highest_snr = snr
                best_spectrum = sp

        if best_spectrum is None:
            return None, None  # Handle case where no valid spectrum was found

        # Perform classification on the best spectrum
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            classification = splat.classifyByStandard(best_spectrum, fit_ranges=[0.87,2.39], verbose=True, all=True)
            
        # Normalize the best spectrum
        best_spectrum.normalize([1.27, 1.28])
        
        output = buffer.getvalue()
        matches = re.findall(r'Type (.*): statistic = (.*?),', output)
        statistics = [(type_, float(statistic)) for type_, statistic in matches]
        best_fit = min(statistics, key=lambda x: x[1])
        chi_squares.append(best_fit[1])

        if best_fit[0] in spectral_types:
            data.append((best_spectrum.name, best_fit[0], best_fit[1], highest_snr))
            spectral_type = classification[0]
            fluxes = [np.mean(best_spectrum.flux.value[np.where((best_spectrum.wave.value >= wl[0]) & (best_spectrum.wave.value <= wl[1]))[0]]) for wl in wav_003]
            flux_means.append([best_spectrum.name, spectral_type] + fluxes)
    
    except Exception as e:
        print(f"An error occurred with {name}. Error: {e}. Skipping this spectrum.")

    return data, flux_means


# +
import io
import contextlib
import re
import numpy as np
import splat

sp = splat.Spectrum(file='spex-prism_twa8b.fits',instrument='SPEX-PRISM',name='TWA 8B (M6.0gamma Std)')
splat.initiateStandards(vlg=True)
splat.STDS_VLG_SPEX['M6.0gamma'] = sp

def process_spectrum(name):
    # Initialize an empty list to hold the data
    data = []
    flux_means = []

    try:
        spectra_list = splat.getSpectrum(name=name)
        if not spectra_list:
            return None, None  # Return or handle empty spectra_list as needed

        best_spectrum = None
        highest_snr = -1
        chi_squares = []

        # Iterate through the spectra to find the one with the highest SNR
        for sp in spectra_list:
            snr = sp.snr  # Assuming SNR can be accessed directly; adjust if necessary
            if snr > highest_snr:
                highest_snr = snr
                best_spectrum = sp

        if best_spectrum is None:
            return None, None  # Handle case where no valid spectrum was found

        # Perform classification on the best spectrum
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            classification = splat.classifyByStandard(best_spectrum, fit_ranges=[0.87,2.39], verbose=True, all=True)
            
        # Normalize the best spectrum
        best_spectrum.normalize([1.27, 1.28])
        
        output = buffer.getvalue()
        matches = re.findall(r'Type (.*): statistic = (.*?),', output)
        statistics = [(type_, float(statistic)) for type_, statistic in matches]
        best_fit = min(statistics, key=lambda x: x[1])
        chi_squares.append(best_fit[1])

        if best_fit[0] in spectral_types:
            data.append((best_spectrum.name, best_fit[0], best_fit[1], highest_snr))
            spectral_type = classification[0]
            fluxes = [np.mean(best_spectrum.flux.value[np.where((best_spectrum.wave.value >= wl[0]) & (best_spectrum.wave.value <= wl[1]))[0]]) for wl in wav_002]
            flux_means.append([best_spectrum.name, spectral_type] + fluxes)
    
    except Exception as e:
        print(f"An error occurred with {name}. Error: {e}. Skipping this spectrum.")

    return data, flux_means
