# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:15:08 2021

@author: katierobinson
"""

# Time control
start_year = 2015
end_year = 2018
timestep = 3      # Should match timestep of input data. Units = hours.

# Input parameters
params_file = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Tuning/Tuning_Params_Initial10k.csv' # Path to csv file containing MF, aice, asnow.
sim = 666                                             # Corresponds to a row of the params file.
R2S = 1.0                                                  # Rain-to-snow threshold (degree C).
Glacier_ID = 'KRH'                                         # Identifier for the glacier. All inputs/outputs will also have the same ID.

# Input data
Model_functions = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel' # Path to folder where model functions script is.
Precip_inputs = 'D:/BiasCorrected_files/KRH'                                 # Path to folder where downscaled & bias-corrected data is.
Temp_inputs = 'D:/BiasCorrected_files/KRH'
Solar_inputs = 'D:/Solar_files/KRH'

# Input geometry
Easting_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/KRH_Xgrid.txt' # Paths to text files defining Easting/Northing coords of every model gridcell
Northing_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/KRH_Ygrid.txt'
Sfc_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/KRH_SfcType.txt'      # Path to text file where 1 = off-glacier, 0 = on-glacier, NaN = not in the domain.

# Debris parameters
debris_parameterization = 'Sub-debris melt scaling' # Debris treatment. Options are (1) None (no debris), (2) 'Boolean debris' (set a_ice = 0 in debris cells), or (3) 'Sub-debris melt scaling' (thickness dependent scaling)
debris_map = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/DebrisThickness/KRH_debrismap.txt' # Path to text file. If debris parameterization is (1) file should be same as Sfc; (2) debris cells are finite, else NaN, (3) cells w debris have value = debris thickness, else NaN
cleanice_melt = 2.0277           # Observed clean ice melt (m w.e).
peak_melt = 2.1717               # Observed peak melt (m w.e.).
peak_melt_thickness = 0.006      # Thickness at which peak melt occurs (m).
transition_thickness = 0.019     # Thickness at which melt equals clean ice melt (m).
b0 = 11.0349260206858            # b0 and k are site-specific params from Rounce et al. (2021).
k = 1.98717418666925

# Outputs
OUTPUT_PATH = 'D:/TuningOutputs/Tuning_Nov8_TESTING_ALL_OUTPUTMODES'                   # Path to directory where model outputs should be saved
SaveMBonly = True                        # If true, model saves mass balance only. If false, model saves Ice Melt, Snow Melt, Refreezing, Mass Balance