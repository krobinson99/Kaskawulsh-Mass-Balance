# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:15:08 2021

@author: katierobinson
"""

#time_control
time_step = 10800     #timestep - keep in units of seconds so that this is more flexible?
start_year = 2006
end_year = 2018
start_day = 1
#domain
glacier_id = 'kaskonly'       #points to a specific glacier outline --> dont change this: it s the same for the kaskawulsh inputs and catchment inputs
Considering_Catchment = False #are we simulating the whole catchment? (True) or just the Kask glacier (False)

#physics/dynamics
debris = True  #turn debris on or off
debris_treatment = 'Boolean' # or 'Variable Thickness'
debris_thickness_map = #path to map 
transition_thickness = 0.04 #transition between melt-enhancing and melt-reducing debris, in meters

Rain_to_snow = 1   #set the rain to snow threshold (line ~132)
Refreezing = True
Temp_shift = True # do you want to change the entire temperature array up or down by a uniform constant? True = yes 
temp_shift_factor = -1.7 # what is the temperature shift. + is an increase in temp, - is a decrease
Bias_CorrectionT = True #are you using bias corrected temp files as input (True) or not (False)
Bias_CorrectionP = True #are you using bias corrected temp files as input (True) or not (False)

Tuning = True
param_total = 1000 #how many parameter combinations to generate for the tuning process

#file_names
#Input_path = .../.../… #tell the model where to find the input files
params_filename = 'final_params_deb.csv' #tell the model which file contains the parameters we want to use for this run
#where are the downscaled/bias corrected inputs stored
#T_inputs = 'D:\\Katie\\Mass Balance Model\\MassBalanceModel_KatiesVersion\\Final_runs'
T_inputs = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\Kaskonly_R2S=1'
P_inputs = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\Kaskonly_R2S=1'
SR_inputs = 'F:\\Mass Balance Model\\Kaskonly_Downscaled_NoBC'
#File_glacier_in = ?? #need to change so that debris/non-debris case is more clear
Output_path = 'F:\\Mass Balance Model\\OUTPUTS\\Diagnostic\\Tshift=-1.7' #tell the model where to put the output netcdf files
File_sufix = ".nc"
