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

#parameterizations
#debris
debris = True  #turn debris on or off
debris_treatment = 'Boolean' # or 'Variable Thickness'
cleaniceM = 2.0277000000000003 #observed dirty ice melt #2.9200183038347 (DR) #m w.e.
peakM = 2.1717000000000004 #observed peak melt #6.62306028549649 (DR) #m w.e.
peakM_thickness = 0.006 # m (reference value)
peakthickness_uncertainty = 0.003

transition_thickness = 0.019 #m (reference value)
transitionthickness_uncertainty = 0.007

b0 = 11.0349260206858 #b0 and k are params fitted to KW from Rounce et al. (2021)
k = 1.98717418666925

if debris_treatment == 'Boolean':
    debris_thickness_map = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel\debmap_boolean.npy' #path to map 
elif debris_treatment == 'Variable Thickness': 
    debris_thickness_map = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel\debmap_variableh.npy' #path to map 
#transition_thickness = 0.04 #transition between melt-enhancing and melt-reducing debris, in meters

Rain_to_snow = 1   #set the rain to snow threshold (line ~132)
Refreezing = True
Temp_shift = False # do you want to change the entire temperature array up or down by a uniform constant? True = yes 
temp_shift_factor = 0 # what is the temperature shift. + is an increase in temp, - is a decrease
Bias_CorrectionT = True #are you using bias corrected temp files as input (True) or not (False)
Bias_CorrectionP = True #are you using bias corrected temp files as input (True) or not (False)

Tuning = False
param_total = 3 #how many parameter combinations to generate for the tuning process

#file_names
#Input_path = .../.../… #tell the model where to find the input files
params_filename = 'final_params_deb.csv' #tell the model which file contains the parameters we want to use for this run
#where are the downscaled/bias corrected inputs stored
#T_inputs = 'D:\\Katie\\Mass Balance Model\\MassBalanceModel_KatiesVersion\\Final_runs'
T_inputs = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\Kaskonly_R2S=1'
P_inputs = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\Kaskonly_R2S=1'
SR_inputs = 'F:\\Mass Balance Model\\Kaskonly_Downscaled_NoBC'
#File_glacier_in = ?? #need to change so that debris/non-debris case is more clear
Output_path = 'D:\TuningOutputs\Test_0628' #tell the model where to put the output netcdf files
File_sufix = ".nc"
