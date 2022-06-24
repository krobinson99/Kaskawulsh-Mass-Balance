# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 13:27:44 2022

tuning script revamped!! Now with spatially distributed ELA :) 

@author: agribbon
"""

# new features to include in the tuning script: 
# 1. tuning script should take outputs from the mass balance model itself, since I will be making lots of
#       edits to the model soon and don;t want to have to also change the tuning script
# WHAT are the tuning criteria:
    #a) geodetic mass balance from EITHER 1977-2007 or 2007-2018
    #b) spatially + temporally distruted ELAs
    
# procedure:
#    1. add option to MBM to run as TUNING or run as regular
#           - this will decide if the model uses specified input parameters OR randomly generates 1000
#           - if tuning is specified in MBM run: save ONLY mass balance outputs as .npy array not NETCDF: to save time!
#    2. run MBM with tuning = TRUE --> save outputs to a specific folder 
#    3. calculate geodetic mass balance for the calibration period for each of the 1000 simulations
#    4. only retain the sims (files saved by sim number) that pass the geodetic mass balance
#           - save txt file with list of sim numbers that should be thrown out
#           - use bash rm * commands to delete the non-passing files (echo $listval, rm file_{$listval}) (something like that)
#    5. for each time frame that an ELA is dilineated for --> compare each of the passing sims to it
#           - Q: do I need to add up to MB to that point? OR just look at the single timestep
#                - I think I'd need to add everything up to that point since MB is (change in mass)/dt
#                - add up the MB for that hydrological balance year? iestarting Oct 1 (?)
#                - is snowline the same as equilibrium line in the model? YES IF you add the MB timesteps
#                   starting on Oct 1 (beginning of hydrological year.) because then summed MB will only be < 0 where no snowpack exists (ie total abl exceeded total acc)
#                - Q for Gwenn: confirm that I should sum MB starting at hydrological year and not starting at model beginning. 
#    6. add non-passing sims to the rm list before running bash script
#    7. for all the sims that PASS both stages: save the geodetic MB and metrics about how close the ELA was to see how strict the tuning criteria are being
#           - during the testing stages, save the MB and ELAs from ALL 1000 runs, we can use them to decide how strict to make the criteria! (ie sort all results into a bar chart)
#    8. save the params that correspond the to passing sims in a new list: will be used to run the model! 
#    
#
#
    
    
