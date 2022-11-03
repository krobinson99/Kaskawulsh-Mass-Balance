# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 10:37:34 2022

Generate melt curves (melt vs debris thickenss) that can be converted into 
melt enhancement factors for the whole domain (using a function).
The melt curve should be generated using values of peak-melt thickness and
transition thickness from a guassian distribution.
(see \Kaskawulsh-Mass-Balance/FieldWork2022/ablationstakes.py)

@author: katierobinson
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import os
import sys
from scipy.interpolate import splrep
from scipy.interpolate import splev
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import regridXY_something



#1. generate the original MF curve with the original params (linear to 2 cm, the eq(1) from 2 cm onwards)

#params
peakM = 6.62306028549649 # from hdopt_params.csv (melt_mwea_2cm)
b0 = 11.0349260206858 #from hdopt_params.csv (b0)
k = 1.98717418666925 # from hdopt_params.csv
cleaniceM = 2.9200183038347
peakthickness = 0.02 #m

debrisparams ='F:/Mass Balance Model/Kaskawulsh-Mass-Balance/FieldWork2022/debrisparameters.csv'

df = pd.read_csv(debrisparams)
    #df['Tributary'].astype('Int64')
arr = np.array(df)

refcase = arr[1]
peakthickness_ref = (refcase[1])/100
peakthickness_ref_unc = (refcase[2])/100
transition_ref = (refcase[3])/100
transition_ref_unc = (refcase[4])/100 #converted to units = m

refcase_data ='F:/Mass Balance Model/Kaskawulsh-Mass-Balance/FieldWork2022/pddaverage_debris.csv'

df = pd.read_csv(refcase_data)
arr2 = np.array(df)
ref_debthickness = (arr2[:,0][:-1])/100
ref_debthickness[1] = 0.000001
ref_melt = arr2[:,1][:-1]
dirtyicemelt_obs = ref_melt[1]
peakmelt_obs = ref_melt[2]

def generate_meltfactors(peakM,peakthickness,h_transition,b0 = 11.0349260206858,k = 1.98717418666925):
    '''
    h_transition = transition thickness (default is 13cm)
    In DR's original curve there was no h_transition specified but it just endedup being 13cm
    to get the original curve, use setcleaniceM = 2.9200183038347 (from hdopt_params.csv)
    '''
    cleaniceM = b0/(1+(k*b0*h_transition))
    
    debristhickness = np.linspace(0,2,10000)
    melt = []
    for h in debristhickness:
        if h <= peakthickness:
            M = (((peakM-cleaniceM)/(peakthickness))*h) + cleaniceM #melt is linear between the clean ice melt and the peak melt
            melt.append(M)
        else:
            M = b0/(1+(k*b0*h)) #after the peak melt the curve follows eq(1) from Rounce et al. 2021
            melt.append(M)
        # melt should be capped at the "peak melt"
    for i in range(0,len(melt)):
        if melt[i] > peakM:
            melt[i] = peakM
        else:
            pass
    
    
    meltfactors = np.array(melt)/cleaniceM
    
    return debristhickness, melt, meltfactors, cleaniceM
        
debrisDR, meltDR, meltfactorsDR, cleaniceDR = generate_meltfactors(6.62306028549649,0.02,0.13,11.0349260206858,1.98717418666925)
        
plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plt.plot(debrisDR,meltDR)
plt.margins(x=0)
#plt.xlim(0,0.05)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt (m w.e.)')
plt.title('Original melt curve (DR)')
plt.subplot(1,2,2)
plt.plot(debrisDR,meltfactorsDR)
plt.margins(x=0)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt factor')
plt.title('Original melt curve (DR)')
plt.tight_layout()
#plt.savefig('Meltcurve_DR_original.png',bbox_inches = 'tight')

#2. Now try it with the same parameters EXCEPT change the peak melt THICKNESS and transition thickness to the value determined in ablationstakes.py

debris_refvals, melt_refvals, meltfacts_refvals, cleanice_refvals = generate_meltfactors(6.62306028549649,peakthickness_ref,transition_ref,11.0349260206858,1.98717418666925)


plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plt.plot(debris_refvals,melt_refvals)
plt.margins(x=0)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt (m w.e.)')
plt.title('DR equation with field vals: \n Thickness at peak melt = 0.6 cm \n Transition thickness = 1.9 cm')
plt.subplot(1,2,2)
plt.plot(debris_refvals,meltfacts_refvals)
plt.margins(x=0)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt factor')
plt.title('Reference Values \n Thickness at peak melt = 0.6 cm \n Transition thickness = 1.9 cm')
plt.tight_layout()
#plt.savefig('Meltcurve_referenceparams.png',bbox_inches = 'tight')


#Now plot what the range of debris melt factor curves would like like for 2 scenarios:
    # a) we test the range of param values using the DR curve.
    # b) we test the range of param values using our own curve (from field data)
    
#modify this func. from ablationstakes.py so that there's always a knot at the specified peak thickness 
def splrepcurve(peakmelt_thickness,transition_thickness,varypeak,k=3,debris=ref_debthickness,melt=ref_melt):
    if varypeak == True:
        knots = np.array([peakmelt_thickness]) #,transition_thickness])
    else:
        knots = np.array([transition_thickness])
    tck = splrep(debris, melt, k=k,s=0.01,t=knots)
    xnew = np.linspace(0, np.max(debris)+0.01, 100)
    ynew = splev(xnew,tck)
    
    return xnew, ynew, tck

#2 rows, 2 cols: top row = using DR equation: changing peak melt thickness vs changing transition thickness
#bottom row: using my own curves from field obs. Q: how to make them changeable based on the given peak melt and trans thickness.

plt.figure(figsize=(12,12))
plt.subplot(2,2,1)
for i in np.linspace(peakthickness_ref-peakthickness_ref_unc,peakthickness_ref+peakthickness_ref_unc,15):
    debris_refvals, melt_refvals, meltfacts_refvals, cleanice_refvals = generate_meltfactors(peakM,i,transition_ref,b0,k)
    plt.plot(debris_refvals,meltfacts_refvals)
plt.margins(x=0)
plt.xlim(0,1)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt factors')
plt.title('DR Equation with varying thickness at peak melt')

plt.subplot(2,2,2)
for i in np.linspace(transition_ref-transition_ref_unc,transition_ref+transition_ref_unc,15):
    debris_refvals, melt_refvals, meltfacts_refvals, cleanice_refvals = generate_meltfactors(peakM,peakthickness_ref,i,b0,k)
    plt.plot(debris_refvals,meltfacts_refvals)
plt.margins(x=0)
plt.xlim(0,1)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt factors')
plt.title('DR Equation with varying transition thickness')
#plt.savefig('DRmeltcurves_withrefvals.png',bbox_inches = 'tight')

plt.subplot(2,2,3)
#vary peak thickness but keep transition thickness the same 
cim = ref_melt[0]
for i in np.linspace(peakthickness_ref-peakthickness_ref_unc,peakthickness_ref+peakthickness_ref_unc,15):
    debris, melt, tck = splrepcurve(i,transition_ref,varypeak=True)
    meltfacts = melt/cim
    plt.plot(debris,meltfacts)
plt.margins(x=0)
plt.xlim(0.003,0.009)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt factors')
plt.title('Field data (spline) with varying thickness at peak melt')

plt.subplot(2,2,4)
#vary peak thickness but keep transition thickness the same 
cim = ref_melt[0]
for i in np.linspace(transition_ref-transition_ref_unc,transition_ref+transition_ref_unc,15):
    debris, melt, tck = splrepcurve(peakthickness_ref,i,varypeak=False)
    meltfacts = melt/cim
    plt.plot(debris,meltfacts)
plt.margins(x=0)
plt.xlim(0.012,0.026)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt factors')
plt.title('Field data (spline) with varying transition thickness')
#plt.savefig('Allpossiblemeltcurves.png',bbox_inches = 'tight')

###############################################################################
# TRY DEFINING A CURVE THAT IS PARTIALLY LINEAR (CLEAN ICE TO PEAK TO TRANSITION)
###### AND PARTIALLY EQ (1) FROM ROUNCE ET AL. (> TRANSITION THICKNESS) #######
###############################################################################

def generate_hybrid_curve(cleaniceM,peakM,peakM_thickness,transition_thickness,b0 = 11.0349260206858,k = 1.98717418666925):
    '''
    generate a curve describing the relationship between melt and debris thickness
    that is a combo of the field data and eq(1) from Rounce et al. 2021
    '''
    debristhickness_list = [0,peakM_thickness,transition_thickness]
    melt_list = [cleaniceM,peakM,cleaniceM]
    
    #figure out where to start the second part of the curve
    h_curvestart = (b0-cleaniceM)/(cleaniceM*k*b0)
    
    for h in np.linspace(h_curvestart,3,100):
        M = b0/(1+(k*b0*h))
        melt_list.append(M)
        debristhickness_list.append(h)
        
    #close the discontinuity by shifting the second piece of the curve
     # so that it starts at the transition thickness and cleanice melt
     
    #1. calculate the difference between the clean ice melt and the melt at
        # the transition thickness
        
    debristhickness = np.array(debristhickness_list)
    melt = np.array(melt_list)
    
    #verticalshift = cleaniceM - melt[3]
    #melt[3:] = melt[3:] + verticalshift
    
    horizontal_shift = transition_thickness - h_curvestart
    debristhickness[3:] = debristhickness[3:] + horizontal_shift
    
    meltfactors = melt/cleaniceM
    
    return debristhickness, melt, meltfactors

#shift horizontally instead of vertically?
    #vertical shift (take the curve at the transition thickness and move it down) leads to negative melt values
    # horizontal shift (take the curv at the clean ice melt and move it left)

deb, m, mf = generate_hybrid_curve(cleaniceM,peakM,peakthickness_ref,transition_ref)    
deb2, m2, mf2 = generate_hybrid_curve(dirtyicemelt_obs,peakmelt_obs,peakthickness_ref,transition_ref)    

plt.figure(figsize=(12,4))
plt.subplot(1,2,1)
plt.plot(deb,mf)
plt.xlim(-0.05,2)
plt.ylim(0,2.3)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt factors')
plt.margins(x=1)
plt.title('Clean ice melt & peak melt = Rounce model values',fontsize=9)

plt.subplot(1,2,2)
plt.plot(deb2,mf2)
plt.xlim(-0.05,2)
plt.ylim(0,1.4)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt factors')
plt.margins(x=1)
plt.title('Clean ice melt & peak melt = Observed (field work) values',fontsize=9)
#plt.savefig('Hybridmeltcurves.png',bbox_inches = 'tight')

plt.figure(figsize=(5,5))
plt.plot(debrisDR,meltDR/meltDR[0],linewidth=3,color='k')
plt.plot(deb2,mf2,linewidth=3,color='blue')
plt.xlim(-0.05,2)
plt.xlabel('Debris Thickness (m)',fontsize=14)
plt.ylabel('Melt Factor',fontsize=14)
plt.margins(x=1)
plt.legend(['Rounce et al. (2021) melt curve','Melt curve using site-specific values'],fontsize=12)
plt.tight_layout()
#plt.savefig('meltcurves_comparison.png',bbox_inches = 'tight')


#plot 4 hybrid curves: top row: using DR values for peakM and cleaniceM, bottom row: using field values for peakM and cleaniceM
# left column: changing peak melt thickness, right column: changing transition thickness
plt.figure(figsize=(10,10))
plt.subplot(2,2,1)
legend = []
for i in np.linspace(peakthickness_ref-peakthickness_ref_unc,peakthickness_ref+peakthickness_ref_unc,8):
    deb, m, mf = generate_hybrid_curve(cleaniceM,peakM,i,transition_ref)   
    plt.plot(deb,mf)
    legend.append(str(np.round(i*100,2)) + ' cm')
plt.legend(legend)
plt.margins(x=0)
plt.xlim(0,0.05)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt factors')
plt.title('Varying thickness at peak melt \n Clean ice melt & peak melt = Rounce model values',fontsize=9)

plt.subplot(2,2,2)
legend = []
for i in np.linspace(transition_ref-transition_ref_unc,transition_ref+transition_ref_unc,8):
    deb, m, mf = generate_hybrid_curve(cleaniceM,peakM,peakthickness_ref,i)   
    plt.plot(deb,mf)
    legend.append(str(np.round(i*100,2)) + ' cm')
plt.legend(legend)
plt.margins(x=0)
plt.xlim(0,0.05)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt factors')
plt.title('Varying transition thickness \n Clean ice melt & peak melt = Rounce model values',fontsize=9)

plt.subplot(2,2,3)
legend = []
for i in np.linspace(peakthickness_ref-peakthickness_ref_unc,peakthickness_ref+peakthickness_ref_unc,8):
    deb, m, mf = generate_hybrid_curve(dirtyicemelt_obs,peakmelt_obs,i,transition_ref)   
    plt.plot(deb,mf)
    legend.append(str(np.round(i*100,2)) + ' cm')
plt.legend(legend)
plt.margins(x=0)
plt.xlim(0,0.05)
plt.ylim(0.8,1.2)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt factors')
plt.title('Varying thickness at peak melt \n Clean ice melt & peak melt = Observed (field work) values',fontsize=9)

plt.subplot(2,2,4)
legend = []
for i in np.linspace(transition_ref-transition_ref_unc,transition_ref+transition_ref_unc,8):
    deb, m, mf = generate_hybrid_curve(dirtyicemelt_obs,peakmelt_obs,peakthickness_ref,i)   
    plt.plot(deb,mf)
    legend.append(str(np.round(i*100,2)) + ' cm')
plt.legend(legend)
plt.margins(x=0)
plt.xlim(0,0.05)
plt.ylim(0.8,1.2)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt factors')
plt.title('Varying transition thickness \n Clean ice melt & peak melt = Observed (field work) values',fontsize=9)
plt.tight_layout()
#plt.savefig('Hybridmeltcurves_ranges.png',bbox_inches = 'tight')


###############################################################################
# WRITE FUNCTIONS TO GENERATE DEBRIS PARAMS AND CALCULATE THE MELT FACTORS MASK
###############################################################################
# (work on the function here but eventually add to the model functions script)
File_glacier_in = os.path.join('F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel','Kaskonly_deb.txt')
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
Ix = glacier[:,3] 
Iy = glacier[:,4] 
Ih = glacier[:,2] 
debris_array = glacier[:,6]

#-------Turn vectors into 3D gridded inputs--------------------

Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
nanlocs3 = np.where(np.isnan(Zgrid))

#Setup debris mask for use in radiation parameters
debris_grid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, debris_array)
debris_m = np.zeros(debris_grid.shape)
debris_m[np.where(debris_grid > 100)] = 1
debris_m[np.where(debris_grid <= 100)] = np.nan
debris_m[nanlocs3] = np.nan

plt.figure(figsize=(8,5))
#plt.contourf(debristhickness_array, cmap = 'PuOr', levels = np.round(np.linspace(0,3,10),1))
plt.contourf(np.flipud(debris_m), cmap = 'PuOr', levels = np.linspace(0,1,3))
#plt.contourf((debris_m[:180,100:]), cmap = 'PuOr', levels = np.round(np.linspace(0,1,3),1))
legend = plt.colorbar()
legend.ax.set_ylabel('Debris True / False', rotation=270,fontsize=14,labelpad=20)
plt.tight_layout()
plt.title('Kaskawulsh Boolean Debris Map',fontsize=14)
#plt.savefig(os.path.join(Path2files,'KW_booleandebris.png'),bbox_inches = 'tight')

    
def generate_meltfactors(debrismask,cleaniceM,peakM,peakM_thickness,transition_thickness,b0 = 11.0349260206858,k = 1.98717418666925):
    '''
    the final version of this function should always be in the Model Functions script
    #function should take in:
        # debris mask,
        # cleanicemelt and peak melt value
        # transition thickness and thickness at peak melt 
        # Rounce et al. model params for the kaskawulsh
    #function should return:
        # 2D melt factors map with same shape as input debris mask. should be 1 in non debris areas
    so that it doesnt change the melt of cleanice when multiplied with the melt array
    '''
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    #make an empty array of same shape as debris mask
    #loop through each gridcell in debris mask
    #if nan (ie no debris) place a 1 in the meltfactors array
    #sort the debris cells by thickness and calculate corresponding MF, then place in melt array
    
    meltfactors = np.empty(debrismask.shape)
    for x in range(len(debrismask)):
        for y in range(len(debrismask[0])):
            #print(debrismask[x,y])
            if np.isnan(debrismask[x,y]):
                meltfactors[x,y] = 1
            elif debrismask[x,y] <= peakM_thickness:
                Melt = (((peakM - cleaniceM)/peakM_thickness)*debrismask[x,y]) + cleaniceM
                Meltfact = Melt/cleaniceM
                meltfactors[x,y] = Meltfact
                #print(str(debrismask[x,y]) + ' ,,, ' + str(Meltfact))
                #calculate the melt factor on a linear line between clean ice melt and peak melt:
            elif (debrismask[x,y] > peakM_thickness) and (debrismask[x,y] <= transition_thickness):
                #print(debrismask[x,y])
                yintercept = peakM - (((cleaniceM-peakM)/(transition_thickness-peakM_thickness))*peakM_thickness)
                Melt = (((cleaniceM-peakM)/(transition_thickness-peakM_thickness))*debrismask[x,y]) + yintercept
                Meltfact = Melt/cleaniceM
                meltfactors[x,y] = Meltfact
            elif debrismask[x,y] > transition_thickness:
                #determine how much this part of the curve needs to be shifted over to meet the linear parts
                h_curvestart = (b0-cleaniceM)/(cleaniceM*k*b0)
                horizontal_shift = transition_thickness - h_curvestart
                Melt = b0/(1+(k*b0*(debrismask[x,y]-horizontal_shift))) #calculate the melt on the original curve(take away the shift) so that you're actually getting the melt at the desired thickness
                Meltfact = Melt/cleaniceM
                meltfactors[x,y] = Meltfact
            else:
                print('debris thickness not accounted for in if statements')

    return meltfactors 
    
    
    
test = generate_meltfactors(debris_m,cleaniceM,peakM,peakthickness_ref,transition_ref)
    
# test it with the actual rounce et al. debris map
debmap = np.load('F:\Mass Balance Model\Kaskawulsh-Mass-Balance\DebrisThickness\debmap_variableh.npy')

meltfacts_test = generate_meltfactors(debmap,dirtyicemelt_obs,peakmelt_obs,peakthickness_ref,transition_ref)
meltfacts_test[nanlocs3] = np.nan
#meltfacts_test[np.where(meltfacts_test == 0)] = np.nan

plt.figure(figsize=(8,5))
#plt.contourf(debristhickness_array, cmap = 'PuOr', levels = np.round(np.linspace(0,3,10),1))
#norm = mcolors.TwoSlopeNorm(vmin=0., vcenter=1.0, vmax=1.2)
plt.contourf(np.flipud(meltfacts_test[:180,100:]),  cmap = 'coolwarm', levels = np.linspace(0,2,41))
#plt.contourf((debris_m[:180,100:]), cmap = 'PuOr', levels = np.round(np.linspace(0,1,3),1))
legend = plt.colorbar()
legend.ax.set_ylabel('Melt Factor', rotation=270,fontsize=14,labelpad=20)
xticks = [0,50,100,150,200]
xlabels = []
for x in xticks:
    km = int(x*0.2)
    xlabels.append(km)
plt.xticks(ticks=xticks,labels=xlabels)
yticks = [0,50,100,150]
ylabels = []
for y in yticks:
    km = int(y*0.2)
    ylabels.append(km)
plt.yticks(ticks=yticks,labels=ylabels)
#roundedX = Xgrid+18
#roundedY = Ygrid+15
#plt.xticks(ticks=[0,50,100,150,200],labels=[int(roundedX[0][100:][0]),int(roundedX[0][100:][50]),int(roundedX[0][100:][100]),int(roundedX[0][100:][150]),int(roundedX[0][100:][200])],fontsize=12,rotation=20)
#plt.yticks(ticks=[0,40,80,120,160],labels=[int(roundedY[:,0][:180][0]),int(roundedY[:,0][:180][40]),int(roundedY[:,0][:180][80]),int(roundedY[:,0][:180][120]),int(roundedY[:,0][:180][160])],fontsize=12,rotation=20)
#plt.xlabel('Easting',fontsize=12)
#plt.ylabel('Northing',fontsize=12)
plt.xlabel('distance (km)',fontsize=12)
plt.ylabel('distance (km)',fontsize=12)
plt.title('Site-Specific Sub-Debris Melt Factors',fontsize=14)
plt.tight_layout()
plt.savefig('New_Meltfactors.png',bbox_inches = 'tight')

deb = []
mf = []
for i in range(0,len(meltfacts_test)):
    for j in range(0,len(meltfacts_test[0])):
        if np.isnan(debmap[i,j]):
            pass
        else:
            deb.append(debmap[i,j])
            mf.append(meltfacts_test[i,j])
            
plt.figure(figsize=(6,5))
plt.scatter(deb,mf,s=2)
plt.xlim(-0.01,0.4)        

oldKWmeltfacts = generate_meltfactors(debmap,cleaniceM,peakM,0.02,0.13)
oldKWmeltfacts[nanlocs3] = np.nan
#oldKWmeltfacts[np.where(oldKWmeltfacts == 0)] = np.nan

plt.figure(figsize=(8,5))
#plt.contourf(debristhickness_array, cmap = 'PuOr', levels = np.round(np.linspace(0,3,10),1))
#norm = mcolors.TwoSlopeNorm(vmin=0., vcenter=1.0, vmax=1.2)
plt.contourf(np.flipud(oldKWmeltfacts[:180,100:]),  cmap = 'coolwarm',levels = np.linspace(0,2,41)) #levels = np.round(np.linspace(0,2,9),1))
#plt.contourf((debris_m[:180,100:]), cmap = 'PuOr', levels = np.round(np.linspace(0,1,3),1))
legend = plt.colorbar()
legend.ax.set_ylabel('Melt Factor', rotation=270,fontsize=14,labelpad=20)
xticks = [0,50,100,150,200]
xlabels = []
for x in xticks:
    km = int(x*0.2)
    xlabels.append(km)
plt.xticks(ticks=xticks,labels=xlabels)
yticks = [0,50,100,150]
ylabels = []
for y in yticks:
    km = int(y*0.2)
    ylabels.append(km)
plt.yticks(ticks=yticks,labels=ylabels)
plt.xlabel('distance (km)',fontsize=12)
plt.ylabel('distance (km)',fontsize=12)
plt.title('Rounce et al. (2021) Sub-Debris Melt Factors',fontsize=14)
plt.tight_layout()
plt.savefig('Old_Meltfactors.png',bbox_inches = 'tight')

#calculate how much of the debris covered area is reducing melt
gridcellarea = 0.2*0.2 # km
alldebcells = len((np.where(debmap >= 0))[0])
alldebarea = alldebcells * gridcellarea

totalglaciercells = len((np.where(oldKWmeltfacts >= 0))[0])
totalglacierarea = totalglaciercells * gridcellarea


thindebcells = len((np.where(debmap <= transition_ref))[0])
thindebarea = thindebcells * gridcellarea

enhancedmeltarea = thindebarea/alldebarea*100

enhancedmeltarea_gl = thindebarea/totalglacierarea*100
