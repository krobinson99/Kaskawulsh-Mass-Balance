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
import pandas as pd
from scipy.interpolate import splrep
from scipy.interpolate import splev

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
    
    horizontal_shift = transition_thickness - debristhickness[3]
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


#plot 4 hybrid curves: top row: using DR values for peakM and cleaniceM, bottom row: using field values for peakM and cleaniceM
# left column: changing peak melt thickness, right column: changing transition thickness
plt.figure(figsize=(10,10))
plt.subplot(2,2,1)
legend = []
for i in np.linspace(peakthickness_ref-peakthickness_ref_unc,peakthickness_ref+peakthickness_ref_unc,8):
    print(i*100)
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

    
    
    
    
    
    
    
    
    
