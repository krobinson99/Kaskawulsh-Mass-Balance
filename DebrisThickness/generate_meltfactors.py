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
plt.savefig('Meltcurve_DR_original.png',bbox_inches = 'tight')

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
plt.savefig('Meltcurve_referenceparams.png',bbox_inches = 'tight')


#Now plot what the range of debris melt factor curves would like like for 2 scenarios:
    # a) we test the range of param values using the DR curve.
    # b) we test the range of param values using our own curve (from field data)
    

#2 rows, 2 cols: top row = using DR equation: changing peak melt thickness vs changing transition thickness
#bottom row: using my own curves from field obs. Q: how to make them changeable based on the given peak melt and trans thickness.

plt.figure(figsize=(12,12))
plt.subplot(2,2,1)
for i in np.linspace(peakthickness_ref-peakthickness_ref_unc,peakthickness_ref+peakthickness_ref_unc,15):
    debris_refvals, melt_refvals, meltfacts_refvals, cleanice_refvals = generate_meltfactors(peakM,i,transition_ref,b0,k)
    plt.plot(debris_refvals,meltfacts_refvals)
plt.margins(x=0)
plt.xlim(0,2)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt (m w.e.)')
plt.title('DR Equation with varying thickness at peak melt')

plt.subplot(2,2,2)
for i in np.linspace(transition_ref-transition_ref_unc,transition_ref+transition_ref_unc,15):
    debris_refvals, melt_refvals, meltfacts_refvals, cleanice_refvals = generate_meltfactors(peakM,peakthickness_ref,i,b0,k)
    plt.plot(debris_refvals,meltfacts_refvals)
plt.margins(x=0)
plt.xlim(0,2)
plt.xlabel('debris thickness (m)')
plt.ylabel('melt (m w.e.)')
plt.title('DR Equation with varying transition thickness')

#modify this func. from ablationstakes.py so that there's always a knot at the specified peak thickness 
def slprepcurve(debris,k,melt=height_change*(p_ice/1000)):
    knots = np.array([peakmelt_thickness,transition thickness])
    tck = splrep(debris[cluster==1], melt[cluster==1],k=k,s=0.01,t=knots)
    xnew = np.linspace(0, np.max(debris[cluster==1])+0.01, 100)
    ynew = splev(xnew,tck)
    
    return xnew, ynew, tck



