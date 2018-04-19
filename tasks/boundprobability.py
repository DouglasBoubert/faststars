""" Calculate the probability that each star is bound
"""
import re

import numpy as np
from astropy.coordinates import SkyCoord as coord
import astropy.units as un
from math import sin, cos, sqrt
from galpy.potential import MWPotential2014, vesc
from scipy.stats import beta

from astrocats.catalog.photometry import PHOTOMETRY
from astrocats.catalog.utils import is_number, pbar, single_spaces, uniq_cdl
from astrocats.catalog.correlation import CORRELATION
from astrocats.catalog.quantity import QUANTITY
from ..faststars import FASTSTARS
from ..utils import name_clean, kinematic_sampler, solarmotion_sampler
mas_to_deg = 1e-3/3600.0

##### Coordinate transformation   
epoch=2000.0
k=4.74057
if epoch == 2000.0:
    cube= 122.932/180.*np.pi
    dec_ngp= 27.12825/180.*np.pi
    ra_ngp= 192.85948/180.*np.pi
elif epoch == 1950.0:
    cube= 123./180.*np.pi
    dec_ngp= 27.4/180.*np.pi
    ra_ngp= 192.25/180.*np.pi
    
T_1 = np.matrix([ [ np.cos(cube) , np.sin(cube) , 0. ], [ np.sin(cube) , -np.cos(cube) , 0. ], [ 0. , 0., 1.] ])
T_2 = np.matrix([ [ -np.sin(dec_ngp) , 0. , np.cos(dec_ngp) ], [ 0. , -1., 0.], [ np.cos(dec_ngp), 0., np.sin(dec_ngp) ] ])
T_3 = np.matrix([ [ np.cos(ra_ngp) , np.sin(ra_ngp) , 0. ], [ np.sin(ra_ngp) , -np.cos(ra_ngp) , 0. ], [ 0. , 0., 1.] ])
T=T_1*T_2*T_3

errorifmissing = 0.5
confidenceintervalcoverage = 0.6827
    
def best_parameter(PARAMETER):
    nPARAMETER = len(PARAMETER)
    errPARAMETER = 9999.9*np.ones(nPARAMETER)
    for i in range(nPARAMETER):
        if QUANTITY.DERIVED == True:
            continue
        elif QUANTITY.E_VALUE in PARAMETER[i]:
            errPARAMETER[i] = float(PARAMETER[i][QUANTITY.E_VALUE])
        elif (QUANTITY.E_LOWER_VALUE in PARAMETER[i] and QUANTITY.E_UPPER_VALUE in PARAMETER[i]):
            errPARAMETER[i] = max(float(PARAMETER[i][QUANTITY.E_LOWER_VALUE]),float(PARAMETER[i][QUANTITY.E_UPPER_VALUE]))
    return np.argmin(errPARAMETER)
    

def do_boundprobability(catalog):
    
    task_str = catalog.get_current_task_str()
    keys = list(catalog.entries.keys())

    for oname in pbar(keys, task_str):
        # Some events may be merged in cleanup process, skip them if
        # non-existent.
        try:
            name = catalog.add_entry(oname)
        except Exception:
            catalog.log.warning(
                '"{}" was not found, suggests merge occurred in cleanup '
                'process.'.format(oname))
            continue

        if (FASTSTARS.RA not in catalog.entries[name] or
                FASTSTARS.DEC not in catalog.entries[name] or (FASTSTARS.LUM_DIST not in catalog.entries[name] and FASTSTARS.PARALLAX not in catalog.entries[name])):
            # Don't need to worry about not having any velocities. Can't be in the catalogue if it doesn't have velocities!
            continue
        else:
            # Obtain highest precision position
            best_position_i = best_parameter(catalog.entries[name][FASTSTARS.RA])
            Mradec=str(catalog.entries[name][FASTSTARS.RA][best_position_i]['value'])+str(catalog.entries[name][FASTSTARS.DEC][best_position_i]['value'])
            c=coord(Mradec,unit=(un.hourangle, un.deg),frame='icrs')
            
            # Obtain galactic coordinates
            ra = c.ra.rad
            dec = c.dec.rad
            l = c.galactic.l.rad
            b = c.galactic.b.rad
            
            # Storage for means, errors and correlation
            kine_n = 4
            kine_values = np.zeros(kine_n)
            kine_errors = 9999.9*np.ones(kine_n)
            kine_corr = np.diag(np.ones(kine_n))
            kine_L = 1.35
            kine_k = 3.0
            
            # Flags for having proper motions and radial velocities
            havepropermotions = False
            havevelocities = False

            # Do we have a parallax?
            if FASTSTARS.PARALLAX in catalog.entries[name]:
                # Yes. 
                kine_values[0] = float(catalog.entries[name][FASTSTARS.PARALLAX][0]['value'])
                kine_errors[0] = float(catalog.entries[name][FASTSTARS.PARALLAX][0]['e_value'])
                
                # Is the parallax correlated with the proper motions?
                if QUANTITY.CORRELATIONS in catalog.entries[name][FASTSTARS.PARALLAX]:
                    # Yes. Loop through correlations.
                    parallax_corr = catalog.entries[name][FASTSTARS.PARALLAX][0]['correlations']
                    n_corr = len(parallax_corr)
                    for i in range(n_corr):
                        if parallax_corr[i]['quantity'] == 'propermotionra':
                            kine_corr[0,1] = kine_corr[1,0] = float(parallax_corr[i]['value'])
                        if parallax_corr[i]['quantity'] == 'propermotiondec':
                            kine_corr[0,2] = kine_corr[2,0] = float(parallax_corr[i]['value'])
            
            # Do we have proper motions?
            if (FASTSTARS.PROPER_MOTION_RA in catalog.entries[name] and FASTSTARS.PROPER_MOTION_DEC in catalog.entries[name]):
                # Yes.
                havepropermotions = True
                kine_values[1] = float(catalog.entries[name][FASTSTARS.PROPER_MOTION_RA][0]['value'])
                kine_values[2] = float(catalog.entries[name][FASTSTARS.PROPER_MOTION_DEC][0]['value'])
                
                # Do the proper motions have errors?
                if (QUANTITY.E_VALUE in catalog.entries[name][FASTSTARS.PROPER_MOTION_RA][0] and QUANTITY.E_VALUE in catalog.entries[name][FASTSTARS.PROPER_MOTION_DEC][0]):
                    kine_errors[1] = float(catalog.entries[name][FASTSTARS.PROPER_MOTION_RA][0]['e_value'])
                    kine_errors[2] = float(catalog.entries[name][FASTSTARS.PROPER_MOTION_DEC][0]['e_value'])
                else:
                    kine_errors[1] = errorifmissing * kine_values[1]
                    kine_errors[2] = errorifmissing * kine_values[2]
                  
                # Are the proper motions correlated?
                if (QUANTITY.CORRELATIONS in catalog.entries[name][FASTSTARS.PROPER_MOTION_RA][0] and QUANTITY.CORRELATIONS in catalog.entries[name][FASTSTARS.PROPER_MOTION_DEC][0]):
                    # Yes. Loop through correlations.
                    propermotionra_corr = catalog.entries[name][FASTSTARS.PROPER_MOTION_RA][0]['correlations']
                    n_corr = len(propermotionra_corr)
                    for i in range(n_corr):
                        if propermotionra_corr[i]['quantity'] == 'propermotiondec':
                            kine_corr[1,2] = kine_corr[2,1] = float(propermotionra_corr[i]['value'])
                        
            # Do we have a radial velocity?
            if FASTSTARS.VELOCITY in catalog.entries[name]:
                # Yes.
                havevelocities = True
                
                # What is the best radial velocity?
                best_velocity_i = best_parameter(catalog.entries[name][FASTSTARS.VELOCITY])
                kine_values[3] = float(catalog.entries[name][FASTSTARS.VELOCITY][best_velocity_i]['value'])
                
                # Does the radial velocity have an error?
                if QUANTITY.E_VALUE in catalog.entries[name][FASTSTARS.VELOCITY][best_velocity_i]:
                    kine_errors[3] = float(catalog.entries[name][FASTSTARS.VELOCITY][best_velocity_i]['e_value'])
                else:
                    kine_errors[3] = errorifmissing * kine_values[3]
            
            # Do we have a previous photometric distance?
            if FASTSTARS.LUM_DIST in catalog.entries[name]:
                # Yes.
                
                # What is the best photometric distance?
                best_lumdist_i = best_parameter(catalog.entries[name][FASTSTARS.LUM_DIST])
                kine_mu = float(catalog.entries[name][FASTSTARS.LUM_DIST][best_lumdist_i]['value'])
                
                # Does the photometric distance have an error?
                if QUANTITY.E_VALUE in catalog.entries[name][FASTSTARS.LUM_DIST][best_lumdist_i]:
                    kine_sigma = float(catalog.entries[name][FASTSTARS.LUM_DIST][best_lumdist_i]['e_value'])
                else:
                    kine_sigma = errorifmissing * kine_mu
                    
                kine_L = kine_sigma**2/kine_mu
                kine_k = kine_mu/kine_L
                
            # Sample
            kine_samples = kinematic_sampler(name,kine_values,kine_errors,kine_corr,KINE_L=kine_L,KINE_K=kine_k)
            kine_samples_corrected = np.copy(kine_samples)
            kine_samples_corrected[:,1] *= k*kine_samples[:,0]
            kine_samples_corrected[:,2] *= k*kine_samples[:,0]
            kine_samples_n = kine_samples.shape[0]

            # Calculate solar reflex correction
            cosl = cos(l)
            sinl = sin(l)
            cosb = cos(b)
            sinb = sin(b)
            cosdec = cos(dec)
            sindec = sin(dec)
            cosra = cos(ra)
            sinra = sin(ra)
            A = np.array([[cosra*cosdec, -sinra, -cosra*sindec], [sinra*cosdec, cosra, -sinra*sindec], [sindec, 0., cosdec]])
            B = np.dot(T, A) # T*A
            invB = np.linalg.inv(B)
            kine_samples_solar = solarmotion_sampler(name,kine_samples_n)
            kine_vhel = kine_samples_solar[:,2:] # samples in UVW_solar
            kine_vhel[:,1] += kine_samples_solar[:,1] # add the vdisk component
            solarreflex = np.einsum('ij,kj->ki',invB, kine_vhel)
            kine_samples_corrected[:,1:] -= solarreflex
            
            # If had no proper motions, then assume had exactly solar reflex
            if havepropermotions == False:
                kine_samples_corrected[:,1] = 0.0
                kine_samples_corrected[:,2] = 0.0
            
            # If had no radial velocity, then assume had exactly solar reflex
            if havevelocities == False:
                kine_samples_corrected[:,3] = 0.0
                
            # Calculate galactic velocity
            kine_vgrf = np.sqrt(kine_samples_corrected[:,1]**2+kine_samples_corrected[:,2]**2+kine_samples_corrected[:,3]**2)
            
            # Calculate escape velocity
            kine_galrad = np.sqrt(kine_samples_solar[:,0]**2+(kine_samples_corrected[:,0]*cosb)**2-2.*kine_samples_solar[:,0]*kine_samples_corrected[:,0]*cosb*cosl)
            kine_vesc = vesc(MWPotential2014,kine_galrad*un.kpc)*kine_samples_solar[:,1]
            
            # Bound probability
            kine_bound_n = np.where( kine_vesc > kine_vgrf )[0].shape[0]
            kine_betaa = kine_bound_n+0.5
            kine_betab = kine_samples_n-kine_bound_n+0.5
            kine_bound_percentiles = [beta.ppf(confidenceintervalcoverage/2.0,kine_betaa,kine_betab),beta.ppf(0.5,kine_betaa,kine_betab),beta.ppf(1.0-confidenceintervalcoverage/2.0,kine_betaa,kine_betab)]
            
            # Debugging
            #print(kine_values)
            #print(kine_errors)
            #print(kine_corr)
            #print(kine_mu)
            #print(kine_L,kine_k)
            #print(kine_galrad)
            #print(kine_samples_corrected[:,0])
            
            # Store all outcomes.
            source = catalog.entries[name].add_self_source()
            boundprobability_upperlimit = (havepropermotions == False or havevelocities == False)
            catalog.entries[name].add_quantity(FASTSTARS.BOUND_PROBABILITY, str(kine_bound_percentiles[1]), e_lower_value=str(kine_bound_percentiles[0]), e_upper_value=str(kine_bound_percentiles[2]), upperlimit = boundprobability_upperlimit, source=source, derived=True)
            catalog.entries[name].add_quantity(FASTSTARS.ESCAPE_VELOCITY, str(kine_vesc.mean()), e_value=str(kine_vesc.std()), u_value='km/s', source=source, derived=True)
            catalog.entries[name].add_quantity(FASTSTARS.VELOCITY, str(kine_vgrf.mean()), e_value=str(kine_vgrf.std()), u_value='km/s', lowerlimit = boundprobability_upperlimit, source=source, derived=True, kind='galactocentric')
            catalog.entries[name].add_quantity(FASTSTARS.LUM_DIST, str(kine_samples_corrected[:,0].mean()), e_value=str(kine_samples_corrected[:,0].std()), u_value='kpc', source=source, derived=True)
    catalog.journal_entries()
    
    
    return
