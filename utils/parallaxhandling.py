"""
"""

import numpy as np
import emcee
from math import log

__all__ = ['parallax_to_distance','kinematic_sampler','solarmotion_sampler']

meanvdisk = 238.0
sigmavdisk = 9.0
meanRsolar = 8.27
sigmaRsolar = 0.29
meanUsolar = 11.1
sigmaUsolar = np.sqrt(0.75**2+1.0**2)
meanVsolar = 12.24
sigmaVsolar = np.sqrt(0.47**2+2.0**2)
meanWsolar = 7.25
sigmaWsolar = np.sqrt(0.37**2+0.5**2)

def parallax_to_distance(NAME,PARALLAX,PARALLAX_ERROR,PHOTD=None,PHOTD_ERROR=None):
    # Constants
    L = 1.35
    nstep_burn = 100
    nstep = 50
    ndim, nwalkers = 1, 100
    N = nstep*nwalkers

    # Set random number generator based on name
    seed = NAME
    def craft_seed(SEED):
        nSEED = len(SEED)
        # The modulo ensures we don't overrun.
        indSEED = range(2*(nSEED+3))
        return np.prod(np.array([ord(SEED[i%nSEED]) for i in indSEED])) % 4294967296
    init_random_state = np.random.RandomState(seed=craft_seed(seed))
    samp_random_state = np.random.RandomState(seed=craft_seed(seed[::-1]))

    # Branch depending on whether we are using exponentially decreasing or photometric distance prior
    if PHOTD == None:
        # No photometric distance estimate, using Astraatmadja & Bailer-Jones (2016) prior
        lnprob_divisor = 1./(2.*PARALLAX_ERROR**2.)
        posterior = lambda x: -(1./x-PARALLAX)**2.*lnprob_divisor+2.*log(x)-x/L
        coeff = [1./L, -2., PARALLAX/(2.*PARALLAX_ERROR**2.), -1./(2.*PARALLAX_ERROR**2.)]
    else:
        # Photometric distance estimate used as prior
        lnprob_parallax_divisor = 1./(2.*PARALLAX_ERROR**2.)
        lnprob_photd_divisor = 1./(2.*PHOTD_ERROR**2.)
        posterior = lambda x: -(1./x-PARALLAX)**2.*lnprob_parallax_divisor-(x-PHOTD)**2.*lnprob_photd_divisor
        coeff = [1., -PHOTD, 0., PARALLAX*(PHOTD_ERROR/PARALLAX_ERROR)**2., -(PHOTD_ERROR/PARALLAX_ERROR)**2.]
        
    def lnprob(x):
        if x > 0:
            return posterior(x)
        else:
            return -np.inf
        
    # Identify maximum-likelihood location for initial positions
    r = np.roots(coeff)
    real_r = r.real[abs(r.imag)<1e-5]
    p0 = [max(real_r)+1e-1*init_random_state.randn(1) for i in range(nwalkers)]
    
    # Sample with burn-in
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
    pos, prob, state = sampler.run_mcmc(p0, nstep_burn, rstate0=samp_random_state.get_state())
    sampler.reset()
    sampler.run_mcmc(pos, nstep, rstate0=samp_random_state.get_state())
    distance=sampler.flatchain
    return np.median(distance), np.std(distance)
    
def kinematic_sampler(NAME,KINE_VALUES,KINE_ERRORS,KINE_CORR,KINE_L=1.35,KINE_K=3.0):
    # KINE_VALUES = (PARALLAX,PMRA,PMDEC,VRAD)

    # Constants
    nstep_burn = 100
    nstep = 100
    ndim, nwalkers = 4, 100
    N = nstep*nwalkers

    # Set random number generator based on name
    seed = NAME
    def craft_seed(SEED):
        nSEED = len(SEED)
        # The modulo ensures we don't overrun.
        indSEED = range(2*(nSEED+3))
        return np.prod(np.array([ord(SEED[i%nSEED]) for i in indSEED])) % 4294967296
    init_random_state = np.random.RandomState(seed=craft_seed(seed))
    samp_random_state = np.random.RandomState(seed=craft_seed(seed[::-1]))

    # Form covariance matrix and invert
    KINE_DIAG = np.diag(KINE_ERRORS)
    KINE_COV = np.dot(KINE_DIAG,np.dot(KINE_CORR,KINE_DIAG))
    KINE_COV_INV = np.linalg.inv(KINE_COV)

    # Find optimal start positions
    coeff = [1.0/KINE_L, -(KINE_K-1.0), KINE_VALUES[0]/KINE_ERRORS[0]**2., -1.0/KINE_ERRORS[0]**2.]
    r = np.roots(coeff)
    real_r = r.real[abs(r.imag)<1e-5]
    p0 = [np.array([real_r[0],KINE_VALUES[1],KINE_VALUES[2],KINE_VALUES[3]])+1e-1*init_random_state.randn(ndim) for i in range(nwalkers)]
    
    # Define log-posterior
    def lnprob(x):
        y = np.copy(x)
        y[0] = 1.0/y[0]
        if x[0] > 0:
            return -0.5*np.dot((y-KINE_VALUES).T,np.dot(KINE_COV_INV,(y-KINE_VALUES)))+(KINE_K-1.0)*log(x[0])-x[0]/KINE_L
        else:
            return -np.inf
    
    # Sample with burn-in
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
    pos, prob, state = sampler.run_mcmc(p0, nstep_burn, rstate0=samp_random_state.get_state())
    sampler.reset()
    sampler.run_mcmc(pos, nstep, rstate0=samp_random_state.get_state())
    
    return sampler.flatchain
    
def solarmotion_sampler(NAME,NSAMP):
    # Set random number generator based on name
    seed = NAME
    def craft_seed(SEED):
        nSEED = len(SEED)
        # The modulo ensures we don't overrun.
        indSEED = range(2*(nSEED+3))
        return np.prod(np.array([ord(SEED[i%nSEED]) for i in indSEED])) % 4294967296
    np.random.seed(seed=craft_seed(seed))
    
    # Sample from normal distributions
    RSOLAR = np.random.normal(meanRsolar,sigmaRsolar,NSAMP)
    VDISK  = np.random.normal(meanvdisk, sigmavdisk, NSAMP)
    USOLAR = np.random.normal(meanUsolar,sigmaUsolar,NSAMP)
    VSOLAR = np.random.normal(meanVsolar,sigmaVsolar,NSAMP)
    WSOLAR = np.random.normal(meanWsolar,sigmaWsolar,NSAMP)
    
    # Form of array is Rsolar,vdisk,Usolar,Vsolar,Wsolar
    SAMPLES_SOLAR = np.vstack([RSOLAR,VDISK,USOLAR,VSOLAR,WSOLAR]).T
    return SAMPLES_SOLAR
