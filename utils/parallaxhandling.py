"""
"""

import numpy as np
import emcee
from math import log

__all__ = ['parallax_to_distance']

# Constants
L = 1.35
nstep_burn = 100
nstep = 50
ndim, nwalkers = 1, 100
N = nstep*nwalkers

def parallax_to_distance(NAME,PARALLAX,PARALLAX_ERROR,PHOTD=None,PHOTD_ERROR=None):
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
    sampler.run_mcmc(pos, nstep)
    distance=sampler.flatchain
    return np.median(distance), np.std(distance)
