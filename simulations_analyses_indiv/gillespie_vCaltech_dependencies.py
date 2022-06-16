#Gillespie code is adapted from code here: http://be150.caltech.edu/2020/content/lessons/12_stochastic_simulation.html

import multiprocessing
import tqdm

import numpy as np
import scipy.stats as st
import numba

import biocircuits

# INTERACTIVE plotting modules
import bokeh.io
import bokeh.plotting

bokeh.io.output_notebook()

# for regulation:
def MM_TF_link(TF_mean, Hill_coef, max_effect):
    EC_50 = np.power(max_effect * TF_mean**Hill_coef - TF_mean**Hill_coef, 1.0/Hill_coef)
    def TF_mult_factor(tf_val):
        return max_effect * (tf_val**Hill_coef) / ((EC_50**Hill_coef) + (tf_val**Hill_coef))
    return TF_mult_factor
# MM_TF_link(60, 2, 16)(70)

def simple_propensity(propensities, population, t, kon, koff, kon2, koff2, beta_m, beta_p, beta_m2, gamma_m, gamma_p, gamma_m2):
    """Updates an array of propensities given a set of parameters (ie txn rate, decay rate),
    and an array of populations (key entities, ie TF mRNA, TF protein).
    """
    # Unpack population
    s, m, p, m2, s2 = population # bursting STATE, TF mRNA, TF protein, Target mRNA

    #incorporate Hill function for regulation
    kon2_adj = kon2*MM_TF_link(60, 4, 16)(p) #the first num should be mean(TF(P))
    # print(kon2, kon2_adj)
    
    # Update propensities
    # betas = production; gammas = decay
    propensities[0] = (1 - s)*kon ## TF switching to burst ON
    propensities[1] = s*koff ## TF switching to burst OFF
    propensities[2] = s*beta_m  # Make TF mRNA transcript (transcription rate * whether TF is bursting)
    propensities[3] = gamma_m * m # degrade TF mRNA
    propensities[4] = beta_p * m # Make TF protein
    propensities[5] = gamma_p * p  # Degrade TF protein
    propensities[6] = s*beta_m2 # Make target mRNA ## INCORPORATE BURST SWITCHING & HILL FUNCTION HERE
    propensities[7] = gamma_m2 * m2 # degrade Target mRNA
    propensities[8] = (1 - s2)*kon2_adj ## Target switching to burst ON
    propensities[9] = s2*koff2 ## Target switching to burst OFF
    
def sample_discrete(probs): #samples which reaction to run next
    """Randomly sample an index with probability given by probs."""
    # Generate random number
    q = np.random.rand()

    # Find index
    i = 0
    p_sum = 0.0
    while p_sum < q:
        p_sum += probs[i]
        i += 1
    return i - 1

def gillespie_draw(propensity_func, propensities, population, t, args=()):
    """
    Draws a reaction and the time it took to do that reaction.

    Parameters
    ----------
    propensity_func : function
        Function with call signature propensity_func(population, t, *args)
        used for computing propensities. This function must return
        an array of propensities.
    propensities : ndarray
        Propensities for each reaction as a 1D Numpy array.
    population : ndarray
        Current population of particles (key entities of interest, ie TF(RNA), TF(P))
    t : float
        Value of the current time.
    args : tuple, default ()
        Arguments to be passed to `propensity_func`.

    Returns
    -------
    rxn : int
        Index of reaction that occured.
    time : float
        Time it took for the reaction to occur.
    """
    # Compute propensities
    propensity_func(propensities, population, t, *args)

    # Sum of propensities
    props_sum = propensities.sum()

    # Compute next time
    time = np.random.exponential(1.0 / props_sum) ## exponentially distributed state-switching

    # Compute discrete probabilities of each reaction
    rxn_probs = propensities / props_sum

    # Draw reaction from this distribution
    rxn = sample_discrete(rxn_probs) #sample_discrete(rxn_probs) or sample_discrete_scipy(rxn_probs)

    return rxn, time

def gillespie_ssa(propensity_func, update, population_0, time_points, args=()):
    """
    Uses the Gillespie stochastic simulation algorithm to sample
    from probability distribution of particle counts over time.

    Parameters
    ----------
    propensity_func : function
        Function of the form f(params, t, population) that takes the current
        population of particle counts and returns an array of propensities
        for each reaction.
    update : ndarray, shape (num_reactions, num_chemical_species)
        Entry i, j gives the change in particle counts of species j
        for chemical reaction i.
    population_0 : array_like, shape (num_chemical_species)
        Array of initial populations of all chemical species.
    time_points : array_like, shape (num_time_points,)
        Array of points in time for which to sample the probability
        distribution.
    args : tuple, default ()
        The set of parameters to be passed to propensity_func.

    Returns
    -------
    sample : ndarray, shape (num_time_points, num_chemical_species)
        Entry i, j is the count of chemical species j at time
        time_points[i].
    """

    # Initialize output
    pop_out = np.empty((len(time_points), update.shape[1]), dtype=np.int)

    # Initialize and perform simulation
    i_time = 1
    i = 0
    t = time_points[0]
    population = population_0.copy()
    pop_out[0, :] = population
    propensities = np.zeros(update.shape[0])
    while i < len(time_points):
        while t < time_points[i_time]:
            # draw the event and time step
            event, dt = gillespie_draw(propensity_func, propensities, population, t, args)

            # Update the population
            population_previous = population.copy()
            population += update[event, :]

            # Increment time
            t += dt

        # Update the index
        i = np.searchsorted(time_points > t, True)

        # Update the population
        pop_out[i_time : min(i, len(time_points))] = population_previous

        # Increment index
        i_time = i

    return pop_out


