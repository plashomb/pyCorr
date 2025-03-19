# -*- coding: utf-8 -*-
"""Metropolis Algorithm for estimating correlation functions. 

This module computes correlation functions for the quantum mechanical harmonic oscillator 
and estimates the energy difference between the ground and first excited state. 

Set the global parameters for the simulation, and run by typing, 

.. code-block:: console 


  $ source ./venv/bin/activate 
  (venv) $ python3 src/Lattice_QCD.py

Be sure to activate the virtual environment first. 

Global Variables: 
    :a (float64): The lattice spacing size
    :N (int): Number of lattice sites
    :N_COR (int): Correlation sample spacing  
    :N_CF (int): Number of configurations/paths
    :EPS (float): Metropolis sweep step size
    :BIN_S (int): Bin size
    :LEN_G (int): Length of binned G 
    :BTS_N (int): Number of bootstrapping copies 


"""

import numpy as np
from numpy.typing import NDArray, DTypeLike 
import random as r
import matplotlib.pyplot as plt 
from math import exp

import logging  


"""
Declaring initial parameters and constants
"""

log = logging.getLogger(__name__)
logging.basicConfig(format='%(asctime)s %(levelname)s:  %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

log.info('Initializing parameters')

a = 1/2                         # Lattice spacing size
N = 20                          # Number of lattice sites
N_COR = 20                      # Correlation sample spacing  
N_CF = 1000                     # Number of configurations/paths
EPS = 1.4                       # Metropolis sweep step size
BIN_S = 4                       # Bin size
LEN_G = int(N_CF/BIN_S)         # Length of binned G 
BTS_N = 100                     # Number of bootstrapping copies 



def update(x: NDArray[np.float64]) -> NDArray[np.float64]:
    r"""Performs Monte Carlo update

    Randomly generates a new configuration x using the Metropolis algorithm
    and updates it if it meets the correct criteria.

    Note: 
        This uses a for-loop via range(), so it will be slow.

    Todo: 
        Need to provide much better documentatoin of this portion, and fix the issues with 
        tuning EPS. 
    
    Args: 
        x (NDArray[np.float64]): Configuration :math:`x` before Monte Carlo update 
    
    Returns: 
        NDArray[np.float64]: The configuration :math:`x` after performing a Monte Carlo update

    """
    times_accepted = 0
    for j in range(0, N): 
        old_x = x[j]
        old_Sj = S(j,x)
        x[j] += np.random.uniform(-EPS,EPS)
        dS = S(j,x) - old_Sj
        if dS > 0 and exp(-dS) < np.random.uniform(0,1): 
            x[j] = old_x
            times_accepted += 1
    # Uncomment for tuning of EPS
    acceptance_rate = times_accepted/N * 100 
    #if (acceptance_rate < 40 or acceptance_rate > 60):
    #    log.warning(f"Accept/reject ratio: {times_accepted/N * 100}%")
#    else: 
#        log.info(f"Accept/reject ratio: {times_accepted/N * 100}%")
    return x


def S(j,x): 
    r"""Computes the action S for the quantum mechanical harmonic oscillator 

    Defines the action :math:`S` for the harmonic oscillator. The discretization 
    used is given by the following: 

    .. math:: 

        S[x] = \frac{a}{2}x^2_j + \frac{x_j - x_{j+1} - x_{j-1}}{a} x_j 

    This forms the first-order discretization of the action. 

    Args:  
        j (int): The lattice site :math:`j` 
        x (NDArray[np.float64]): The configuration :math:`x`

    Returns:  
        NDArray[np.float64]: The action :math:`S(j,x)`
    """             
    jp = (j+1)%N                                
    jm = (j-1)%N
    return a*x[j]**2/2 + x[j]*(x[j]-x[jp]-x[jm])/a  


def compute_G(x: NDArray[np.float64],t: np.float64) -> np.float64:  
    r"""Computes the correlation function 

    Computes the estimator for the correlation function 
    :math:`\langle G(x,t) \rangle` for time slice :math:`t` 
    for one configuration :math:`x(t)`. Specifically, it computes, 

    .. math::
        
        G(t) = \frac{1}{N} \sum_j \langle \langle x(t_j + t) x(t_j) \rangle \rangle

    In other words, G(t) will be a vector where :math:`t` indicates how many time-slices 
    of separation between :math:`x(t_j)` and :math:`x(t_j + t)` appearing in the summation.
    This calculation is performed on only a single path or "configuration" :math:`x^{(\alpha)}`.

    Todo: 
        This can be made more pythonic without using the for loops. 

    Args:
        x (NDArray[np.float64]): The configuration :math:`x(t)`
        t (np.float64): The timeslice `t`

    Returns: 
        np.float64: The correlation function :math:`\langle G(x,t) \rangle`. 
        
    """
    g = 0; 
    for j in range(0, N):
        g = g + x[j]*x[(j+t)%N]
    return g/N 


def compute_G_cubic(x,n):
    r"""Computes the correlation function for cubic source 

    Computes the correlation function for a time slice using cubic sources and sinks. 
    """
    g = 0; 
    for j in range(0, N):
        g = g + (x[j]**3)*(x[(j+n)%N]**3)
    return g/N 


def MCpaths(x,G): 
    r"""Performs Monte Carlo to generate the correlation function 

    Paths are randomly generated using :py:func:`update` until each configuration :math:`x^{(k)}`
    sampled is such that, :math:`Pr[x] \propto e^{-S[x]}`. Once this is the case, the correlation 
    function is computed for a path and at every timeslice. A few more updates are performed, and  
    another sampling of the correlation function is computed. These in-between updates are done 
    to ensure that no two configurations are related to each other. This allows us to form a 
    2D array for the correlation function.
   
    
    Todo: 
        * This is the bulk of the calculation, so should look into using Cython to speed up loops.
        * This needs a lot better documentation, and needs to be fixed. Currently takes G and x as 
          arguments, but never uses them.
        * If they need to take a "default" value for x and G, should provide that in function signature. 

    Args: 
        x (NDArray[np.float64]): 
    
    """

    log.info("Starting MC for correlation function generation...")     

    x = np.zeros(N) #Initializing to zero performs a "cold start" 

    # Performs thermalization of the array 
    # by repeatedly running update() on the  
    # configurations x until thermalized 

    log.info("Thermalizing array...") 
    for j in range(0, 5*N_COR):
        update(x)
        if ((j+1) % 10 == 0): 
            log.info(f"Performed {j+1} updates") 
    log.info("Array thermalized.")

    log.info("Calculating correlation functions...")
    for alpha in range(0, N_CF): #Calculates G[alpha][n]
        for j in range(0, N_COR): # Perform N_COR updates between each sampling  
            update(x)
        for n in range(0, N): 
            G[alpha, n] = compute_G(x,n) # Sample correlation functions  
        if ((alpha+1) % 100 == 0):
            log.info(f"Sampled {alpha+1} correlation functions")  

    log.info("Monte Carlo complete! Wrapping up...") 


def binning(G):
    r"""Performs binning on the correlation function 

    Takes the "ensemble" of configuration functions :math:`G^{(1)}, G^{(2)}, \ldots`
    and divides it up into binned averages as,  

    .. math:: 

      \begin{matrix} 
      \bar{G}^{(1)} \equiv \frac{G^{(1)} + G^{(2)} + G^{(3)} + G^{(4)}}{4}\\
      \bar{G}^{(2)} \equiv \frac{G^{(5)} + G^{(6)} + G^{(7)} + G^{(8)}}{4}\\
      \bar{G}^{(3)} \equiv \frac{G^{(9)} + G^{(10)} + G^{(11)} + G^{(12)}}{4}\\
      \vdots 
      \end{matrix}

    Notice that if binsize=N, then since the dimensions of G is N_cf*N 
    and for the binned version it is (1/binsize)*N_cf*N, then if binsize = N_cf, 
    then one simply obtains a matrix of size N which just corresponds to the 
    average of each time slice over all alpha configurations.

    Args: 
        G (NDArray[np.float64]): The correlation function ensemble 

    Returns:
        NDArray[np.float64]: The binned correlation function 
    """

    G_binned = []
    
    for i in range(0, len(G), BIN_S):
        sumb_G = 0 
        
        for j in range(BIN_S):
            sumb_G += G[i+j]
        G_binned.append(sumb_G/BIN_S)
    
    return G_binned


def bootstrap(G):
    r"""Performs bootstrap sampling of the binned correlation function  
    
        Bootstrap is a type of resampling that takes the correlation function
        values :math:`\{G_0, G_1, G_2, \ldots\}` and produces a set of new 
        values based on randomly picking from the full set, with possible 
        duplicates. This produces a whole new set. This is repeated as many 
        times as required in order to obtain errorbar estimates from the new 
        samples.

    Todo: 
        The current implementation is not very pythonic. Could simply do it in a single 
        list comprehension
    """

    G_bootstrap = []
    for i in range(LEN_G):
        alpha = int(r.uniform(0, LEN_G))
        G_bootstrap.append(G[alpha])
    return G_bootstrap        


def ave_t_G(G):
    r"""Computes the configuration-averaged correlation function 

    For each time slice, the correlation function :math:`G` is averaged over all configurations
    and the result is stored in the array t_aveG. This produces the Monte Carlo estimage of 
    :math:`langle G(t) rangle` over some discrete set of time slices. Or, in other words, 
    we have :math:`{G(t_0), G(t_1), ldots}`

    :return: Configuration-averaged G(t)
    """
    t_aveG = []

    for n in range(0, N):
        aveG = 0
        for alpha in range(LEN_G):
            aveG += G[alpha][n]
        aveG /= LEN_G 
        t_aveG.append(aveG)
    return t_aveG            


def stdCalc_dE(ave_tot_G, btsp):
    r"""
    Calculates the difference excitation energy delta E(t) for each time slice t
    and the uncertainty in G and delta E(t) and stores each as an array

    At large time slices, it's expected that: 

    .. math:: 

        \log\Big( \frac{G(t)}{G(t+a)} \Big) \rightarrow (E_1 - E_0)a,

    meaning the difference in energy can be computed from the ratio of correlated functions
    at adjacent time-slices. Specifically, this gives estimates as, 

    .. math:: 
        \Delta E = \frac{1}{a} \log\Big( \frac{G_{n}}{G_{n+1}} \Big) 

    The theory prediction is that this energy difference for the quantum mechanical harmonic 
    oscillator between the ground and first excited state should be unity. 

    Args: 
        ave_tot_G (NDArray[np.float64]): The configuration-average correlation functions 
        btsp (NDArray[np.float64]): The bootstrapped copies 

    Returns: 
        [np.float64]: The energy difference, its uncertainty, and the correlation function error bars

    """

    unc_dE = []
    unc_G = []
    dE = []  
    
    for n in range(N):
        err_1 = []
        err_2 = []
        dE.append((1/a)*np.log(np.abs(ave_tot_G[n]/ave_tot_G[(n+1)%N])))
        
        for beta in range(BTS_N):
            err_1.append((1/a)*(np.log(np.abs(btsp[beta][n] / btsp[beta][(n + 1)%N]))))
            err_2.append(btsp[beta][n])
            
        unc_dE.append(np.std(err_1))
        unc_G.append(np.std(err_2))
    
    return unc_dE, unc_G, dE 


def plot_energy(dE, unc_dE):
    r"""Plots the theoretical and numerical energy differences and their uncertainties

    Args: 
        dE (NDArray[np.float64]): The numerical estimate of energy differences :math:`\Delta E` 
        unc_dE (NDArray[np.float64]): The uncertainties of the energy differences 
    """ 
   
    t_slices = np.arange(0.0, 10.0, a) # Break up time slices by lattice spacing 'a' 
 
    # Plot settings 
    plt.figure('Energy Difference') 
    plt.ylabel('$\Delta$E(t)')
    plt.xlabel('t')
    plt.ylim(0.5, 1.5)
    plt.xlim(-.5, 2.9)
    #plt.legend() 

    # Theoretical Data 
    plt.plot((-.5, 6), (1, 1), lw = 1, color = 'Black') 

    # Experimental Data 
    plt.errorbar(t_slices, dE, yerr = unc_dE, fmt = 'o', markersize = 4, capsize = 3)
    #plt.plot(range(0, N), dE, 'bo', markersize = 4)


def plot_corr(ave_tot_G, unc_G): 
    r"""Plots the configuration-average correlation function and uncertainties

    Args: 
        ave_tot_G (NDArray[np.float64]): The configuration-averaged correlation function 
        unc_G (NDArray[np.float64]): The correlation function uncertainties
    """ 
    t = np.arange(0.0, 10.0, 0.01)  # Time slices 
    exp_t = (1/2)*np.exp(-t)        # Theoretical value

    t_slice = np.arange(0.0, 10.0, a) # Break up time slices by lattice spacing 'a'  

    # Plot settings  
    plt.figure('Correlation Function')
    plt.xlim(-.5, 10)
    plt.xlabel('t')
    plt.ylabel('G(t)')
    #plt.legend(['']) 

    # Theoretical Data
    #plt.plot(t, exp_t, lw=1, color = 'Black')

    # Experimental Data 
    #plt.plot(t_slice, ave_tot_G, 'ro', markersize = 4)
    plt.errorbar(t_slice, ave_tot_G, yerr=unc_G, fmt ='-o', markersize = 4, capsize = 3)




def main(): 
    r"""
    Declaration of Main function 
    """
    
    # Defining path x as an array and the correlation function as a 2D array
    x = np.zeros(N);
    G_cor = np.zeros((N_CF, N)); 
    btsp = []
    ave_tot_G = []

    # Creates the correlation functions for each time slice and configuration
    MCpaths(x, G_cor)
   
    # Use bootstrapping to create extra data for extracting statistics 
    binnedG = binning(G_cor)
    for _ in range(BTS_N):    
        btsp.append(ave_t_G(bootstrap(binnedG)))
        
    # Compute averages and uncertainties over configurations  
    ave_tot_G = ave_t_G(binnedG)                        # Averages
    unc_dE, unc_G, dE = stdCalc_dE(ave_tot_G, btsp)     # Uncertainties 

    # Plot results 
    plot_energy(dE, unc_dE)
    plot_corr(ave_tot_G, unc_G) 
    
    plt.show() 

if __name__ == "__main__": 
    main()
        
        
        
        
