# -*- coding: utf-8 -*-
"""
Lattice_QCD.py

This is a Metropolis Monte Carlo algorithm for Lattice QCD calculations
"""

import numpy as np
import random as r
import matplotlib.pyplot as plt 
from math import exp


"""
Declaring initial parameters and constants
"""

a = 1/2                         # Lattice spacing size
N = 20                          # Number of lattice sites
N_COR = 20                      # Correlation sample spacing  
N_CF = 10000                     # Number of configurations/paths
EPS = 1.4                       # Metropolis sweep step size
EPS_T = 0                       # A parameter for tuning EPS 
BIN_S = 4                       # Bin size
LEN_G = int(N_CF/BIN_S)         # Length of binned G 
BTS_N = 100                     # Number of bootstrapping copies 


def update(x): 
    """
    Monte Carlo update

    Randomly generates a new configuration x using the Metropolis algorithm
    and updates it if it meets the correct criteria. 
    :param Configuration x before MC update 
    :return: Configuration x after MC update
    """
    EPS_T = 0
    for j in range(0, N): 
        old_x = x[j]
        old_Sj = S(j,x)
        x[j] += np.random.uniform(-EPS,EPS)
        dS = S(j,x) - old_Sj
        if dS > 0 and exp(-dS) < np.random.uniform(0,1): 
            x[j] = old_x
            EPS_T += 1
    # Uncomment for tuning of EPS
    #print('Accept/reject ratio: {0}%'.format(EPS_T/N * 100))
    return x


def S(j,x): 
    """ 
    Defining the action S for harmonic oscillator 
    :param Lattice site j and configuration x
    :return: Returns the action S(j,x)
    """             
    jp = (j+1)%N                                
    jm = (j-1)%N
    return a*x[j]**2/2 + x[j]*(x[j]-x[jp]-x[jm])/a  


def compute_G(x,t):
    """
    Computes the correlation function for time slice t for one configuration x 
    """
    g = 0; 
    for j in range(0, N):
        g = g + x[j]*x[(j+t)%N]
    return g/N 


def compute_G_cubic(x,n):
    """
    Computes the correlation function for a time slice using cubic sources and sinks. 
    """
    g = 0; 
    for j in range(0, N):
        g = g + (x[j]**3)*(x[(j+n)%N]**3)
    return g/N 


def MCpaths(x,G): 
    """
    Paths are randomly generated using the update(x) and the correlation functions
    are calculated for each path and time slice and stored in a 2D array G
    """

    print("Starting MC for correlation function generation...\n")     

    for j in range(0,N):                        #Empties array to all zero  
        x[j]=0;                                 #Sometimes called a "cold start" 

    print("Thermalizing array...\n") 
    for j in range(0, 5*N_COR):                 #Thermalizes array
        update(x)
        if ((j+1) % 10 == 0): 
            print("Performed {0} updates".format(j+1)) 
    print("\nArray thermalized.\n")

    print("Calculating correlation functions...\n")  
    for alpha in range(0, N_CF):                #Calculates G[alpha][n]
        for j in range(0, N_COR):               #Perform N_COR updates before sampling  
            update(x)
        for n in range(0, N):                   #Sample correlation functions
            G[alpha, n] = compute_G(x,n) 
        if ((alpha+1) % 100 == 0):
            print("Sampled {0} correlation functions".format(alpha+1))  

    print("\nMC complete! Wrapping up...\n") 


def binning(G):
    """
    Creates binned copies of G of a specified binsize. 

    Notice that if binsize=N, then since the dimensions of G is N_cf*N 
    and for the binned version it is (1/binsize)*N_cf*N, then if binsize = N_cf, 
    then one simply obtains a matrix of size N which just corresponds to the 
    average of each time slice over all alpha configurations. 
    """

    G_binned = []
    
    for i in range(0, len(G), BIN_S):
        sumb_G = 0 
        
        for j in range(BIN_S):
            sumb_G += G[i+j]
        G_binned.append(sumb_G/BIN_S)
    
    return G_binned


def bootstrap(G):
    """
    Creates bootstrapped copies of the binned G. 
    """
    
    G_bootstrap = []
    for i in range(LEN_G):
        alpha = int(r.uniform(0, LEN_G))
        G_bootstrap.append(G[alpha])
    return G_bootstrap        


def ave_t_G(G):
    """
    For each time slice, calculates and averages G over all configurations
    and stores in an array t_aveG
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
    """
    Calculates the difference excitation energy delta E(t) for each time slice t
    and the uncertainty in G and delta E(t) and stores each as an array
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
    """
    Plots the theoretical and numerical energy differences
    """ 
   
    t_slices = np.arange(0.0, 10.0, a) # Break up time slices by lattice spacing 'a' 
 
    # Plot settings 
    plt.figure('Energy Difference') 
    plt.ylabel('$\Delta$E(t)')
    plt.xlabel('t')
    plt.ylim(0.5, 1.5)
    plt.xlim(-.5, 2.9)
    plt.legend() 

    # Theoretical Data 
    plt.plot((-.5, 6), (1, 1), lw = 1, color = 'Black') 

    # Experimental Data 
    plt.errorbar(t_slices, dE, yerr = unc_dE, fmt = 'o', markersize = 4, capsize = 3)
    #plt.plot(range(0, N), dE, 'bo', markersize = 4)


    plt.show()        


def plot_corr(ave_tot_G, unc_G): 
    """
    Plots the theoretical and numerical correlation functions
    """ 
    t = np.arange(0.0, 10.0, 0.01)  # Time slices 
    exp_t = (1/2)*np.exp(-t)        # Theoretical value

    t_slice = np.arange(0.0, 10.0, a) # Break up time slices by lattice spacing 'a'  

    # Plot settings  
    plt.figure('Correlation Function')
    plt.xlim(-.5, 10)
    plt.xlabel('t')
    plt.ylabel('G(t)')
    plt.legend(['']) 

    # Theoretical Data
    #plt.plot(t, exp_t, lw=1, color = 'Black')

    # Experimental Data 
    #plt.plot(t_slice, ave_tot_G, 'ro', markersize = 4)
    plt.errorbar(t_slice, ave_tot_G, yerr=unc_G, fmt ='-o', markersize = 4, capsize = 3)


    plt.show()



def main(): 
    """
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
    
   
if __name__ == "__main__": 
    main()
        
        
        
        
