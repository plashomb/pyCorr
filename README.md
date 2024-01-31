# pyCorr 
A Python script for computing correlation functions and the energy difference between the 
ground state and first excited state of a quantum-mechanical harmonic oscillator. 


## Motivation 
This was written as my own choice for a final project for an introductory, elementary particle physics course that I took. 

I wrote it as a way to get a better sense of how lattice calculations are made for graduate research in Lattice QCD. 

## Purpose 

In the path integral formulation of quantum 
mechanics and quantum field theory, many 
physical observables can be extracted from 
correlation functions. In quantum mechanics, 
these come in the following form: 

$$G(t) \equiv \langle \langle x(t_2) x(t_1) \rangle \rangle = \frac{ \int \mathcal{D} x(t) x(t_2) x(t_2) e^{-S[x]}}{\int \mathcal{D} x(t) e^{e-S[x]}} $$


These kinds of integrals can be performed using 
Monte Carlo integration using algorithms like 
the Metropolis algorithm. This involves breaking 
up time into discrete "slices." 

If we then insert a complete set of energy eigenstates and look at large time separations, we find that this simplifies to, 

$$G(t) \xrightarrow{t\text{ large}} |\langle E_0 | \tilde{x} | E_1 \rangle|^2 e^{- (E_1 - E_0)}.$$

Here, $E_0$ and $E_1$ are the ground state 
energy and first excited state of the quantum 
mechanical harmonic oscillator. 

If we then compute $\log(G(t)/G(t+a))$, where $a$ is the lattice spacing, we find that, 

$$\log(G(t)/G(t+a)) \rightarrow (E_1 - E_0) a$$

In other words, computing this correlation function $G(t)$ gives us a way of being able to compute the difference in energy of the ground state and first excited state of the harmonic oscillator. 

To compute the correlation function, we use Monte Carlo integration on the integral.

To do this, we use the Metropolis algorithm to randomly generate paths $x$.
In doing so, we are effectively sampling 
the paths from the probability distribution $e^{-S[x]}$, which acts like a "weight" in the integral.
We can then simply average "samplings" of the integral, i.e. 

$$G(t) = \frac{1}{N} \sum_j \langle \langle x(t_j + t) x(t_j) \rangle \rangle$$
where $N$ is the number of lattice sites.

This leaves us with a simple procedure: 

1. "Thermalize" the lattice by updating it many times 
2. Wait $N_{CR}$ times and then sample the quantity we want
3. Repeat 2. Until we have enough samples 
4. Average the samples
5. Fit the averaged samples to expected functions of desired
parameters
6. Extract desired parameters 

In this program, the correlation function $G(t)$
is determined for several time slices and 
plotted. Furthermore, the energy difference 
between the first and second excited state 
of the harmonic oscillator is determined and 
compared to the theoretical value of 1, working 
in units of $\hbar \omega = 1$.

## References 
This program was *heavily* inspired by the following works:  

- [Lattice QCD for Novices](https://arxiv.org/pdf/hep-lat/0506036.pdf)
- [Monte Carlo Applications and Lattice QCD](https://digitalcommons.lib.uconn.edu/cgi/viewcontent.cgi?article=1502&context=srhonors_theses)