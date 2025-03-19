Monte Carlo Correlation Functions
=================================

Here is an explanation of the theory. 

.. Todo:: 

   Need to look into renaming the functions to use snake-case and have more useful names. 



The basic idea is that we have "paths" that are described by :math:`x(t)`. This 
simply describes a set of x-values at specific times. If we were to think of this 
as a set of discrete time-slices, we might have :math:`\{x(t_0), x(t_1), \ldots, x(t_N)\}`. 

We call this a "path" because it describes the position (of something) as a function of time.

Ignoring the physics, we usually want to compute things like, 

.. math:: 

   G(t) \equiv \langle \langle x(t_2) x(t_1) \rangle \rangle \equiv \frac{\int \mathcal{D} x(t) x(t_2) x(t_1) e^{-S[x]}}{\int \mathcal{D} x(t) e^{-S[x]}}



Instead, we simply approximate the correlation function by taking some 
finite set of samples of configurations and then average: 

.. math:: 
  
  G \approx \frac{1}{N_{cfg}} \sum_{i=1}^{N_{cfg}} G(x[i])

The more samples we have, the better the approximation. Unfortunately, this 
tends to produce very noisy signals, in practice, and requires the use of 
variance reduction techniques to improve.

In order to do this, however, we need to pick the method for setting up the configurations. 

To do this, we'll use a method call the Metropolis Algorithm. This involves starting with a 
configuration :math:`x(t)`, calculating the action :math:`S[x]` for that configuration, 
making a _random_ update to the configuration, seeing how the action changed, and then deciding 
to either keep the change or reject it based on a probability proportional to :math:`e^{-S[x]}`.

The basic outline of the Monte Carlo calculation of :math:`\langle \langle \Gamma[x] \rangle \rangle` for some function :math:`\Gamma[x]` depending on a path :math:`x` can be summarized as: 

1. Initialize the path (e.g. set every :math:`x_j` to zero)
2. Update the path :math:`x` using the Metropolis Algorithm until it's thermalized (about :math:`5N_{cor}-10N_{cor}`)
3. Compute :math:`\Gamma[x]` from the thermalized configuration :math:`x`
4. Update the path :math:`N_{cor}` times so that the path :math:`x` is statistically independent from the one previously used to meausre :math:`\Gamma[x]`.
5. Repeat 4 until enough samplings of :math:`\Gamma[x]` have been accumulated. 
