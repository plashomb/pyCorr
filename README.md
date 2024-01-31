# Python Correlation Function 
A Python script for computing correlation functions and energy differences for a quantum-mechanical harmonic oscillator. 


## Motivation 
This was written as my own choice for a final project for an introductory, elementary particle physics course that I took. 

I wrote it as a way to get a better sense of how lattice calculations are made for graduate research in Lattice QCD. 

## Theory 
This project is an example of Monte Carlo integration. In this case, the integral is the path integral of quantum mechanics.

Before we discuss this, first we should look at a simple 
example of averages. Averages, as you probably know, look 
something like this, 

$$\bar{x} = \frac{\sum_i^N x_i}{\sum_i^N1} = \frac{x_1 + x_2 + ... + x_N}{N},$$
where $x_1, x_2, ..., x_N$ are each of the numbers to average, $\bar{x}$ is the average, and $N$ is the total 
number of things to average. 
This example is "discrete", which is to say, that there 
are a set number of variables to average. 

We can, however, have averages of functions of continuous 
variables defined within a specific region. For instance, 

$$\bar{f} = \frac{\int_a^b dx f(x)}{\int_a^b dx} = \frac{\int_a^b dx f(x)}
{b-a}.$$

Recall from the Riemann sum, we have that, 

$$\int_a^b f(x) dx = \lim_{||\Delta x||\rightarrow 0} \sum_{i=1}^N f(x_i) \Delta x_i$$

Note that this is **__not__** always true of every 
function $f(x)$. In the cases we will look at, however,
it will be. 

This can be used to approximate this contunous integral by, 

$$\bar{f} \approx \frac{1}{b-a}\frac{1}{N}\sum_i^N f(x_i) (b - a) = \frac{1}{N}\sum_i^N f(x_i)$$
where, here, we take a large number of values of $x$, 
we evaluate the function $f(x)$ at those points, and then
we average the results. If we use a lot of values $x_i$,
then we would hope to obtain a better and better 
approximation of the integral. Not that, of course, the 
values $x_i$ are, of course, within the range $[a, b]$, 
in this case.

Notice that these points $x_i$ are each selected at equal
distance from each other in the Riemann sum. 
In a Monte Carlo simulation, on the other hand, we sample 
these values $x_i$ *_randomly_*. 


This is fine if every value is equiprobably, or equally 
likely to occur in our average. However, we might expect
that in certain cases, some values are more likely to 
occur than others will be. In cases like this, we need to 
weight the values different from each other using a 
"weighted average." For the discrete case, we have, 

$$\bar{x} = \frac{\sum_i^N w_i x_i}{\sum_i^N w_i} = \frac{w_1 x_1 + w_2 x_2 + ... + w_N x_N}{w_1 + w_2 + ... + w_N},$$
Notice this time that we don't simply have $N$ in the denominator, this time. 

For the continuous case, we have, 

$$\bar{f} = \frac{\int_a^b dx f(x) g(x)}{\int_a^b g(x)dx}.$$

How would we do a Monte Carlo integration this time? In 
this case, we can't simply approximate the integral the 
way we did before because now there is an extra $g(x)$ 
weight inside of the integral. 

The trick is to treat this weight like a probability 
distribution, and to sample our values $x_i$ based on 
this probability distribution. In other words, 

$$P[x_i] \propto g(x) $$





What about the following integral? 

$$ \langle \langle \Gamma[x] \rangle \rangle = \frac{ \int \mathcal{D} x(t) \Gamma[x] e^{-S[x]}}{\int \mathcal{D} x(t) e^{e-S[x]}} $$

While this may look daunting, notice that we have the same 
sort of situation as the weighted average of the function 
$f$ from before. In this case, the function $g(x)$ is 
replaced with the functional $e^{-S[x]}$ that depends on 
the action. 

In other words, the exponential is weighting the paths 
and certain paths will contribute more to the integral 
than others will. 

What we do is form whta is called a "markov chain" until 
the lattice is "thermalized." 

Sampling from these "configurations" is then equivalent to 
sampling $x(t)$ where, 

$$P[x(t)] \propto e^{-S[x]}.$$

How do we determine these configurations, though? We do 
this using an algorithm called the Metropolis Hastings 
algorithm. By following this algorithm, we can generate 
a set of configurations that 

If we then generate $N_{CFG}$ of these configurations, 
$\{x_1(t), x_2(t), ..., x_{N_{CFG}}(t)\}$, then we can use these in our averaging of the quantity we want to get, 

$$\langle \langle \Gamma[x] \rangle \rangle = \frac{1}{N_{CFG}} \sum_{i=1}^{N_{CFG}} \Gamma[x_i]$$


This integral is summing up all of the paths, 


In the theory, there is an important quantity called a correlation function, represented as $G(t)$ which can be written as, 





$$G(t) \equiv \langle \langle x(t_2) x(t_1) \rangle \rangle = \frac{ \int \mathcal{D} x(t) x(t_2) x(t_2) e^{-S[x]}}{\int \mathcal{D} x(t) e^{e-S[x]}} $$



This quantity as an important property. Suppose we look at the ground state and first excited state. If you look at large time seperations, you can find that, 

$$G(t) \xrightarrow{t\text{ large}} |\langle E_0 | \tilde{x} | E_1 \rangle|^2 e^{- (E_1 - E_0)}$$

If we then compute $\log(G(t)/G(t+a))$, where $a$ is the lattice spacing, we find that, 

$$\log(G(t)/G(t+a)) \rightarrow (E_1 - E_0) a$$

In other words, computing this correlation function $G(t)$ gives us a way of being able to compute the difference in energy of the ground state and first excited state of the harmonic oscillator. 

To compute the correlation function, we use Monte Carlo integration on the integral.

To do this, we use an algorithm called the Metropolis-Hastings algorithm to randomly generate paths $x$.
By using Metropolis-Hastings, we are effectively sampling 
the paths from the probability distribution $e^{-S[x]}$, which acts like a "weight" in the integral.
We can then simply average "samplings" of the integral, i.e. 

$$G(t) = \frac{1}{N} \sum_j \langle \langle x(t_j + t) x(t_j) \rangle \rangle$$
where $N$ is the number of lattice sites.

In other words, $G(t)$ in the numerical integration works out to be just a 2D array of numbers. For instance, let's say we pick $t=0$. Then, we have an array, 

```python
G_0 = sum([x(0)*x(0), x(1)*x(1), ..., x(N-1)*x(N-1)])
G_0 = G_0 / N
```
This is where $\langle \langle x(t_2)x(t_1) \rangle \rangle$ is computed from the same time slice $t_2 = t_1$. 

For one time slice away, we can get, 

```python
G_1 = sum([x(0+1)*x(0), x(1+1)*x(1), ..., x(N-1+1)*x(N-1)])
G_1 = G_1 / N 
```
Then, the correlation function $G(t)$ is just the array, 

```python 
G = [G_0, G_1, G_2, ..., G_N]
```

This leaves us with a simple procedure: 

1. "Thermalize" the lattice by updating it many times 
2. Wait $N_{CR}$ times and then sample the quantity we want
3. Repeat 2. Until we have enough samples 
4. Average the samples
5. Fit the averaged samples to expected functions of desired
parameters
6. Extract desired parameters 



## References 
This program was heavily inspired by the following works:  

- [Monte Carlo Applications and Lattice QCD](https://digitalcommons.lib.uconn.edu/cgi/viewcontent.cgi?article=1502&context=srhonors_theses)

- [Lattice QCD for Novices](https://arxiv.org/pdf/hep-lat/0506036.pdf)