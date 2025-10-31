# growthEstimation: R Package for Fisheries Growth and Mortality Analysis

growthEstimation is an R package for analyzing age-at-length and length-frequency data from fisheries surveys. Its primary purpose is to simultaneously estimate von Bertalanffy growth parameters and mortality rates, while correcting for the significant biases introduced by length-stratified sampling.

## **The Problem: Bias in Length-Stratified Samples**

Age data, is often collected using a length-stratified sampling design (i.e., sampling a specific number of fish from each length class.
This sampling method creates a well-known bias:

* At any given length, fish that grow faster are younger.  
* Fish that grow slower are older.  
* Because slower-growing fish take longer to reach a specific length class, they are exposed to natural mortality for a longer period.  
* Therefore, slower-growing individuals are underrepresented in the sample compared to faster-growing individuals of the same length.

Any attempt to deduce growth parameters from this data without simultaneously accounting for mortality will be **highly biased**, typically resulting in an overestimation of the average growth rate.

## **The Model: Stochastic Growth and Mortality**

The growthEstimation package solves this problem by modeling the underlying population dynamics of a cohort. The model is built on the principle that individual growth is stochastic.

1. **Stochastic Growth (SDE):** The growth of a single individual is not deterministic. We model it using a stochastic differential equation (SDE). The drift (mean growth) is described by the von Bertalanffy growth function, and the stochastic component is modeled as a Wiener process (Brownian motion) with a size-dependent diffusion rate.  
2. **Cohort Dynamics (Fokker-Planck Equation):** The evolution of the entire cohort's size distribution over time, subject to both stochastic growth and size-dependent mortality, is described by an advection-diffusion-reaction PDE. This package numerically solves this PDE to find the size distribution of the cohort as a function of age.  
3. **Spawning and Survey Timing:** The model does not assume fish are born on a single day. It incorporates a continuous spawning season (e.g., modeled by a von Mises distribution) to simulate a realistic population structure.  
4. **True Age vs. Observed Annuli:** The model maps a fish's **true age** (a continuous variable) to the **observed annuli** (a discrete count). This mapping function explicitly accounts for the fish's birth date (from the spawning model), the survey date, and the assumed date of annulus formation.

By convolving the cohort's dynamic evolution with the spawning distribution, we can predict the expected distribution of annuli for any length class on any given survey date. The model parameters (growth, diffusion, mortality) are then estimated by maximizing the likelihood of observing the empirical data.

### **Key Features**

* **Simultaneous Estimation:** Estimates von Bertalanffy growth parameters ($L_\infty$, $k$, $t_0$) and mortality rates ($M$) at the same time.  
* **Bias Correction:** Explicitly models the interplay between stochastic growth, mortality, and length-stratified sampling to remove the inherent bias.  
* **Stochastic Growth:** Uses an SDE to model individual growth variability, providing a more realistic representation of the population.  
* **Population Dynamics:** Solves the advection-diffusion-reaction PDE to model the cohort's size distribution over time.  
* **Realistic Timing:** Accounts for seasonal spawning patterns and discrete survey dates.  
* **Annuli Mapping:** Includes a robust function to translate true age to observed annuli, which is critical for likelihood calculations.


For detailed examples and a full walkthrough, please see "[Getting started](articles/growthEstimation.html)".

## **Vignettes**

This package provides two key vignettes to explain the underlying methodology:

1. **[Modelling of age-at-size data](articles/age_simulation.html)** This document explains the complete statistical model, including the cohort evolution, the incorporation of the spawning season, and the crucial mapping from a fish's true age to its observed annuli. It also details the multinomial likelihood function used for parameter estimation.  
2. **[Numerical scheme for solving PDE](articles/numerical_scheme.html)** This document details the numerical methods used to solve the Fokker-Planck equation. It describes the finite volume discretisation, the use of the Thomas algorithm for solving the resulting tridiagonal system, and a stability analysis that proves the scheme is unconditionally stable.
