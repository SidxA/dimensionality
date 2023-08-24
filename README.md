repository structure

main functionalities
struct of dimension reduction
running of dimension reduction
processing of the results
saving the results
loading the results
plotting motivation
plotting example
plotting heatmaps
plotting workflow


todo
move preprocessing to toolkit
readme

# Code Descriptions

This README provides an overview of the content in each of the Julia files included in this repository.

## 01_basic_functionality.jl

- includes essential functions and utilities that serve as the foundation for other scripts in this project

## 02_dimensionreduction_structs.jl

- defines custom data structures and types required for dimension reduction algorithms
- provides the necessary abstractions and data representations for efficient data processing and analysis

## 03_data_preprocessing.jl

- focuses on data preprocessing techniques
- lowpass filter at different frequencies

## 04_RUN_dimensionreduction.jl

- main script for executing dimension reduction algorithms
- run `export JULIA_NUM_THREADS=1` in bash and start julia with desired number of cores to specify CPU usage
- the dimensionality reduction is performed 1 timeseries per core using the same W values to re-write the big matrices to save RAM
- each individual calculation outputs a SSA and NLSA jld2 file with the modes,etc. inside into the specified `outdir`
- just run `/opt/julia-1.9.0/bin/julia --threads 64 04_RUN_dimensionreduction.jl` for using 64 cores
- with the specified parameters changed inside
  - `savedirname = "/net/scratch/lschulz/data/time_series.jld2"` the file containing the `[N,spots,variables]` different values for different filters
  - `wholedata = SharedArray{Float32}(load(savedirname)["data_x"])`with `x` in `[raw,f3,f4,f6]` for corresponding dataset
  - this file contains 18 spots and 16 variables selected for beeing the longest measurement periods of the variables of interest
  - the desired different `W` need to be specified in the list at the parallel loop at the end of the document: this is run at 7 a

## 05_processing.jl

- contains data processing routines, applied after dimension reduction
- includes feature extraction, seasonal cycle 

## 06_load_figuredata.jl

- loading and preparing data for generating figures
- it is important if julia is run externally to run it by `GKSwstype=nul /opt/julia-1.9.0/bin/julia` in order to prevent it from trying to display images directly
- 
## 07_fig1_motivation.jl

- code for generating Figure 1
- illustrates the motivation and background information: need for higher harmonic terms then the fundamental in order to resolve complex climatic time-series

## 08_fig2_results.jl

- contains the code for generating Figure 2
- 3 examples of dimension reduction to extract the seasonal cycle in different spots and with different vegetation related flux data
- `mode_figure_flags` builds individual overview panel with timeseries (with quality flags), spectrum, mode shapes, mode spectral content
    - e.g.
`spot = 2
vari = 2
F = Figure()
p = local_parameters(spots[spot],vars[vari],outdir)
mode_figure_flags(F,p,"test",flags[:,spot,vari],data_tensor)
save(savedirname,F)`
- `large_mode_figure_flags` combines 3 selected examples into a large figure

## 09_fig3_heatmaps.jl

- includes the code for generating Figure 3
- `heatmap_numbers` produces heatmap with numbers, focused solely on differences in extracted number of harmonics between different noise filters
- `characterization_heatmap` also includes little symbolls for noise content, signal sinusoidality, signal entropy (NN sample entropy) and quality flag artefacts

## 10_workflow_figures.jl

- pictures used in the workflow diagram

# DATA


time series raw data and lowpass filtered dataand description of locations,variables,ecosystem types, qc quality flags
`"/net/scratch/lschulz/data/time_series.jld2"`

individual dimension-reduction results
`outdir_raw = "/net/scratch/lschulz/data/dimensionreduction/raw/"
outdir_f3 = "/net/scratch/lschulz/data/dimensionreduction/lowpass3/"
outdir_f4 = "/net/scratch/lschulz/data/dimensionreduction/lowpass4/"
outdir_f6 = "/net/scratch/lschulz/data/dimensionreduction/lowpass6/"`

seasonal cycle prior to the filter of only considering intact mode pairs (for old heatmap)
`"/net/scratch/lschulz/data/seasonal_cycle.jld2"`


data characteristics 
`"/net/scratch/lschulz/data/data_characteristics_without_f3.jld2"`

new seasonal cycle information based on the pair filtering (for new heatmap)
`"/net/scratch/lschulz/data/trends.jld2"`



# Theory

This repository contains code for the time-series analysis of climatic data gathered from the FLUXNET network. The analysis is based on two dimensionality-reduction methods: Singular Spectrum Analysis (SSA) and Nonlinear Laplacian Spectral Analysis (NLSA).

## Introduction

The FLUXNET network is a global network of eddy covariance flux tower sites that provide high-frequency measurements of carbon, water, and energy fluxes between the terrestrial biosphere and the atmosphere. The FLUXNET dataset is a valuable resource for studying the dynamics of terrestrial ecosystems and their response to environmental drivers.

Time-series analysis is a powerful tool for understanding the complex dynamics of the earth system. Dimensionality-reduction methods such as SSA and NLSA are useful for identifying patterns and trends in time-series data. SSA is a data-driven method for decomposing time-series data into a small number of underlying components, while NLSA is a nonlinear extension of traditional spectral analysis methods that can capture complex nonlinear dynamics in time-series data.

## Methodology

The code in this repository is written in Julia 1.8.5 and implements both SSA and NLSA for the analysis of FLUXNET data. The data is delay-embedded, and then the dimensionality-reduction methods are applied to the resulting time-series.

For SSA, the code utilizes the PCA class from the MultivariateStats library to perform the decomposition of time-series data into a set of underlying components. The resulting components can be used to identify trends, oscillations, and other patterns in the data.

For NLSA, the code utilizes the DiffusionMap class from the ManifoldLearning library. The approximated Laplace-Beltrami-operator is then computed.
Its spectrum contains important properties of the underlying dynamical system, which can be used to identify the dominant modes of variability in the data.

    
## Conclusion

One important property of both SSA and NLSA is that the reduced dimensions obtained from these methods are orthogonal to each other. This means that each dimension carries unique information about the dynamics of the system, with minimal redundancy between dimensions.
This orthogonality allows for a spread of information across different time scales, and can provide insight into the underlying mechanisms driving the observed patterns in the data.
It additionally constrains the information quality retrieved by the individual modes.
Here, the resulting harmonics of the fundamental frequency of seasonal behavior/trends are identified among the timescales to reconstruct the seasonality.
This seasonal information is then used to compare different vegetation related variables such as Gross Primary Production with environmental variables
among different locations and ecosystems.


In summary, the use of SSA and NLSA can help us better understand the complex dynamics of the earth system by identifying patterns and trends in time-series data. Phase and amplitude related insights about the seasonal behavior of different variables are gained by collecting the correct orthogonal modes that distinguish the seasonal behavior from residual information.

## References

- Ghil, M., Allen, M.R., Dettinger, M.D., Ide, K., Kondrashov, D., Mann, M.E., Robertson, A.W., Saunders, A., Tian, Y., Varadi, F., and Yiou, P. (2002). Advanced spectral methods for climatic time series. Reviews of Geophysics, 40(1), pp. 1003. https://doi.org/10.1029/2000RG000092
- Vautard, R., Yiou, P., and Ghil, M. (1992). Singular-spectrum analysis: A toolkit for short, noisy chaotic signals. Physica D: Nonlinear Phenomena, 58(1-4), pp. 95-126. https://doi.org/10.1016/0167-2789(92)90103-T
- Nadal, M.E., Harlim, J., and Majda, A.J. (2015). Nonlinear Laplacian spectral analysis for time series with intermittency and low-frequency variability. Proceedings of the National Academy of Sciences, 112(13), pp. 3918-3923. https://doi.org/10.1073/pnas.1419266112
- Coifman, R.R. and Lafon, S. (2006). Diffusion maps. Applied and Computational Harmonic Analysis, 21(1), pp. 5-30. https://doi.org/10.1016/j.acha.2006.04.006
Lahiri, S.N., Santhanam, M.S., and Chattopadhyay, J. (2015). Singular spectrum analysis and its variants: A review. Wiley Interdisciplinary Reviews: Computational Statistics, 7(4), pp. 268-283. https://doi.org/10.1002/wics.1351
- MultivariateStats.jl: Provides a collection of multivariate statistical methods for data analysis. The PCA class is part of this library. Reference: https://multivariatestatsjl.readthedocs.io/en/latest/
- ManifoldLearning.jl: Provides a collection of manifold learning methods for high-dimensional data analysis. The DiffusionMap class is part of this library. Reference: https://manifoldlearningjl.readthedocs.io/en/latest/
