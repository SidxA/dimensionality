# Code Descriptions

This README provides an overview of the content in each of the Julia files included in this repository.

## 01_basic_functionality.jl

This file contains the implementation of basic functionality in Julia. It includes essential functions and utilities that serve as the foundation for other scripts in this project.

## 02_dimensionreduction_structs.jl

The `02_dimensionreduction_structs.jl` file defines custom data structures and types required for dimension reduction algorithms. It provides the necessary abstractions and data representations for efficient data processing and analysis.

## 03_data_preprocessing.jl

The `03_data_preprocessing.jl` file focuses on data preprocessing techniques. It contains functions for cleaning, transforming, and normalizing raw data before it is used in subsequent analysis steps.

## 04_RUN_dimensionreduction.jl

In the `04_RUN_dimensionreduction.jl` file, you will find the main script for executing dimension reduction algorithms. It combines the functionality from previous files and performs dimensionality reduction on the input data, generating reduced-dimensional representations.

## 05_processing.jl

The `05_processing.jl` file contains additional data processing routines that are applied after dimension reduction. This includes feature extraction, outlier detection, or other post-processing steps that are specific to the project's requirements.

## 06_load_figuredata.jl

The `06_load_figuredata.jl` file focuses on loading and preparing data specifically required for generating figures or visualizations. It handles data loading, preprocessing, and formatting to facilitate the creation of informative and visually appealing figures.

## 07_fig1_motivation.jl

In the `07_fig1_motivation.jl` file, you will find the code for generating Figure 1. This figure illustrates the motivation or background information related to the project. It may include visualizations, charts, or other graphical representations.

## 08_fig2_results.jl

The `08_fig2_results.jl` file contains the code for generating Figure 2. This figure presents the results of the analysis performed in the project. It may visualize the performance, accuracy, or any other relevant metrics related to the dimension reduction algorithms.

## 09_fig3_heatmaps.jl

The `09_fig3_heatmaps.jl` file includes the code for generating Figure 3. This figure focuses on presenting heatmaps or similar visualizations, showcasing specific patterns, distributions, or correlations discovered during the analysis.

## 99_old_analysis.jl

The `99_old_analysis.jl` file contains old or deprecated code that was used in previous iterations of the project. It is kept for historical purposes but is no longer actively used in the current implementation.

Please refer to the individual Julia files for more detailed comments and code explanations.


# data

- time series main information loading by `load("/net/scratch/lschulz/data/time_series.jld2")`, holding
    - `data_raw, data_f3, data_f4, data_f6` raw and lowpass-filtered time series
    - `flags` quality flags for selected variables
    - `flag_variables` variables for wjhich the flags are selected
    - `spotslist` 18 fluxnet description names marking the measurement sites
    - `IGBP_list` corresponding 18 classified ecosystems
    - `variables_names` chosen names for the 16 selected variables
    - `variables_original_names` original fluxnet variable names
- the main analysis of the individual time series is performed using the `local_parameters` function that returns a dirty list of the extracted variables of interest
- the computed seasonal cycles are loaded by `load("/net/scratch/lschulz/data/seasonal_cycle.jld2")`, holding
    - `ssa_h_x` : 18x16 matrix of retrieved number of harmonics for SSA
    - `nlsa_h_x` : 18x16 matrix of retrieved number of harmonics for NLSA
    - `ssa_trends_x` : Nx18x16 tensor of the calculated seasonal cycles for SSA
    - `nlsa_trends_x` : Nx18x16 tensor of the calculated seasonal cycles for NLSA
- with `x` out of `[raw,3,4,6]` for the corresponding (un)-filtered time series
- 
# dimensionality reduction methods

- run `export JULIA_NUM_THREADS=1` in bash and start julia with desired number of cores to specify CPU usage
- the dimensionality reduction is performed 1 timeseries per core using the same W values to re-write the big matrices to save RAM
- each individual calculation outputs a SSA and NLSA jld2 file with the modes,etc. inside
- main ingredients contained in `fundamentals.jl`
- just run `/opt/julia-1.9.0/bin/julia --threads 64 iterate_blockdata.jl` for using 64 cores
- with the specified parameters changed inside
  - `savedirname = "/net/scratch/lschulz/data/time_series.jld2"` the file containing the `[N,spots,variables]` different values for different filters
  - `wholedata = SharedArray{Float32}(load(savedirname)["data_x"])`with `x` in `[raw,f3,f4,f6]` for corresponding dataset
  - this file contains 18 spots and 16 variables selected for beeing the longest measurement periods of the variables of interest
  - `outdir = "/net/scratch/lschulz/fluxfullset_midwithnee_lowpass4/"` the directory getting created that the individual files are saved in
  - `the desired different `W` need to be specified in the list at the parallel loop at the end of the document: thsi is run at 7

# analysis

- inside the `analysis.jl` file, individual functions perform the calculations and plotting the create the figures
- it is important if julia is run externally to run it by `GKSwstype=nul /opt/julia-1.9.0/bin/julia` in order to prevent it from trying to display images directly
- the created figures can be saved into specified directory `dir` and looked at there
- for submission change julia Makie fonts by `set_theme!(fonts=(
    regular="Latin Modern Roman",
    bold = "Latin Modern Roman Bold",
    italic = "Latin Modern Roman Italic",
    bold_italic = "Latin Modern Roman Bold Italic",))`


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
