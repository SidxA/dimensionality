# Overview

This repository contains the julia code used for the publication LINK.
dimensionality reduction / time series decomposition
with the aim of accurately identifying the saisonal trend of the following variables

SWC, TS, TA, SW_IN, NEE, RECO, GPP,

of 01.01.2005-31.12.2020 FLUXNET data of sites

CH-Dav, DE-Tha, FI-Hyy, IT-Lav, RU-Fyo, BE-Vie, CH-Lae, DE-Hai, DK-Sor

# Methods

Using Diffusion Maps and Singular Spectrum Analysis, the feasability of extracting a seasonal cycle of sufficient accuracy - e.g. harmonic order of at least 2 - is ivestigated.
The effects of noise are investigated by a FFT-based lowpass filter.

# Repository structure

- main functionalities such as required functions and imported packages are in `01_basic_functionality.jl`
- the struct used to perform multiple calculations in parallel are in `02_dimensionreduction_structs.jl`
- lowpass filter for the datasets are in `03_data_preprocessing.jl`
- executions of the calculations over the datasets are in `03_04_RUN_dimensionreduction.jl` 
- processing of the results such as identifying harmonic modes and combining them are in `05_processing.jl`, same as, saving the results
- loading the results and preparing them for plotting is done in `06_load_figuredata.jl`
- `07_fig1_motivation.jl' holds the code for the creation of the motivation figure
- `08_fig2_results.jl' holds the code for the creation of the individual example/scenario figures
- `09_fig1_heatmaps.jl' holds the code for the creation of the showcasing overview figure
- `10_workflow_figures.jl' holds the code for the creation of parts of the workflow schema

# Code quickstart: full functionality of analyzing a single time series
Everything is tailored for Float32 for FLUXNET precision and computation speed.

    # load prerequisites
    include("REPOSITORY-POSITION/dimensionality/01_basic_functionality.jl")
    # example series
    time_series = load("example.jld2")["example"]
    # number of measurements ( 5114 daily measurements ~ 14 a length)
    N = length(tme_series)
    # embedding length (longest full year period: 7a)
    W = Int(floor(365.25 .* 7))
    # arbitrary number of extracted modes, here set to 48
    k = 48
    # creating a fixed struct for dimensionsality reduction holding the big matrices of the calculations
    # Set the number of samples per year
    sampleyear = 365
    d = matrixholder(N, W, k)
    # execute computation 1 (example without preprocessing)
      # Process the data (centralization)
      put_data_raw(d, time_series)
      # Compute epsilon
      put_epsilon(d)
      # perform nonlinear dimensionality reduction (NLSA), yielding EOF basis set
      put_EOF_diff(d)
      # project the data on the EOF basis set
      calculate(d)
      # Save the nonlinear results
      saving_directory = "unfiltered_nlsa_results"
      save_results(d, saving_directory)
      # Perform singular spectrum analysis (SSA) using the same struct
      put_EOF_ssa(d)
      calculate(d)
      # Save the linear results
      saving_directory = "unfiltered_ssa_results"
      save_results(d, saving_directory)
    # preprocess the example series by lowpass filter with a cutoff frequency in 1/a
    # 6 / year is a reasonable cutoff frequency for information content about the seasonal cycle
    cutoff_frequency = 6
    time_series_lowpassed = lowpass_filter(time_series, cutoff_frequency)
    # execute computation 2
    # save results 2
    # processing of the combined results
    # creation of the results figure

# References

- Ghil, M., Allen, M.R., Dettinger, M.D., Ide, K., Kondrashov, D., Mann, M.E., Robertson, A.W., Saunders, A., Tian, Y., Varadi, F., and Yiou, P. (2002). Advanced spectral methods for climatic time series. Reviews of Geophysics, 40(1), pp. 1003. https://doi.org/10.1029/2000RG000092
- Vautard, R., Yiou, P., and Ghil, M. (1992). Singular-spectrum analysis: A toolkit for short, noisy chaotic signals. Physica D: Nonlinear Phenomena, 58(1-4), pp. 95-126. https://doi.org/10.1016/0167-2789(92)90103-T
- Nadal, M.E., Harlim, J., and Majda, A.J. (2015). Nonlinear Laplacian spectral analysis for time series with intermittency and low-frequency variability. Proceedings of the National Academy of Sciences, 112(13), pp. 3918-3923. https://doi.org/10.1073/pnas.1419266112
- Coifman, R.R. and Lafon, S. (2006). Diffusion maps. Applied and Computational Harmonic Analysis, 21(1), pp. 5-30. https://doi.org/10.1016/j.acha.2006.04.006
Lahiri, S.N., Santhanam, M.S., and Chattopadhyay, J. (2015). Singular spectrum analysis and its variants: A review. Wiley Interdisciplinary Reviews: Computational Statistics, 7(4), pp. 268-283. https://doi.org/10.1002/wics.1351
- MultivariateStats.jl: Provides a collection of multivariate statistical methods for data analysis. The PCA class is part of this library. Reference: https://multivariatestatsjl.readthedocs.io/en/latest/
- ManifoldLearning.jl: Provides a collection of manifold learning methods for high-dimensional data analysis. The DiffusionMap class is part of this library. Reference: https://manifoldlearningjl.readthedocs.io/en/latest/
