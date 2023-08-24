## Overview

This repository contains the `julia 1.9.0` code used for the publication LINK MISSING. 
Dimensionality reduction (time series decomposition) is used for the purpose of
accurately identifying the saisonal trend of the following variables

SWC (soil water content), TS (temperature soil), TA (temperature air), SW_IN (shortwave radiation incoming), NEE (net ecosystem exchange), RECO (respiration of ecosystem), GPP (gross primary production), 

of FLUXNET daily measurement series within 01.01.2005-31.12.2020 of sites 

CH-Dav, DE-Tha, FI-Hyy, IT-Lav, RU-Fyo, BE-Vie, CH-Lae, DE-Hai, DK-Sor.

## Methods

Using Diffusion Maps and Singular Spectrum Analysis, the feasability of extracting a seasonal cycle of sufficient accuracy - e.g. harmonic order of at least 2 - is ivestigated.
The effects of noise are investigated by a FFT-based lowpass filter.

## Repository structure

- main functionalities such as required functions and imported packages are in `01_basic_functionality.jl`
- the struct used to perform multiple calculations in parallel are in `02_dimensionreduction_structs.jl`
- lowpass filter for the datasets are in `03_data_preprocessing.jl`
- executions of the calculations over the datasets are in `03_04_RUN_dimensionreduction.jl` 
- processing of the results such as identifying harmonic modes and combining them are in `05_processing.jl`, same as, saving the results
- loading the results and preparing them for plotting is done in `06_load_figuredata.jl`
- `07_fig1_motivation.jl` holds the code for the creation of the motivation figure
- `08_fig2_results.jl` holds the code for the creation of the individual example/scenario figures
- `09_fig1_heatmaps.jl` holds the code for the creation of the showcasing overview figure
- `10_workflow_figures.jl` holds the code for the creation of parts of the workflow schema

## Code quickstart: analyzing a single time series raw and preprocessed

Everything is tailored for Float32 for FLUXNET precision and computation speed.
This rather lengthy example consists of the subfunctions with main functionality for clarity.
The following code contains all reasoning steps required for the understanding of the algorithmic structure.
All numerical results of the publication can be obtained by employing this procedure sequentially on the individual timeseries.
The material in the repository handles the dataset at once and tailors specific figures.

Dimensionality reduction
```julia
# load prerequisites
include("REPOSITORY-POSITION/dimensionality/01_basic_functionality.jl")
# example GPP series
time_series = load("example.jld2")["example"]
# number of measurements ( 5114 daily measurements ~ 14 a length)
N = length(time_series)
# embedding length (longest full year period: 7a)
W = Int(floor(365.25 .* 7))
# arbitrary number of extracted modes, here set to 48
k = 48
# creating a fixed struct for dimensionsality reduction holding the big matrices of the calculations
include("/net/home/lschulz/dimensionality/02_dimensionreduction_structs.jl")
d = matrixholder(N, W, k)
# Set the number of samples per year
sampleyear = 365
# execute computation 1 (example without preprocessing)
  # Process the data (centralization), inject into the struct
  put_data_raw(d, time_series)
  # Compute epsilon
  put_epsilon(d)
  # perform nonlinear dimensionality reduction (NLSA), yielding EOF basis set
  put_EOF_diff(d)
  # project the data on the EOF basis set
  calculate(d)
  # Save the nonlinear results (name of saving directory)
  save_results(d, "unfiltered_nlsa_results")
  # Perform singular spectrum analysis (SSA) using the same struct
  put_EOF_ssa(d)
  calculate(d)
  # Save the linear results (name of saving directory)
  save_results(d, "unfiltered_ssa_results")
# preprocess the example series by lowpass filter with a cutoff frequency off 6/a
include("/net/home/lschulz/dimensionality/03_data_preprocessing.jl")
cutoff_frequency = 6
time_series_lowpassed = lowpass_filter(time_series, cutoff_frequency)
# execute computation 2
  # Process the data (centralization), inject into the struct
  put_data_raw(d, time_series_lowpassed)
  # Compute epsilon
  put_epsilon(d)
  # perform nonlinear dimensionality reduction (NLSA), yielding EOF basis set
  put_EOF_diff(d)
  # project the data on the EOF basis set
  calculate(d)
  # Save the nonlinear results (name of saving directory)
  save_results(d, "filtered_nlsa_results")
  # Perform singular spectrum analysis (SSA) using the same struct
  put_EOF_ssa(d)
  calculate(d)
  # Save the linear results (name of saving directory)
  save_results(d, "filtered_ssa_results")
```

Processing of the combined results
```julia
# load in the data from the computations and investigate the EOS basis set
include("/net/home/lschulz/dimensionality/05_processing.jl")
p_unfiltered = extract_from_directory("unfiltered_ssa_results","unfiltered_nlsa_results")
p_filtered = extract_from_directory("filtered_ssa_results","filtered_nlsa_results")
# p holds         spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa
# from these, we only require the identified frequencies of individual modes to decide wether they belong to the harmonic strcture of the fundamental frequency
```
find out wether there is at least second harmonic order of seasonal cycle present 
(with a fixed threshold)
```julia

threshold = 0.15

# unfiltered data
freq_list_ssa = p_unfiltered[30]
freq_list_nlsa = p_unfiltered[31]
#two modes per frequency required: booleans - wether they are present
ssa_firstorder = length(findall(x->abs(x-1.0)<=threshold,freq_list_ssa)) >= 2
ssa_secondorder = length(findall(x->abs(x-2.0)<=threshold,freq_list_ssa)) >= 2
nlsa_firstorder = length(findall(x->abs(x-1.0)<=threshold,freq_list_nlsa)) >= 2
nlsa_secondorder = length(findall(x->abs(x-2.0)<=threshold,freq_list_nlsa)) >= 2

# filtered data
freq_list_ssa = p_filtered[30]
freq_list_nlsa = p_filtered[31]
#two modes per frequency required: booleans - wether they are present
ssa_firstorder = length(findall(x->abs(x-1.0)<=threshold,freq_list_ssa)) >= 2
ssa_secondorder = length(findall(x->abs(x-2.0)<=threshold,freq_list_ssa)) >= 2
nlsa_firstorder = length(findall(x->abs(x-1.0)<=threshold,freq_list_nlsa)) >= 2
nlsa_secondorder = length(findall(x->abs(x-2.0)<=threshold,freq_list_nlsa)) >= 2
```
## References

- Ghil, M., Allen, M.R., Dettinger, M.D., Ide, K., Kondrashov, D., Mann, M.E., Robertson, A.W., Saunders, A., Tian, Y., Varadi, F., and Yiou, P. (2002). Advanced spectral methods for climatic time series. Reviews of Geophysics, 40(1), pp. 1003. https://doi.org/10.1029/2000RG000092
- Vautard, R., Yiou, P., and Ghil, M. (1992). Singular-spectrum analysis: A toolkit for short, noisy chaotic signals. Physica D: Nonlinear Phenomena, 58(1-4), pp. 95-126. https://doi.org/10.1016/0167-2789(92)90103-T
- Nadal, M.E., Harlim, J., and Majda, A.J. (2015). Nonlinear Laplacian spectral analysis for time series with intermittency and low-frequency variability. Proceedings of the National Academy of Sciences, 112(13), pp. 3918-3923. https://doi.org/10.1073/pnas.1419266112
- Coifman, R.R. and Lafon, S. (2006). Diffusion maps. Applied and Computational Harmonic Analysis, 21(1), pp. 5-30. https://doi.org/10.1016/j.acha.2006.04.006
Lahiri, S.N., Santhanam, M.S., and Chattopadhyay, J. (2015). Singular spectrum analysis and its variants: A review. Wiley Interdisciplinary Reviews: Computational Statistics, 7(4), pp. 268-283. https://doi.org/10.1002/wics.1351
- MultivariateStats.jl: Provides a collection of multivariate statistical methods for data analysis. The PCA class is part of this library. Reference: https://multivariatestatsjl.readthedocs.io/en/latest/
- ManifoldLearning.jl: Provides a collection of manifold learning methods for high-dimensional data analysis. The DiffusionMap class is part of this library. Reference: https://manifoldlearningjl.readthedocs.io/en/latest/
