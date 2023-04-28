# Time-Series Analysis of Climatic Data from the FLUXNET Network

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
