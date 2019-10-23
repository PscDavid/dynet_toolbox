# dynet_toolbox

The complete collection of functions and scripts for estimating and simulating
time-varying Multivariate Autoregressive processes (tv-MVAR)
by means of Kalman filtering and the Self-Tuning Optimized Kalman filter (STOK)

The toolbox includes:
- One demo (please refer to the file dynet_demo01.m for a brief tutorial)
- Four folders:
    - 'statespace'
        - dynet_SSM_KF.m implements the Kalman filter for state-space modeling of
        physiological time series
        - dynet_SSM_STOK.m implements the STOK algorithm
        filter with self-tuning memory and least-squares reconstruction

    - 'connectivity'
        - dynet_ar2pdc.m estimates the tv PDC from tv-AR coefficients
        - dynet_connplot.m displays connectivity matrices (function of time and
          frequency) for each combination of signals
        - dynet_parpsd.m estimates the AR coefficients in the frequency domain
        and the parametric power spectral density of the input signals

    - 'simulation'
        - dynet_sim.m is the simulation framework for tv-MVAR generated
        surrogate time series
        - review.m displays the 1) structural adjacency matrix, 2) the
        functional adjacency matrix, 3)surrogate time-series in the time domain,
        4) the power spectral density of surrogate time-series

    - 'utilities' contains all the invoked functions to let all the above listed
    functions/scripts to properly work
