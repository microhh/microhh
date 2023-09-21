# Point source dispersion

## Case description
Demonstration case for dispersion from point sources, and the option to have non-periodic lateral boundary conditions (BCs) for scalars. Scalar `s2` uses periodic BCs, and scalar `s1` open BCs where mass can be advected out of the domain.

### Mass conservation
With the doubly-periodic BCs used for `s2`, this scalar should conserve mass (initial + emitted mass), and with its open BCs, scalar `s1` should lose mass relative to `s2`. The Python script `check_mass_conservation.py` checks this. With a scalar emission of `1 kg s-1` the results are as expected: 

    In [1]: run check_mass_conservation.py  
    Time =     0, scalar "s1": correct mass = 0.00 kg,     integral field = 0.00 kg  
    Time =     0, scalar "s2": correct mass = 0.00 kg,     integral field = 0.00 kg  
    Time =  3600, scalar "s1": correct mass = 3600.00 kg,  integral field = 2486.29 kg  
    Time =  3600, scalar "s2": correct mass = 3600.00 kg,  integral field = 3600.39 kg  
    Time =  7200, scalar "s1": correct mass = 7200.00 kg,  integral field = 2664.75 kg  
    Time =  7200, scalar "s2": correct mass = 7200.00 kg,  integral field = 7200.38 kg  
    Time = 10800, scalar "s1": correct mass = 10800.00 kg, integral field = 2875.91 kg  
    Time = 10800, scalar "s2": correct mass = 10800.00 kg, integral field = 10800.37 kg

### Gaussian fit of plume
The Python script `fit_gaussian_curve.py` reads the vertically integrated `s1` values, time averages them, and fits a Gaussian curve at several distances from the plume. As you can see, LES produces nice Gaussian shaped plumes:

![gaussian_fit](https://github.com/julietbravo/microhh/assets/1067637/96d5c32e-5750-4621-b76a-66475d248427)
