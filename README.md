# TOV solver 
by Marina Berbel Palomeque

This program computes the solution to the TOV equations using one of the different equation of state (EoS) model:
* a piecewise polytropic (PP) EoS (Read et al 2009, Phys. Rev. D, 79)
* a tabulated EoS 
* a thermodynamically adaptive slope piecewise polytropic (T-ASPP) EoS (Berbel and Serna, Phys. Rev. D, 2023)

The TOV equations can be found in any book covering stellar structure, such as Neutron Stars 1 Equation of state and structure, P.Haensel, A.Y. Potekhin and D. G. Yakovlev, Springer, 2007

For details about the T-ASPP EoS see the reference * Berbel and Serna, Phys. Rev. D, 2023 *

NOTE: you are free to use this code, but you should cite the following papers if you make use of it for your publications:

1) Berbel and Serna, Phys. Rev. D, 2023

(citation to paper 1 is required only if you use the T-ASPP EOS)

Copyright 2023, Marina Berbel, All rights reserved

## COMPILATION INSTRUCTIONS

There is a Makefile distributed with the source code, simply change the name of the C compiler, if you use a different one, and type make.

## USING THE CODE

The executable is named solver and runs without any parameters. The details of the EoS are specified in a par file, eos.par.
The instructions for the type of computation desired are specified in tov.par. Examples of both type of par files are included in the folder.

It can compute the properties of a single star for a given central density: mass, radius, metric and tidal deformability.
It can also construct a sequence of stars for a range of central densities.

It includes the possibility of computing the precision of the neutron star models. The product of the baryon chemical potential mu_b with the exponential of the metric phi must be constant along the star radius for the star to be in perfect hydrostatic equilibrium. The program computes the maximum relative deviation of m_b*exp(phi) from a constant value. Tests show that at most, it is of a few fractions of percent. For PP and T-ASPP EoSs the models are very accurate, while a larger precision for tabulated EoSs would require a different interpolation technique.

If using a tabulated EoS, the program expects the following columns:
1-rest-mass density (g/cm3) 2-Pressure (g/cm3) 3-internal energy

The code saves some useful information in the following files:

* XXprofile_Y.txt: this file contains the radial profile of a
single star computed with EoS model XX and central density Y if saveRadialProfile is set True in 
the tov.par.  If set to False, mass, radius and tidal deformability will be printed in the console.

* XXstarSeq.txt: this file contains the central density, mass, radius and tidal deformability of the 
stars in the sequence computed.

* A temporal temp_file.txt is created during execution if checkAccuracy is set to True in the tov.par.
The program deletes it after execution.

