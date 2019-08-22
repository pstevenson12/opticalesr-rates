# opticalesr-rates

Python functions and example notebooks for simulating optical ESR experiments

## FineStructureHam()

This class generates the functions necessary to calculate the ground and excited states of NV- and isomorphic systems. The form of the operators, and the notation is based on the equations "Quantum Optics with Single Spins" - Lee Bassett https://arxiv.org/abs/1908.05566.

The parameters used here are:\
D<sub>{gs}</sub>,D<sub>{es}</sub>: Ground and excited state zero-field splittings. Generally a combination of spin-orbit and spin-spin interactions\
g<sub>{gs}</sub>,g<sub>{es}</sub>: g-tensors for the ground and excited states. Should be a list with parallel and perpendicular terms\
&lambda;<sub>z:</sub> Axial spin orbit coupling in the excited state\
&Delta;<sub>1</sub>: Axial spin-spin term; does not cause mixing of spin states\
&Delta;<sub>2</sub>: Non-axial spin-spin term; does cause mixing of spin states, should be zero in inversion symmetric systems \
B: Magnetic field vector in Gauss \
&delta;: Strain interaction terms \
&alpha;: Strain angle; 0 is along axial dimension

All units are in GHz unless otherwise specified

The GetHelp() and GetHamHelp() functions in this class will provide short text outputs with various parameter definitions.

Calculating an optical spectrum with this class has five steps. 
1. Create a dictionary with the parameters for the Hamiltonian (see above for definitions). 
2. Update the class parameters with this using the .SetParams() function. 
3. Build the class Hamiltonian with .BuildHam(). 
4. Calculate the eigenvalues and eigenstates with .CalcEvals(). 
5. Calculate the spectrum with .CalcSpectrum()

There are several checks in the class that will stop you skipping a step. In this version, only selection rules for linear polarization are implemented. 
