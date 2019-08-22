# opticalesr-rates

Python functions and example notebooks for simulating optical ESR experiments

## FineStructureHam()

This class generates the functions necessary to calculate the ground and excited states of NV- and isomorphic systems. The form of the operators, and the notation is based on the equations "Quantum Optics with Single Spins" - Lee Bassett https://arxiv.org/abs/1908.05566.

$ \sum_{\forall i}{x_i^{2}} $

The GetHelp() and GetHamHelp() functions in this class will provide short text outputs with various parameter definitions.

Calculating an optical spectrum with this class has five steps. 
1. Create a dictionary with the parameters for the Hamiltonian (see above for definitions). 
2. Update the class parameters with this using the .SetParams() function. 
3. Build the class Hamiltonian with .BuildHam(). 
4. Calculate the eigenvalues and eigenstates with .CalcEvals(). 
5. Calculate the spectrum with .CalcSpectrum()

There are several checks in the class that will stop you skipping a step. In this version, only selection rules for linear polarization are implemented. 
