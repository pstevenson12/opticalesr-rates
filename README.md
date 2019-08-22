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

## FiveLevelModel()

This class generates the functions required to numerically integrate the rate equations for a five-level system, with a specific view towards simulating color-center like system. (Two ground states, two excited states and a shelving state are assumed in the form of the parameter inputs.)

The states are labelled and ordered: A<sub>G</sub>, B<sub>G</sub>, A<sub>E</sub>, B<sub>E</sub>, S

All rates are assumed zero by default - the use <b>must</b> input some parameters before the rest of the functions will allow themselves to be used. 

The rates incorporated into this version are: \
k<sub>r,A</sub>, k<sub>r,B</sub> - the radiative rates of A and B, respectively \
&phi; - branching ratio of the radiative decay between A and B (i.e. what fraction A-A, A-B)
k<sub>up,A</sub>, k<sub>up,B</sub> - optical excitation rates of A (B) spin-conserving transition \
k<sub>up,AB</sub>, k<sub>up,BA</sub> - optical excitation rates of the non-spin-conserving transitions. This is related to the branching ratio, &phi;, but these two parameters do not have a simple relation here, because k<sub>up,AB</sub> also depends on the detuning of the optical excitation from a particular transition. \
Spin-lattice relaxation rate - ground state T<sub>1</sub> \
ISC - intersystem crossing rate from excited state to shelving state. Currently, both excited states are assumed to have the same ISC rate. This will be updated in a later version. \
ShG - Shelving state to ground state rate (inverse shelving state lifetime) \
&phi;<sub>s</sub> - branching ratio of shelving state relaxation between A<sub>G</sub>, B<sub>G</sub>. 0 is completely into A, 1 is equally into A, B.

Calculating the kinetic traces has X steps
1. Create a dictionary with the parameters for the ratematrix (see above for definitions). 
2. Update the class parameters with this using the .SetRateDict() function. 
3. Build the class Hamiltonian with .BuildRateMatrix(). 
4. Build a parameter dictionary and update the class parameters with this using the .SetSolverParameters() function. 
5. Use the ODE solver with .SolveRateEquations(). 

The output is stored as .Results in the class.
