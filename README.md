# dynamicblade
This MATLAB script numerically simulates the dynamics of a flexible blade under oscillatory flow.  The blade is modeled as an inextensible, linear elastic beam undergoing finite (i.e., nonlinear deformation).  The fluid forces are modeled based on the well-known Morrison force formulation.

This is a finite difference solver that is 2nd order accurate in space. The blade-normal force balance is solved explicitly to yield the bending angle, theta, while the blade-tangential force balance is solved implicitly to yield tension.

For further details, see:
1. M Luhar and HM Nepf (2016) Wave-induced dynamics of flexible blades, Journal of Fluids and Structures, 61: 20-41
2. M Luhar (2012) Analytical and experimental studies of plant-flow interaction at multiple scales, PhD thesis, MIT

The publications above used PIV-measured velocity fields to force the blade, here we use a simple sinusoidal velocity field. 

---
List of Scripts and Functions

1) dynamicblade.m:  This is the central script.  The dimensionless inputs (Cauchy Number; Buoyancy Parameter; Keulegan-Carptenter Number; Drag, Added Mass, Friction Coefficients etc.) are specific in this file, as are the discretization parameters (grid size, time step).  The script updates blade postures and plots them every few time steps. Blade posture, tension, and forces are also saved at each time step.  

2) fdmatrix.m: This generates the finite differencing matrices, which are used in dynamicblade.m

3) fdcoefs.m: This generate finite difference coefficients, which are used in fdmatrix.m

---
Developed by M Luhar* and HM Nepf+
* Department of Aerospace and Mechanical Engineering, USC
+ Department of Civil and Environmental Engineering, MIT

Written by M Luhar (luhar@usc.edu, mluhar@alum.mit.edu, mluhar@cantab.net).

SEE ATTACHED LICENSING AGREEMENT
