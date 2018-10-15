Micro FEM
---------

This library provides finite element models for three types of systems on
a uniform 2D rectangular mesh. The purpose of the library is for easy 
integartion into strucutural optimization algorithms for cantilever like 
structures. 


The first system is for the modal analysis of a Mindlin plate, the second is 
for the modal analysis of a laminate plate structure containing a piezoelectric 
layer, and the third is a solver for poisson's equation that can be used to 
analyse many different physical phenomena. All three have a predetermined 
fixed boundary condition along the left edge of the rectangular mesh. 


To install the package, change directory to the folder containing  

    > pip install .


Examples of the functionality provided by this library are contained in the
examples folder.


Extensions
----------

- Create an object that encapsulates all the analysis from each type of
    solver.
- Avoid rank one numpy arrays and add asserts to check these arrays don't 
    exist.
- Add damping models to the plate and laminate FEM.
- Add the possibility of multiple piezoelectric patches to the laminate FEM.
- Add static solver with applied voltages and forces to plate and laminate FEM.
