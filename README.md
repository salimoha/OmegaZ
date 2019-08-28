# $$\delta$$-DOGS($$\Omega_Z$$)

This repo contains code accompaning the paper, 	[Design of IMEXRK time integration schemes via Delaunay-based derivative-free optimization with nonconvex constraints and grid-based acceleration)]. It includes code for running the several benchmark optimization problems with nonconvex objective function and nonconvex (or even disconnected) feasible domain defined with a set of nonconvex constraint functions, including the derivation of a new, low-storage, high-accuracy, Implicit/Explicit Runge-Kutta (IMEXRK) time integration scheme.

### Dependencies
This code requires the following:
* [SNOPT](https://ccom.ucsd.edu/~optimizers/solvers/snopt/) if the target value is unkown (constant-K approach).
* MATLAB software
* [MATLAB Optimization Toolbox](https://www.mathworks.com/products/optimization.html)


### Usage
To run the code, first go to the main repository and in the MATLAB command line run:

```
>> addTopath
```
Afterwards, all the dependances are going to be installed in your MATLAB path and you can run the test examples. 

### Contact
To ask questions or report issues, please open an issue on the [issues tracker](https://github.com/salimoha/deltaDOGS_OmegaZ/issues).



# deltaDOGS Omega Z

This Delaunay-based Derivative-free Optimization via Global Surrogates family of algorithms, and is used to identify a new, low-storage, high-accuracy, Implicit/Explicit Runge-Kutta (IMEXRK) time integration scheme for the stiff ODEs arising in high performance computing applications, like the simulation of turbulence. 


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

