# Simulation of pure gauge fields in SU(3) on a two-dimensional lattice 

This MATLAB implementation enables the simulation of pure gauge fields in SU(3) on a two-dimensional lattice. 
The implementation aims at providing a comprehensive code that allows for an easy implementation and testing of new numerical integration schemes for the molecular dynamics step in the hybrid Monte Carlo algorithm. 

## Demands on the numerical integration schemes
As the links are elements of the matrix Lie group SU(3), the numerical integration scheme has to satisfy the closure property, i.e., the numerical approximations should also be situated on the Lie group manifold.
Furthermore, the integrator must be time-reversible and volume-preserving to satisfy the detailed balance condition.

## Description 
The current version of the code provides an implementation of the following numerical integration schemes:
- BAB, the velocity version of the Störmer-Verlet / leapfrog method, also known as the Strang splitting [https://doi.org/10.1137/0705041]
- ABABA, the position version of the second-order minimum-norm scheme [https://doi.org/10.1137/0916010 , https://doi.org/10.1016/S0010-4655(02)00754-3]
- BABAB, the velocity version of the second-order minimum-norm scheme [https://doi.org/10.1137/0916010 , https://doi.org/10.1016/S0010-4655(02)00754-3]
- BADAB, five-stage Hessian-free force-gradient integrator of convergence order four [https://doi.org/10.48550/arXiv.2403.10370]
- ABADABA, seven-stage Hessian-free force-gradient integrator of convergence order four [https://doi.org/10.48550/arXiv.2403.10370]
- OMF4 (BABABABABAB), eleven-stage fourth-order minimum norm scheme [https://doi.org/10.1016/S0010-4655(02)00754-3] 
- Yoshida, Yoshida's triple-jump composition scheme of order four using BAB as the base scheme [https://doi.org/10.1016/0375-9601(90)90092-3]
- Suzuki, Suzuki's fractals composition scheme of order four using BAB as the base scheme [https://doi.org/10.1016/0375-9601(90)90962-N]
- advancedComposition, advanced composition scheme of order six using BAB as the base scheme [https://doi.org/10.1090/S0025-5718-97-00873-9]
In general, the methods are using the matrix exponential (which can be efficiently implemented for SU(3) using the approach proposed in [https://doi.org/10.1140/epja/s10050-022-00816-5]). As this code has also been used to investigate the usage of the modified Cayley transform for SU(3) [https://doi.org/10.48550/arXiv.2406.11337], this code also allows to choose the modified Cayley transform instead of the matrix exponential. For integrators up to second order and higher-order methods based on composition techniques keep the desired convergence order. For direct decomposition algorithms (e.g. OMF4) and Hessian-free force-gradient integrators (e.g. BADAB), one observes an order reduction to convergence order two as these algorithms are derived based on the Baker-Campbell-Hausdorff formula.

## Usage
The main function is `Hybrid_Monte_Carlo.m`. This function applies the HMC algorithm fore pure gauge field simulations on a L x L lattice at specified inverse temperature beta. In the molecular dynamics step, it applies a volume-preserving and time-reversible integrator using a constant step size h. For the integrators, different local parameterizations (either the matrix exponential or the modified Cayley transform) can be used to map from the Lie algebra su(3) to the Lie group SU(3). All in all, the function computes a specified number of trajectories of length tau.

The function call reads 

U = Hybrid_Monte_Carlo(L, method_string, param_beta, h, local_param, runname, ntraj, ncheckpoints, tau,U_start)

where the output U will be the final link field. As input parameters, we have 

- L  : lattice size (2D lattice of size L x L)
- method_string : string describing the method used for the MD step. The following methods are implemented: ABABA, ABADABA, advancedComposition, BAB, BABAB, BADAB, OMF4, Suzuki, Yoshida.
- param_beta : inverse temperatur beta
- h : step size used in the MD step
- local_param : local parameterization used to compute links in SU(3). The following local parameterizations are implemented: exponential_map (matrix exponential for SU(3)), caymod (modified Cayley transform). 
- runname : name of the run to avoid overwriting files. The results are then stored in results/runname/
- ncheckpoints : after ncheckpoints trajectories, an intermediate result will be stored.
- tau : trajectory length. 
- U_start : optional argument to start with a link field from a previous run. Otherwise, the simulation starts with a random link field.

After each trajectory, a summary of the trajectory is written into the output file. Furthermore, any checkpoint creates a .mat file that also contains the data to perform some further analysis (e.g. to compute the variance of Delta H in order to evaluate the performance of the integrators)

As an example, there is a file `test_script.m` performing a simulation.

## Author
The implementation of this release has been done by Kevin Schäfers. 
If a general bug in the code is discovered, please send a report to schaefers@math.uni-wuppertal.de 
Furthermore, any kind of suggestions and contributions to this project is highly welcomed.

## Acknowledgements
The present Matlab implementation started as a generalization of a Matlab implementation to perform pure gauge field simulations on a two-dimensional lattice in SU(2). I am grateful that Michèle Wandelt shared her Matlab implementation for SU(2) gauge field simulations that served as a starting point for this project.

## License
The software may be used under the terms of the GNU General Public Licence (GPL).
