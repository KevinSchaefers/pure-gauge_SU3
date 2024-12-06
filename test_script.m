% example script
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------

clear; clc;
addpath(genpath(pwd))

L = 32; % lattice of size 32 x 32
param_beta = 2.0; % inverse temperature
tau = 2.0; % trajectory length
ntraj = 1000; ncheckpoints=ntraj; % 1000 trajectories and only one checkpoint file at the end of the simulation
local_param = 'exponential_map'; % exponential map as the local parameterization 
% (alternatively, one could choose the modified Cayley transform by sett local_param = 'caymod')
h = 0.1; % step size (make sure tau/h is an integer)
method = 'BADAB'; % in this example, we use the five-stage fourth-order Hessian-free force-gradient integrator BADAB.
runname = 'testrun'; % the output file will contain 'testrun' in the name

U = Hybrid_Monte_Carlo(L, method, param_beta, h, local_param, runname, ntraj, ncheckpoints, tau); % we start from a random link field (no U_start)