function writeHeader(result, method, seed, h, traj, local_param)
% WRITEHEADER creates the header of the output file
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: WRITEHEADER(result, method, seed, h, traj, local_param)
%--------------------------------------------------------------------------
% input:    result - output file 
%           method - used integration method
%           seed - seed of the random number generators
%           h - step size
%           traj - trajectory length
%           local_param - the local parameterization used to map from the Lie
%               algebra to the Lie group
% output :  -
%--------------------------------------------------------------------------

global beta;


fprintf(result, '%s %s :\n', '%', method);
fprintf(result,'%s%f\t%s%f\t%s%f\t%s%d\t%s%g\n','% beta = ', beta, 'h = ', h, ...
        'trajectory = ', traj, 'seed = ', seed);
fprintf(result,'%s%s\n%s%s\n', '% local parameterization = ', local_param);
fprintf(result,'%s\n','%  1: total number of trajectories');
fprintf(result,'%s\n','%  2: number of accepted trajectories');
fprintf(result,'%s\n','%  3: total energy');
fprintf(result,'%s\n','%  4: action');
fprintf(result,'%s\n','%  5: kinetic energy');
fprintf(result,'%s\n','%  6: difference between old and new (total) energy');
fprintf(result,'%s\n','%  7: exp(difference between old and new (total) energy)');
fprintf(result,'%s\n','%  8: random number');
fprintf(result,'%s\n','%  9: mean plaquette');
fprintf(result,'%s\n','% 10: accepted (0=false, 1=true)');
fprintf(result,'%s\n','% 11: acceptance rate');