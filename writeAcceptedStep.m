function writeAcceptedStep(result, loop, trajectory, H, S, eKin,...
                deltaH,exp_deltaH,r,meanPlaquette, ...
                accept, acceptance_rate)
% WRITEACCEPTEDSTEP writes the intermediate output into the result file
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: WRITEACCEPTEDSTEP(result, loop, trajectory, H, S, eKin, deltaH, 
%       exp_deltaH, r, meanPlaquette, accept, acceptance_rate)
%--------------------------------------------------------------------------
% input: result            - output file
%        loop              - number of the loop
%        trajectory        - length of accepted trajectories
%        H                 - new Hamiltonian
%        S                 - new action
%        eKin              - new kinetic energy
%        deltaH            - difference of old and new Hamiltonian
%        exp_deltaH        - exponential function of deltaH
%        r                 - random number
%        meanPlaquette     - mean value of the new plaquette
%        accept            - accepted? 0 = no, 1 = yes
%        acceptance rate   - cumulative acceptance rate
% output :  -
%--------------------------------------------------------------------------

fprintf(result,'%4d\t',loop);
fprintf(result,'%4d\t',trajectory);
fprintf(result,'%1.10e\t',H);
fprintf(result,'%1.10e\t',S);
fprintf(result,'%1.10e\t',eKin);
fprintf(result,'%1.10e\t',deltaH);
fprintf(result,'%1.10e\t',exp_deltaH);
fprintf(result,'%1.6f\t',r);
fprintf(result,'%1.10e\t',meanPlaquette);
fprintf(result,'%1d\t',accept);
fprintf(result,'%3.2f',acceptance_rate);
fprintf(result,'\n');