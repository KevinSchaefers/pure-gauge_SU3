function [accepted, exp_deltaH, r] = Metropolis_step(H_old, H_new)
% METROPOLIS_STEP performs the accept/reject step inside the 
% Hybrid Monte Carlo (HMC) algorithm.
% The new configuration is accepted with probability min(1,exp(H_old-H_new))
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: U = METROPOLIS_STEP(H_old,H_new)
%--------------------------------------------------------------------------
% input:    H_old - H(p_0,q_0), the value of the Hamiltonian at the
%           beginning of the trajectory
%           H_new - H(p_1,q_1), the value of the Hamiltonian at the end of 
%           the trajectory
% output :  accepted - 1 if step is accepted, 0 otherwise
%           exp_deltaH - value of exp(H_old - H_new)
%           r - the random number used in the acceptance step. r=0 if no
%           random number has been used, i.e., if exp_deltaH>=1.
%--------------------------------------------------------------------------
    
    accepted = 0;       % step is not yet accepted
    r = 0;              % random number (only needed in case energy has increased)
    
    exp_deltaH = exp(H_old - H_new); 
    
    if exp_deltaH >= 1
        % energy has decreased
        accepted = 1;
    else
        % energy has increased -> Monte Carlo step
        r = rand;
        if r < exp_deltaH
            accepted = 1;
        end
    end
end