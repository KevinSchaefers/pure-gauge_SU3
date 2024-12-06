function [U,iP] = advancedComposition(U, iP, local_param, h, traj)
% advancedComposition computes one trajectory of the molecular dynamics 
% (MD) step using an advanced composition scheme of convergence order
% p=6 using BAB as the base scheme. The composition scheme has been
% introduced in <a
% href="matlab:web('https://doi.org/10.1090/S0025-5718-97-00873-9')">[Kahan and Li 1997]</a>.
% and uses the composition weights gamma_1 = 0.78451361047755726382,  
%   gamma_2 = 0.23557321335935813368, gamma_3 = -1.1776799841788710069,
%   gamma_4 = 1.3151863206839112189, and gamma_5 = gamma_3,
%   gamma_6 = gamma_2, gamma_7 = gamma_1 
% so that the overall composition method is a splitting method with 15
% stages.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: [U,iP] = advancedComposition(U, iP, local_param, h, traj)
%--------------------------------------------------------------------------
% input:    U               : links
%           iP              : conjugate momenta
%           local_param     : local parameterization to map into the
%                               Lie group
%           h               : step size
%           traj            : trajectory length
% output:   U               : updated links
%           iP              : updated momenta
%--------------------------------------------------------------------------
    
    % symmetric coefficients of the composition
    gamma_1 = 0.78451361047755726382;  
    gamma_2 = 0.23557321335935813368;
    gamma_3 = -1.1776799841788710069;
    gamma_4 = 1.3151863206839112189;
    
    % first step (momentum update)
    iP = iP + gamma_1/2 * h * g(U);

    nit = abs(traj/h); % nit time steps of the integrator 
    for it=1:nit
        % link update
        omega = get_omega(local_param,gamma_1*h*iP); % omega_1 = gamma_1 h dPsi_{omega}^{-1}(iP)
        psi_iP = feval(local_param, omega); % applying the local parameterization 
        % to the Lie algebra element omega, resulting in Lie group element psi_iP
        U = matMultField(psi_iP,U); % multiplication of Lie group elements
        
        iP = iP + (gamma_1+gamma_2) * h/2 * g(U); % momentum update
        
        omega = get_omega(local_param,gamma_2*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);
        
        iP = iP + (gamma_2 + gamma_3)/2 * h * g(U);
        
        omega = get_omega(local_param,gamma_3*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);
        
        iP = iP + (gamma_3 + gamma_4)/2 * h * g(U);

        omega = get_omega(local_param,gamma_4*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        iP = iP + (gamma_3 + gamma_4)/2 * h * g(U);

        omega = get_omega(local_param,gamma_3*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        iP = iP + (gamma_2 + gamma_3)/2 * h * g(U);

        omega = get_omega(local_param,gamma_2*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        iP = iP + (gamma_1 + gamma_2)/2 * h * g(U);

        omega = get_omega(local_param,gamma_1*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        if it < nit 
            % one can merge the last momentum update of step it and the first
            % momentum update of step it+1 
            iP = iP + gamma_1 * h * g(U);
        else
            % in the last step of the trajectory, there is no step it+1 
            % anymore, i.e., we use step size gamma_1/2 * h
            iP = iP + gamma_1/2 * h * g(U);
        end
    end
end