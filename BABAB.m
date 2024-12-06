function [U,iP] = BABAB(U, iP, local_param, h, traj)
% BABAB computes one trajectory of the molecular dynamics (MD) step using 
% the integrator BABAB, the second-order min.-norm scheme proposed by 
% Omelyan et al. which is a five-stage velocity-version of convergence 
% order two. The integrator reads
%   Phi_h = exp(lambda*h*B)exp(h/2 * A)exp((1-2*lambda)*h*B)exp(h/2 * A)exp(lambda*h*B)
% where A denotes a link update and B denotes a momentum update. 
% The used integrator coefficient is lambda = 0.19318333. For more
% details on the integrator, see <a href="matlab:web('https://doi.org/10.1016/S0010-4655(02)00754-3')">[Omelyan et al. 2003]</a>.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: [U,iP] = BABAB(U, iP, local_param, h, traj)
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
    
    lambda = 0.19318333; % integrator coefficient 

    % first step (momentum update)
    iP = iP + lambda*h*g(U);   % momentum update

    nit = traj/abs(h);  % nit time steps of the integrator 
    for it=1:nit
        % link update
        omega = get_omega(local_param,0.5*h*iP); % omega_1 = h/2 * dPsi_{omega}^{-1}(iP)
        psi_iP = feval(local_param, omega); % applying the local parameterization 
            % to the Lie algebra element omega, resulting in Lie group element psi_iP
        U = matMultField(psi_iP,U); % multiplication of Lie group elements

        iP = iP + (1-2*lambda)*h*g(U);

        omega = get_omega(local_param,0.5*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        if it < nit 
            % one can merge the last momentum update of step it and the first
            % link update of step it+1 
            iP = iP + 2*lambda*h*g(U);
        else 
            % in the last step of the trajectory, there is no step it+1 
            % anymore, i.e., we use step size lambda*h
            iP = iP + lambda*h*g(U);
        end
    end
end