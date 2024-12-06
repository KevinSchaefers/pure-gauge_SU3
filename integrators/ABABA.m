function [U,iP] = ABABA(U, iP, local_param, h, traj)
% ABABA computes one trajectory of the molecular dynamics (MD) step using 
% the integrator ABABA, the second-order min.-norm scheme proposed by 
% Omelyan et al. which is a five-stage position-version of convergence 
% order two. The integrator reads
%   Phi_h = exp(lambda*h*A)exp(h/2 * B)exp((1-2*lambda)*h*A)exp(h/2 * B)exp(lambda*h*A)
% where A denotes a link update and B denotes a momentum update. 
% The used integrator coefficient is lambda = 0.19318333. For more
% details on the integrator, see <a href="matlab:web('https://doi.org/10.1016/S0010-4655(02)00754-3')">[Omelyan et al. 2003]</a>.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: [U,iP] = ABABA(U, iP, local_param, h, traj)
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
    
    % integrator coefficient
    lambda = 0.19318333;

    % first half-step (link update)
    omega = get_omega(local_param,lambda*h*iP); 
    psi_iP = feval(local_param, omega); % applying the local parameterization 
        % to the Lie algebra element omega, resulting in Lie group element psi_iP
    U = matMultField(psi_iP,U); % multiplication of Lie group elements

    % iteration
    nit = traj/abs(h);  % nit time steps of the integrator 
    for it=1:nit
        iP = iP + 0.5*h*g(U);   % momentum update

        omega = get_omega(local_param,(1-2*lambda)*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        iP = iP + 0.5*h*g(U);

        if it < nit 
            % one can merge the last link update of step it and the first
            % link update of step it+1 
            omega = get_omega(local_param,2*lambda*h*iP);
            psi_iP = feval(local_param, omega);
            U = matMultField(psi_iP,U);
        else 
            % in the last step of the trajectory, there is no step it+1 
            % anymore, i.e., we use step size lambda*h
            omega = get_omega(local_param,lambda*h*iP);
            psi_iP = feval(local_param, omega);
            U = matMultField(psi_iP,U);
        end
    end
end