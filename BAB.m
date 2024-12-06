function [U,iP] = BAB(U, iP, local_param, h, traj)
% BAB computes one trajectory of the molecular dynamics (MD) step using 
% the integrator BAB, the velocity version of the Stoermer-Verlet/leapfrog
% method of convergence order two.
% The integrator reads
%   Phi_h = exp(h/2 * B)exp(h*A)exp(h/2 * B)
% where A denotes a link update and B denotes a momentum update. 
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: [U,iP] = BAB(U, iP, local_param, h, traj)
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

    % first step (momentum update)
    iP = iP + 0.5*h*g(U);   % momentum update

    nit = traj/abs(h);  % nit time steps of the integrator 
    for it=1:nit
        % link update
        omega = get_omega(local_param,h*iP); % omega_1 = lambda h dPsi_{omega}^{-1}(iP)
        psi_iP = feval(local_param, omega); % applying the local parameterization 
            % to the Lie algebra element omega, resulting in Lie group element psi_iP
        U = matMultField(psi_iP,U); % multiplication of Lie group elements

        if it < nit 
            % one can merge the last momentum update of step it and the first
            % link update of step it+1 
            iP = iP + h*g(U);
        else 
            % in the last step of the trajectory, there is no step it+1 
            % anymore, i.e., we use step size lambda*h
            iP = iP + 0.5*h*g(U);
        end
    end
end