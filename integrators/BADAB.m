function [U,iP] = BADAB(U, iP, local_param, h, traj)
% BADAB computes one trajectory of the molecular dynamics (MD) step using 
% the integrator BADAB, a five-stage Hessian-free force-gradient 
% integrator of convergence order four that has been introuduced in 
% <a href="matlab:web('')">[Yin and Mawhinney 2011]</a>.
% The integrator reads
%   Phi_h = exp(h/6 * B)exp(h/2 * A)exp((2/3)*h*D(h,2/3,1/72))exp(h/2 * A)exp(h/6 * B)
% where A denotes a link update, B denotes a momentum update, and D denotes
% a Hessian-free force-gradient step. For more details on Hessian-free
% force-gradient integrators, see <a href="matlab:web('https://doi.org/10.48550/arXiv.2403.10370')">[Schaefers et al. 2024]</a>.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: [U,iP] = BADAB(U, iP, local_param, h, traj)
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
    iP = iP + h/6 * g(U); % momentum update

    % iteration
    nit = traj/abs(h); % nit time steps of the integrator
    for it=1:nit
        % link update
        omega = get_omega(local_param,0.5*h*iP); % omega_1 = 0.5 * h dPsi_{omega}^{-1}(iP)
        psi_iP = feval(local_param, omega); % applying the local parameterization 
            % to the Lie algebra element omega, resulting in Lie group element psi_iP
        U = matMultField(psi_iP,U); % multiplication of Lie group elements

        U_tmp = matMultField(feval(local_param,get_omega(local_param,h^2/24 * g(U))),U); % temporary link-update
        iP = iP + 2*h/3 * g(U_tmp); % Hessian-free force-gradient step using the temporary link field U_tmp

        omega = get_omega(local_param,0.5*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        if it < nit
            % one can merge the last momentum update of step it and the first
            % link update of step it+1 
            iP = iP + 1/3 * h * g(U);
        else 
            % in the last step of the trajectory, there is no step it+1 
            % anymore, i.e., we use step size h/6
            iP = iP + 1/6 * h * g(U);
        end
    end
end