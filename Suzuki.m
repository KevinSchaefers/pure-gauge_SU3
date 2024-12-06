function [U,iP] = Suzuki(U, iP, local_param, h, traj)
% Suzuki computes one trajectory of the molecular dynamics (MD) step using 
% Suzuki's fractals and BAB as the base method, resulting in a eleven-stage 
% velocity-version of convergence order four. The composition scheme has been
% introduced in <a
% href="matlab:web('https://doi.org/10.1016/0375-9601(90)90962-N')">[Suzuki 1990]</a>.
% and uses the composition weights 
% gamma_1 = gamma_2 = gamma_4 = gamma_5 = 1/(4-nthroot(4,3)),
% gamma_3 = 1 - 4*gamma_1.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: [U,iP] = Suzuki(U, iP, local_param, h, traj)
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
    gamma_1 = 1/(4-nthroot(4,3)); % = gamma_2 = gamma_4 = gamma_5
    gamma_3 = 1-4*gamma_1;
    
    % first step (momentum update)
    iP = iP + gamma_1 * 0.5 * h * g(U); % momentum update

    nit = traj/abs(h); % nit time steps of the integrator 
    for it=1:nit
        % link update
        omega = get_omega(local_param,gamma_1*h*iP); % omega_1 = gamma_1 h dPsi_{omega}^{-1}(iP)
        psi_iP = feval(local_param, omega); % applying the local parameterization 
        % to the Lie algebra element omega, resulting in Lie group element psi_iP
        U = matMultField(psi_iP,U); % multiplication of Lie group elements
        
        iP = iP + gamma_1 * h * g(U);
        
        omega = get_omega(local_param,gamma_1*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);
        
        iP = iP + (gamma_1 + gamma_3)/2 * h * g(U);
        
        omega = get_omega(local_param,gamma_3*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);
        
        iP = iP + (gamma_1 + gamma_3)/2 * h * g(U);

        omega = get_omega(local_param,gamma_1*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        iP = iP + gamma_1 * h * g(U);

        omega = get_omega(local_param,gamma_1*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        if it < nit
            % one can merge the last momentum update of step it and the first
            % link update of step it+1 
            iP = iP + gamma_1 * h * g(U);
        else
            % in the last step of the trajectory, there is no step it+1 
            % anymore, i.e., we use step size 0.5*gamma_1*h
            iP = iP + gamma_1/2 * h * g(U);
        end
    end
end