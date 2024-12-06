function [U,iP] = ABADABA(U, iP, local_param, h, traj)
% ABADABA computes one trajectory of the molecular dynamics (MD) step using 
% the integrator ABADABA, a seven-stage Hessian-free force-gradient 
% integrator of convergence order four that has been introuduced in 
% <a href="matlab:web('https://doi.org/10.48550/arXiv.2403.10370')">[Schaefers et al. 2024]</a>.
% The integrator reads
%   Phi_h = exp(a1*h*A)exp(b1*h*B)exp((0.5-a1)*h*A)exp((1-2*b1)*h*D(h,1-2*b1,c2))exp((0.5-a1)*h*A)exp(b1*h*B)exp(a1*h*A)
% where A denotes a link update, B denotes a momentum update, and D denotes
% a Hessian-free force-gradient step. The used integrator coefficients
% are a1 = 0.089775972994422, b1 = 0.247597680043986, 
% c2 = 0.006911440413815.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: [U,iP] = ABADABA(U, iP, local_param, h, traj)
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

    % integrator coefficients
    a1 = 0.089775972994422;
    b1 = 0.247597680043986;
    c2 = 0.006911440413815;

    % first step (link update)
    omega = get_omega(local_param,a1*h*iP); % omega_1 = a1 * h dPsi_{omega}^{-1}(iP)
    psi_iP = feval(local_param, omega); % applying the local parameterization 
        % to the Lie algebra element omega, resulting in Lie group element psi_iP
    U = matMultField(psi_iP,U); % multiplication of Lie group elements

    nit = traj/abs(h); % nit time steps of the integrator 
    for it=1:nit
        iP = iP + b1 * h * g(U); % momentum update

        omega = get_omega(local_param,(0.5-a1)*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        U_tmp = matMultField(feval(local_param,get_omega(local_param,(2*c2)/(1-2*b1) * h^2 * g(U))),U); % temporary link-update
        iP = iP + (1-2*b1)*h * g(U_tmp); % Hessian-free force-gradient step using the temporary link field U_tmp

        omega = get_omega(local_param,(0.5-a1)*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        iP = iP + b1 * h * g(U);
        
        if it < nit
            % one can merge the last link update of step it and the first
            % link update of step it+1 
            omega = get_omega(local_param,2*a1*h*iP);
            psi_iP = feval(local_param, omega);
            U = matMultField(psi_iP,U);
        else
            % in the last step of the trajectory, there is no step it+1 
            % anymore, i.e., we use step size a1*h
            omega = get_omega(local_param,a1*h*iP);
            psi_iP = feval(local_param, omega);
            U = matMultField(psi_iP,U);
        end
    end
end