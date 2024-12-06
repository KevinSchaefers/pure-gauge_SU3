function [U,iP] = OMF4(U, iP, local_param, h, traj)
% OMF4 computes one trajectory of the molecular dynamics (MD) step using 
% the integrator OMF4, the fourth-order min.-norm scheme proposed by Omelyan et 
% al. which is a eleven-stage velocity-version. 
% The integrator reads
%   Phi_h = exp(b1*h*B)exp(a2*h*A)exp(b2*h*B)exp(a3*h*A)exp((0.5-b1-b2)*h*B)
%               exp((1-2(a2+a3))*h*A)exp((0.5-b1-b2)*h*B)exp(a3*h*A)
%                   exp(b2*h*B)exp(a2*h*A)exp(b1*h*B)
% where A denotes a link update and B denotes a momentum update. 
% The used integrator coefficients are 
% a2 = 0.253978510841060,
% a3 = -0.032302867652700,
% b1 = 0.083983152628767,
% b2 = 0.682236533571909.
% For more details on the integrator, see <a href="matlab:web('https://doi.org/10.1016/S0010-4655(02)00754-3')">[Omelyan et al. 2003]</a>.
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

    % integrator coefficients
    a2 = 0.253978510841060;
    a3 = -0.032302867652700;
    b1 = 0.083983152628767;
    b2 = 0.682236533571909;

    % first step (momentum update)
    iP = iP + b1 * h * g(U); % momentum update

    nit = traj/abs(h); % nit time steps of the integrator 
    for it=1:nit
        % link update 
        omega = get_omega(local_param,a2*h*iP); % omega_1 = a2 * h dPsi_{omega}^{-1}(iP)
        psi_iP = feval(local_param, omega); % applying the local parameterization 
        % to the Lie algebra element omega, resulting in Lie group element psi_iP
        U = matMultField(psi_iP,U); % multiplication of Lie group elements

        iP = iP + b2*h * g(U);

        omega = get_omega(local_param,a3*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        iP = iP + (0.5 - (b1+b2))*h * g(U);

        omega = get_omega(local_param,(1-2*(a2+a3))*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        iP = iP + (0.5 - (b1+b2))*h * g(U);

        omega = get_omega(local_param,a3*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        iP = iP + b2*h * g(U);

        omega = get_omega(local_param,a2*h*iP);
        psi_iP = feval(local_param, omega);
        U = matMultField(psi_iP,U);

        if it < nit
            % one can merge the last link update of step it and the first
            % link update of step it+1 
            iP = iP + 2*b1 * h * g(U);
        else 
            % in the last step of the trajectory, there is no step it+1 
            % anymore, i.e., we use step size b1*h
            iP = iP + b1 * h * g(U);
        end
    end
end