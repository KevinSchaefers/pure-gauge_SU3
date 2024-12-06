function [U,iP] = integration_step(U,iP,method_string,h,traj,local_param)
% INTEGRATION_STEP calls the preferred integrator using the preferred local
% parameterization to perform the molecular dynamics (MD) step.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: [U,iP] = INTEGRATION_STEP(U,iP,method_string,h,traj,local_param)
%--------------------------------------------------------------------------
% input:    U             :     links
%           iP            :     conjugate momenta
%           method_string :     string of the preferred integration
%                               method
%           h             :     step size
%           traj          :     trajectory length
%           local_param   :     string of the preferred local
%                               parameterization     
% output :  updated link field U and conjugate momenta iP
%--------------------------------------------------------------------------
    
    switch method_string
        case 'BAB'
            % leapfrog / Stoermer Verlet method (velocity version)
            [U,iP] = BAB(U,iP,local_param,h,traj);
        case 'BABAB'
            % second-order min.-norm scheme velocity version 
            % [Omelyan et al. 2003]
            [U,iP] = BABAB(U,iP,local_param,h,traj);
        case 'ABABA'
            % second-order min.-norm scheme position version 
            % [Omelyan et al. 2003]
            [U,iP] = ABABA(U,iP,local_param,h,traj);
        case 'OMF4'
            % fourth-order min.-norm scheme (11-stage velocity version)
            % [Omelyan et al. 2003] <- ID: BABABABABAB
            [U,iP] = OMF4(U,iP,local_param,h,traj);
        % Hessian-free force-gradient integrators
        case 'BADAB'
            % five-stage Hessian-free force-gradient integrator of order
            % four [Schaefers et al. 2024]
            [U,iP] = BADAB(U,iP,local_param,h,traj);
        case 'ABADABA'
            % seven-stage Hessian-free force-gradient integrator of order
            % four [Schaefers et al. 2024]
            [U,iP] = ABADABA(U,iP,local_param,h,traj);
        % composition techniques based on BAB
        case 'Yoshida'
            % Yoshida's triple-jump applied to BAB, resulting in a
            % seven-stage integrator of order four
            [U,iP] = Yoshida(U,iP,local_param,h,traj);
        case 'Suzuki'
            % Suzuki's fractals applied to BAB, resulting in a eleven-stage
            % integrator of order four
            [U,iP] = Suzuki(U,iP,local_param,h,traj);
        case 'advancedComposition'
            % advanced comoposition applied to BAB, resulting in a
            % 15-stage sixth-order integrator
            [U,iP] = advancedComposition(U,iP,local_param,h,traj);
        otherwise 
            error('Preferred method unknown!');
    end
end