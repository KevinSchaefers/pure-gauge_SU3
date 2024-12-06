function omega = get_omega(local_param,A)
% GET_OMEGA computes d_\Psi_{\Omega=0}^{-1}(A) where \Psi denotes the local
% parameterization used to map from su(3) to SU(3).
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: omega = GET_OMEGA(local_param,A)
%--------------------------------------------------------------------------
% input:    local_param - the local parameterization used to map from su(3)
%               to SU(3)
%           A - Lie algebra elements in su(3) (usually A = iP in this code)
% output :  omega - equals d_\Psi_{\Omega=0}^{-1}(A) where \Psi is the 
%               local parameterization     
%--------------------------------------------------------------------------
    switch local_param
        case {'exponential_map'} % matrix exponential 
            omega = A;
        case {'caymod'} % modified Cayley transform
            omega = 0.5*A;
        otherwise
            error('no valid local parameterization!');
    end
end