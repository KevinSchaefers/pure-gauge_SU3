function E_kin = kinetic_energy(iP)
% KINETIC_ENERGY computes the kinetic energy 
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: E_kin = KINETIC_ENERGY(iP)
%--------------------------------------------------------------------------
% input:    iP              : conjugate momenta
% output:   E_kin           : kinetic energy
%--------------------------------------------------------------------------

    E_kin = sum(sum(iP.^2));
end