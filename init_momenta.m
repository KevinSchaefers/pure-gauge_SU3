function iP = init_momenta(nlinks)
% INIT_MOMENTA initializes the conjugated momenta. P is traceless and
% Hermitian, iP is traceless and anti-Hermitian and thus in the Lie algebra
% su(3).
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: iP = INIT_MOMENTA(nlinks)
%--------------------------------------------------------------------------
% input:    nlinks - number of momenta that have to be initialized       
% output :  iP - stores the initialized momenta iP as linear combination of
% the basis matrices, i.e., it is a 2D array of size nlinks x 8 
%--------------------------------------------------------------------------
                                
    c = 1/sqrt(2);
    iP = c*randn(nlinks,8);
end