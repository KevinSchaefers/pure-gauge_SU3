function [S, meanPlaquette] = action(U)
% ACTION computes the gauge action of the link field U
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: [S,meanPlaquette] = ACTION(U)
%--------------------------------------------------------------------------
% input: link field U
% output : gauge action S
%          meanPlaquette
%--------------------------------------------------------------------------
    
    global lattice_size beta hop;
    % lattice_size; % number of the plaquettes
    % beta;         % factor
    % hop;          % field with indices of adjacent link matrices
    % U;
    
    N = 3; % dimension of the link matrices: N x N
    
    % plaquette
    w = matMultField(U(hop(1:lattice_size,1),:),U(hop(1:lattice_size,2),:));
    w = matMultField(w,conjugateTransposeField(U(hop(1:lattice_size,3),:)));
    w = matMultField(w,conjugateTransposeField(U(hop(1:lattice_size,4),:)));
    
    plaq = sum(real(w(:,1))+real(w(:,5))+real(w(:,9)));
    
    c = -beta/N;
    S = c * plaq + beta * lattice_size; % action      
    meanPlaquette = plaq/lattice_size;
end