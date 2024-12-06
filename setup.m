function setup(L,D)
% SETUP initializes the lattice and introduces the basis matrices of su(3) 
% (the Gell-Mann matrices multiplied by the imaginary unit i) as global 
% variables.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: SETUP(L,D)
%--------------------------------------------------------------------------
% input:    L - number of lattice points in one direction
%           D - dimension of the lattice (currently, only D=2 is supported)
% output :  none (lattice and Gell-Mann matrices are global)
%--------------------------------------------------------------------------

    lattice(L,D);               % computes hop including the adjacent links/momenta
                                % dimension is 2*L^D x 7
    
    setGellmann();              % sets the Gell-Mann-matrices as global variables                        
end

function lattice(L,D)
% LATTICE initializes the lattice and calculates the indices of adjacent
% links/momenta
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: LATTICE(L,D)
%--------------------------------------------------------------------------
% input:    L - number of lattice points in one direction
%           D - dimension of the lattice (currently, only D=2 is supported)
% output :  none - the indices of adjacent links/momenta are stored in the
%           global variable 'hop'
%--------------------------------------------------------------------------
    
    global lattice_size;        % size of the lattice
    nlinks = 2*lattice_size;    % number of links
    global hop;
    hop = zeros(nlinks,7);      % stores the adjacent links / momenta
    
    switch D
        case {1}
            error('dimension too small!')
        case {2}
            disp(strcat('Lattice in 2 dimensions : ', num2str(L), 'x', num2str(L)));
            
            % vertically adjacent plaquettes:
            % hop:       ___3___
            %           |       |
            %           4       2
            %           |___1___|
            %           |       |
            %           7       5   
            %           |___6___|
            for j=1:L
                for k=1:L
                    index = L*(j-1) + k;
                    
                    hop(index,1) = index;                                           % link itself 
                    hop(index,2) = L*(j-1) + mod(k,L) + lattice_size+1;             % right top
                    hop(index,3) = L*mod(j,L) + k;                                  % top
                    hop(index,4) = L*(j-1) + k + lattice_size;                      % left top
                    hop(index,5) = L*mod(j-2+L,L) + mod(k,L) + lattice_size+1;      % right bottom
                    hop(index,6) = L*mod(j-2+L,L) + k;                              % bottom
                    hop(index,7) = L*mod(j-2+L,L) + mod(k-1+L,L) + lattice_size+1;  % left bottom
                end
            end
            
            % horizontally adjacent plaquettes:
            % hop:      ___2_______5___
            %          |       |       |
            %          3       1       6
            %          |___4___|___7___|
            for j=1:L
                for k=1:L
                    index = L*(j-1) + k + lattice_size;
                    
                    hop(index,1) = index;                                           % link itself
                    hop(index,2) = L*mod(j,L) + mod(k-2+L,L) + 1;                   % left top
                    hop(index,3) = L*(j-1) + mod(k-2+L,L) + lattice_size + 1;       % left
                    hop(index,4) = L*(j-1) + mod(k-2+L,L) + 1;                      % left bottom
                    hop(index,5) = L*mod(j,L) + k;                                  % right top
                    hop(index,6) = L*(j-1) + mod(k,L) + lattice_size + 1;           % right
                    hop(index,7) = L*(j-1) + k;                                     % right bottom 
                end
            end
        case {3}
            error('3 dimensions - not yet implemented');
        case {4} 
            error('4 dimensions - not yet implemented');
        otherwise
            error('dimension is too large!');
    end

end

function setGellmann()
% SETGELLMANN sets the eight Gell-Mann matrices (multiplied by the 
% imaginary unit i) as the basis matrice sof the Lie algebra su(3).
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: SETGELLMANN()
%--------------------------------------------------------------------------
% input:    none
% output :  none (basis matrices B1,...,B8 are global)
%--------------------------------------------------------------------------
    
    global B1;
    global B2;
    global B3;
    global B4;
    global B5;
    global B6;
    global B7;
    global B8;
    
    i = sqrt(-1);
    
    B1 = i * [0,1,0;
          1,0,0;
          0,0,0];
    B2 = i * [0,-i,0;
          i, 0,0;
          0, 0,0];
    B3 = i * [1, 0,0;
          0,-1,0;
          0, 0,0];
    B4 = i * [0,0,1;
          0,0,0;
          1,0,0];
    B5 = i * [0,0,-i;
          0,0, 0;
          i,0, 0];
    B6 = i * [0,0,0;
          0,0,1;
          0,1,0];
    B7 = i * [0,0, 0;
          0,0,-i;
          0,i, 0];
    B8 = i * [1/sqrt(3),         0,          0;
                  0, 1/sqrt(3),          0;
                  0,         0, -2/sqrt(3)];
end