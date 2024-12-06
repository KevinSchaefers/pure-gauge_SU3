function detA = detField(A) 
% DETFIELD computes the determinant of all matrices stored in the matrix A 
% (matrix of size nlinks x 9).
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: detA = DETFIELD(A)
%--------------------------------------------------------------------------
% input:    A - matrix of size nlinks x 9 where each row represents a 3x3
%           matrix.  
% output :  detA - vector of size nlinks x 1 where detA(i) is the 
%           determinant of the matrix that is stored in A(i,:).
%--------------------------------------------------------------------------
    detA = A(:,1) .* (A(:,5) .* A(:,9) - A(:,6) .* A(:,8)) ... 
        + A(:,2) .* (A(:,6) .* A(:,7) - A(:,4) .* A(:,9)) ...
        + A(:,3) .* (A(:,4) .* A(:,8) - A(:,5) .* A(:,7));
end