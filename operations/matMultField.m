function C = matMultField(A,B)
% MATMULTFIELD computes the matrix products C(i) = A(i)*B(i) of 3x3 matrices
% where the matrices are stored row-wise in A, B and C, respectively.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: C = MATMULTFIELD(A,B)
%--------------------------------------------------------------------------
% input:    A - matrix of size nlinks x 9 where each row represents a 3x3
%           matrix. 
%           B - matrix of size nlinks x 9 where each row represents a 3x3
%           matrix. 
% output :  C - matrix of size nlinks x 9 where each row represents a 3x3
%           matrix. C(i,:) is the result of the matrix-matrix
%           multiplication of A(i,:) and B(i,:).
%--------------------------------------------------------------------------
    
    C = zeros(size(A)); % C = zeros(nlinks,9)
    C(:,1) = A(:,1) .* B(:,1) + A(:,2) .* B(:,4) + A(:,3) .* B(:,7);   % first
    C(:,2) = A(:,1) .* B(:,2) + A(:,2) .* B(:,5) + A(:,3) .* B(:,8);   % ...
    C(:,3) = A(:,1) .* B(:,3) + A(:,2) .* B(:,6) + A(:,3) .* B(:,9);   % row
    C(:,4) = A(:,4) .* B(:,1) + A(:,5) .* B(:,4) + A(:,6) .* B(:,7);   % second
    C(:,5) = A(:,4) .* B(:,2) + A(:,5) .* B(:,5) + A(:,6) .* B(:,8);   % ...
    C(:,6) = A(:,4) .* B(:,3) + A(:,5) .* B(:,6) + A(:,6) .* B(:,9);   % row
    C(:,7) = A(:,7) .* B(:,1) + A(:,8) .* B(:,4) + A(:,9) .* B(:,7);   % third
    C(:,8) = A(:,7) .* B(:,2) + A(:,8) .* B(:,5) + A(:,9) .* B(:,8);   % ...
    C(:,9) = A(:,7) .* B(:,3) + A(:,8) .* B(:,6) + A(:,9) .* B(:,9);   % row
end