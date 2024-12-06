function invA = inv_3x3_field(A)
% INV_3X3_FIELD computes the inverse of all 3x3-matrices that are stored
% row-wise in the matrix A of size N x 9.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: invA = INV_3X3_FIELD(A)
%--------------------------------------------------------------------------
% input:    A - matrix of size nlinks x 9 where each row represents a 3x3
%           matrix.  
% output :  invA - matrix of size nlinks x 9 where each row invA(i,:) 
%           represents the inverse of A(i,:)
%--------------------------------------------------------------------------
    
    % formula for inverse of 3x3 matrix based on det and adj
    invA = [A(:,5).*A(:,9) - A(:,6).*A(:,8), ...
            A(:,3).*A(:,8) - A(:,2).*A(:,9), ...
            A(:,2).*A(:,6) - A(:,3).*A(:,5), ...
            A(:,6).*A(:,7) - A(:,4).*A(:,9), ...
            A(:,1).*A(:,9) - A(:,3).*A(:,7), ...
            A(:,3).*A(:,4) - A(:,1).*A(:,6), ...
            A(:,4).*A(:,8) - A(:,5).*A(:,7), ...
            A(:,2).*A(:,7) - A(:,1).*A(:,8), ...
            A(:,1).*A(:,5) - A(:,2).*A(:,4)];

    invA = invA./(A(:,1).*invA(:,1) + A(:,2) .*invA(:,4) + A(:,3).*invA(:,7));
    % invA has size nlinks x 9
end