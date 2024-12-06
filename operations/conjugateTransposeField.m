function res = conjugateTransposeField(A)
% CONJUGATETRANSPOSEFIELD applies the dagger operator (complex conjugation 
% + transpose) to all matrices stored in the matrix A 
% (matrix of size nlinks x 9).
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: res = CONJUGATETRANSPOSEFIELD(A)
%--------------------------------------------------------------------------
% input:    A - matrix of size N x 9 where each row represents a 3x3
%           matrix.  
% output :  res - matrix of size N x 9 where res(i,:) represents 
%           conj(transpose(A(i,:))).
%--------------------------------------------------------------------------
    tmp = conj(A);
    res = zeros(size(A));
    res(:,1) = tmp(:,1);
    res(:,2) = tmp(:,4);
    res(:,3) = tmp(:,7);
    res(:,4) = tmp(:,2);
    res(:,5) = tmp(:,5);
    res(:,6) = tmp(:,8);
    res(:,7) = tmp(:,3);
    res(:,8) = tmp(:,6);
    res(:,9) = tmp(:,9);
end