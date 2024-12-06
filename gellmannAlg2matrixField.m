function B = gellmannAlg2matrixField(A)
% GELLMANNALG2MATRIXFIELD transforms any representation of a Lie algebra
% element in su(3) as element in R^8 into its matrix representation, stored
% row-wise as a vector in R^9.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: B = GELLMANNALG2MATRIXFIELD(A)
%--------------------------------------------------------------------------
% input:    A - matrix of size nlinks x 8 where each row is a representation of
%           a Lie algebra element in su(3) as a vector in R^8
% output :  B - matrix of size nlinks x 9 where each row B(i,:) is the matrix
%           representation of the Lie algebra element A(i,:), stored
%           row-wise as a vector in R^9.
%--------------------------------------------------------------------------
    
    nlinks = size(A, 1);
    B = zeros(nlinks, 9);
    
    i = sqrt(-1);

    % Define intermediate variables for better readability
    term1 = i * A(:,1);
    term3 = i * A(:,3);
    term4 = i * A(:,4);
    term6 = i * A(:,6);
    term8 = (i/sqrt(3)) * A(:,8);
    
    % Compute matrix elements
    B(:,1) = term3 + term8;
    B(:,2) = A(:,2) + term1;
    B(:,3) = A(:,5) + term4;
    B(:,4) = -conj(B(:,2)); %B(:,4) = -A(:,2) + i*A(:,1);
    B(:,5) = -term3 + term8;
    B(:,6) = A(:,7) + term6;
    B(:,7) = -conj(B(:,3)); %B(:,7) = -A(:,5) + i*A(:,4);
    B(:,8) = -conj(B(:,6)); %B(:,8) = -A(:,7) + i*A(:,6);
    B(:,9) = -2*term8;
end