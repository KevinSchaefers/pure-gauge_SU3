function cayA = caymod(A)
% CAYMOD applies the modified Cayley transform to map elements from su(3)
% into SU(3). For details on the modified Cayley transform, see <a
% href="matlab:web('https://arxiv.org/abs/2406.11337')">[Schaefers Peardon Guenther 2024]</a>.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: cayA = CAYMOD(A)
%--------------------------------------------------------------------------
% input:    A               : array of size nlinks x 8 storing where each
%   row corresponds to a single element in su(3) that is uniquely defined
%   as a linear combination of the basis matrices. The coefficients of this
%   linear combination are stored in A.
% output:   cayA            : Lie group elements in SU(3), stored as a 2D
%   array of size nlinks x 9.
%--------------------------------------------------------------------------
    
    B = gellmannAlg2matrixField(A);    % transform coefficients A into the corresponding matrix in the Lie algebra

    gamma_inv = -0.5*sum(A.^2,2)./imag_det_su3(A);
    gamma_inv(isinf(gamma_inv)) = 0;

    sin_theta = -0.5 * (gamma_inv -sign(gamma_inv).*sqrt(gamma_inv.^2 + 1));
    cos_theta = sqrt(1 - sin_theta.^2);

    inv_term = (1i.*sin_theta - cos_theta).*B;
    inv_term(:,1:4:9) = inv_term(:,1:4:9) + 1;
    inv_term = inv_3x3_field(inv_term);

    B = (cos_theta + 1i.*sin_theta).*B;
    B(:,1:4:9) = B(:,1:4:9) + 1;

    cayA = matMultField(inv_term,B);
end

function res = imag_det_su3(A)
    tmp = A.^2;
    res = A(:,3).* (-tmp(:,4) - tmp(:,5) + tmp(:,6) + tmp(:,7)) ...
        - 2*(A(:,1).*A(:,4) + A(:,2).*A(:,5)).*A(:,6) ...
        + 2*(A(:,2).*A(:,4) - A(:,1).*A(:,5)).*A(:,7) ...
        + sqrt(3)/3 * (-2*tmp(:,1)-2*tmp(:,2)-2*tmp(:,3)+tmp(:,4)...
        +tmp(:,5)+tmp(:,6)+tmp(:,7)+(2/3).*tmp(:,8)).*A(:,8);
end