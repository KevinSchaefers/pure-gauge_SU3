function expA = exponential_map(A)
% EXPONENTIAL_MAP applies the matrix exponential to map elements from su(3)
% into SU(3). For details on the computation of the matrix exponential for
% SU(3), see <a
% href="matlab:web('https://link.springer.com/article/10.1140/epja/s10050-022-00816-5')">[Kaiser 2022]</a>.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: expA = EXPONENTIAL_MAP(A)
%--------------------------------------------------------------------------
% input:    A               : array of size nlinks x 8 storing where each
%   row corresponds to a single element in su(3) that is uniquely defined
%   as a linear combination of the basis matrices. The coefficients of this
%   linear combination are stored in A.
% output:   expA            : Lie group elements in SU(3), stored as a 2D
%   array of size nlinks x 9.
%--------------------------------------------------------------------------

    nlinks = size(A,1);
    norm_val = sqrt(sum(A.^2,2));
    Sigma = A./(norm_val);
    tmp = Sigma.^2;

    det_Sigma = -(Sigma(:,3).* (-tmp(:,4) - tmp(:,5) + tmp(:,6) + tmp(:,7)) ...
        - 2*(Sigma(:,1).*Sigma(:,4) + Sigma(:,2).*Sigma(:,5)).*Sigma(:,6) ...
        + 2*(Sigma(:,2).*Sigma(:,4) - Sigma(:,1).*Sigma(:,5)).*Sigma(:,7) ...
        + sqrt(3)/3 * (-2*tmp(:,1)-2*tmp(:,2)-2*tmp(:,3)+tmp(:,4)...
        +tmp(:,5)+tmp(:,6)+tmp(:,7)+ (2/3).*tmp(:,8)).*Sigma(:,8));

    Psi = 1/3 * acos((3*sqrt(3)/2) .* det_Sigma);
    cos_Psi = cos(Psi);
    sin_Psi = sin(Psi);
    z1 = 2/sqrt(3) * cos_Psi; z1_squared = z1.^2;
    z2 = -sin_Psi - cos_Psi/sqrt(3); z2_squared = z2.^2;
    z3 = -(z1+z2); z3_squared = z1_squared + z2_squared + 2*z1.*z2;

    B = gellmannAlg2matrixField(Sigma)./1i;
    B_squared = matMultField(B,B);
    Id = repelem(reshape(eye(3),[1,9]),nlinks,1);

    expA = ( exp(1i*z1 .* norm_val)./(3*z1_squared - 1) ) .* ( (z1_squared - 1) .*  Id + z1.*B + B_squared) ...
        + ( exp(1i*z2 .* norm_val)./(3*z2_squared - 1) ) .* ( (z2_squared - 1) .*  Id + z2.*B + B_squared) ...
        + ( exp(1i*z3 .* norm_val)./(3*z3_squared - 1) ) .* ( (z3_squared - 1) .*  Id + z3.*B + B_squared);
end