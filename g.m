function y = g(U)
% G computes the gauge force, i.e., the sum of two stables belonging to one 
% link (all) in the Gell-Mann basis. 
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: y = G(U)
%--------------------------------------------------------------------------
% input: link field U
% output : gauge force y
%--------------------------------------------------------------------------
    global beta;
    % beta  : coupling constant 

    staple = sum_of_staples(U); % computes the sum of staples
    w_TA = TA(U,staple);    % applies the traceless and anti-Hermitian operator
    
    y = -beta/3 * w_TA;
end

function v = sum_of_staples(U)
% SUM_OF_STAPLES computes the sum of the two staples belonging to one link
% (all).
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: v = sum_of_staples(U)
%--------------------------------------------------------------------------
% input: link field U
% output : v - sum of staples
%--------------------------------------------------------------------------
    
  global hop; %field of indices of adjacent link matrices

  nlinks = size(U,1);
  n=nlinks/2;

  v = zeros(size(U));

  % sum of the 2 staples belonging to one horizontal link
  v1 = matMultField( U(hop(1:n,2),:) , conjugateTransposeField(U(hop(1:n,3),:)));
  v1 = matMultField(v1,conjugateTransposeField(U(hop(1:n,4),:)));
  v2 = matMultField(conjugateTransposeField(U(hop(1:n,5),:)),...
      conjugateTransposeField(U(hop(1:n,6),:)));
  v2 = matMultField(v2,U(hop(1:n,7),:));
  v(1:n,:)=v1+v2;
  
  % sum of the 2 staples belonging to one vertical link
  v1 = matMultField(conjugateTransposeField(U(hop(n+1:nlinks,2),:)),...
      conjugateTransposeField(U(hop(n+1:nlinks,3),:)));
  v1 = matMultField(v1,U(hop(n+1:nlinks,4),:));
  v2 = matMultField(U(hop(n+1:nlinks,5),:),...
      conjugateTransposeField(U(hop(n+1:nlinks,6),:)));
  v2 = matMultField(v2,conjugateTransposeField(U(hop(n+1:nlinks,7),:)));
  v(n+1:nlinks,:)=v1+v2;
end

function w_algebra = TA(U, staple)
% TA applies the traceless and anti-Hermitian operator { }_{TA}.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: w_algebra = TA(U,staple)
%--------------------------------------------------------------------------
% input:    U - link field
%           staple - staples (result from sum_of_staples()) 
% output :  w_algebra - result after application of the traceless and
%           anti-Hermitian operator that yields elements in su(3)
%--------------------------------------------------------------------------
    N=3;
    nlinks = size(U,1);
    w_algebra = zeros(nlinks,8);

    M = matMultField(U,staple);
    M = M - conjugateTransposeField(M);
    traceField = M(:,1) + M(:,5) + M(:,9);
    tmp = 0.5 * M;
    tmp(:,1) = tmp(:,1) - 1/(2*N) * traceField;
    tmp(:,5) = tmp(:,5) - 1/(2*N) * traceField;
    tmp(:,9) = tmp(:,9) - 1/(2*N) * traceField;

    w_algebra(:,1) = imag(tmp(:,2));
    w_algebra(:,2) = real(tmp(:,2));
    w_algebra(:,3) = imag(tmp(:,1))+ 0.5*imag(tmp(:,9));
    w_algebra(:,4) = imag(tmp(:,3));
    w_algebra(:,5) = real(tmp(:,3));
    w_algebra(:,6) = imag(tmp(:,6));
    w_algebra(:,7) = real(tmp(:,6));
    w_algebra(:,8) = (-sqrt(3)/2)*imag(tmp(:,9));
end