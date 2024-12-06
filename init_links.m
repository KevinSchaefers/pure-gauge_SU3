function U = init_links(starthotcold,nlinks)
% INIT_LINKS initializes the link field.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: U = INIT_LINKS(starthotcold,nlinks)
%--------------------------------------------------------------------------
% input:    starthotcold - 0 = cold start, otherwise hot start
%           nlinks - number of links that have to be initialized        
% output :  initialized link field U 
%--------------------------------------------------------------------------
    
    if starthotcold == 0
        % cold start -> all links are initialized as the identitiy
        % matrix
        U = zeros(nlinks,9);
        U(:,1) = 1;
        U(:,5) = 1;
        U(:,9) = 1;
    else
        c = 1/sqrt(2);
        tmp = c*randn(nlinks,8);
        U = exponential_map(tmp);
    end
end