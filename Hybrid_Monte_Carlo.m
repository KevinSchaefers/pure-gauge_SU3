function U = Hybrid_Monte_Carlo(L, method_string, param_beta, h, local_param, runname, ntraj, ncheckpoints, tau, U_start)
% HYBRID_MONTE_CARLO applies the HMC algorithm for pure gauge field
% simulations on a L x L lattice at specified inverse temperature param_beta.
% In the molecular dynamics (MD) step, it applies a volume-preserving and
% time-reversible integrator using a constant step size h that is defined 
% via method_string in order to satisfy the detailed balance condition. 
% For the integrators, different local parameterizations local_param can be
% used to obtain links in SU(3). All in all, the function computes ntraj
% trajectories of length tau.
%--------------------------------------------------------------------------
% Kevin Schaefers (v1, 2024)
%--------------------------------------------------------------------------
% call: U = HYBRID_MONTE_CARLO(L,method_string,param_beta,h,local_param,...
%               ostring,ntraj,ncheckpoints,tau,U_start)
%--------------------------------------------------------------------------
% input:    L - lattice size (2D lattice of size L x L)
%           method_string - string describing the method used for the MD
%           step. The following methods are implemented: ABABA, ABADABA,
%           advancedComposition, BAB, BABAB, BADAB, OMF4, Suzuki, Yoshida.
%           param_beta - inverse temperatur beta
%           h - step size used in the MD step
%           local_param - local parameterization used to compute links in
%           SU(3). The following local parameterizations are implemented:
%           exponential_map (matrix exponential for SU(3)), caymod 
%           (modified Cayley transform). For details on the computation of
%           the matrix exponential, see <a
%           href="matlab:web('https://link.springer.com/article/10.1140/epja/s10050-022-00816-5')">[Kaiser 2022]</a>.
%           For details on the modified Cayley transform for SU(3), see 
%           <a href="matlab:web('https://arxiv.org/abs/2406.11337')">[Schaefers Peardon Guenther 2024]</a>.
%           runname - name of the run to avoid overwriting files. The
%           results are then stored in results/runname/
%           ncheckpoints - after ncheckpoints trajectories, an intermediate
%           result will be stored.
%           tau - trajectory length. 
%           U_start - optional argument to start with a link field from a
%           previous run. Otherwise, the simulation starts with a random
%           link field.
% output :  U - link field
%--------------------------------------------------------------------------      

    % Initialization
    global nlinks;
    d = 2;                      % dimension of the lattice (currently, only 2D lattices are supported)
    global beta;
    beta = param_beta;
    global lattice_size;
    lattice_size = L^d;         % sizeof the lattice
    nlinks = 2*lattice_size;    % number of links
    setup(L,d);                 % general setup
    
    seed = 12345;               % seed for the random number generators
    rand('seed', seed);
    randn('seed',seed);
    
    if nargin < 10
        starthotcold = 1;           % hot start = 1, cold start = 0
        U = init_links(starthotcold,nlinks);   % initializes the links U
    else 
        % start with the link field from a previous run
        if size(U_start,1) == nlinks 
            U = U_start;
        else 
            error('The passed link field has another size than L x L!');
        end
    end
    

    if ~exist("results/"+runname, 'dir')
       mkdir("results/"+runname);
    end
    
    % declaring names for the output/checkpoint files
    out_file = strcat('results/',runname,'/SU(3)_',method_string, '_', num2str(L,'%1d'), 'x', ...
           num2str(L,'%1d'), '_beta',num2str(param_beta,'%1.1f'),'_h', ...
           num2str(h,'%1.4f'), '_', local_param, ...
           '.txt');

    checkpoint_file = strcat('results/',runname,'/SU(3)_',method_string, '_', num2str(L,'%1d'), 'x', ...
           num2str(L,'%1d'), '_beta',num2str(param_beta,'%1.1f'),'_h', ...
           num2str(h,'%1.4f'), '_', local_param, ...
           '.checkpoint.mat');
       
    disp(out_file);
    result = fopen(out_file, 'w+');
    
    % set variables
    acc_traj = 0;       % computed trajectories (accepted)      
    nloops = 0;         % number of computed trajectories (all)
    
    writeHeader(result,method_string,seed,h,tau,local_param);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [S, meanPlaquette] = action(U);     % initial action
    
    tic;    % measure the time
    startctime = cputime;
    while nloops < ntraj     % loop over the trajectories
        % compute new momenta iP
        iP = init_momenta(nlinks);
        E_kin = kinetic_energy(iP); % kinetic energy
        
        H = E_kin + S;              % Hamiltonian H = total energy
        
        % MOLECULAR DYNAMICS (MD) STEP
        U_old = U;                  % save U in case of rejected step
        [U,iP] = integration_step(U,iP,method_string, h, tau, local_param); 
        % integration step with the chosen method method_string and 
        % local parameterization local_param
        
        [S_new, meanPlaquette_new] = action(U); % action after the MD step
        E_kin_new = kinetic_energy(iP); % energy after the MD step
        H_new = E_kin_new + S_new;      % total energy after the MD step
        
        % HYBRID MONTE CARLO (HMC) ACCEPTANCE STEP
        deltaH = H - H_new; % difference in the Hamiltonian (H_old - H_new)
        [accepted, exp_deltaH, r] = Metropolis_step(H, H_new);
        
        if accepted == 1    % step accepted
            acc_traj = acc_traj + 1;
            
            % overwrite the values by the new ones
            S = S_new;
            meanPlaquette = meanPlaquette_new; 
            E_kin = E_kin_new;
            H = H_new;
        else % step rejected, we keep the old values S and U_old
            U = U_old;
        end
        
        nloops = nloops + 1;
        acceptance_rate = 100*acc_traj/nloops;
        
        % output:
        writeAcceptedStep(result, nloops, acc_traj, H, S, E_kin, deltaH, ...
            exp_deltaH, r, meanPlaquette, accepted, acceptance_rate);
        
        % write a checkpoint...
        if (mod(nloops,ncheckpoints)==0) % write every ncheckpoints steps a checkpoint
            wtime = toc;
            ctime = cputime - startctime;

            disp('loops:');
            disp(nloops);

            fprintf(result, '%s %1.8e s\n','% wallclocktime = ',wtime);
            fprintf(result, '%s %1.8e s\n','% cputime = ',ctime);

            seed_rand = rand('seed');
            seed_randn = randn('seed');
            save(checkpoint_file, 'U', 'U_old', 'meanPlaquette', 'E_kin', 'S', 'H', ...
                'nloops', 'seed_rand', 'seed_randn', 'accepted', 'acc_traj', ...
                'tau', 'wtime', 'ctime');
        end
    end

    wtime = toc;
    ctime = cputime - startctime;
    fprintf(result, '%s %1.8e s\n','% wallclocktime = ',wtime);
    fprintf(result, '%s %1.8e s\n','% cputime = ',ctime);
end

