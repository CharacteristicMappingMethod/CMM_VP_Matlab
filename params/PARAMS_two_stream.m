n=2^8; % 4097
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.Nmap           = n; % grid dimension of the submaps n x n
params.Nfine        =2^8; % grid for measuring and saves of f
params.Nsampling    =2^9; % upsampling grid of veloctiy
params.Nplotting    =1024; % plotting grid
params.nv           = n;
params.Lv           = 5*pi;
params.detTol       = 1e-2; % incomp. threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condition
params.l            = 1;
params.k            = 0.2;
params.eps          = 5e-2;
params.v0           = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain size
params.Lx           = 2*pi/params.k;
params.L = [params.Lx, params.Lv*2];                                        % domain size
params.dom = [0, 0, params.Lx, 2*params.Lv];                                % domain boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.CFL          = 0.9;
params.T_end        = 40;
params.dt           = 0.002;
params.iplot        = 100; % plot every iplot time steps
params.ihist        = 100;
params.ilog         = 10;
params.dt_hist      = 1;
params.dt_log       = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.bump_transition_width = 0.1*params.Lv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options:
params.case            = 'two_stream';
params.filter          = 'gauss';