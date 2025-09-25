n=2^8; % 4097
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.Nmap           = n; % grid dimension of the submaps n x n
params.Nfine        =2^8; % grid for measuring and saves of f
params.Nsampling    =2^10; % upsampling grid of veloctiy
params.Nplotting    =2^8; % plotting grid
params.nv           = n;
params.Lv           = 6;
params.detTol       = 1e-2; % incomp. threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condition
params.l            = 1;
params.kDr            = 0.26;
params.eps          = 5e-2;
params.v0           = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain size
params.Lx           = 2*pi/params.kDr;
params.L = [params.Lx, params.Lv*2];                                        % domain size
params.dom = [0, 0, params.Lx, 2*params.Lv];                                % domain boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.CFL          = 1;
params.T_end        = 1000;
params.dt0           = 0.25;
params.dt_constant  = 0.25;
params.iplot        = 100; % plot every iplot time steps
%params.ihist        = 100;
%params.ilog         = 10;
params.dt_hist      = 100;
params.dt_log       = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.bump_transition_width = 0.1*params.Lv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options:
params.case            = 'keen_waves';
params.filter          = 'gauss';
params.do_plot=1;

% Keen waves parameters from paper
params.wDr = 0.37;                          % Drive frequency
params.aDr = 0.2;                           % Drive amplitude (canonical)
params.TDr = 100;                           % Drive duration

% Drive timing parameters
params.t0 = 0;                              % Start time
params.tL = 69;                             % Left ramp start
params.twL = 20;                            % Left ramp width
params.twR = 20;                            % Right ramp width
params.tR = 207 + params.TDr;               % Right ramp start


% Initial condition: spatially uniform Maxwellian

% External drive field: E_Pond(x,t) = a_Dr * k_Dr * a(t) * sin(k_Dr * x - w_Dr * t)
g = @(t) 0.5 * (tanh((t - params.tL) / params.twL) - tanh((t - params.tR) / params.twR));
a = @(t) (g(t) - g(params.t0)) / (1 - g(params.t0));
params.E_ext = @(x, t) params.aDr * params.kDr * a(t) .* sin(params.kDr * x - params.wDr * t);