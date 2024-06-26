n=2^8; % 4097
params.eqname = "BGK";
params.Nmap         =n;
params.Nfine        =n;
params.Nsampling    =n;
params.Nplotting    =n;
params.nv           = n;
params.Lv           = 25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condition
params.l            = 1;
params.k            = 0.5;
params.eps          = 5e-2;
params.v0           = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain size
params.Lx           = 5;
params.L = [params.Lx, params.Lv*2];                                        % domain size
params.dom = [0, 0, params.Lx, 2*params.Lv];                                % domain boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.CFL          = 0.9;
params.T_end        = 0.12;
params.dt           = 0.005;
params.iplot        = 100; % plot every iplot time steps
params.ihist        = 100; % only used if dt_hist not set
params.ilog         = 10;  %only used if dt_log not set
params.dt_hist      = 0.01;
params.dt_log       = 0.001;
params.do_plot = 1;
params.dt_constant = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% source term
params.use_sourceterm = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.ideal_gas_constant = 8.3145;  % J *mol-1 K^-1 ideale gas constant
params.relaxation_time=0.01; % particle responds time 1/tau (f-fM) 0.00001
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.bump_transition_width = 0.1*params.Lv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options:
params.case          = 'sod_shock_tube';
params.filter_name          = 'gauss';
params.myplot =  @plot_sod;
params.genVel = @gen_Velocity_Vlasov_Boltzmann;
params.mymeasure = @measure_VP;