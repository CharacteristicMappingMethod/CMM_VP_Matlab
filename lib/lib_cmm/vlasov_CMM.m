function [params,f,XB,YB] = vlasov_main(params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main function of the Vlasov Poisson CMM solver
% Authors: P.Krah and  Xi-Yuan Yin
% Date: 27.07.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------
% Data management
%----------------------------------------
if(~exist('clim', 'builtin'))
    clim = @(varargin) caxis(varargin{:});
end
if isfield(params,"dat_dir")
    dat_dir = params.dat_dir;
else
    dat_dir = "./data/";
end

if isfield(params,"load_existing")
    load_existing = params.load_existing;
else
    load_existing = 0;
end

if ~isfield(params, 'do_plot')
    do_plot = 0;
else
    do_plot = params.do_plot;
end


%----------------------------------------
% PARAMETER Initialization
%----------------------------------------
f0 = inicond(params);                 % initial condition as a executable function
T_end = params.T_end;                 % maximal time
dom = params.dom;                     % domain boundaries
CFL = 1;

% filter
if params.filter == "truncate"
    genKFOp = @(R) setupKelvinFilteredOp_truncate(R);
else
    ptrunc = 8; % suggested p=10 for drop from 1e-1 to 1e-8 within 20% of R
    genKFOp = @(R) setupKelvinFilteredOp_gauss_2p(R, ptrunc);
end

%----------------------------------------
% GRID Initialization
%----------------------------------------
mgsz = params.Nmap * [1 1];         % grid points of the map
vgsz = params.Nsampling *[1 1];     % grid of velocity (here same as sampling grid)
sgsz = params.Nsampling * [1 1];    % grid for sampling f
pgsz = params.Nplotting * [1 1];    % grid for plotting f
fgsz = params.Nfine * [1 1];        % fine grid for saving;

params.grid = set_grids(params,{mgsz,vgsz,sgsz,pgsz,fgsz}, ...
    ["map","velocity","sampling","plotting","fine"]);

NsavesF = 5;    % number of large saves
RPlot = 128;

%R = 24; % Mapgrid/2
R = params.Nmap/2;
if isfield(params, 'dt_constant') && isfield(params, 'dt0')
    dt0 = params.dt0;
    if isfield(params, 'ilog') params.dt_log = params.ilog *dt0; end
    if isfield(params, 'ihist') params.dt_hist = params.ihist *dt0; end
else
    dt0 = 0.5*max(dom([3 4])-dom([1 2]))/max(mgsz)/CFL;
    params.dt0 = dt0;
end

if ~isfield(params, 'detTol')
    display("default remap threshold used")
    detTol = 1e-3;
else
    detTol = params.detTol;
end
dtThres = 5e-3*dt0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD the data if already exists!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir = dat_dir + params.case + "_Map" + params.Nmap+ "_Fine" + params.Nfine + "_Nt"+ceil(params.T_end/dt0);
mkdir(dir);
fname = dir + "/params.mat";
if exist(fname,"file") && load_existing
    dat = load(fname);
    params = dat.params;
    f = params.fhist_fine(:,:,end);

    XB = params.map_end(:,:,1);
    YB = params.map_end(:,:,2);
    f  = params.f_end;
    display("loading simulation results: " + fname)
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xgf, ygf] = genMeshPoints(dom, fgsz);
f = f0(xgf, ygf);

params.time_hist(1) = 0;
params.fhist_fine(:,:,1) =  f;

clear xgf ygf f B

evalVOp = @(BMList, sgsz) genVOp(params,f0, dom, mgsz, genKFOp(R), BMList, sgsz);
evalVOp_plot = @(BMList, sgsz) genVOp(params,f0, dom, mgsz, genKFOp(RPlot), BMList, sgsz);


% plotting grid also need to be periodic, if really want to plot periodic functions properly, wrap
% all the plotting functions
[xgp, ygp] = genMeshPoints(dom, pgsz);

BMjet = zeros(prod(mgsz), 4, 2);
BMList = [];


rmtime = [];



res = [15 15];


t = 0;
nmaps = 0;
nstep = 1;



[u1h, u2h, f] = evalVOp({BMjet}, params.grid.plotting);
fh = fft2(f);
u1 = ifft2(u1h, 'symmetric');
u2 = ifft2(u2h, 'symmetric');



[u1h, u2h, f, XB, YB, JB, Efield] = evalVOp([BMList; {BMjet}], params.grid.sampling);
params.time_log(1) = t;
il = 1;
[params.Mass(il),params.Momentum(il),params.Epot(il),params.Ekin(il),params.Etot(il),params.L2norm(il),~,~,params.Emodes(il,:)] = measure(f,Efield,params.grid.sampling);
params.dM_rel(il) = abs(params.Mass(il) - params.Mass(1))/abs(params.Mass(1));
params.dE_rel(il) = abs(params.Etot(il) - params.Etot(1))/abs(params.Etot(1));
params.dL2_rel(il) = abs(params.L2norm(il) - params.L2norm(1))/abs(params.L2norm(1));
params.nmaps(il) = nmaps;


[u1h, u2h, f, XB, YB, JB, Efield] = evalVOp([BMList; {BMjet}], params.grid.fine);

[Mf,Pf,~,~,Etotf,L2Normf] = measure(f,Efield,params.grid.fine);

%% plotting first frame

nxs = ceil(pgsz(2)/24);
nys = ceil(pgsz(1)/24);


if (~isfield(params, 'dt_log'))
    params.dt_log = params.ilog * dt0;
end
if (~isfield(params, 'dt_hist'))
    params.dt_hist = params.ihist * dt0;
end
% --------------
dt_hist = params.dt_hist;
dt_log = params.dt_log;
% --------------

if do_plot
    fig = figure(6);
    fig.Position = [100 100 1800 1200];
end
%% start sim
[BMjet, u1j_m, u1j_mm, u2j_m, u2j_mm, t, tm ,tmm] =  initEulerB0_L3(params,BMjet, evalVOp, mgsz, vgsz, params.grid.sampling);

ifbreak = false;
hist_log ='';
elapse_time = 0;
elapse_tot = 0;
ik = 0;
il = 0;

% after the initEuler it might happen that t has already run several
% dt_logs
% so we have to set the counter accordingly
it_log = floor(t/dt_log)+1;
it_hist = floor(t/dt_hist)+1;

while (norm(t-params.T_end)>1e-12)
    il = il + 1;
    ik = ik + 1;
    %---------------------------------------------------------
    % evaluate the poisson equation to compute the velocities
    %---------------------------------------------------------
    tic();
    [u1h, u2h, f, XB, YB, JB, Efield] = evalVOp([BMList; {BMjet}], params.grid.sampling);
    u1j = dataF2H(params,u1h, sgsz, vgsz);
    u2j = dataF2H(params,u2h, sgsz, vgsz);

    % measure time for poisson equation solve
    elapse_time(1,ik) = toc();

    % -----------------
    % make a time step
    % -----------------
    dt = compute_time_step(params, f,u1j ,u2j , t);

    u1j_f = defExt_L3(u1j, u1j_m, u1j_mm, t, tm, tmm);
    u2j_f = defExt_L3(u2j, u2j_m, u2j_mm, t, tm, tmm);

    [bmjet, bM] = RK4_jetstages(params, vgsz, u1j_f, u2j_f, t+dt, -dt, mgsz);
    BMjet = HMapCompose(params, mgsz, BMjet, mgsz, bmjet);

    % measure time for time integration
    elapse_time(2,ik) = toc()-elapse_time(1,ik);

    %%%%%%%%%%%%%%
    % Remapping ?
    %%%%%%%%%%%%%%
    Mdet = (1+BMjet(:,2,1)).*(1+BMjet(:,3,2)) - BMjet(:,3,1).*BMjet(:,2,2);
    maxLogDet = max(abs(log(Mdet)));

    if maxLogDet > detTol

        BMList = [BMList; {BMjet}];
        rmtime = [rmtime; t];

        BMjet = 0.*BMjet;
        nmaps = nmaps +1;

    end
    % measure time for remapping
    elapse_time(3,ik) = toc()-elapse_time(2,ik);
    elapse_time(4,ik) = toc();

    u1j_mm = u1j_m;
    u2j_mm = u2j_m;
    tmm =tm;
    u1j_m = u1j;
    u2j_m = u2j;
    tm = t;
    t = t+dt;
    %%%%%%%%%%%%
    % SAVING
    %%%%%%%%%%%%
    if (abs(t - (it_hist)*dt_hist)<1e-12)
        it_hist = it_hist + 1;
        [u1h, u2h, f, XB, YB, JB, Efield] = evalVOp([BMList; {BMjet}], params.grid.fine);
        [params.rk,params.spectr_fft(it_hist,:)] = spectrum(f);

        params.time_hist(it_hist) = t;
        params.backmap(it_hist).BMList = BMList;
        params.backmap(it_hist).BMjet = BMjet;
        params.fhist_fine(:,:,it_hist) =  f;
        hist_log =' (config saved)';
    end
    %%%%%%%%%%%%
    % LOGGING
    %%%%%%%%%%%%
    if (abs(t - (it_log)*dt_log)<1e-12 || abs(t - (it_hist)*dt_hist)<1e-12)
        [~, ~, f, ~, ~, ~, Efield] = evalVOp([BMList; {BMjet}], params.grid.sampling);
        it_log = it_log + 1;
        params.time_log(it_log) = t;
        [params.Mass(it_log),params.Momentum(it_log),params.Epot(it_log),params.Ekin(it_log),params.Etot(it_log),params.L2norm(it_log),~,~,params.Emodes(it_log,:)] = measure(f,Efield,params.grid.sampling);
        params.dM_rel(it_log) = abs(params.Mass(it_log) - params.Mass(1))/abs(params.Mass(1));
        params.dE_rel(it_log) = abs(params.Etot(it_log) - params.Etot(1))/abs(params.Etot(1));
        params.dL2_rel(it_log) = abs(params.L2norm(it_log) - params.L2norm(1))/abs(params.L2norm(1));
        params.nmaps(it_log) = nmaps;
        params.tcpu_op(it_log) = sum(elapse_time(1,1:ik))/ik;    % average time
        params.tcpu_rk3(it_log) = sum(elapse_time(2,1:ik))/ik;   % average time
        params.tcpu_remap(it_log) = sum(elapse_time(3,1:ik))/ik; % average time
        average_time_per_iter = toc()/ik;
        elapse_tot = elapse_tot + toc();
        params.elapse_time_tot(it_log) = elapse_tot;       % total elapsed time

        ik = 0; % reset ik
        disp(['t = ' num2str(t) ' tcpu(/iter) = ' num2str(toc()) ' (' num2str(average_time_per_iter) ') '...
            '; nmaps = ' num2str(nmaps) ' dM = ' num2str(params.dM_rel(it_log)) ...
            ' dP = ' num2str(params.Momentum(it_log)) ' dE = ' num2str(params.dE_rel(it_log)) ...
            ' dL2 = ' num2str(params.dL2_rel(it_log)) hist_log ]);
        hist_log ='';
    end
    %%%%%%%%%%%%
    % PLOTTING
    %%%%%%%%%%%%
    if (mod(il,params.iplot)==0)
        BMjet_t = HMapCompose(params, mgsz, BMjet, mgsz, bmjet);
        [u1h, u2h, f, XB, YB, JB,Efield] = evalVOp_plot([BMList; {BMjet_t}], params.grid.plotting);

        fh = fft2(f);

        u1 = ifft2(u1h, 'symmetric');
        u2 = ifft2(u2h, 'symmetric');

        clf;
        subplot(2,4,1:2)
        pc(xgp,ygp,f)
        title("$f(x,v,t)$")
        colorbar
        subplot(2,4,3:4)
        U1 = u1;
        U2 = u2;
        imu = pcolor(xgp, ygp, sqrt(U1.^2 + U2.^2));
        title("$f(x,v,t)$")
        shading flat
        imu.FaceColor = 'interp';
        set(gca, 'YDir', 'normal');
        hold on

        uqv = quiver(xgp(1:nys:end, 1:nxs:end), ygp(1:nys:end, 1:nxs:end), u1(1:nys:end, 1:nxs:end), u2(1:nys:end, 1:nxs:end), 'r');
        subplot(2,4,5)
        semilogy(params.time_log,params.dM_rel,'-+')
        xlabel("time $t$")
        title('dM')
        subplot(2,4,6)
        semilogy(params.time_log,params.Momentum,':x')
        title('dP')
        xlabel("time $t$")
        subplot(2,4,7)
        semilogy(params.time_log,params.dE_rel,'--.')
        title('dE')
        xlabel("time $t$")
        subplot(2,4,8)
        semilogy(params.time_log,params.dL2_rel,'--.')
        title('dL2')
        xlabel("time $t$")
        pause(0.0001)
    end

    if any(isnan(BMjet), "all")
        warning("have nan of inf values, skip to finish");
        ifbreak = true;
        nframe = cFrame;
        break;
    end


    if ifbreak
        break;
    end

end

if do_plot
    saveas(fig, datapath+'/dist_sim_stats.eps');
end
if (norm(t-params.T_end)>1e-12)
    warning("\n!!!!!!!!!!!!!!!!!!!!!\nEnd time not reached!\n!!!!!!!!!!!!!!!!!!!!!\n");
    keyboard
end

BMjet_t = HMapCompose(params, mgsz, BMjet, mgsz, bmjet);
[u1h, u2h, f, XB, YB, JB,Efield] = evalVOp([BMList; {BMjet}], params.grid.fine);

[Mf1,Pf1,~,~,Etotf1,L2Norm1] = measure(f,Efield,params.grid.fine);


params.dMtot = abs(Mf-Mf1)/abs(Mf);
params.dPtot = abs(Pf1);
params.dEtot = abs(Etotf-Etotf1)/abs(Etotf);
params.dL2tot = abs(L2Normf-L2Norm1)/abs(L2Normf);
params.map_end(:,:,1) = XB;
params.map_end(:,:,2) = YB;
params.f_end = f;
save('-v7.3', dir +'/params.mat', 'params');
%[params.Mass(il),params.Momentum(il),params.Epot(il),params.Ekin(il),params.Etot(il)] = measure(f,Efield,params.sample_grid.dx,params.sample_grid.dv,params.sample_grid.Vgrid);

end

function [u1h, u2h, f, M1g, M2g, Jac, Efield] = genVOp(params,f0, dom, mgsz, deKFOp, BMList, grid)

xgs = grid.Xgrid;
vgs = grid.Vgrid;


[M1g, M2g, Jac] = fullMapEval(dom, mgsz, BMList, xgs, vgs);
f = f0(M1g, M2g);

[dphi_xh,kx,kv] =vPoisson(f, deKFOp, grid.size, dom, grid.dv);
v_periodic = grid.v_periodic;

dphi_xh = dphi_xh + compute_external_Efield(params, params.grids(1).x, params.time + dt);
dphi_x = reshape(ifft(dphi_xh,"symmetric"),1,[]);
[u2, u1] = meshgrid(dphi_x,v_periodic);

u1h = fft2(-u1);
u2h = fft2(-u2);
Efield = - dphi_x;
end
