

function [fhist,params] =  vlasov1d(params)

if ~isfield(params, 'dt_log')
    params.dt_log = params.ilog * params.dt;
end
if ~isfield(params, 'dt_hist')
    params.dt_hist = params.ihist * params.dt;
end
if ~isfield(params, 'eqname')
    params.eqname = "VP";
end

% convert params:
params.y            = params.v;
params.Ly           = 2*params.Lv;
params.L = [params.Lx, params.Lv*2];
params.ny           = params.nv;

time = 0;
it   = 1;
it_hist = 1;
it_log   = 1;
% initialize grid, wavenumbers, etc
params = geometry_wavenumbers(params);
params.dv = params.dy;
params.V = params.Y;

% compute velocity field for periodic continuation
params.kv = params.ky;
[params.v_periodic,params.grid.Sigma] = velocity_periodicfication_spectral(params,params.bump_transition_width);
[~,params.VP] = meshgrid(params.x,params.v_periodic);
% initial condition
[phi,f] = inicond_vlasov(params);

% compute inital electric field
fk = fft2(f);
tic();
[~, Efield] = nonlinear_vlasov(params,time,fk);
params.tcpu_op =toc();

% meassure first time
params.time_log(it_log) = time;
[params.Mass(it_log),params.Momentum(it_log),params.Epot(it_log),params.Ekin(it_log),params.Etot(it_log),params.L2norm(it_log),params.spectr_fft(it_log,:),rk] = measure(f',Efield, params.grid);
params.dM_rel(it_log) = abs(params.Mass(it_log) - params.Mass(1))/abs(params.Mass(1));
params.dP_rel(it_log) = abs(params.Momentum(it_log));
params.dE_rel(it_log) = abs(params.Etot(it_log) - params.Etot(1))/abs(params.Etot(1));
params.dL2_rel(it_log) = abs(params.L2norm(it_log) - params.L2norm(1))/abs(params.L2norm(1));
params.elapsed_time_per_iter(it_log)=0;
params.fhist(:,:,it_hist) = f;
params.time_hist(it_hist) = time;

figure(4)
set(gcf,'Position',[100 100 1800 1200])

% --------------
dt_hist = params.dt_hist;
it_hist = 1;
dt_log = params.dt_log;
it_log = 1;
% --------------

hist_log = "";
ie = 0;
dt_mean = 0;
tic();
while (time < params.T_end)
    % checks cfl criteria and condition for logging and backupsaves
    dt = give_time_step(params, Efield, time);

    % advance fluid in time
    [fk, Efield] = RK3(params,time,dt,fk);
    %[fk, Efield] = euler(params,time,dt,fk);
    %[fk, Efield] = make_step(params,time,dt,fk);
    % update the time step
    time = time + dt;

    % calculate the mean dt until next logging
    dt_mean = dt_mean + dt;

    ie = ie+1;
    it = it+1;
    if (abs(time - (it_hist)*dt_hist)<1e-15)
        it_hist = it_hist+1;
        f = cofitxy(fk);
        params.fhist(:,:,it_hist) = f;
        params.time_hist(it_hist) = time;
        hist_log = " (config saved)";
    end

    % we also save the logs when we backup
    if (abs(time - (it_log)*dt_log)<1e-15 || abs(time - (it_hist)*dt_hist)<1e-15)
        it_log = it_log + 1;
        params.time_log(it_log) = time;
        f = cofitxy(fk);
        params.elapsed_time_tot(it_log) = toc();

        params.dt_mean(it_log) = dt_mean/ie;
        params.elapsed_time_per_iter(it_log) =(params.elapsed_time_tot(it_log)-params.elapsed_time_tot(it_log-1))/ie;
        params.tcpu_rk3(it_log) = params.elapsed_time_per_iter(it_log) - 3*params.tcpu_op;
        [params.Mass(it_log),params.Momentum(it_log),params.Epot(it_log),params.Ekin(it_log),params.Etot(it_log),params.L2norm(it_log),params.spectr_fft(it_log,:),rk] = measure(f',Efield, params.grid);
        params.dM_rel(it_log) = abs(params.Mass(it_log) - params.Mass(1))/abs(params.Mass(1));
        params.dP_rel(it_log) = abs(params.Momentum(it_log));
        params.dE_rel(it_log) = abs(params.Etot(it_log) - params.Etot(1))/abs(params.Etot(1));
        params.dL2_rel(it_log) = abs(params.L2norm(it_log) - params.L2norm(1))/abs(params.L2norm(1));
        fprintf("it: %d time: %.2f dt: %.1e tcpu/iter: %.1e dE: %.1e dP: %.1e dM: %.1e dL2: %.1e" ...
            +hist_log+ "\n",it, time,dt_mean/ie, params.elapsed_time_per_iter(it_log), params.dE_rel(it_log), params.dP_rel(it_log), params.dM_rel(it_log), params.dL2_rel(it_log));
        hist_log = "";
        ie = 0;
        dt_mean = 0;
    end
    % plot if time to do so
    if (mod(it,params.iplot)==0)
        clf;
        if params.eqname=="BGK"
             plot_BGK(params,fk,Efield)
        else
             plot_VP(params,fk,Efield)
        end
    end
end
params.freqs = rk;
fhist = params.fhist;
end




function [fk_new, Efield] = RK2(params,time,dt,fk)

% compute non-linear terms
[nlk, Efield] = nonlinear_vlasov(params,time,fk);

% integrating factor:
vis = 1;%cal_vis_2d(dt);
fk_new = vis.*(fk + dt*nlk );
fk_new = dealias(params,fk_new);

% compute non-linear terms at new time level
[nlk2, Efield] = nonlinear_vlasov(params,time,fk_new);

% advance to final velocity
fk_new = vis.*fk + 0.5*dt*(nlk.*vis + nlk2);

% Dealiasing
fk_new = dealias(params,fk_new);
end

function [fk_new, Efield] = RK3(params,time,dt,fk)

% Compute non-linear terms at the beginning
[nlk1, Efield] = nonlinear_vlasov(params, time, fk);

% Calculate intermediate stage values
k1 = dt * nlk1;

k2_temp = nonlinear_vlasov(params, time + 0.5 * dt, fk + 0.5 * k1);
k2 = dt * k2_temp;

k3_temp = nonlinear_vlasov(params, time + dt, fk - k1 + 2 * k2);
k3 = dt * k3_temp;

% Update distribution function using weighted average of stages
fk_new = fk + (1/6) * (k1 + 4 * k2 + k3);

% Dealiasing
fk_new = dealias(params, fk_new);
end


function [fk_new, Efield] = euler(params,time,dt,fk)

% Compute non-linear terms at the beginning
[nlk1, Efield] = nonlinear_vlasov(params, time, fk);

% Calculate intermediate stage values
k1 = dt * nlk1;

fk_new = fk + k1;

% Dealiasing
fk_new = dealias(params, fk_new);
end

% advance fluid in time
function   [fknew, Efield] = make_step(params,time,dt,fk)

    tau = params.relaxation_time*0.5;
    time_new = time + dt;
    % advance fluid in time
    Nsteps = 5;
    
    for n = 1:Nsteps
        fkm = fk;
    [fk, Efield] =euler(params,time,tau,fkm);
    time = time + tau;
    end
    assert(time_new-time>0,"N*tau > dt")
    fknew = fk + (time_new - time)*(fk-fkm)/tau;
    fk_new = dealias(params, fk_new);
end



function [nlk, Efield]=nonlinear_vlasov(params,time,fk)
%% computes the non-linear terms + penalization in Fourier space.
% you can also add explicit laplacian (explicit diffusion) if you
% want to use higher order time stepping (maybe RK4) without
% integrating factor

f = cofitxy(fk);

% init source term
ST = 0;
if params.eqname == "BGK"
        Efield = zeros(size(f,1),1);
    dphi_x = Efield;

    Rconst = params.ideal_gas_constant; % J *mol-1 K^-1 ideale gas constant
    fm = calc_maxwellian_spectral(f,Rconst,params.V);

    tau = params.relaxation_time;
    ST = 1/tau * (fm- f); % source term
else
    phi_k = give_potential(params,f);
    dphi_x =ifft(1i*params.kx'.*phi_k,'symmetric');
    Efield = - dphi_x;
end

if params.v_periodic
    [u(:,:,1), u(:,:,2)] = meshgrid(params.v_periodic,dphi_x);
else
    [u(:,:,1), u(:,:,2)] = meshgrid(params.v,dphi_x);
end


nlk = fft2(u(:,:,1).*cofitxy(1i*params.Kx.* fk) + u(:,:,2).*cofitxy(1i*params.Ky.*fk) + ST);

% delete aliased modes
nlk = dealias(params,nlk);

% emergency brake
if (max(max(max(abs(u))))>1e3)
    error('diverged..')
end
end