function [BMjet, u1j_m, u1j_mm, u2j_m, u2j_mm, t, tm ,tmm] =  initEulerB0_L3(params,BMjet, evalVOp, mgsz, vgsz, samplinggrid)

sgsz = samplinggrid.size;
t = 0;
dom = params.dom;

[u1h, u2h, f, M1g, M2g, Jac] = evalVOp({BMjet}, samplinggrid,t);


u1j = dataF2H(params,u1h, sgsz, vgsz);
u2j = dataF2H(params,u2h, sgsz, vgsz);

dt = compute_time_step(params, f, u1j ,u2j , t);

u1j_f = @(t) u1j;
u2j_f = @(t) u2j;


bmjet = RK4_jetstages(params, vgsz, u1j_f, u2j_f, t+dt, -dt, mgsz);

BMjet_t = HMapCompose(params, mgsz, BMjet, mgsz, bmjet);


u1j_m = u1j;
u2j_m = u2j;
tm = t;
t = t+dt;


% next step eval an improve the Euler step, update map
[u1h, u2h, f, M1g, M2g, Jac] = evalVOp({BMjet_t}, samplinggrid,t);

u1j = dataF2H(params,u1h, sgsz, vgsz);
u2j = dataF2H(params,u2h, sgsz, vgsz);

u1j_f = defExt_L2(u1j, u1j_m, t, tm);
u2j_f = defExt_L2(u2j, u2j_m, t, tm);

u1j_mm =u1j;
u2j_mm =u2j;

bmjet = RK4_jetstages(params, vgsz, u1j_f, u2j_f, t, -dt, mgsz);

BMjet = HMapCompose(params, mgsz, BMjet, mgsz, bmjet);


% update velocity again with new map and move to next step
[u1h, u2h, f, M1g, M2g, Jac] = evalVOp({BMjet}, samplinggrid,t);

u1j = dataF2H(params, u1h, sgsz, vgsz);
u2j = dataF2H(params, u2h, sgsz, vgsz);

% Heuns on velocity
u1j = 0.5*(u1j + u1j_mm);
u2j = 0.5*(u2j + u2j_mm);

u1j_f = defExt_L2(u1j, u1j_m, t, tm);
u2j_f = defExt_L2(u2j, u2j_m, t, tm);

bmjet = RK4_jetstages(params, vgsz, u1j_f, u2j_f, t+dt, -dt, mgsz);

BMjet = HMapCompose(params, mgsz, BMjet, mgsz, bmjet);


u1j_mm = u1j_m;
u2j_mm = u2j_m;
tmm =tm;
u1j_m = u1j;
u2j_m = u2j;
tm = t;
t = t+dt;

return