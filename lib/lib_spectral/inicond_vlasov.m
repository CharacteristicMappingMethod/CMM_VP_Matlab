
function [phi,f] = inicond_vlasov(params)


switch params.inicond
    case 'landau_damping'
        f = 1/params.l*(1+params.eps*cos(params.k*params.X))/sqrt(2*pi).*exp(-params.V.^2/2);
        phi = cofitxy(give_potential(params,f));  % note that for vlassov we return the distribution function not the velocity vecto
    case 'two_stream'
        f = 1/params.l*(1+params.eps*cos(params.k*params.X))/(2*sqrt(2*pi)).*(exp(-(params.V-params.v0).^2/2)+exp(-(params.V+params.v0).^2/2));
        phi = cofitxy(give_potential(params,f));
    case "sod_shock_tube"
        Rconst = params.ideal_gas_constant;
        % left side:
        rho = 1;
        Umacro = 0;
        p = 1;
        T = p./rho;
        fl = @(x,v) rho./(2*pi*Rconst*T).^0.5 .* exp(-((v)-Umacro ).^2./(2*Rconst*T));

        % right side:
        rho = 0.125;
        Umacro = 0;
        p = 0.1;
        Lx = params.L(1);
        T = p./rho;
        fr = @(x,v) rho./(2*pi*Rconst*T).^0.5 .* exp(-((v)-Umacro ).^2./(2*Rconst*T));
        dx = params.dx;
        bump = @(x) 0.5 + 0.5 * tanh(x/(2*dx));
        %bump = @(x)
        f0 = @(x,v) (1-bump(abs(x-0.5*Lx)-0.2*Lx)).*fl(x,v) + bump(abs(x-0.5*Lx)-0.2*Lx).* fr(x,v);
        f = f0(params.X,params.V);
        phi = zeros(size(params.X));
    otherwise
        error('Initial condition unkown...')
end
end