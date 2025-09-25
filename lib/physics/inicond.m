function [f0] = inicond(params)


Lv = params.L(2);
eps = params.eps;

if params.case == "landau_damping"
    k = params.k;
    f0 = @(x,v) (1 + eps*cos(x*k)) ./sqrt(2*pi).*exp(-(v-Lv/2).^2/2);
elseif params.case == "two_stream"
    k = params.k;
    v0 = params.v0;
    f0 = @(x,v) (1 + eps*cos(k*x)) ./(2*sqrt(2*pi)).*(exp(-(v-Lv/2-v0).^2/2)+exp(-(v-Lv/2+v0).^2/2));
elseif params.case =="keen_waves"
    f0 = @(x, v) exp(-(v-Lv/2).^2/2) / sqrt(2*pi);
else
    display("Case: "+params.case+ " does not exist!")
end


end