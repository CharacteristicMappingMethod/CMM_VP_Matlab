close all
clear all
clc

addpath(genpath('./lib/'),genpath('./params/'))
INIT_DEFAULTS
params.ilog = 1;

%---------------------------
%% here, all case-specific options are hidden:
%---------------------------
%PARAMS_non_lin_landau_damping;
%PARAMS_landau_damping;
%PARAMS_two_stream;
PARAMS_BGK_sod_shock_tube;

params = convert_cmm2fourier(params);
%% start simulation
tic()
[fhist,params] = vlasov1d(params);
toc()
folder = "./data/"+params.case+"/Kn"+num2str(params.relaxation_time)+"_N"+num2str(params.Nfine);
mkdir(folder)
save(folder+"/params","params")
save(folder+"/fhist","fhist")