clc;
clear all;
tic

addpath(genpath('./../..'))
load('./../combined_posteriors.mat');
load('./../../data/triplicate_data.mat');
color_ofthe_fit = [1 0 0]*0.5;
color_ofthe_fill = [0.95 0 0];
transparency = 0.25;
linewidth = 2;

tvec = 0:0.05:15.75; % for better viz
%% new parameters

pars2 = update_pars(pars1,pars_from_dist(chain_stored4(5001:end,:)),mcmcpars);
pars2.epsilon = ones(1,10);
pars2.prob = [0 0 0 0 0]';

   pars2.phi = 1.0e-07 *[         0    0.6000         0         0      0;
                             0.18    0.8000    0.25        0         0;
                                  0         0    0.900         0         0;
                                  0         0         0    0.6   0.1285;
                                  0         0         0    0.6   0.22];


pars2.tau =[0 3 0 0 0;
            1.7 2.7 2.3 0 0;
            0 0 2 0 0 ;
            0 0 0 1.8 4.7;
            0 0 0 2.3 2];
pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);

% r
pars2.r = [0.18,0.25,0.3,0.68,0.52]' ;


% beta

pars2.beta = [0  1.7231         0         0         0;
  200.7512  205.9496  100.1492         0         0;
         0         0   20.7017         0         0;
         0         0         0  522.0549  60.2599;
         0         0         0  485.1209  50.9918];




pars2.Dc = 5e6;
pars2.Dc2 = 6.13e6;
pars2.Dc3 =12e6;
pars2.Dc4 = 19.3e5;
pars2.Dc5 = 16.0e5;


pars2.NE = 200*pars2.M; %otherwise this will take a very long time

max_NE = round(max(max(pars2.NE)));
model = SEIVD_diff_NE_diff_debris_abs(5,5,max_NE);
model.host_growth = 0;
model.viral_decay = 0;
model.viral_adsorb = 0;
model.lysis_reset = 0;
model.debris_inhib = 2;
model.debris_inhib2 = 2;
model.debris_inhib3 = 2;
model.debris_inhib4 = 2;
model.debris_inhib5 = 2;

model.diff_beta = 0;


model.name = 'SEIVD-diffabs';

seed = 3500;

while exist('revised'+string(seed)+'.mat','file') == 2
    seed = seed + 1;
end




[t_median,S_median,V_median,D_median,I_median,E_median] =  simulate_ode(model,pars2,tvec,pars2.S0,pars2.V0); % mcmc parameter set

%% inference part

mcmcoptions.nsimu = 10000;
transient_id = 30;
lambda = 0;
include_pars = {'beta','phi','tau','r','Dc','Dc2','Dc3','Dc4','Dc5'};

% can vectorize for all the 10 time series,
% then ssfun needsto be a 10 element vector -- now set to sum.
% variance controls width of deviation -- so acceptable rate
% N0 is the var of the error variance -- so helps to tune that -- so
% changes the sticking behavior -- more means less accurate -- so explore
% sigma of error chain more. 

mcmcmodel.sigma2 = 200; % (initial) error variance from residuals of the lsq fit
mcmcmodel.S20 = mcmcmodel.sigma2;
mcmcmodel.N0 =50;
mcmcoptions.updatesigma = 1;



mcmcoptions.method  = 'dram';
mcmcpars = mcmcpars_setup_new(pars2,pars2,include_pars,flags,model);


mcmcparam = mcmcpars2param(mcmcpars);

%% priors and initial starting positions
% priors are truncated gaussians
% name -- start -- min -- max -- mean -- std


theta_optimized(1:9) = [1.81,118.2,60.5,20,8.246,99.2,437,93,413];
theta_optimized(10:18) = log([5.14e-8,1.45e-8,9.57e-8,2.5e-8,1.227e-7,1.56e-7, 6.46e-8,1.31e-7,8.02e-8 ])/log(10);
theta_optimized(19:27)  = [1.749,1.89,2.19,1,1.9,1.47,2.14,1.42,1.97]; 
theta_optimized(28:32) = [0.19,0.245,0.22,0.28,0.25];
theta_optimized(33:37) = log([5e6 6.13e6 12e6  19.3e5 16.0e5])/log(10);


std = [2.02,68.8,50,20,21.9,43.8,79.5,50,86, 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,   0.56,0.18,0.19,2,0.4,0.13,0.13,0.12,0.09,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5]

%%% initial values set
for i = 1:37
mcmcparam{1,i}{1,2} = theta_optimized(i); %initial value
mcmcparam{1,i}{1,5} = theta_optimized(i); %mean
mcmcparam{1,i}{1,6} = std(i);
end



%%% setting up means of priors.



%% mcmc run

mcmcmodel.ssfun = @(theta,data) ssfun(theta,data,pars2,mcmcpars,model,lambda); 
[mcmcresults, chain, s2chain]= mcmcrun(mcmcmodel,data,mcmcparam,mcmcoptions);

%%
pars_afterinf = update_pars(pars2,median(chain(transient_id:end, :)),mcmcpars);
[t_after,S_after,V_after,D_after,I_after,E_after] =  simulate_ode(model,pars_afterinf,tvec,pars2.S0,pars2.V0); % mcmc parameter set

save('from_os_prior_16.mat')
