clear all;
clc;

%just some settings 

% do not change
flags.phi_entire_matrix = 0;
flags.ssfun_normalized = 0;
flags.tau_mult = 1;
flags.tau_new = 0;


%% add all files to path
addpath(genpath('./..'));

%% Settings for running the model
NH = 6;
NV = 6;



load('data/parameters_example','pars'); % parameters without nans
pars1 = pars;
load('data/parameters'); % true parameter set with nans


%% new parameters;
pars1.r = 0.1+ rand(NH,1)*0.2 ;
pars1.M = randi([0 1], NH,NV);
pars1.M = [1 0 0 0 0 0;
    1 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 1 0;
    0 0 0 1 1 0;
    0 0 0 0 0 1];
pars1.beta = pars1.M.*(100+200*rand(NH,NV));
pars1.tau = pars1.M.*(1+5*rand(NH,NV));
pars1.phi = pars1.M.*(0.5e-9+1e-8*rand(NH,NV));
pars1.epsilon = ones(1,NH+NV);
pars1.eta = zeros(NH,NV);
pars1.NE = 100*pars1.M;
pars1.a = 1*rand(NH,NH);
pars1.m = 0.5*ones(NV,1);

pars1.S0 = 3e6*ones(NH,1);
pars1.V0 = 1e6*ones(NV,1);
pars1.NH  = NH;
pars1.NV = NV;

% finally fix the model with given number of boxes.
max_NE = round(max(max(pars1.NE)));
model = SEIVD_diff_NE_diff_debris_abs(NH,NV,max_NE);
model.host_growth = 0;
model.viral_decay = 1;
model.viral_adsorb = 0;
model.lysis_reset = 0;
model.debris_inhib = 0;
model.debris_inhib2 = 0;
model.debris_inhib3 = 0;
model.debris_inhib4 = 0;
model.debris_inhib5 = 0;



% this was used to modify latent periods to be longer, is not used, ignore
% for code review - the multiplication factor was set to 1 anyway.
pars1.tau = pars1.tau*flags.tau_mult;

% a change of variable was done to test inference priors, not needed for code review
pars1.eta(pars1.tau>0) = 1./pars1.tau(pars1.tau>0);

tvec = 0:0.05:10; % for better viz
[t1,S1,V1,D1] = simulate_ode(model,pars1,tvec,pars1.S0,pars1.V0); % initial parameter set

%%

figure(1)
for i=1:model.NH
plot(t1,S1(:,i),'-',LineWidth=2); hold on;%microbes
xlabel('Time (hours)');ylabel('Host (cells/ml)');set(gca,'Yscale','log');
%title(['Points are from averaged experimental datasets, line is from simulated SEIV model' ...
    %' (N_E=10)'])
end

figure(2)
for i=1:model.NV
plot(t1,V1(:,i),'-',LineWidth=2); hold on; %virus
xlabel('Time (hours)');ylabel('Phage (cells/ml)');set(gca,'Yscale','log');
%title('Points are from averaged experimental datasets, line is from simulated SEIV model (N_E=10)')
end



