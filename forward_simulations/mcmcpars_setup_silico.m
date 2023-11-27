function mcmcpars = mcmcpars_setup_new(pars, pars_example, include_pars, flags,model)
% default values - sub-index (which elements for each parameter), limits, priors, and starting values
% default values are calculated from values within pars and pars_example

% growth rate r
tmp_mcmcpars.r.log = 0;
tmp_mcmcpars.r.subid = 1:length(pars.r);
tmp_mcmcpars.r.lims = [0 1];
tmp_mcmcpars.r.startval = pars.r;
tmp_mcmcpars.r.priormu = pars.r;
tmp_mcmcpars.r.priorstd = ones(1,length(pars.r))*0.2;



% burst size beta
tmp_mcmcpars.beta.log = 0;
tmp_mcmcpars.beta.subid = find(pars.M);
tmp_mcmcpars.beta.lims = [0 1000];
tmp_mcmcpars.beta.startval = pars_example.beta(find(pars.M));
tmp_mcmcpars.beta.priormu = pars.beta(find(pars.M));
tmp_mcmcpars.beta.priorstd = ones(1,length(tmp_mcmcpars.beta.priormu))*500;

% adsorption rate phi
    tmp_mcmcpars.phi.log = 1;
    tmp_mcmcpars.phi.subid = find(pars.M);
    tmp_mcmcpars.phi.lims = [1e-12 1e-5];
    tmp_mcmcpars.phi.startval = pars_example.phi(find(pars.M));
    tmp_mcmcpars.phi.priormu = pars.phi(find(pars.M));
    tmp_mcmcpars.phi.priorstd = ones(1,length(tmp_mcmcpars.beta.priormu))*2;


% inverse latent period eta
tmp_mcmcpars.eta.log = 0;
tmp_mcmcpars.eta.subid = find(pars.M);

if flags.tau_new == 0
tmp_mcmcpars.eta.lims = [0 10];
elseif flags.tau_new == 1
  tmp_mcmcpars.eta.lims = [0.4 10];
end

tmp_mcmcpars.eta.startval = pars_example.eta(find(pars.M));
tmp_mcmcpars.eta.priormu = pars.eta(find(pars.M));
tmp_mcmcpars.eta.priorstd = std(pars.eta(find(pars.M)),'omitnan')*ones(1,length(find(pars.M)));



% latent period
tmp_mcmcpars.tau.log = 0;
tmp_mcmcpars.tau.subid = find(pars.M);
tmp_mcmcpars.tau.lims = [0.1 10];
tmp_mcmcpars.tau.startval = 1./pars_example.eta(find(pars.M));
tmp_mcmcpars.tau.priormu = tmp_mcmcpars.tau.startval;
tmp_mcmcpars.tau.priorstd = ones(1,length(tmp_mcmcpars.tau.priormu))*10;

% number of compartments
% will do later


tmp_mcmcpars.m.log = 1;
tmp_mcmcpars.m.subid = 1:length(pars.m);
tmp_mcmcpars.m.lims = [1e-5 1]; % 1000days - 1hr
tmp_mcmcpars.m.startval = pars.m;
tmp_mcmcpars.m.priormu = pars.m;
tmp_mcmcpars.m.priorstd =1*ones(1,length(pars.m));

% only include these parameters
if ~exist('include_pars','var') || isempty(include_pars) % everything by default
   include_pars = fieldnames(tmp_mcmcpars);
end
for i = 1:length(include_pars)
   mcmcpars.(include_pars{i}) = tmp_mcmcpars.(include_pars{i}); 
end

end


%% post processing and plots (subpar version now)

