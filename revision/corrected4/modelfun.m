function  y = modelfun(data,theta,  pars2,mcmcpars,model )
%MODELFUN Summary of this function goes here
%   Detailed explanation goes here

pars_samples = update_pars(pars2,theta,mcmcpars);


%%%%%  this was not included before.
pars_samples.eta(pars_samples.tau>0) = 1./pars_samples.tau(pars_samples.tau>0);

[t, S, V, D, I,E] = simulate_ode(model, pars_samples, data.xdata, pars2.S0, pars2.V0);
y = [S,V];

end

