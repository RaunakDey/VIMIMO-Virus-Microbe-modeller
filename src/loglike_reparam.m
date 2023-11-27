function [ss, err_vec, err_id, err_labels, err_legends] = loglike_reparam(theta, data, pars, mcmcpars, model, lambda, theta_start, theta_std)
% log likelihood function
% for updating the pars struct: pars, theta, mcmcpars
% for simulating the ode: model, pars, S0, V0, and time vector from data
% for computing ss: ydata from data and model output from ode

if ~exist('lambda','var') || isempty(lambda)
    lambda = 0;
end

% grab qpcr data
tvec = data.xdata;
ydata = data.ydata;



theta = (theta + theta_start).*theta_std;

% simulate ode with updated pars
if ~isempty(theta)
    pars = update_pars(pars,theta,mcmcpars);
end

%%%%%  this was not included before.
pars.eta(pars.tau>0) = 1./pars.tau(pars.tau>0);


S0 = pars.S0;
V0 = pars.V0;
[~,S,V] = simulate_ode(model,pars,tvec,S0,V0);
ymodel = [S V];
ymodel(ymodel<0) = 0; % get rid of negative densities

% compute error for each timepoint
log_error_square = (log(ymodel)-log(ydata)).^2;

% total error, not normalized
ss = sum(log_error_square(~isnan(log_error_square)));

% total error with lasso on phi matrix
ss = ss + lambda*norm(pars.phi,1);



%% more detailed error quantities

% host vs phage
les_host = log_error_square(:,data.id.S);
les_phage = log_error_square(:,data.id.V);

% totals
err = ss;

% by channel
err_channel_host = sum(les_host,'omitnan');
err_channel_phage = sum(les_phage,'omitnan');

% phi matrix
err_phi = lambda*norm(pars.phi,1);

% combine
err_vec = [err err_channel_host err_channel_phage err_phi];
err_id = [1 repmat(2,[1 length(err_channel_host)]) repmat(3,[1 length(err_channel_phage)]) 4];
err_labels = ["total","host channels","phage channels","phi L1 norm"];
err_legends = {{},...
    cellstr(data.labels.host(~all(isnan(les_host)))),...
    cellstr(data.labels.phage(~all(isnan(les_phage)))),...
    {}};

end