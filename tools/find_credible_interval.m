function [S_min,S_max,V_min,V_max,S_median,V_median] = find_credible_interval(chain,transient_id_new,mcmcpars,confidence_limit,model, pars_without_nan, s2chain)
%FIND_CONFIDENCE_INTERVALS : by Raunak Dey: gives confidence limits, max + min; inputs:
%chain, %of intervals, 
%   Detailed explanation goes here

chain_to_be_used = chain(transient_id_new:end,:);
[length_chain_to_be_used,~]=size(chain_to_be_used);



tvec = 0:0.05:15.75; % for better viz


for i = 1:length_chain_to_be_used
    i
    %pars_from_dist = @(chain_to_be_used) chain_to_be_used;
    pars_used = update_pars(pars_without_nan,chain_to_be_used(i,:),mcmcpars);
    [~,S_recon(:,:,i),V_recon(:,:,i),D_recon(:,:,i)] = simulate_ode(model,pars_used,tvec,pars_used.S0,pars_used.V0); % mcmc parameter set
    
    %S_recon(:,:,i) = exp(log(S_recon(:,:,i)) + randn(size(S_recon(:,:,i))) * sqrt(s2chain(i))  );
    %V_recon(:,:,i) = exp(log(V_recon(:,:,i)) + randn(size(V_recon(:,:,i))) * sqrt(s2chain(i))  );
    
    S_recon(:,:,i) = S_recon(:,:,i) + randn(size(S_recon(:,:,i))) * sqrt(s2chain(i))  ;
    V_recon(:,:,i) = V_recon(:,:,i) + randn(size(V_recon(:,:,i))) * sqrt(s2chain(i))  ;
    

end



for i = 1:model.NH
    S_reconstructed_i = reshape(S_recon(:,i,:),length(tvec),length_chain_to_be_used);
    S_min(:,i) = quantile(S_reconstructed_i',(1-confidence_limit)/2);
    S_max(:,i) = quantile(S_reconstructed_i',1-(1-confidence_limit)/2);
    S_median(:,i) = median(S_reconstructed_i');
    clear S_reconstructed_i;
end



for i=1:model.NV

    V_reconstructed_i = reshape(V_recon(:,i,:),length(tvec),length_chain_to_be_used);
    V_min(:,i) = quantile(V_reconstructed_i',(1-confidence_limit)/2);
    V_max(:,i) = quantile(V_reconstructed_i',1-(1-confidence_limit)/2);
    V_median(:,i) = median(V_reconstructed_i');
    clear V_reconstructed_i;

end


end





