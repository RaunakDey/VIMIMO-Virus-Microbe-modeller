clc;
clear all;

%% load data
load('./v14-4.mat');

%% resample
cutoff = 250;
chain_resampled = chain(1,:);

for i = 30000:100:40000
    i
    if loglikefun(chain(i,:),data,pars2,mcmcpars,model,0)< cutoff
        chain_resampled(end+1,:) = chain(i,:);
    end
end
