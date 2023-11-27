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
transient_id = 300;
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
%mcmcoptions.adaptint = 1000; %100 is default

% proposal kernel
% mcmcoptions.qcov = zeros(37,37);
% for i = 1:9
% mcmcoptions.qcov(i,i) = 1;
% end
% 
% for i = 10:18
% mcmcoptions.qcov(i,i) = 0.01;
% end
% 
% for i = 19:27
% mcmcoptions.qcov(i,i) = 0.01;
% end
% 
% for i = 28:32
% mcmcoptions.qcov(i,i) = 0.01;
% end
% 
% for i = 33:37
% mcmcoptions.qcov(i,i) = 0.01;
% end



mcmcoptions.method  = 'dram';
mcmcpars = mcmcpars_setup_new(pars2,pars2,include_pars,flags,model);
mcmcparam = mcmcpars2param(mcmcpars);

%% priors and initial starting positions
% priors are truncated gaussians
% name -- start -- min -- max -- mean -- std


load('./v14-4.mat','theta_optimized','mcmcresults');

%%% initial values set
for i = 1:37
mcmcparam{1,i}{1,2} = theta_optimized(i);
end

for i = 1:9
mcmcparam{1,i}{1,6} = 5000;
mcmcresults.qcov(i,i) = 10;
mcmcresults.qcov2(i,i) = 10;
mcmcresults.RDR{1,1}(i,i) = 10*mcmcresults.RDR{1,1}(i,i);
end


mcmcmodel.S20 = 50;
mcmcmodel.sigma2 = 50;
mcmcmodel.N0 = 10;
mcmcresults.S20 = 50;
mcmcresults.N0 = 10;

%%% setting up means of priors.




%first save
save("v14-6.mat")

%% mcmc run

mcmcmodel.ssfun = @(theta,data) ssfun(theta,data,pars2,mcmcpars,model,lambda); 
[mcmcresults, chain, s2chain]= mcmcrun(mcmcmodel,data,mcmcparam,mcmcoptions,mcmcresults);


pars_afterinf = update_pars(pars2,median(chain(transient_id:end, :)),mcmcpars);
[t_after,S_after,V_after,D_after,I_after,E_after] =  simulate_ode(model,pars_afterinf,tvec,pars2.S0,pars2.V0); % mcmc parameter set

%second save
save("v14-6.mat")
%% new plots
% hf4 = figure;
% subplot(2,5,1)
% errorbar(time/60,mean(1e3*host1'),std(1e3*host1'),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
% set(gca, 'YScale', 'log');
% set(gca,'fontname','times')  % Set it to times
% ylim([1e5 1e8]);
%     xlim([0 16]);
%     xticks([0 2 4 6 8 10 12 14 16]);
%     set(gca,'FontSize',20);
%     axis('square');
%     yticks([1e5 1e6 1e7 1e8]);
% ylabel({'Host density';'(cells/ml)'});
% title('CBA 4','FontSize',18);
%     plot(t_after,S_after(:,1),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
% 
% 
% 
% 
% subplot(2,5,2)
% errorbar(time/60,mean(1e3*host2'),std(1e3*host2'),'o','MarkerSize',8,  'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255] );hold on;
% set(gca, 'YScale', 'log');
% set(gca,'fontname','times')  % Set it to times
% ylim([1e5 1e8]);
%     xlim([0 16]);
%    xticks([0 2 4 6 8 10 12 14 16]);
%     set(gca,'FontSize',20);
%     axis('square');
%     yticks([1e5 1e6 1e7 1e8]);
%     title('CBA 18','FontSize',18);
%    
%     plot(t_after,S_after(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
% 
% 
% 
% 
% subplot(2,5,3)
% errorbar(time/60,mean(1e3*host3'),std(1e3*host3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
% set(gca, 'YScale', 'log');
% set(gca,'fontname','times')  % Set it to times
% ylim([1e5 1e8]);
%     xlim([0 16]);
%     xticks([0 2 4 6 8 10 12 14 16]);
%     set(gca,'FontSize',20);
%     axis('square');
%     yticks([1e5 1e6 1e7 1e8]);
% title('CBA 38','FontSize',18);
% 
%     plot(t_after,S_after(:,3),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
% 
% 
% 
% subplot(2,5,4)
% errorbar(time/60,mean(1e3*host4'),std(1e3*host4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
% set(gca, 'YScale', 'log');set(gca,'FontSize',20)
% set(gca,'fontname','times')  % Set it to times
% ylim([1e5 1e8]);
%     xlim([0 16]);
%   xticks([0 2 4 6 8 10 12 14 16]);
%   axis('square');
%     yticks([1e5 1e6 1e7 1e8]);
%     title('PSA H100','FontSize',18);
%     
%     plot(t_after,S_after(:,4),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
% 
% 
% 
% 
% subplot(2,5,5)
% errorbar(time/60,mean(1e3*host5'),std(1e3*host5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
% set(gca, 'YScale', 'log');set(gca,'FontSize',20);
% set(gca,'fontname','times')  % Set it to times
% ylim([1e5 1e8]);
%     xlim([0 16]);
%     xticks([0 2 4 6 8 10 12 14 16]);
%     axis('square');
%     yticks([1e5 1e6 1e7 1e8]);
%     title('PSA 13-15','FontSize',18);
%     
%     plot(t_after,S_after(:,5),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
% 
% %xlabel("Time (hours)");
% %ylabel("Host density (cell/ml)");
% 
% 
% 
% 
% 
% 
% subplot(2,5,6)
% errorbar(time/60,mean(1e3*virus1'),std(1e3*virus1'),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
% set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
% ylim([1e4 1e11]);
%     xlim([0 16]);
%    xticks([0 2 4 6 8 10 12 14 16]);
%     set(gca,'FontSize',20);
%     axis('square');
%     yticks([1e4 1e6 1e8 1e10]);
% ylabel({'Phage density';'(virions/ml)'});
% title('\phi18:2','FontSize',18);
%     
%     plot(t_after,V_after(:,1),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
% 
% subplot(2,5,7)
% errorbar(time/60,mean(1e3*virus2'),std(1e3*virus2'),'o','MarkerSize',8,  'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255] );hold on;
% set(gca, 'YScale', 'log');
% set(gca,'fontname','times')  % Set it to times
% ylim([1e4 1e11]);
%     xlim([0 16]);
%   xticks([0 2 4 6 8 10 12 14 16]);
%     set(gca,'FontSize',20);
%     axis('square');
%     yticks([1e4 1e6 1e8 1e10]);
%     title('\phi18:3','FontSize',18);
%     
%     plot(t_after,V_after(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
% 
% 
% 
% 
% subplot(2,5,8)
% errorbar(time/60,mean(1e3*virus3'),std(1e3*virus3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
% set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
% ylim([1e4 1e11]);
%     xlim([0 16]);
%   xticks([0 2 4 6 8 10 12 14 16]);
%     set(gca,'FontSize',20);
%     axis('square');
%    yticks([1e4 1e6 1e8 1e10]);
%    title('\phi38:1','FontSize',18);
%     
%     plot(t_after,V_after(:,3),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
% 
% subplot(2,5,9)
% errorbar(time/60,mean(1e3*virus4'),std(1e3*virus4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
% set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
% ylim([1e4 1e11]);
%     xlim([0 16]);
%    xticks([0 2 4 6 8 10 12 14 16]);
%     set(gca,'FontSize',20);
%     axis('square');
%    yticks([1e4 1e6 1e8 1e10]);
%    title('PSA HP1','FontSize',18);
%     
%     plot(t_after,V_after(:,4),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
% 
% 
% 
% 
% subplot(2,5,10)
% errorbar(time/60,mean(1e3*virus5'),std(1e3*virus5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
% set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
% ylim([1e4 1e11]);
%     xlim([0 16]);
%  xticks([0 2 4 6 8 10 12 14 16]);
%     set(gca,'FontSize',20);
%     axis('square');
%   yticks([1e4 1e6 1e8 1e10]);
%   title('PSA HS6','FontSize',18);
%     
%     plot(t_after,V_after(:,5),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
%     %legend('Data','95% confidence interval','Bayesian fit');
%     %legend('Box','off');
% 
% han=axes(hf4,'visible','off'); 
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% set(gca,'FontSize',20);
% set(gca,'fontname','times')  % Set it to times
% xlabel("Time (hours)");
% 
