clear all;
clc;

%% SEIVD
x=1:9;



load('./v14-4.mat');

% change the log transformed part of the chain to original
chain(:,10:18) = exp(2.303.*chain(:,10:18));

%start_new = 7000;
start_new = 30001;
beta(:,1) = median(chain(start_new:end,1:9)) ;
r(:,1) = median(chain(start_new:end,28:32)) ;
phi(:,1) = median(chain(start_new:end,10:18)) ;
tau(:,1) = median(chain(start_new:end,19:27)) ;

error_beta_seivd = std(chain(start_new:end,1:9)) ;
error_r_seivd = std(chain(start_new:end,28:32)) ;
error_phi_seivd = std(chain(start_new:end,10:18)) ;
error_tau_seivd = std(chain(start_new:end,19:27)) ;


%% one-steps

% one steps conventional 
beta_onestep_conventional = [0.94,91.2,25.7,0,10.4,44.3,319,54.2,357];
tau_onestep_conventional = [1,1.25,1.5,0,1,0.833,1.5,0.667,1.167];
phi_onestep_conventional = [1.83e-7,1.88e-8,1.33e-7,0,9.9e-8,1.79e-7,7.6e-8,1.87e-7,1e-7];

%one step bayesian (these are corrected for the indices)
beta_onestep = [118.2,1.18,60.5,-200,8.246,99.2,93,437,413];
tau_onestep = [1.89,1.749,2.19,0,1.9,1.47,1.42,2.14,1.97];
phi_onestep = [ 1.45e-8,5.14e-8, 9.57e-8,0,1.227e-7,1.56e-7, 1.31e-7, 6.46e-8, 8.02e-8 ];

beta_onestep_error = [68.8,2.02,50,0,21.9,43.8,50,79.5,86];
tau_onestep_error = [0.18,0.56,0.19,0,0.4,0.13,0.12,0.13,0.09];
phi_onestep_error = [9.57e-9, 1.92e-8,1.39e-8,0,1.68e-8, 2.29e-8,1.99e-8,1.75e-8,1.7e-8];

r_onestep_conventional = [0.19,0.245,0.22,0.28,0.25];



% these are the old ones with wrong indices serial numbers
% beta_onestep = [1.81,118.2,60.5,-200,8.246,99.2,437,93,413];
% tau_onestep = [1.749,1.89,2.19,0,1.9,1.47,2.14,1.42,1.97];
% phi_onestep = [5.14e-8,1.45e-8,9.57e-8,0,1.227e-7,1.56e-7, 6.46e-8,1.31e-7,8.02e-8 ];
% 
% beta_onestep_error = [2.02,68.8,50,0,21.9,43.8,79.5,50,86];
% tau_onestep_error = [0.56,0.18,0.19,0,0.4,0.13,0.13,0.12,0.09] ;
% phi_onestep_error = [1.92e-8,9.57e-9,1.39e-8,0,1.68e-8, 2.29e-8,1.75e-8,1.99e-8,1.7e-8];
% 
% r_onestep_conventional = [0.19,0.245,0.22,0.28,0.25];








boxes_onestep = [77,74,85,0,55,83,98,87,109];
boxes_onestep_error = [41,39,37,0,40,34,35,35,31];



%%  one step vs community

gap = 0.1;
market_size_given = 10;
%color_green = [70/255,210/255,130/255];
color_green = [171,193,157]./255;

figure
subplot(2,2,1)
errorbar(x-gap,beta,error_beta_seivd,'ro','MarkerSize',market_size_given,'MarkerFaceColor','auto','Marker','diamond','LineWidth',2);
hold on;
errorbar(x+gap,beta_onestep,beta_onestep_error,'o','MarkerSize',market_size_given,'MarkerFaceColor',color_green,'Color',color_green,'LineWidth',2);
set(gca,'FontSize',20);
xticks(1:9);
ylabel({'Burst sizes' ;'\beta (virions/cell)'});
ylim([0 600]);
yticks(0:100:600);
xline(5.5,'--k',LineWidth=2.5);
%legend('Inferred parameters from scaled-up SEIVD model (Fig 4)','Inferred parameters from scaled-up SEIV model (Fig 3)','Inferred parameters from pairwise SEIV model (Fig 2)');
 %legend('Parameters inferred with scaled-up SEIVD community model','Parameters inferred with scaled-up SEIV community model','Parameters from pairwise SEIV model');
 %legend('Box','off');

subplot(2,2,2)
errorbar(x-gap,tau,error_tau_seivd,'ro','MarkerSize',market_size_given,'MarkerFaceColor','auto','Marker','diamond','LineWidth',2);
hold on;
errorbar(x+gap,tau_onestep,tau_onestep_error,'o','MarkerSize',market_size_given,'MarkerFaceColor',color_green,'Color',color_green,'LineWidth',2);
set(gca,'FontSize',20);
xticks(1:9);
ylabel({'Latent periods'; ' \tau (hr)'});
ylim([0.2 6]);
xline(5.5,'--k',LineWidth=2.5);




% subplot(2,2,3)
% errorbar(x,exp(2.3026*phi(:,1)),exp(2.3026*error_phi_seivd),'ko','MarkerSize',8);
% hold on;
% errorbar(x,exp(2.3026*phi(:,2)),exp(2.3026*error_phi_seiv),'ro','MarkerSize',8);
% errorbar(x,phi_onestep,phi_onestep_error,'bo','MarkerSize',8);
% set(gca,'FontSize',20);
% xticks(1:9);
% ylabel('\phi');

%error_phi_seivd_rescaled = (exp(2.3026*(phi + error_phi_seivd')) - exp(2.3026*(phi - error_phi_seivd')))/2;

% 
subplot(2,2,3)
errorbar(x-gap,phi,error_phi_seivd,'ro','MarkerSize',market_size_given,'MarkerFaceColor','auto','Marker','diamond','LineWidth',2);
hold on;
errorbar(x+gap,phi_onestep, phi_onestep_error,'o','MarkerSize',market_size_given,'MarkerFaceColor',color_green,'Color',color_green,'LineWidth',2);
set(gca,'FontSize',20);
xticks(1:9);
set(gca, 'YScale', 'log');
ylabel('\phi (ml/hr)');
%yticks([1e-10 1e-9 1e-8 1e-7]);
xline(5.5,'--k',LineWidth=2.5);
%ylim([1e-10 10e-7]);


% subplot(2,2,3)
% errorbar(x-gap,phi(:,1),error_phi_seivd,'ro','MarkerSize',market_size_given,'MarkerFaceColor','auto','Marker','diamond');
% hold on;
% errorbar(x,phi(:,2),error_phi_seiv,'ko','MarkerSize',market_size_given,'MarkerFaceColor','k','Marker','square');
% plot(x+gap,log(phi_onestep)/2.3026,'o','MarkerSize',market_size_given,'MarkerFaceColor',color_green,'Color',color_green);
% set(gca,'FontSize',20);
% xticks(1:9);
% ylabel({'Log-adsorption rates','log_{10}\phi (ml/hr)'});
% ylim([-10 -6.5]);
% xline(5.5,'--k',LineWidth=2.5);


subplot(2,2,4)
set(gca, 'YScale', 'linear')
errorbar((1:5)-gap ,r,error_r_seivd,'ro','MarkerSize',market_size_given,'MarkerFaceColor','auto','Marker','diamond','LineWidth',2);
hold on;
plot((1:5)+gap,r_onestep_conventional,'o','MarkerSize',market_size_given,'MarkerFaceColor',color_green,'Color',color_green,'LineWidth',2);
set(gca,'FontSize',20);
xticks(1:5);
xticklabels({'CBA 4','CBA 18','CBA 38','PSA H100','PSA 13-15'});
set(gca,'XTickLabelRotation',90)
ylabel({'Growth rates','r (cells/hr)'});
xline(3.5,'--k',LineWidth=2.5);

