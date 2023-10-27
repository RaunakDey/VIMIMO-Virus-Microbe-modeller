clc;
clear all;

confidence_run = false;

% colors
blue = [0 0.4470 0.7410];
yellow = [ 0.9290 0.6940 0.1250];
purple = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];

addpath(genpath('./..'));



seed = 20031;
load('./../results/SEIVD-diff-all-seed'+string(seed)+'.mat')
chain_stored = chain(1:10000,:);

seed2 = 20029;
load('./../results/SEIVD-diff-all-seed'+string(seed2)+'.mat')
chain_stored2 = chain;

seed3 = 20027;
load('./../results/SEIVD-diff-all-seed'+string(seed3)+'.mat')
chain_stored3 = chain;


seed4 = 20030;
load('./../results/SEIVD-diff-all-seed'+string(seed4)+'.mat')
chain_stored4 = chain;

%% add the beta part
%14,17,18
% 31 29 30 27

% 32, 30, 31, 35

load('./../results/SEIVD-onlybeta-seed20032.mat','chain')
chain_stored(:,6:14) = chain;

load('./../results/SEIVD-onlybeta-seed20030.mat','chain')
chain_stored2(:,6:14) = chain;

load('./../results/SEIVD-onlybeta-seed20031.mat','chain')
chain_stored3(:,6:14) = chain;

load('./../results/SEIVD-onlybeta-seed20030.mat','chain')
chain_stored4(:,6:14) = chain;


%% add the tau part
% do same thing as tau
%31,32,33

% load('./../results/SEIVD-beta-tau-seed20043.mat','chain');
% chain_stored(:,34:42) = chain(:,10:18);
% 
% load('./../results/SEIVD-beta-tau-seed20047.mat','chain');
% chain_stored2(:,34:42) = chain(:,10:18);
% 
% load('./../results/SEIVD-beta-tau-seed20045.mat','chain');
% chain_stored3(:,34:42) = chain(:,10:18);
% 
% load('./../results/SEIVD-beta-tau-seed20046.mat','chain');
% chain_stored4(:,34:42) = chain(:,10:18);


load('./../results/SEIVD-onlytau-seed20039.mat','chain');
chain_stored(:,34:42) = chain;

load('./../results/SEIVD-onlytau-seed20040.mat','chain');
chain_stored2(:,34:42) = chain;

load('./../results/SEIVD-onlytau-seed20041.mat','chain');
chain_stored3(:,34:42) = chain;

load('./../results/SEIVD-onlytau-seed20042.mat','chain');
chain_stored4(:,34:42) = chain;

%% add the phi part
load('./../results/SEIVD-phi-seed20031.mat','chain');
chain_stored(:,15:23) = chain;

load('./../results/SEIVD-phi-seed20032.mat','chain');
chain_stored2(:,15:23) = chain;

load('./../results/SEIVD-phi-seed20033.mat','chain');
chain_stored3(:,15:23) = chain;

load('./../results/SEIVD-phi-seed20034.mat','chain');
chain_stored4(:,15:23) = chain;


%% priors

priors.r.mean = [0.19 0.24 0.24 0.28 0.25];
priors.r.std = 0.02245*ones(1,5);

priors.beta.mean = [1.5 75 29.6 70 37.5 102.4 437 93.0 413.7];
priors.beta.std = 300*ones(1,9);

priors.tau.mean = [1.36 2.07 1.4 1 1.8 1.5 2.1 1.42 2.0];
priors.tau.std = 10*ones(1,9);

priors.phi.mean = log([5.3 1.2 14 10 4.97 16 6.5 13 8]*1e-8)/log(10);
priors.phi.std = 4.32*ones(1,9);

title_names = ['\phi18:2--CBA4','\phi18:3--CBA18','\phi18:2--CBA18','\phi38:1--CBA18','\phi38:1--CBA38','PSA HP1--PSA H100','PSA HS6--PSA H100','PSA HP1--PSA 13-15','PSA HS6--PSA 13-15'];
%% plots

figure(1)
for i = 1:9
subplot(9,2,2*i-1)
plot(chain_stored(5000:10000,i+5),Color=blue);hold on;
xlim([0 5000]);
if i == 9
    xlabel('trace positions')
end
%title(title_names(i))

ylabel('\beta')



subplot(9,2,2*i)
histogram(chain_stored(5001:10000,i+5),'DisplayStyle','stairs','NumBins',25,'Normalization','pdf','LineWidth',1,'EdgeColor',blue); hold on;
plot(gaussian(0:1000,priors.beta.mean(i),priors.beta.std(i)),'LineWidth',2,'Color','k');
xline(priors.beta.mean(i),Color='r',LineWidth=2);
%title(title_names(i))

xlim([0 1000])
end
xlabel('\beta ')


figure(2)
for i = 1:9
subplot(9,2,2*i-1)
plot(chain_stored(5000:10000,i+33),Color=blue); hold on;
xlim([0 5000]);

subplot(9,2,2*i)
histogram(chain_stored(5001:10000,i+33),'DisplayStyle','stairs','NumBins',100,'Normalization','pdf','LineWidth',1,'EdgeColor',blue);
hold on;
plot(-1:20,gaussian(-1:20,priors.tau.mean(i),priors.tau.std(i)),'LineWidth',2,'Color','k')
xline(priors.tau.mean(i),Color='r',LineWidth=2);

xlim([0 20])

end
xlabel('\tau (hrs)')

figure(3)
for i = 1:9
subplot(9,2,2*i-1)
plot(chain_stored(5000:10000,i+14),Color=blue); hold on;
xlim([0 5000]);

subplot(9,2,2*i)
histogram(chain_stored(5001:10000,i+14),'DisplayStyle','stairs','NumBins',50,'Normalization','pdf','LineWidth',1,'EdgeColor',blue);
hold on;
plot(-15:0.1:-5,gaussian(-15:0.1:-5,priors.phi.mean(i),priors.phi.std(i)),'LineWidth',2,'Color','k');
xline(priors.phi.mean(i),Color='r',LineWidth=2);
set(gca, 'YScale', 'log');
xlim([-15 -5]);

end
xlabel('log(\phi) ml/hrs')



figure(4)
for i = 1:5
subplot(5,2,2*i-1)
plot(chain_stored(5000:10000,i),Color=blue); hold on;
xlim([0 5000]);

subplot(5,2,2*i)
histogram(chain_stored(5001:10000,i),'DisplayStyle','stairs','NumBins',50,'Normalization','pdf','LineWidth',1,'EdgeColor',blue);
hold on;

plot(0:0.01:0.8,gaussian(0:0.01:0.8,priors.r.mean(i),priors.r.std(i)),'LineWidth',2,'Color','k')
xline(priors.r.mean(i),Color='r',LineWidth=2);
xlim([0 0.8]);

end
xlabel('r (hrs^{-1})')


%% second chain

figure(1)
for i = 1:9
subplot(9,2,2*i-1)
plot(chain_stored2(5000:10000,i+5),Color=yellow); hold on;
xlim([0 5000]);

subplot(9,2,2*i)
histogram(chain_stored2(5001:10000,i+5),'DisplayStyle','stairs','NumBins',25,'Normalization','pdf','LineWidth',1,'EdgeColor',yellow); hold on;

xlim([0 1000])
end
xlabel('\beta ')


figure(2)
for i = 1:9
subplot(9,2,2*i-1)
plot(chain_stored2(5000:10000,i+33),Color=yellow); hold on;
xlim([0 5000]);

subplot(9,2,2*i)
histogram(chain_stored2(5001:10000,i+33),'DisplayStyle','stairs','NumBins',50,'Normalization','pdf','LineWidth',1,'EdgeColor',yellow);

xlim([0 20])

end
xlabel('\tau (hrs)')

figure(3)
for i = 1:9
subplot(9,2,2*i-1)
plot(chain_stored2(5000:10000,i+14),Color=yellow); hold on;
xlim([0 5000]);

subplot(9,2,2*i)
histogram(chain_stored2(5001:10000,i+14),'DisplayStyle','stairs','NumBins',50,'Normalization','pdf','LineWidth',1,'EdgeColor',yellow);
hold on;
set(gca, 'YScale', 'log');
xlim([-15 -5]);

end
xlabel('log(\phi) ml/hrs')



figure(4)
for i = 1:5
subplot(5,2,2*i-1)
plot(chain_stored2(5000:10000,i),Color=yellow); hold on;
xlim([0 5000]);

subplot(5,2,2*i)
histogram(chain_stored2(5001:10000,i),'DisplayStyle','stairs','NumBins',50,'Normalization','pdf','LineWidth',1,'EdgeColor',yellow);
hold on;
xlim([0 0.8]);

end
xlabel('r (hrs^{-1})')


%% third chain

figure(1)
for i = 1:9
subplot(9,2,2*i-1)
plot(chain_stored3(5000:10000,i+5),Color=purple); hold on;
xlim([0 5000]);

subplot(9,2,2*i)
histogram(chain_stored3(5001:10000,i+5),'DisplayStyle','stairs','NumBins',25,'Normalization','pdf','LineWidth',1,'EdgeColor',purple); hold on;

xlim([0 1000])
end
xlabel('\beta ')


figure(2)
for i = 1:9
subplot(9,2,2*i-1)
plot(chain_stored3(5000:10000,i+33),Color=purple); hold on;
xlim([0 5000]);

subplot(9,2,2*i)
histogram(chain_stored3(5001:10000,i+33),'DisplayStyle','stairs','NumBins',50,'Normalization','pdf','LineWidth',1,'EdgeColor',purple);

xlim([0 20])

end
xlabel('\tau (hrs)')

figure(3)
for i = 1:9
subplot(9,2,2*i-1)
plot(chain_stored3(5000:10000,i+14),Color=purple); hold on;
xlim([0 5000]);

subplot(9,2,2*i)
histogram(chain_stored3(5001:10000,i+14),'DisplayStyle','stairs','NumBins',50,'Normalization','pdf','LineWidth',1,'EdgeColor',purple);
hold on;
set(gca, 'YScale', 'log');
xlim([-10 -7]);

end
xlabel('log(\phi) ml/hrs')



figure(4)
for i = 1:5
subplot(5,2,2*i-1)
plot(chain_stored3(5000:10000,i),Color=purple); hold on;
xlim([0 5000]);

subplot(5,2,2*i)
histogram(chain_stored3(5001:10000,i),'DisplayStyle','stairs','NumBins',50,'Normalization','pdf','LineWidth',1,'EdgeColor',purple);
hold on;
xlim([0 0.8]);

end
xlabel('r (hrs^{-1})')

%% fourth

figure(1)
for i = 1:9
subplot(9,2,2*i-1)
plot(chain_stored4(5000:10000,i+5),Color=green); hold on;
xlim([0 5000]);

subplot(9,2,2*i)
histogram(chain_stored4(5001:10000,i+5),'DisplayStyle','stairs','NumBins',25,'Normalization','pdf','LineWidth',1,'EdgeColor',green); hold on;

xlim([0 1000])
end
xlabel('\beta ')


figure(2)
for i = 1:9
subplot(9,2,2*i-1)
plot(chain_stored4(5000:10000,i+33),Color=green); hold on;
xlim([0 5000]);

subplot(9,2,2*i)
histogram(chain_stored4(5001:10000,i+33),'DisplayStyle','stairs','NumBins',50,'Normalization','pdf','LineWidth',1,'EdgeColor',green);

xlim([0 20])

end
xlabel('\tau (hrs)')

figure(3)
for i = 1:9
subplot(9,2,2*i-1)
plot(chain_stored4(5000:10000,i+14),Color=green); hold on;
xlim([0 5000]);

subplot(9,2,2*i)
histogram(chain_stored4(5001:10000,i+14),'DisplayStyle','stairs','NumBins',50,'Normalization','pdf','LineWidth',1,'EdgeColor',green);
hold on;
set(gca, 'YScale', 'log');
xlim([-12 -6]);

end
xlabel('log(\phi) ml/hrs')



figure(4)
for i = 1:5
subplot(5,2,2*i-1)
plot(chain_stored4(5000:10000,i),Color=green); hold on;
xlim([0 5000]);

subplot(5,2,2*i)
histogram(chain_stored4(5001:10000,i),'DisplayStyle','stairs','NumBins',50,'Normalization','pdf','LineWidth',1,'EdgeColor',green);
hold on;
xlim([0 0.8]);

end
xlabel('r (hrs^{-1})')


%%
if confidence_run == true

pars_from_dist = @(chain) median(chain);
pars2 = update_pars(pars2,pars_from_dist(chain_stored4(5000:end,:)),mcmcpars);
tvec = 0:0.05:15.75; % for better viz
[t2,S2,V2,D2] = simulate_ode(model,pars2,tvec,pars2.S0,pars2.V0); % mcmc parameter set


[S_min,S_max,V_min,V_max,S_median,V_median] = find_confidence_interval_looped(chain_stored4,5000,mcmcpars,0.95,model, pars2);


load('./../data/triplicate_data.mat');

 time_2 = [t2', fliplr(t2')];
linewidth = 2;

color_ofthe_fit = [1 0 0]*0.5;
color_ofthe_fill = [0.95 0 0];
transparency = 0.25;


hf4 = figure;
subplot(2,5,1)
errorbar(time/60,mean(1e3*host1'),std(1e3*host1'),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times
ylim([1e6 1e7]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
ylabel({'Host density';'(cells/ml)'});
title('CBA 4','FontSize',18);
inBetween = [S_min(:,1)', fliplr(S_max(:,1)')];
    fill(time_2, inBetween, color_ofthe_fill,'FaceAlpha',transparency ,'LineStyle','none');
    plot(t2,S_median(:,1),'--','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,2)
errorbar(time/60,mean(1e3*host2'),std(1e3*host2'),'o','MarkerSize',8,  'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255] );hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 18','FontSize',18);
    inBetween = [S_min(:,2)', fliplr(S_max(:,2)')];
    fill(time_2, inBetween, color_ofthe_fill,'FaceAlpha',transparency ,'LineStyle','none');
    plot(t2,S_median(:,2),'--','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,3)
errorbar(time/60,mean(1e3*host3'),std(1e3*host3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times
ylim([5e5 2e7]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
title('CBA 38','FontSize',18);
inBetween = [S_min(:,3)', fliplr(S_max(:,3)')];
    fill(time_2, inBetween, color_ofthe_fill,'FaceAlpha',transparency ,'LineStyle','none');
    plot(t2,S_median(:,3),'--','Color',color_ofthe_fit,'LineWidth',linewidth);



subplot(2,5,4)
errorbar(time/60,mean(1e3*host4'),std(1e3*host4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA H100','FontSize',18);
    inBetween = [S_min(:,4)', fliplr(S_max(:,4)')];
    fill(time_2, inBetween, color_ofthe_fill,'FaceAlpha',transparency ,'LineStyle','none');
    plot(t2,S_median(:,4),'--','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,5)
errorbar(time/60,mean(1e3*host5'),std(1e3*host5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e6 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA 13-15','FontSize',18);
    inBetween = [S_min(:,5)', fliplr(S_max(:,5)')];
    fill(time_2, inBetween, color_ofthe_fill,'FaceAlpha',transparency ,'LineStyle','none');
    plot(t2,S_median(:,5),'--','Color',color_ofthe_fit,'LineWidth',linewidth);

%xlabel("Time (hours)");
%ylabel("Host density (cell/ml)");






subplot(2,5,6)
errorbar(time/60,mean(1e3*virus1'),std(1e3*virus1'),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e10]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e4 1e6 1e8 1e10]);
ylabel({'Phage density';'(virions/ml)'});
title('\phi18:2','FontSize',18);
    inBetween = [V_min(:,1)', fliplr(V_max(:,1)')];
    fill(time_2, inBetween, color_ofthe_fill,'FaceAlpha',transparency ,'LineStyle','none');
    plot(t2,V_median(:,1),'--','Color',color_ofthe_fit,'LineWidth',linewidth);

subplot(2,5,7)
errorbar(time/60,mean(1e3*virus2'),std(1e3*virus2'),'o','MarkerSize',8,  'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255] );hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e10]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e4 1e6 1e8 1e10]);
    title('\phi18:3','FontSize',18);
    inBetween = [V_min(:,2)', fliplr(V_max(:,2)')];
    fill(time_2, inBetween, color_ofthe_fill,'FaceAlpha',transparency ,'LineStyle','none');
    plot(t2,V_median(:,2),'--','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,8)
errorbar(time/60,mean(1e3*virus3'),std(1e3*virus3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e9]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('\phi38:1','FontSize',18);
    inBetween = [V_min(:,3)', fliplr(V_max(:,3)')];
    fill(time_2, inBetween, color_ofthe_fill,'FaceAlpha',transparency ,'LineStyle','none');
    plot(t2,V_median(:,3),'--','Color',color_ofthe_fit,'LineWidth',linewidth);

subplot(2,5,9)
errorbar(time/60,mean(1e3*virus4'),std(1e3*virus4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('PSA HP1','FontSize',18);
    inBetween = [V_min(:,4)', fliplr(V_max(:,4)')];
    fill(time_2, inBetween, color_ofthe_fill,'FaceAlpha',transparency ,'LineStyle','none');
    plot(t2,V_median(:,4),'--','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,10)
errorbar(time/60,mean(1e3*virus5'),std(1e3*virus5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e6 1e10]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('PSA HS6','FontSize',18);
    inBetween = [V_min(:,5)', fliplr(V_max(:,5)')];
    fill(time_2, inBetween, color_ofthe_fill,'FaceAlpha',transparency ,'LineStyle','none');
    plot(t2,V_median(:,5),'--','Color',color_ofthe_fit,'LineWidth',linewidth);
    %legend('Data','95% confidence interval','Bayesian fit');
    %legend('Box','off');

han=axes(hf4,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
xlabel("Time (hours)");

end
