%load('./../corrected3/revised80002.mat');
%close all;

%load('./../corrected3/revised90015.mat');

%load('./../corrected4/from_os_prior_20000.mat');
%load('./../corrected4/v14-4.mat');
load('./../corrected4/v17.mat');


close all;



transparency = 0.1;
color = [70/255,130/255,180/255];

figure(1)
subplot(2,5,6)
plot(data.xdata,data.ydata(:,6),'o','MarkerEdgeColor','k','MarkerFaceColor',color);hold on;
set(gca, 'YScale', 'log');
subplot(2,5,7)
plot(data.xdata,data.ydata(:,7),'o','MarkerEdgeColor','k','MarkerFaceColor',color);hold on;
set(gca, 'YScale', 'log');
subplot(2,5,8)
plot(data.xdata,data.ydata(:,8),'o','MarkerEdgeColor','k','MarkerFaceColor',color);hold on;
set(gca, 'YScale', 'log');
subplot(2,5,9)
plot(data.xdata,data.ydata(:,9),'o','MarkerEdgeColor','k','MarkerFaceColor',color);hold on;
set(gca, 'YScale', 'log');
subplot(2,5,10)
plot(data.xdata,data.ydata(:,10),'o','MarkerEdgeColor','k','MarkerFaceColor',color);hold on;
set(gca, 'YScale', 'log');


subplot(2,5,1)
plot(data.xdata,data.ydata(:,1),'o','MarkerEdgeColor','k','MarkerFaceColor',color);hold on;
set(gca, 'YScale', 'log');
subplot(2,5,2)
plot(data.xdata,data.ydata(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor',color);hold on;
set(gca, 'YScale', 'log');
subplot(2,5,3)
plot(data.xdata,data.ydata(:,3),'o','MarkerEdgeColor','k','MarkerFaceColor',color);hold on;
set(gca, 'YScale', 'log');
subplot(2,5,4)
plot(data.xdata,data.ydata(:,4),'o','MarkerEdgeColor','k','MarkerFaceColor',color);hold on;
set(gca, 'YScale', 'log');
subplot(2,5,5)
plot(data.xdata,data.ydata(:,5),'o','MarkerEdgeColor','k','MarkerFaceColor',color);hold on;
set(gca, 'YScale', 'log');



    
  %%

for i = 100:100:5000
 pars_samples = update_pars(pars2,chain(i,:),mcmcpars);
 pars_samples.beta2 = pars_samples.beta;

[t3,S3,V3,~] = simulate_ode(model,pars_samples,tvec,pars2.S0,pars2.V0);
    figure(1)
    subplot(2,5,6)
    patchline(t3,V3(:,1),'edgecolor',[0 0 0],'linewidth',1,'edgealpha',transparency);hold on;
  
    subplot(2,5,7)
    patchline(t3,V3(:,2),'edgecolor',[0 0 0],'linewidth',1,'edgealpha',transparency);hold on;

    subplot(2,5,8)
    patchline(t3,V3(:,3),'edgecolor',[0 0 0],'linewidth',1,'edgealpha',transparency);hold on;

    subplot(2,5,9)
    patchline(t3,V3(:,4),'edgecolor',[0 0 0],'linewidth',1,'edgealpha',transparency);hold on;

    subplot(2,5,10)
    patchline(t3,V3(:,5),'edgecolor',[0 0 0],'linewidth',1,'edgealpha',transparency);hold on;

    
    subplot(2,5,1)
    patchline(t3,S3(:,1),'edgecolor',[0 0 0],'linewidth',1,'edgealpha',transparency);hold on;

    subplot(2,5,2)
    patchline(t3,S3(:,2),'edgecolor',[0 0 0],'linewidth',1,'edgealpha',transparency);hold on;
    
    subplot(2,5,3)
    patchline(t3,S3(:,3),'edgecolor',[0 0 0],'linewidth',1,'edgealpha',transparency);hold on;
    
    subplot(2,5,4)
    patchline(t3,S3(:,4),'edgecolor',[0 0 0],'linewidth',1,'edgealpha',transparency);hold on;
    
    subplot(2,5,5)
    patchline(t3,S3(:,5),'edgecolor',[0 0 0],'linewidth',1,'edgealpha',transparency);hold on;


end

%% touch up

figure(1)

subplot(2,5,1)
hold on;
set(gca,'FontSize',14);
title('CBA 4');
ylabel('Host density (cells/ml)')
xticks(0:2:16);
ylim([1e5 1e8]);
yticks([1e5, 1e6, 1e7, 1e8]);

subplot(2,5,2)
hold on;
set(gca,'FontSize',14);
title('CBA 18');
xticks(0:2:16);
ylim([1e5 1e8]);
yticks([1e5, 1e6, 1e7, 1e8]);


subplot(2,5,3)
hold on;
set(gca,'FontSize',14);
title('CBA 38');
xticks(0:2:16);
ylim([1e5 1e8]);
yticks([1e5, 1e6, 1e7, 1e8]);

subplot(2,5,4)
hold on;
set(gca,'FontSize',14);
title('PSA H100');
xticks(0:2:16);
ylim([1e5 1e8]);
yticks([1e5, 1e6, 1e7, 1e8]);

subplot(2,5,5)
hold on;
set(gca,'FontSize',14);
title('PSA 13-15');
xticks(0:2:16);
ylim([1e5 1e8]);
yticks([1e5, 1e6, 1e7, 1e8]);


subplot(2,5,6)
hold on;
set(gca,'FontSize',14);
title('\phi18:2');
ylabel('Phage density (virions/ml)')
xticks(0:2:16);
ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]);




subplot(2,5,7)
hold on;
set(gca,'FontSize',14);
title('\phi18:3');
xticks(0:2:16);
ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]);



subplot(2,5,8)
hold on;
set(gca,'FontSize',14);
title('\phi38:1');
xticks(0:2:16);
ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]);
xlabel('Time (hrs)')


subplot(2,5,9)
hold on;
set(gca,'FontSize',14);
title('PSA HP1');
xticks(0:2:16);
ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]);


subplot(2,5,10)
hold on;
set(gca,'FontSize',14);
title('PSA HS6');
xticks(0:2:16);
ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]);

    