clear all;
clc;

load("synthetic_Data6.mat");

% sample data




data.xdata = tvec(1:10:end);
data.ydata = S1(1:10:end,:);
data.ydata(:,NH+1:NV+NH) = V1(1:10:end,:);


data.ydata = data.ydata + randn(size(data.ydata)).*data.ydata*0.1;


data.id.S = 1:NH;
data.id.V = NH+1:NH+NV;
data.labels.phage(1:NV) = '.';
data.labels.phage_short(1:NV) = '.';
data.labels.host(1:NH) = '.';
data.labels.host_short(1:NH) = '.';



figure(1)
for i=1:NH
    subplot(2,NH,i)
plot(data.xdata,data.ydata(:,i),'Marker','o',MarkerFaceColor=[0.3 0.5 0.6],LineStyle='none');
set(gca,'YScale','log');
hold on;



subplot(2,NH,i+NH)
plot(data.xdata,data.ydata(:,i+NH),'Marker','o',MarkerFaceColor=[0.3 0.5 0.6],LineStyle='none');
set(gca,'YScale','log');
hold on;
end

%%
mcmcoptions.nsimu = 1000;
transient_id = 300;
lambda = 0;
include_pars = {'beta','phi','tau','r'};

 ssfun = @loglikefun;
% can vectorize for all the 10 time series,
% then ssfun needs to be a 10 element vector -- now set to sum.
mcmcmodel.S20 = 2.6;
mcmcmodel.N0 = 100;
mcmcoptions.updatesigma = 1;

mcmcoptions.method  = 'dram';
mcmcpars = mcmcpars_setup_silico(pars1,pars1,include_pars,flags,model);
mcmcparam = mcmcpars2param(mcmcpars);
mcmcmodel.ssfun = @(theta,data) ssfun(theta,data,pars1,mcmcpars,model,lambda); 
[mcmcresults, chain, s2chain]= mcmcrun(mcmcmodel,data,mcmcparam,mcmcoptions);


%% 

save('result-7a.mat');

%% post processing

% for i = 300:1000
%     err(i-299) = loglikefun(chain(i,:),data,pars2,mcmcpars,model,0);
% end
% 
% figure(3)
% for i=1:NH
%     subplot(2,NH,i)
% plot(data.xdata,data.ydata(:,i),'Marker','o',MarkerFaceColor=[0.3 0.5 0.9],MarkerEdgeColor ='k',LineStyle='none');
% set(gca,'YScale','log');
% hold on;
% 
% 
% 
% subplot(2,NH,i+NH)
% plot(data.xdata,data.ydata(:,i+NH),'Marker','o',MarkerFaceColor=[0.3 0.5 0.9],MarkerEdgeColor ='k', LineStyle='none');
% set(gca,'YScale','log');
% hold on;
% end
% 
for i = 300:100:10000
    i
    err = loglikefun(chain(i,:),data,pars1,mcmcpars,model,0);
    if err<80
        pars2 = update_pars(pars1,chain(i,:),mcmcpars);
        [t2,S2,V2,D2] = simulate_ode(model,pars2,tvec,pars2.S0,pars2.V0); % initial parameter set
        for j=1:NH
            subplot(2,NH,j)
            patchline(t2,S2(:,j),'edgecolor',[0.5 0.5 0.5],'edgealpha',0.2);hold on;
            subplot(2,NH,j+NH)
            patchline(t2,V2(:,j),'edgecolor',[0.5 0.5 0.5],'edgealpha',0.2);hold on;
        end

    else

    end
   
end
% 
% %% posteriors
% 
figure(10)
for i = 1:9
    subplot(2,5,i)
    smoothHistogram(chain(300:1000,i),10); hold on;
    xlabel('Burst size (virions/cell)');
    hold on;
    plot(0:1:1000,gaussian(0:1:1000,mcmcparam{1,i}{1,5},mcmcparam{1,i}{1,6}),'-k',LineWidth=2 );
    xlim([0 500]);
    xline(chain(1,i),'-r',LineWidth=2)
    set(gca,'FontSize',14);
    ylabel('PDF');
end
    


figure(11)
for i = 1:9
    subplot(2,5,i)
    smoothHistogram(chain(300:1000,i+9),10); hold on;
    xlabel('log \phi (ml/cell/hr)');
    hold on;
    plot(-10:0.1:-7.5,gaussian(-10:0.1:-7.5,mcmcparam{1,i+9}{1,5},mcmcparam{1,i+9}{1,6}),'-k',LineWidth=2 );
    xline(chain(1,i+9),'-r',LineWidth=2)
    %ylim([0 1e-2])
    set(gca,'FontSize',14);
    ylabel('PDF');
end
    

figure(12)
for i = 1:9
    subplot(2,5,i)
    smoothHistogram(chain(300:1000,i+18),5); hold on;
    xlabel('Latent period (hrs)');
    hold on;
    plot(0:0.1:10,gaussian(0:0.1:10,mcmcparam{1,i+18}{1,5},mcmcparam{1,i+18}{1,6}),'-k',LineWidth=2 );
    xline(chain(1,i+18),'-r',LineWidth=2)
    set(gca,'FontSize',14);
    
    ylabel('PDF');
end

figure(13)
for i = 1:6
    subplot(1,6,i)
    smoothHistogram(chain(300:1000,i+27),5); hold on;
    xlabel('Host growth rate (cells/hr)');
    hold on;
    plot(0:0.01:0.5,gaussian(0:0.01:0.5,mcmcparam{1,i+27}{1,5},mcmcparam{1,i+27}{1,6}),'-k',LineWidth=2 );
    xline(chain(1,i+27),'-r',LineWidth=2)
    xlim([0 0.5]);
    set(gca,'FontSize',14);
    ylabel('PDF');
end


