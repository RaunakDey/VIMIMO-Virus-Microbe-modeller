close all;
clear all;
clc;

%% import datasets as matlab file

%load('./from_os_prior_20000.mat');

%load('./v14-4.mat')

load('./v17.mat');

%% add paths
addpath(genpath('./../../'));


%% construct the prior

% all of them are gaussian priors in log or linear scales
mu_prior = mcmcresults.prior(:,1);
sigma_prior = mcmcresults.prior(:,2);

transient_id = 100;

figure(1)

for i=1:9
subplot(2,5,i)
xaxis = linspace( mu_prior(i) - 3*sigma_prior(i),  mu_prior(i) + 3*sigma_prior(i) ,100);
plot(xaxis, ...
    gaussian(xaxis,mu_prior(i),sigma_prior(i)), ...
    "Color",'k',LineWidth=2);
hold on;
smoothHistogram(chain(transient_id:end,i),6,'b');
xlim([min(chain(transient_id:end,i)) max(chain(transient_id:end,i)) ]);
xlabel('\beta (virions/cell)');
ylabel('PDF');
end



figure(2)
for i=1:9
subplot(2,5,i)
i=i+9;
xaxis = linspace( mu_prior(i) - 3*sigma_prior(i),  mu_prior(i) + 3*sigma_prior(i) ,100);
plot(xaxis, ...
    gaussian(xaxis,mu_prior(i),sigma_prior(i)), ...
    "Color",'k',LineWidth=2);
hold on;
smoothHistogram(chain(transient_id:end,i),6,'b');
xlabel('\phi (ml/hr/cell)');
ylabel('PDF');
xlim([min(chain(transient_id:end,i)) max(chain(transient_id:end,i)) ]);
end




figure(3)

for i=1:9
subplot(2,5,i)
i=i+18;
xaxis = linspace( mu_prior(i) - 3*sigma_prior(i),  mu_prior(i) + 3*sigma_prior(i) ,100);
plot(xaxis, ...
    gaussian(xaxis,mu_prior(i),sigma_prior(i)), ...
    "Color",'k',LineWidth=2);
hold on;
smoothHistogram(chain(transient_id:end,i),6,'b');
xlabel('\tau (hr)');
ylabel('PDF');
xlim([min(chain(transient_id:end,i)) max(chain(transient_id:end,i)) ]);
end




figure(4)

for i=1:5
subplot(2,5,i)
i=i+27;
xaxis = linspace( mu_prior(i) - 3*sigma_prior(i),  mu_prior(i) + 3*sigma_prior(i) ,100);
plot(xaxis, ...
    gaussian(xaxis,mu_prior(i),sigma_prior(i)), ...
    "Color",'k',LineWidth=2);
hold on;
smoothHistogram(chain(transient_id:end,i),6,'b');
xlabel('r (/hr)');
ylabel('PDF');
xlim([min(chain(transient_id:end,i)) max(chain(transient_id:end,i)) ]);
end

for i=1:5
subplot(2,5,i+5)
i=i+32;
xaxis = linspace( mu_prior(i) - 3*sigma_prior(i),  mu_prior(i) + 3*sigma_prior(i) ,100);
plot(xaxis, ...
    gaussian(xaxis,mu_prior(i),sigma_prior(i)), ...
    "Color",'k',LineWidth=2);
hold on;
smoothHistogram(chain(transient_id:end,i),6,'b');
xlabel('Dc');
ylabel('PDF');
xlim([min(chain(transient_id:end,i)) max(chain(transient_id:end,i)) ]);
end

