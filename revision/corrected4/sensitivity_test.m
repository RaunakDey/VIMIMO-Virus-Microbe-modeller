clc;
clear all;
tic

load('./v14-4.mat');

%% all converted to linear scale


theta_start = chain(1,:);

%convert it back to linear scale
theta_start(10:18) = 10.^chain(1,10:18);
theta_start(33:37) = 10.^chain(1,33:37);

% when simulate do not do log scale back conversion, so changing the
% properties.

mcmcpars_linear = mcmcpars;
mcmcpars_linear.phi.log =  0;
mcmcpars_linear.Dc.log =  0;
mcmcpars_linear.Dc2.log =  0;
mcmcpars_linear.Dc3.log =  0;
mcmcpars_linear.Dc4.log =  0;
mcmcpars_linear.Dc5.log =  0;


for i = 1:37
theta = theta_start;
i

dummy = theta_start(i) + (-0.3:0.1:0.3)*theta_start(i);
for j = 1:7
    theta(i) = dummy(j);
    error(i,j) = ssfun(theta,data,pars2,mcmcpars_linear,model,0);
end


end

%%
n = 280;
sigma_ll = 10;
LL = -n/2*log(2*pi*sigma_ll^2) - error/(sigma_ll^2);

%% optimal sigma_ll

x = 0:0.1:100;
plot(-100*log(2*pi*x) - 200./x,'-o');
%% Likelihood visualization
figure(2)
for i = 1:9
subplot(2,5,i)
plot(theta_start(i) + (-0.3:0.1:0.3)*theta_start(i),LL(i,:),'-bo',LineWidth=2);
xlabel('\beta');
ylabel('LL')
set(gca,'FontSize',20)
end


figure(3)
for i = 10:18
subplot(2,5,i-9)
plot(theta_start(i) + (-0.3:0.1:0.3)*theta_start(i),LL(i,:),'-bo',LineWidth=2);
xlabel('\phi');
ylabel('LL');
set(gca,'FontSize',14)
end


figure(4)
for i = 19:27
subplot(2,5,i-18)
plot(theta_start(i) + (-0.3:0.1:0.3)*theta_start(i),LL(i,:),'-bo',LineWidth=2);
xlabel('\tau');
ylabel('LL')
set(gca,'FontSize',20)
end

%% compute prior std

% I am going to assume gaussian, which is terribly wrong!
% but should work for providing a rough estimate.

for i=1:36
std_est(i) = 1/sqrt( 2*log(mean([error(i,1),error(i,7)])/error(i,4)) / (0.6*theta_start(i))^2)

end

%%

std_est = abs(std_est);

for i = 1:37
 i

theta_right = theta_start;
theta_left = theta_start;
theta_right(i) = theta_start(i) + std_est(i);
theta_left(i) = theta_start(i) - std_est(i);
err_right(i) =  ssfun(theta_right,data,pars2,mcmcpars_linear,model,0);
err_left(i) = ssfun(theta_left,data,pars2,mcmcpars_linear,model,0);
end



