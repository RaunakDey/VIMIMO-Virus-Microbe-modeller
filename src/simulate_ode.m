function [t, S, V, D, I,E] = simulate_ode(model, pars, tvec, S0, V0)

% set up initial conditions
y0 = model.zeros();
y0(model.id.S) = S0;
y0(model.id.V) = V0;

host_growth = model.host_growth ;
viral_decay = model.viral_decay ;
viral_adsorb = model.viral_adsorb ;
lysis_reset = model.lysis_reset;
debris_inhib = model.debris_inhib;
debris_inhib2 = model.debris_inhib2;
debris_inhib3 = model.debris_inhib3;
debris_inhib4 = model.debris_inhib4;
debris_inhib5 = model.debris_inhib5;



diff_beta = model.diff_beta;
NE  = round(max(max(pars.NE)));

if model.name == 'SEIV'+string(model.NE)
    model = SEIV_diff_NE(model.NH,model.NV,NE);
elseif strcmp(model.name,'SEIVD-diffabs-diffbeta')
    model = SEIVD_diff_NE_diff_debris_abs_diffbeta(model.NH,model.NV,NE);
elseif strcmp(model.name,'SEIVD-diffabs')
    model = SEIVD_diff_NE_diff_debris_abs(model.NH,model.NV,NE);
else
    model = SEIVD_diff_NE_diff_debris(model.NH,model.NV,NE);

end

model.host_growth = host_growth;
model.viral_decay = viral_decay;
model.viral_adsorb = viral_adsorb;
model.lysis_reset = lysis_reset;
model.debris_inhib = debris_inhib;
model.debris_inhib2 = debris_inhib2;
model.debris_inhib3 = debris_inhib3;
model.debris_inhib4 = debris_inhib4;
model.debris_inhib5 = debris_inhib5;


model.diff_beta = diff_beta;



% integrate
reltol = 1e-8;
opts = odeset('RelTol',reltol);
ode = @(t,y) model.ode(t,y,pars);
[t, y] = ode45(ode,tvec,y0,opts);

% output timeseries + compute total host and virus density
S = model.sum_hosts(y);
V = model.sum_viruses(y);
D = y(:,model.id.D); % debris
I = y(:,model.id.I); % infections. % should be a 5*5 so 25 objects.
E = y(:,model.id.E);

% measurement bias
if ~all(isnan(pars.epsilon))
    S = S.*pars.epsilon(1:model.NH);
    V = V.*pars.epsilon(model.NH+(1:model.NV));
end

end