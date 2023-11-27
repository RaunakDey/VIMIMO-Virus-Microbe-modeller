classdef ode_funs
    % every model is a subclass of ode_funs
    
    % ode_funs contains variations on particular parts of the core
    % virus-host differential equations (see properties)
    
	% to use a function handle inside an ode function of an ode_funs subclass, do:
	%   fun = obj.fun               %(this "unpacking" line is necessary)
	%   dX = fun(inputs) + ...
    
    properties
        host_growth = 0; % 0 = exponetial growth; 1 = global carrying capacity w/ niche competition
        viral_decay = 0; % 0 = no decay; 1 = exponential decay
        viral_adsorb = 0; % 0 = viruses only adsorb to their hosts; 1 = viruses may adsorb to any host
        lysis_reset = 0; % 0 = no lysis reset; 1 = adsorption by new virsues resets lytic cycle from I to first E class
        debris_inhib = 0; % 0 = no debris inhibition; 1 = build-up of dead cells inhibits infection (requires defining pars.Dc)
        debris_inhib2 = 0;
        debris_inhib3 = 0;
        debris_inhib4 = 0;
        debris_inhib5 = 0;

        diff_beta = 0;
        %NH = 5; %% have to remove later -- this is default value
    end
    
    methods
        
        function mystr = odestr(obj)
            
            mystr = sprintf('%d%d%d%d%d',obj.host_growth,obj.viral_decay,obj.viral_adsorb,obj.lysis_reset,obj.debris_inhib);
        end
        
        function handle = viral_growth(obj)
            if obj.diff_beta == 0
                handle = @(pars,t, Imat, OH, etaeff)  (pars.beta.*etaeff.*Imat)'* OH;
            else
                handle = @(pars,t, Imat, OH, etaeff)  ((pars.beta + (pars.beta2 - pars.beta).*heaviside(t-2.6*pars.tau).*pars.M) .*etaeff.*Imat)'* OH;
            end
        end



        function handle = exposed_transition_fun(obj)
            % output = dEmat2 (NH x NV x NE-1)
            % use depends on model structure, user does not interact with this ftn
            if obj.NE==1 % NE exists for every model subclass (if this throws an error, check subclass constructor)
               handle = @(etaeff,Emat) [];
           else
               handle = @(etaeff,Emat) etaeff.*Emat(:,:,1:end-1) - etaeff.*Emat(:,:,2:end);
           end
        end

        function handle = host_growth_fun(obj)
            % output = component of dS (NH x 1)
            if obj.host_growth==0
                handle = @(pars,S,N) pars.r.*S;
            else
                handle = @(pars,S,N) pars.r.*S.*(1-pars.a*N/pars.K);
            end
        end
        
        function handle = viral_decay_fun(obj)
            % output = component of dV (NV x 1)
            % also compatible with W state variable, then output = component of dW (NV x 1)
            if obj.viral_decay==0
                handle = @(pars,V) 0;
            else
                handle = @(pars,V) pars.m.*V;
            end
        end
         
        function handle = viral_adsorb_fun(obj)
            % output = adsorption matrix (NH x NV)
            % viral_adsorb_fun is only used in dV equation; assuming adsorption by non-specific viruses have no affect on hosts
            if obj.viral_adsorb==0
                handle = @(pars) pars.M.*pars.phi;
            else
                handle = @(pars) pars.phi;
            end
        end       

        function handle = lysis_reset_fun(obj)
            % output = component of dImat (to be subtracted) and dEmat (to be added) (NH x NV)
            if obj.lysis_reset==0
                handle = @(pars,Imat,OH,V) 0;
            else 
                viral_adsorb_fun = obj.viral_adsorb_fun;
                handle = @(pars,Imat,OH,V) pars.epsilon_reset.*viral_adsorb_fun(pars).*Imat.*(OH*V');
            end
        end
        
        function handle = viral_debris_interaction(obj)
            if obj.debris_inhib == 2
                handle = @(pars,V,D) (pars.prob.*(pars.phi*V))*D; 
            else
                handle = @(pars,V,D) 0;
            end
        end

        function handle = debris_inhib_fun(obj)
            % output = scalar to be multiplied with state variable S
            if obj.debris_inhib==0
                handle = @(pars,D) 1;
            elseif obj.debris_inhib == 1
                handle = @(pars,D) 1/(1+D/pars.Dc)^2;
            elseif obj.debris_inhib == 2
                handle = @(pars,D) 1./(1+(D./pars.Dc).^2); %exponent corrected
                %handle = @(pars,D) (1+0/(1+(pars.Dc/D)^2));
            elseif obj.debris_inhib == 3
                handle = @(pars,D) 1/(1+(pars.Dc/D)^2); %exponent corrected and reciprocal

            end
        end


        function handle = debris_inhib_fun_second(obj)
            if obj.debris_inhib2==0
                handle = @(pars,D) 1;
            elseif obj.debris_inhib2 == 1
                handle = @(pars,D) 1/(1+D/pars.Dc2)^2;
            elseif obj.debris_inhib2 == 2
                handle = @(pars,D) 1/(1+(D/pars.Dc2)^2); %exponent corrected
                %handle = @(pars,D) (1+0/(1+(pars.Dc2/D)^2));
            elseif obj.debris_inhib2 == 3
                handle = @(pars,D) 1/(1+(pars.Dc2/D)^2); %exponent corrected and reciprocal

            end
        end


        function handle = debris_inhib_fun_third(obj)
            if obj.debris_inhib3==0
                handle = @(pars,D) 1;
            elseif obj.debris_inhib3 == 1
                handle = @(pars,D) 1/(1+D/pars.Dc3)^2;
            elseif obj.debris_inhib3 == 2
                handle = @(pars,D) 1/(1+(D/pars.Dc3)^2); %exponent corrected
                %handle = @(pars,D) (1+0/(1+(pars.Dc3/D)^2));
            elseif obj.debris_inhib3 == 3
                handle = @(pars,D) 1/(1+(pars.Dc3/D)^2); %exponent corrected and reciprocal

            end
        end
  

        function handle = debris_inhib_fun_fourth(obj)
            if obj.debris_inhib4==0
                handle = @(pars,D) 1;
            elseif obj.debris_inhib4 == 1
                handle = @(pars,D) 1/(1+D/pars.Dc4)^2;
            elseif obj.debris_inhib4 == 2
                handle = @(pars,D) 1/(1+(D/pars.Dc4)^2); %exponent corrected
            elseif obj.debris_inhib4 == 3
                handle = @(pars,D) 1/(1+(pars.Dc4/D)^2); %exponent corrected and reciprocal

            end
        end

        function handle = debris_inhib_fun_fifth(obj)
            if obj.debris_inhib5==0
                handle = @(pars,D) 1;
            elseif obj.debris_inhib5 == 1
                handle = @(pars,D) 1/(1+D/pars.Dc5)^2;
            elseif obj.debris_inhib5 == 2
                handle = @(pars,D) 1/(1+(D/pars.Dc5)^2); %exponent corrected
            elseif obj.debris_inhib5 == 3
                handle = @(pars,D) 1/(1+(pars.Dc5/D)^2); %exponent corrected and reciprocal

            end
        end

    end
    


end




