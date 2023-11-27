classdef SEIVD_diff_NE_diff_debris_abs_diffbeta < ode_funs
    
    properties
        name;
        statevars = ["S","E","I","V"];
        NH;
        NV;
        NE;
        id;
    end
    
    methods
        
        function obj = SEIVD_diff_NE_diff_debris_abs_diffbeta(NH,NV,NE)
            obj.name = sprintf('SEIV%d',NE);
            obj.NH = NH;
            obj.NV = NV;
            obj.NE = NE;
            id.S = 1:NH;
            id.E = (1:NH*NV*NE)+id.S(end);
            id.I = (1:NH*NV)+id.E(end);
            id.V = (1:NV)+id.I(end);
            id.D = 1+id.V(end); % debris
            id.Emat = reshape(id.E,[NH NV NE]);
            id.Imat = reshape(id.I,[NH NV]);
            obj.id = id;
        end

        function vec = zeros(obj)
            %vec = zeros(obj.id.V(end),1);
            vec = zeros(obj.id.D,1);
        end
        
        function S = sum_hosts(obj,y) % all host cells including infected
            S = y(:,obj.id.S);
            for i = 1:obj.NH % could use bsxfun or arrayfun here?
               S(:,i) = S(:,i) + sum(y(:,obj.id.Emat(i,:)),2) + sum(y(:,obj.id.Imat(i,:)),2); 
            end
        end
        
        function V = sum_viruses(obj,y) % free phage only
           V = y(:,obj.id.V); 
        end
        
        function dydt = ode(obj,t,y,pars)
            % so we have the variable t here already.

            OH = ones(obj.NH,1);
            OV = ones(obj.NV,1);
            S = y(obj.id.S);
            Emat = y(obj.id.Emat);
            Imat = y(obj.id.Imat);
            V = y(obj.id.V);
            D = y(obj.id.D);
            N = S+sum(Emat,3)*OV+Imat*OV;
            etaeff = pars.eta*(pars.NE+1);

            exposed_transition_fun = obj.exposed_transition_fun;
            host_growth_fun = obj.host_growth_fun;
            viral_decay_fun = obj.viral_decay_fun;
            viral_adsorb_fun = obj.viral_adsorb_fun; % only use in dV
            lysis_reset_fun = obj.lysis_reset_fun;
            viral_debris_interaction = obj.viral_debris_interaction;
            diffbeta_function = obj.diffbeta_function;

            debris_inhib_fun = obj.debris_inhib_fun;
            debris_inhib_fun_second = obj.debris_inhib_fun_second;
            debris_inhib_fun_third = obj.debris_inhib_fun_third;
            debris_inhib_fun_fourth = obj.debris_inhib_fun_fourth;
            debris_inhib_fun_fifth = obj.debris_inhib_fun_fifth;


                  
            Sdeb(1) =  y(1)*debris_inhib_fun(pars,D); % debris
            Sdeb(2) =  y(2)*debris_inhib_fun_second(pars,D);
            Sdeb(3) =  y(3)*debris_inhib_fun_third(pars,D);

             %Sdeb(1) =  y(1);
             %Sdeb(2) =  y(2);
             %Sdeb(3) =  y(3);

            Sdeb(4) =  y(4)*debris_inhib_fun_fourth(pars,D);
            Sdeb(5) =  y(5)*debris_inhib_fun_fifth(pars,D);
            
            %Sdeb = y;
            

            Ndeb(1) =  N(1)*debris_inhib_fun(pars,D); % debris
            Ndeb(2) =  N(2)*debris_inhib_fun_second(pars,D);
            Ndeb(3) =  N(3)*debris_inhib_fun_third(pars,D);
            
%             Ndeb(1) =  N(1);
%             Ndeb(2) =  N(2);
%             Ndeb(3) =  N(3);

            Ndeb(4) =  N(4)*debris_inhib_fun_fourth(pars,D);
            Ndeb(5) =  N(5)*debris_inhib_fun_fifth(pars,D);
            
                       
            %Sdeb = S.*debris_inhib_fun(pars,D); % debris
            Sdeb = Sdeb';
            
            dS = host_growth_fun(pars,S,N) - Sdeb.*((pars.M.*pars.phi)*V);
            dEmat = (pars.M.*pars.phi).*(Sdeb*V') - etaeff.*Emat(:,:,1) + lysis_reset_fun(pars,Imat,OH,V);
            dEmat2 = exposed_transition_fun(etaeff,Emat);
            for i = 1:obj.NH
                for j = 1:obj.NV
                    if (pars.NE(i,j) ~= obj.NE  && pars.NE(i,j) ~=0) 
                        Emat(i,j,pars.NE(i,j)+1:end) = Emat(i,j,pars.NE(i,j));
                    end
                end
            end
            
            dImat = etaeff.*Emat(:,:,end) - etaeff.*Imat - lysis_reset_fun(pars,Imat,OH,V);
            dV = (diffbeta_function.*etaeff.*Imat)'*OH - V.*(viral_adsorb_fun(pars)'*Ndeb') - viral_decay_fun(pars,V) - viral_debris_interaction(pars,V,D);
            dD = sum(etaeff(:).*Imat(:)); % sum across all pairs for net lysis rate
            
            dydt = [dS; dEmat(:); dEmat2(:); dImat(:); dV; dD];
            
        end
    end
end

