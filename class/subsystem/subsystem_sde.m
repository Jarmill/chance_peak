classdef subsystem_sde <subsystem_interface
    %SUBSYSTEM A subsystem dx=f(t, x)dt + g(t, x)dw of a stochastic 
    %differential equation where there is no additional uncertainty
    
    properties
        varnames = {'t', 'x'};
        
        %sde properties
        g = 0;
        jump = [];
        lam_handle = @(d) 0;
        DISCRETE_TIME=0;
        Tmax = 1;
    end

    properties(Access = protected)
        %use private properties for numeric function evaluations
        %object-oriented matlab is slow with public variables 
        %
        %TODO: check if protected is also slow
        g_;
    end
    
    methods
        function obj = subsystem_sde(loc_supp, dyn, sys_id, loc_id)
            %SUBSYSTEM_SDE Construct an instance of this class
            %dyn: structure with fields 'f', 'g', 'jump'
            %the 'jump' part will allow for Levy processes (Poisson jumps),
            %and will be implemented later. Right now, f and g are the
            %important parts
            
            %process input
            if nargin < 3
                sys_id = 1;
            end
            
            if nargin < 4
                loc_id = [];
            end
            
            
            obj@subsystem_interface(loc_supp, dyn.f, sys_id, loc_id);
            
            if isfield(dyn, 'g')
                obj.g = dyn.g;
            end
            
            if isfield(dyn, 'jump')
                obj.jump = dyn.jump;
            end
            
            if isfield(obj.vars, 'y')
                obj.vars = rmfield(obj.vars, 'y');
            end

            obj.g_ = dyn.g;
            
            obj.DISCRETE_TIME = loc_supp.DISCRETE_TIME;
            
            obj.lam_handle = loc_supp.lam_handle;
            obj.Tmax = loc_supp.Tmax;
        end
        
        %% Constraints                       
        function Ay = cons_liou(obj, d)
            %CONS_LIOU Liouville Equation includes an affine combination of
            %Lie derivatives (continuous systems only)
            
            %no box inputs, simple to perform
            
            if obj.DISCRETE_TIME
                
                %terrible code, hack. fix later, put this in measures
                v = mmon([obj.vars.t; obj.vars.x], d);
                if ~isempty(obj.vars.t)                
                    Rv = subs(v, [obj.vars.t; obj.vars.x], [obj.vars.t+1/obj.Tmax; obj.f]);
                    Rv_int = integrate_var(Rv, obj.vars.lam, obj.lam_handle);
    %             f_curr = obj.var_sub(vars_old, f_old);
    %             Rv = subs(v, [obj.vars.t; obj.vars.x], ...
    %                 [obj.vars.t+t_shift; f_curr]);


                    Rv_int_subs = subs_vars(Rv_int, [obj.vars.t; obj.vars.x], [obj.meas_occ.vars.t; obj.meas_occ.vars.x]);
                    v_subs = subs_vars(v, [obj.vars.t; obj.vars.x], [obj.meas_occ.vars.t; obj.meas_occ.vars.x]);
    %             Rv_int = integrate_var(Rv, [obj.vars.lam], lam_handle);

                    %divide by Tmax to get sensible integrals
                    %elapsed time = 1 (when compressed to 1/Tmax)
                    Ay = (mom(Rv_int_subs)-mom(v_subs));

                    Ay = Ay *(1/obj.Tmax);
                else
                    Rv = subs(v, [obj.vars.x], [obj.f]);
                    Rv_int = integrate_var(Rv, obj.vars.lam, obj.lam_handle);
                    
                    Rv_int_subs = subs_vars(Rv_int, [obj.vars.x], [ obj.meas_occ.vars.x]);
                    v_subs = subs_vars(v, [obj.vars.x], [obj.meas_occ.vars.x]);
    %             Rv_int = integrate_var(Rv, [obj.vars.lam], lam_handle);

                    %divide by Tmax to get sensible integrals
                    %elapsed time = 1 (when compressed to 1/Tmax)
                    Ay = (mom(Rv_int_subs)-mom(v_subs));
                end
                
                
%                 Ay_push = obj.meas_occ.mom_push_integrate(d, obj.vars, ...
%                     obj.f, obj.lam_handle, 1/obj.Tmax);
%                 Ay_orig = obj.meas_occ.mom_monom(d);
%                 
%                 Ay = Ay_push - Ay_orig;
                    
            else
                Ay_f = obj.meas_occ.mom_lie(d, obj.get_vars, obj.f);
                Ay_g = obj.meas_occ.mom_hess(d, obj.get_vars, obj.g);

                %Lv = v_t + f' grad v + 0.5 g' hess v g
                Ay = Ay_f + 0.5*Ay_g;
            end
        
            
        end       

        %% sampler

        function g_out = g_eval(obj, data)
            %data: [t, x, th, w, b] as required            
            g_out = eval(obj.g_, obj.get_vars(), data);
        end
        
        
        %% Dual Recovery 
        function obj = dual_process(obj, v, zeta)
            %DUAL_PROCESS store dual functions and compute nonnegative
            %functions for this subsystem
            
            %auxiliary function v            
            obj.dual.v = v;
            
            vx = diff(v, obj.vars.x);
            vxx = diff(vx', obj.vars.x);
            
            obj.dual.Lv = vx*obj.f + obj.g'*vxx+obj.g;
            
            if ~isempty(obj.vars.t)
                obj.dual.Lv = obj.dual.Lv + diff(v, obj.vars.t);
            end
            
                                   
            obj.nn_ = obj.dual.nn;
        end    
    end
end

