classdef subsystem_sde <subsystem_interface
    %SUBSYSTEM A subsystem dx=f(t, x)dt + g(t, x)dw of a stochastic 
    %differential equation where there is no additional uncertainty
    
    properties
        varnames = {'t', 'x'};
        
        %sde properties
        g = 0;
        jump = [];
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
        end
        
        %% Constraints                       
        function Ay = cons_liou(obj, d)
            %CONS_LIOU Liouville Equation includes an affine combination of
            %Lie derivatives (continuous systems only)
            
            %no box inputs, simple to perform
             Ay_f = obj.meas_occ.mom_lie(d, obj.get_vars, obj.f);
             Ay_g = obj.meas_occ.mom_hess(d, obj.get_vars, obj.g);
             
             %Lv = v_t + f' grad v + 0.5 g' hess v g
             Ay = Ay_f + 0.5*Ay_g;
           
            
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

