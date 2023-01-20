classdef location_sde < location_interface
    %LOCATION_SDE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        varnames = {'t', 'x'};
        g = [];
        jump = [];
    end
    
    methods
        function obj = location_sde(loc_supp, dyn, objective, id)
            %LOCATION_SDE Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3             
                %by default, no objective
                objective = [];                
            end
            
            if nargin < 4
                id = [];            
            end

            if ~iscell(dyn.f)
                dyn.f = {dyn.f};
            end

            if ~iscell(dyn.g)
                dyn.g= {dyn.g};
            end

            obj@location_interface(loc_supp, dyn.f, objective, id);
           
            
            if ~loc_supp.TIME_INDEP % && obj.supp.SCALE_TIME
                %scale for time if time is a variable
                Tmax = loc_supp.Tmax;
                for i = 1:length(dyn.f)
                    dyn.f{i} = Tmax * dyn.f{i};
                    dyn.g{i} = Tmax * dyn.g{i};
                end
                loc_supp.Tmax = 1;
            end
            %only a single SDE system

            obj.sys = cell(length(dyn.f), 1);
            for i = 1:length(dyn.f)
                dyn_curr = struct('f', dyn.f{i}, 'g', dyn.g{i});
                obj.sys{i} = subsystem_sde(loc_supp, dyn_curr, i, id);
            end
        end
        
        %% Constraints
        
        function [objective, cons_eq, cons_ineq, len_dual] = all_cons(obj, d)
            %LOC_CONS all constraints involving solely location
            %does not the objective anymore
            %
            %Output:
            %   cons_eq: equality constraints
            %   cons_ineq: inequality constraints (objective)
            
            %gather all constraints together
            liou = obj.liou_con(d);
            len_liou = length(liou);
           
            [objective, cons_ineq, cons_eq] = obj.objective_con();
            
            cons_ineq = [];

            %package up the output
            len_dual = struct;
            len_dual.v = len_liou;
%             len_dual.beta = length(cons_ineq);
            
            %ensure this is the correct sign
%             cons_eq = [-liou]==0;   
            cons_eq = [(-liou)==0; cons_eq];
        end      
        
        function [len_out] = len_eq_cons(obj)
            %LEN_EQ_CONS Number of equality constraints strictly in this
            %location            
            len_out = obj.len_dual.v;
        end
        
        function vars_out = get_vars_end(obj)
            %GET_VARS_END variables at endpoint measures
            %   initial and terminal, without time-dependent
            vars_out = [obj.vars.t; obj.vars.x];
        end
        
        
        
        function [obj_max, obj_con_ineq, obj_con_eq] = objective_con(obj, objective)
            %OBJECTIVE_CON deal with the objective, which may be maximin
            %                        
            %Inputs:
            %   objective:  the objective that should be maximized
            %
            %Outputs:
            %   obj_max:        a 2x1 vector with [objective; objective^2]
            %                   used for the chance-peak SOCP
            %   obj_con_ineq:   inequality constraints
            %   obj_con_eq:     equality constraints
            %
            %
            %TODO: This should maybe go in the manager
            %The current implementation is only for peak estimation
            
            %TODO: include support for putting objectives on initial and
            %occupation measures as well as the terminal measure
            if nargin == 1
                objective = obj.objective;
            end
                                    
            obj_con_eq = [];
            obj_con_ineq = [];
            
            %extra terms for the chance-peak
            %p(x) + r*sqrt(p(x)^2 - p(x))
            obj_subs_2 = 0;
            
            var_end = obj.var_index(obj.vars, {'t', 'x'});
            if isempty(objective)
                obj_max = [0; 0];
            elseif length(objective) == 1    
                obj_subs = obj.term.var_sub_mom(var_end, objective);                
%                 obj_max = (obj_subs);                            
                
                obj_subs_2 = obj.term.var_sub_mom(var_end, objective^2);
                
                obj_max = [obj_subs; obj_subs_2];
            else
                obj_subs = obj.term.var_sub_mom(var_end, objective);
                q_name = ['q_', num2str(obj.id)];
                mpol(q_name, 1, 1);
                q = eval(q_name);
                muq = meas(q);
                obj.cost_q = q;
                
                obj_max = mom(q);
                obj_con_eq = [mass(q) == 1];
                obj_con_ineq= (mom(q) <= obj_subs);
%                 obj_subs_2 = mom(q^2);
                
                obj_max = [mom(q); mom(q^2)];
            end            
            
            %call the chance-peak function
            

        end

        
        function supp_con_out = supp_con(obj)
            %SUPP_CON get support constraints of measures
            
            
            %terminal measure support 
            if ~isempty(obj.term)
                term_supp =  obj.term.supp();
            else
                term_supp =  [];
            end
            
            %initial measure support 
            if ~isempty(obj.init)
                init_supp =  obj.init.supp();
            else
                init_supp =  [];
            end
            
            %subsystem measure support
            sys_supp = [];
            for i = 1:length(obj.sys)
                sys_supp = [sys_supp; obj.sys{i}.get_supp()];
            end
            
            supp_con_out = [init_supp;
                            term_supp;
                            sys_supp];
        end
        
        
        function obj = dual_process(obj, d, rec_eq, rec_ineq, gamma)
             %DUAL_PROCESS turn the dual variables from solution into 
             %polynomials and interpretable quantities
             %
             %Input:
             %  d:          2*order, order of auxiliary polynomials
             %  rec_eq:     dual variables from equality constraints
             %  rec_ineq:   dual variables from inequality constraints
             %  gamma:      objective value (as a dual variable)
             
             %Need to fix this for the SDE case
             
             %numeric quantities
             obj.dual.solved = 1;
             
             obj.dual.beta = rec_ineq;
             obj.dual.gamma = gamma;
             
             %process the polynomials
             
             %auxiliary function v
             v_coeff = rec_eq(1:obj.len_dual.v);
             monom = mmon(obj.get_vars_end, 0, d);
             obj.dual.v = v_coeff'*monom;
             
             %iterate through all subsystems
             monom_all = mmon(obj.get_vars_end(), 0, d);
             
             %TODO: confirm that all abscont relations have the same length
             len_monom_all = length(monom_all);

             %ship off dual variables for processing in subsystem 
             obj.sys{1} = obj.sys{1}.dual_process(obj.dual.v);       


            %nonnegativity of location (not subsystems)
            %initial measure
            if isempty(obj.init)
                nn_init = 0;
            else
                nn_init = obj.dual.gamma - obj.dual.v;
            end
            
            %terminal measure
            if isempty(obj.term)
                nn_term = 0;
            else
                nn_term = obj.dual.v - obj.dual.beta'*obj.objective;
            end
            obj.dual.nn = [nn_init; nn_term];
             
        end     
        
        
        %% Sampling ???
    end
end

