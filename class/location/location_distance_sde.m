classdef location_distance_sde < location_sde & location_distance_interface
    %LOCATION_DISTANCE_SDE 
    %perform distance estimation with an SDE
    

    methods
        function obj = location_distance_sde(unsafe_supp, dyn, objective, id)
            %LOCATION_SDE Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3             
                %by default, no objective
                objective = [];                
            end
            
            if nargin < 4
                id = [];            
            end
            
            objective = 0;
            
            obj@location_sde(unsafe_supp, dyn, objective, id);
            obj@location_distance_interface(unsafe_supp, id);
            
%             obj.varnames = [obj.varnames, {'y'}];
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
           
            [objective, cons_eq, cons_ineq, len_dual] = all_cons@location_sde(obj, d);
            
            marg = obj.marg_wass_con(d);
            
            cons_eq = [cons_eq; marg==0];

            %package up the output
            len_dual.w = length(marg);                       
        end      
        
        function [len_out] = len_eq_cons(obj)
            %LEN_EQ_CONS Number of equality constraints strictly in this
            %location            
            len_out = obj.len_dual.v;
        end
        
        function [obj_max, obj_con_ineq, obj_con_eq] = objective_con(obj, d, objective)
            %OBJECTIVE_CON deal with the distance objective
            %bad code, violates Don't Repeat Yourself principles
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
            if nargin == 2
                objective = obj.get_objective();
            end
                                    
            obj_con_eq = [];
            obj_con_ineq = [];
            
            %extra terms for the chance-peak
            %p(x) + r*sqrt(p(x)^2 - p(x))
            obj_subs_2 = 0;
            
%             var_dist = obj.var_index(obj.vars, {'x', 'y'});
            var_dist = [obj.vars.x; obj.vars.y];
            if isempty(objective)
                obj_max = [0; 0];
            elseif length(objective) == 1    
                obj_subs = obj.wass{1}.mom_objective(objective, var_dist);                
%                 obj_max = (obj_subs);                            
                
                obj_subs_2 = obj.wass{1}.mom_objective(objective^2, var_dist);
                
                obj_max = [obj_subs; obj_subs_2];
            else
                obj_subs = obj.wass{1}.mom_objective(objective, var_dist);
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

        
        function obj_out = get_objective(obj)
            %return the objective to maximize (this is a chance-peak
            %maximizing SDE)
            obj_out = -obj.dist;
        end
        
        function supp_con_out = supp_con(obj)
            %SUPP_CON get support constraints of measures
            
            
            %terminal measure support 
            supp_con_out = supp_con@location_sde(obj);
            wass_con = supp_con@location_distance_interface(obj);
            
            supp_con_out = [supp_con_out; wass_con];
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
             monom = mmon(obj.get_vars, 0, d);
             obj.dual.v = v_coeff'*monom;
             
             %iterate through all subsystems
             monom_all = mmon(obj.get_vars(), 0, d);
             
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
        
        function [optimal, mom_out, corner] = recover(obj, tol)
            %recovery: do this later
           if nargin < 2
               tol = 0;
           end
           optimal = 0;
           mom_out = [];
           corner = [];
        end
        %% Sampling ???
    end
end

