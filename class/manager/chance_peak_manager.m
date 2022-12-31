classdef chance_peak_manager < manager_interface
    %CHANCE_PEAK_MANAGER manager for probabalistic (chance-based) peak
    %estimation. Next step: add SDE dynamics
    
%     properties
%         Property1
%     end
    
    methods
        function obj = chance_peak_manager(chance_supp, dyn, objective)
            %CHANCE_PEAK_MANAGER Construct an instance of this class
            if nargin == 2
                objective = [];
            end

            loc_curr = location_sde(chance_supp, dyn, objective);
            
            obj@manager_interface(loc_curr);
                                 
        end
        
        function [objective, mom_con, supp_con, len_dual] = cons(obj,d, Tmax)
            %PEAK_CONS formulate support and measure constraints for peak
            %program at degree d
            %Input:
            %   d:      Monomials involved in relaxation (2*order)
            %   Tmax:   Maximum time (only when time-independent)
            %
            %Output:
            %   objective:  target to maximize  (@mom)
            %   mom_con:    moment constraints  (@momcon)
            %   supp_con:   support constraints (@supcon)
          

            supp_con = obj.loc.supp_con();       %support constraint                 
                        
            %gather all constraints in each location
            %(loop over locations)
            [obj_max, cons_eq, cons_ineq, len_dual] = ...
                obj.loc.all_cons(d);
            
            %objective relations (specific care needs to be taken for
            %chance-peak)
            [objective, cons_eq_obj] = obj.objective_process(obj_max);
            
            cons_eq = [cons_eq; cons_eq_obj];
            
            %finalize moment constraints
            
            %mass of initial measure sums to one
            mass_init_con = (obj.loc.mass_init() - 1 == 0);
            if islogical(mass_init_con)
                mass_init_con = [];
            end
            
            if obj.loc.supp.TIME_INDEP
                mass_occ_con = (obj.loc.mass_occ() <= Tmax);
            else
                mass_occ_con = [];
            end

            %time independent: mass of sum of occupation measures are less
            %than Tmax. Implement/fix this?
            %TODO: get this done
                
%             len_liou = length(liou_con);

            mom_con = [cons_eq; cons_ineq; mass_init_con; mass_occ_con];            
        end    
        
        function [objective, cons_eq_obj] = objective_process(obj, obj_max)
            %OBJECTIVE_PROCESS implement the SOC constraints for the
            %chance-peak bound
            
            p = obj_max(1);
            p2 = obj_max(2);
            r = obj.loc.supp.get_r_const();
            
            if r > 0
                %there is some funky chance-peak business going on
                %need to add the SOC constraint now
                
                %sup z*r + p
                %norm([1-p^2, 2z, 2p], 2) <= 1+p^2
                
                %SDP form
                %
                %[1+p^2,  1-p^2,   2z,  2p;
                % 1-p^2,      1,    0,   0;
                %    2z,      0     1,   0;
                %    2p,      0,    0,   1]  is PSD;
                
                mpol('s', 3)
                mus = meas(s);
                
                svector = reshape(mom(s*s'), [], 1);
                seye = reshape(eye(3), [], 1);
                cons_eq_obj = [(mass(mus)==1+p2);
                                mom(s(1)) == 1-p2;
                                mom(s(3)) == 2*p;
                                svector==seye];
                            
                objective = obj_max(1) + mom(s(2))*(r/2);
            else
                %the mean is being maximized
                %no funky business, this remains an LP
                objective = obj_max(1);
                cons_eq_obj = [];
            end
            
        end
        
        function obj = dual_process(obj, d, dual_rec, len_dual)
            %TODO: need to perform dual recovery at some point
%             obj = obj;
        end
        

    end
end

