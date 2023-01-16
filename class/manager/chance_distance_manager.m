classdef chance_distance_manager < manager_interface & chance_manager_interface
    %CHANCE_DISTANCE_MANAGER manager for probabalistic (chance-based) 
    %distance estimation
    
    methods
        function obj = chance_distance_manager(chance_supp, dyn, objective)
            %CHANCE_PEAK_MANAGER Construct an instance of this class
            if nargin == 2
                objective = [];
            end

            loc_curr = location_distance_sde(chance_supp, dyn, objective);            
            obj@manager_interface(loc_curr);
            obj@chance_manager_interface();  
            
            obj.MAXIMIZE = 1;
        end
        
        function [sol, obj] = run(obj, order, Tmax)
            %take the negative of the distance
            if nargin < 3
                Tmax = 1;
            end
            [sol, obj] = run@manager_interface(obj, order, Tmax);
            
            sol.obj_rec = -sol.obj_rec;
        end
    end
end

