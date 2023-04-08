classdef chance_peak_manager < manager_interface & chance_manager_interface
    %CHANCE_PEAK_MANAGER manager for probabalistic (chance-based) peak
    %estimation. Next step: add Levy dynamics
    
    methods
        function obj = chance_peak_manager(chance_supp, dyn, objective)
            %CHANCE_PEAK_MANAGER Construct an instance of this class
            if nargin == 2
                objective = [];
            end

            loc_curr = location_sde(chance_supp, dyn, objective);            
            obj@manager_interface(loc_curr);
            obj@chance_manager_interface();
                                 
        end
    end
end

