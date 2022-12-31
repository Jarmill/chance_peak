classdef chance_peak_manager < peak_manager
    %CHANCE_PEAK_MANAGER manager for probabalistic (chance-based) peak
    %estimation. Next step: add SDE dynamics
    
%     properties
%         Property1
%     end
    
    methods
        function obj = chance_peak_manager(chance_supp, f, objective)
            %CHANCE_PEAK_MANAGER Construct an instance of this class
            if nargin == 2
                objective = 0;
            end

            obj@peak_manager(chance_supp, f, objective);
                                 
        end
        
        function obj = dual_process(obj, d, dual_rec, len_dual)
            %TODO: need to perform dual recovery at some point
            obj = dual_process@peak_manager(d, dual_rec, len_dual)
        end
        

    end
end

