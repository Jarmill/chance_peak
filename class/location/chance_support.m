classdef chance_support < loc_support
    %CHANCE_SUPPORT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        epsilon; %probability bound
                 %peak is achieved at probability (1-epsilon)
                 
        bound_vp = 0; %either Cantelli or VP bounds allowed (in this implementation)
    end
    
    methods
        function obj = chance_support(vars, epsilon, loc_ref)
            %CHANCE_SUPPORT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                epsilon = 0.01;
            end
            if nargin < 3
                loc_ref = [];
            end
            obj@loc_support(vars, loc_ref);
            obj.epsilon = epsilon;
        end
        
        function r = get_r_const(obj, epsilon)
            %get the tail-bound constant (Cantelli or VP)
            if nargin < 2
                epsilon = obj.epsilon;
            end
            
            if obj.bound_vp
                r = sqrt(4/(9*epsilon) - 1);
            else
                r = sqrt(1/(epsilon) - 1);
            end
            
        end
    end
end

