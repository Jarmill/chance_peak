classdef chance_distance_support < loc_support & chance_support_interface & unsafe_support_interface
    %CHANCE_distance_SUPPORT Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = chance_distance_support(vars, epsilon, loc_ref)
            %CHANCE_SUPPORT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                epsilon = 0.1;
            end
            if nargin < 3
                loc_ref = [];
            end
            obj@loc_support(vars, loc_ref);
            obj@chance_support_interface(epsilon);
            obj@unsafe_support_interface(vars, loc_ref);
        end
       
    end
end

