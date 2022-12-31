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
            
            obj@location_interface(loc_supp, dyn, objective, id);
           
            %only a single SDE system
            obj.sys = subsystem_sde(loc_supp, dyn, 1, id);
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

