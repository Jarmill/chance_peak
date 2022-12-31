classdef chance_support < loc_support
    %CHANCE_SUPPORT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        epsilon; %probability bound
                 %peak is achieved at probability (1-epsilon)

        bound_type = 'cantelli'; %either mean, Cantelli or VP bounds allowed (in this implementation)
        %'mean':  bound the max mean of the distribution, no fancy math
        %'cantelli': Cantelli bound on value-at-risk
        %'vp': Vysochanskij-Petunin bound on value-at-risk (requires unimodal
        %distributions and epsilon < 1/6)
    end
    
    methods
        function obj = chance_support(vars, epsilon, loc_ref)
            %CHANCE_SUPPORT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                epsilon = 0.1;
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
            
            if strcmp(obj.bound_type, 'vp')
                r = sqrt(4/(9*epsilon) - 1);
            elseif strcmp(obj.bound_type, 'cantelli')
                r = sqrt(1/(epsilon) - 1);
            else %default to mean (effectively epsilon = 0.5)
                r=0;
            end
            
        end
    end
end

