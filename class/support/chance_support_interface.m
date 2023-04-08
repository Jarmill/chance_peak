classdef (Abstract) chance_support_interface
    %CHANCE_SUPPORT_INTERFACE Allow for a location support to contain a
    %chance-constraint
    %   Detailed explanation goes here
    
    properties
        epsilon; %probability bound
                 %peak is achieved at probability (1-epsilon)

        bound_type = 'cantelli'; %either mean, Cantelli or VP bounds allowed (in this implementation)
        %'mean':  bound the max mean of the distribution, no fancy math
        %'cantelli': Cantelli bound on value-at-risk
        %'vp': Vysochanskij-Petunin bound on value-at-risk (requires unimodal
        %distributions and epsilon < 1/6)
        
        p_supp = []; %bounds on p inside X (if needed)
    end
    
    methods
        function obj = chance_support_interface(epsilon)
            %CHANCE_SUPPORT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 1
                epsilon = 0.1;
            end
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

