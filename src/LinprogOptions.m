classdef LinprogOptions
    properties
        
        Algorithm;
        Display;
        ConstraintTolerance = 1e-06;
        MaxIterations = 200;
        OptimalityTolerance = 1e-06;
        MaxTime = Inf;
        
    end
    
    methods
        
        function obj = LinprogOptions(Algorithm,Display)
            obj.Algorithm = Algorithm;
            obj.Display = Display;
        end
        
    end

end