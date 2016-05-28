classdef (Abstract) FormationFlying
    
    properties (Abstract)
        time
        mu
        t0
        numPeriod
    end
    
    methods (Abstract)
        makeTimeVector(obj)
%         propagateModel(obj)
%         plotTrajectory(obj)
    end
    
    methods
        function obj = FormationFlying(t0,numPeriod,mu)
            obj.t0 = t0;
            obj.numPeriod = numPeriod;
            obj.mu = mu;
%             obj.time = makeTimeVector();
        end
    end
end