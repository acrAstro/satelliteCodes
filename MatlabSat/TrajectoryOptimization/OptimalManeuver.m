classdef OptimalManeuver < handle
    
    properties
        X
        U
        T
        samples
        t0
        dt
        tf
        umax
        B
        motionModel
        descriptor
    end
    
    methods
        function obj = OptimalManeuver(initStruct)
            obj.descriptor = initStruct.descriptor;
            obj.motionModel = eval(strcat(obj.descriptor,'(initStruct)'));
            obj.t0 = initStruct.timeParams.t0;
            obj.dt = initStruct.timeParams.dt;
            obj.tf = initStruct.timeParams.tf;
            if isempty(initStruct.maneuverParams{1})
                obj.samples = [];
            else
                obj.samples = initStruct.maneuverParams{1};
            end
            obj.B = initStruct.maneuverParams{2};
            obj.umax = initStruct.maneuverParams{3};
            obj.motionModel.makeTimeVector();
            obj.motionModel.makeDiscreteMatrices();
        end
    end
    
end