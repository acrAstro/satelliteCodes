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
        model
        motion
    end
    
    methods
        function obj = OptimalManeuver(params,initStruct)
            obj.model = params.Model;
            switch obj.model
                case 'GA'
                    obj.motion = GimAlfriendSTM(initStruct)
                case 'HCW'
                case 'LERM'
                case 'SS'
            end
            obj.samples = params.timeParams(1);
            obj.t0 = params.timeParams(2);
            obj.dt = params.timeParams(3);
            obj.tf = params.timeParams(4);
            obj.umax = params.inputParams(1);
            obj.B = params.inputParams(2);
            
        end
    end
    
end