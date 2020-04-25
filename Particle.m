classdef Particle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Weight
        X
        Y
        Lm
        LmP
        T_x
        T_y
    end
    
    methods
        function obj = Particle(N_PARTICLES,N_LM,LM_SIZE)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.Weight = 1/N_PARTICLES;
            obj.X = 0;
            obj.Y = 0;
            obj.Lm = zeros(N_LM, LM_SIZE);
            obj.LmP = [];
            for i=1:N_LM
                obj.LmP = [obj.LmP,zeros(LM_SIZE, LM_SIZE)];
            end
            obj.T_x = [];
            obj.T_y = [];
            
            
        end
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

