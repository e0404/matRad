classdef matRad_SequencingIons < matRad_SequencingBase
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Constant)
        name = 'Particle IMPT Scanning Sequencing';
        shortName = 'SequencingParticle';
        possibleRadiationModes = {'protons','helium','carbon'};
    end 

    methods
        function obj = matRad_SequencingIons(inputArg1,inputArg2)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end