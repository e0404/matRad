classdef matRad_SequencingBase
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Constant, Abstract)
        name;                       %Descriptive Name
        shortName;                  %Short name for referencing
        possibleRadiationModes;     %Possible radiation modes for the respective StfGenerator
    end

    properties (Access = public)
        radiationMode;              %Radiation Mode
        visBool = true              % vis bool
    end

    methods
        function this = matRad_SequencingBase()

        end
    end
end