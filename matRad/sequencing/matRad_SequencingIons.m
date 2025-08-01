classdef matRad_SequencingIons < matRad_SequencingBase
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Constant)
        name = 'Particle IMPT Scanning Sequencing';
        shortName = 'SequencingParticle';
        possibleRadiationModes = {'protons','helium','carbon'};
    end 

    methods  (Static)
        function [available,msg] = isAvailable(pln,machine)
                % see superclass for information            
                       
                if nargin < 2
                    machine = matRad_loadMachine(pln);
                end
                %checkBasic
                available = isfield(machine,'meta') && isfield(machine,'data');
    
                available = available && any(isfield(machine.meta,{'machine','radiationMode'}));
        
                if ~available
                    msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';                
                else
                    msg = [];
                end
   
                 %check modality
                 checkModality = any(strcmp(matRad_SequencingIons.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_SequencingIons.possibleRadiationModes, pln.radiationMode));
                    
                 %Sanity check compatibility
                 if checkModality
                    checkModality = strcmp(machine.meta.radiationMode,pln.radiationMode);
                 end
        
                 available  = available && checkModality;
       
        end
    end
end