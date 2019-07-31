classdef MatRad_MCsquareBaseData
%MatRad_MCsquareBaseData Maps the matRad base data to MCsquare base data /
%phase space file
%
%
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        machine         %matRad base data machine struct
        bdl_path = ''   %stores path to generated file
        nozzleToIso     %Nozzle to Isocenter Distance
        smx             %Scanning magnet X to isocenter Distance
        smy             %Scanning magnet y to isocenter Distance
        dataTable       %Optical beam parameters
    end
    
    methods
        function obj = MatRad_MCsquareBaseData(machine,focusIx)
            %MatRad_MCsquareBaseData Construct an instance of the MCsquare
            %Base data format using a focus index
            
            obj.machine = machine;
            
            if isfield(machine.meta,'BAMStoIsoDist')
                obj.nozzleToIso = machine.meta.BAMStoIsoDist;
            else
                warning('No information on BAMS to isocenter distance. Using generic value of 500mm');
                obj.nozzleToIso = 500;
            end
            
            obj.smx = machine.meta.SAD;
            obj.smy = machine.meta.SAD;
            
            dataTable.NominalEnergy = [machine.data(:).energy]';
            rowSz = size(dataTable.NominalEnergy);
            
            dataTable.MeanEnergy = dataTable.NominalEnergy; %Needs a correct value, just a dummy
            dataTable.EnergySpread = 0.8*ones(rowSz); %Was assumed during creation of the generic base data set
            
            dataTable.ProtonsMU = ones(rowSz)*1e6; %Needs a correct value
            
            %Calculate width & divergence at iso in sigma
            fwhmIso = machine.meta.LUT_bxWidthminFWHM(2,focusIx);
            sigmaIso = 0.5* fwhmIso / sqrt(2*log(2));
            
            spotDiv = sigmaIso / machine.meta.SAD;
            
            spotDiv = 1e-6*ones(rowSz);
            
            sigmaIso = arrayfun(@(s) s.sigma(1),[machine.data(:).initFocus]);
            sigmaIso = sigmaIso';
            
            if ~isfield(machine.data,'sigma')
                dataTable.Weight1 = ones(rowSz);%1 - [machine.data(:).weight]';
                dataTable.SpotSize1x = sigmaIso;
                dataTable.Divergence1x = spotDiv;
                dataTable.Correlation1x = ones(rowSz)*0;
                dataTable.SpotSize1y = sigmaIso;
                dataTable.Divergence1y = spotDiv;
                dataTable.Correlation1y = ones(rowSz)*0;
                
                dataTable.Weight2 = zeros(rowSz);%[machine.data(:).weight]';
                dataTable.SpotSize2x = sigmaIso;
                dataTable.Divergence2x = spotDiv;
                dataTable.Correlation2x = ones(rowSz)*0;
                dataTable.SpotSize2y = sigmaIso;
                dataTable.Divergence2y = spotDiv;
                dataTable.Correlation2y = ones(rowSz)*0;
            else
                dataTable.Weight1 = ones(rowSz);
                dataTable.SpotSize1x = sigmaIso;
                dataTable.Divergence1x = spotDiv;
                dataTable.Correlation1x = ones(rowSz)*0;
                dataTable.SpotSize1y = sigmaIso;
                dataTable.Divergence1y = spotDiv;
                dataTable.Correlation1y = ones(rowSz)*0;
                
                dataTable.Weight2 = ones(rowSz)*0;
                dataTable.SpotSize2x = sigmaIso;
                dataTable.Divergence2x = spotDiv;
                dataTable.Correlation2x = ones(rowSz)*0;
                dataTable.SpotSize2y = sigmaIso;
                dataTable.Divergence2y = spotDiv;
                dataTable.Correlation2y = ones(rowSz)*0;
            end
            
            obj.dataTable = struct2table(dataTable);
        end
        
        
        function [obj] = writeToBDLfile(obj,filepath)
            %writeToBDLfile write the base data to file "filepath"
            
            machine = obj.machine;
            
            try
                
                fileID = fopen(filepath,'w');
                
                %Header
                %fprintf(fileID,'--matRad: Beam Model for machine %s (%s)--\n',machine.meta.machine,machine.meta.dataType);
                fprintf(fileID,'--UPenn beam model (double gaussian)--\n');
                fprintf(fileID,'# %s\n',machine.meta.description);
                fprintf(fileID,'# created by %s on %s\n\n',machine.meta.created_by,machine.meta.created_on);
                
                fprintf(fileID,'Nozzle exit to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n',obj.nozzleToIso);
                
                fprintf(fileID,'SMX to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n',obj.smx);
                
                fprintf(fileID,'SMY to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n',obj.smy);
                
                fprintf(fileID,'Beam parameters\n%d energies\n\n',size(obj.dataTable,1));

                
                fclose(fileID);               
                
                        
                writetable(obj.dataTable,[filepath '.tmp'],'Delimiter','\t','FileType','text');
                
                fileID = fopen([filepath '.tmp'],'r');
                tableTxt = fread(fileID);
                fclose(fileID);
                
                fileID = fopen(filepath,'a');
                fwrite(fileID,tableTxt);
                fclose(fileID);
                
                delete([filepath '.tmp']);
                
                obj.bdl_path = filepath;
                
            catch MException
                error(MException.message);
            end
        end
    end
end

