classdef MatRad_MCsquareBaseData
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        machine
        bdl_path = ''
        nozzleToIso
        smx
        smy
        dataTable
    end
    
    methods
        function obj = MatRad_MCsquareBaseData(machine,focusIx)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
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
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
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

