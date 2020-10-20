classdef MatRad_MCsquareBaseData < MatRad_MCemittanceBaseData
    % MatRad_MCsquareBaseData class for calculating MCsquare base data and
    % writing it to a .txt file, for MCsquare to use
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = public)
        function obj = MatRad_MCsquareBaseData(machine,stf)
            %Call MatRad_MCemmitanceBaseData constructor
            obj = obj@MatRad_MCemittanceBaseData(machine,stf);             
        end
                
        function obj = writeMCsquareData(obj,filepath)
            %function that writes a data file containing Monte Carlo base
            %data for a simulation with MCsquare
            
            %look up focus indices
            focusIndex = obj.selectedFocus(obj.energyIndex);
            
            %save mcData acording to used focus index in selectedData
            selectedData = [];
            for i = 1:numel(focusIndex)
                selectedData = [selectedData, obj.monteCarloData(focusIndex(i), i)];
            end
            
            %remove field not needed for MCsquare base data
            selectedData = rmfield(selectedData, 'FWHMatIso');
            
            %write MCsqaure data base file
            try
                
                fileID = fopen(filepath,'w');
                
                %Header
                %fprintf(fileID,'--matRad: Beam Model for machine %s (%s)--\n',machine.meta.machine,machine.meta.dataType);
                fprintf(fileID,'--UPenn beam model (double gaussian)--\n');
                fprintf(fileID,'# %s\n', obj.machine.meta.description);
                fprintf(fileID,'# created by %s on %s\n\n', obj.machine.meta.created_by, obj.machine.meta.created_on);
                
                fprintf(fileID,'Nozzle exit to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n',obj.nozzleToIso);
                
                fprintf(fileID,'SMX to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n',obj.smx);
                
                fprintf(fileID,'SMY to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n',obj.smy);
                
                for i = 1:length(obj.rangeShifters)
                    raShi = obj.rangeShifters(i);
                    fprintf(fileID,'Range Shifter parameters\n');
                    fprintf(fileID,'RS_ID = %d\n',raShi.ID);
                    fprintf(fileID,'RS_type = binary\n');
                    fprintf(fileID,'RS_material = 64\n'); %Material ID Hardcoded for now (PMMA)
                    fprintf(fileID,'RS_density = 1.19\n'); %Maetiral density Hardcoded for now (PMMA)
                    fprintf(fileID,'RS_WET = %f\n\n',raShi.eqThickness);
                end
                    
                
                fprintf(fileID,'Beam parameters\n%d energies\n\n',size(selectedData,2));
                
                fn = fieldnames(selectedData);
                for names = 1:size(fn,1)
                    fprintf(fileID, fn{names});
                    fprintf(fileID, '\t');
                end
                fprintf(fileID, '\n');
                
                for k = 1:size(selectedData,2)
                    for m = 1:numel(fn)
                        fprintf(fileID, '%g', selectedData(k).(fn{m}));
                        fprintf(fileID, '\t');
                    end
                    fprintf(fileID, '\n');
                end
                
                fclose(fileID);
                
                obj.bdl_path = filepath;
                
            catch MException
                error(MException.message);
            end
        end        
    end    
end

