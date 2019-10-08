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
            
            SAD = machine.meta.SAD;
            
            obj.smx = SAD;
            obj.smy = SAD;
            
            dataTable.NominalEnergy = [machine.data(:).energy]';
            rowSz = size(dataTable.NominalEnergy);
            
            for i = 1:numel(machine.data)
                %find position of bragg peak
                [maxV, maxIndex] = max(machine.data(i).Z);
                
                %find position of 80% dose fall off position
                [~,ix] = min(abs(machine.data(i).Z(maxIndex:end) - 0.8*maxV));
                
                %interpolate exact position of 80% dose fall off
                r80(i) = interp1(machine.data(i).Z(maxIndex+ix-5:maxIndex+ix+5),machine.data(i).depths(maxIndex+ix-5:maxIndex+ix+5),0.8*maxV);
                r80(i) = r80(i) + machine.data(i).offset;
                
                
                %find and interpolate exact positions of 50% dose range
                [~,ix] = min(abs(machine.data(i).Z(maxIndex:end) - 0.5*maxV));
                d50_r = interp1(machine.data(i).Z(maxIndex+ix-5:maxIndex+ix+5),machine.data(i).depths(maxIndex+ix-5:maxIndex+ix+5),0.5*maxV);
                
                [~,ix] = min(abs(machine.data(i).Z(1:maxIndex) - 0.5*maxV));
                d50_l = interp1(machine.data(i).Z(ix-5:ix+5),machine.data(i).depths(ix-5:ix+5),0.5*maxV);
                
                %calculate FWHM of bragg peak
                w50 = abs(d50_r - d50_l);
                
                %calculate energy straggling according to "The physics of
                %chared paticle therapy, M. Bangert
                fullSigSq = (w50 / 6.14)^2;
                sigRangeStragSq = (0.012*r80(i)^0.95)^2;
                
                sigEnergy(i) = sqrt(fullSigSq - sigRangeStragSq);

            end

            
            
            %dataTable.MeanEnergy = arrayfun(@(r) exp(3.464048 + 0.561372013*log(r/10) - 0.004900892*log(r/10)^2+0.001684756748*log(r/10)^3),r80'); %Energy calculated with formula in MCSquare documentation
            dataTable.MeanEnergy =  (r80./0.022).^(1/1.77)';
            dataTable.EnergySpread = sigEnergy'; %energy and energy spread as described in mcSquare documentation
            
            dataTable.ProtonsMU = ones(rowSz)*1e6; %Doesn't seem to be implemented in MCsquare despite in BDL file?
            
            %fwhmIso = machine.meta.LUT_bxWidthminFWHM(2,focusIx);
            
            %interpolate sigma^2 at isocenter
            sigma0 = arrayfun(@(s) interp1(s.dist(focusIx,:),s.sigma(focusIx,:),SAD),[machine.data(:).initFocus]);
            
            
            %use curve fitting for finding spot sizes
            for i = 1:numel(machine.data)
                sigmaxSq = machine.data(i).initFocus.sigma(focusIx,:).^2;
                z = machine.data(i).initFocus.dist(focusIx,:) - SAD;
                func = @(a, x) sigma0(i)^2 - 2 *a*sigma0(i) + a^2*x.^2;
                ft = fittype(func);
                myFit = fit(z', sigmaxSq', ft);
                y = machine.data(i).initFocus.sigma(focusIx,:).^2;%- sigma0(i).^2;
                spotDiv(i) = sqrt(myFit.a);
                sigmaIso(i) = sqrt(sigma0(i)^2 + myFit.a * obj.nozzleToIso^2);
                sigmaCorr(i) = (sigma0(i) - spotDiv(i) *obj.nozzleToIso)/ func(spotDiv(i), obj.nozzleToIso);
            end

            
            sigmaIso = sigmaIso';
            spotDiv(spotDiv == 0) = 1e-20;
            spotDiv = spotDiv';
            sigmaCorr = sigmaCorr';
            
            if ~isfield(machine.data,'sigma')
                dataTable.Weight1 = ones(rowSz);
                dataTable.SpotSize1x = sigmaIso;
                dataTable.Divergence1x = spotDiv;
                dataTable.Correlation1x = sigmaCorr;
                dataTable.SpotSize1y = sigmaIso;
                dataTable.Divergence1y = spotDiv;
                dataTable.Correlation1y = sigmaCorr;
                
                dataTable.Weight2 = zeros(rowSz);
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

