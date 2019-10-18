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
        function obj = MatRad_MCsquareBaseData(machine,stf,pln)
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
            
            %select needed energies and according focus indices            
            plannedEnergies     = [stf.ray(:).energy];
            focusIndex          = [stf.ray(:).focusIx];
            [~, ind]            = unique(plannedEnergies);
            plannedEnergies     = plannedEnergies(ind);
            focusIndex          = focusIndex(ind);
            [~ ,energyIndex, ~] = intersect([machine.data(:).energy],plannedEnergies);
            
            
            dataTable.NominalEnergy = plannedEnergies';
            rowSz = size(dataTable.NominalEnergy);
            meanEnergy       = zeros(rowSz);
            sigEnergy = zeros(rowSz);
            
            %problem with calculating sigmaEnegry?
            problemSig = false;
            
            count = 1;
            for i = energyIndex'
                
                %look up whether MonteCarlo data have already been
                %calculated, if so do not recalculate
                if isfield(machine.data(i),'mcData')
                    if (isempty(machine.data(i).mcData) == 0)
                        meanEnergy(count) = machine.data(i).mcData.MeanEnergy;
                        sigEnergy(count)  = machine.data(i).mcData.EnergySpread;
                        count = count + 1;
                        continue;
                    end
                end
                                   
                %interpolate range at 80% dose after peak.
                [maxV, maxI] = max(machine.data(i).Z);
                [~, r80ind] = min(abs(machine.data(i).Z(maxI:end) - 0.8 * maxV));
                r80ind = r80ind - 1;
                r80 = interp1(machine.data(i).Z(maxI + r80ind - 1:maxI + r80ind + 1), ...
                                 machine.data(i).depths(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV) ...
                               + machine.data(i).offset;
          
                %find FWHM w50 of bragg peak
                [~, d50rInd] = min(abs(machine.data(i).Z(maxI:end) - 0.5 * maxV));
                d50rInd = d50rInd - 1;
                d50_r = interp1(machine.data(i).Z(maxI + d50rInd - 1:maxI + d50rInd + 1), ...
                                machine.data(i).depths(maxI + d50rInd - 1:maxI + d50rInd + 1), 0.5 * maxV);

                if (machine.data(i).Z(1) < 0.4 * maxV)
                    [~, d50lInd] = min(abs(machine.data(i).Z(1:maxI) - 0.5*maxV));
                    d50_l = interp1(machine.data(i).Z(d50lInd - 1:d50lInd + 1), ...
                                    machine.data(i).depths(d50lInd - 1:d50lInd + 1), 0.5 * maxV);
                    w50 = d50_r - d50_l;

                else
                    problemSig = true;
                    w50 = (d50_r - machine.data(i).depths(maxI)) * 2;
                end
                
                %calculate energy straggling using formulae from paper "An
                %analytical approximation of the Bragg curve for 
                %therapeuticproton beams" by T. Bortfeld
                fullSigSq = (w50 / 6.14)^2;
                sigRangeStragSq = (0.012*r80)^2;
                
                if((fullSigSq - sigRangeStragSq) > 0)
                    sigEnergy(count) = sqrt(fullSigSq - sigRangeStragSq);
                else
                    sigEnergy(count) = 0.6; %set according to MCsquare documentation
                    problemSig = true;
                end
                
                meanEnergy(count) = exp(3.464048 + 0.561372013*log(r80/10) - 0.004900892*log(r80/10)^2+0.001684756748*log(r80/10)^3);
                            
                count = count + 1;
            end
            
            if problemSig
                warning('Calculation of FWHM of bragg peak in base data not possible! Using simple approximation for energy spread');
            end
            
            divergence  = zeros(rowSz);
            correlation = zeros(rowSz);
            spotsize    = zeros(rowSz);
            
            count = 1;
            for i = energyIndex'
                
                %look up whether MonteCarlo data have already been
                %calculated, if so do not recalculate
                if isfield(machine.data(i),'mcData')
                    if (isempty(machine.data(i).mcData) == 0)
                        spotsize(count)    = machine.data(i).mcData.SpotSize1x;
                        divergence(count)  = machine.data(i).mcData.Divergence1x;
                        correlation(count) = machine.data(i).mcData.Correlation1x;
                        count = count + 1;
                        continue;
                    end
                end
                
                %calculate geometric distances and extrapolate spot size at nozzle
                SAD = machine.meta.SAD;
                z     = -(machine.data(i).initFocus.dist(focusIndex(count),:) - SAD);
                sigmaSq = machine.data(i).initFocus.sigma(focusIndex(count),:).^2;

                %fit Courant-Synder equation to data using ipopt
                sigmaNull = sqrt(interp1(z,sigmaSq,0));
                             
                qRes = @(rho, sigmaT) (sigmaSq -  (sigmaNull^2 - 2*sigmaNull*rho*sigmaT.*z + sigmaT^2.*z.^2));

                funcs.objective = @(x) sum(qRes(x(1), x(2)).^2);
                funcs.gradient  = @(x) [  2 * sum(qRes(x(1), x(2)) .* (2 * sigmaNull * x(2) * z)); 
                                          2 * sum(qRes(x(1), x(2)) .* (2 * sigmaNull * x(1) * z  - 2 * x(2) * z.^2))];     
                            
                options.lb = [-0.99, -Inf];
                options.ub = [ 0.99,  Inf];
                
                options.ipopt.hessian_approximation = 'limited-memory';
                options.ipopt.limited_memory_update_type = 'bfgs';
                options.ipopt.print_level = 1;
                
                start = [0.9; 0.1];
                [result, ~] = ipopt (start, funcs, options);
                rho    = result(1);
                sigmaT = result(2);
                
                %calculate divergence, spotsize and correlation at nozzle
                DivergenceAtNozzle = sigmaT;
                SpotsizeAtNozzle = sqrt(sigmaNull^2 - 2 * rho * sigmaNull * sigmaT * obj.nozzleToIso + sigmaT^2 * obj.nozzleToIso^2);
                CorrelationAtNozzle = (rho * sigmaNull - sigmaT * obj.nozzleToIso) / SpotsizeAtNozzle;

                spotsize(count)    = SpotsizeAtNozzle;
                divergence(count)  = DivergenceAtNozzle;
                correlation(count) = CorrelationAtNozzle;
                count = count + 1;
            end
            
                dataTable.MeanEnergy = meanEnergy;
                dataTable.EnergySpread = sigEnergy;
                dataTable.ProtonsMU = ones(rowSz)*1e6;
           
                dataTable.Weight1 = ones(rowSz);
                dataTable.SpotSize1x = spotsize;
                dataTable.Divergence1x = divergence;
                dataTable.Correlation1x = correlation;
                dataTable.SpotSize1y = spotsize;
                dataTable.Divergence1y = divergence;
                dataTable.Correlation1y =  correlation;
                
                dataTable.Weight2 = zeros(rowSz);
                dataTable.SpotSize2x = zeros(rowSz);
                dataTable.Divergence2x = zeros(rowSz);
                dataTable.Correlation2x = zeros(rowSz);
                dataTable.SpotSize2y = zeros(rowSz);
                dataTable.Divergence2y = zeros(rowSz);
                dataTable.Correlation2y = zeros(rowSz);
            
                count = 1;
                newEntry = false;
                for i = energyIndex'
                    
                    %look up whether MonteCarlo data have already been
                    %calculated, if so do not recalculate
                    if isfield(machine.data(i),'mcData')
                        if (isempty(machine.data(i).mcData) == 0)
                            count = count + 1;
                            continue;
                        end
                    end
                
                    %write needed new entries in machine data
                    machine.data(i).mcData.MeanEnergy    = dataTable.MeanEnergy(count); 
                    machine.data(i).mcData.EnergySpread  = dataTable.EnergySpread(count); 
                    machine.data(i).mcData.ProtonsMU     = dataTable.ProtonsMU(count);
                    
                    machine.data(i).mcData.Weight1       = dataTable.Weight1(count);
                    machine.data(i).mcData.SpotSize1x    = dataTable.SpotSize1x(count);
                    machine.data(i).mcData.Divergence1x  = dataTable.Divergence1x(count);
                    machine.data(i).mcData.Correlation1x = dataTable.Correlation1x(count);
                    machine.data(i).mcData.SpotSize1y    = dataTable.SpotSize1y(count);
                    machine.data(i).mcData.Divergence1y  = dataTable.Divergence1y(count);
                    machine.data(i).mcData.Correlation1y = dataTable.Correlation1y(count);

                    machine.data(i).mcData.Weight2       = dataTable.Weight2(count);
                    machine.data(i).mcData.SpotSize2x    = dataTable.SpotSize2x(count);
                    machine.data(i).mcData.Divergence2x  = dataTable.Divergence2x(count);
                    machine.data(i).mcData.Correlation2x = dataTable.Correlation2x(count);
                    machine.data(i).mcData.SpotSize2y    = dataTable.SpotSize2y(count);
                    machine.data(i).mcData.Divergence2y  = dataTable.Divergence2y(count);
                    machine.data(i).mcData.Correlation2y = dataTable.Correlation2y(count);
                    
                    newEntry = true;

                    count = count +1;
                end
            
            obj.dataTable = struct2table(dataTable);
            
            %save new base data file if new entries were needed and
            %calculated
            if newEntry
                 save(strcat('../../',pln.radiationMode, '_', pln.machine, '.mat'),'machine');
                 fprintf(['New MonteCarlo data entries have been saved in', ' ', pln.radiationMode, '_', pln.machine, '.mat \n']);
            end
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

