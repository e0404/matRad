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
        function obj = MatRad_MCsquareBaseData(machine,stf)
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
            r80       = zeros(rowSz);
            sigEnergy = zeros(rowSz);
            
            %problem with calculating sigmaEnegry?
            problemSig = false;
            
            count = 1;
            for i = energyIndex'
                %interpolate range at 80% dose after peak.
                [maxV, maxI] = max(machine.data(i).Z);
                [~, r80ind] = min(abs(machine.data(i).Z(maxI:end) - 0.8 * maxV));
                r80ind = r80ind - 1;
                r80(count) = interp1(machine.data(i).Z(maxI + r80ind - 1:maxI + r80ind + 1), ...
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
                sigRangeStragSq = (0.012*r80(count))^2;
                
                if((fullSigSq - sigRangeStragSq) > 0)
                    sigEnergy(count) = sqrt(fullSigSq - sigRangeStragSq);
                else
                    sigEnergy(count) = 0.6; %set according to MCsquare documentation
                    problemSig = true;
                end
                count = count + 1;
            end
            
            if problemSig
                warning('Calculation of FWHM of bragg peak in base data not possible! Using simple approximation for energy spread');
            end
            
            dataTable.MeanEnergy = arrayfun(@(r) exp(3.464048 + 0.561372013*log(r/10) - 0.004900892*log(r/10)^2+0.001684756748*log(r/10)^3),r80')'; %Energy calculated with formula in MCSquare documentation
            %dataTable.MeanEnergy = dataTable.MeanEnergy';
            dataTable.EnergySpread = sigEnergy; %Energy spread as described in mcSquare documentation
            
            dataTable.ProtonsMU = ones(rowSz)*1e6; %Doesn't seem to be implemented in MCsquare despite in BDL file?
                       

            dataBDL.spotSize = zeros(rowSz);
            dataBDL.divergence = zeros(rowSz);
            dataBDL.correlation = zeros(rowSz);

            
            count = 1;
            for i = energyIndex'
                
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


                dataBDL.spotsize(count)    = SpotsizeAtNozzle;
                dataBDL.divergence(count)  = DivergenceAtNozzle;
                dataBDL.correlation(count) = CorrelationAtNozzle;
                count = count + 1;
            end
           
                dataTable.Weight1 = ones(rowSz);
                dataTable.SpotSize1x = dataBDL.spotsize';
                dataTable.Divergence1x = dataBDL.divergence;
                dataTable.Correlation1x = dataBDL.correlation;
                dataTable.SpotSize1y = dataBDL.spotsize';
                dataTable.Divergence1y = dataBDL.divergence;
                dataTable.Correlation1y =  dataBDL.correlation;
                
                dataTable.Weight2 = zeros(rowSz);
                dataTable.SpotSize2x = zeros(rowSz);
                dataTable.Divergence2x = zeros(rowSz);
                dataTable.Correlation2x = zeros(rowSz);
                dataTable.SpotSize2y = zeros(rowSz);
                dataTable.Divergence2y = zeros(rowSz);
                dataTable.Correlation2y = zeros(rowSz);
            
            
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

