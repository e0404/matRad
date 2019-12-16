classdef MatRad_MCsquareBaseData
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_MCsquareBaseData Maps the matRad base data to MCsquare base data /
% phase space file
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
        mcSquareData    %MCsquare Phase space data struct
        selectedFocus   %array containing selected focus indices per energy
    end
    
    properties (SetAccess = private)
        stfCompressed   
        problemSigma
        dataTable       %Optical beam parameter table used for BDL generation
        energyIndex     %Indices of calculated energies
    end
    
    methods
        function obj = MatRad_MCsquareBaseData(machine,stf)
            %MatRad_MCsquareBaseData Construct an instance of the MCsquare
            %Base data format using a focus index
            
            %stfCompressed states whether mcSquareData are calculated for
            %all energies (false) or only for energies which exist in given
            %stf. If function is called without stf stfCompressed = false.
            if nargin < 2
                obj.stfCompressed = false;
            else
                obj.stfCompressed = true;
            end
            
            obj.machine = machine;
            obj.problemSigma = false;
            obj.selectedFocus = ones(numel(machine.data),1) * NaN;
            
            if isfield(machine.meta,'BAMStoIsoDist')
                obj.nozzleToIso = machine.meta.BAMStoIsoDist;
            else
                warning('No information on BAMS to isocenter distance. Using generic value of 500mm');
                obj.nozzleToIso = 500;
            end
            
            SAD = machine.meta.SAD;
            
            obj.smx = SAD;
            obj.smy = SAD;
            
            obj.mcSquareData = [];
            
            %select needed energies and according focus indices by using stf         
            if obj.stfCompressed
                tmp = [stf(:).ray];
                plannedEnergies     = [tmp.energy];
                focusIndex          = [tmp.focusIx];
                [~, ind]            = unique(plannedEnergies);
                plannedEnergies     = plannedEnergies(ind);
                focusIndex          = focusIndex(ind);
                [~ ,obj.energyIndex, ~] = intersect([machine.data(:).energy],plannedEnergies);
            
            %if no stf was refered all energies are chosen, while setting
            %the focus index for all energies to preliminary 1
            else
                plannedEnergies = [machine.data(:).energy];
                focusIndex = ones(size(plannedEnergies));
                [~ ,obj.energyIndex, ~] = intersect([machine.data(:).energy],plannedEnergies);
            end
            
            obj.selectedFocus(obj.energyIndex) = focusIndex;
             
            count = 1;
            for ii = 1:numel(obj.energyIndex)
                
                i = obj.energyIndex(ii);
                
                %look up whether MonteCarlo data are already present in 
                %machine file , if so do not recalculate
                if isfield(machine.data(i),'mcSquareData')
                    if (isempty(machine.data(i).mcSquareData) == 0)
                        obj.mcSquareData = [obj.mcSquareData, machine.data(i).mcSquareData];
                        count = count + 1;
                        continue;
                    end
                end
                
                
                %calculate mcSquareData for given energy and every focus
                %index
                data = [];
                energyData = obj.fitPhaseSpaceForEnergy(i);
                
                for j = 1:size(machine.data(i).initFocus.sigma,1)
                    
                    tmp = energyData;
                    opticsData = obj.fitBeamOpticsForEnergy(i, j);
                    
                    f = fieldnames(opticsData);
                    for a = 1:length(f)
                        tmp.(f{a}) = opticsData.(f{a});
                    end

                data = [data; tmp];
                end
                
                obj.mcSquareData = [obj.mcSquareData, data];
                
                count = count + 1;
            end
            
            %throw out warning if there was a problem in calculating the
            %width of the Bragg peak in obj.fitBeamOpticsForEnergy
            if obj.problemSigma
                warning('Calculation of FWHM of bragg peak in base data not possible! Using simple approximation for energy spread');
            end           
        end
        
        function mcDataEnergy = fitPhaseSpaceForEnergy(obj,energyIx)
            %function to calculate mean energy and energy spread used by 
            %mcSquare for given energy
            
            i = energyIx;
            
            mcDataEnergy.NominalEnergy = obj.machine.data(i).energy;
                       
            newDepths = linspace(0,obj.machine.data(i).depths(end),numel(obj.machine.data(i).depths) * 100);
            newDose   = interp1(obj.machine.data(i).depths, obj.machine.data(i).Z, newDepths, 'spline');       

            %find FWHM w50 of bragg peak
            [maxV, maxI] = max(newDose);
            [~, r80ind] = min(abs(newDose(maxI:end) - 0.8 * maxV));
            r80ind = r80ind - 1;
            r80 = interp1(newDose(maxI + r80ind - 1:maxI + r80ind + 1), ...
                             newDepths(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV) ...
                           + obj.machine.data(i).offset;


            [~, d50rInd] = min(abs(newDose(maxI:end) - 0.5 * maxV));
            d50rInd = d50rInd - 1;
            d50_r = interp1(newDose(maxI + d50rInd - 1:maxI + d50rInd + 1), ...
                                    newDepths(maxI + d50rInd - 1:maxI + d50rInd + 1), 0.5 * maxV);

%             if (newDose(1) < 0.4 * maxV)
                [~, d50lInd] = min(abs(newDose(1:maxI) - 0.5*maxV));
                d50_l = interp1(newDose(d50lInd - 1:d50lInd + 1), ...
                                newDepths(d50lInd - 1:d50lInd + 1), 0.5 * maxV);
%                 w50 = d50_r - d50_l;
%             %if width left of peak can be determined use twice the width to
%             %the right and throw out a warning after calculation
%             else
%                 w50 = (d50_r - newDepths(maxI)) * 2;
%                 obj.problemSigma = true;
%             end
%             fwhmBragg = w50;


            %calculate mean energy according to the mcSquare documentation
            %using the 80% dose range
%             mcDataEnergy.MeanEnergy = exp(3.464048 + 0.561372013*log(r80/10) - 0.004900892*log(r80/10)^2+0.001684756748*log(r80/10)^3); 
            mcDataEnergy.MeanEnergy = (r80 / 10 / 0.0022)^(1/1.77);

            %calculate energy straggling using formulae from paper "An
            %analytical approximation of the Bragg curve for 
            %therapeuticproton beams" by T. Bortfeld
%             fullSigSq =  2.4281 * (fwhmBragg  / 6.14)^2;
            fullSigSq1 = ((d50_r - newDepths(maxI)) / 0.62);
            fullSigSq = ((newDepths(maxI) - d50_l) / 5.52);

%             mcDataEnergy.topasFit = fullSigSq;
%             mcDataEnergy.r80 = r80;
            sigRangeStragSq = (0.012 * r80^0.935)^2; %Theoretical Formula
%             mcDataEnergy.sigRF = sigRangeStragSq;
%             sigRangeStragSq = (-3.809e-07 * r80^3  + 0.0008254  * r80^2 + ...
%                                -0.003203  * r80    +  0.1212); %TOPAS fit  
%             sigRangeStragSq = (-1.774e-08 * r80^3  +  4.969e-05  * r80^2 + ...
%                                8.524e-05  * r80    + -0.003357); %TOPAS fit
%             sigRangeStragSq = 4.941e-05 * r80^1.905;
               
% p1 =  -3.809e-07  (-8.797e-07, 1.179e-07)
%        p2 =   0.0008254  (0.0005713, 0.001079)
%        p3 =   -0.003203  (-0.04217, 0.03576)
%        p4 =      0.1212  (-1.618, 1.86)          

       
       
       

            %calculate Energy straggling using total range straggling,
            %catch error when sqrt gives imaginary results, then set energy
            %straggling to zero
            if((fullSigSq - sigRangeStragSq) >= 0)
                mcDataEnergy.EnergySpread = sqrt(fullSigSq - sigRangeStragSq);
            else
                mcDataEnergy.EnergySpread = 0;                 
                obj.problemSigma = true;
            end            
        end
            
        function mcDataOptics = fitBeamOpticsForEnergy(obj,energyIx, focusIndex)
            %function to calculate beam optics used by mcSquare for given
            %energy
            
            i = energyIx;        

            %calculate geometric distances and extrapolate spot size at nozzle
            SAD = obj.machine.meta.SAD;
            z     = -(obj.machine.data(i).initFocus.dist(focusIndex,:) - SAD);
            sigmaSq = obj.machine.data(i).initFocus.sigma(focusIndex,:).^2;

            %fit Courant-Synder equation to data using ipopt, formulae
            %given in mcSquare documentation
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
            DivergenceAtNozzle  = sigmaT;
            SpotsizeAtNozzle    = sqrt(sigmaNull^2 - 2 * rho * sigmaNull * sigmaT * obj.nozzleToIso + sigmaT^2 * obj.nozzleToIso^2);
            CorrelationAtNozzle = (rho * sigmaNull - sigmaT * obj.nozzleToIso) / SpotsizeAtNozzle;

            %save calcuated beam optics data in mcData
            mcDataOptics.ProtonsMU     = 1e6;

            mcDataOptics.Weight1       = 1;
            mcDataOptics.SpotSize1x    = SpotsizeAtNozzle;
            mcDataOptics.Divergence1x  = DivergenceAtNozzle;
            mcDataOptics.Correlation1x = CorrelationAtNozzle;
            mcDataOptics.SpotSize1y    = SpotsizeAtNozzle;
            mcDataOptics.Divergence1y  = DivergenceAtNozzle;
            mcDataOptics.Correlation1y = CorrelationAtNozzle;

            mcDataOptics.Weight2       = 0;
            mcDataOptics.SpotSize2x    = 0;
            mcDataOptics.Divergence2x  = 0;
            mcDataOptics.Correlation2x = 0;
            mcDataOptics.SpotSize2y    = 0;
            mcDataOptics.Divergence2y  = 0;
            mcDataOptics.Correlation2y = 0;
        end
                          
        function obj = writeToBDLfile(obj,filepath)
            %writeToBDLfile write the base data to file "filepath"
            
            
            %look up focus indices
            focusIndex = obj.selectedFocus(obj.energyIndex);
            
            %save mcData acording to used focus index in dataTable
            selectedData = [];
            for i = 1:numel(focusIndex)
                
                selectedData = [selectedData, obj.mcSquareData(focusIndex(i), i)];
            end
                                    
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
          
        function obj = saveMatradMachine(obj,name)
            %save previously calculated mcSquareData in new baseData file
            %with given name
            
            machine = obj.machine;
            [~ ,energyIndex, ~] = intersect([obj.machine.data(:).energy], [obj.mcSquareData(:).NominalEnergy]);
            
            machineName = [obj.machine.meta.radiationMode, '_', name];
            
            count = 1;
            for i = energyIndex'
               
                machine.data(i).mcSquareData = obj.mcSquareData(:,count);
                
                count = count + 1;
            end
            
            save(strcat('../../', machineName, '.mat'),'machine');
        end
   end
end

