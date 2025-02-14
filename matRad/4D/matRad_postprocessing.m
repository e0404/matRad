function resultGUI = matRad_postprocessing(resultGUI, dij, pln, cst, stf)
% matRad postprosseing function accounting for
%       minimum number of particles per spot
%       minimum number of particles per iso-energy slice
%   
% call
%   resultGUI =  matRad_postprocessing(resultGUI, dij, pln, cst, stf)

% input
%   resultGUI   struct containing optimized fluence vector
%   dij:        matRad dij struct
%   pln:        matRad pln struct
%   cst:        matRad cst struct
%   stf:        matRad stf struct
%
% output
%   resultGUI:  new w and doses in resultGUI
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

round2 = @(a,b)round(a*10^b)/10^b;
 
if strcmp(pln.radiationMode,'protons')
    Imin = 500000/1e6; % intensity per spot
    minNrParticlesIES = 25000000; % intensity per energy slice
elseif strcmp(pln.radiationMode,'carbon')
     Imin = 15000/1e6;   
     minNrParticlesIES = 0;
else
    matRad_cfg.dispError('postprocessing only implemented for proton and carbon ion therapy')
end

% remember old solution
resultGUI.optW = resultGUI.w;

% Get Biological Model
if ~isfield(pln,'bioModel')
    pln.bioModel = 'none';
end
if ~isa(pln.bioModel,'matRad_BiologicalModel')
    pln.bioModel = matRad_BiologicalModel.validate(pln.bioModel,pln.radiationMode);
end

% Optimization Quantity
if ~isfield(pln,'propOpt') || ~isfield(pln.propOpt, 'quantityOpt') || isempty(pln.propOpt.quantityOpt)
    pln.propOpt.quantityOpt = pln.bioModel.defaultReportQuantity;
    matRad_cfg.dispWarning('quantityOpt was not provided, using quantity suggested by biological model: %s',pln.propOpt.quantityOpt);    
end

if isequal(pln.propOpt.quantityOpt,'physicalDose')
    resultGUI.optDose = resultGUI.physicalDose;
else
    resultGUI.optRBExDose = resultGUI.RBExDose;
end

% manipulate fluence vector
% set intensities to zero if below threshold/2
resultGUI.w(resultGUI.w<Imin/2) = 0;
% set intensities to thresholds if below threshold and above threshold/2
resultGUI.w(resultGUI.w<Imin & resultGUI.w>=Imin/2) = Imin;   

% recalculate cubes!
calcCubes = matRad_calcCubes(resultGUI.w,dij,1);

% compare dose
if isequal(pln.propOpt.quantityOpt,'physicalDose')
    resultGUI.physicalDose = calcCubes.physicalDose;
    relIntDoseDif = (1-sum(resultGUI.physicalDose(:))/sum(resultGUI.optDose(:)))*100;
else
    resultGUI.RBExDose = calcCubes.RBExDose;
    relIntDoseDif = (1-sum(resultGUI.RBExDose(:))/sum(resultGUI.optRBExDose(:)))*100;
end

if relIntDoseDif ~= 0
    matRad_cfg.dispInfo('Relative difference in integral dose after deleting spots: %f %%\n',relIntDoseDif);
end

%% delete IES with less than XXX particles
if(minNrParticlesIES ~= 0)   
    
    % Find IES values
    for i = 1:pln.propStf.numOfBeams 
        iesArray = [];
        for j = 1:stf(i).numOfRays
            iesArray = unique([iesArray stf(i).ray(j).energy]);
        end  
       
        for iesArrayIx = 1:length(iesArray) 
            
            iesEnergy = iesArray(iesArrayIx);
            numParticlesIES = 0;
            bixelsIES = [];
    
            for j = 1:stf(i).numOfRays 
                % find index of used energy (round to keV for numerical reasons
                bixelNb = find(round2(stf(i).ray(j).energy,4) == round2(iesEnergy,4));
                
                if length(bixelNb)==1 % one IES found
                    
                    bixelIndex = find(dij.beamNum==i & dij.rayNum==j & dij.bixelNum==bixelNb);

                    numParticles = round(1e6*resultGUI.w(bixelIndex));
                    
                    % check whether there are (enough) particles for beam delivery
                    if (numParticles >= Imin)
                        numParticlesIES = numParticlesIES + numParticles;
                        bixelsIES = [bixelsIES bixelIndex];
                    end
                end
            end % ray
            
            if(numParticlesIES < minNrParticlesIES && ~isempty(bixelsIES))  % not enough particles in IES, all spots are deleted
                matRad_cfg.dispInfo("IES %f in beam %d deleted\n", iesEnergy, i);
                resultGUI.w(bixelsIES) = 0;
            end              
                
        end % IES
    end % beam
   
    % recalculate cubes!
    calcCubes = matRad_calcCubes(resultGUI.w,dij,1);

    % compare dose
    if isequal(pln.bioModel.model,'none')
        resultGUI.physicalDose = calcCubes.physicalDose;
        relIntDoseDif = (1-sum(resultGUI.physicalDose(:))/sum(resultGUI.optDose(:)))*100;
    else
        resultGUI.RBExDose = calcCubes.RBExDose;
        relIntDoseDif = (1-sum(resultGUI.RBExDose(:))/sum(resultGUI.optRBExDose(:)))*100;
    end

    if relIntDoseDif ~= 0
        matRad_cfg.dispInfo('Relative difference in integral dose after deleting IES: %f %%\n',relIntDoseDif);
    end

end                        
  
