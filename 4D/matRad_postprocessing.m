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
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

round2 = @(a,b)round(a*10^b)/10^b;
 
if strcmp(pln.radiationMode,'protons')
    Imin = 500000/1e6; % intensity per spot
    minNrParticlesIES = 25000000; % intensity per energy slice
elseif strcmp(pln.radiationMode,'carbon')
     Imin = 15000/1e6;   
     minNrParticlesIES = 0;
else
    error('postprocessing only implemented for proton and carbon ion therapy')
end

% remember old solution
resultGUI.optW = resultGUI.w;
if isequal(pln.bioParam.model,'none')
    resultGUI.optDose = resultGUI.physicalDose;
else
    resultGUI.optRBExD = resultGUI.RBExD;
end

% manipulate fluence vector
% set intensities to zero if below threshold/2
resultGUI.w(resultGUI.w<Imin/2) = 0;
% set intensities to thresholds if below threshold and above threshold/2
resultGUI.w(resultGUI.w<Imin & resultGUI.w>=Imin/2) = Imin;   

% recalculate cubes!
calcCubes = matRad_calcCubes(resultGUI.w,dij,1);

% compare dose
if isequal(pln.bioParam.model,'none')
    resultGUI.physicalDose = calcCubes.physicalDose;
    relIntDoseDif = (1-sum(resultGUI.physicalDose(:))/sum(resultGUI.optDose(:)))*100;
else
    resultGUI.RBExD = calcCubes.RBExD;
    relIntDoseDif = (1-sum(resultGUI.RBExD(:))/sum(resultGUI.optRBExD(:)))*100;
end

if relIntDoseDif ~= 0
    fprintf(['Relative difference in integral dose after deleting spots: ' num2str(relIntDoseDif) '%%\n']);
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
                fprintf(['IES ' num2str(iesEnergy) ' in beam ' num2str(i) ' deleted\n']);
                resultGUI.w(bixelsIES) = 0;
            end              
                
        end % IES
    end % beam
   
    % recalculate cubes!
    calcCubes = matRad_calcCubes(resultGUI.w,dij,1);

    % compare dose
    if isequal(pln.bioParam.model,'none')
        resultGUI.physicalDose = calcCubes.physicalDose;
        relIntDoseDif = (1-sum(resultGUI.physicalDose(:))/sum(resultGUI.optDose(:)))*100;
    else
        resultGUI.RBExD = calcCubes.RBExD;
        relIntDoseDif = (1-sum(resultGUI.RBExD(:))/sum(resultGUI.optRBExD(:)))*100;
    end

    if relIntDoseDif ~= 0
        fprintf(['Relative difference in integral dose after deleting IES: ' num2str(relIntDoseDif) '%%\n']);
    end

end                        
  
