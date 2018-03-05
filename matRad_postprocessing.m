function resultGUI = matRad_postprocessing(resultGUI, dij, pln, cst, stf)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad postprosseing function accounting for
%       minimum number of particles per spot - necessary for lmdout
%       minimum number of particles per IES - not necessary
%   
% 
% call
%   resultGUI =  matRad_postprocessing(resultGUI, dij, pln, cst, stf)

% input
%   resultGUI   struct containing optimized fluence vector
%   dij:        matRad dij struct
%   pln:        matRad pln struct
%   cst:        for calc_cubes
%   stf:        necessary for IES delete

% output
%   resultGUI:  new w and doses in resultGUI
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
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
    Imin = 500000/1e6;  %for protons
    minNrParticlesIES = 25000000;    %for protons
elseif strcmp(pln.radiationMode,'carbon')
     Imin = 15000/1e6;   
     minNrParticlesIES = 0;
else
    warning('postprocessing only for proton and carbon therapy')
end


%%manipulate fluence vector
w = resultGUI.w;
resultGUI.optW = w;
lw = length(w);
for i = 1:lw
    if(w(i) < Imin/2)
        w(i) = 0;
    elseif(w(i) > Imin/2 && w(i) < Imin)
         w(i) = Imin;   
    end
end

%%calc dose (nicht backprojection da options nicht zur Verf�gung steht,
resultGUI.w = w;

if isequal(pln.bioParam.model,'none')
    resultGUI.optDose = resultGUI.physicalDose;
else
    resultGUI.optRBExDose = resultGUI.RBExDose;
    if  strcmp(pln.bioParam.model,'constRBE') && strcmp(pln.radiationMode,'protons')
        % check if a constant RBE is defined - if not use 1.1
        if ~isfield(dij,'RBE')
            dij.RBE = 1.1;
        end
    end  
end

CalcCubes = matRad_calcCubes(w,dij,cst,1);

if isequal(pln.bioParam.model,'none')
    %%calc difference to optimized dose (not necessary, can be deleted)
    resultGUI.physicalDose = CalcCubes.physicalDose;
    relIntDoseDif = (1-sum(resultGUI.physicalDose(:))/sum(resultGUI.optDose(:)))*100;
else
    resultGUI.RBExDose = CalcCubes.RBExDose;
    relIntDoseDif = (1-sum(resultGUI.RBExDose(:))/sum(resultGUI.optRBExDose(:)))*100;
end

if relIntDoseDif ~= 0
fprintf(['Relative difference in integral dose after deleting spots: ' num2str(relIntDoseDif) '%%\n']);
end

%% delete IES with less than XXX particles
if(minNrParticlesIES ~= 0)   
    
    % Find IES values
    for beamNb = 1:pln.propStf.numOfBeams 
        iesArray = [];
        for rayNb=1:stf(beamNb).numOfRays
            iesArray = unique([iesArray stf(beamNb).ray(rayNb).energy]);
        end  
       
        for iesArrayIx=1:length(iesArray) 
            iesEnergy = iesArray(iesArrayIx);
            NrParticlesIES = 0;
            counter = 1;
    
            for rayNb=1:stf(beamNb).numOfRays 
                % find index of used energy (round to keV for numerical reasons
                bixelNb = find(round2(stf(beamNb).ray(rayNb).energy,4) == round2(iesEnergy,4)==1);
                %energyIx = find(round2(iesEnergy,4) == round2(availableEnergies,4)==1);

                if length(bixelNb)==1 % one IES found
                    bixelIndex = find([dij.beamNum==beamNb & dij.rayNum==rayNb & dij.bixelNum==bixelNb]==1);

                     voxel_nbParticles = resultGUI.w(bixelIndex);
                     voxel_nbParticles = round(1e6*voxel_nbParticles);

                    % check whether there are (enough) particles for beam delivery
                    if (voxel_nbParticles>=Imin)
                        NrParticlesIES = NrParticlesIES + voxel_nbParticles;
                        tmpBixelIndex(counter) = bixelIndex;
                        counter = counter +1;
                    end
                end
            end %ray
            if(NrParticlesIES < minNrParticlesIES && counter > 1)  % not enough particles in IES, all spots are deleted
                fprintf(['IES ' num2str(iesEnergy) ' in beam ' num2str(beamNb) ' deleted\n']);
                while(counter > 1)
                    counter = counter -1;
                    w(tmpBixelIndex(counter)) = 0;
                end              
            end 
                
        end %IES
    end %beam
   
    
resultGUI.w = w;

if isequal(pln.bioParam.model,'none') 
    resultGUI.OptDose = resultGUI.physicalDose;
else
    resultGUI.OptRBExDose = resultGUI.RBExDose;
end
CalcCubes = matRad_calcCubes(w,dij,cst,1);

if isequal(pln.bioParam.model,'none')
%%calc difference to optimized dose (not necessary, can be deleted)
resultGUI.physicalDose = CalcCubes.physicalDose;
relIntDoseDif = (1-sum(resultGUI.physicalDose(:))/sum(resultGUI.OptDose(:)))*100;

else
    resultGUI.RBExDose = CalcCubes.RBExDose;
    relIntDoseDif = (1-sum(resultGUI.RBExDose(:))/sum(resultGUI.OptRBExDose(:)))*100;

end

if relIntDoseDif ~= 0
fprintf(['Relative difference in integral dose after deleting IES: ' num2str(relIntDoseDif) '%%\n']);
end

end
                        
  
