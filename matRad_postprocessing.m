function resultGUI = matRad_postprocessing(resultGUI, dij, pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad postprosseing function accounting for
%       minimum number of particles per spot
%       minimum number of particles per IES???
%       (scan path)???
%       ...
% 
% call
%   resultGUI = matRad_fluenceOptimization(resultGUI, dij,cst,pln)
%
% input
%   resultGUI   struct containing optimized fluence vector
%   dij:        matRad dij struct
%   pln:        matRad pln struct
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%
% References
%   -
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


% issue warning if biological optimization impossible
if sum(strcmp(pln.bioOptimization,{'effect','RBExD'}))>0 && (~isfield(dij,'mAlphaDose') || ~isfield(dij,'mSqrtBetaDose'))
    warndlg('Alpha and beta matrices for effect based and RBE optimization not available - physical optimization is carried out instead.');
    pln.bioOptimization = 'none';
end

if ~isdeployed % only if _not_ running as standalone
    
    % add path for optimization functions
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(matRadRootDir,'optimization'))
    
    % get handle to Matlab command window
    mde         = com.mathworks.mde.desk.MLDesktop.getInstance;
    cw          = mde.getClient('Command Window');
    xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);
    h_cw        = handle(xCmdWndView,'CallbackProperties');

    % set Key Pressed Callback of Matlab command window
    set(h_cw, 'KeyPressedCallback', @matRad_CWKeyPressedCallback);

end

if ~strcmp(pln.radiationMode,'protons') && ~strcmp(pln.radiationMode,'carbon')
  error('HITXML plan for this radiationMode not supported!');
end
Imin = 500000/1e6;  %for protons
minNrParticlesIES = 25000000;    %for protons
if strcmp(pln.radiationMode,'carbon')
     Imin = 15000/1e6;   
     minNrParticlesIES = 0;
end



%%manipulate fluence vector
w = resultGUI.w;
%Imin = pln.minNrParticles/1e6;
lw = length(w);
for i = 1:lw
    if(w(i) < Imin/2)
        w(i) = 0;
    elseif(w(i) > Imin/2 && w(i) < Imin)
        %w(i) = Imin;    WIEDER ÄNDERN; SO KONSISTENT MIT LMDOUT
         w(i) = 0;
    end
end

%%calc dose (nicht backprojection da options nicht zur Verfügung steht,
%%nicht calc cubes da sonst resultGUI.physicalDose überschrieben wird
resultGUI.finalDose = reshape(dij.physicalDose{1}*w,dij.dimensions);
%d = matRad_backProjection(w,dij,'none');

%resultGUI.finalDose = reshape(d{1},dij.dimensions);
resultGUI.finalw = w;

%%calc difference to optimized dose (not necessary, can be deleted)
relIntDoseDif = (1-sum(resultGUI.physicalDose(:))/sum(resultGUI.finalDose(:)))*100;

fprintf(['Relative difference in integral dose after deleting spots: ' num2str(relIntDoseDif) '%%\n']);


%% delete IES with less than XXX particles
if(minNrParticlesIES ~= 0)
    
    % get data from workspace
    stf       = evalin('base','stf');
    % helper function for energy selection
    round2 = @(a,b)round(a*10^b)/10^b;
    
    % Find IES values
    for beamNb = 1:pln.numOfBeams 
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

                     voxel_nbParticles = resultGUI.finalw(bixelIndex);
                     voxel_nbParticles = round(1e6*voxel_nbParticles);

                    % check whether there are (enough) particles for beam delivery
                    if (voxel_nbParticles>pln.minNrParticles)
                        NrParticlesIES = NrParticlesIES + voxel_nbParticles;
                        tmpBixelIndex(counter) = bixelIndex;
                        counter = counter +1;
                    end
                end
            end %ray
            if(NrParticlesIES < minNrParticlesIES && counter > 1)  % not enough particles in IES, all spots are deleted
                fprintf(['IES ' num2str(iesEnergy) ' deleted\n']);
                while(counter > 1)
                    counter = counter -1;
                    w(tmpBixelIndex(counter)) = 0;
                end              
            end 
                
        end %IES
    end %beam
   
    
%%calc dose
%d = matRad_backProjection(w,dij,'none');
resultGUI.finalDose = reshape(dij.physicalDose{1}*w,dij.dimensions);
%resultGUI.finalDose = reshape(d{1},dij.dimensions);
resultGUI.finalw = w;

%%calc difference to optimized dose (not necessary, can be deleted)
relIntDoseDif = (1-sum(resultGUI.physicalDose(:))/sum(resultGUI.finalDose(:)))*100;

fprintf(['Relative difference in integral dose after deleting IES: ' num2str(relIntDoseDif) '%%\n']);

maxDosefinal = max(resultGUI.finalDose(:));
maxDiff = max(resultGUI.physicalDose(:) - resultGUI.finalDose(:));
    
fprintf(['Absolute difference: ' num2str(maxDiff) 'Gy. Max dose: ' num2str(maxDosefinal) 'Gy\n']);
end
                        
  
