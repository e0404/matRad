function [resultGUI] = matRad_SFUDoptimization(pln, cst, dij, ct, stf)
% Calculation of single field uniform dose (SFUD) optimization
% If provided the dij matrix is used for optimisation, otherwise single
% beam dijs are calculated (memory saving).
%
% call
%   [resultGUI] = matRad_SFUDoptimization(pln, cst, dij)
%   or
%   [resultGUI] = matRad_SFUDoptimization(pln, cst, [], ct, stf)
%
%
% input
%   pln:         matRad pln struct
%   cst:         matRad cst struct
%   dij:         matRad dij struct (optional)
%                  
%   ct:          matRad ct struct (optional, only needed if no dij provided)
%   stf:         matRad stf struct (optional, only if needed no dij provided)
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%   (info:      struct containing information about optimization)
%
% References
%   [1]    https://ro-journal.biomedcentral.com/articles/10.1186/s13014-016-0705-8
%
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

% adjust cst for single beams
sb_cst = cst;
for i=1:size(sb_cst,1)
    for j = 1:size(sb_cst{i,6},1)
        % biological dose splitting for carbon
        if strcmp(pln.bioOptimization, 'LEMIV_effect') || ...
                        strcmp(pln.bioOptimization, 'LEMIV_RBExD')
            ab = sb_cst{i,5}.alphaX / sb_cst{i,5}.betaX;
            % dose per fraction
            fx_dose = sb_cst{i,6}(j).dose / pln.numOfFractions;
            % calculate dose per beam per fraction according to [1]
            fx_dose = -0.5*ab +sqrt( 0.25*ab^2 + ...
                fx_dose/pln.numOfBeams *(fx_dose + ab));
            % calculate pseudo total Dose per Beam
            sb_cst{i,6}(j).dose = fx_dose * pln.numOfFractions;
            
        % physical dose splitting
        else
            sb_cst{i,6}(j).dose = sb_cst{i,6}(j).dose/pln.numOfBeams;
        end
    end
end

if ~isempty(dij)
    % calculate dij with total dij being present
        
    % initialise total weight vector
    wTot = zeros(dij.totalNumOfBixels,1);

    for i = 1:pln.numOfBeams
        % columns in total dij for single beam
        sb_col = find(dij.beamNum == i);
        % construct dij for single beam
        sb_dij.numOfBeams = 1;
        sb_dij.doseGrid = dij.doseGrid;
        sb_dij.ctGrid = dij.ctGrid;
        sb_dij.numOfRaysPerBeam = dij.numOfRaysPerBeam(i);
        sb_dij.totalNumOfRays = sb_dij.numOfRaysPerBeam;
        sb_dij.totalNumOfBixels = size(sb_col, 1);
        sb_dij.numOfScenarios = dij.numOfScenarios;
        sb_dij.bixelNum = dij.bixelNum(sb_col);
        sb_dij.rayNum = dij.rayNum(sb_col);
        sb_dij.beamNum = dij.beamNum(sb_col);
        sb_dij.physicalDose{1} = dij.physicalDose{1}(:, sb_col);
        if isfield(dij, 'RBE')
            sb_dij.RBE = dij.RBE;
        end
        if isfield(dij, 'mLETDose')
            sb_dij.mLETDose = dij.mLETDose(:, sb_col);
        end
        if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
            sb_dij.mAlphaDose{1} = dij.mAlphaDose{1}(:, sb_col);
            sb_dij.mSqrtBetaDose{1} = dij.mSqrtBetaDose{1}(:, sb_col);
        end

        % adjust pln to one beam only
        sb_pln = pln;
        sb_pln.gantryAngles = pln.gantryAngles(i);
        sb_pln.couchAngles = pln.couchAngles(i);
        sb_pln.numOfBeams = 1;
        sb_pln.isoCenter = pln.isoCenter(i,:);
        
        % optimize single beam
        sb_resultGUI = matRad_fluenceOptimization(sb_dij,sb_cst,sb_pln);    

        % merge single beam weights into total weight vector
        wTot(sb_col) = sb_resultGUI.w;
    end

    % calculate dose
    resultGUI = matRad_calcCubes(wTot,dij,cst,1);

else
    % calculate SFUD without total dij
    
    % initialise total weight vector
    wTot = [];

    for i = 1:pln.numOfBeams
        fprintf(['optimizing beam ' num2str(i) '...\n']);
        % single beam stf
        sb_stf = stf(i);

        % adjust pln to one beam only
        sb_pln = pln;
        sb_pln.isoCenter = pln.isoCenter(i,:);
        sb_pln.numOfBeams = 1;
        sb_pln.gantryAngles = pln.gantryAngles(i);
        sb_pln.couchAngles = pln.couchAngles(i);

        % calculate single beam dij
        sb_dij = matRad_calcParticleDose(ct,sb_stf,sb_pln,sb_cst,false);

        % optimize single beam
        sb_resultGUI = matRad_fluenceOptimization(sb_dij,sb_cst,sb_pln);    

        % merge single beam weights into total weight vector
        wTot = cat(1,wTot,sb_resultGUI.w);

    end

    fprintf('Calculate total dose...\n');
    % calculate dose
    resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,wTot);
    resultGUI.w = wTot;
end
        
end %eof

