function pln = matRad_plnConsistency(pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_plnConsistency checks the consistency of the optimized quantities 
% for the multiple modalities; all optimization quantities should be the
% same!
%
% call
%   pln = matRad_plnConsistency(pln, param)
%
% input
%   pln(:):          pln stuct for multiple modalities  
% output
%   pln struct
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();
AvailableQuantitiesForOpt = {'physicalDose', 'effect', 'RBExD'};
if numel(pln) == 1
    matRad_cfg.dispInfo(['single pln is by definition consistent with itself... \n']);
elseif numel(pln) > 1
    maxh = [];
    % find level in the heirarchy (physical, RBExD, Effect)
    for i = 1:numel(pln)
        switch pln(i).bioParam.quantityOpt
            case AvailableQuantitiesForOpt(1)   % physicalDose
                maxh = [maxh 1];
            case AvailableQuantitiesForOpt(2)   % effect
                maxh = [maxh 2];
            case AvailableQuantitiesForOpt(3)   % RBExD
                maxh = [maxh 3];
        end
    end
    

    % for plnIdx = 1:numel(pln)
    % 
    %         % Avoid case in which varRBE for protons and constRBE for
    %         % photons (set model none for photons, so avoids the instance of constRBE backProjection later on. Make sure that only case in which there is constRBE_BP is when protons have constRBE))
    % 
    %         % In any case if RBExD with photons, the other modalities get a
    %         % varRBE prjection. RBExD wins over all, so it is set as a
    %         % quantity for optimization. And then you will need a model for
    %         % protons and it will be MCN (or constRBE)
    %         bioModelExeptThisPln = [pln.bioParam];
    %         bioModelExeptThisPln(plnIdx) = [];
    % 
    %          if strcmp(pln(plnIdx).radiationMode, 'photons') && (strcmp(pln(plnIdx).bioParam.model, 'constRBE')) && strcmp(pln(plnIdx).bioParam.quantityOpt, 'RBExD') %if ohoton plan and constRBE model
    %             if ~all(strcmp({bioModelExeptThisPln.model}, 'constRBE'))
    %                 pln(plnIdx).bioParam = matRad_bioModel(pln(plnIdx).radiationMode, 'RBExD', 'none'); % switch off the constRBE model for photons is it is the only one and we have RBExD. Means the other models have varRBE (?) 
    %             end
    %          end
    % end

    if numel(unique(maxh)) > 1
        maxh = max(maxh);
        %This will probably clash with new class based bioModel implementation
        quantityOpt = AvailableQuantitiesForOpt{maxh};

        matRad_cfg.dispInfo(['make plan consistent. optimized quantity chosen by hirarchy is ' quantityOpt '. \n' ]);

        for i = 1:numel(pln)
            pln(i).bioParam = matRad_bioModel(pln(i).radiationMode, quantityOpt, pln(i).bioParam.model);
        end
    elseif numel(unique(maxh)) == 1

        % case in which both plans are
        matRad_cfg.dispInfo(['pln is consistent... \n']);
    end
end
