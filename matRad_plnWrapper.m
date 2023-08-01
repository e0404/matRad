function plnJO = matRad_plnWrapper(pln)
% Combines the arbitrary number of input plans into a single plan for the
% different modalities.
% 
% call
%   plnJO = matRad_plnWrapper(pln)
%
% input
%   pln:       array of pln structure for the different modalities (if any)
%
% output
%   plnJO:      synthetic overarching pln stuct for Joint Opt 
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
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
nPlans = length(pln);
if nPlans>0
    
    allFields = arrayfun(@(modality) fieldnames(pln(modality)), [1:nPlans], 'UniformOutput',false);
    allFields = unique(cat(1,allFields{:}));
    plnJO = cell2struct(cell(1,length(allFields)), allFields, 2);
   
    %First run plnConsistency, it is the same as current MixModality, small
    %corrections are introduced to handle output generation with matRad_cfg
    pln = matRad_plnConsistency(pln); %to be reviewed
   
    % Save now the plns. This will have now the updated pln.bioParam
    originalPlans = pln;
   
    %Define the Pln properties
    plnJO.numOfFractions  = sum([pln(:).numOfFractions]); %Use total number of fractions
    plnJO.radiationMode   = 'MixMod';
    plnJO.machine         = 'MixMod';
    plnJO.propStf         = [pln(:).propStf];
    plnJO.numOfModalities = nPlans;

    %%%%%%% propDoseCalc %%%%%
    allFields = arrayfun(@(modality) fieldnames(pln(modality).propDoseCalc), [1:plnJO.numOfModalities], 'UniformOutput',false);
    allFields = unique(cat(1,allFields{:}));
    plnJO.propDoseCalc(1:plnJO.numOfModalities) = cell2struct(cell(1,length(allFields)), allFields, 2);

     for modalityIdx=1:plnJO.numOfModalities
        currModalityFields = fieldnames(pln(modalityIdx).propDoseCalc);
        for fieldIdx=1:length(currModalityFields)
           plnJO.propDoseCalc(modalityIdx).(currModalityFields{fieldIdx}) = pln(modalityIdx).propDoseCalc.(currModalityFields{fieldIdx});
        end
    end

    %%%%%%% propStf %%%%%
    allFields = arrayfun(@(modality) fieldnames(pln(modality).propStf), [1:plnJO.numOfModalities], 'UniformOutput',false);
    allFields = unique(cat(1,allFields{:}));
    plnJO.propStf(1:plnJO.numOfModalities) = cell2struct(cell(1,length(allFields)), allFields, 2);

    for modalityIdx=1:plnJO.numOfModalities
        currModalityFields = fieldnames(pln(modalityIdx).propStf);
        for fieldIdx=1:length(currModalityFields)
           plnJO.propStf(modalityIdx).(currModalityFields{fieldIdx}) = pln(modalityIdx).propStf.(currModalityFields{fieldIdx});
        end
    end

    %%%%%%% propOpt %%%%%
    allFields = arrayfun(@(modality) fieldnames(pln(modality).propOpt), [1:plnJO.numOfModalities], 'UniformOutput',false);
    allFields = unique(cat(1,allFields{:}));
    plnJO.propOpt = cell2struct(cell(1,length(allFields)), allFields, 2);
    
    
    for modalityIdx=1:plnJO.numOfModalities

        currModalityFields = fieldnames(pln(modalityIdx).propOpt);
        for fieldIdx=1:length(currModalityFields)
            if ~any(strcmp(currModalityFields{fieldIdx}, {'spatioTemp', 'STfractions', 'STscenarios'}))
                plnJO.propOpt.(currModalityFields{fieldIdx}) = pln(modalityIdx).propOpt.(currModalityFields{fieldIdx});
            end
        end
    end

    if ~isfield(plnJO.propOpt, {'spatioTemp'})
        plnJO.propOpt.spatioTemp = zeros(1,plnJO.numOfModalities);
        plnJO.propOpt.STfractions = num2cell(pln(:).numOfFractions);
        plnJO.propOpt.STscenarios = ones(1,plnJO.numOfModalities);
    end

    for modalityIdx=1:plnJO.numOfModalities
        if pln(modalityIdx).propOpt.spatioTemp
            plnJO.propOpt.spatioTemp(modalityIdx) = 1;

            if isfield(pln(modalityIdx).propOpt, 'STscenarios') && pln(modalityIdx).propOpt.STscenarios>1
                
                if isfield(pln(modalityIdx).propOpt, 'STfractions') && sum(pln(modalityIdx).propOpt.STfractions) == pln(modalityIdx).numOfFractions
                    
                    plnJO.propOpt.STscenarios(modalityIdx) = pln(modalityIdx).propOpt.STscenarios;
                    
                    if length(pln(modalityIdx).propOpt.STfractions) == pln(modalityIdx).propOpt.STscenarios
                        plnJO.propOpt.STfractions(modalityIdx) = {pln(modalityIdx).propOpt.STfractions};
                        matRad_cfg.dispInfo(['STfractionation for modality ', num2str(modalityIdx), ' correctly set with ', num2str(plnJO.propOpt.STscenarios(modalityIdx)), ' scenarios and [', num2str(plnJO.propOpt.STfractions{modalityIdx}), '] fractions\n']);
                    else
                        matRad_cfg.dispError(['Number of fractions for modality ', num2str(modalityIdx), ' needs to be set for every scenario']);

                    end
                else
                    matRad_cfg.dispError(['Fractionation scheme of modality ',num2str(modalityIdx), ' must match the total number of fractions']);
                end
            else
                plnJO.propOpt.spatioTemp(modalityIdx) = 0;
                plnJO.propOpt.STscenarios(modalityIdx) = 1;
                plnJO.propOpt.STfractions(modalityIdx) = {pln(modalityIdx).numOfFractions};
                matRad_cfg.dispWarning(['STfractionation required for modality ',num2str(modalityIdx), ' but no fractionation scheme specified. Disabling STfractionation.\n']);
            end
        else
            plnJO.propOpt.spatioTemp(modalityIdx) = 0;
            plnJO.propOpt.STscenarios = [plnJO.propOpt.STscenarios, 1];
            plnJO.propOpt.STfractions = [plnJO.propOpt.STfractions, {pln(modalityIdx).numOfFractions}];
        end
    end
    
    % Feed the first bio model quantity, they are consistent
    plnJO.bioParam = matRad_bioModel(plnJO.radiationMode,pln(1).bioParam.quantityOpt, plnJO.radiationMode);
    plnJO.originalPlans = originalPlans;
    plnJO.multScen = [pln(:).multScen];

       
   
else
   %Do nothing
   plnJO = pln;
end


end