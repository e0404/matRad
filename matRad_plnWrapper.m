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

    %%% Disable STfractionation fo the time being
    for modalityIdx=1:plnJO.numOfModalities
        if (isfield(pln(modalityIdx).propOpt, 'spatioTemp')) && pln(modalityIdx).propOpt.spatioTemp
            matRad.dispWarning('Sorry no spatioTemporal avalability yet');
        end

    end
    plnJO.propOpt.spatioTemp = zeros(1,plnJO.numOfModalities);
    plnJO.propOpt.STfractions = {pln.numOfFractions};
    plnJO.propOpt.STscenarios = ones(1,plnJO.numOfModalities);

    % Feed the first bio model quantity, they are consistent
    plnJO.bioParam = matRad_bioModel(plnJO.radiationMode,pln(1).bioParam.quantityOpt, plnJO.radiationMode);
    plnJO.originalPlans = originalPlans;
    plnJO.multScen = [pln(:).multScen];


   % for k=1:length(currentFields)
   %    if isempty(getfield(plnJO,currentFields{1,k})) && isstruct(pln(1).(currentFields{1,k}))
   %      %For all pln fields that are structures, check that the number of
   %      %fields is the same
   %      plnfld = matRad_fieldConsistency(pln,currentFields(1,k));
   %      if strcmp(currentFields{1,k},'propDoseCalc')
   %          plnJO.(currentFields{1,k}) = [plnfld(:,1)];                         % check for consistency and singletonnness of propDoseCalc and PropOpt
   %      elseif strcmp(currentFields{1,k},'propOpt')
   %          plnJO.(currentFields{1,k}) = [plnfld(:,1)];   
   % 
   %      else
   %          plnJO.(currentFields{1,k}) = [plnfld(:)];                         % check for consistency and singletonnness of propDoseCalc and PropOpt
   %      end
   % 
   % 
   %    elseif isempty(getfield(plnJO,currentFields{1,k})) && ~isstruct(pln(1).(currentFields{1,k}))
   %       %If the field is a class, just keep it
   %       plnJO.(currentFields{1,k}) = [pln(:).(currentFields{1,k})];
   %    end
   % end
   % plnJO.bioParam = matRad_bioModel(plnJO.radiationMode,pln(1).bioParam.quantityOpt, plnJO.radiationMode);
   % %Save the original plans as well
   % plnJO.originalPlans = originalPlans;
   % plnJO.numOfModalities = nPlans;

   % SORT OUT THE ST Fields


   % if isfield(plnJO.propOpt,'spatioTemp')
   % 
   %     for i = 1: plnJO.numOfModalities
   %          plnJO.propOpt.spatioTemp(i)  =    pln(i).propOpt.spatioTemp ;
   %     end 
   %     if isfield(plnJO.propOpt,'STscenarios') && sum(plnJO.propOpt.spatioTemp)>=1
   %         for i = 1: plnJO.numOfModalities
   %              plnJO.propOpt.STscenarios(i)  =    pln(i).propOpt.STscenarios;
   %         end
   %     else
   %          plnJO.propOpt.STscenarios = ones( nPlans, 1);
   %     end 
   %     if isfield(plnJO.propOpt,'STfractions')
   %         for i = 1: plnJO.numOfModalities
   %              plnJO.propOpt.STfractions{i}  =    pln(i).propOpt.STfractions;
   %         end
   %     else
   %         for i = 1: plnJO.numOfModalities
   %             plnJO.propOpt.STfractions{i} = floor(pln(i).numOfFractions/plnJO.propOpt.STscenarios(i))* ones(plnJO.propOpt.STscenarios(i),1);
   %             plnJO.propOpt.STfractions{i} = plnJO.propOpt.STfractions{i}(end) + mod(pln(i).numOfFractions,plnJO.propOpt.STscenarios(i));
   %         end
   %     end 
   % else 
   %     plnJO.propOpt.spatioTemp = zeros(plnJO.numOfModalities,1);
   %     plnJO.propOpt.STscenarios = ones(plnJO.numOfModalities,1);
   %     for i = 1: plnJO.numOfModalities
   %          plnJO.propOpt.STfractions{i} = floor(pln(i).numOfFractions/pln(i).propOpt.STscenarios)* ones(pln(i).propOpt.STscenarios,1);
   %          plnJO.propOpt.STfractions{i} = plnJO.propOpt.STfractions{i}(end) + mod(pln(i).numOfFractions,plnJO.propOpt.STscenarios(i));
   %     end
   % end
       
   
else
   %Do nothing
   plnJO = pln;
end


end