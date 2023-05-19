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

nPlans = length(pln);
if nPlans>0
   originalPlans = pln;
   %Initialize Pln struct
   currentFields = fieldnames(pln(1));
   currentFields(:,2) = cell(size(currentFields)); %why the two cells 
   currentFields = currentFields';
   plnJO = struct(currentFields{:});
   
   %First run plnConsistency, it is the same as current MixModality, small
   %corrections are introduced to handle output generation with matRad_cfg
   pln = matRad_plnConsistency(pln); %to be reviewed
   
   %Define the Pln properties
   plnJO.numOfFractions = sum([pln(:).numOfFractions]); %Use total number of fractions
   plnJO.radiationMode = 'MixMod';
   plnJO.machine       = 'MixMod';
   plnJO.propStf       = [pln(:).propStf];
  
   for k=1:length(currentFields)
      if isempty(getfield(plnJO,currentFields{1,k})) && isstruct(pln(1).(currentFields{1,k}))
        %For all pln fields that are structures, check that the number of
        %fields is the same
        plnfld = matRad_fieldConsistency(pln,currentFields(1,k));
        if strcmp(currentFields{1,k},'propDoseCalc')
            plnJO.(currentFields{1,k}) = [plnfld(:,1)];                         % check for consistency and singletonnness of propDoseCalc and PropOpt
        elseif strcmp(currentFields{1,k},'propOpt')
            plnJO.(currentFields{1,k}) = [plnfld(:,1)];   
            
        else
            plnJO.(currentFields{1,k}) = [plnfld(:)];                         % check for consistency and singletonnness of propDoseCalc and PropOpt
        end
     
       
      elseif isempty(getfield(plnJO,currentFields{1,k})) && ~isstruct(pln(1).(currentFields{1,k}))
         %If the field is a class, just keep it
         plnJO.(currentFields{1,k}) = [pln(:).(currentFields{1,k})];
      end
   end
   plnJO.bioParam = matRad_bioModel(plnJO.radiationMode,pln(1).bioParam.quantityOpt, plnJO.radiationMode);
   %Save the original plans as well
   plnJO.originalPlans = originalPlans;
   plnJO.numOfModalities = nPlans;

   % SORT OUT THE ST Fields


   if isfield(plnJO.propOpt,'spatioTemp')

       for i = 1: plnJO.numOfModalities
            plnJO.propOpt.spatioTemp(i)  =    pln(i).propOpt.spatioTemp ;
       end 
       if isfield(plnJO.propOpt,'STscenarios') && sum(plnJO.propOpt.spatioTemp)>=1
           for i = 1: plnJO.numOfModalities
                plnJO.propOpt.STscenarios(i)  =    pln(i).propOpt.STscenarios;
           end
       else
            plnJO.propOpt.STscenarios = ones( nPlans, 1);
       end 
       if isfield(plnJO.propOpt,'STfractions')
           for i = 1: plnJO.numOfModalities
                plnJO.propOpt.STfractions{i}  =    pln(i).propOpt.STfractions;
           end
       else
           for i = 1: plnJO.numOfModalities
               plnJO.propOpt.STfractions{i} = floor(pln(i).numOfFractions/plnJO.propOpt.STscenarios(i))* ones(plnJO.propOpt.STscenarios(i),1);
               plnJO.propOpt.STfractions{i} = plnJO.propOpt.STfractions{i}(end) + mod(pln(i).numOfFractions,plnJO.propOpt.STscenarios(i));
           end
       end 
   else 
       plnJO.propOpt.spatioTemp = zeros(plnJO.numOfModalities,1);
       plnJO.propOpt.STscenarios = ones(plnJO.numOfModalities,1);
       for i = 1: plnJO.numOfModalities
            plnJO.propOpt.STfractions{i} = floor(pln(i).numOfFractions/pln(i).propOpt.STscenarios)* ones(pln(i).propOpt.STscenarios,1);
            plnJO.propOpt.STfractions{i} = plnJO.propOpt.STfractions{i}(end) + mod(pln(i).numOfFractions,plnJO.propOpt.STscenarios(i));
       end
   end
       
   
else
   %Do nothing
   plnJO = pln;
end


end