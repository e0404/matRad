function Pln = matRad_plnWrapper(pln)
% Combines the arbitrary number of input plans into a single plan for the
% different modalities.
% 
% call
%   Pln = matRad_plnWrapper(pln)
%
% input
%   pln:       array of pln structure for the different modalities (if any)
%
% output
%   Pln
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end


nPlans = length(pln);
if nPlans>1
   OriginalPlans = pln;
   %Initialize Pln struct
   CurrentFields = fieldnames(pln(1));
   CurrentFields(:,2) = cell(size(CurrentFields));
   CurrentFields = CurrentFields';
   Pln = struct(CurrentFields{:});
   
   %First run plnConsistency, it is the same as current MixModality, small
   %corrections are introduced to handle output generation with matRad_cfg
   pln = matRad_plnConsistency(pln); %to be reviewed
   
   %Define the Pln properties
   Pln.numOfFractions = sum([pln(:).numOfFractions]); %Use total number of fractions
   Pln.radiationMode = 'MixMod';
   Pln.machine       = 'MixMod';
   Pln.propStf       = [pln(:).propStf];
  
   for k=1:length(CurrentFields)
      if isempty(getfield(Pln,CurrentFields{1,k})) && isstruct(pln(1).(CurrentFields{1,k}))
        %For all pln fields that are structures, check that the number of
        %fields is the same
        pln = matRad_fieldConsistency(pln,CurrentFields(1,k));
        Pln.(CurrentFields{1,k}) = [pln(:).(CurrentFields{1,k})];
      elseif isempty(getfield(Pln,CurrentFields{1,k})) && ~isstruct(pln(1).(CurrentFields{1,k}))
         %If the field is a class, just keep it
         Pln.(CurrentFields{1,k}) = [pln(:).(CurrentFields{1,k})];
      end
   end

   %Save the original plans as well
   Pln.OriginalPlans = OriginalPlans;
   Pln.numOfModalities = nPlans;
else
   %Do nothing
   Pln = pln;
end


end