function pln = matRad_fieldConsistency(pln, field)
%Checks that the given field of the pln struct array contains stuctures
%that are consistent, i.e contain trhe same fields
% 
% call
%   pln = matRad_fieldConsistency(pln, field)
%
% input
%   pln:       array of pln structure for the different modalities
%   field:     specific pln field to be addressed
%
% output
%   pln
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end
fields = [];
   for k=1:size(pln,2)
      fields = [fields; fieldnames(pln(k).(field{:}))];
   end
   
   Totfields = unique(fields);
   for k=1:size(pln,2)
      IsPlanField = isfield(pln(k).(field{:}),Totfields);

      for m=1:size(IsPlanField,1)
       if ~IsPlanField(m)
          pln(k).(field{:}) = setfield(pln(k).(field{:}),Totfields{m},[]);
       end
      end
   end
end