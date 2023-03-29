function out = matRad_fieldConsistency(pln, field)
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
fieldcontainer = [];
   for k=1:size(pln,2)
      fieldcontainer = [fieldcontainer; fieldnames(pln(k).(field{:}))];
   end
   
   totfields = unique(fieldcontainer);
   for k=1:size(pln,2)
      isPlanField = isfield(pln(k).(field{:}),totfields);

      for m=1:size(isPlanField,1)
       if ~isPlanField(m)
          pln(k).(field{:}) = setfield(pln(k).(field{:}),totfields{m},[]);
       end
      end
   end
 out = [pln().(field{:})];
end