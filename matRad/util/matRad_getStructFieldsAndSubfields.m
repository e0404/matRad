function fieldStruct = matRad_getStructFieldsAndSubfields(s)
% matRad_getSubfolders: 
% Recursively loops into the struct fields and outputs fields and subfields
% call
%    fieldStruct = matRad_getStructFieldsAndSubfields(s)
%    input:
%       s:              input structure
%   output:
%       fieldStruct:   cell array containing subfields of s
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the
% help edit
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fieldStruct = fieldnames(s);
currLevelFields  = fieldStruct;
nFields  = numel(fieldStruct);

for subIdx=1:nFields
    currField = s.(currLevelFields{subIdx});
    if isstruct(currField)
        subFields = matRad_getStructFieldsAndSubfields(currField);
        fieldStruct = [fieldStruct; cellfun(@(sF) [fieldStruct{subIdx}, '.', sF], subFields, 'UniformOutput',false)];
    end
end

end