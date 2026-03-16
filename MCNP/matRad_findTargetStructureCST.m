function cstTargetIndex = matRad_findTargetStructureCST(cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB function to find the index of the radiotherapy target structure.
%
%
% call
%   cstTargetIndex = matRad_findTargetStructureCST(cst)
%
% input
%   cst
%
% output:
%   cstTargetIndex as array
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 01/2024
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cstTargetIndex = [];

for counter =1:size(cst,1)
    if strcmp(cst{counter,3}, 'TARGET')
        cstTargetIndex = [cstTargetIndex, counter];
    end
end

if isempty(cstTargetIndex)
    error('No target structure contoured or structure not found! Note: Target structure has to be set in matRad.')
end