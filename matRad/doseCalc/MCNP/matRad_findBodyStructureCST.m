function cstBodyIndex = matRad_findBodyStructureCST(cst, bodyStructureName)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB function to find the index of the radiotherapy structure body.
%
% Note: Body structure is the only structure that has to be contoured.
%
% call
%   cstBodyIndex = matRad_findBodyStructureCST(cst)
%
% input
%   cst
%   bodyStructureName
%
% output:
%   cstBodyIndex
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 05/2019
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

cstBodyIndex = 1;

while ~strcmpi(cst{cstBodyIndex,2}, bodyStructureName)
    cstBodyIndex = cstBodyIndex +1;
    if cstBodyIndex > size(cst,1)
        matRad_cfg.dispError('No body structure contoured or structure not found! Note: Body structure has to be named BODY (case insensitive).');
    end
end