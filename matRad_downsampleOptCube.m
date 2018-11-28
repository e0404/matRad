function [dij, cst] = matRad_downsampleOptCube(cst, dij ) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to downsample the optimization cube
% checker board method
% call
%   dij, cst] = matRad_downsampleOptCube(cst, dij, k ) 
%
% input
%   dij:     dose influence matrix
%   cst:     matRad cst struct
%   
%
% output
%   dij:     downsample dij matrices
%   cst:     cst with relevant samples
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


% checker board  3d mask 
check_mask =[];
for i=1: dij.dimensions(1)
    for j=1:dij.dimensions(2)
        for k=1:dij.dimensions(3)
            check_mask(i,j,k) = rem(sum([i j k]),2);
        end 
    end 
end
chid = find(check_mask);
% select idx from each structure therefore MUST have body
% contour
select_idx = [];
for i= 1: size(cst,1)
   if ~(numel(cst{i,4}{1})*prod([dij.resolution.x dij.resolution.y dij.resolution.z])/1000 < 100)            %only to do if volume of VOI < 100cm^2
       [~,vx] = intersect(cst{i,4}{1}, chid);
       select_idx = [select_idx; cst{i,4}{1}(vx)];
   else
       select_idx = [select_idx; cst{i,4}{1}];
   end 
   
end 

select_idx = unique(select_idx);



% altering dij
dij.numOfVoxels = numel(select_idx);
if isfield(dij,'alphaX') && isfield(dij,'betaX')
    dij.alphaX = dij.alphaX(select_idx);
    dij.betaX = dij.betaX(select_idx);
    dij.alphaXspare = dij.alphaXspare(select_idx);
    dij.betaXspare = dij.betaXspare(select_idx);
    dij.abX   = dij.abX(select_idx);
end

for i = 1: numel(dij.physicalDose) 
   dij.physicalDose{i}  =   dij.physicalDose{i}(select_idx,:);
end 

if isfield(dij, 'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
    for i = 1: numel(dij.mAlphaDose) 
        dij.mAlphaDose{i}    = dij.mAlphaDose{i}(select_idx,:);
        dij.mSqrtBetaDose{i} = dij.mSqrtBetaDose{i}(select_idx,:);
        
    end 
    dij.mAlphaDoseSpare{1} = dij.mAlphaDoseSpare{1}(select_idx,:);
    dij.mSqrtBetaDoseSpare{1} = dij.mSqrtBetaDoseSpare{1}(select_idx,:);
end 

% new map for cst
for i= 1: size(cst,1)
    [~,cst{i,4}{1}] = intersect(select_idx, cst{i,4}{1});   
end 

end 