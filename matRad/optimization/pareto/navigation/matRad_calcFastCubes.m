function cubes = matRad_calcFastCubes(w,dij,pln)
% matRad computation of cube for plan optimization quantity
%
% call
%   resultGUI = matRad_calcFastCubes(w,dij,pln)
%
% input
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   pln:     matRad pln struct
%
% output
%   cubes:   Array storing cubes of optimization quantity
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch pln.bioParam.quantityOpt
        case 'physicalDose'
            cubes = reshape(dij.physicalDose{1}*w,dij.doseGrid.dimensions);
        case 'RBExD'
            if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
                %TODO                
                effect = (dij.mAlphaDose{1} * w + (dij.mSqrtBetaDose{1} * w).^2);
                ix = dij.bx{1}~=0 & effect(:) > 0;
                cubes = zeros(dij.doseGrid.dimensions);  
                cubes(ix) = (sqrt(dij.ax{1}(ix).^2 + 4 .* dij.bx{1}(ix).* effect(ix)) - dij.ax{1}(ix))./(2.*dij.bx{1}(ix));

            else
                cubes = reshape(dij.physicalDose{1}*w,dij.doseGrid.dimensions)*dij.RBE;
            end
        case 'effect'
            cubes = reshape(full(dij.mAlphaDose{1} * w + (dij.mSqrtBetaDose{1} * w).^2),dij.doseGrid.dimensions);
    end
    cubes = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
                                     cubes, ...
                                     dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);
end