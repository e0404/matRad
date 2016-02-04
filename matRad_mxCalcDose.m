function d = matRad_mxCalcDose(dij,w,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix based dose calculation and for biological base data also RBE,
% alpha, effect, etc.
% 
% call
%   d = matRad_mxCalcDose(dij,w)
%
% input
%   dij:    dose influence matrix
%   w:      bixel weight vector
%   cst:    matRad cst struct
%
% output
%   d:      struct containing the dose cube. If alpha and beta cubes are 
%           also available, rbe cube, RBE-weighted dose etc. will also be
%           calculated
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the Eclipse Public License 1.0 (EPL-1.0).
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.
%
% You should have received a copy of the EPL-1.0 in the file license.txt
% along with matRad. If not, see <http://opensource.org/licenses/EPL-1.0>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d.w = w;
d.physicalDose = reshape(dij.physicalDose*w,dij.dimensions);

if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
    
    fprintf('Calculating alpha/beta/effect/cube...');
   
    a_x = zeros(size(d.physicalDose));
    b_x = zeros(size(d.physicalDose));
    for  i = 1:size(cst,1)
        % Only take OAR or target VOI.
        if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') 
            a_x(cst{i,4}) = cst{i,5}.alphaX;
            b_x(cst{i,4}) = cst{i,5}.betaX;
        end
    end
    
    d.effect = (dij.mAlphaDose*w+(dij.mSqrtBetaDose*w).^2);
    d.effect = reshape(d.effect,dij.dimensions);
    
    d.RBEWeightedDose = ((sqrt(a_x.^2 + 4 .* b_x .* d.effect) - a_x)./(2.*b_x));
    d.RBE = d.RBEWeightedDose./d.physicalDose;
    
    d.alpha = (dij.mAlphaDose.*spfun(@(x)1./x,dij.physicalDose)) * d.w;
    d.alpha = reshape(d.alpha,dij.dimensions);
    d.beta = ( (dij.mSqrtBetaDose.*spfun(@(x)1./x,dij.physicalDose)) * d.w ).^2;
    d.beta = reshape(d.beta,dij.dimensions);
    
    fprintf(' done!\n')
end

