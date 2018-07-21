function [vMeanRBExD,vStdRBExD ] = matRad_effectToRBExDose(vMeanEffect, vStdEffect ,vAlpha_x,vBeta_x)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  matRad function to propagate the first and second moment of the biological effect
%  to first and second moment of the RBE-weighted dose
%
%  call
%   [vMeanRBExD,vStdRBExD ] = matRad_effectToRBExDose(vMeanEffect, vStdEffect ,vAlpha_x,vBeta_x)
%
%
% input
%   vMeanEffect:            first raw moment of the biological effect
%   vStdEffect:             second central moment of the biological effect
%   vAlpha_x                   reference radio-sensitivity parameter of photons
%   vBeta_x                    reference radio-sensitivity parameter of photons
%
% output
%   vMeanRBExD:           matRad's bioParam structure containing information
%   vStdRBExD                   about the choosen biological model
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vMeanRBExD = zeros(size(vMeanEffect));
vStdRBExD  = zeros(size(vMeanEffect));

gamma = vAlpha_x./(2*vBeta_x);
ix    = gamma > 0 & ~isinf(gamma) ;
x0    = vMeanEffect;
A0    =  sqrt(vBeta_x(ix).^-1 .* x0(ix) + gamma(ix).^2) - gamma(ix);
A1    =  ( 1 * vBeta_x(ix).^-1 )./( 2*sqrt(vBeta_x(ix).^-1 .* x0(ix) + gamma(ix).^2));           % 1 
A2    = -( 1 * vBeta_x(ix).^-2 )./( 4*((vBeta_x(ix).^-1    .* x0(ix) + gamma(ix).^2).^(3/2))) ;  % 2
A4    = -(15 * vBeta_x(ix).^-4 )./(16*((vBeta_x(ix).^-1    .* x0(ix) + gamma(ix).^2).^(7/2))) ;  % 24
    
vMeanRBExD(ix) =   A0 + (1/2)*A2.*vStdEffect(ix).^2 + (1/24)*A4.*vStdEffect(ix).^4;
vStdRBExD(ix)  =   sqrt(A1.^2 .* vStdEffect(ix).^2  +  (0.5 * A2.^2 .* vStdEffect(ix).^4));


end

