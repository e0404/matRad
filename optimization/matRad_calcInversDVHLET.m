function LET = matRad_calcInversDVHLET(volume,LETVec)
% matRad inverse DVH LET (LET Volume Histogram) calculation
% 
% call
%   LET = matRad_calcInversDVHLET(volume,LETVec)
%
% input
%   volume:     rel volume of structure
%   LETVec:    LET vector of specific structure
%
% output
%   LET:       LET that corresponds to rel volume
%
% References
%   -
%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sort dose values
LETPoints = sort(LETVec,'descend');

ix = max([1 ceil(volume*numel(LETPoints))]);

LET = LETPoints(ix);

end
