function [allMatch, msg] = matRad_compareDijStf(stf,dij)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% call
%   matching = matRad_comparePlnDijStf(pln,stf,dij)
%
% input
%   dij:                        matRad dij struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%
% output
%
%   allMatch:                   flag is true if they all match
%   matching:                   message to display
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allMatch=true;
msg = [];

    %% compare number of rays per beam in dij and stf
    stf_RaysPerBeam=[stf.numOfRays];
    if numel(stf_RaysPerBeam) ~= numel(dij.numOfRaysPerBeam) ... % different size
            || ~isempty(find(stf_RaysPerBeam-dij.numOfRaysPerBeam,1)) % different values
        msg= 'Number of rays do not match';
        allMatch=false;
        return
    end
    stf_gantryAngles=[stf.gantryAngle];
    if dij.numOfBeams ~= numel(stf_gantryAngles)
        msg= 'Number of beams do not match';
        allMatch=false;
        return
    end

end