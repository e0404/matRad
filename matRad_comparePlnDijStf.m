function [allMatch, msg] = matRad_comparePlnDijStf(pln,stf,dij)
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


%% compare gantry angles in dij, stf and pln
stf_gantryAngles=[stf.gantryAngle];
if dij.numOfBeams ~= numel(stf_gantryAngles) || dij.numOfBeams ~= numel(pln.propStf.gantryAngles) ... % different size
        || ~isempty(find(stf_gantryAngles-pln.propStf.gantryAngles, 1))  % values in stf and pln do not match
    allMatch=false;
    msg= 'Gantry angles do not match';
    return
end

%% compare couch angles in stf and pln
stf_couchAngles=[stf.couchAngle];
if numel(stf_couchAngles) ~= numel(pln.propStf.couchAngles) ... % different size
        || ~isempty(find(stf_couchAngles-pln.propStf.couchAngles, 1))  % values in stf and pln do not match
    allMatch=false;
    msg= 'Couch angles do not match';
    return
end

%% compare Bixel width in stf and pln
if stf(1).bixelWidth ~= pln.propStf.bixelWidth
    allMatch=false;
    msg= 'Bixel width does not match';
    return
end

%% compare radiation mode in stf and pln
if ~strcmp(stf(1).radiationMode, pln.radiationMode)
    allMatch=false;
    msg= 'Radiation mode does not match';
    return
end

%% compare isocenter in stf and pln for each gantry angle
for i = 1:numel(pln.propStf.gantryAngles)
    if ~isempty(find(stf(i).isoCenter - pln.propStf.isoCenter(i,:) ,1))
        allMatch=false;
        msg= 'Isocenter does not match';
        return
    end
end

%% compare number of rays per beam in dij and stf
stf_RaysPerBeam=[stf.numOfRays];
if numel(stf_RaysPerBeam) ~= numel(dij.numOfRaysPerBeam) ... % different size
        || ~isempty(find(stf_RaysPerBeam-dij.numOfRaysPerBeam,1)) % different values
    msg= 'Number of rays does not match';
    allMatch=false;
    return
end

end