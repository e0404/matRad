function [allMatch, msg] = matRad_comparePlnStf(pln,stf)
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allMatch=true;
msg = [];

%% check if steering information is available in plan from the begining
if ~isfield(pln,'propStf')
    allMatch=false;
    msg= 'No steering information in plan';
    return
end

%% compare number of gantry angles, but ignore if numOfBeams not set
if isfield(pln.propStf,'numOfBeams') && pln.propStf.numOfBeams ~= numel(stf) 
        msg= 'Number of beams do not match';
        allMatch=false;
        return
end

%% compare gantry angles in  stf and pln
stfGantryAngles=[stf.gantryAngle];
if ~isfield(pln.propStf,'gantryAngles') || numel(stfGantryAngles) ~= numel(pln.propStf.gantryAngles) ... % different size
        || ~isempty(find(stfGantryAngles-pln.propStf.gantryAngles, 1))  % values in stf and pln do not match  % values in stf and pln do not match
    allMatch=false;
    msg= 'Gantry angles do not match';
    return
end

%% compare couch angles in stf and pln
stfCouchAngles=[stf.couchAngle];
if ~isfield(pln.propStf,'couchAngles') || numel(stfCouchAngles) ~= numel(pln.propStf.couchAngles) ... % different size
        || ~isempty(find(stfCouchAngles-pln.propStf.couchAngles, 1))  % values in stf and pln do not match
    allMatch=false;
    msg= 'Couch angles do not match';
    return
end

%% compare Bixel width in stf and pln
bixelMatch = false;
if isfield(pln.propStf,'bixelWidth') && isfield(stf(1),'bixelWidth')
    if isnumeric(stf(1).bixelWidth) && isequal(stf(1).bixelWidth,pln.propStf.bixelWidth)
        bixelMatch = true;
    elseif ischar(stf(1).bixelWidth) && strcmp(stf(1).bixelWidth,'field')
        bixelMatch = true;
    end
end

if ~bixelMatch
    allMatch=false;
    msg= 'Bixel width does not match';
    return
end

%% compare radiation mode in stf and pln
if ~isfield(pln,'radiationMode') || ~strcmp(stf(1).radiationMode, pln.radiationMode)
    allMatch=false;
    msg= 'Radiation mode does not match';
    return
end

%% compare isocenter in stf and pln for each gantry angle
for i = 1:numel(pln.propStf.gantryAngles)
    if size(pln.propStf.isoCenter,1) == 1
        isoCenter = repmat(pln.propStf.isoCenter,numel(stf),1);
    else
        isoCenter = pln.propStf.isoCenter;
    end
    if size(isoCenter,1) ~= numel(stf) || any(stf(i).isoCenter - isoCenter(i,:) ~= 0)
        allMatch=false;
        msg= 'Isocenters do not match';
        return
    end
end

end