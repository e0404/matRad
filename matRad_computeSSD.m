function stf = matRad_computeSSD(stf,ct,mode)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad SSD calculation
% 
% call
%   stf = matRad_computeSSD(stf,ct,mode)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   mode:           optional parameter specifying how to handle multiple
%                   cubes to compute one SSD
% output
%   stf:            matRad steering information struct
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    mode = 'first';
end

% set density threshold for SSD computation
densityThreshold = 0.05;

if strcmp(mode,'first')
    
    for i = 1:size(stf,2)
        SSD = cell(1,stf(i).numOfRays);
        for j = 1:stf(i).numOfRays
            [alpha,~,rho,~,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                                 ct.resolution, ...
                                 stf(i).sourcePoint, ...
                                 stf(i).ray(j).targetPoint, ...
                                 {ct.cube{1}});
            ixSSD = find(rho{1} > densityThreshold,1,'first');

            
            if isempty(ixSSD)
                warning('ray does not hit patient. Trying to fix afterwards...');
            elseif ixSSD(1) == 1
                warning('Surface for SSD calculation starts directly in first voxel of CT\n');
            end

            % calculate SSD
            SSD{j} = double(2 * stf(i).SAD * alpha(ixSSD));
            stf(i).ray(j).SSD = SSD{j};            
        end
        
        % try to fix SSD by using SSD of closest neighbouring ray
        SSDnotSet = find(cellfun('isempty',SSD));
        if ~isempty(SSDnotSet)
            rayPos_bev = reshape([stf(i).ray(:).rayPos_bev]',[3 stf(i).numOfRays])';
            for j = SSDnotSet
                stf(i).ray(j).SSD =  matRad_closestNeighbourSSD(rayPos_bev, SSD, rayPos_bev(j,:));
            end
        end
    end
else
    error('mode not defined for SSD calculation');
end


% default setting only use first cube
function bestSSD = matRad_closestNeighbourSSD(rayPos, SSD, currPos)
    vDistances = sum((rayPos - currPos).^2,2);
    [~, vIdx]   = sort(vDistances);
    for ix = vIdx'
        bestSSD = SSD{ix};
        % if SSD has been found, bestSSD is not empty
        if ~any(isempty(bestSSD))
            break
        end
    end
    if any(isempty(bestSSD))
        error('Could not fix SSD calculation.');
    end
end






end
