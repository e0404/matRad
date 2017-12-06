function stf = matRad_computeSSD(ct,stf,pln,mode)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad SSD calculation
% 
% call
%   stf = matRad_computeSSD(stf,ct,mode)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   pln:            matRad plan meta information struct
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

function bestSSD = closestNeighbourSSD(rayPos, SSD, currPos)
    distances = sum((rayPos - currPos).^2,2);
    [~, idx] = sort(distances);
    for a = 1:numel(idx)
        bestSSD = SSD{idx(a)};
        % if SSD has been found, bestSSD is not empty
        if ~any(isempty(bestSSD))
            break
        end
    end
    if any(isempty(bestSSD))
        error('Could not fix SSD calculation.');
    end
end


% default setting only use first cube
if nargin < 4
    mode = 'first';
end

% set density threshold for SSD computation
densityThreshold = 0.05;

if strcmp(mode,'first')
    
    for CtScen = 1:pln.multScen.numOfCtScen
       for i = 1:size(stf,2)
           SSD = cell(1,stf(i).numOfRays);
           for j = 1:stf(i).numOfRays
               [alpha,~,rho,~,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                                    ct.resolution, ...
                                    stf(i).sourcePoint, ...
                                    stf(i).ray(j).targetPoint, ...
                                    {ct.cube{CtScen}});
               ixSSD = find(rho{1} > densityThreshold,1,'first');

            if isempty(ixSSD)
                warning('Ray is off patient. Trying to fix afterwards...');

            elseif ixSSD(1) == 1
                warning('Surface for SSD calculation starts directly in first voxel of CT\n');
            end
            
            SSD{j} = double(2 * stf(i).SAD * alpha(ixSSD));
            stf(i).ray(j).SSD{CtScen} = SSD{j};            
            end
           
           % bestSSD if not already set
           SSDnotSet = find(cellfun('isempty',SSD));
           if ~isempty(SSDnotSet)
               rayPos_bev = reshape([stf(i).ray(:).rayPos_bev],[],3);
               for j = SSDnotSet
                    pos = rayPos_bev(j,:);
                    stf(i).ray(j).SSD{CtScen} = closestNeighbourSSD(rayPos_bev, SSD, pos);
               end
           end
       end
    end

else
    error('mode not defined for SSD calculation');
end

end

