function stf = matRad_computeSSD(stf,ct,varargin)
% matRad SSD calculation
% 
% call
%   stf = matRad_computeSSD(stf,ct)
%   stf = matRad_computeSSD(stf,ct,Name,Value)
%
% input
%   ct:     ct cube
%   stf:    matRad steering information struct
%   
%   Optional Name/Value Properties:
%   mode:   optional parameter specifying how to handle multiple
%           cubes to compute one SSD. Only 'first' isimplemented
%
%   densityThreshold: value determining the skin threshold.
%
% output
%   stf:    matRad steering information struct
%
% References
%   -
%
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


matRad_cfg = MatRad_Config.instance();

%Parse arguments
p = inputParser();
%We leave the required parameters to avoid duplicates
%p.addRequired('stf',@isstruct);
%p.addRequired('ct',@isstruct);
p.addParameter('mode','first');
p.addParameter('densityThreshold',matRad_cfg.propDoseCalc.defaultSsdDensityThreshold,@(x) isnumeric(x) && isscalar(x));
p.addParameter('showWarning',true,@(x) islogical(x) && isscalar(x));

p.parse(varargin{:});

mode             = p.Results.mode;
densityThreshold = p.Results.densityThreshold;
boolShowWarning  = p.Results.showWarning;



if strcmp(mode,'first')
    
    for i = 1:size(stf,2)
        SSD = cell(1,stf(i).numOfRays);
        for j = 1:stf(i).numOfRays
            [alpha,~,rho,d12,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                                 ct.resolution, ...
                                 stf(i).sourcePoint, ...
                                 stf(i).ray(j).targetPoint, ...
                                 {ct.cube{1}});
            ixSSD = find(rho{1} > densityThreshold,1,'first');

            if boolShowWarning
                if isempty(ixSSD)
                    matRad_cfg.dispWarning('ray does not hit patient. Trying to fix afterwards...');
                    boolShowWarning = false;
                elseif ixSSD(1) == 1
                    matRad_cfg.dispWarning('Surface for SSD calculation starts directly in first voxel of CT!');
                    boolShowWarning = false;
                end
            end
            
            % calculate SSD
            SSD{j} = double(d12* alpha(ixSSD));
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
    matRad_cfg.dispError('mode not defined for SSD calculation');
end

end

% default setting only use first cube
function bestSSD = matRad_closestNeighbourSSD(rayPos, SSD, currPos)
    vDistances = sum((rayPos - repmat(currPos,size(rayPos,1),1)).^2,2);
    [~, vIdx]   = sort(vDistances);
    for ix = vIdx'
        bestSSD = SSD{ix};
        % if SSD has been found, bestSSD is not empty
        if ~any(isempty(bestSSD))
            break
        end
    end
    if any(isempty(bestSSD))
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Could not fix SSD calculation.');
    end
end

