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

% default setting only use first cube
if nargin < 3
    mode = 'first';
end

% set density threshold for SSD computation
densityThreshold = 0.05;

if strcmp(mode,'first')
    
    for i = 1:size(stf,2)
        for j = 1:stf(i).numOfRays
            [alpha,~,rho,~,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                                 ct.resolution, ...
                                 stf(i).sourcePoint, ...
                                 stf(i).ray(j).targetPoint, ...
                                 {ct.cube{1}});
            ixSSD = find(rho{1} > densityThreshold,1,'first');

            if isempty(ixSSD)== 1
                warning('Surface for SSD calculation starts directly in first voxel of CT\n');
            end

            % calculate SSD
            stf(i).ray(j).SSD = double(2 * stf(i).SAD * alpha(ixSSD));

        end
    end
    
    % include rangeshifter data if not yet available 
    for i = 1:size(stf,2)
        for j = 1:stf(i).numOfRays
            if ~isfield(stf(i).ray(j), 'rashiEqThickness') || ~isfield(stf(i).ray(j), 'rashiIsoCenterDistance') ...
                    || isempty(stf(i).ray(j).rashiEqThickness) || isempty(stf(i).ray(j).rashiIsoCenterDistance)
                for k = 1:numel(stf(i).ray(j).energy)
                    stf(i).ray(j).rangeShifter(k).ID = 0;
                    stf(i).ray(j).rangeShifter(k).eqThickness = 0;
                    stf(i).ray(j).rangeShifter(k).sourceRashiDistance = 0;
                end
            end
        end
    end

else
    
    error('mode not defined for SSD calculation');
    
end
