function jacobStruct = matRad_getJacobStruct(dij,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: jacobian structure function for inverse planning supporting max dose
% constraint, min dose constraint, min max dose constraint, min mean, max
% min, min max mean constraint, min EUD constraint, max EUDconstraint, 
% min max EUD constraint, exact DVH constraint, max DVH constraint, 
% min DVH constraint 
% 
% call
%   jacobStruct = matRad_getJacobStruct(dij,cst)
%
% input
%   dij: dose influence matrix
%   cst: matRad cst struct
%
% output
%   jacobStruct: jacobian of constraint function
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% check consitency of constraints
for i = 1:size(cst,1)    
    for j = 1:numel(cst{i,6})
        if isequal(cst{i,6}(j).type, 'max mean dose constraint') || ...
           isequal(cst{i,6}(j).type, 'min mean dose constraint') || ...
           isequal(cst{i,6}(j).type, 'min max mean dose constraint')

                % -> no other max, min or min max mean dose constraint
                for k = j+1:numel(cst{i,6})
                    if isequal(cst{i,6}(k).type, 'max mean dose constraint') || ...
                       isequal(cst{i,6}(k).type, 'min mean dose constraint') || ...
                       isequal(cst{i,6}(k).type, 'min max mean dose constraint')
                            error('Simultatenous definition of min, max and or min max mean dose constraint\n');
                    end
                end
        elseif isequal(cst{i,6}(j).type, 'max EUD constraint') || ...
               isequal(cst{i,6}(j).type, 'min EUD constraint') || ...
               isequal(cst{i,6}(j).type, 'min max EUD constraint')

                % -> no other max, min or min max mean dose constraint
                for k = j+1:numel(cst{i,6})
                    if isequal(cst{i,6}(k).type, 'max EUD constraint') || ...
                       isequal(cst{i,6}(k).type, 'min EUD constraint') || ...
                       isequal(cst{i,6}(k).type, 'min max EUD constraint')
                            error('Simultatenous definition of min, max and or min max EUD constraint\n');
                    end
                end

        elseif isequal(cst{i,6}(j).type, 'exact DVH constraint') ||...
               isequal(cst{i,6}(j).type, 'max DVH constraint') ||...
               isequal(cst{i,6}(j).type, 'min DVH constraint')

            % -> no other DVH constraint
            for k = j+1:numel(cst{i,6})
                if (isequal(cst{i,6}(k).type, 'exact DVH constraint') && isequal(cst{i,6}(j).dose,cst{i,6}(k).dose)) || ...
                   (isequal(cst{i,6}(k).type, 'exact DVH constraint') && isequal(cst{i,6}(j).volume,cst{i,6}(k).volume)) || ... 
                   (isequal(cst{i,6}(k).type, 'max DVH constraint')   && isequal(cst{i,6}(j).dose,cst{i,6}(k).dose)) || ...
                   (isequal(cst{i,6}(k).type, 'max DVH constraint')   && isequal(cst{i,6}(j).volume,cst{i,6}(k).volume)) || ... 
                   (isequal(cst{i,6}(k).type, 'min DVH constraint')   && isequal(cst{i,6}(j).dose,cst{i,6}(k).dose)) || ...
                   (isequal(cst{i,6}(k).type, 'min DVH constraint')   && isequal(cst{i,6}(j).volume,cst{i,6}(k).volume))

                        error('Simultatenous definition of DVH constraint\n');
                end
            end    
        end
    end
end

% Initializes constraints
jacobStruct = sparse([]);

% loop over all scenarios
for i = 1:dij.numOfScenarios

    % compute objective function for every VOI.
    for  j = 1:size(cst,1)

        % Only take OAR or target VOI.
        if ~isempty(cst{j,4}) && ( isequal(cst{j,3},'OAR') || isequal(cst{j,3},'TARGET') )

            % loop over the number of constraints for the current VOI
            for k = 1:numel(cst{j,6})
                
                % only in the nominal case or for robust optimization
                if i == 1 || strcmp(cst{j,6}(k).robustness,'probabilistic') || ...
                             strcmp(cst{j,6}(k).robustness,'voxel-wise worst case')

                    if isequal(cst{j,6}(k).type, 'max dose constraint') || ...
                       isequal(cst{j,6}(k).type, 'min dose constraint') || ...
                       isequal(cst{j,6}(k).type, 'max mean dose constraint') || ...
                       isequal(cst{j,6}(k).type, 'min mean dose constraint') || ...
                       isequal(cst{j,6}(k).type, 'min max mean dose constraint') || ...
                       isequal(cst{j,6}(k).type, 'max EUD constraint') || ...
                       isequal(cst{j,6}(k).type, 'min EUD constraint') || ...
                       isequal(cst{j,6}(k).type, 'min max EUD constraint') || ...
                       isequal(cst{j,6}(k).type, 'exact DVH constraint') || ...
                       isequal(cst{j,6}(k).type, 'max DVH constraint') || ... 
                       isequal(cst{j,6}(k).type, 'min DVH constraint')

                       jacobStruct = [jacobStruct; spones(mean(dij.physicalDose{i}(cst{j,4},:)))];
                    
                    end

                end

            end

        end
        
    end
   
end