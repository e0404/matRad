function hessianStruct = matRad_getHessianStruct(dij,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: hessian structure function for inverse planning 
% supporting squared underdosage, squared overdosage, squared deviation, 
% mean dose objectives, EUD objectives, DVH objectives, (exact) max dose 
% constraint, (exact) min dose constraint, min mean dose constraint, max 
% mean dose constraint, min EUD constraint, max EUD constraint, max DVH 
% constraint, min DVH constraint 
% 
% call
%   hessianStruct = matRad_getHessianStruct(dij,cst)
%
% input
%   dij: dose influence matrix
%   cst: matRad cst struct
%
% output
%   hessianStruct: hessian matrix (sparse, lower triangular matrix)
%
% References
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

% check consitency of constraints
for i = 1:size(cst,1)    
    for j = 1:numel(cst{i,6})
        if isequal(cst{i,6}(j).type, 'max mean dose constraint') || ...
           isequal(cst{i,6}(j).type, 'min mean dose constraint')

                % -> no other max, min mean dose constraint
                for k = j+1:numel(cst{i,6})
                    if isequal(cst{i,6}(k).type, 'max mean dose constraint') || ...
                       isequal(cst{i,6}(k).type, 'min mean dose constraint') 
                            error('Simultatenous definition of min mean and max mean dose constraint\n');
                    end
                end
        elseif isequal(cst{i,6}(j).type, 'max EUD constraint') || ...
               isequal(cst{i,6}(j).type, 'min EUD constraint') 

                % -> no other max, min mean dose constraint
                for k = j+1:numel(cst{i,6})
                    if isequal(cst{i,6}(k).type, 'max EUD constraint') || ...
                       isequal(cst{i,6}(k).type, 'min EUD constraint')
                            error('Simultatenous definition of min EUD and max EUD constraint\n');
                    end
                end

        elseif isequal(cst{i,6}(j).type, 'max DVH constraint') ||...
               isequal(cst{i,6}(j).type, 'min DVH constraint')

            % -> no other DVH constraint
            for k = j+1:numel(cst{i,6})
                if (isequal(cst{i,6}(k).type, 'max DVH constraint') && isequal(cst{i,6}(k).dose,cst{i,6}(j).dose)) || ...
                   (isequal(cst{i,6}(k).type, 'max DVH constraint') && isequal(cst{i,6}(k).volume,cst{i,6}(j).volume)) || ... 
                   (isequal(cst{i,6}(k).type, 'min DVH constraint') && isequal(cst{i,6}(k).dose,cst{i,6}(j).dose)) || ...
                   (isequal(cst{i,6}(k).type, 'min DVH constraint') && isequal(cst{i,6}(k).volume,cst{i,6}(j).volume))

                        error('Simultatenous definition of DVH constraint\n');
                end
            end    
        end
    end
end

% initialize hessian structure (zero, only linear constraints/objectives)
hessianStruct = sparse(dij.totalNumOfBixels, dij.totalNumOfBixels);

% check whether there are non-linear constraints/objectives
for i = 1:size(cst,1)

    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints/objectives for the current VOI
        for j = 1:numel(cst{i,6})

            % if conventional opt: just check for nominal dose
            if strcmp(cst{i,6}(j).robustness,'none')

               if ismember(cst{i,6}(j).type, {'square underdosing','square overdosing','square deviation', 'EUD',...
                        'min dose constraint','max dose constraint',...
                        'min EUD constraint','max EUD constraint',...
                        'max DVH constraint','min DVH constraint',...
                        'max DVH objective' ,'min DVH objective'})

                   hessianStruct = sparse(tril(ones(dij.totalNumOfBixels)));
                   
                   % return once non-linear objective/constraint has been encountered
                   return

                end

            end

        end

    end

end
  