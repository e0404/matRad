function [cl,cu] = matRad_getConstBounds(cst,numOfScenarios,type)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT get constraint bounds function
% 
% call
%   [cl,cu] = matRad_getConstBounds(cst,type)
%
% input
%   cst:  matRad cst struct
%   numOfScenarios: number of scenarios considered during optimization
%   type: type of optimizaiton; either 'none','effect' or 'RBExD'
%
% output
%   cl: lower bounds on constraints
%   cu: lower bounds on constraints
%
% References
%
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

% Initialize bounds
cl = [];
cu = [];

% loop over all scenarios
for i = 1:numOfScenarios
    
    % compute objective function for every VOI.
    for  j = 1:size(cst,1)

        % Only take OAR or target VOI.
        if ~isempty(cst{j,4}) && ( isequal(cst{j,3},'OAR') || isequal(cst{j,3},'TARGET') )

            % loop over the number of constraints for the current VOI
            for k = 1:numel(cst{j,6})
               
                % only in the nominal case or for robust optimization
                if i == 1 || strcmp(cst{j,6}(k).robustness,'probabilistic') || ...
                             strcmp(cst{j,6}(k).robustness,'voxel-wise worst case')

                    if isequal(type,'none') || isequal(type,'RBExD') 
                        param = cst{j,6}(k).dose;
                    elseif isequal(type,'effect')
                        param = cst{j,5}.alphaX .* cst{j,6}(k).dose + cst{j,5}.betaX .* cst{j,6}(k).dose.^2;
                    end
                    
                    if isequal(cst{j,6}(k).type, 'max dose constraint') 
             
                        cl = [cl;-inf];
                        cu = [cu;param];

                    elseif isequal(cst{j,6}(k).type, 'min dose constraint') 

                        cl = [cl;param];
                        cu = [cu;inf];

                    elseif isequal(cst{j,6}(k).type, 'min mean dose constraint') 
                    
                        cl = [cl;param];
                        cu = [cu;inf];

                    elseif isequal(cst{j,6}(k).type, 'max mean dose constraint') 

                        cl = [cl;-inf];
                        cu = [cu;param];

                    elseif isequal(cst{j,6}(k).type, 'min max mean dose constraint') 

                        cl = [cl;param(1)];
                        cu = [cu;param(2)];

                    elseif isequal(cst{j,6}(k).type, 'min EUD constraint') 

                        cl = [cl;param];
                        cu = [cu;inf];

                    elseif isequal(cst{j,6}(k).type, 'max EUD constraint') 

                        cl = [cl;-inf];
                        cu = [cu;param];

                    elseif isequal(cst{j,6}(k).type, 'min max EUD constraint') 

                        cl = [cl;param(1)];
                        cu = [cu;param(2)];

                    elseif isequal(cst{j,6}(k).type, 'exact DVH constraint')

                        cl = [cl;cst{j,6}(k).volume/100];
                        cu = [cu;cst{j,6}(k).volume/100];

                    elseif isequal(cst{j,6}(k).type, 'max DVH constraint') 

                        cl = [cl;-inf];
                        cu = [cu;cst{j,6}(k).volume/100];

                        % alternative constraint calculation 1/4 %                
                        % cl = [cl;-inf];
                        % cu = [cu;0];
                        % alternative constraint calculation 1/4 %

                    elseif isequal(cst{j,6}(k).type, 'min DVH constraint') 

                        cl = [cl;cst{j,6}(k).volume/100];
                        cu = [cu;inf];

                        % alternative constraint calculation 2/4 %                
                        % cl = [cl;-inf];
                        % cu = [cu;0];
                        % alternative constraint calculation 2/4 %

                    end % constraint switch
                    
                end % if nominal scenario or robust opt
                
            end % over all objectives of structure

        end % if structure not empty and target or oar
        
    end % over all structures
   
end % over all scenarios

