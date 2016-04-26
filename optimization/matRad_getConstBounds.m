function [cl,cu] = matRad_getConstBounds(cst,type)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT get constraint bounds function
% 
% call
%   [cl,cu] = matRad_getConstBounds(cst,type)
%
% input
%   cst:  matRad cst struct
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

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
        
        % check consitency of constraints
        for j = 1:numel(cst{i,6})
            if isequal(cst{i,6}(j).type, 'max dose constraint') || ...
               isequal(cst{i,6}(j).type, 'min dose constraint') || ...
               isequal(cst{i,6}(j).type, 'min max dose constraint')
           
                % -> no other max, min or min max dose constraint
                for k = j+1:numel(cst{i,6})
                    if isequal(cst{i,6}(k).type, 'max dose constraint') || ...
                       isequal(cst{i,6}(k).type, 'min dose constraint') || ...
                       isequal(cst{i,6}(k).type, 'min max dose constraint')
                            error('Simultatenous definition of min, max and or min max dose constraint\n');
                    end
                end
            elseif isequal(cst{i,6}(j).type, 'max mean dose constraint') || ...
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
                    
        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            if isequal(type,'none') || isequal(type,'RBExD') 
                Param = cst{i,6}(j).dose;
            elseif isequal(type,'effect')
                Param = cst{i,5}.alphaX .* cst{i,6}(j).dose + cst{i,5}.betaX .* cst{i,6}(j).dose.^2;
            end
            
            if isequal(cst{i,6}(j).type, 'max dose constraint') 
             
                cl = [cl;-inf];
                cu = [cu;Param];
                
            elseif isequal(cst{i,6}(j).type, 'min dose constraint') 
               
                cl = [cl;Param];
                cu = [cu;inf];
                
            elseif isequal(cst{i,6}(j).type, 'min mean dose constraint') 
    
                cl = [cl;Param];
                cu = [cu;inf];
                
            elseif isequal(cst{i,6}(j).type, 'max mean dose constraint') 
    
                cl = [cl;-inf];
                cu = [cu;Param];
                
            elseif isequal(cst{i,6}(j).type, 'min max mean dose constraint') 
    
                cl = [cl;Param(1)];
                cu = [cu;Param(2)];
                
            elseif isequal(cst{i,6}(j).type, 'min EUD constraint') 
    
                cl = [cl;Param];
                cu = [cu;inf];
                
            elseif isequal(cst{i,6}(j).type, 'max EUD constraint') 
    
                cl = [cl;-inf];
                cu = [cu;Param];
                
            elseif isequal(cst{i,6}(j).type, 'min max EUD constraint') 
    
                cl = [cl;Param(1)];
                cu = [cu;Param(2)];
                
            elseif isequal(cst{i,6}(j).type, 'exact DVH constraint')
                
                cl = [cl;cst{i,6}(j).volume/100];
                cu = [cu;cst{i,6}(j).volume/100];
                
            elseif isequal(cst{i,6}(j).type, 'max DVH constraint') 
                
                cl = [cl;-inf];
                cu = [cu;cst{i,6}(j).volume/100];
                
                % alternative constraint calculation 1/4 %                
                % cl = [cl;-inf];
                % cu = [cu;0];
                % alternative constraint calculation 1/4 %

            elseif isequal(cst{i,6}(j).type, 'min DVH constraint') 

                cl = [cl;cst{i,6}(j).volume/100];
                cu = [cu;inf];

                % alternative constraint calculation 2/4 %                
                % cl = [cl;-inf];
                % cu = [cu;0];
                % alternative constraint calculation 2/4 %
                
            end
      
        end
        
    end
end
   
end