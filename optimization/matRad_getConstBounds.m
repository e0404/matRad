function [cl,cu] = matRad_getConstBounds(constraint,param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT get constraint bounds function
% 
% call
%   [cl,cu] = matRad_getConstBounds(constraint,param)
%
% input
%   constraint: matRad constraint struct
%   param:      reference parameter
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
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isequal(constraint.type, 'max dose constraint') 

    cl = -inf;
    cu = param;

elseif isequal(constraint.type, 'min dose constraint') 

    cl = param;
    cu = inf;

elseif isequal(constraint.type, 'min mean dose constraint') 

    cl = param;
    cu = inf;

elseif isequal(constraint.type, 'max mean dose constraint') 

    cl = -inf;
    cu = param;

elseif isequal(constraint.type, 'min EUD constraint') 

    cl = param;
    cu = inf;

elseif isequal(constraint.type, 'max EUD constraint') 

    cl = -inf;
    cu = param;

elseif isequal(constraint.type, 'max DVH constraint') 

    cl = -inf;
    cu = constraint.volume/100;

    % alternative constraint calculation 1/4 %                
    % cl = [cl;-inf];
    % cu = [cu;0];
    % alternative constraint calculation 1/4 %

elseif isequal(constraint.type, 'min DVH constraint') 

    cl = constraint.volume/100;
    cu = inf;

    % alternative constraint calculation 2/4 %                
    % cl = [cl;-inf];
    % cu = [cu;0];
    % alternative constraint calculation 2/4 %
        
end % constraint switch
