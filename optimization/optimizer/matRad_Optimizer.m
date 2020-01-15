classdef (Abstract) matRad_Optimizer < handle
% matRad_Optimizer This is the superclass for all optimizer to be used
% within matRad
%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties (Abstract)
        options %options struct
        wResult
        resultInfo
    end
    
   
    %These should be abstract methods, however Octave can't parse them. As soon 
    %as Octave is able to do this, they should be made abstract again 
    methods %(Abstract)        
        function obj = optimize(obj,w0,optiProb,dij,cst)
          error('Function needs to be implemented!');
        end
        
        function [msg,statusflag] = GetStatus(obj)
          error('Function needs to be implemented!');
        end
    end
end

