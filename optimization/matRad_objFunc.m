function f = matRad_objFunc(w,dij,cst,type)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: objective function for inverse planning supporting mean dose
% objectives, EUD objectives, squared overdosage, squared underdosage,
% squared deviation and DVH objectives
% 
% call
%    f = matRad_objFunc(w,dij,cst,type)
%
% input
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   cst:  matRad cst struct
%   type: type of optimizaiton; either 'none','effect' or 'RBExD'
%
% output
%   f: objective function value
%
% References
%   [1] http://www.sciencedirect.com/science/article/pii/S0958394701000577
%   [2] http://www.sciencedirect.com/science/article/pii/S0360301601025858
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

% get current dose / effect / RBExDose vector
d = matRad_backProjection(w,dij,type);

% Initializes f
f = 0;

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
    
        % get dose vector in current VOI
        d_i = d(cst{i,4});
                
        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            % reference dose/effect/RBExDose
             if isempty(strfind(cst{i,6}(j).type,'constraint'))
                if (~isequal(cst{i,6}(j).type, 'mean') && ~isequal(cst{i,6}(j).type, 'EUD')) &&...
                        isequal(type,'effect') 

                    d_ref = dij.ax(cst{i,4}).*cst{i,6}(j).dose + dij.bx(cst{i,4})*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
             end
            
            if isequal(cst{i,6}(j).type, 'square underdosing') 
               
                if ~isequal(cst{i,3},'OAR')
                    % underdose : dose minus prefered dose
                    underdose = d_i - d_ref;

                    % apply positive operator
                    underdose(underdose>0) = 0;

                    % calculate objective function
                    f = f + (cst{i,6}(j).penalty/size(cst{i,4},1))*(underdose'*underdose);

                else
                    disp(['square underdosing constraint for ' cst{i,2} ' will be skipped'])
                end
                
            elseif isequal(cst{i,6}(j).type, 'square overdosing')
                
                    % overdose : dose minus prefered dose
                    overdose = d_i - d_ref;

                    % apply positive operator
                    overdose(overdose<0) = 0;

                    % calculate objective function
                    f = f + (cst{i,6}(j).penalty/size(cst{i,4},1))*(overdose'*overdose);
                
           elseif isequal(cst{i,6}(j).type, 'square deviation')
               
               if ~isequal(cst{i,3},'OAR')
                    % deviation : dose minus prefered dose
                    deviation = d_i - d_ref;

                    % claculate objective function
                    f = f + (cst{i,6}(j).penalty/size(cst{i,4},1))*(deviation'*deviation);

                else
                    disp(['square deviation constraint for ' cst{i,2} ' will be skipped'])
                end
                
            elseif isequal(cst{i,6}(j).type, 'mean')              
                
                if ~isequal(cst{i,3},'TARGET')
                    % calculate objective function
                    f = f + (cst{i,6}(j).penalty/size(cst{i,4},1))*sum(d_i);

                else
                    disp(['mean constraint for ' cst{i,2} ' will be skipped'])
                end
                
            elseif isequal(cst{i,6}(j).type, 'EUD') 
               
               if ~isequal(cst{i,3},'TARGET')
                    % get exponent for EUD
                    exponent = cst{i,6}(j).EUD;

                    % calculate objective function and delta
                    if sum(d_i.^exponent)>0

                        f = f + cst{i,6}(j).penalty * nthroot((1/size(cst{i,4},1)) * sum(d_i.^exponent),exponent);

                    end

               else
                    disp(['EUD constraint for ' cst{i,2} ' will be skipped'])
               end
                
            elseif isequal(cst{i,6}(j).type, 'max DVH objective') ||...
                   isequal(cst{i,6}(j).type, 'min DVH objective')
                
                % get reference Volume
                refVol = cst{i,6}(j).volume/100;
                
                % calc deviation
                deviation = d_i - d_ref;
                
                % calc d_ref2: V(d_ref2) = refVol
                d_ref2 = matRad_calcInversDVH(refVol,d_i);
                
                % apply lower and upper dose limits
                if isequal(cst{i,6}(j).type, 'max DVH objective')
                     deviation(d_i < d_ref | d_i > d_ref2) = 0;
                elseif isequal(cst{i,6}(j).type, 'min DVH objective')
                     deviation(d_i > d_ref | d_i < d_ref2) = 0;
                end

                % claculate objective function
                f = f + (cst{i,6}(j).penalty/size(cst{i,4},1))*(deviation'*deviation);
                
            end
      
        end
        
    end
end
   
end