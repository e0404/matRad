function f = matRad_objFuncWrapper(w,dij,cst,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT objective function wrapper
% 
% call
%   f = matRad_objFuncWrapper(w,dij,cst,type)
%
% input
%   w:       beamlet/ pencil beam weight vector
%   dij:     matRad dose influence struct
%   cst:     matRad cst struct
%   options: option struct defining the type of optimization
%
% output
%   f: objective function value
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

% get current dose / effect / RBExDose vector
[d,d_exp,Omega] = matRad_backProjection(w,dij,cst,options);

% Initialize f
f = 0;

% if composite worst case optimization is used then create a vector for book keeping
for i = 1:size(cst,1)
  for j = 1:numel(cst{i,6})
      if strcmp(cst{i,6}(j).robustness,'COWC')
         f_COWC = zeros(options.numOfScen,1);break;
      end
  end
end

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})

            % only perform gradient computations for objectives
            if isempty(strfind(cst{i,6}(j).type,'constraint'))

                % compute reference
                if (~isequal(cst{i,6}(j).type, 'mean') && ~isequal(cst{i,6}(j).type, 'EUD')) &&...
                    isequal(options.quantityOpt,'effect') 

                    d_ref = cst{i,5}.alphaX*cst{i,6}(j).dose + cst{i,5}.betaX*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
                
               % if conventional opt: just sum objectives of nominal dose
                if strcmp(cst{i,6}(j).robustness,'none')

                    d_i = d{1}(cst{i,4}{1});

                    f = f + matRad_objFunc(d_i,cst{i,6}(j),d_ref);
                    
                
                % if prob opt: sum up expectation value of objectives
                elseif strcmp(cst{i,6}(j).robustness,'PROB')
                    
                        d_i = d_exp{1}(cst{i,4}{1});
                        
                        f   = f +  matRad_objFunc(d_i,cst{i,6}(j),d_ref);
                        
                        % only one variance term per VOI
                        if j == 1
                            f = f + w' * Omega{i};
                        end
                    
                % if voxel-wise worst case or voxel-wise conformitiy (only for target structures)
                elseif strcmp(cst{i,6}(j).robustness,'VWWC') || strcmp(cst{i,6}(j).robustness,'VWWC_CONF')

                    % prepare min/max dose vector
                    if ~exist('d_tmp','var')
                         d_tmp = [d{:}];
                    end
                    
                    d_Scen = d_tmp(cst{i,4}{1},:);
                    d_max = max(d_Scen,[],2);
                    d_min = min(d_Scen,[],2);
                         
                    if isequal(cst{i,3},'OAR')
                        d_i = d_max;
                    elseif isequal(cst{i,3},'TARGET')
                        d_i = d_min;
                    end

                    if strcmp(cst{i,6}(j).robustness,'VWWC')
                        f = f + matRad_objFunc(d_i,cst{i,6}(j),d_ref);
                    elseif strcmp(cst{i,6}(j).robustness,'VWWC_CONF') && isequal(cst{i,6}(j).type, 'square overdosing')
                        f = f + matRad_objFunc(d_max,cst{i,6}(j),d_ref);
                    elseif strcmp(cst{i,6}(j).robustness,'VWWC_CONF') && isequal(cst{i,6}(j).type, 'square underdosing')
                        f = f + matRad_objFunc(d_min,cst{i,6}(j),d_ref);     
                    end
                    
                % composite worst case consideres ovarall the worst objective function value       
                elseif strcmp(cst{i,6}(j).robustness,'COWC')
                   
                     for ixScen = 1:options.numOfScen

                        d_i = d{ixScen}(cst{i,4}{1});
         
                        f_COWC(ixScen) = f_COWC(ixScen) + matRad_objFunc(d_i,cst{i,6}(j),d_ref);
                        
                     end   
            
                % objective-wise worst case consideres the worst individual objective function value        
                elseif strcmp(cst{i,6}(j).robustness,'OWC')
                    
                     matRad_dispToConsole(['not yet implemented \n'],param,'error');
                 
                end
       
            end
            
        end
    
    end
end

% extract the worst total objective function value across all scenarios
if exist('f_COWC','var')
   f = f + max(f_COWC);
end
