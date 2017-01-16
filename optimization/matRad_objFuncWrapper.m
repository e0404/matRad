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

global fScaling

% get current dose / effect / RBExDose vector
d = matRad_backProjection(w,dij,options);

% Initialize f
f = 0;

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
                    isequal(options.bioOpt,'LEMIV_effect') 

                    d_ref = cst{i,5}.alphaX*cst{i,6}(j).dose + cst{i,5}.betaX*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
                
               % if conventional opt: just sum objectives of nominal dose
                if strcmp(cst{i,6}(j).robustness,'none')

                    d_i = d{1}(cst{i,4}{1});

                    f = f + matRad_objFunc(d_i,cst{i,6}(j),d_ref);

                % if prob opt: sum up expectation value of objectives
                elseif strcmp(cst{i,6}(j).robustness,'probabilistic')

                    for k = 1:dij.numOfScenarios

                        d_i = d{k}(cst{i,4}{1});

                        f = f + dij.probOfScenarios(k) * matRad_objFunc(d_i,cst{i,6}(j),d_ref);

                    end

                % if voxel-wise worst case: sum up objective of min/max dose
                elseif strcmp(cst{i,6}(j).robustness,'WC')

                    % prepare min/max dose vector we have chosen voxel-wise worst case
                    if ~exist('d_max','var')
                        d_max = max([d{:}],[],2);
                        d_min = min([d{:}],[],2);
                    end

                    if isequal(cst{i,3},'OAR')
                        d_i = d_max(cst{i,4}{1});
                    elseif isequal(cst{i,3},'TARGET')
                        d_i = d_min(cst{i,4}{1});
                    end

                    f = f + matRad_objFunc(d_i,cst{i,6}(j),d_ref);
                    
                % if coveraged based opt    
                elseif strcmp(cst{i,6}(j).robustness,'coverage')
                    
                    % calc invers DCH
                    Q_ref  = cst{i,6}(j).coverage/100;
                    V_ref  = cst{i,6}(j).volume/100;
                    d_ref2 = matRad_calcInversDCH(V_ref,Q_ref,d,dij,cst(i,:)); 
                    
                    if dij.numOfScenarios > 1
                        
                        for k = 1:dij.numOfScenarios
                            
                            % get VOI dose in current scenario
                            d_i = d{k}(cst{i,4}{1});

                            % get voxel dependent weigthing
                            voxelWeighting = 1; 
                            
                            % calculate dose deviations from d_ref
                            fTmp(k) = matRad_objFunc(d_i,cst{i,6}(j),d_ref,d_ref2,voxelWeighting);
                            
                        end
                        
                        % claculate objective function
                        f = f + sum(dij.ScenProb.*fTmp);
                        
                    else
                        
                        % get VOI ScenUnion dose of nominal scneario
                        cstLogical = strcmp(cst(:,2),[cst{i,2},' ScenUnion']);
                        d_i        = d{1}(cst{cstLogical,5}.voxelID);
                        
                        % get voxel dependent weigthing
                        voxelWeighting = 5*cst{cstLogical,5}.voxelProb;
                        
                        % calculate objective function
                        f = f + matRad_objFunc(d_i,cst{i,6}(j),d_ref,d_ref2,voxelWeighting);
                        
                    end
                    
                end
            
            end
       
        end
            
    end
    
end

% apply objective scaling
f = fScaling*f;
