function c = matRad_constFuncWrapper(w,dij,cst,type)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT constraint function wrapper
% 
% call
%   c = matRad_constFuncWrapper(w,dij,cst,type)
%
% input
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   cst:  matRad cst struct
%   type: type of optimizaiton; either 'none','effect' or 'RBExD'
%
% output
%   c: constraint function value
%
% Reference
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

global cScaling
global CONSTRAINT
global matRad_iteration

% get current dose / effect / RBExDose vector
d = matRad_backProjection(w,dij,type);

% Initializes constraints
c = [];

% compute objective function for every VOI.
for  i = 1:size(cst,1)

    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            % only perform computations for constraints
            if ~isempty(strfind(cst{i,6}(j).type,'constraint'))
                
                % compute reference
                if (~isequal(cst{i,6}(j).type, 'max dose constraint') && ~isequal(cst{i,6}(j).type, 'min dose constraint') &&...
                    ~isequal(cst{i,6}(j).type, 'min mean dose constraint') && ~isequal(cst{i,6}(j).type, 'max mean dose constraint') &&...
                    ~isequal(cst{i,6}(j).type, 'min max mean dose constraint') && ~isequal(cst{i,6}(j).type, 'min EUD constraint') &&...
                    ~isequal(cst{i,6}(j).type, 'max EUD constraint') && ~isequal(cst{i,6}(j).type, 'min max EUD constraint')) &&...
                    isequal(type,'effect')
                     
                    d_ref = cst{i,5}.alphaX*cst{i,6}(j).dose + cst{i,5}.betaX*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end

                % if conventional opt: just add constraints of nominal dose
                if strcmp(cst{i,6}(j).robustness,'none')

                    d_i = d{1}(cst{i,4}{1});

                    c = [c; matRad_constFunc(d_i,cst{i,6}(j),d_ref)];

                % if prob opt or voxel-wise worst case: add constraints of all dose scenarios
                elseif strcmp(cst{i,6}(j).robustness,'probabilistic') || strcmp(cst{i,6}(j).robustness,'voxel-wise worst case')
                    
                    for k = 1:dij.numOfScenarios
                        
                        d_i = d{k}(cst{i,4}{1});
                        
                        c = [c; matRad_constFunc(d_i,cst{i,6}(j),d_ref)];
                        
                    end
                    
                % if coveraged based opt   
                elseif strcmp(cst{i,6}(j).robustness,'coverage')
                    
                    if isequal(cst{i,6}(j).type, 'max DCH Area constraint') || ...
                       isequal(cst{i,6}(j).type, 'min DCH Area constraint')
                   
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
                                cTmp(k) = matRad_constFunc(d_i,cst{i,6}(j),d_ref,d_ref2,voxelWeighting);

                            end

                            % claculate constraint function
                            c = [c; sum(dij.ScenProb.*cTmp)];

                        else

                            % get VOI ScenUnion dose of nominal scneario
                            cstLogical = strcmp(cst(:,2),[cst{i,2},' ScenUnion']);
                            d_i        = d{1}(cst{cstLogical,5}.voxelID);

                            % get voxel dependent weigthing
                            voxelWeighting = 1; 

                            % claculate constraint function
                            c = [c; matRad_constFunc(d_i,cst{i,6}(j),d_ref,d_ref2,voxelWeighting)];
                        end

                    elseif isequal(cst{i,6}(j).type, 'max DCH Theta constraint') || ...
                           isequal(cst{i,6}(j).type, 'min DCH Theta constraint')
                       
                        if dij.numOfScenarios > 1
                            
                            for k = 1:dij.numOfScenarios

                                % get VOI dose in current scenario
                                d_i = d{k}(cst{i,4}{1});

                                % calculate volumes
                                volume_pi(k) = matRad_constFunc(d_i,cst{i,6}(j),d_ref);
                                
                            end
                            
                            % get scenario probabilities
                            scenProb = dij.ScenProb;
                            
                        else
                            
                            for k = 1:cst{i,5}.VOIShift.ncase
                                
                                % get VOI dose in current scenario
                                if isequal(cst{i,5}.VOIShift.shiftType,'rounded')
                                    d_i = d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.roundedShift.idxShift(k));
                                elseif isequal(cst{i,5}.VOIShift.shiftType,'linInterp')
                                    error('linInterp in DCH Theta constraint not implemented yet')
                                end                                
                                    
                                % calculate volumes
                                volume_pi(k) = matRad_constFunc(d_i,cst{i,6}(j),d_ref);

                            end
                            
                            % get scenario probabilities
                            scenProb = 1/cst{i,5}.VOIShift.ncase;  % assume equiprobable scenarios
                            
                        end                        
                        
                        % calculate constraint function
                        c = [c; sum(scenProb.*(volume_pi >= cst{i,6}(j).volume/100))];                                            
                       
                    end % if Area or Theta constrained is used

                end % if we are in the nominal sceario, rob opt or COP
            
            end % if type is constraint

        end % over all defined constraints & objectives

    end % if structure not empty and oar or target

end % over all structures

% save unscaled constraint
CONSTRAINT(:,matRad_iteration+1) = c;

% apply constraint scaling
c = cScaling.*c;
 
