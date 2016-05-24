function f = matRad_objFuncWrapper(w,dij,cst,type)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT objective function wrapper
% 
% call
%   f = matRad_objFuncWrapper(w,dij,cst,type)
%
% input
%   w:    beamlet/ pencil beam weight vector
%   dij:  matRad dose influence struct
%   cst:  matRad cst struct
%   type: switch to indicate biological optimization
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

global matRad_voxelWeighting;
global matRad_iteration;

% initialize voxel calc Flag
[matRad_voxelWeighting{:,2}] = deal(true);

% get current dose / effect / RBExDose vector
d = matRad_backProjection(w,dij,type);

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
                    isequal(type,'effect') 

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
                elseif strcmp(cst{i,6}(j).robustness,'voxel-wise worst case')

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
                    
                    d_i = [];
                    
                    % get cst index of VOI that corresponds to VOI ring
                    cstidx = find(strcmp(cst(:,2),cst{i,2}(1:end-4)));
         
                    if dij.numOfScenarios > 1
                        % get dose of VOI that corresponds to VOI ring
                        for k = 1:dij.numOfScenarios
                            d_i{k} = d{k}(cst{cstidx,4}{1});
                        end
                        
                        % calc invers DCH of VOI
                        refQ   = cst{i,6}(j).coverage/100;
                        refVol = cst{i,6}(j).volume/100;
                        d_ref2 = matRad_calcInversDCH(refVol,refQ,d_i,dij.numOfScenarios);
                        
                    else
                        idxNom = 1:dij.numOfVoxels;
                        for k = 1:100%size(cst{i,5}.shift_vox,2)
                            idxShift      = cst{i,5}.shift_vox(2,k) + cst{i,5}.shift_vox(1,k)*dij.dimensions(1) + cst{i,5}.shift_vox(3,k)*dij.dimensions(2)*dij.dimensions(1);
                            idx           = idxNom - idxShift;
                            %idx           = mod(idxNom - idxShift,dij.numOfVoxels + 1);
                            %idx(idx == 0) = 1;
                            idx           = idx(cst{cstidx,4}{1});
                            d_i{k}        = d{1}(idx);
                        end
                        
                        % calc invers DCH of VOI
                        refQ   = cst{i,6}(j).coverage/100;
                        refVol = cst{i,6}(j).volume/100;
                        d_ref2 = matRad_calcInversDCH(refVol,refQ,d_i,100);
                        
                    end
                    
                    % get dose of VOI ring
                    d_i = d{1}(cst{i,4}{1});
                    
                    % get voxel dependetn weigthing
                    if matRad_iteration < 0
                        voxelWeighting = 1;
                        
                    else    
                        if isequal(cst{i,5}.voxelWeightingType,'heurWeighting')
                            matRad_calcVoxelWeighting(i,j,cst,d_i,d_ref,d_ref2)
                            voxelWeighting = matRad_voxelWeighting{i,1};

                        elseif isequal(cst{i,5}.voxelWeightingType,'probWeighting')
                            voxelWeighting = 5*cst{i,5}.voxelProb;    
                        end
                        
                    end
    
                    f = f + matRad_objFunc(d_i,cst{i,6}(j),d_ref,d_ref2,voxelWeighting);

                end
            
            end
       
        end
            
    end
    
end
