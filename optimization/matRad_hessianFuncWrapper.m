function hessian = matRad_hessianFuncWrapper(w,sigma,lambda,dij,cst,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: hessian function for inverse planning supporting 
% squared underdosage, squared overdosage, squared deviation, mean dose 
% objectives, EUD objectives, DVH objectives, (exact) max dose constraint, 
% (exact) min dose constraint, min mean dose constraint, max mean dose 
% constraint, min EUD constraint, max EUD constraint, max DVH constraint,
% min DVH constraint 
% 
% call
%   hessian = matRad_hessianFuncWrapper(w,sigma,lambda,dij,cst,options)
%
% input
%   w:      bixel weight vector
%   sigma:  scalar factor on the objective
%   lambda: Lagrange multipliers for the constraints
%   dij:    dose influence matrix
%   cst:    matRad cst struct
%   options: option struct defining the type of optimization
%
% output
%   hessian: hessian of the Langrangian at the current point
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

% get current dose / effect / RBExDose vector
d = matRad_backProjection(w,dij,options);

% % initialize projection matrices and id containers
% DoseProjection          = sparse([]);
% mAlphaDoseProjection    = sparse([]);
% mSqrtBetaDoseProjection = sparse([]);
% voxelID                 = [];
% constraintID            = 0;
% scenID                  = [];
% scenID2                 = [];

% initialize hessian diagonal
hessianDiag = sparse(dij.numOfVoxels,1);
constraintCounter = 0;

% compute hessian diagonal function for every VOI.
for i = 1:size(cst,1)

    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of objective and constraints for the current VOI
        for j = 1:numel(cst{i,6})

            % discriminate between objectives and constraints
            if isempty(strfind(cst{i,6}(j).type,'constraint'))
                
                % compute reference
                if (~isequal(cst{i,6}(j).type, 'mean') && ~isequal(cst{i,6}(j).type, 'EUD')) &&...
                    isequal(options.bioOpt,'LEMIV_effect') 

                    d_ref = cst{i,5}.alphaX*cst{i,6}(j).dose + cst{i,5}.betaX*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
                
                % if conventional opt: just add objectives of nominal dose
                if strcmp(cst{i,6}(j).robustness,'none')
                    
                    d_i = d{1}(cst{i,4}{1});              
                    hessianDiag = hessianDiag + sigma * matRad_hessianFunc(dij,d_i,cst{i,6}(j),cst{i,4}{1},d_ref);                    
                else
                    error('robust exact optimization not supported yet!!!')
                end
                
            else                
                
                % compute reference
                if (~isequal(cst{i,6}(j).type, 'max dose constraint')      && ~isequal(cst{i,6}(j).type, 'min dose constraint')          &&...
                    ~isequal(cst{i,6}(j).type, 'max mean dose constraint') && ~isequal(cst{i,6}(j).type, 'min mean dose constraint') && ...
                    ~isequal(cst{i,6}(j).type, 'min EUD constraint')       && ~isequal(cst{i,6}(j).type, 'max EUD constraint'))           && ...
                    ~isequal(cst{i,6}(j).type, 'max dose constraint (exact)')      && ~isequal(cst{i,6}(j).type, 'min dose constraint (exact)') &&...
                    isequal(options.bioOpt,'LEMIV_effect')
                     
                    d_ref = cst{i,5}.alphaX*cst{i,6}(j).dose + cst{i,5}.betaX*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
                
                % if conventional opt: just add constraints of nominal dose
                if strcmp(cst{i,6}(j).robustness,'none')

                    if isequal(cst{i,6}(j).type, 'max dose constraint (exact)') || isequal(cst{i,6}(j).type, 'min dose constraint (exact)')
                        % skip calculation hessian diagonal, zero by definition
                        constraintCounter = constraintCounter + size(cst{i,4}{1},1);
                    else
                        
                        d_i = d{1}(cst{i,4}{1});
                        constraintCounter = constraintCounter + 1;                        
                        hessianDiag = hessianDiag + lambda(constraintCounter) * matRad_hessianFunc(dij,d_i,cst{i,6}(j),cst{i,4}{1},d_ref);    
                    end                    
                    
                    if isequal(options.bioOpt,'LEMIV_effect')

                        error('biological optimization not supported yet for exact optimization!!!')

%                        mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];
%                        mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
%                                                   sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*jacobVec,dij.numOfVoxels,numel(cst{i,4}{1}))];
%                        voxelID                 = [voxelID ;cst{i,4}{1}];
%                        constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];

                    elseif isequal(options.bioOpt,'LEMIV_RBExD')

                        error('biological optimization not supported yet for exact optimization!!!')
                                        
%                        scaledEffect = (dij.gamma(cst{i,4}{1}) + d_i);
% 
%                        delta = jacobVec./(2*dij.bx(cst{i,4}{1}).*scaledEffect);
% 
%                        mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,delta,dij.numOfVoxels,1)];
%                        mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
%                                                   sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*delta,dij.numOfVoxels,numel(cst{i,4}{1}))];
%                        voxelID                 = [voxelID ;cst{i,4}{1}];
%                        constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];

                    end
                    
                else
                    error('robust exact optimization not supported yet!!!')
                end
            end
        end
    end
end

% construct hessian matrix from diagonal
hessian = matRad_constructHessian(hessianDiag,dij,options);

end
