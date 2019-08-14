function [cl,cu] = matRad_getConstraintBounds(optiProb,cst)
% matRad IPOPT get constraint bounds wrapper function
% 
% call
%   [cl,cu] = matRad_getConstraintBounds(cst,options)
%
% input
%   cst:            matRad cst struct
%
% output
%   cl: lower bounds on constraints
%   cu: lower bounds on constraints
%
% References
%   -
%
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


BPtype = class(optiProb.BP);
isEffectBP = strcmp(BPtype,'matRad_EffectProjection');

% Initialize bounds
cl = [];
cu = [];

% compute objective function for every VOI.
for  i = 1:size(cst,1)

    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            optiFunc = cst{i,6}{j};
            
            % only perform computations for constraints
            %if ~isempty(strfind(cst{i,6}{j}.type,'constraint'))
            if isa(optiFunc,'DoseConstraints.matRad_DoseConstraint')
                
                
                if isEffectBP
                    doses = optiFunc.getDoseParameters();
                
                    effect = cst{i,5}.alphaX*doses + cst{i,5}.betaX*doses.^2;
                    
                    optiFunc = optiFunc.setDoseParameters(effect);
                end

                    
                 cl = [cl;optiFunc.lowerBounds(numel(cst{i,4}{1}))];
                 cu = [cu;optiFunc.upperBounds(numel(cst{i,4}{1}))];
                    
                %end
            end

        end % over all objectives of structure

    end % if structure not empty and target or oar

end % over all structures
   

