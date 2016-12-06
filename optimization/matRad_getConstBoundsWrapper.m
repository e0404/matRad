function [cl,cu] = matRad_getConstBoundsWrapper(cst,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT get constraint bounds wrapper function
% 
% call
%   [cl,cu] = matRad_getConstBounds(cst,options)
%
% input
%   cst:            matRad cst struct
%   options:        option struct defining the type of optimization
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


% Initialize bounds
cl = [];
cu = [];

% compute objective function for every VOI.
for  i = 1:size(cst,1)

    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            % only perform computations for constraints
            if ~isempty(strfind(cst{i,6}(j).type,'constraint'))

                if isequal(options.bioOpt,'none') || isequal(options.ID,'protons_const_RBExD') ||  isequal(options.bioOpt,'LEMIV_RBExD')
                    param = cst{i,6}(j).dose;
                elseif isequal(options.bioOpt,'LEMIV_effect')
                    param = cst{i,5}.alphaX .* cst{i,6}(j).dose + cst{i,5}.betaX .* cst{i,6}(j).dose.^2;
                end

                if strcmp(cst{i,6}(j).robustness,'none')

                    [clTmp,cuTmp] = matRad_getConstBounds(cst{i,6}(j),param);
                    
                    cl = [cl;clTmp];
                    cu = [cu;cuTmp];
                    
                end
            end

        end % over all objectives of structure

    end % if structure not empty and target or oar

end % over all structures
   

