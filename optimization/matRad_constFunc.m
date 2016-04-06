function c = matRad_constFunc(w,dij,cst,type)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: constraint function for inverse planning supporting max dose
% constraint, min dose constraint, min max dose constraint, min mean, max
% min, min max mean constraint, min EUD constraint, max EUDconstraint, 
% min max EUD constraint, exact DVH constraint, max DVH constraint, 
% min DVH constraint 
% 
% call
%   c = matRad_constFunc(w,dij,cst,type)
%
% input
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   cst:  matRad cst struct
%   type: type of optimizaiton; either 'none','effect' or 'RBExD'
%
% output
%   c: value of constraints
%
% Reference
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

% Initializes constraints
c = [];

% loop over all scenarios
for i = 1:dij.numOfScenarios

    % compute objective function for every VOI.
    for  j = 1:size(cst,1)

        % Only take OAR or target VOI.
        if ~isempty(cst{j,4}) && ( isequal(cst{j,3},'OAR') || isequal(cst{j,3},'TARGET') )

            % get dose vector in current VOI
            d_i = d{i}(cst{j,4});

            % loop over the number of constraints for the current VOI
            for k = 1:numel(cst{j,6})
                
                % only in the nominal case or for robust optimization
                if i == 1 || strcmp(cst{j,6}(k).robustness,'probabilistic') || ...
                             strcmp(cst{j,6}(k).robustness,'voxel-wise worst case')


                    if isequal(cst{i,6}(j).type, 'max dose constraint')

                        epsilon = 1e-3;
                        d_i_max = max(d_i);
                
                        c = [c;d_i_max + epsilon * log( sum(exp((d_i - d_i_max)/epsilon)) )];
            
                    elseif isequal(cst{i,6}(j).type, 'min dose constraint')

                        epsilon = 1e-3;
                        d_i_min = min(d_i);

                        c = [c;d_i_min - epsilon * log( sum(exp((d_i_min - d_i)/epsilon)) )];

                    elseif isequal(cst{j,6}(k).type, 'min mean dose constraint') || ...
                           isequal(cst{j,6}(k).type, 'max mean dose constraint') || ...
                           isequal(cst{j,6}(k).type, 'min max mean dose constraint')

                        c = [c;mean(d_i)];

                    elseif isequal(cst{j,6}(k).type, 'min EUD constraint') || ...
                           isequal(cst{j,6}(k).type, 'max EUD constraint') || ...
                           isequal(cst{j,6}(k).type, 'min max EUD constraint')

                        exponent = cst{j,6}(k).EUD;

                        c = [c;mean(d_i.^exponent)^(1/exponent)];

                    elseif isequal(cst{j,6}(k).type, 'exact DVH constraint') || ...
                           isequal(cst{j,6}(k).type, 'max DVH constraint') || ... 
                           isequal(cst{j,6}(k).type, 'min DVH constraint')

                        % reference dose/effect/RBExDose
                        if isequal(type,'effect')
                            d_ref = dij.ax(cst{j,4}).*cst{j,6}(k).dose + dij.bx(cst{j,4})*cst{j,6}(k).dose^2;
                        else
                            d_ref = cst{j,6}(k).dose;
                        end

                        c = [c;sum(d_i >= d_ref)/size(cst{j,4},1)];

                        % alternative constraint calculation 3/4 %
                        % % get reference Volume
                        % refVol = cst{i,6}(j).volume/100;
                        % 
                        % % calc deviation
                        % deviation = d_i - d_ref;
                        % 
                        % % calc d_ref2: V(d_ref2) = refVol
                        % d_ref2 = matRad_calcInversDVH(refVol,d_i);
                        % 
                        % % apply lower and upper dose limits
                        % if isequal(cst{i,6}(j).type, 'max DVH constraint')
                        %    deviation(d_i < d_ref | d_i > d_ref2) = 0;
                        % elseif isequal(cst{i,6}(j).type, 'min DVH constraint')
                        %    deviation(d_i > d_ref | d_i < d_ref2) = 0;
                        % end
                        % 
                        % %c = [c;sum(deviation)];                              % linear deviation
                        % %c = [c;deviation'*deviation];                        % square devioation
                        % c = [c;(1/size(cst{i,4},1))*(deviation'*deviation)]; % square deviation with normalization
                        % %c = [c;(deviation).^2'*(deviation).^2];               % squared square devioation
                        % alternative constraint calculation 3/4 %
                    
                    end % switch for different constraints
                    
                end % if we are in the nominal sceario or rob opt
                
            end % over all defined constraints & objectives
            
        end % if structure not empty and oar or target
        
    end % over all structures
    
end % over all scenarios

