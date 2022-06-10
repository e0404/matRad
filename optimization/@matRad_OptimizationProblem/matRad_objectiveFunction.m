function f = matRad_objectiveFunction(optiProb,w,dij,cst)
% matRad IPOPT objective function wrapper
%
% call
%   f = matRad_objectiveFuncWrapper(optiProb,w,dij,cst)
%
% input
%   optiProb: matRad optimization problem
%   w:        beamlet/ pencil beam weight vector
%   dij:      matRad dose influence struct
%   cst:      matRad cst struct
%   scenario: index of dij scenario to consider (optional: default 1)
%
% output
%   f: objective function value
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

matRad_cfg = MatRad_Config.instance();

% get current dose / effect / RBExDose vector
optiProb.BP.compute(dij,w);
d = optiProb.BP.GetResult();

% get probabilistic quantities (nearly no overhead if empty)
[dExp,dOmega] = optiProb.BP.GetResultProb();

% get the used scenarios
useScen  = optiProb.BP.scenarios;
scenProb = optiProb.BP.scenarioProb;

% retrieve matching 4D scenarios
fullScen = cell(ndims(d),1);
[fullScen{:}] = ind2sub(size(d),useScen);
contourScen = fullScen{1};

% initialize f
f = 0;

% required for COWC opt
f_COWC = zeros(numel(useScen),1);

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
        
        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            objective = cst{i,6}{j};
            
            % only perform gradient computations for objectiveectives
            if isa(objective,'DoseObjectives.matRad_DoseObjective')
                
                % rescale dose parameters to biological optimization quantity if required
                objective = optiProb.BP.setBiologicalDosePrescriptions(objective,cst{i,5}.alphaX,cst{i,5}.betaX);
                
                % retrieve the robustness type
                robustness = objective.robustness;
                
                switch robustness
                    case 'none' % if conventional opt: just sum objectives of nominal dose
                        for ixScen = useScen
                            d_i = d{ixScen}(cst{i,4}{useScen(1)});
                            f = f + objective.computeDoseObjectiveFunction(d_i);
                        end
                        
                    case 'STOCH' % if prob opt: sum up expectation value of objectives
                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);
                            
                            d_i = d{ixScen}(cst{i,4}{ixContour});
                            
                            f   = f + scenProb(s) * objective.computeDoseObjectiveFunction(d_i);
                            
                        end
                        
                    case 'PROB' % if prob opt: sum up expectation value of objectives
                        
                        d_i = dExp{1}(cst{i,4}{1});
                        
                        f   = f +  objective.computeDoseObjectiveFunction(d_i);
                        
                        p = objective.penalty/numel(cst{i,4}{1});
                        
                        % only one variance term per VOI
                        if j == 1
                            f = f + p * w' * dOmega{i,1};
                        end
                        
                    case 'VWWC'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                        contourIx = unique(contourScen);
                        if ~isscalar(contourIx)
                            % voxels need to be tracked through the 4D CT,
                            % not yet implemented
                            matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                        end
                        
                        % prepare min/max dose vector
                        if ~exist('d_tmp','var')
                            d_tmp = [d{useScen}];
                        end
                        
                        d_Scen = d_tmp(cst{i,4}{contourIx},:);
                        
                        d_max = max(d_Scen,[],2);
                        d_min = min(d_Scen,[],2);
                        
                        if isequal(cst{i,3},'OAR')
                            d_i = d_max;
                        elseif isequal(cst{i,3},'TARGET')
                            d_i = d_min;
                        end
                        
                        f = f + objective.computeDoseObjectiveFunction(d_i);
                        
                    case 'VWWC_INV'  %inverse voxel-wise conformitiy - consider the maximum and minimum dose in the target and optimize the dose conformity
                        contourIx = unique(contourScen);
                        if ~isscalar(contourIx)
                            % voxels need to be tracked through the 4D CT,
                            % not yet implemented
                            matRad_cfg.dispWarning('4D inverted VWWC optimization is currently not supported');
                        end
                        
                        % prepare min/max dose vector
                        if ~exist('d_tmp','var')
                            d_tmp = [d{useScen}];
                        end
                        
                        d_Scen = d_tmp(cst{i,4}{contourIx},:);
                        d_max = max(d_Scen,[],2);
                        d_min = min(d_Scen,[],2);
                        
                        if isequal(cst{i,3},'OAR')
                            d_i = d_min;
                        elseif isequal(cst{i,3},'TARGET')
                            d_i = d_max;
                        end
                        
                        f = f + objective.computeDoseObjectiveFunction(d_i);
                        
                    case 'COWC'  % composite worst case consideres ovarall the worst objective function value
                        
                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);
                            
                            d_i = d{ixScen}(cst{i,4}{ixContour});
                            
                            f_COWC(s) = f_COWC(s) + objective.computeDoseObjectiveFunction(d_i);
                        end
                        
                    case 'OWC'   % objective-wise worst case considers the worst individual objective function value
                        
                        f_OWC = zeros(numel(useScen),1);
                        
                        for s = 1:numel(useScen)
                            ixScen    = useScen(s);
                            ixContour = contourScen(s);
                            
                            d_i = d{ixScen}(cst{i,4}{ixContour});
                            f_OWC(s) = objective.computeDoseObjectiveFunction(d_i);
                        end
                        
                        % compute the maximum objective function value
                        switch optiProb.useMaxApprox
                            case 'logsumexp'
                                fMax = optiProb.logSumExp(f_OWC);
                            case 'pnorm'
                                fMax = optiProb.pNorm(f_OWC,numel(useScen));
                            case 'none'
                                fMax = max(f_OWC);
                            case 'otherwise'
                                matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
                                fMax = max(f_OWC);
                        end
                        f = f + fMax;
                        
                    otherwise
                        matRad_cfg.dispError('Robustness setting %s not supported!',objective.robustness);
                        
                end  %robustness type
            end  % objective check
        end %objective loop
    end %empty check
end %cst structure loop


%Handling the maximum of the composite worst case part
fMax = max(f_COWC(:));
if fMax > 0
    switch optiProb.useMaxApprox
        case 'logsumexp'
            fMax = optiProb.logSumExp(f_COWC);
        case 'pnorm'
            fMax = optiProb.pNorm(f_COWC,numel(useScen));
        case 'cheapCOWC'
            fScenProb = zeros(size(dij.physicalDose));
            fScenProb(useScen) = scenProb;
            f_COWC_tmp = zeros(size(dij.physicalDose));
            f_COWC_tmp(useScen) = f_COWC;
            fMax = optiProb.cheapCOWC(f_COWC_tmp,fScenProb);
        case 'none'
            fMax = max(f_COWC);
        case 'otherwise'
            matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
            fMax = max(f_COWC);
    end
end
%Sum up max of composite worst case part
f = f + fMax;
