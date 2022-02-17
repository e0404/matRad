function c = matRad_constraintFunctions(optiProb,w,dij,cst)
% matRad IPOPT callback: constraint function for inverse planning
% supporting max dose constraint, min dose constraint, min mean dose constraint,
% max mean dose constraint, min EUD constraint, max EUD constraint,
% max DVH constraint, min DVH constraint
%
% call
%   c = matRad_constraintFunctions(optiProb,w,dij,cst)
%
% input
%   optiProb:   option struct defining the type of optimization
%   w:          bixel weight vector
%   dij:        dose influence matrix
%   cst:        matRad cst struct
%
% output
%   c:          value of constraints
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


% get current dose / effect / RBExDose vector
%d = matRad_backProjection(w,dij,optiProb);
optiProb.BP = optiProb.BP.compute(dij,w);
d = optiProb.BP.GetResult();

if ~isempty(optiProb.BP_LET)
    optiProb.BP_LET = optiProb.BP_LET.compute(dij,w);
    LET = optiProb.BP_LET.GetResult();
end

if ~isempty(optiProb.BP_DADRfixed)
    optiProb.BP_DADRfixed = optiProb.BP_DADRfixed.compute(dij,w);
    DADR = optiProb.BP_DADRfixed.GetResult();
end

% Initializes constraints
c = [];

% compute objective function for every VOI.
for  i = 1:size(cst,1)

    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})

            constraint = cst{i,6}{j};

            % only perform computations for constraints
            % if ~isempty(strfind(obj.type,'constraint'))
            if isa(constraint,'DoseConstraints.matRad_DoseConstraint')

                % if we have effect optimization, temporarily replace doses with effect
                % Maybe we should put some switch into the classes for that
                if (~isequal(constraint.name, 'max dose constraint')      && ~isequal(constraint.name, 'min dose constraint')      &&...
                        ~isequal(constraint.name, 'max mean dose constraint') && ~isequal(constraint.name, 'min mean dose constraint') && ...
                        ~isequal(constraint.name, 'min EUD constraint')       && ~isequal(constraint.name, 'max EUD constraint'))     && ...
                        (isa(optiProb.BP,'matRad_EffectProjection') && ~isa(optiProb.BP,'matRad_VariableRBEProjection'))

                    doses = constraint.getDoseParameters();

                    effect = cst{i,5}.alphaX*doses + cst{i,5}.betaX*doses.^2;

                    constraint = constraint.setDoseParameters(effect);
                end

                d_i = d{1}(cst{i,4}{1});
                c = [c; constraint.computeDoseConstraintFunction(d_i)];


            end % if we are a constraint

            if isa(constraint,'DADRConstraints.matRad_DADRConstraint')
                % if conventional opt: just sum objectiveectives of nominal LET
                %if strcmp(cst{i,6}{j}.robustness,'none')

                DADR_i = DADR{1}(cst{i,4}{1});
                c = [c; constraint.computeDADRConstraintFunction(DADR_i)];

                %end
            end

        end % over all defined constraints & objectives

    end % if structure not empty and oar or target

end % over all structures
