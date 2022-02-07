function jacob = matRad_constraintJacobian(optiProb,wCombined,dij,cst)
% matRad IPOPT callback: jacobian function for inverse planning
% supporting max dose constraint, min dose constraint, min mean dose constraint,
% max mean dose constraint, min EUD constraint, max EUD constraint, max DVH
% constraint, min DVH constraint
%
% call
%   jacob = matRad_jacobFunc(optiProb,w,dij,cst)
%
% input
%   optiProb: option struct defining the type of optimization
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   cst:  matRad cst struct
%
% output
%   jacob: jacobian of constraint function
%
% References
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
%d = matRad_backProjection(w,dij,optiProb);
%d = optiProb.matRad_backProjection(w,dij);
optiProb.BP = optiProb.BP.compute(dij,wCombined);
dTilda = optiProb.BP.GetResult();
d{1} = dTilda{1}(1:dij.doseGrid.numOfVoxels);
DADR{1} = dTilda{1}(1+dij.doseGrid.numOfVoxels:end);
indexPB = size(dij.physicalDose{1},2);
w = wCombined(1:indexPB);
I = wCombined(indexPB+1:end);

% initialize jacobian (only single scenario supported in optimization)
jacob = sparse([]);

% initialize projection matrices and id containers
DoseProjection{1}          = sparse([]);
DADRProjection{1}          = sparse([]);
mAlphaDoseProjection{1}    = sparse([]);
mSqrtBetaDoseProjection{1} = sparse([]);
voxelID                     = [];
constraintID                = [];
scenID2                     = [];

% compute objective function for every VOIfor i = 1:size(cst,1)
for i = 1:size(cst,1)

    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') ||  isequal(cst{i,3},'TARGET') )


        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})

            obj = cst{i,6}{j}; %Get the Optimization Object

            % only perform computations for constraints
            %if ~isempty(strfind(obj.type,'constraint'))
            if isa(obj,'DoseConstraints.matRad_DoseConstraint')

                % compute reference
                if (~isequal(obj.name, 'max dose constraint')      && ~isequal(obj.name, 'min dose constraint')      &&...
                        ~isequal(obj.name, 'max mean dose constraint') && ~isequal(obj.name, 'min mean dose constraint') && ...
                        ~isequal(obj.name, 'min EUD constraint')       && ~isequal(obj.name, 'max EUD constraint'))      && ...
                        (isa(optiProb.BP,'matRad_EffectProjection') && ~isa(optiProb.BP,'matRad_VariableRBEProjection'))

                    doses = obj.getDoseParameters();

                    effect = cst{i,5}.alphaX*doses + cst{i,5}.betaX*doses.^2;

                    obj = obj.setDoseParameters(effect);
                end

                % if conventional opt: just add constraints of nominal dose
                %if strcmp(cst{i,6}{j}.robustness,'none')

                d_i = d{1}(cst{i,4}{1});

                %jacobVec =  matRad_jacobFunc(d_i,cst{i,6}{j},d_ref);
                jacobSub = obj.computeDoseConstraintJacobian(d_i);

                nConst = size(jacobSub,2);

                %tic;

                %Iterate through columns of the sub-jacobian
                %TODO: Maybe this could all be function of the projection
                %Objects???
                if isa(optiProb.BP,'matRad_DADRProjection') % The only option we need ?
                    startIx = size(DoseProjection{1},2) + 1;
                    %First append the Projection matrix with sparse zeros
                    DoseProjection{1}          = [DoseProjection{1},sparse(dij.doseGrid.numOfVoxels,nConst)];

                    %Now directly write the jacobian in there
                    DoseProjection{1}(cst{i,4}{1},startIx:end) = jacobSub;
                elseif isa(optiProb.BP,'matRad_DoseProjection') && ~isempty(jacobSub) || isa(optiProb.BP,'matRad_ConstantRBEProjection')

                    startIx = size(DoseProjection{1},2) + 1;
                    %First append the Projection matrix with sparse zeros
                    DoseProjection{1}          = [DoseProjection{1},sparse(dij.doseGrid.numOfVoxels,nConst)];

                    %Now directly write the jacobian in there
                    DoseProjection{1}(cst{i,4}{1},startIx:end) = jacobSub;

                elseif isa(optiProb.BP,'matRad_EffectProjection') && ~isempty(jacobSub)

                    if isa(optiProb.BP,'matRad_VariableRBEProjection')
                        scaledEffect = (dij.gamma(cst{i,4}{1}) + d_i);
                        jacobSub = jacobSub./(2*dij.bx(cst{i,4}{1}) .* scaledEffect);
                    end

                    startIx = size(mAlphaDoseProjection{1},2) + 1;

                    %First append the alphaDose matrix with sparse
                    %zeros then insert
                    mAlphaDoseProjection{1}    = [mAlphaDoseProjection{1},sparse(dij.doseGrid.numOfVoxels,nConst)];
                    mAlphaDoseProjection{1}(cst{i,4}{1},startIx:end) = jacobSub;

                    %The betadose has a different structure due to the
                    %quadratic transformation, but in principle the
                    %same as above
                    mSqrtBetaDoseProjection{1} =  [mSqrtBetaDoseProjection{1}, sparse(repmat(cst{i,4}{1},nConst,1),repmat(1:numel(cst{i,4}{1}),1,nConst),2*reshape(jacobSub',[],1),dij.doseGrid.numOfVoxels,nConst*numel(cst{i,4}{1}))];

                    if isempty(constraintID)
                        newID = 1;
                    else
                        newID = constraintID(end)+1;
                    end

                    voxelID = [voxelID;repmat(cst{i,4}{1},nConst,1)];                         %Keep track of voxels for organizing the sqrt(beta)Dose projection later
                    constraintID = [constraintID, ...
                        reshape(ones(numel(cst{i,4}{1}),1)*[newID:newID+nConst-1],[1 nConst*numel(cst{i,4}{1})])];  %Keep track of constraints for organizing the sqrt(beta)Dose projection later

                end


            elseif isa(obj,'DADRConstraints.matRad_DADRConstraint')

                d_i = DADR{1}(cst{i,4}{1});

                %jacobVec =  matRad_jacobFunc(d_i,cst{i,6}{j},d_ref);
                jacobSub = obj.computeDoseConstraintJacobian(d_i);

                nConst = size(jacobSub,2);
                % Assuming we only have isa(optiProb.BP,'matRad_FLASHProjection')
                startIx = size(DADRProjection{1},2) + 1;
                %First append the Projection matrix with sparse zeros
                DADRProjection{1}          = [DADRProjection{1},sparse(dij.doseGrid.numOfVoxels,nConst)];

                %Now directly write the jacobian in there
                DADRProjection{1}(cst{i,4}{1},startIx:end) = jacobSub;

            end

        end

    end

end

% enter if statement also for protons using a constant RBE
if isa(optiProb.BP,'matRad_DoseProjection')% Eventually changed

    if ~isempty(DoseProjection{1})
        jacob = DoseProjection{1}' * dij.physicalDose{1};
    end

elseif isa(optiProb.BP,'matRad_DADRProjection')
    if ~isempty(DoseProjection{1})
        jacobDose = DoseProjection{1}' * dij.physicalDose{1};
        jacob = [jacob; [jacobDose zeros(size(jacobDose))]];
    end

    if ~isempty(DADRProjection{1})
        d = dij.physicalDose{1}*w;
        wFinal=length(w);
        % Computation of the DADR
        num = zeros(size(d));
        DADR = zeros([size(d,1) 1]);
        for j =1:wFinal
            num = num+ dij.physicalDose{1}(:,j).*w(j).*dij.physicalDose{1}(:,j)*I(j)*10^4;
        end
        DADR(dij.ixDose,1) = num(dij.ixDose)./d(dij.ixDose);

        % Computation of dGamma/dI
        tempVar= 10^4*dij.physicalDose{1}.^2.*w';
        IGrad_helper=zeros([length(d) length(w)]);
        IGrad_helper(dij.ixDose,:) = (tempVar(dij.ixDose,:)./d(dij.ixDose));

        jacobDADRI = sparse(DADRProjection{1}' * IGrad_helper);


        % Computation of dGamma/dw
        tempVar2= 10^4*dij.physicalDose{1}.^2.*I';
        IGrad_helper2=zeros([length(d) length(w)]);
        IGrad_helper2(dij.ixDose,:) = (tempVar2(dij.ixDose,:)./d(dij.ixDose));

        jacobDADRw = DADRProjection{1}' * IGrad_helper2;

        IGrad_helper3=zeros([length(d) length(w)]);
        IGrad_helper3(dij.ixDose,:) = ...
            (dij.physicalDose{1}(dij.ixDose,:)./d(dij.ixDose)).*DADR(dij.ixDose);


        jacobDADRw = sparse(jacobDADRw-DADRProjection{1}' * IGrad_helper3);

        jacob = [jacob; [jacobDADRw jacobDADRI ]];



        %         index = find(d);
        %         IGrad_helper = zeros(size(dij.physicalDose{1}));
        %         IGrad_helper(index,:)= (d(index)./d(index).^2).*(dij.physicalDose{1}(index,:).^2.*w');
        %
        %         jacobDADR = sparse(DADRProjection{1}' * IGrad_helper);
        %         jacob = [jacob; [zeros(1,indexPB) jacobDADR ]];
        %
        tempVar= dij.physicalDose{1}.^2.*w';
        IGrad_helper = zeros([length(d) length(w)]);
        IGrad_helper(dij.ixDose,:) = (tempVar(dij.ixDose,:)./d(dij.ixDose));
        %         grad = DADRProjection{1}';



        %         IGrad_helper = (1./d).*(10^6*(dij.physicalDose{1}.^2.*w'));
        %         [C,RT] = rmmissing(0*IGrad_helper );
        %          clean = (1-RT);
        %          B = IGrad_helper(clean>0,:);
        %          grad = DADRProjection{1}';
        %        jacobDADR = sparse((grad(:,clean>0) * B));


        %jacob = [jacob; [zeros(size(jacobDADR)) jacobDADR ]];


    end

elseif isa(optiProb.BP,'matRad_ConstantRBEProjection')

    if ~isempty(DoseProjection{1})
        jacob = DoseProjection{1}' * dij.RBE * dij.physicalDose{1};
    end

elseif isa(optiProb.BP,'matRad_EffectProjection')

    if ~isempty(mSqrtBetaDoseProjection{1}) && ~isempty(mAlphaDoseProjection{1})
        mSqrtBetaDoseProjection{1} = mSqrtBetaDoseProjection{1}' * dij.mSqrtBetaDose{1} * w;
        mSqrtBetaDoseProjection{1} = sparse(voxelID,constraintID,mSqrtBetaDoseProjection{1},...
            size(mAlphaDoseProjection{1},1),size(mAlphaDoseProjection{1},2));

        jacob   = mAlphaDoseProjection{1}' * dij.mAlphaDose{1} +...
            mSqrtBetaDoseProjection{1}' * dij.mSqrtBetaDose{1};

    end
end

end