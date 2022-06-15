function resultGUI = matRad_calcCubes(w,dij,scenNum)
% matRad computation of all cubes for the resultGUI struct
% which is used as result container and for visualization in matRad's GUI
%
% call
%   resultGUI = matRad_calcCubes(w,dij)
%   resultGUI = matRad_calcCubes(w,dij,scenNum)
%
% input
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   scenNum: optional: number of scenario to calculated (default 1)
%
% output
%   resultGUI: matRad result struct
%
% References
%   -
%
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

if nargin < 3
    scenNum = 1;
end

resultGUI.w = w;

% get bixel - beam correspondence
for i = 1:dij.numOfBeams
    beamInfo(i).suffix = ['_beam', num2str(i)];
    beamInfo(i).logIx  = (dij.beamNum == i);
end
beamInfo(dij.numOfBeams+1).suffix = '';
beamInfo(dij.numOfBeams+1).logIx  = true(size(resultGUI.w,1),1);


%% Physical Dose
doseFields = {'physicalDose','doseToWater'};
doseQuantities = {'','_std','_batchStd'};
% compute physical dose for all beams individually and together
for j = 1:length(doseFields)
    for k = 1:length(doseQuantities)
        % Check if combination is a field in dij, otherwise skip
        if isfield(dij,[doseFields{j} doseQuantities{k}])
            % Handle standard deviation fields and add quadratically
            if contains(doseQuantities,'std','IgnoreCase',true)
                for i = 1:length(beamInfo)
                    resultGUI.([doseFields{j}, doseQuantities{k}, beamInfo(i).suffix]) = sqrt(reshape(full(dij.([doseFields{j} doseQuantities{k}]){scenNum}.^2 * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions));
                    resultGUI.([doseFields{j}, doseQuantities{k}, beamInfo(i).suffix])(isnan(resultGUI.([doseFields{j}, doseQuantities{k}, beamInfo(i).suffix]))) = 0;
                end
            % Handle normal fields as usual
            else
                for i = 1:length(beamInfo)
                    resultGUI.([doseFields{j}, doseQuantities{k}, beamInfo(i).suffix]) = reshape(full(dij.([doseFields{j} doseQuantities{k}]){scenNum} * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions);
                end
            end
        end
    end
end


%% LET
% consider LET
if isfield(dij,'mLETDose')
    for i = 1:length(beamInfo)
        LETDoseCube                                 = reshape(full(dij.mLETDose{scenNum} * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions);
        resultGUI.(['LET', beamInfo(i).suffix])     = zeros(dij.doseGrid.dimensions);
        ix                                          = resultGUI.(['physicalDose', beamInfo(i).suffix]) > 0;
        resultGUI.(['LET', beamInfo(i).suffix])(ix) = LETDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix);
    end
end


%% RBE weighted dose
% consider RBE for protons and skip varRBE calculation
if isfield(dij,'RBE') && isscalar(dij.RBE)
    for i = 1:length(beamInfo)
        resultGUI.(['RBExD', beamInfo(i).suffix]) = resultGUI.(['physicalDose', beamInfo(i).suffix]) * dij.RBE;
    end
elseif any(contains(fieldnames(dij),'alpha','IgnoreCase',true))
    % Load RBE models if MonteCarlo was calculated for multiple models
    if isfield(dij,'RBE_models')
        RBE_model = cell(1,length(dij.RBE_models)+1);
        for i = 1:length(dij.RBE_models)
            RBE_model{i+1} = ['_' dij.RBE_models{i}];
        end
    else
        RBE_model = {''};
    end

    % consider biological optimization
    wBeam = (resultGUI.w .* beamInfo(i).logIx);
    ix = dij.bx(:,scenNum)~=0 & resultGUI.(['physicalDose', beamInfo(i).suffix])(:) > 0;

    % Loop through RBE models
    for j = 1:length(RBE_model)
        % Check if combination is a field in dij, otherwise skip
        if isfield(dij,['mAlphaDose' RBE_model{j}])
            for i = 1:length(beamInfo)
                % Calculate effect from alpha- and sqrtBetaDose
                resultGUI.(['effect', RBE_model{j}, beamInfo(i).suffix])       = full(dij.(['mAlphaDose' RBE_model{j}]){scenNum} * wBeam + (dij.(['mSqrtBetaDose' RBE_model{j}]){scenNum} * wBeam).^2);
                resultGUI.(['effect', RBE_model{j}, beamInfo(i).suffix])       = reshape(resultGUI.(['effect', RBE_model{j}, beamInfo(i).suffix]),dij.doseGrid.dimensions);

                % Calculate RBExD from the effect
                resultGUI.(['RBExD', RBE_model{j}, beamInfo(i).suffix])        = zeros(size(resultGUI.(['effect', RBE_model{j}, beamInfo(i).suffix])));
                resultGUI.(['RBExD', RBE_model{j}, beamInfo(i).suffix])(ix)    = (sqrt(dij.ax(ix).^2 + 4 .* dij.bx(ix) .* resultGUI.(['effect', RBE_model{j}, beamInfo(i).suffix])(ix)) - dij.ax(ix))./(2.*dij.bx(ix));

                % Divide RBExD with the physicalDose to get the plain RBE cube
                resultGUI.(['RBE', RBE_model{j}, beamInfo(i).suffix])          = resultGUI.(['RBExD', RBE_model{j}, beamInfo(i).suffix])./resultGUI.(['physicalDose', beamInfo(i).suffix]);

                % Initialize alpha/beta cubes
                resultGUI.(['alpha', RBE_model{j}, beamInfo(i).suffix])        = zeros(dij.doseGrid.dimensions);
                resultGUI.(['beta',  RBE_model{j}, beamInfo(i).suffix])        = zeros(dij.doseGrid.dimensions);
                resultGUI.(['alphaDoseCube', RBE_model{j}, beamInfo(i).suffix])        = zeros(dij.doseGrid.dimensions);
                resultGUI.(['SqrtBetaDoseCube',  RBE_model{j}, beamInfo(i).suffix])        = zeros(dij.doseGrid.dimensions);

                % Calculate alpha and weighted alphaDose
                AlphaDoseCube                                    = full(dij.(['mAlphaDose' RBE_model{j}]){scenNum} * wBeam);
                resultGUI.(['alpha', RBE_model{j}, beamInfo(i).suffix])(ix)    = AlphaDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix);
                resultGUI.(['alphaDoseCube', RBE_model{j}, beamInfo(i).suffix])(ix)   = AlphaDoseCube(ix);

                % Calculate beta and weighted sqrtBetaDose
                SqrtBetaDoseCube                                 = full(dij.(['mSqrtBetaDose' RBE_model{j}]){scenNum} * wBeam);
                resultGUI.(['beta', RBE_model{j}, beamInfo(i).suffix])(ix)     = (SqrtBetaDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix)).^2;
                resultGUI.(['SqrtBetaDoseCube', RBE_model{j}, beamInfo(i).suffix])(ix)   = SqrtBetaDoseCube(ix);
            end
        end
    end
end

%% Final processing
% Add non-processed MC tallies
% Note that the tallies are already computed per beam and altogether
if isfield(dij,'MC_tallies')
    for f = 1:numel(dij.MC_tallies)
        tally = dij.MC_tallies{f};
        % skip tallies processed above
        if ~isfield(resultGUI,tally) && ~contains(tally,'std','IgnoreCase',true)
            tallyCut = strsplit(tally,'_');
            if size(tallyCut,2) > 1
                beamNum = str2num(cell2mat(regexp(tally,'\d','Match')));
                resultGUI.(tally) = reshape(full(dij.(tallyCut{1}){scenNum}(:,beamNum)),dij.doseGrid.dimensions);
            else
                resultGUI.(tally) = reshape(full(dij.(tally){scenNum} * wBeam),dij.doseGrid.dimensions);
            end
        end
    end
end

% Remove suffix for RBExD if there's only one available
if any(contains(fieldnames(dij),'alpha','IgnoreCase',true)) && isfield(dij,'RBE_models') && length(dij.RBE_models) == 1
    % Get fieldnames that include the specified RBE model
    fnames = fieldnames(resultGUI);
    fnames = fnames(contains(fnames,dij.RBE_models{1},"IgnoreCase",true));

    % Rename fields and remove model specifier if there's only one
    for f = 1:length(fnames)
        resultGUI.(erase(fnames{f},['_',dij.RBE_models{1}])) = resultGUI.(fnames{1});
    end

    % Remove old fields
    resultGUI = rmfield(resultGUI,fnames);
end

% group similar fields together
resultGUI = orderfields(resultGUI);

% interpolation if dose grid does not match ct grid
if any(dij.ctGrid.dimensions~=dij.doseGrid.dimensions)
    myFields = fieldnames(resultGUI);
    for i = 1:numel(myFields)
        if numel(resultGUI.(myFields{i})) == dij.doseGrid.numOfVoxels

            % interpolate!
            resultGUI.(myFields{i}) = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
                resultGUI.(myFields{i}), ...
                dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);

        end
    end
end

end

