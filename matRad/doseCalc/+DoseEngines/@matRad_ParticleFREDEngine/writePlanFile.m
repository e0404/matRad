function writePlanFile(this, fName, stf)
% FRED helper to write data to plan.inp file
% call
%   writePlanFile(fName, stf)
% 
% input
%   fName: string specifying the file path and name for saving the data.
%   stf:   Fred stf struct
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matRad_cfg = MatRad_Config.instance();

fID = fopen(fName, 'w');

try
    totalNumOfBixels = sum([stf.totalNumOfBixels]);
    if this.calcDoseDirect
        simulatedPrimariesPerBixel = max([1, floor(this.numHistoriesDirect/totalNumOfBixels)]);
    else
        simulatedPrimariesPerBixel = this.numHistoriesPerBeamlet;
    end
    
    fprintf(fID, 'nprim = %i\n', simulatedPrimariesPerBixel);

    layerCounter = 0;
    bixelCounter = 0;

    % loop oever the fields
    for i=1:numel(stf)

        %Loop over energy layers
        for j=1:numel(stf(i).energies)

            fprintf(fID, '#Bixels Field%i, Layer%i\n', i-1,layerCounter+j-1);

            % Print bixel info (ID, Position, Direction, Weight)
            for k=1:stf(i).energyLayer(j).nBixels
                currBixel.beamletID = num2str(bixelCounter+k-1);
                currBixel.P         = arrayfun(@(idx) num2str(idx, '%2.3f'), [stf(i).energyLayer(j).rayPosY(k),stf(i).energyLayer(j).rayPosX(k),0], 'UniformOutput', false);
                currBixel.v         = arrayfun(@(idx) num2str(idx, '%2.5f'), [stf(i).energyLayer(j).rayDivY(k),stf(i).energyLayer(j).rayDivX(k),1], 'UniformOutput', false);
                currBixel.w         = num2str(stf(i).energyLayer(j).numOfPrimaries(k), '%2.7f');

                printStructToDictionary(fID, currBixel, ['S', num2str(bixelCounter+k-1)],2);
                
            end
                
            currLayer.Energy   = num2str(stf(i).energies(j));
            currLayer.Espread  = num2str(stf(i).energySpreadFWHMMev(j));

            % Select source model
            switch this.sourceModel
                case 'gaussian'
                    currLayer.FWHM     = num2str(stf(i).FWHMs(j));
                case 'emittance'
                    currLayer.emittanceX  = num2str(stf(i).emittanceX(j), '%1.10f');
                    currLayer.twissAlphaX = num2str(stf(i).twissAlphaX(j),'%1.10f');
                    currLayer.twissBetaX  = num2str(stf(i).twissBetaX(j), '%1.10f');
                case 'sigmaSqrModel'
                    currLayer.sSQr_a = num2str(stf(i).sSQr_a(j));
                    currLayer.sSQr_b = num2str(stf(i).sSQr_b(j));
                    currLayer.sSQr_c = num2str(stf(i).sSQr_c(j));
            end

            % Specify the beamlets in current layer
            currLayer.beamlets = arrayfun(@(idx) ['S', num2str(idx)], bixelCounter:stf(i).energyLayer(j).nBixels+bixelCounter-1, 'UniformOutput', false);

            %Print layer
            printStructToDictionary(fID, currLayer, ['L', num2str(layerCounter+j-1)],1);
            fprintf(fID, '\n');
            bixelCounter = bixelCounter + stf(i).energyLayer(j).nBixels;
        end
    
        % Estimate field dimension
        fieldLim = max(abs([stf(i).energyLayer.rayPosX,stf(i).energyLayer.rayPosY])) + 10*max([stf(i).FWHMs]);

        %Write field parameters
        currF.fieldNumber = i-1;
        currF.GA          = num2str(stf(i).gantryAngle);
        currF.CA          = num2str(stf(i).couchAngle);
        currF.ISO         = arrayfun(@num2str, stf(i).isoCenter, 'UniformOutput', false);
        currF.dim         = arrayfun(@num2str, [fieldLim, fieldLim, 0.1], 'UniformOutput', false);
        currF.Layers      = arrayfun(@(idx) ['L', num2str(idx)], layerCounter:numel(stf(i).energies)+layerCounter-1, 'UniformOutput', false);

        layerCounter = layerCounter + numel(stf(i).energies);

        printStructToDictionary(fID, currF, ['F', num2str(i-1)]);
        fprintf(fID, '\n');
    end

    %% Build plan

    % BAMsToIso is the same for all fields
    plan.SAD    = stf(1).BAMStoIsoDist;
    plan.Fields = arrayfun(@(i) ['F', num2str(i)], 0:numel(stf)-1, 'UniformOutput', false);
            
    printStructToDictionary(fID, plan, 'plan');
    
catch
    matRad_cfg.dispError('Failed to write plan file');
end

fclose(fID);
end

function printStructToDictionary(fID, S, sName, indentTabs)
% Helper function to convert struct fields into FRED specific python dictionary 
% call
%   printStructToDictionary(fID, S, sName, indentTabs)
% 
% input
%   fID:        ID of file to write
%   S:          struct to convert
%   sName:      variable name of the printed dictionary
%   indentTabs: optional, adds indentation to the left of pruinted code
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('indentTabs', 'var') || isempty(indentTabs)
    indentTabs = 0;
end

indentString = repmat('\t',1,indentTabs);


fprintf(fID, indentString);

% variable defintion
fprintf(fID, 'def: %s = {', sName);

% Get all fields to print
sFields = fieldnames(S);


for sFieldIdx =1:numel(sFields)

    currField = sFields{sFieldIdx};

    % write field name
    fprintf(fID, '''%s'': ',currField);

    % write one or multiple values
    if ~iscell(S.(currField))
        fprintf(fID, '%s', num2str(S.(currField)));
    else
        fprintf(fID, '[');
        for elementIdx=1:numel(S.(currField))-1
            fprintf(fID, '%s, ', S.(currField){elementIdx});
        end
        fprintf(fID, '%s]', S.(currField){end});

    end

    if sFieldIdx ~= numel(sFields)
        fprintf(fID, ', ');
    end
end

% Close dictionary defintion
fprintf(fID, '}\n');
end