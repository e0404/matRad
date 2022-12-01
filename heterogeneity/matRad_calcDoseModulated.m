function [resultGUI,pln] = matRad_calcDoseModulated(ct,stf,pln,cst,weights)

% Instance of matRad configuration class
matRad_cfg =  MatRad_Config.instance();

% load appropriate configs from pln or from class
if isfield(pln,'propMC')
    pln = matRad_cfg.getDefaultClass(pln,'propMC');
end
pln = matRad_cfg.getDefaultClass(pln,'propHeterogeneity');
pln.propHeterogeneity.calcHetero = false;

% load default parameters in case they haven't been set yet
pln = matRad_cfg.getDefaultProperties(pln,'propDoseCalc');

% Get number of samples from class
samples = pln.propHeterogeneity.sampling.numOfSamples;

% Output current mode
matRad_cfg.dispInfo(['Calculating modulated dose in mode: ' pln.propHeterogeneity.sampling.mode]);

switch pln.propHeterogeneity.sampling.mode
    case 'TOPAS'
        %         pln.propMC.materialConverter.densityCorrection = 'Schneider_TOPAS'; %'rspHLUT','Schneider_TOPAS','Schneider_matRad'
        pln.propMC.materialConverter.HUSection = 'advanced'; %'default','advanced'
        pln.propMC.materialConverter.HUToMaterial = 'default'; %'default','simpleLung','advanced'
        if strcmp(pln.propHeterogeneity.sampling.method,'poisson')
            pln.propMC.materialConverter.addSection = 'poisson';
        else
            pln.propMC.materialConverter.addSection = 'sampledDensities'; %'none','lung','poisson','sampledDensities' (the last 2 only with modulation)
        end

        numHistories = pln.propMC.numHistories;
        calcExternal = pln.propMC.externalCalculation;
    case 'MCsquare'
        calcExternal = false;
        numHistories = pln.propMC.numHistories;
        pln.propMC.materialConverter.addSection = 'sampledDensities';
    case 'matRad'
        calcExternal = false;
    otherwise
        matRad_cfg.dispError('No sampling mode other than TOPAS and matRad implemented');
end

% parallelComputationTOPAS =1;
% if parallelComputationTOPAS
%
%     for i = 1:samples
%         ct_mod{i} = pln.propHeterogeneity.modulateDensity(ct,cst,pln,Pmod,modulation);
%     end
%
%     pln.propMC.proton_engine = 'TOPAS';
%         if strcmp(modulation,'poisson')
%             pln.propMC.materialConverter = 'HUToWaterSchneider_mod';
%         else
%             if ~isprop(pln.propMC,'materialConverter')
%                 pln.propMC.materialConverter = 'HUToWaterSchneider';
%             end
%         end
%
%         resultGUI_mod = matRad_calcDoseDirectMC(ct_mod,stf,pln,cst,weights,mode{2}/samples,mode{3});
%
%         if ~mode{3}
%             %     resultGUI.(['physicalDose',num2str(s)]) = resultGUI.(['physicalDose',num2str(s)]) + resultGUI_mod.physicalDose/s;
%             if strcmp(pln.bioParam.quantityOpt,'RBExD')
%                 resultGUI.RBExD = resultGUI.RBExD + resultGUI_mod.RBExD/samples;
%             end
%             resultGUI.physicalDose = resultGUI.physicalDose + resultGUI_mod.physicalDose/samples;
%             std{i} = resultGUI_mod.std;
%         end
%
%
%
% else

% Turn info messages off for modulation
logLevel = matRad_cfg.logLevel;
matRad_cfg.logLevel = 2;

% Perform resampling to dose grid if necessary (modulation is performed on the resampled grid)
switch pln.propHeterogeneity.sampling.mode
    case 'TOPAS'
        pln.propMC.engine = 'TOPAS';
        [ctR,cstR,stfR] = matRad_resampleCTtoGrid(ct,cst,pln,stf);
    case 'MCsquare'
        pln.propMC.engine = 'MCsquare';
        ctR = ct;
        cstR = cst;
        stfR = stf;
    case 'matRad'
        ctR = ct;
        cstR = cst;
end

% Initialize waitbar if calculation locally
if ~calcExternal
    figureWait = waitbar(0,['calculate modulated dose for ' pln.propHeterogeneity.sampling.mode ': Sample ' num2str(1) '/' num2str(samples)]);
    % prevent closure of waitbar and show busy state
    set(figureWait,'pointer','watch');
end

% set this flag so that the modulated cube is not overwritten in matRad_calcDoseInit
pln.propDoseCalc.useGivenEqDensityCube = true;

% Allocate empty resultGUI and space for individual physical doses to calculate their standard deviation
resultGUI = struct;
data = cell(samples,1);

for i = 1:samples
    % Modulate density of ct cube
    ct_mod = pln.propHeterogeneity.modulateDensity(ctR,cstR,pln);

    % Save number of samples in modulated CT (e.g. used for TOPAS folder generation)
    ct_mod.sampleIdx = i;

    % Switch between the different modes of calculation currently implemented
    % WARNING: Implementation of MCsquare is currently not finished
    switch pln.propHeterogeneity.sampling.mode
        case {'TOPAS','MCsquare'}
            % Set TOPAS parameters
            pln.propMC.numOfRuns = 1;
            pln.propHeterogeneity.sampling.numHistories = numHistories/samples;

            % Calculate dose with modulated CT
            resultGUI_mod = matRad_calcDoseDirectMC(ct_mod,stfR,pln,cstR,weights);

            if ~calcExternal
                % Accumulate averaged results
                resultGUI = pln.propHeterogeneity.accumulateOverSamples(resultGUI,resultGUI_mod,samples);

                % Save individual standard deviation
                if isfield(resultGUI_mod,'physicalDose_std')
                    resultGUI.physicalDose_std_individual{i} = resultGUI_mod.physicalDose_std;
                end
            end
        case 'matRad'
            % Calculate dose with modulated CT
            resultGUI_mod = matRad_calcDoseDirect(ct_mod,stf,pln,cstR,weights);

            % Accumulate averaged results
            resultGUI = pln.propHeterogeneity.accumulateOverSamples(resultGUI,resultGUI_mod,samples);
    end

    % Save individual physical doses to calculate standard deviation
    if ~calcExternal
        data{i} = resultGUI_mod.physicalDose;
    end

    % Update waitbar if available
    if exist('figureWait','var') && ishandle(figureWait)
        waitbar(i/samples,figureWait,['calculate modulated dose for ' pln.propHeterogeneity.sampling.mode ': Sample ' num2str(i+1) '/' num2str(samples)]);
    end

end

if ~calcExternal
    % Calculate standard deviation between samples
    resultGUI.physicalDose_std = pln.propHeterogeneity.calcSampleStd(data,resultGUI.physicalDose);
end

% Close Waitbar
if exist('figureWait','var') && ishandle(figureWait)
    delete(figureWait);
end

% Change loglevel back to default;
matRad_cfg.logLevel = logLevel;

end

