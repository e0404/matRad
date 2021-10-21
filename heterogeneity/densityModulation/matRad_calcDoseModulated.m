function resultGUI = matRad_calcDoseModulated(ct,stf,pln,cst,param,weights,samples,Pmod,modulation,continuous)

global matRad_cfg;
matRad_cfg =  MatRad_Config.instance();

if nargin < 10
    continuous = true;
end
if nargin < 9
    modulation = 'binomial';
end
if nargin < 7
    samples = 50;
end
if nargin < 8
    Pmod = 800;
end

if iscell(param)
    pln.propHeterogeneity.mode = param{1};
    histories = param{2};
    calcOpenstack = param{3};
    
    if strcmp(pln.propHeterogeneity.mode,'TOPAS')
        
        pln.propHeterogeneity.mode = 'TOPAS';
        pln.propMC.materialConverter.densityCorrection.mode = 'TOPAS2'; %'default','TOPAS1','TOPAS2'
        pln.propMC.materialConverter.HUSection = 'advanced'; %'default','advanced'
        pln.propMC.materialConverter.HUToMaterial = 'advanced'; %'default','simpleLung','advanced'
        if strcmp(modulation,'poisson')
            pln.propMC.materialConverter.densityCorrection.addSection = 'poisson';
        else
            pln.propMC.materialConverter.densityCorrection.addSection = 'sampledDensities'; %'none','lung','poisson','sampledDensities' (the last 2 only with modulation)
        end
    else
        matRad_cfg.dispError('No sampling mode other than TOPAS implemented');
    end
else
    pln.propHeterogeneity.mode = 'matRad';
end

% parallelComputationTOPAS =1;
% if parallelComputationTOPAS
%
%     for i = 1:samples
%         ct_mod{i} = matRad_modulateDensity(ct,cst,pln,Pmod,modulation);
%     end
%
%     pln.propMC.proton_engine = 'TOPAS';
%         if strcmp(modulation,'poisson')
%             pln.propMC.materialConverter = 'HUToWaterSchneider_mod';
%         else
%             if ~isfield(pln.propMC,'materialConverter')
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

% Turn info and warning messages off for modulation
logLevel = matRad_cfg.logLevel;
matRad_cfg.logLevel = 3;
    
if iscell(param) %case TOPAS
    [ctR,cstR] = matRad_resampleTopasGrid(ct,cst,pln,stf);
else %case matRad
    ctR = ct;
    cstR = cst;
end
resultGUI.physicalDose = zeros(ct.cubeDim);
if strcmp(pln.bioParam.quantityOpt,'RBExD')
    resultGUI.RBExD = zeros(ct.cubeDim);
end

% set this flag so that the modulated cube is not overwritten in matRad_calcDoseInit
pln.propDoseCalc.useGivenEqDensityCube = true;

for i = 1:samples
    fprintf('Dose calculation for CT %i/%i \n',i,samples)
    ct_mod = matRad_modulateDensity(ctR,cstR,pln,Pmod,modulation,continuous);
    ct_mod.sampleIdx = i;
    %%  
    if iscell(param) %case TOPAS
        if strcmp(pln.propHeterogeneity.mode,'TOPAS')
            pln.propMC.proton_engine = 'TOPAS';
            pln.propMC.numOfRuns = 1;
            resultGUI_mod = matRad_calcDoseDirectMC(ct_mod,stf,pln,cstR,weights,histories/samples,calcOpenstack);
            
            if ~calcOpenstack
                %     resultGUI.(['physicalDose',num2str(s)]) = resultGUI.(['physicalDose',num2str(s)]) + resultGUI_mod.physicalDose/s;
                if strcmp(pln.bioParam.quantityOpt,'RBExD')
                    resultGUI.RBExD = resultGUI.RBExD + resultGUI_mod.RBExD/samples;
                end
                resultGUI.physicalDose = resultGUI.physicalDose + resultGUI_mod.physicalDose/samples;
                std{i} = resultGUI_mod.physicalDose_std;
            end
        else
            error('Make sure you selected the correct environment for modulation.')
        end
    else %case matRad
        resultGUI_mod = matRad_calcDoseDirect(ct_mod,stf,pln,cstR,weights);
        
        if strcmp(pln.bioParam.quantityOpt,'RBExD')
            resultGUI.RBExD = resultGUI.RBExD + resultGUI_mod.RBExD/samples;
        end
        resultGUI.physicalDose = resultGUI.physicalDose + resultGUI_mod.physicalDose/samples;
    end
    data{i} = resultGUI_mod.physicalDose;
    
end

% Calculate Standard deviation
topasMeanDiff = 0;
for k = 1:samples
    topasMeanDiff = topasMeanDiff + (data{k} - resultGUI.physicalDose).^2;
end
topasVarMean = topasMeanDiff./(samples - 1)./samples;
topasStdMean = sqrt(topasVarMean);

topasStdSum = topasStdMean * samples;
topasVarSum = topasStdSum.^2;

resultGUI.physicalDose_std = sqrt(topasVarSum);

%Change loglevel back to default;
matRad_cfg.logLevel = logLevel;

end

