function resultGUI = matRad_calcDoseModulated(ct,stf,pln,cst,mode,weights,samples,Pmod,modulation,continuous)

global matRad_cfg;
matRad_cfg =  MatRad_Config.instance();

if nargin < 9
    modulation = 'binomial';
    continuous = false;
elseif nargin < 10
    continuous = false;
end

if iscell(mode)
    pln.propHeterogeneity.mode = 'TOPAS';
    if strcmp(modulation,'poisson')
        pln.propMC.materialConverter = 'HUToWaterSchneider_mod';
    else
        if ~isfield(pln.propMC,'materialConverter')
            pln.propMC.materialConverter = 'HUToWaterSchneider_custom';
        end
    end
    pln.propMC.materialConverter = 'HUToWaterSchneider_custom';
else
    pln.propHeterogeneity.mode = 'matRad';
end

resultGUI.physicalDose = zeros(ct.cubeDim);
if strcmp(pln.bioParam.quantityOpt,'RBExD')
    resultGUI.RBExD = zeros(ct.cubeDim);
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
matRad_cfg.logLevel = 1;

for i = 1:samples
    fprintf('Dose calculation for CT %i/%i \n',i,samples)
    ct_mod = matRad_modulateDensity(ct,cst,pln,Pmod,modulation,continuous);
    ct_mod.sampleIdx = i;
    
    if iscell(mode) %case TOPAS
        if strcmp(mode{1},'TOPAS')
            pln.propMC.proton_engine = 'TOPAS';
            pln.propMC.numOfRuns = 1;
            resultGUI_mod = matRad_calcDoseDirectMC(ct_mod,stf,pln,cst,weights,mode{2}/samples,mode{3});
            
            if ~mode{3}
                %     resultGUI.(['physicalDose',num2str(s)]) = resultGUI.(['physicalDose',num2str(s)]) + resultGUI_mod.physicalDose/s;
                if strcmp(pln.bioParam.quantityOpt,'RBExD')
                    resultGUI.RBExD = resultGUI.RBExD + resultGUI_mod.RBExD/samples;
                end
                resultGUI.physicalDose = resultGUI.physicalDose + resultGUI_mod.physicalDose/samples;
                std{i} = resultGUI_mod.physicalDose_std;
            end
        else
            error('error')
        end
    else %case matRad
        resultGUI_mod = matRad_calcDoseDirect(ct_mod,stf,pln,cst,weights);
        
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

