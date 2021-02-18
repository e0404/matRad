function resultGUI = matRad_calcDoseModulated(ct,stf,pln,cst,mode,weights,samples,Pmod,modulation)

if nargin < 9
    modulation = 'binomial';
end

resultGUI.physicalDose = zeros(ct.cubeDim);
if strcmp(pln.bioParam.quantityOpt,'RBExD')
    resultGUI.RBExD = zeros(ct.cubeDim);
end

for i = 1:samples
    ct_mod = matRad_modulateDensity(ct,cst,pln,Pmod,modulation);
    
    if iscell(mode) %case TOPAS
        if strcmp(mode{1},'TOPAS')
            pln.propMC.proton_engine = 'TOPAS';
            if strcmp(modulation,'poisson')
                pln.propMC.materialConverter = 'HUToWaterSchneider_mod';
            else
                pln.propMC.materialConverter = 'HUToWaterSchneider';
            end
            resultGUI_mod = matRad_calcDoseDirectMC(ct_mod,stf,pln,cst,weights,mode{2},mode{3});
        else
            error('error')
        end
    else %case matRad
        resultGUI_mod = matRad_calcDoseDirect(ct_mod,stf,pln,cst,weights);
    end
    
    %     resultGUI.(['physicalDose',num2str(s)]) = resultGUI.(['physicalDose',num2str(s)]) + resultGUI_mod.physicalDose/s;
    if strcmp(pln.bioParam.quantityOpt,'RBExD')
        resultGUI.RBExD = resultGUI.RBExD + resultGUI_mod.RBExD/samples;
    end
    resultGUI.physicalDose = resultGUI.physicalDose + resultGUI_mod.physicalDose/samples;
end

end

