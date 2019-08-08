function pln = matRad_checkHeterogeneity(pln,cst,param)

% Setting default parameters for heterogeneity correction
if isfield(pln,'heterogeneity')
    if ~isfield(pln.heterogeneity,'useDoseCurves') || isempty(pln.heterogeneity.useDoseCurves)
        pln.heterogeneity.useDoseCurves = false;
    end
    if ~isfield(pln.heterogeneity,'calcHetero') || isempty(pln.heterogeneity.calcHetero)
        pln.heterogeneity.calcHetero = false;
    end
    if ~isfield(pln.heterogeneity,'type') || isempty(pln.heterogeneity.type)
        if pln.heterogeneity.calcHetero
            pln.heterogeneity.type = 'complete';
            warning('Heterogeneity correction wanted but no correction type specified, default was selected.'); 
        else
            pln.heterogeneity.type = 'none';
        end
    end
else
     pln.heterogeneity.useDoseCurves = false;
     pln.heterogeneity.calcHetero = false; 
     pln.heterogeneity.type = 'none';
end

if pln.heterogeneity.calcHetero
    matRad_dispToConsole('Heterogeneity correction enabled. \n',param,'info');
    heteroCST = false;
    for i = 1:length(cst(:,1))
        if isfield(cst{i,5},'heterogeneityCorrection')
            heteroCST = true;
            break
        end
    end
    if ~heteroCST
       warning('Heterogeneity correction enabled but no usable data in cst.'); 
    end
else
    matRad_dispToConsole('Heterogeneity correction disabled. \n',param,'info');
end
end