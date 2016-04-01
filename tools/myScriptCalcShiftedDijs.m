%% calculate a couple of shifted dijs

shifts = [ 5  0  0; ...
          -5  0  0; ...
           0  5  0; ...
           0 -5  0; ...
           0  0  5; ...
           0  0 -5]/2;
       
       
dijShifts                  = dij;
dijShifts.physicalDose     = [];
dijShifts.physicalDose{1}  = dij.physicalDose;
if isfield(dij,'mAlphaDose')
    dijShifts.mAlphaDose       = [];
    dijShifts.mAlphaDose{1}    = dij.mAlphaDose;
    dijShifts.mSqrtBetaDose    = [];
    dijShifts.mSqrtBetaDose{1} = dij.mSqrtBetaDose;
end

for i = 1:size(shifts,1)
    
    stfTmp = matRad_shiftStf(stf,ct,shifts(i,:));
    
    dijTmp = matRad_calcParticleDose(ct,stfTmp,pln,cst);
    
    dijShifts.physicalDose{i+1} = dijTmp.physicalDose;
    
    if isfield(dijTmp,'mAlphaDose')
        dijShifts.mAlphaDose{i+1}    = dijTmp.mAlphaDose;
        dijShifts.mSqrtBetaDose{i+1} = dijTmp.mSqrtBetaDose;
    end
    
end

dijShifts.numOfScenarios = numel(dijShifts.physicalDose);

dijShifts.probOfScenarios = ones(dijShifts.numOfScenarios,1)/dijShifts.numOfScenarios;

clear dijTmp;

%% adjust cst

for i = 1:size(cst,1)
   
    for j = 1:numel(cst{i,6})
        
        % options 'none', 'probabilistic', 'voxel-wise worst case', 'objective-wise
        % worst case'
        
        cst{i,6}(j).robustness = 'voxel-wise worst case';
        %cst{i,6}(j).robustness = 'none';
        %cst{i,6}(j).robustness = 'probabilistic';
        
    end
    
end

clear i j;

%%

resultGUI.physicalDose = reshape(dijShifts.physicalDose{2}*resultGUI.w,dij.dimensions);

