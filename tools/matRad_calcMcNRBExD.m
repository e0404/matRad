function RBExDose  = matRad_calcMcNRBExD(dij, cst, resultGUI)

dimensions = [size(resultGUI.physicalDose,1) size(resultGUI.physicalDose,2) size(resultGUI.physicalDose,3)]; 

AlphaDose = zeros(dimensions);
SqrtBetaDose = zeros(dimensions);

AlphaDose = reshape(dij.mAlphaDose{1} * resultGUI.w, dij.dimensions);
SqrtBetaDose = reshape(dij.mSqrtBetaDose{1} * resultGUI.w, dij.dimensions);
  
a_x = zeros(size(resultGUI.physicalDose));
b_x = zeros(size(resultGUI.physicalDose));

for i = 1:size(cst,1)
    if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') 
            a_x(cst{i,4}{1}) = cst{i,5}.alphaX;
            b_x(cst{i,4}{1}) = cst{i,5}.betaX;
    end    
end

% only compute where we have biologically defined tissue
ix = a_x~=0; 
    
Effect = AlphaDose+SqrtBetaDose.^2;   
    
RBExDose     = zeros(dimensions);
RBExDose(ix) = ((sqrt(a_x(ix).^2 + 4 .* b_x(ix) .* Effect(ix)) - a_x(ix))./(2.*b_x(ix)));
           
    
% only compute where we have finite dose
ix = resultGUI.physicalDose~=0; 
    
RBE     = zeros(dimensions);
RBE(ix) = RBExDose(ix)./resultGUI.physicalDose(ix);
