function [bixelDose] = matRad_calcParticleDoseBixelAPM(radDepths, latDistsX, latDistsZ, sigmaIni_sq, baseData,FlagBioOpt,vTissueIndex, sError)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the three dij components  (two lateral, one depth
% component) which are needed for robust treatment planning within the APM
% framework
% 
% call
%   dose = matRad_calcParticleDoseBixelAPM(radDepths,latDistsX, latDistsY, baseData, gaussBaseData, uct, BixelNum)
%
% input
%   radDepths:      radiological depths
%   latDistsX:      lateral x distance in BEV from central ray
%   latDistsZ:      lateral z distance in BEV from central ray
%   baseData:       base data required for particle dose calculation (contains the weights, widhts and means of)
%   BixelNum:       indicates the current bixel number
%
% output
%   Lx_ij:          particle dose at specified locations as linear vector
%   LY_ij:          particle dose at specified locations as linear vector
%   Z_ij:           particle dose at specified locations as linear vector
%
% References
%   [1] http://www.ncbi.nlm.nih.gov/pubmed/23877218
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function handle for calculating lateral dose
Gauss    = @(x,mu,SqSigma) 1./(sqrt(2*pi.*SqSigma)).*exp(-((x - mu).^2./(2.*SqSigma)));

% function handle for calculating depth doses                          
sumGauss = @(x,mu,SqSigma,w) ((1./sqrt(2*pi*ones(numel(x),1) * SqSigma') .* ...
                              exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) * SqSigma' ))) * w);

% add potential offset
depths   = baseData.depths + baseData.offset;

% convert from MeV cm^2/g per primary to Gy mm^2 per 1e6 primaries
conversionFactor = 1.6021766208e-02;

if ~isfield(baseData,'sigma')
    
    % interpolate depth dose, sigmas, and weights    
    X = matRad_interp1(depths,[baseData.sigma1 baseData.weight baseData.sigma2],radDepths);
    
    sigmaSq_Narr  = (X(:,1)).^2 + sigmaIni_sq;
    sigmaSq_Bro   = (X(:,3)).^2 + sigmaIni_sq;
    w             = (X(:,2));
    
else
    % interpolate depth dose and sigma
    X = matRad_interp1(depths,baseData.sigma,radDepths);

    %compute lateral sigma
    sigmaSq = X.^2 + sigmaIni_sq;
end


%% calculate lateral sigma in x direction %%%
latMeanOffsetX = 0;
latMeanOffsetZ = 0;
shiftErrorLat  = sError.latRndSq + sError.latSysSq;

if ~isfield(baseData,'sigma')
                          
    SqSigmaX_Narr = sigmaSq_Narr + shiftErrorLat;
    SqSigmaZ_Narr = sigmaSq_Narr + shiftErrorLat;
    SqSigmaX_Bro  = sigmaSq_Bro  + shiftErrorLat;
    SqSigmaZ_Bro  = sigmaSq_Bro  + shiftErrorLat;
    
    % only save the lateral dose of the narrow component
    NarrX = sqrt(1-w) .* Gauss(latDistsX,0,SqSigmaX_Narr);
    NarrZ = sqrt(1-w) .* Gauss(latDistsZ,0,SqSigmaZ_Narr);
    BroX  = sqrt(w)   .* Gauss(latDistsX,0,SqSigmaX_Bro);
    BroZ  = sqrt(w)   .* Gauss(latDistsZ,0,SqSigmaZ_Bro);
    
    bixelDose.L_ij  = NarrX .* NarrZ + BroX .* BroZ;
       
    radDist         = sqrt(latDistsX.^2 + latDistsZ.^2);
    bixelDose.Lx_ij = sqrt(((1-w) .* exp((-(radDist.^2))./(2*SqSigmaX_Narr))./(sqrt(SqSigmaX_Narr).*sqrt(SqSigmaZ_Narr).*2*pi)) +...
                              (w  .* exp((-(radDist.^2))./(2*SqSigmaX_Bro)) ./(sqrt(SqSigmaX_Bro) .*sqrt(SqSigmaZ_Bro) .*2*pi)));
    bixelDose.Lz_ij = sqrt(((1-w) .* exp((-(radDist.^2))./(2*SqSigmaZ_Narr))./(sqrt(SqSigmaX_Narr).*sqrt(SqSigmaZ_Narr).*2*pi)) +...
                              (w  .* exp((-(radDist.^2))./(2*SqSigmaZ_Bro)) ./(sqrt(SqSigmaX_Bro) .*sqrt(SqSigmaZ_Bro) .*2*pi))); 

else

    SqSigmaX = sigmaSq + shiftErrorLat;
    SqSigmaZ = sigmaSq + shiftErrorLat;
    
    bixelDose.Lx_ij = Gauss(latDistsX,latMeanOffsetX,SqSigmaX);
    bixelDose.Lz_ij = Gauss(latDistsZ,latMeanOffsetZ,SqSigmaZ);
    bixelDose.L_ij  = bixelDose.Lx_ij .*  bixelDose.Lz_ij;
    
end

%% calculate sigma in range direction 
SqSigmaRangeOffset = sError.rangeRndSq + sError.rangeSysSq;


bixelDose.Z_ij = conversionFactor * baseData.LatCutOff.CompFac * sumGauss(radDepths - baseData.offset,baseData.Z.mean,...
                                              (baseData.Z.width).^2 + SqSigmaRangeOffset,baseData.Z.weight);
                                                                                                               
if FlagBioOpt
   
    tissueClasses = unique(vTissueIndex);
    bixelDose.Z_Aij = zeros(numel(radDepths),1);
    bixelDose.Z_Bij = zeros(numel(radDepths),1);
    
    for i = 1:numel(tissueClasses)
       ix = vTissueIndex == tissueClasses(i);
       bixelDose.Z_Aij(ix)  = conversionFactor * baseData.LatCutOff.CompFac * ...
                              sumGauss(radDepths(ix)-baseData.offset,baseData.alphaDose(tissueClasses(i)).mean,...
                             (baseData.alphaDose(tissueClasses(i)).width).^2 + SqSigmaRangeOffset,baseData.alphaDose(tissueClasses(i)).weight);
             
       bixelDose.Z_Bij(ix)  = conversionFactor * baseData.LatCutOff.CompFac * ...
                              sumGauss(radDepths(ix)-baseData.offset,baseData.SqrtBetaDose(tissueClasses(i)).mean,...
                             (baseData.SqrtBetaDose(tissueClasses(i)).width).^2  + SqSigmaRangeOffset,baseData.SqrtBetaDose(tissueClasses(i)).weight)';  
    end
       
end



if ~isfield(baseData,'sigma')
    bixelDose.physDose =  bixelDose.L_ij  .* bixelDose.Z_ij;
else
    bixelDose.physDose =  bixelDose.Lx_ij .* bixelDose.Lz_ij .* bixelDose.Z_ij;
end

                                        
end


