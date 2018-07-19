function [PSI_y_corr_ijm, PSI_y_uncorr_ijm,LinIx,LinIxInv] =  matRad_calcSecRangeMomFast(amplificationFactor, numBixels,numRays,vEnergyIx,vBeamIx,vCandidates,Z,mMean,mWeight,mWidth,mDepth,LUT,...
   randError,sysError,mRandCovRadDepth, mSysCovRadDepth,mBioCov,relUCTbio, mCovSpot) %#codegen
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the second raw moment for depth dose profile.
% This function is optimized for ray wise block correlation and should not
% be used for other correlatin assumptions. This function utilizes the
% following corr assumption for range erros:
% All spots on the same ray are perfectly correlated whereas
% spots from different rays are uncorrelated
%
% call
%   [PSI_y_corr_ijm, PSI_y_uncorr_ijm,LinIx,LinIxInv] =  matRad_calcSecRangeMomFast(amplificationFactor, numBixels,numRays,vEnergyIx,vBeamIx,vCandidates,Z,mMean,mWeight,mWidth,mDepth,LUT,...
%                                                                           randError,sysError,mRandCovRadDepth, mSysCovRadDepth)
%
% input
%   amplificationFactor:      factor to will be multiplied with the resulting second central moment (converstion factor, const RBE etc.)
%   numBixels:                number of bixels
%   numRays:                  number of rays
%   vEnergyIx:                energy index vector for all spots
%   vBeamIx:                  beam bixel look up vector
%   vCandidates:              correlated spots
%   Z:                        radiological depth of voxel i of beam s
%   mMean:                    base data - mean of gaussian components
%   mWeight:                  base data - weight of gaussian components
%   mWidth:                   base data - width of gaussian components
%   mDepth:                   base data - depth values
%   LUT:                      base data - look up table to use only a reduced number of gaussian components
%   randError:                random range error
%   sysError:                 systematic range error
%   mRandCovRadDepth:         random correlation between bixels
%   mSysCovRadDepth:          systematic correlation between bixels
%   mBioCov:                  correlation for bio uncertainties - usually its the same structure as for range uncertainties - ray-wise correlation
%   relUCTbio                 relative uncertainty of weights 0.2 equals a 20 % uncertainty
%   mCovSpot                  covariance matrix of gaussian weights
%
%
% output
%   PSI_y_corr_ijm:           second central moment of the depth dose in y assuming correlation (single fraction term)
%   PSI_y_uncorr_ijm:         second central moment of the lateral dose in y assuming no random correlation (only systematic)
%   LinIx:                    linear indices
%   LinIxInv:                 inverted linear indices for symmetric filling
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 Wieser Hans-Peter
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


assert(isa(numBixels,'double'));
assert(isa(numRays,'double'));
assert(isa(vEnergyIx,'double'));
assert(isa(vBeamIx,'double'));
assert(isa(vCandidates,'double'));
assert(isa(Z,'double'));
assert(isa(mMean,'double'));
assert(isa(mWeight,'double'));
assert(isa(mWidth,'double'));
assert(isa(mDepth,'double'));
assert(isa(LUT,'logical'));
assert(isa(randError,'double'));
assert(isa(sysError,'double'));
assert(isa(mRandCovRadDepth,'double'));
assert(isa(mSysCovRadDepth,'double'));
assert(isa(mBioCov,'double'));
assert(isa(amplificationFactor,'double'));
assert(isa(relUCTbio,'double'));
assert(isa(mCovSpot,'double'));

upperBound     = 25000;
lowerBound     = 15;
MaxNumDepth    = 550;
MaxNumEnergies = 300;

coder.varsize('numBixels',[1 1],[0 0]);
coder.varsize('numRays',[1 1],[0 0]);
coder.varsize('vEnergyIx',[1 upperBound],[0 1]);
coder.varsize('vBeamIx',[1 upperBound],[0 1]);
coder.varsize('vCandidates',[1 upperBound],[0 1]);
coder.varsize('Z',[1 lowerBound],[0 1]);
coder.varsize('mMean',[lowerBound MaxNumEnergies],[1 1]);
coder.varsize('mWeight',[lowerBound MaxNumEnergies],[1 1]);
coder.varsize('mWidth',[lowerBound MaxNumEnergies],[1 1]);
coder.varsize('mDepth',[MaxNumDepth MaxNumEnergies],[1 1]);
coder.varsize('LUT',[15  MaxNumDepth MaxNumEnergies],[1 1 1]);
coder.varsize('randError',[1 1],[0 0]);
coder.varsize('sysError',[1 1],[0 0]);
coder.varsize('mRandCovRadDepth',[upperBound upperBound],[1 1]);
coder.varsize('mSysCovRadDepth',[upperBound upperBound],[1 1]);
coder.varsize('mBioCov',[upperBound upperBound],[1 1]);
coder.varsize('amplificationFactor',[1 upperBound],[0 1]);
coder.varsize('relUCTbio',[1 1],[0 0]);
coder.varsize('mCovSpot',[lowerBound upperBound],[1 1]);

numComp = size(mWeight,1);

estimatedElements = max([nnz(mRandCovRadDepth) nnz(mSysCovRadDepth)  nnz(mBioCov)]);
PSI_y_corr_ijm    = zeros(estimatedElements,1);
PSI_y_uncorr_ijm  = zeros(estimatedElements,1);
LinIx             = zeros(estimatedElements,1);
LinIxInv          = zeros(estimatedElements,1);

ixSave = 1;

for jj = 1:numel(vCandidates)
   
   

   j = vCandidates(jj);
   
   bixelIxDepthBioPhys = find((mBioCov(jj,jj:end))) + jj - 1;
   
   bixelIxDepthPhys    = find((mRandCovRadDepth(jj,jj:end)> 0 | ...
                               mSysCovRadDepth(jj,jj:end) > 0 )) + jj - 1;
                
   bixelIxDepthBio     = bixelIxDepthBioPhys(~ismember(bixelIxDepthBioPhys,bixelIxDepthPhys));               
      
   bixelIxDepth        = unique([bixelIxDepthPhys bixelIxDepthBioPhys]);  
   PSI_y_corr_ijmTmp   = zeros(numel(bixelIxDepth),1);
   PSI_y_uncorr_ijmTmp = zeros(numel(bixelIxDepth),1);
   
   
   lengthBx     = numel(bixelIxDepthBio);
   vBeamIx_j    = vBeamIx(j);
   vEnergyIx_j  = vEnergyIx(j);
   ampliFac     = amplificationFactor(bixelIxDepthBio);
   
   ixLut_j      = find(mDepth(:,vEnergyIx_j) > Z(1,vBeamIx_j),1,'first');  % find depth index
   vIxLUT_j     = logical(squeeze(LUT(:,ixLut_j,vEnergyIx_j)));
   
   Dev_j        = Z(1,vBeamIx_j) - mMean(vIxLUT_j,vEnergyIx_j);
   w_j          = mWeight(vIxLUT_j,vEnergyIx_j)';
   
   vBeamIx_m    = vBeamIx(vCandidates(bixelIxDepthBio));
   vEnergyIx_m  = vEnergyIx(vCandidates(bixelIxDepthBio));
   mIxLUT_m     = false(numComp,lengthBx);
   
   [~,ixLut_m]  = min(abs(bsxfun(@minus,mDepth(:,vEnergyIx_m),Z(1,vBeamIx_m))),[],1);
   
   cntLUT = 1;
   for p = 1:lengthBx
      mIxLUT_m(:,cntLUT) = LUT(:,ixLut_m(cntLUT),vEnergyIx_m(cntLUT));
      cntLUT             = cntLUT + 1;
   end

   SigmaJsq = mWidth(vIxLUT_j,vEnergyIx_j) + mRandCovRadDepth(jj,jj) + mSysCovRadDepth(jj,jj);
   
   Z_j    = exp(-((Dev_j).^2./(2.*SigmaJsq)))./(sqrt(2*pi.*SigmaJsq));
   
   % calculate second raw moments only for biological uncertainties 
   for mm = 1:lengthBx
       
      energyIx_m = vEnergyIx_m(mm);
      vIxLUT_m   = mIxLUT_m(:,mm);
      
      Dev_m    = Z(1,vBeamIx_m(mm)) - mMean(vIxLUT_m,energyIx_m);
      w_m      = mWeight(vIxLUT_m,energyIx_m);
      SigmaMsq = mWidth(vIxLUT_m,energyIx_m) + mRandCovRadDepth(mm,mm) + mSysCovRadDepth(mm,mm);
      
      Z_m  = exp(-((Dev_m).^2./(2.*SigmaMsq)))./(sqrt(2*pi.*SigmaMsq));
      Z_jm = (Z_j * Z_m');
      mW_CovAlphaDose_jm = mCovSpot(vIxLUT_j,vIxLUT_m) .* ((w_j' * relUCTbio) * (w_m' * relUCTbio));
      mBioOffset         = mW_CovAlphaDose_jm  + (w_j' * w_m');
      vHelpJ             = ones(numel(w_j),1)';
      vHelpM             = ones(numel(w_m),1);
         
      logIx = (bixelIxDepthBio(mm)==bixelIxDepth);
      PSI_y_corr_ijmTmp(logIx)   = amplificationFactor(jj) * (vHelpJ * (Z_jm .* mBioOffset) * vHelpM); 
      PSI_y_uncorr_ijmTmp(logIx) = PSI_y_corr_ijmTmp(logIx); 
      
   end
   
   % calculate second raw moments for biological and physical uncertainties 
   lengthBx     = numel(bixelIxDepthPhys);
   vBeamIx_j    = vBeamIx(j);
   vEnergyIx_j  = vEnergyIx(j);
   ampliFac     = amplificationFactor(bixelIxDepthPhys);
   
   ixLut_j      = find(mDepth(:,vEnergyIx_j) > Z(1,vBeamIx_j),1,'first');  % find depth index
   vIxLUT_j     = logical(squeeze(LUT(:,ixLut_j,vEnergyIx_j)));
   
   Dev_j        = Z(1,vBeamIx_j) - mMean(vIxLUT_j,vEnergyIx_j);
   w_j          = mWeight(vIxLUT_j,vEnergyIx_j)';
   LinInd_MM    = (bixelIxDepthPhys + (bixelIxDepthPhys-1)*numel(vCandidates))';
  
   vBeamIx_m    = vBeamIx(vCandidates(bixelIxDepthPhys));
   vEnergyIx_m  = vEnergyIx(vCandidates(bixelIxDepthPhys));
   mIxLUT_m     = false(numComp,lengthBx);
   
   [~,ixLut_m] = min(abs(bsxfun(@minus,mDepth(:,vEnergyIx_m),Z(1,vBeamIx_m))),[],1);
   
   cntLUT = 1;
   for p = 1:lengthBx
      mIxLUT_m(:,cntLUT)    = LUT(:,ixLut_m(cntLUT),vEnergyIx_m(cntLUT));
      cntLUT                = cntLUT + 1;
   end

   if randError > 0
      randCovRadDepth     = mRandCovRadDepth(jj,jj);
      vRandCovRadDepth    = mRandCovRadDepth(jj,bixelIxDepthPhys)';
      vRandCovRadDepth22  = mRandCovRadDepth(LinInd_MM);
   else
      randCovRadDepth = 0; vRandCovRadDepth = 0; vRandCovRadDepth22 = 0;
   end
   
   if sysError > 0
      sysCovRadDepth      = mSysCovRadDepth(jj,jj);
      vSysCovRadDepth     = mSysCovRadDepth(jj,bixelIxDepthPhys)';
      vSysCovRadDepth22   = mSysCovRadDepth(LinInd_MM);
   else
      sysCovRadDepth = 0; vSysCovRadDepth = 0; vSysCovRadDepth22   = 0;
   end
   
   
   for mm = 1:lengthBx
      
      energyIx_m           = vEnergyIx_m(mm);
      SqSigmaRandUnCorr_jm = [0; 0; 0; 0];
      vIxLUT_m             = mIxLUT_m(:,mm);
      
      if randError > 0
         SqSigmaRandCorr_jm   = [randCovRadDepth; vRandCovRadDepth(mm); vRandCovRadDepth(mm); vRandCovRadDepth22(mm)];
         SqSigmaRandUnCorr_jm = [randCovRadDepth; 0;  0; vRandCovRadDepth22(mm)];
      else
         SqSigmaRandCorr_jm = [0; 0; 0; 0];
      end
      
      if sysError > 0
         SqSigmaSys_jm = [sysCovRadDepth; vSysCovRadDepth(mm);  vSysCovRadDepth(mm); vSysCovRadDepth22(mm)];
      else
         SqSigmaSys_jm = [0; 0; 0; 0];
      end
      
      Dev_m       = Z(1,vBeamIx_m(mm)) - mMean(vIxLUT_m,energyIx_m);
      w_m         = mWeight(vIxLUT_m,energyIx_m);
      SigmaCorr   = SqSigmaSys_jm + SqSigmaRandCorr_jm;
      SigmaUnCorr = SqSigmaSys_jm + SqSigmaRandUnCorr_jm;
      
      vLaSi11 = mWidth(vIxLUT_j,vEnergyIx_j) + SigmaCorr(1,1);
      vLaSi22 = mWidth(vIxLUT_m,energyIx_m)  + SigmaCorr(4,1);
      vLaSi12 = SigmaCorr(2,1);
      vLaSi21 = SigmaCorr(3,1);
      
      mW_CovAlphaDose_jm = mCovSpot(vIxLUT_j,vIxLUT_m) .* ((w_j' * relUCTbio) * (w_m' * relUCTbio));
      mBioOffset         = mW_CovAlphaDose_jm  + (w_j' * w_m');
      vHelpJ             = ones(numel(w_j),1)';
      vHelpM             = ones(numel(w_m),1);
         
      Det      = vLaSi11*vLaSi22' - (vLaSi12*vLaSi21');
      FracDet  = 1./(2*pi*sqrt(Det));
      ExpTerm  = FracDet .* exp(-.5./Det.*((Dev_j.^2)*vLaSi22' - bsxfun(@times,(Dev_j*Dev_m'),(vLaSi12+vLaSi21)) + vLaSi11*(Dev_m.^2)'));
      PSI_y_corr_ijmTmp(mm) = amplificationFactor(jj) * (vHelpJ * (ExpTerm .* mBioOffset) * vHelpM);
      
      vLaSi12 = SigmaUnCorr(2,1);
      vLaSi21 = SigmaUnCorr(3,1);
      
      Det      = vLaSi11*vLaSi22' - (vLaSi12*vLaSi21');
      FracDet  = 1./(2*pi*sqrt(Det));
      ExpTerm  = FracDet .* exp(-.5./Det.*((Dev_j.^2)*vLaSi22' - bsxfun(@times,(Dev_j*Dev_m'),(vLaSi12+vLaSi21)) + vLaSi11*(Dev_m.^2)'));
      PSI_y_uncorr_ijmTmp(mm) = amplificationFactor(jj) * (vHelpJ * (ExpTerm .* mBioOffset) * vHelpM);
      
   end
   
   lengthBx = numel(bixelIxDepth);
   mIdx     = vCandidates(bixelIxDepth)';
   jIdx     = j.*ones(length(mIdx),1);
   
   PSI_y_corr_ijm(ixSave:ixSave+lengthBx -1)   = PSI_y_corr_ijmTmp;
   PSI_y_uncorr_ijm(ixSave:ixSave+lengthBx -1) = PSI_y_uncorr_ijmTmp;
   LinIx(ixSave:ixSave+lengthBx -1)            = jIdx + (mIdx-1)*numBixels;
   LinIxInv(ixSave:ixSave+lengthBx -1)         = mIdx + (jIdx-1)*numBixels;
   
   ixSave = ixSave + lengthBx;
     
end

PSI_y_corr_ijm    = PSI_y_corr_ijm(1:ixSave -1);
PSI_y_uncorr_ijm  = PSI_y_uncorr_ijm(1:ixSave -1);
LinIx             = LinIx(1:ixSave -1);
LinIxInv          = LinIxInv(1:ixSave -1);


