function [PSI_x_corr_ijm, PSI_x_uncorr_ijm, PSI_z_corr_ijm, PSI_z_uncorr_ijm, subIx] =  ...
                              matRad_calcSecLatMomFast(vCandidates,numBeams,numBixels,vBixelindexBeamOffset,vBixelindexBeam,...
                                                       Lx,Lz,latSqSigma,randError,sysError) %#codegen                                                  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the second raw moment of the lateral dose profile.
% This function is optimized for beam wise block correlation and should not
% be used for other correlatin assumptions. This function utilizes the
% following corr assumption for setup erros:
% All spots from the same beam direction are perfectly correlated whereas
% spots from different beam directions are uncorrelated.
% 
% call
% [PSI_x_corr_ijm, PSI_x_uncorr_ijm, PSI_z_corr_ijm, PSI_z_uncorr_ijm, subIx] =  ...
%                               matRad_calcSecLatMomFast(vCandidates,numBeams,numBixels,vBixelindexBeamOffset,vBixelindexBeam,...
%                                                        Lx,Lz,latSqSigma,randError,sysError,randCovLateral,sysCovLateral)%
% input
%   vCandidates:              lateral correlated bixel indices 
%   numBeams:                 number of beams
%   numBixels:                number of bixels
%   vBixelindexBeamOffset:    start and end bixel indices of each beam direction
%   vBixelindexBeam:          bixel beam look up table
%   Lx:                       lateral distance in x direction of voxel i to spot j
%   Lz:                       lateral distance in z direction of voxel i to spot j
%   latSqSigma:               beam spread (sigma of lateral beam profile in the patient considering the inital beam width
%   randError:                random setup error
%   sysError:                 systematic setup error
%   randCovLateral:           (for further developments to include arbritrary correlation assumptions) 
%   sysCovLateral:            (for further developments to include arbritrary correlation assumptions) 
%
% output
%   PSI_x_corr_ijm:            second central moment of the lateral dose in x assuming correlation (single fraction term) 
%   PSI_x_uncorr_ijm:          second central moment of the lateral dose in x assuming  no random correlation (only systematic) 
%   PSI_z_corr_ijm:            second central moment of the lateral dose in x assuming correlation (single fraction term) 
%   PSI_z_uncorr_ijm:          second central moment of the lateral dose in x assuming no random correlation correlation (only systematic) 
%   subIx                      linear indices
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
                                                    
                                                                                                                                                    
assert(isa(vCandidates,'double'));
assert(isa(numBixels,'double'));
assert(isa(vBixelindexBeamOffset,'double'));  
assert(isa(vBixelindexBeam,'double'));   
assert(isa(Lx,'double'));  
assert(isa(Lz,'double')); 
assert(isa(latSqSigma,'double'));  
assert(isa(randError,'double'));  
assert(isa(sysError,'double'));  
assert(isa(numBeams,'double'));  
% assert(isa(randCovLateral,'double'));  
% assert(isa(sysCovLateral,'double'));  

maxNumBixel   = 30000;
uppperBound2  = 10000;

coder.varsize('vCandidates',[1 maxNumBixel],[0 1]);
coder.varsize('numBixels',[1 1]);
coder.varsize('numBeams',[1 1]);
coder.varsize('vBixelindexBeamOffset',[1 uppperBound2],[0 1]);
coder.varsize('vBixelindexBeam',[maxNumBixel 1],[1 0]);
coder.varsize('Lx',[1 maxNumBixel],[0 1]);
coder.varsize('Lz',[1 maxNumBixel],[0 1]);
coder.varsize('latSqSigma',[uppperBound2 1],[1 0]);
coder.varsize('randError',[1 1],[0 0]);
coder.varsize('sysError',[1 1],[0 0]);
% coder.varsize('randCovLateral',[maxNumBixel 1],[1 0]);
% coder.varsize('sysCovLateral', [maxNumBixel 1],[1 0]);

numCandidates        =  numel(vCandidates);
ixHelper             =  (1:1:numCandidates)'; 
estimatedElemLateral = round(3*((numCandidates^2)./numBeams));
PSI_x_uncorr_ijm              = zeros(estimatedElemLateral,1);
PSI_x_corr_ijm                = zeros(estimatedElemLateral,1);
PSI_z_uncorr_ijm              = zeros(estimatedElemLateral,1);
PSI_z_corr_ijm                = zeros(estimatedElemLateral,1);
 
subIx = zeros(estimatedElemLateral,1); 
index = 1; cntSpot = 1;

for  i = 1:numCandidates

      j = vCandidates(1,i);
      upperBeamIdx = vBixelindexBeamOffset(vBixelindexBeam(j)+1);
      lowerBeamIdx = vBixelindexBeamOffset(vBixelindexBeam(j));

      bLatSqSig                = false(1,numCandidates);
      ixBixelLatFull           = vCandidates(cntSpot:end);
      ixx                      = (ixBixelLatFull <= upperBeamIdx & ixBixelLatFull >= lowerBeamIdx);
      ixBixelLat               = ixBixelLatFull(ixx);
      bLatSqSig(1,cntSpot:end) = ixx;

      if ~isempty(ixBixelLat)

         mIdx              = false(numCandidates,1); 
         mIdx(cntSpot:end) = ixx;
         mIdx              = ixHelper(mIdx);
         currSize          = length(mIdx);
         jIdx              = cntSpot.*ones(currSize,1);

         LinIdx    = jIdx + (mIdx-1)*numCandidates;
         LinIdxInv = mIdx + (jIdx-1)*numCandidates; 
         
         idxSpotBeamLUT = vBixelindexBeam(j) == vBixelindexBeam(ixBixelLat);

         vLaSi11       = (randError^2 + sysError^2) + latSqSigma(cntSpot);     
         vLaSi12       = (randError^2 + sysError^2) * idxSpotBeamLUT;
         vLaSi12uncorr = (sysError^2  * idxSpotBeamLUT);
         vLaSi22       =  sysError^2 + randError^2 + latSqSigma(bLatSqSig); 

         % x-direction
         Dev_ij = Lx(cntSpot);  
         Dev_im = Lx(mIdx)'; 
         
         DetCorr     = abs((vLaSi11 .* vLaSi22) - (vLaSi12 .* vLaSi12));
         FracDetCorr = (1./(2*pi.*real(sqrt(DetCorr))));
         
         PSI_x_corr_ijmTmp     = FracDetCorr .* exp(-.5 .* ((Dev_ij.*(vLaSi22./DetCorr) + Dev_im.*(-vLaSi12./DetCorr)).* Dev_ij +...
                                                           (Dev_ij.*(-vLaSi12./DetCorr) + Dev_im.*(vLaSi11./DetCorr)).* Dev_im));
         
         DetUnCorr     = abs((vLaSi11 .* vLaSi22) - (vLaSi12uncorr .* vLaSi12uncorr));
         FracDetUnCorr = (1./(2*pi.*real(sqrt(DetUnCorr))));
         
         PSI_x_uncorr_ijmTmp       = FracDetUnCorr .* exp(-.5 .* ((Dev_ij.*(vLaSi22./DetUnCorr) + Dev_im.*(-vLaSi12uncorr./DetUnCorr)).* Dev_ij +...
                                                                  (Dev_ij.*(-vLaSi12uncorr./DetUnCorr) + Dev_im.*(vLaSi11./DetUnCorr)).* Dev_im));

         % z-direction
         Dev_ij = Lz(cntSpot);  
         Dev_im = Lz(mIdx)'; 
         
         ExpTerm     = -.5 .* ((Dev_ij.*(vLaSi22./DetCorr) + Dev_im.*(-vLaSi12./DetCorr)).* Dev_ij +...
                               (Dev_ij.*(-vLaSi12./DetCorr) + Dev_im.*(vLaSi11./DetCorr)).* Dev_im);

         PSI_z_corr_ijmTmp    = FracDetCorr .* exp(ExpTerm);

         
         ExpTerm     = -.5 .* ((Dev_ij.*(vLaSi22./DetUnCorr) + Dev_im.*(-vLaSi12uncorr./DetUnCorr)).* Dev_ij +...
                               (Dev_ij.*(-vLaSi12uncorr./DetUnCorr) + Dev_im.*(vLaSi11./DetUnCorr)).* Dev_im);

         PSI_z_uncorr_ijmTmp  = FracDetUnCorr .* exp(ExpTerm);
         
         if numel(PSI_x_corr_ijmTmp) ~= 1
             offset                    = 2;
             vRange                    = index:index+numel(PSI_x_corr_ijmTmp)*2-offset;
             PSI_x_corr_ijm(vRange')   = [PSI_x_corr_ijmTmp; PSI_x_corr_ijmTmp(offset:end)];
             PSI_x_uncorr_ijm(vRange') = [PSI_x_uncorr_ijmTmp; PSI_x_uncorr_ijmTmp(offset:end)];
             PSI_z_corr_ijm(vRange')   = [PSI_z_corr_ijmTmp; PSI_z_corr_ijmTmp(offset:end)];
             PSI_z_uncorr_ijm(vRange') = [PSI_z_uncorr_ijmTmp; PSI_z_uncorr_ijmTmp(offset:end)];
             subIx(vRange)    = [LinIdx; LinIdxInv(offset:end)];
         else
             vRange           = index;
             PSI_x_corr_ijm(vRange)    = PSI_x_corr_ijmTmp;
             PSI_x_uncorr_ijm(vRange)  = PSI_x_uncorr_ijmTmp;
             PSI_z_corr_ijm(vRange)    = PSI_z_corr_ijmTmp;
             PSI_z_uncorr_ijm(vRange)  = PSI_z_uncorr_ijmTmp;
             subIx(vRange)    = LinIdx;
         end

         index         = index + numel(vRange);
         
      end

      cntSpot      = cntSpot + 1;    

end

PSI_x_corr_ijm    = PSI_x_corr_ijm(1:index-1,1);
PSI_x_uncorr_ijm  = PSI_x_uncorr_ijm(1:index-1,1);

PSI_z_corr_ijm    = PSI_z_corr_ijm(1:index-1,1);
PSI_z_uncorr_ijm  = PSI_z_uncorr_ijm(1:index-1,1);

subIx    = subIx(1:index-1,1); 

end



