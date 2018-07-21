function  [cst,resultGUI] =  matRad_calcParticleVarParallel(ct,cst,stf,pln,dij,resultGUI,ix,param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad's  seriel variance calculation function to compute based on APM
% the second central moment
%
% call
%   [cst,resultGUI] =  matRad_calcParticleVarSeriel(cst,stf,pln,dij,resultGUI,ix)
%
% input
%   cst:              matRad critical structure struct cst
%   stf:              matRad steering information struct
%   pln:              matRad plan meta information struct
%   dij:              matRad dij dose influence struct
%   resultGUI:        matRads resultGUI struct
%   ix                voxel indices for which variance should be calcualted
%
% output
%   cst:              matRad critical structure struct cst
%   resultGUI:        matRads resultGUI struct
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath([fileparts(mfilename('fullpath')) filesep 'MEX']);

SGGauss  = @(x,mu,SqSigma) 1./(sqrt(2*pi.*SqSigma)).*exp(-((x - mu).^2./(2.*SqSigma)));
DBGauss  = @(x,mu,SqSigma1,SqSigma2,W) ((1-W).*(1./(sqrt(2*pi.*SqSigma1)).*exp(-((x - mu).^2./(2.*SqSigma1))))) + ...
   ((W).*(1./(sqrt(2*pi.*SqSigma2)).*exp(-((x - mu).^2./(2.*SqSigma2))))) ;

sumGauss = @(x,mu,SqSigma,w) ((1./sqrt(2*pi*ones(numel(x),1) * SqSigma') .* ...
   exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) * SqSigma' ))) * w);

USE_MEX_LAT   = false;
USE_MEX_DEPTH = false;

if exist(['matRad_calcSecLatMom_mex.' mexext],'file') == 3
   USE_MEX_LAT = true;
end

if exist(['matRad_calcSecRangeMom_mex.' mexext],'file') == 3
   USE_MEX_DEPTH = true;
end

% check which depth dose component should be used for variance calculation
if ~pln.bioParam.bioOpt
   depthComp = 'Z';
else
   depthComp = 'alphaDose';
end

% calculate rED or rSP from HU
ct = matRad_calcWaterEqD(ct, pln, param);

% check if robust pencil beam weights are available
robSuffix = 'Rob';
if param.CALC_OMEGA
   robSuffix = '';
end

% find target voxel indices
omegaVoxel = [];
if pln.probOpt.omegaTarget
   for i=1:size(cst,1)
      if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
         omegaVoxel = [omegaVoxel;vertcat(cst{i,4}{:})];
      end
   end
else
   omegaVoxel = ix;
end

%% load base data set and compute reduced components
load([pln.radiationMode '_' pln.machine '.mat']);
SAD     = machine.meta.SAD;
[mMean, mWeight, mWidth, mLatSigma, mLUT, mDepth, vOffset, vCutOff] = matRad_getBaseDataAPM(ct,stf,pln,dij.vTissueIndex,depthComp);

% define biological uncertainties
if ~pln.bioParam.bioOpt
   relUCTbio = 0.0;
else
   relUCTbio = pln.multScen.relBioUCT;
end

mCovSpot = pln.multScen.mCovSpot;

doubleGauss = false;
if ~isfield(machine.data(1),'sigma')
   doubleGauss = true;
end

clear  machine;

% define some variables
stdSingleFrac        = zeros(dij.dimensions);
stdTotFrac           = zeros(dij.dimensions);
numBixels            = dij.totalNumOfBixels;
numBeams             = length([pln.propStf.gantryAngles]);
numRays              = dij.totalNumOfRays;
cubeDim              = size(resultGUI.(['physicalDose' robSuffix]));
numFrac              = pln.numOfFractions;
isoCenter            = [stf(1).isoCenter];
CTres                = [ct.resolution.x ct.resolution.y ct.resolution.z];
w                    = resultGUI.(['w' robSuffix]);
conversionFactorSq   = 1.6021766208e-02^2;
vBeamNum             = dij.beamNum;
vSigmaIni_sq         = dij.sigmaIni_sq';
vEnergyIx            = dij.energyIx';
includeLateralUTC    = (sum(strcmp(pln.probOpt.InputUCT,{'phys','biophys'})) > 0 || (strcmp(robSuffix,'Rob') && strcmp(pln.probOpt.InputUCT,{'bio'})));

% set overlap priorities
[cst,voiIndexCube] = matRad_setOverlapPriorities(cst,dij.dimensions);

% loop over all affected VOIs
% check how many and which structures are affected by std calculation
numOmega = 0;
lutVOI = 0;
for i = 1:size(cst,1)
   for j = 1:size(cst{i,6},1)
      if ~isempty(cst{i,6}(j))
         numOmega            = numOmega + 1;
         lutVOI(numOmega) = i; %#ok<AGROW>
      end
   end
end

totInfo = whos;
fprintf([' current memory consumption ' num2str((sum([totInfo.bytes]))/1e9) ' GB \n'])
dijInfo = whos('dij');
plnInfo = whos('pln');
fprintf(['dij: ' num2str(dijInfo.bytes/1e9) ' GB; pln: ' num2str(plnInfo.bytes/1e9) ' GB \n']);fprintf('');


for k = 1:numOmega
   
   mOmega          = sparse(numBixels,numBixels);
   [boolMember,~]  = ismember(ix,cst{lutVOI(k),4}{1});
   ixSub           = ix(boolMember);
   [~,ixMap]       = ismember(ixSub,ix);
   
   % initialize depth dose vector
   if ~pln.bioParam.bioOpt
      dZ         = dij.dZ{1}(ixSub,:);                           % depth doseh component including the the conversion factor
      doseExp    = resultGUI.([pln.bioParam.quantityVis 'Exp' robSuffix])(ixSub);
   else
      dZ         = dij.dZa{1}(ixSub,:);                          % alpha dose depth component
      dZb        = dij.dZb{1}(ixSub,:);                          % sqrt beta  dose depth component
      doseExp    = resultGUI.(['effectExp' robSuffix])(ixSub);
   end
   
   SpotLUT          = dij.SpotLUT{1}(ixSub,:);
   radDepth         = dij.radDepth{1}(ixSub,:);
   stdSingleFracTmp = zeros(numel(ixMap),1);
   stdTotFracTmp    = zeros(numel(ixMap),1);
   
   RBEsq = 1;
   if strcmp(pln.bioParam.model,'constRBE')
      dZ    = dZ * pln.bioParam.RBE;
      RBEsq = pln.bioParam.RBE^2;
   end
   
   %% display complexity
   fprintf(['in total are ' num2str(length(ixSub)) ' voxels affected \n']);fprintf('');
  
   % loop over all voxels and calculate the secand central moment
   %splitapply https://de.mathworks.com/matlabcentral/answers/368003-inevitable-broadcast-variable-in-parfor
   
   % running parfor with an additional parameter NumWorker to switch
   % between a serieal and parallel computation does not work. Moreover,
   % parfor with 0 works takes in simplied example twice as long as a
   % simple for loop
   for i = 1:numel(ixSub)
      
      cntVoxel = i;
      
      tmpCandidates = round(SpotLUT(cntVoxel,:));                   % for the prostate case there was one single non-integer value
      vCandidates   = full(tmpCandidates(tmpCandidates~=0));        % maybe directly store these values
      numCandidates = length(vCandidates);
      numCorr       = numCandidates^2;
      
      % in case of no correlation jump to next voxel
      if isempty(vCandidates)
         cntVoxel = cntVoxel + 1;
         continue
      end
      
      % compute lateral dose in x and z
      [xDist,zDist] = matRad_getLateralDistVoxelSpots(ixSub(i),cubeDim,CTres,isoCenter,stf,pln,vCandidates,SAD);
      
      radDepthVoxel = single(full(radDepth(cntVoxel,vBeamNum(vCandidates))));
      sqSigmaInI    = vSigmaIni_sq(vCandidates);
      subEnergyIx   = vEnergyIx(vCandidates);
      
      if doubleGauss
         latSqSigma = zeros(numCandidates,3);
      else
         latSqSigma = zeros(numCandidates,1);
      end
      
      % calc sigmas
      for jj = 1:numCandidates
         if doubleGauss
            latSqSigma(jj,:) = matRad_interp1(mDepth(:,subEnergyIx(jj)), mLatSigma(:,subEnergyIx(jj),:),radDepthVoxel(jj));
            
            latSqSigma(jj,1) = latSqSigma(jj,1).^2 + sqSigmaInI(jj);
            latSqSigma(jj,3) = latSqSigma(jj,3).^2 + sqSigmaInI(jj);
            vW               = latSqSigma(:,2);
         else
            latSqSigma(jj) = matRad_interp1(mDepth(:,subEnergyIx(jj)), mLatSigma(:,subEnergyIx(jj)),radDepthVoxel(jj)).^2  + sqSigmaInI(jj);
         end
      end
      
      % get lateral errors to calcualte the expected lateral dose
      shiftErrorLatExp = pln.multScen.shiftSDrnd.^2 + pln.multScen.shiftSDsys.^2;
      
      % calc lateral dose
      if doubleGauss
         
         if includeLateralUTC
            SqSigmaX_Narr = latSqSigma(:,1) + shiftErrorLatExp;
            SqSigmaZ_Narr = latSqSigma(:,1) + shiftErrorLatExp;
            SqSigmaX_Bro  = latSqSigma(:,3) + shiftErrorLatExp;
            SqSigmaZ_Bro  = latSqSigma(:,3) + shiftErrorLatExp;
         else
            SqSigmaX_Narr = latSqSigma(:,1);
            SqSigmaZ_Narr = latSqSigma(:,1);
            SqSigmaX_Bro  = latSqSigma(:,3);
            SqSigmaZ_Bro  = latSqSigma(:,3);
         end
         
        latSqSigmaError = SqSigmaX_Narr;

%         rSq = sqrt(xDist.^2 + zDist.^2);
%         tmp_dLx_full = sqrt(((1-vW) .* exp((-(rSq.^2))./(2*SqSigmaX_Narr))./(sqrt(SqSigmaX_Narr).*sqrt(SqSigmaZ_Narr).*2*pi)) +...
%                                (vW  .* exp((-(rSq.^2))./(2*SqSigmaX_Bro)) ./(sqrt(SqSigmaX_Bro) .*sqrt(SqSigmaZ_Bro) .*2*pi)));
%         tmp_dLz_full = sqrt(((1-vW) .* exp((-(rSq.^2))./(2*SqSigmaZ_Narr))./(sqrt(SqSigmaX_Narr).*sqrt(SqSigmaZ_Narr).*2*pi)) +...
%                                (vW  .* exp((-(rSq.^2))./(2*SqSigmaZ_Bro)) ./(sqrt(SqSigmaX_Bro) .*sqrt(SqSigmaZ_Bro) .*2*pi)));
 
        dLx = SGGauss(xDist,0,SqSigmaX_Narr);
        dLz = SGGauss(zDist,0,SqSigmaZ_Narr);
         
      else
         if includeLateralUTC
            latSqSigmaError = latSqSigma + shiftErrorLatExp;
         else
            latSqSigmaError = latSqSigma;
         end
         
         dLx = SGGauss(xDist,0,latSqSigmaError);
         dLz = SGGauss(zDist,0,latSqSigmaError);
      end
      
      vdLxCorr    = reshape(dLx*dLx',[numCorr 1]);     % create squared lateral dose in x vector
      vdLxUnCorr  = vdLxCorr;
      vdLzCorr    = reshape(dLz*dLz',[numCorr 1]) ;    % create squared lateral dose in z vector
      vdLzUnCorr  = vdLzCorr;
      
      PSI_Z_Corr   = (dZ(cntVoxel,:)' * dZ(cntVoxel,:));
      PSI_Z_Uncorr = PSI_Z_Corr;
      
      %% start lateral correlation computation
      
      % get indices
      vI    = repmat(vCandidates',[numCandidates 1]);
      vJ    = zeros(numCorr,1); loopIx = 1;
      LinIx = zeros(numCandidates);
      
      for h = 1:numCandidates
         vJ(loopIx:loopIx+numCandidates-1) = vCandidates(h);
         loopIx                            = loopIx + numCandidates;
         LinIx(:,h)                        = vCandidates' .* ones(numCandidates,1) + (vCandidates(h)'-1)*numBixels ;
      end
      
      linIxDiag = 1:numBixels+1:numBixels^2;
      
      if (pln.multScen.shiftSDrnd > 0 || pln.multScen.shiftSDsys > 0) && includeLateralUTC
         if USE_MEX_LAT
            [vPSI_XcorrTmp, vPSI_XuncorrTmp, vPSI_ZcorrTmp, vPSI_ZuncorrTmp, LinIdx] = ...
               matRad_calcSecLatMomFast_mex(...
               vCandidates,numBeams,numBixels,pln.multScen.bixelIndexBeamOffset,pln.multScen.bixelBeamLUT,xDist',zDist',latSqSigma(:,1),...
               pln.multScen.shiftSDrnd,pln.multScen.shiftSDsys);
         else
            [vPSI_XcorrTmp, vPSI_XuncorrTmp, vPSI_ZcorrTmp, vPSI_ZuncorrTmp, LinIdx] = ...
               matRad_calcSecLatMomFast(...
               vCandidates,numBeams,numBixels,pln.multScen.bixelIndexBeamOffset,pln.multScen.bixelBeamLUT,xDist',zDist',latSqSigma(:,1),...
               pln.multScen.shiftSDrnd,pln.multScen.shiftSDsys);
         end
         
      else
         if USE_MEX_LAT
            [vPSI_XcorrTmp, vPSI_XuncorrTmp, vPSI_ZcorrTmp, vPSI_ZuncorrTmp, LinIdx] = ...
               matRad_calcSecLatMomFast_mex(...
               vCandidates,numBeams,numBixels,pln.multScen.bixelIndexBeamOffset,pln.multScen.bixelBeamLUT,xDist',zDist',latSqSigma(:,1),0,0);
         else
            [vPSI_XcorrTmp, vPSI_XuncorrTmp, vPSI_ZcorrTmp, vPSI_ZuncorrTmp, LinIdx] = ...
               matRad_calcSecLatMomFast(...
               vCandidates,numBeams,numBixels,pln.multScen.bixelIndexBeamOffset,pln.multScen.bixelBeamLUT,xDist',zDist',latSqSigma(:,1),0,0);
         end
         
      end
      
       % create  lateral correlation matrices for  voxel i
      vdLxCorr(LinIdx)    = vPSI_XcorrTmp;
      PSI_Lx_Corr         = sparse(vI,vJ,vdLxCorr,numBixels,numBixels);
      vdLxUnCorr(LinIdx)  = vPSI_XuncorrTmp;
      PSI_Lx_Uncorr       = sparse(vI,vJ,vdLxUnCorr,numBixels,numBixels);
      
      
      vdLzCorr(LinIdx)    = vPSI_ZcorrTmp;
      PSI_Lz_Corr         = sparse(vI,vJ,vdLzCorr,numBixels,numBixels);
      vdLzUnCorr(LinIdx)  = vPSI_ZuncorrTmp;
      PSI_Lz_Uncorr       = sparse(vI,vJ,vdLzUnCorr,numBixels,numBixels);
      
      if (pln.multScen.rangeSDrnd > 0 || pln.multScen.rangeSDsys > 0) || relUCTbio > 0
         
         numCandidates = numel(vCandidates);
         LinIx         = zeros(numCandidates);
         vCompFac      = zeros(numCandidates,1);
         
         for h = 1:numCandidates
            jIdx        = vCandidates';
            mIdx        = vCandidates(h).*ones(length(jIdx),1);
            LinIx(:,h)  = jIdx + (mIdx-1)*numBixels;
            vCompFac(h) = vCutOff{dij.energyIx(vCandidates(h)),dij.beamNum(vCandidates(h))}.CompFac^2;
         end
         
         randCovRadDepth = full(pln.multScen.mCovRangeRnd(LinIx));
         sysCovRadDepth  = full(pln.multScen.mCovRangeSys(LinIx));
         bioCov          = full(pln.multScen.mCovBio(LinIx));
         
         %% check diagonal double indices
         if USE_MEX_DEPTH
            [PSI_CorrTmp, PSI_UncorrTmp,LinIx,LinIxInv] =  matRad_calcSecRangeMomFast_mex(conversionFactorSq * vCompFac' * RBEsq,...
               numBixels,numRays,vEnergyIx',vBeamNum',vCandidates,full(radDepth(cntVoxel,:)),mMean,mWeight,mWidth...
               ,mDepth,mLUT,pln.multScen.rangeSDrnd,pln.multScen.rangeSDsys,randCovRadDepth,sysCovRadDepth,bioCov,...
               relUCTbio,mCovSpot);
         else
            [PSI_CorrTmp, PSI_UncorrTmp,LinIx,LinIxInv] =  matRad_calcSecRangeMomFast(conversionFactorSq * vCompFac * RBEsq,...
               numBixels,numRays,vEnergyIx',vBeamNum',vCandidates,full(radDepth(cntVoxel,:)),mMean,mWeight,mWidth...
               ,mDepth,mLUT,pln.multScen.rangeSDrnd,pln.multScen.rangeSDsys,randCovRadDepth,sysCovRadDepth,bioCov,...
               relUCTbio,mCovSpot);
         end
         
         PSI_Z_Corr(LinIx)      = PSI_CorrTmp;
         PSI_Z_Corr(LinIxInv)   = PSI_CorrTmp;
         PSI_Z_Uncorr(LinIx)    = PSI_UncorrTmp;
         PSI_Z_Uncorr(LinIxInv) = PSI_UncorrTmp;
         
      end
      
      % PSI_Z_Corr either contains physical quantites or the alpha dose component
      mPSI_Corr    = (PSI_Lx_Corr   .* PSI_Lz_Corr   .* PSI_Z_Corr);
      mPSI_UnCorr  = (PSI_Lx_Uncorr .* PSI_Lz_Uncorr .* PSI_Z_Uncorr);
      PSI_Corr     = w' * (mPSI_Corr)   * w;
      PSI_UnCorr   = w' * (mPSI_UnCorr) * w;
      
      dLxSparse = sparse(ones(numCandidates,1),vCandidates,dLx,1,numBixels);
      dLzSparse = sparse(ones(numCandidates,1),vCandidates,dLz,1,numBixels);
      
      doseExpij = (dLxSparse .* dLzSparse .* dZ(cntVoxel,:));
      
      Var_Corr      = (PSI_Corr   - (doseExpij*w)^2);
      Var_UnCorr    = (PSI_UnCorr - (doseExpij*w)^2);
      Var_FracTot   = (1/numFrac) .* (Var_Corr + (Var_UnCorr * (numFrac-1)));
      
      if sum(isnan(Var_Corr)) > 0
         matRad_dispToConsole('matRad_calcParticleVarSeriel: NaN values in std',[],'warning');
      end
      
      stdTotFracTmp(i) = sqrt(Var_FracTot);
      if imag(stdTotFracTmp(i))
         stdTotFracTmp(i) = 0;
      end
      stdSingleFracTmp(i)    = sqrt(Var_Corr);
      if imag(stdSingleFracTmp(i))
         stdSingleFracTmp(i) = 0;
      end
      
      % fill omega matrix
      if param.CALC_OMEGA && ismember(ixSub(i),omegaVoxel)
         
         doseExpij2 = doseExpij' * doseExpij;      
         PSI  = (1/numFrac) .* (mPSI_Corr + (mPSI_UnCorr * (numFrac-1)));
         mOmega     = mOmega + (PSI - doseExpij2);
         
      end
      
      matRad_progress(cntVoxel,numel(ixSub));
      
   end % eof of voxel loop

   stdTotFrac(ix(ixMap))    = stdTotFracTmp;
   stdSingleFrac(ix(ixMap)) = stdSingleFracTmp;
   
   cst{lutVOI(k),6}.mOmega = mOmega;
   
end   % THIS needs to go to the end  of the voxel loop


if ~pln.bioParam.bioOpt
   
   resultGUI.([pln.bioParam.quantityVis 'StdSingleFrac' robSuffix]) = stdSingleFrac;
   resultGUI.([pln.bioParam.quantityVis 'StdTotFrac' robSuffix])    = stdTotFrac;
   
else
   
   if ~strcmp(pln.bioParam.quantityOpt,'effect')
      matRad_dispToConsole('matRad_calcParticleVarSerial: not supported',[],'error')
   end
   
   resultGUI.([pln.bioParam.quantityOpt 'StdSingleFrac' robSuffix]) = stdSingleFrac;
   resultGUI.([pln.bioParam.quantityOpt 'StdTotFrac' robSuffix])    = stdTotFrac;
   
   [~,vStdRBExDsingleFrac]          = matRad_effectToRBExDose(resultGUI.(['effectExp' robSuffix]), stdSingleFrac ,dij.alphaX,dij.betaX);
   [vMeanRBExD,vStdRBExDtotFrac]    = matRad_effectToRBExDose(resultGUI.(['effectExp' robSuffix]), stdTotFrac    ,dij.alphaX,dij.betaX);
   
   resultGUI.([pln.bioParam.quantityVis 'Exp'  robSuffix])             = reshape(vMeanRBExD,dij.dimensions);
   resultGUI.([pln.bioParam.quantityVis 'StdSingleFrac'  robSuffix])   = reshape(vStdRBExDsingleFrac,dij.dimensions);
   resultGUI.([pln.bioParam.quantityVis 'StdTotFrac'  robSuffix])      = reshape(vStdRBExDtotFrac,dij.dimensions);
   
end



end % eof function

