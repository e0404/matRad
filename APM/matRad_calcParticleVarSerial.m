function  [cst,resultGUI] =  matRad_calcParticleVarSerial(ct,cst,stf,pln,dij,resultGUI,ix,param)
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
numComp = numel(machine.data(1).Z.width);
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

%clear machine;                                  
% define some variables
stdSingleFrac        = zeros(dij.dimensions);
stdTotFrac           = zeros(dij.dimensions);
numBixels            = dij.totalNumOfBixels;
cubeDim              = size(resultGUI.(['physicalDose' robSuffix]));
numFrac              = pln.numOfFractions;
isoCenter            = [stf(1).isoCenter];
CTres                = [ct.resolution.x ct.resolution.y ct.resolution.z];
radDepth             = dij.radDepth{1}(ix,:);
w                    = resultGUI.(['w' robSuffix]);
conversionFactor     = 1.6021766208e-02;   
includeLateralUTC    = (sum(strcmp(pln.probOpt.InputUCT,{'phys','biophys'})) > 0 || (strcmp(robSuffix,'Rob') && strcmp(pln.probOpt.InputUCT,{'bio'})));

% set overlap priorities
[cst,voiIndexCube]   = matRad_setOverlapPriorities(cst,dij.dimensions);

% initialize depth dose vector
dZb        = zeros(numel(ix),1); 
if ~pln.bioParam.bioOpt
   dZ         = dij.dZ{1}(ix,:);                           % depth doseh component including the the conversion factor
   doseExp    = resultGUI.([pln.bioParam.quantityVis 'Exp' robSuffix])(ix);
else 
   dZ         = dij.dZa{1}(ix,:);                          % alpha dose depth component   
   dZb        = dij.dZb{1}(ix,:);                          % sqrt beta  dose depth component                                                
   doseExp    = resultGUI.(['effectExp' robSuffix])(ix);  
end

RBEsq = 1;
if strcmp(pln.bioParam.model,'constRBE')
   dZ    = dZ * pln.bioParam.RBE;
   RBEsq = pln.bioParam.RBE^2;
end

%% display complexity
totInfo = whos;
fprintf([' current memory consumption ' num2str((sum([totInfo.bytes]))/1e9) ' GB \n'])
dijInfo = whos('dij');
plnInfo = whos('pln');
dZInfo  = whos('dZ');
matRad_dispToConsole(['dij: ' num2str(dijInfo.bytes/1e9) ' GB; pln: ' num2str(plnInfo.bytes/1e9) ' GB;' ' dZ: ' num2str(dZInfo.bytes/1e9) ' GB \n'],param,'info')
matRad_dispToConsole(['In total are ' num2str(length(ix)) ' voxels affected \n'],param,'info');

cntVoxel = 1;

% loop over all voxels and calculate the secand central moment
for i = ix'
    
    tmpCandidates = dij.SpotLUT{1}(i,:);
    vCandidates   = full(tmpCandidates(tmpCandidates~=0));                  % maybe directly store these values
    numCandidates = length(vCandidates);
    
    % in case of no correlation jump to next voxel
     if isempty(vCandidates)
        cntVoxel = cntVoxel + 1;
        continue
     end
     
     % compute lateral dose in x and z
     [xDist,zDist] = matRad_getLateralDistVoxelSpots(i,cubeDim,CTres,isoCenter,stf,pln,vCandidates,SAD);
     
     radDepthVoxel = full(radDepth(cntVoxel,dij.beamNum(vCandidates))) + 0;
     sqSigmaInI    = dij.sigmaIni_sq(vCandidates);   
     subEnergyIx   = dij.energyIx(vCandidates);
    
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
 
        tmp_dLx = SGGauss(xDist,0,SqSigmaX_Narr);
        tmp_dLz = SGGauss(zDist,0,SqSigmaZ_Narr);
        
        dLx = sparse(ones(numCandidates,1),vCandidates,tmp_dLx,1,numBixels);
        dLz = sparse(ones(numCandidates,1),vCandidates,tmp_dLz,1,numBixels);
        
     else
        
        if includeLateralUTC
           latSqSigmaError = latSqSigma(:,1) + shiftErrorLatExp;
        else
           latSqSigmaError = latSqSigma(:,1);
        end
        
        dLx = sparse(ones(numCandidates,1),vCandidates,SGGauss(xDist,0,latSqSigmaError),1,numBixels);
        dLz = sparse(ones(numCandidates,1),vCandidates,SGGauss(zDist,0,latSqSigmaError),1,numBixels);
        
     end
     
    PSI_Lx_Corr    = (dLx' * dLx);  
    PSI_Lz_Corr    = (dLz' * dLz);
    PSI_Z_Corr     = (dZ(cntVoxel,:)' * dZ(cntVoxel,:));   
    PSI_Lx_Uncorr  = PSI_Lx_Corr;
    PSI_Lz_Uncorr  = PSI_Lz_Corr;
    PSI_Z_Uncorr   = PSI_Z_Corr;   
     
    
    %% consistency check
    if doubleGauss
       rSq = sqrt(xDist.^2 + zDist.^2);
       tmp_dLx_full = sqrt(((1-vW) .* exp((-(rSq.^2))./(2*SqSigmaX_Narr))./(sqrt(SqSigmaX_Narr).*sqrt(SqSigmaZ_Narr).*2*pi)) + ...
                              (vW  .* exp((-(rSq.^2))./(2*SqSigmaX_Bro)) ./(sqrt(SqSigmaX_Bro) .*sqrt(SqSigmaZ_Bro) .*2*pi)));
       tmp_dLz_full = sqrt(((1-vW) .* exp((-(rSq.^2))./(2*SqSigmaZ_Narr))./(sqrt(SqSigmaX_Narr).*sqrt(SqSigmaZ_Narr).*2*pi)) + ...
                              (vW  .* exp((-(rSq.^2))./(2*SqSigmaZ_Bro)) ./(sqrt(SqSigmaX_Bro) .*sqrt(SqSigmaZ_Bro) .*2*pi)));
       dLxDG = sparse(ones(numCandidates,1),vCandidates,tmp_dLx_full,1,numBixels);
       dLzDG = sparse(ones(numCandidates,1),vCandidates,tmp_dLz_full,1,numBixels);
    end
    
    if ~pln.bioParam.bioOpt
       if ~doubleGauss
          if abs( (full(dLx .* dLz .* dZ(cntVoxel,:)) * w )   - doseExp(cntVoxel) ) > 1e-4
             warning('recalculated dose does not match !');
          end
       else
          if abs( (full(dLxDG .* dLzDG .* dZ(cntVoxel,:)) * w )   - doseExp(cntVoxel) ) > 1e-4
             warning('recalculated dose does not match !');
          end
       end
    else
       if ~doubleGauss
          if abs(   full(  ((dLx .* dLz .* dZ(cntVoxel,:)) * w) + ((dLx .* dLz .* dZb(cntVoxel,:)) * w)^2) - doseExp(cntVoxel)) > 1e-4
             warning('recalculated effect does not match !');
          end
       else
          if abs(   full(  ((dLxDG .* dLzDG .* dZ(cntVoxel,:)) * w) + ((dLxDG .* dLzDG .* dZb(cntVoxel,:)) * w)^2) - doseExp(cntVoxel)) > 1e-4
             warning('recalculated effect does not match !');
          end
       end
    end
    
    % define lateral cut offs  - this is a bit overkill - cause the 
    % cut off is already defined by spots which contribute to a certain
    % voxel. However for the APM calculation it is possible to set an
    % additional cutoff.
    if  pln.probOpt.LatCutOff > 0
        lateralCutOff  = sqrt(pln.probOpt.LatCutOff.*latSqSigmaError);
        bLcutOffX      = abs(xDist) < lateralCutOff;
        bLcutOffZ      = abs(zDist) < lateralCutOff;
    end
    
    % estimate the number of elements in the PSI matrix
    iniGuess         = numBixels*numBixels;
    cntPSIx          = 1; cntPSIz    = 1;  cntPSIrange = 1;
    vPSI_Lx_i        = zeros(iniGuess,1);
    vPSI_Lx_j        = zeros(iniGuess,1);
    vPSI_Lx_lin      = zeros(iniGuess,1);
    vPSI_Lx_CorrV    = zeros(iniGuess,1);
    vPSI_Lx_UncorrV  = zeros(iniGuess,1);
    
    vPSI_Lz_i        = zeros(iniGuess,1);
    vPSI_Lz_j        = zeros(iniGuess,1);
    vPSI_Lz_lin      = zeros(iniGuess,1);
    vPSI_Lz_CorrV    = zeros(iniGuess,1);
    vPSI_Lz_UncorrV  = zeros(iniGuess,1);
    
    vPSI_Z_lin       = zeros(iniGuess,1);
    vPSI_Z_CorrV     = zeros(iniGuess,1);
    vPSI_Z_UncorrV   = zeros(iniGuess,1);
    
    cntSpot = 1;
    
    for  j = vCandidates    % only loop over bixel combinations that contribute dose to voxel i
        
         lowerBeamIdx = pln.multScen.bixelIndexBeamOffset(dij.beamNum(j));
         upperBeamIdx = pln.multScen.bixelIndexBeamOffset(dij.beamNum(j)+1);

         % find correlated spot indices in 
         bLatSqSig                = false(1,numCandidates);
         ixBixelLatFull           = vCandidates(cntSpot:end);                                          % loop through vCandidates in a syemmtric manner
         ixx                      = (ixBixelLatFull < upperBeamIdx & ixBixelLatFull >= lowerBeamIdx);  % check if all spots belong to the same beam direction
        
         ixBixelLat               = ixBixelLatFull(ixx);                                                % obtained reduced correlated spot indices
         bLatSqSig(1,cntSpot:end) = ixx;                                                                % obtain corresponding logical index vector
             
         % apply additional cut off in x direction
         ixBixelLatX = ixBixelLat';
         bLatSqSigX  = bLatSqSig;
         
         if pln.probOpt.LatCutOff > 0
             ixCut       = bLcutOffX(ixx)';
             ixBixelLatX = ixBixelLatX(ixCut);
             bLatSqSigX(bLatSqSigX==true) = ixCut;
         end
         
         if ~isempty(ixBixelLatX)
             
              % obtain indices 
              mIdx      = ixBixelLatX; 
              jIdx      = j.*ones(numel(mIdx),1);
              
              LinIdx    = jIdx + (mIdx-1) * numBixels;
              LinIdxInv = mIdx + (jIdx-1) * numBixels;
              
              Dev_ij    = xDist(cntSpot);          % lateral distance in x from voxel i to central axis of pencil beam j            
              Dev_im    = xDist(bLatSqSigX);       % lateral distance in x from voxel i to central axis of pencil beam m       
              
              ixSpotBeamLUT = dij.beamNum(j) == dij.beamNum(mIdx);

              vLaSi11  =     latSqSigma(cntSpot,1) + pln.multScen.shiftSDrnd ^2 + pln.multScen.shiftSDsys^2 ;  % combines Lambda and Sigma -> LaSi
              vLaSi22  =  latSqSigma(bLatSqSigX,1) + pln.multScen.shiftSDrnd ^2 + pln.multScen.shiftSDsys^2 ;  
              vLaSi12uncorr = (pln.multScen.shiftSDsys^2  * ixSpotBeamLUT);
              
              if (pln.multScen.shiftSDrnd > 0 || pln.multScen.shiftSDsys > 0) && includeLateralUTC
                 vLaSi12  =  (pln.multScen.shiftSDrnd ^2 + pln.multScen.shiftSDsys^2 )* ixSpotBeamLUT;
                 % vLaSi22   = full(pln.multScen.mCovLatRnd(mIdx) + latSqSigma(bLatSqSigX,1));
              else
                 vLaSi12  =  0 * ixSpotBeamLUT;
              end
              
              if USE_MEX_LAT
                 vPSI_CorrTmp   = matRad_calcSecLatMom_mex(vLaSi11,vLaSi22,vLaSi12,vLaSi12,Dev_ij,Dev_im);
                 vPSI_UncorrTmp = matRad_calcSecLatMom_mex(vLaSi11,vLaSi22,vLaSi12uncorr,vLaSi12uncorr,Dev_ij,Dev_im);
              else
                 vPSI_CorrTmp   = matRad_calcSecLatMom(vLaSi11,vLaSi22,vLaSi12,vLaSi12,Dev_ij,Dev_im);
                 vPSI_UncorrTmp = matRad_calcSecLatMom(vLaSi11,vLaSi22,vLaSi12uncorr,vLaSi12uncorr,Dev_ij,Dev_im);
              end
               
              numPSI = numel(vPSI_CorrTmp);            
              
              if numPSI ~= 1 
                  offset = 2;
                  vRangeX                  = (cntPSIx:cntPSIx+numPSI*2-offset)';
                  vPSI_Lx_i(vRangeX)       = [jIdx; mIdx(offset:end)];   % take care of diagonal element, double indices will be added up when sparse matrix is created
                  vPSI_Lx_j(vRangeX)       = [mIdx; jIdx(offset:end)];
                  vPSI_Lx_CorrV(vRangeX)   = [vPSI_CorrTmp;   vPSI_CorrTmp(offset:end)];
                  vPSI_Lx_UncorrV(vRangeX) = [vPSI_UncorrTmp; vPSI_UncorrTmp(offset:end)];
                  vPSI_Lx_lin(vRangeX)     = [LinIdx; LinIdxInv(offset:end)];   
              else
                  vRangeX = cntPSIx;
                  vPSI_Lx_i(vRangeX)       = jIdx;
                  vPSI_Lx_CorrV(vRangeX)   = vPSI_CorrTmp;
                  vPSI_Lx_UncorrV(vRangeX) = vPSI_UncorrTmp;
                  vPSI_Lx_lin(vRangeX)     = LinIdx;
              end
              
              cntPSIx  = cntPSIx  + numel(vRangeX);
              
         end % eof ~isempty(ixBixelLatX)
            
         % apply additional cut off in x direction
         ixBixelLatZ = ixBixelLat';
         bLatSqSigZ  = bLatSqSig;
         
         if pln.probOpt.LatCutOff > 0
             ixCut       = bLcutOffZ(ixx)';
             ixBixelLatZ = ixBixelLatZ(ixCut);
             bLatSqSigZ(bLatSqSigZ==true) = ixCut;
         end
         
         if ~isempty(ixBixelLatZ)
             
             % obtain indices 
              mIdx      = ixBixelLatZ; 
              jIdx      = j.*ones(numel(mIdx),1);
              LinIdx    = jIdx + (mIdx-1)* numBixels;
              LinIdxInv = mIdx + (jIdx-1)* numBixels;
              
              Dev_ij    = zDist(cntSpot);          % lateral distance in z from voxel i to central axis of pencil beam j            
              Dev_im    = zDist(bLatSqSigZ);       % lateral distance in z from voxel i to central axis of pencil beam m       
                   
              ixSpotBeamLUT = dij.beamNum(j) == dij.beamNum(mIdx);

              vLaSi11  =     latSqSigma(cntSpot,1) + pln.multScen.shiftSDrnd ^2 + pln.multScen.shiftSDsys^2 ;  % combines Lambda and Sigma -> LaSi
              vLaSi22  =  latSqSigma(bLatSqSigZ,1) + pln.multScen.shiftSDrnd ^2 + pln.multScen.shiftSDsys^2 ;  
              vLaSi12uncorr = (pln.multScen.shiftSDsys^2  * ixSpotBeamLUT);

              if (pln.multScen.shiftSDrnd > 0 || pln.multScen.shiftSDsys > 0) && includeLateralUTC
                 vLaSi12  =  (pln.multScen.shiftSDrnd ^2 + pln.multScen.shiftSDsys^2 )* ixSpotBeamLUT;
                 % vLaSi22   = full(pln.multScen.mCovLatRnd(mIdx) + latSqSigma(bLatSqSigX,1));
              else
                 vLaSi12  =  0 * ixSpotBeamLUT;
              end
              
              if USE_MEX_LAT
                 vPSI_CorrTmp   = matRad_calcSecLatMom_mex(vLaSi11,vLaSi22,vLaSi12,vLaSi12,Dev_ij,Dev_im);
                 vPSI_UncorrTmp = matRad_calcSecLatMom_mex(vLaSi11,vLaSi22,vLaSi12uncorr,vLaSi12uncorr,Dev_ij,Dev_im);
              else
                 vPSI_CorrTmp   = matRad_calcSecLatMom(vLaSi11,vLaSi22,vLaSi12,vLaSi12,Dev_ij,Dev_im);
                 vPSI_UncorrTmp = matRad_calcSecLatMom(vLaSi11,vLaSi22,vLaSi12uncorr,vLaSi12uncorr,Dev_ij,Dev_im);
              end
              
              numPSI = numel(vPSI_CorrTmp);
              
              if numPSI ~= 1
                  offset = 2;
                  vRangeZ                  = (cntPSIz:cntPSIz+numPSI*2-offset)';
                  vPSI_Lz_i(vRangeZ)       = [jIdx; mIdx(offset:end)];   % take care of diagonal element, double indices will be added during creation of sparse matrix
                  vPSI_Lz_j(vRangeZ)       = [mIdx; jIdx(offset:end)];
                  vPSI_Lz_CorrV(vRangeZ)   = [vPSI_CorrTmp; vPSI_CorrTmp(offset:end)];
                  vPSI_Lz_UncorrV(vRangeZ) = [vPSI_UncorrTmp; vPSI_UncorrTmp(offset:end)];
                  vPSI_Lz_lin(vRangeZ)     = [LinIdx; LinIdxInv(offset:end)];
              else
                  vRangeZ = cntPSIz;
                  vPSI_Lz_i(vRangeZ)       = jIdx;
                  vPSI_Lz_CorrV(vRangeZ)   = vPSI_CorrTmp;
                  vPSI_Lz_UncorrV(vRangeZ) = vPSI_UncorrTmp;
                  vPSI_Lz_lin(vRangeZ)     = LinIdx;
              end
              
             cntPSIz  = cntPSIz  + numel(vRangeZ);
              
         end
        
         % calc second raw moment of range - assuming rand and sys errors
         % have the same correlation structure
         bixelIdxDepth = [];
         if pln.multScen.rangeSDrnd > 0
             bixelIdxDepth = unique([bixelIdxDepth find((pln.multScen.mCovRangeRnd(j,j:vCandidates(end)) > 0 )) + j - 1]);
         end
         
         if pln.multScen.rangeSDsys > 0
             bixelIdxDepth = unique([bixelIdxDepth find((pln.multScen.mCovRangeSys(j,j:vCandidates(end)) > 0 )) + j - 1]);
         end
         
         if relUCTbio > 0
             bixelIdxDepth  = unique([bixelIdxDepth find((pln.multScen.mCovBio(j,j:vCandidates(end)) > 0 )) + j - 1]);
         end
             
         % obtain look up table for reduced components for spot j
          beamIx_j   = dij.beamNum(j);
          energyIx_j = dij.energyIx(j); 
          if pln.probOpt.useReducedComp 
               ixLut_j  = find(mDepth(:,energyIx_j) > radDepth(cntVoxel,beamIx_j),1,'first');  % find depth index  
               vIxLUT_j = logical(squeeze(mLUT(:,ixLut_j,energyIx_j)));
          else
               vIxLUT_j = true(numComp,1);
          end
           
          Dev_j       = radDepth(cntVoxel,beamIx_j) - mMean(vIxLUT_j,energyIx_j);  
          SigmaSq_j   = mWidth(vIxLUT_j,energyIx_j);
          Weight_j    = mWeight(vIxLUT_j,energyIx_j);
          LinInd_MM   = bixelIdxDepth'  + ((bixelIdxDepth-1)*numBixels)';
         
         
          if pln.multScen.rangeSDrnd > 0
                randCovRadDepth     = pln.multScen.mCovRangeRnd(j,j);                           
                vRandCovRadDepth    = full(pln.multScen.mCovRangeRnd(j,bixelIdxDepth))';
                vRandCovRadDepth22  = full(pln.multScen.mCovRangeRnd(LinInd_MM));
          else
                vRandCovRadDepth22 = 0;
          end
             
          
          if pln.multScen.rangeSDsys > 0
                 sysCovRadDepth      = pln.multScen.mCovRangeSys(j,j);
                 vSysCovRadDepth     = full(pln.multScen.mCovRangeSys(j,bixelIdxDepth))';
                 vSysCovRadDepth22   = full(pln.multScen.mCovRangeSys(LinInd_MM));
          else
                 vSysCovRadDepth22   = 0;
          end
         
           cntRange = 1;
           
           % loop over all correlated spot combinations
           if ~isempty(bixelIdxDepth) 
             for m = bixelIdxDepth
              
              SqSigmaRandUnCorr_jm = 0;

              if pln.multScen.rangeSDrnd > 0
                  SqSigmaRandCorr_jm   = [randCovRadDepth vRandCovRadDepth(cntRange) vRandCovRadDepth(cntRange) vRandCovRadDepth22(cntRange)];
                  SqSigmaRandUnCorr_jm = [randCovRadDepth 0  0 vRandCovRadDepth22(cntRange)];
              else
                  SqSigmaRandCorr_jm = 0;
              end

              if pln.multScen.rangeSDsys > 0
                  SqSigmaSys_jm = [sysCovRadDepth vSysCovRadDepth(cntRange)  vSysCovRadDepth(cntRange) vSysCovRadDepth22(cntRange)];
              else
                  SqSigmaSys_jm = 0;
              end
              
              % obtain look up table for reduced components for spot m
              beamIx_m   = dij.beamNum(m); 
              energyIx_m = dij.energyIx(m); 
              
              
              if pln.probOpt.useReducedComp
                  ixLut_m  = find(mDepth(:,energyIx_m) > radDepth(cntVoxel,beamIx_m),1,'first');   %this line of code is still slow
                  vIxLUT_m = (mLUT(:,ixLut_m,energyIx_m));
              else
                  vIxLUT_m = true(numComp,1);
              end
              
              
              Weight_m    = mWeight(vIxLUT_m,energyIx_m);
              Dev_m       = radDepth(cntVoxel,beamIx_m) - mMean(vIxLUT_m,energyIx_m);
              SigmaCorr   = full(SqSigmaRandCorr_jm   + SqSigmaSys_jm);
              SigmaUnCorr = full(SqSigmaRandUnCorr_jm + SqSigmaSys_jm);


              compFac     = vCutOff{energyIx_j,dij.beamNum(j)}.CompFac  * vCutOff{energyIx_m,dij.beamNum(m)}.CompFac;
              
              %mW_CovAlphaDose_jm  = caWeightUCT.mCovalphaDose(j*numComp-(numComp-1):j*numComp,m*numComp-(numComp-1):m*numComp);
              mW_CovAlphaDose_jm = mCovSpot(vIxLUT_j,vIxLUT_m).* ((Weight_j*relUCTbio) *  (Weight_m'*relUCTbio));
                   
              if USE_MEX_DEPTH
                  
                  PSI_CorrTmp   = conversionFactor^2 * compFac * RBEsq *  matRad_calcSecRangeMom_mex(SigmaSq_j + SigmaCorr(1,1),...
                                                             mWidth(vIxLUT_m,energyIx_m) + SigmaCorr(:,4),...
                                                             SigmaCorr(1,2),SigmaCorr(1,3),Dev_j,Dev_m,...
                                                             Weight_j, Weight_m, mW_CovAlphaDose_jm);   

                  PSI_UncorrTmp = conversionFactor^2 * compFac * RBEsq * matRad_calcSecRangeMom_mex(SigmaSq_j + SigmaUnCorr(1,1),...
                                                             mWidth(vIxLUT_m,energyIx_m) + SigmaUnCorr(:,4),...
                                                             SigmaUnCorr(1,2),SigmaUnCorr(1,3),Dev_j,Dev_m,...
                                                             Weight_j, Weight_m, mW_CovAlphaDose_jm);
              else
                  PSI_CorrTmp   = conversionFactor^2 * compFac * RBEsq *  matRad_calcSecRangeMom(SigmaSq_j + SigmaCorr(1,1),...
                                                             mWidth(vIxLUT_m,energyIx_m) + SigmaCorr(:,4),...
                                                             SigmaCorr(1,2),SigmaCorr(1,3),Dev_j,Dev_m,...
                                                             Weight_j, Weight_m, mW_CovAlphaDose_jm);   

                  PSI_UncorrTmp = conversionFactor^2 * compFac * RBEsq * matRad_calcSecRangeMom(SigmaSq_j + SigmaUnCorr(1,1),...
                                                             mWidth(vIxLUT_m,energyIx_m) + SigmaUnCorr(:,4),...
                                                             SigmaUnCorr(1,2),SigmaUnCorr(1,3),Dev_j,Dev_m,...
                                                             Weight_j, Weight_m, mW_CovAlphaDose_jm);
                  
              end
                  
                  
              if j == m
                  LinIdx                      = j + (m-1)* numBixels;
                  vPSI_Z_CorrV(cntPSIrange)   = PSI_CorrTmp;
                  vPSI_Z_UncorrV(cntPSIrange) = PSI_UncorrTmp;
                  vPSI_Z_lin(cntPSIrange)     = LinIdx;
                  cntPSIrange                 = cntPSIrange  + 1;  
              else
                  LinIdx                      = j + (m-1)* numBixels;
                  LinIdxInv                   = m + (j-1)* numBixels;
                  vRangeZ                     = [cntPSIrange:cntPSIrange+1]';
                  vPSI_Z_CorrV(vRangeZ)       = [PSI_CorrTmp; PSI_CorrTmp];
                  vPSI_Z_UncorrV(vRangeZ)     = [PSI_UncorrTmp; PSI_UncorrTmp];
                  vPSI_Z_lin(vRangeZ)         = [LinIdx; LinIdxInv];
                  cntPSIrange                 = cntPSIrange  + 2;  
              end
              
              cntRange     = cntRange + 1;
              
             end
          end
          cntSpot  = cntSpot + 1;
    end
  
    % fill sparse matrix from vectors
    PSI_Lx_Corr(vPSI_Lx_lin(1:cntPSIx-1))     = vPSI_Lx_CorrV(1:cntPSIx-1);
    PSI_Lz_Corr(vPSI_Lz_lin(1:cntPSIz-1))     = vPSI_Lz_CorrV(1:cntPSIz-1);
    PSI_Lx_Uncorr(vPSI_Lx_lin(1:cntPSIx-1))   = vPSI_Lx_UncorrV(1:cntPSIx-1);
    PSI_Lz_Uncorr(vPSI_Lz_lin(1:cntPSIz-1))   = vPSI_Lz_UncorrV(1:cntPSIz-1);
    PSI_Z_Corr(vPSI_Z_lin(1:cntPSIrange-1))   = vPSI_Z_CorrV(1:cntPSIrange-1);
    PSI_Z_Uncorr(vPSI_Z_lin(1:cntPSIrange-1)) = vPSI_Z_UncorrV(1:cntPSIrange-1);
    
    % PSI_Z_Corr either contains physical quantites or the alpha dose component
    mPSI_Corr    = (PSI_Lx_Corr   .* PSI_Lz_Corr   .* PSI_Z_Corr);         
    mPSI_UnCorr  = (PSI_Lx_Uncorr .* PSI_Lz_Uncorr .* PSI_Z_Uncorr); 
    PSI_Corr     = w' * (mPSI_Corr)   * w;
    PSI_UnCorr   = w' * (mPSI_UnCorr) * w; 
    
    doseExpij = (dLx .* dLz .* dZ(cntVoxel,:));
    
%     % consider alpha dose and beta dose    
%     if pln.bioParam.bioOpt
%        % create a helper vector vC
%        mII_Corr   = spalloc(numBixels,numBixels,1);
%        mII_UnCorr = spalloc(numBixels,numBixels,1);
%        mTmpCorr   = PSI_Lx_Corr   .* PSI_Lz_Corr;             
%        mTmpUnCorr = PSI_Lx_Uncorr .* PSI_Lz_Uncorr;
%        dZbTmp     = dZb(cntVoxel,:);
%        %res        = zeros(NumBixels,1);                    % uncomment for alternative solution

%       % only calculate the second term in BxB dimension when omega matrix is filled
%       if param.CALC_OMEGA
% 
%           cnt  = 1;
%           vC                   = (dLx .* dLz .* dZ(cntVoxel,:)) .* w';    % dZ is in this case dZa
%           numOfBixelsContainer = ceil(numBixels/10);
%           CorrCont             = cell(numOfBixelsContainer,1);    % correlated bixel container
%           UnCorrCont           = cell(numOfBixelsContainer,1);    % uncorrelated bixel container
%           
%           for l = 1:numBixels  
% 
%              vTmp  =  vC.*dZbTmp(l);  
%              CorrCont{mod(cnt-1,numOfBixelsContainer)+1,1}    = sum((mTmpCorr(l,:)  .*dZbTmp)' * vTmp,2);
%              UnCorrCont{mod(cnt-1,numOfBixelsContainer)+1,1}  = sum((mTmpUnCorr(l,:).*dZbTmp)' * vTmp,2);
% 
%              % alternative solution res(g) = w' * (sum((mTmpUnCorr(g,:).*dZbTmp)' * vTmp,2));
%              if mod(cnt,numOfBixelsContainer) == 0 || cnt == numBixels
%                  mII_Corr(:,(ceil(cnt/numOfBixelsContainer)-1)*numOfBixelsContainer+1:cnt)   = [CorrCont{1:mod(cnt-1,numOfBixelsContainer)+1,1}];    %#ok<SPRIX> % needed for omega matrix calculation
%                  mII_UnCorr(:,(ceil(cnt/numOfBixelsContainer)-1)*numOfBixelsContainer+1:cnt) = [UnCorrCont{1:mod(cnt-1,numOfBixelsContainer)+1,1}];  %#ok<SPRIX> % needed for omega matrix calculation
%              end
%              cnt = cnt + 1;
%           end
%           II_Corr    = w' * mII_Corr * w;
%           II_UnCorr  = w' * mII_UnCorr * w;  
%        else
%           vC  = (dLx .* dLz .* dZ(cntVoxel,:)) * w;    % dZ is in this case dZa
%           II_Corr = 0; II_UnCorr = 0; 
%           for l = 1:numBixels      
%              II_Corr   = II_Corr   + w(l) * (sum(w.*((( mTmpCorr(l,:)'   * dZbTmp(l)) .* dZbTmp') * vC)));
%              II_UnCorr = II_UnCorr + w(l) * (sum(w.*((( mTmpUnCorr(l,:)' * dZbTmp(l)) .* dZbTmp') * vC)));  
%           end 
% 
%        end  
% 
%        mIII_approxCorr    = (PSI_Lx_Corr  .* PSI_Lz_Corr   .* (dZb(cntVoxel,:)'*dZb(cntVoxel,:))) ;          % needed for omega matrix calculation
%        mIII_approxUnCorr  = (PSI_Lx_Uncorr.* PSI_Lz_Uncorr .* (dZb(cntVoxel,:)'*dZb(cntVoxel,:))) ;          % needed for omega matrix calculation
%        III_approxCorr     =  (w'* (mIII_approxCorr)   * w)^2;
%        III_approxUnCorr   =  (w'* (mIII_approxUnCorr) * w)^2;
% 
%        PSI_Corr      =  PSI_Corr   + 2*II_Corr   + III_approxCorr;
%        PSI_UnCorr    =  PSI_UnCorr + 2*II_UnCorr + III_approxUnCorr;
%
%    end
    
      Var_Corr      = (PSI_Corr   - (doseExpij*w)^2);
      Var_UnCorr    = (PSI_UnCorr - (doseExpij*w)^2);
      Var_FracTot   = (1/numFrac) .* (Var_Corr + (Var_UnCorr * (numFrac-1)));
     
    if sum(isnan(Var_Corr)) > 0 || Var_Corr < 0
        matRad_dispToConsole('matRad_calcParticleVarSeriel: NaN values in std',[],'warning');
    end  

    stdSingleFrac(i) = sqrt(Var_Corr); 
    if imag(stdSingleFrac(i))
       stdSingleFrac(i) = 0;
    end
    stdTotFrac(i) = sqrt(Var_FracTot); 
    if imag(stdTotFrac(i))
       stdTotFrac(i) = 0;
    end
        
    % fill omega matrix
    if param.CALC_OMEGA && ismember(i,omegaVoxel)
        
        doseExpij2 = doseExpij' * doseExpij;      
        PSI  = (1/numFrac) .* (mPSI_Corr + (mPSI_UnCorr * (numFrac-1)));
        cst{voiIndexCube(i),6}.mOmega = cst{voiIndexCube(i),6}.mOmega + (PSI - doseExpij2);
        
    end
      
    matRad_progress(cntVoxel,length(ix));
    cntVoxel = cntVoxel + 1;
    
end % eof of voxel loop


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

