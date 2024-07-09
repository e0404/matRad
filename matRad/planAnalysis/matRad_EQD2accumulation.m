function result = matRad_EQD2accumulation(pln1,ct1,cst1,dose1,prescribedDose1, ...
                                          pln2,ct2,cst2,dose2,prescribedDose2)
                             
% matRad function to accumulate and compare dose and EQD2 for two treatment 
% plans 
%
% call
%   result = matRad_EQD2accumulation(pln1,ct1,cst1,dose1,prescribedDose1, ...
%                                    pln2,ct2,cst2,dose2,prescribedDose2)
%
% input
%   pln1/2:             matRad pln struct
%   ct1/2:              matRad ct struct
%   cst1/2:             matRad cst struct
%   dose1/2:            3D (RBE-weighted) dose cubes
%   prescribedDose1/2:  prescribed doses of the respective dose cubes
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             
                             
% parameters that can be changed to fit to the RT plans, also alpha/beta
% for the EQD_2 calculation can be changed
 
checkImageRegis = true;                                % true or false, shows pictures of the particle and photon CT to check if the image registration went wrong (true shows pictures)

alphaBetaRatio = 2;                                     % alpha/beta value

numOfFractions1 = pln1.numOfFractions;                              % number of photon fractions
numOfFractions2 = pln2.numOfFractions;                               % number of proton fractions

%% image registration and recalculation of the photon dose cube to fit in the particle dose cube
% Rfixed/Rmoving are helping the optimizer
fixed = ct2.cubeHU{1};
Rfixed = imref3d(ct2.cubeDim,ct2.resolution.y,ct2.resolution.x,ct2.resolution.z);
moving = ct1.cubeHU{1};
Rmoving = imref3d(ct1.cubeDim,ct1.resolution.y,ct1.resolution.x,ct1.resolution.z);

% choice of the optimizer for the image registration, for this ct
% registration monomodal is highly recommended
% below are possible changes for the optimization, values are the default
% values; only change is done in the maximum iterations (default 100)

% [optimizer, metric] = imregconfig('multimodal');
[optimizer, metric] = imregconfig('monomodal');

% optimizer parameters for multimodal
% optimizer.GrowthFactor = 1.05;
% optimizer.Epsilon = 1.5e-6;
% optimizer.InitialRadius = 0.0063;
% optimizer.MaximumIterations = 100;

% optimizer changes for monomodal
% optimizer.GradientMagnitudeTolerance = 1.0e-4;
% optimizer.MinimumStepLength = 1.0e-5;
% optimizer.MaximumStepLength = 0.0625;
optimizer.MaximumIterations = 300;
% optimizer.RelaxationFactor = 0.5;

geomtform = imregtform(moving,Rmoving, fixed,Rfixed, 'affine', optimizer, metric,'DisplayOptimization',false);
warpDose1 = imwarp(dose1,Rmoving,geomtform,'OutputView',Rfixed);

%% add photon target to particle cst, deleting structs from the cst without volume information (e.g. reference points)
% cst is the cell where the changes are done and that is use at the end for
% the DVH and the quality indicators

% deleting empty structs in reference cst
cst2(cellfun(@isempty,cst2(:,4)),:) = [];
      
% adding target structs with the help of the geometrical information
% from the image registration
for i = 1 : size(cst1)
        
    if strcmp(cst1{i,3} ,'TARGET')
    
        cube1 = zeros(ct1.cubeDim);
      
        structIndices = cst1{i,4}{1};
        cube1(structIndices) = 1;
        newPhotonCube = imwarp(cube1,Rmoving,geomtform,'OutputView',Rfixed);
        
        % ist das hier zulässig? vergrössern wir hier durch interpolation
        % bei imwarp nicht die volumina?
        newStructIndices = find(newPhotonCube > 0);
        
        cst2{end+1,1}    = size(cst2) + 1;
        cst2{end  ,2}    = [cst1{i,2} '_Photon'];
        cst2{end  ,3}    = 'TARGET';
        cst2{end  ,4}{1} = newStructIndices;
        cst2{end  ,5}    = cst1{i,5};
        cst2{end  ,6}    = cst1{i,6};
        
    end
end

%% EQD_2 calculations (calculates RBExD added dose, photon EQD_2, particle EQD_2, added EQD_2 and RBExD added divided by EQD_2 added [and invers ^-1])

result.totalDose = numOfFractions1 * warpDose1 + numOfFractions2 * dose2;

result.EQD2_1    = warpDose1 * numOfFractions1 .* ((warpDose1 + alphaBetaRatio) / (2 + alphaBetaRatio));
result.EQD2_2    = dose2 * numOfFractions2 .* ((dose2 + alphaBetaRatio) / (2 + alphaBetaRatio));
result.totalEQD2 = result.EQD2_1 + result.EQD2_2;

result.EQD2ratio       = result.totalEQD2 ./ result.totalDose;
result.EQD2ratioInvers = result.totalDose ./ result.totalEQD2;


%% DVH calculation
    
aimedDose = numOfFractions1 * prescribedDose1 + numOfFractions2 * prescribedDose2;
aimedEQD2 = numOfFractions1 * prescribedDose1 .* ((prescribedDose1 + alphaBetaRatio) / (2 + alphaBetaRatio)) + numOfFractions2 * prescribedDose2 .* ((prescribedDose2 + alphaBetaRatio) / (2 + alphaBetaRatio)); 
        
dvh_dose = matRad_calcDVH(cst2,result.totalDose ./ aimedDose);
dvh_EQD2 = matRad_calcDVH(cst2,result.totalEQD2 ./ aimedEQD2);
dvh_EQD2ratio = matRad_calcDVH(cst2,result.EQD2ratio);

figure;
matRad_showDVH(dvh_dose,cst2,pln2,1);
hold on
matRad_showDVH(dvh_EQD2,cst2,pln2,2);
title(['Added dose cubes divided by aimed dose, straight line RBExD added [aimed: ' num2str(aimedDose) 'Gy], dashed line EQD_2 added [aimed : ' num2str(aimedEQD2) 'Gy]']);
xlabel('relative Dose [%]');
legend(cst2{:,2});
hold off


figure;
matRad_showDVH(dvh_EQD2ratio,cst2,pln2,1);
title('Added EQD_2 divided by added RBExD');
xlabel('relative Dose [%]');

%% Qualitiy indicators
fieldNamesResult = fieldnames(result);

for resultGUInumber = 1 : numel(fieldNamesResult)
    qualityIndicators.(fieldNamesResult{resultGUInumber}) = matRad_calcQualityIndicators(cst2,pln2,result.(fieldNamesResult{resultGUInumber}));
end


%% image registration check (shows ct pictures)
if checkImageRegis == true
    movedCT = imwarp(moving,Rmoving,geomtform,'OutputView',Rfixed);
   
    figure,imshowpair(squeeze( fixed(:,50,:)),squeeze( moving(:,50,:)),'Scaling','joint');
    title('fixed vs moving');
    figure,imshowpair(squeeze( fixed(:,50,:)),squeeze( movedCT(:,50,:)),'scaling','joint');
    title('fixed vs moved');

    figure,imshowpair(fixed(:,:,50),moving(:,:,50),'Scaling','joint');
    title('fixed vs moving');
    figure,imshowpair( fixed(:,:,50),movedCT(:,:,50),'scaling','joint');
    title('fixed vs moved');
end