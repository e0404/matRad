% parameters that can be changed to fit to the RT plans, also alpha/beta
% for the EQD_2 calculation can be changed
photonFileName = 'RTPlan_9W5A49NQ_photons.mat';         % name of the photon file      

particleFileName = 'RTPlan_9W5A49NQ_protons.mat';        % name of the particle file
 
checkImageRegis = false;                                % true or false, shows pictures of the particle and photon CT to check if the image registration went wrong (true shows pictures)

alphaBetaRatio = 2;                                     % alpha/beta value

numOfFractionsPhoton = 25;                              % number of photon fractions
numOfFractionsProtons= 5;                               % number of proton fractions
numOfFractionsCarbon = 6;                               % number of carbon fractions

aimedDosePhotonPF = 2;                                  % aimed in the target volume per fraction
aimedDoseProtonsPF = 2;                                 % aimed in the target volume per fraction
aimedDoseCarbonPF = 3;                                  % aimed in the target volume per fraction

%% image registration and recalculation of the photon dose cube to fit in the particle dose cube

load(photonFileName);
ct_photon = ct;
cst_photon = cst;

if isfield(resultGUI, {'addedDose'})
    resultGUI_Photon.physicalDose = resultGUI.addedDose;
elseif isfield(resultGUI, {'physicalDose'})
    resultGUI_Photon.physicalDose = resultGUI.physicalDose;
end

load(particleFileName);
particle = pln.radiationMode;
     
ct_particle = ct;
cst_particle = cst;
resultGUI_particle = resultGUI;

clear resultGUI

%fixed cube is the cube from the particle ct, photon ct is the moving cube
%Rfixed/Rmoving are helping the optimizer
fixed = ct_particle.cubeHU{1};
Rfixed = imref3d(ct_particle.cubeDim,ct_particle.resolution.y,ct_particle.resolution.x,ct_particle.resolution.z);
moving = ct_photon.cubeHU{1};
Rmoving = imref3d(ct_photon.cubeDim,ct_photon.resolution.y,ct_photon.resolution.x,ct_photon.resolution.z);

%choose of the optimizer for the image registration, for this ct
%registration monomodal is highly recommended
%below are possible changes for the optimization, values are the default
%values; only change is done in the maximum iterations (default 100)

%[optimizer, metric] = imregconfig('multimodal');
[optimizer, metric] = imregconfig('monomodal');

%optimizer parameters for multimodal
% optimizer.GrowthFactor = 1.05;
% optimizer.Epsilon = 1.5e-6;
% optimizer.InitialRadius = 0.0063;
% optimizer.MaximumIterations = 100;

%optimizer changes for monomodal
% optimizer.GradientMagnitudeTolerance = 1.0e-4;
% optimizer.MinimumStepLength = 1.0e-5;
% optimizer.MaximumStepLength = 0.0625;
optimizer.MaximumIterations = 300;
% optimizer.RelaxationFactor = 0.5;


geomtform = imregtform(moving,Rmoving, fixed,Rfixed, 'affine', optimizer, metric,'DisplayOptimization',false);
resultGUI.movedPhotonDose = imwarp(resultGUI_Photon.physicalDose,Rmoving,geomtform,'OutputView',Rfixed);

%% add photon target to particle cst, deleting structs from the cst without volume information (e.g. reference points)
% cst is the cell were the changes are done and that is use at the end for
% the DVH and the quality indicators

cst_photonSize = size(cst_photon);
cst_ParticleSize = size(cst_particle);
cst_photonTargetNumber = 0;
cst_shift = 0;

%deleting structs
for cstNumberParticle = 1 : cst_ParticleSize(1)
    if isempty(cst_particle{cstNumberParticle,4}{1})

        cst(cstNumberParticle - cst_shift,:) = [];
        cst_shift = cst_shift + 1;
       
    end
end
        
cst_Size = size(cst);


%adding photon target structs with the help of the geometrical formation
%from the image registration
for cstNumberPhoton = 1 : cst_photonSize(1)
    
    photonCube = zeros(ct_photon.cubeDim);
    
    if strcmp(cst_photon{cstNumberPhoton,3} ,'TARGET')
        
        cst_photonTargetNumber = cst_photonTargetNumber +1;
        
        structIndices = cst_photon{cstNumberPhoton,4}{1};
        photonCube(structIndices) = 1;
        newPhotonCube = imwarp(photonCube,Rmoving,geomtform,'OutputView',Rfixed);
        newStructIndices = find(newPhotonCube > 0);
        
        cst{cst_Size(1)+ cst_photonTargetNumber,1} = cst_Size(1)+ cst_photonTargetNumber;
        cst{cst_Size(1)+ cst_photonTargetNumber,2} = [cst_photon{cstNumberPhoton,2} '_Photon'];
        cst{cst_Size(1)+ cst_photonTargetNumber,3} = 'TARGET';
        cst{cst_Size(1)+ cst_photonTargetNumber,4}{1} = newStructIndices;
        cst{cst_Size(1)+ cst_photonTargetNumber,5} = cst_photon{cstNumberPhoton,5};
        cst{cst_Size(1)+ cst_photonTargetNumber,6} = cst_photon{cstNumberPhoton,6};
        
    end
end

%% EQD_2 calculations (calculates RBExD added dose, photon EQD_2, particle EQD_2, added EQD_2 and RBExD added divided by EQD_2 added [and invers ^-1])

fieldNamesresultGUIparticle = fieldnames(resultGUI_particle);   
fieldNumbresultGUIparticle = numel(fieldNamesresultGUIparticle);         
particleRBExDName = string(fieldNamesresultGUIparticle(fieldNumbresultGUIparticle - 1));

particleRBExD = resultGUI_particle.(particleRBExDName);
photonPhysicalDose = resultGUI.movedPhotonDose;

if max(max(max(particleRBExD))) < 1
    if strcmp(particle,'carbon')
        particleRBExD = numOfFractionsParticle * particleRBExD;
    elseif strcmp(particle,'protons')
        particleRBExD = numOfFractionsParticle * particleRBExD;
    end
end

resultGUI.particleRBExD = particleRBExD;

if strcmp(particle,'carbon')
    numOfFractionsParticle = numOfFractionsCarbon;
    aimedDoseParticlePF = aimedDoseCarbonPF;
elseif strcmp(particle,'protons')
    numOfFractionsParticle = numOfFractionsProtons;
    aimedDoseParticlePF = aimedDoseProtonsPF;
end

resultGUI.addedRBExD = particleRBExD * numOfFractionsParticle + photonPhysicalDose * numOfFractionsPhoton;

resultGUI.photonEQD_2 = photonPhysicalDose * numOfFractionsPhoton .* ((photonPhysicalDose + alphaBetaRatio) / (2 + alphaBetaRatio));
resultGUI.particleEQD_2 = particleRBExD * numOfFractionsParticle .* ((particleRBExD + alphaBetaRatio) / (2 + alphaBetaRatio));
resultGUI.addedEQD_2 = resultGUI.photonEQD_2 + resultGUI.particleEQD_2;

resultGUI.addedEQD_2_addedRBExD_ratio = resultGUI.addedEQD_2 ./ resultGUI.addedRBExD;
resultGUI.addedEQD_2_addedRBExD_ratio_invers = (resultGUI.addedEQD_2 ./ resultGUI.addedRBExD).^(-1);

%% DVH calculation
    
aimedDoseAddedRBExD = numOfFractionsParticle * aimedDoseParticlePF + numOfFractionsPhoton * aimedDosePhotonPF;
aimedDoseAddedEQD_2 =  aimedDoseParticlePF * numOfFractionsParticle .* ((aimedDoseParticlePF + alphaBetaRatio) / (2 + alphaBetaRatio)) + aimedDosePhotonPF * numOfFractionsPhoton .* ((aimedDosePhotonPF + alphaBetaRatio) / (2 + alphaBetaRatio)); 
        
dvh_RBExD = matRad_calcDVH(cst,resultGUI.addedRBExD ./ aimedDoseAddedRBExD);
dvh_EQD_2 = matRad_calcDVH(cst,resultGUI.addedEQD_2 ./ aimedDoseAddedEQD_2);
dvh_ratioRBExD_EQD_2 = matRad_calcDVH(cst,resultGUI.addedEQD_2_addedRBExD_ratio);

figure;
matRad_showDVH(dvh_RBExD,cst,pln,1);
hold on
matRad_showDVH(dvh_EQD_2,cst,pln,2);
title(['Added dose cubes divided by aimed dose, straight line RBExD added [aimed: ' num2str(aimedDoseAddedRBExD) 'Gy], dashed line EQD_2 added [aimed : ' num2str(aimedDoseAddedEQD_2) 'Gy]']);
xlabel('relative Dose [%]');
legend(cst{:,2});
hold off


figure;
matRad_showDVH(dvh_ratioRBExD_EQD_2,cst,pln,1);
title('Added EQD_2 divided by added RBExD');
xlabel('relative Dose [%]');

%% Qualitiy indicators

fieldNamesResultGUI = fieldnames(resultGUI);
resultGUIsize = numel(fieldNamesResultGUI);


for resultGUInumber = 1 : resultGUIsize
    qualityIndicators.(fieldNamesResultGUI{resultGUInumber}) = matRad_calcQualityIndicators(cst,pln,resultGUI.(fieldNamesResultGUI{resultGUInumber}));
end


%% image registration check (shows ct pictures)

if checkImageRegis == true
    movedCT = imwarp(moving,Rmoving,geomtform,'OutputView',Rfixed);

    resultGUI.notmovedPhoton = resultGUI_Photon.physicalDose;
    
    figure,imshowpair(squeeze( fixed(:,50,:)),squeeze( moving(:,50,:)),'Scaling','joint');
    title('fixed vs moving');
    figure,imshowpair(squeeze( fixed(:,50,:)),squeeze( movedCT(:,50,:)),'scaling','joint');
    title('fixed vs moved');

    figure,imshowpair(fixed(:,:,50),moving(:,:,50),'Scaling','joint');
    title('fixed vs moving');
    figure,imshowpair( fixed(:,:,50),movedCT(:,:,50),'scaling','joint');
    title('fixed vs moved');
end