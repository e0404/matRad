% this script creates a testing ct,cst,stf,pln for external radiation
% therapy, that can be red in and used in testing

%% create ct
ct = struct();
ct.cubeDim = [20,10,10];
ct.resolution.x = 10;
ct.resolution.y = 10;
ct.resolution.z = 10;
ct.numOfCtScen = 1;
ct.cubeHU{1} = ones(ct.cubeDim) * -1000;

ct.cubeHU{1}(2:19,2:9,2:9) = 0;
VolHelper = false(ct.cubeDim);
VolHelper(2:19,2:9,2:9) = true;
ixBody = find(VolHelper);
VolHelper = false(ct.cubeDim);
VolHelper(10:11,5:6,5:6) = true;
ixTarget = find(VolHelper);

for i = 1:2
    cst{i,1}                = i-1;
    cst{i,5}.TissueClass    = 1;
    cst{i,5}.alphaX         =  0.1000;
    cst{i,5}.betaX          = 0.0500;
    cst{i,5}.Priority       = i;
    cst{i,5}.Visible        = 1;
    cst{i,5}.visibleColor   = [0 0 0];

end
cst{1,2} = 'Target'; 
cst{1,3} = 'TARGET'; 
cst{1,4}{1} = ixTarget;
cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(800,45));
cst{2,2} = 'Body'; 
cst{2,3} = 'OAR'; 
cst{2,4}{1} = ixBody;
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(400,0));

clear VolHelper ixBody ixTarget i
%% create pln, stf
radModes = ["protons","helium","carbon","VHEE"];
for radMode = radModes
    %radMode = 'carbon'; %protons,helium,carbon;
    
    pln.radiationMode = char(radMode);            
    pln.machine               = 'Generic';                                         
    pln.numOfFractions        = 30;
    pln.propStf.gantryAngles  = [0,180];
    pln.propStf.couchAngles   = zeros(size(pln.propStf.gantryAngles));
    pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
    pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
    
    pln.propStf.longitudinalSpotSpacing = 8;
    pln.propStf.bixelWidth = 10;
    pln.propDoseCalc.doseGrid.resolution = struct('x',10,'y',10,'z',10); %[mm]
    
    %pln.bioModel = matRad_bioModel(pln.radiationMode,'none');
    
    %% Generate Beam Geometry STF
    pln.propStf.addMargin    = false; %to make smaller stf, les bixel
    stf = matRad_generateStf(ct,cst,pln);
    ct = matRad_calcWaterEqD(ct, pln);
    %% Dose Calculation
    dij = matRad_calcDoseInfluence(ct,cst,stf,pln);
    resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij);
    
    %% save basic data
    save([char(radMode) '_testData.mat'],'ct','cst','pln','stf','dij','resultGUI','-v7');
end