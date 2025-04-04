function VmcOptions = matRad_vmcOptions(pln,ct)

%% run options

% number of paralle MC simulations
if isfield(pln.propDoseCalc.vmcOptions,'numOfParMCSim')
    VmcOptions.run.numOfParMCSim    = pln.propDoseCalc.vmcOptions.numOfParMCSim;
else
    VmcOptions.run.numOfParMCSim    = 4;
end
if isunix && VmcOptions.run.numOfParMCSim > 1
    VmcOptions.run.numOfParMCSim = 1;
end

% number of histories per bixel
if isfield(pln.propDoseCalc.vmcOptions,'nCasePerBixel')
    VmcOptions.run.nCasePerBixel    = pln.propDoseCalc.vmcOptions.nCasePerBixel;
else
    VmcOptions.run.nCasePerBixel    = 5000;
end

% relative dose cutoff
VmcOptions.run.relDoseCutoff    = 10^(-3);

% version (Carleton, dkfz, etc.)
VmcOptions.run.version = pln.propDoseCalc.vmcOptions.version;

% set absolute calibration factor
% CALCULATION
% absolute_calibration_factor = 1/D(depth = 100,5mm) -> D(depth = 100,5mm) = 1Gy
% SETUP
% SAD = 1000mm, SCD = 500mm, bixelWidth = 5mm, IC = [240mm,240mm,240mm]
% fieldsize@IC = 105mm x 105mm, phantomsize = 81 x 81 x 81 = 243mm x 243mm x 243mm
% rel_Dose_cutoff = 10^(-3), ncase = 500000/bixel
switch pln.propDoseCalc.vmcOptions.version
    case 'Carleton'
        switch pln.propDoseCalc.vmcOptions.source
            case 'phsp'
                
                d_50mm = 9.351001892810018e-07;
                d_50mm_error = 7.668474434354598e-09;
                
                VmcOptions.run.absCalibrationFactorVmc      = 1./d_50mm;
                VmcOptions.run.absCalibrationFactorVmc_err  = d_50mm_error./(d_50mm.^2);
        end
    case 'dkfz'
        VmcOptions.run.absCalibrationFactorVmc  = 99.818252282632300;
end

%% source

VmcOptions.source.myName       = 'some_source';                                         % name of source
VmcOptions.source.monitorUnits = 1;
switch pln.propDoseCalc.vmcOptions.source
    case 'beamlet'
    VmcOptions.source.spectrum     = fullfile(runsPath,'spectra','var_6MV.spectrum');   % energy spectrum source (only used if no mono-Energy given)
    VmcOptions.source.charge       = 0;                                                 % charge (-1,0,1)
    VmcOptions.source.type         = 'beamlet';
    
    case 'phsp'
    VmcOptions.source.particleType  = 2;
    VmcOptions.source.type          = 'phsp';
end

%% transport parameters

VmcOptions.McParameter.automatic_parameter  = 'yes';                       % if yes, automatic transport parameters are used
VmcOptions.McParameter.spin                 = 0;                           % 0: spin effects ignored; 1: simplistic; 2: full treatment

%% MC control

VmcOptions.McControl.ncase  = VmcOptions.run.nCasePerBixel;                % number of histories
VmcOptions.McControl.nbatch = 10;                                          % number of batches

%% variance reduction

VmcOptions.varianceReduction.repeatHistory      = 0.041;
VmcOptions.varianceReduction.splitPhotons       = 1;   
VmcOptions.varianceReduction.photonSplitFactor  = -80;

%% quasi random numbers

VmcOptions.quasi.base      = 2;                                                 
VmcOptions.quasi.dimension = 60;                                             
VmcOptions.quasi.skip      = 1;

%% geometry
switch pln.propDoseCalc.vmcOptions.version
    case 'Carleton'
        VmcOptions.geometry.XyzGeometry.methodOfInput = 'MMC-PHANTOM';  % input method ('CT-PHANTOM', 'individual', 'groups')
    case 'dkfz'
        VmcOptions.geometry.XyzGeometry.methodOfInput = 'CT-PHANTOM';   % input method ('CT-PHANTOM', 'individual', 'groups')
end
VmcOptions.geometry.dimensions          = ct.cubeDim;
VmcOptions.geometry.XyzGeometry.Ct      = 'CT';                         % name of geometry

%% scoring manager
VmcOptions.scoringOptions.startInGeometry               = 'CT';            % geometry in which partciles start their transport
VmcOptions.scoringOptions.doseOptions.scoreInGeometries = 'CT';            % geometry in which dose is recorded
VmcOptions.scoringOptions.doseOptions.scoreDoseToWater  = 'yes';           % if yes output is dose to water
VmcOptions.scoringOptions.outputOptions.name            = 'CT';            % geometry for which dose output is created (geometry has to be scored)
VmcOptions.scoringOptions.outputOptions.dumpDose        = pln.propDoseCalc.vmcOptions.dumpDose;               % output format (1: format=float, Dose + deltaDose; 2: format=short int, Dose)



end