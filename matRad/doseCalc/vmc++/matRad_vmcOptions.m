function vmcOpts = matRad_vmcOptions(pln, ct)

%% run options

% number of parallel MC simulations
if isfield(pln.propDoseCalc.vmcOptions, 'numOfParMCSim')
    vmcOpts.run.numOfParMCSim    = pln.propDoseCalc.vmcOptions.numOfParMCSim;
else
    vmcOpts.run.numOfParMCSim    = 4;
end
if isunix && vmcOpts.run.numOfParMCSim > 1
    vmcOpts.run.numOfParMCSim = 1;
end

% number of histories per bixel
if isfield(pln.propDoseCalc.vmcOptions, 'nCasePerBixel')
    vmcOpts.run.nCasePerBixel    = pln.propDoseCalc.vmcOptions.nCasePerBixel;
else
    vmcOpts.run.nCasePerBixel    = 5000;
end

% relative dose cutoff
vmcOpts.run.relDoseCutoff    = 10^(-3);

% version (Carleton, dkfz, etc.)
vmcOpts.run.version = pln.propDoseCalc.vmcOptions.version;

% set absolute calibration factor
% CALCULATION
% absolute_calibration_factor = 1/D(depth = 100,5mm) -> D(depth = 100,5mm) = 1Gy
% SETUP
% SAD = 1000mm, SCD = 500mm, bixelWidth = 5mm, IC = [240mm,240mm,240mm]
% fieldsize@IC = 105mm x 105mm, phantomsize = 81 x 81 x 81 = 243mm x 243mm x 243mm
% rel_Dose_cutoff = 10^(-3), ncase = 500000/bixel
vmcOpts.run.absCalibrationFactorVmc_err = 0;
switch pln.propDoseCalc.vmcOptions.version
    case 'Carleton'
        switch pln.propDoseCalc.vmcOptions.source
            case 'phsp'

                d_50mm = 9.351001892810018e-07;
                d_50mm_error = 7.668474434354598e-09;

                vmcOpts.run.absCalibrationFactorVmc      = 1 ./ d_50mm;
                vmcOpts.run.absCalibrationFactorVmc_err  = d_50mm_error ./ (d_50mm.^2);
        end
    case 'dkfz'
        vmcOpts.run.absCalibrationFactorVmc  = 99.818252282632300;
end

%% source

% locate the runs folder of the VMC++ environment for the spectrum file
matRad_cfg = MatRad_Config.instance();
vmcPath = fullfile(matRad_cfg.matRadSrcRoot, 'doseCalc', 'vmc++');
switch pln.propDoseCalc.vmcOptions.version
    case 'Carleton'
        runsPath = fullfile(vmcPath, 'run');
    case 'dkfz'
        runsPath = fullfile(vmcPath, 'runs');
end

vmcOpts.source.myName       = 'some_source';                                         % name of source
vmcOpts.source.monitorUnits = 1;
switch pln.propDoseCalc.vmcOptions.source
    case 'beamlet'
        % energy spectrum source (only used if no mono-Energy given)
        vmcOpts.source.spectrum     = fullfile(runsPath, 'spectra', 'var_6MV.spectrum');
        vmcOpts.source.charge       = 0;                                                 % charge (-1,0,1)
        vmcOpts.source.type         = 'beamlet';

    case 'phsp'
        vmcOpts.source.particleType  = 2;
        vmcOpts.source.type          = 'phsp';
end

%% transport parameters

vmcOpts.McParameter.automatic_parameter  = 'yes';                       % if yes, automatic transport parameters are used
vmcOpts.McParameter.spin                 = 0;                           % 0: spin effects ignored; 1: simplistic; 2: full treatment

%% MC control

vmcOpts.McControl.ncase  = vmcOpts.run.nCasePerBixel;                % number of histories
vmcOpts.McControl.nbatch = 10;                                          % number of batches

%% variance reduction

vmcOpts.varianceReduction.repeatHistory      = 0.041;
vmcOpts.varianceReduction.splitPhotons       = 1;
vmcOpts.varianceReduction.photonSplitFactor  = -80;

%% quasi random numbers

vmcOpts.quasi.base      = 2;
vmcOpts.quasi.dimension = 60;
vmcOpts.quasi.skip      = 1;

%% geometry
switch pln.propDoseCalc.vmcOptions.version
    case 'Carleton'
        vmcOpts.geometry.XyzGeometry.methodOfInput = 'MMC-PHANTOM';  % input method ('CT-PHANTOM', 'individual', 'groups')
    case 'dkfz'
        vmcOpts.geometry.XyzGeometry.methodOfInput = 'CT-PHANTOM';   % input method ('CT-PHANTOM', 'individual', 'groups')
end
vmcOpts.geometry.dimensions          = ct.cubeDim;
vmcOpts.geometry.XyzGeometry.Ct      = 'CT';                         % name of geometry

%% scoring manager
vmcOpts.scoringOptions.startInGeometry               = 'CT';            % geometry in which particles start their transport
vmcOpts.scoringOptions.doseOptions.scoreInGeometries = 'CT';            % geometry in which dose is recorded
vmcOpts.scoringOptions.doseOptions.scoreDoseToWater  = 'yes';           % if yes output is dose to water
vmcOpts.scoringOptions.outputOptions.name            = 'CT';            % geometry for which dose output is created (geometry has to be scored)
% output format (1: format=float, Dose + deltaDose; 2: format=short int, Dose)
vmcOpts.scoringOptions.outputOptions.dumpDose        = pln.propDoseCalc.vmcOptions.dumpDose;

end
