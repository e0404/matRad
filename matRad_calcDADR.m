function  resultGUI = matRad_calcDADR(dij,resultGUI,I)
%
% function to compute the dose-averaged dose-rate for each voxel of the dose grid
%
% input
%   dij:        matRad dij struct, dij dij.physicalDose
%   resultGUI:  matRad resultGUI struct, resultGUI.w
%
% Output
% dose_ratei:           matRad dose_ratei array: average dose-rate for each voxel i
%
% References
%   [1] https://pubmed.ncbi.nlm.nih.gov/33068294/
%                                          
%                                          
matRad_cfg = MatRad_Config.instance();
w = resultGUI.w;

if nargin < 3 && ~isfield(dij,'fixedCurrent')
    defaultIntensity = 300; %nA         3e9/1e6; %default intensity noramlized for 1e6 particles
    I = ones(length(w), 1)*defaultIntensity; 
    matRad_cfg.dispInfo('Using default Intensity of %g nA for each beamlet!\n',defaultIntensity);
elseif nargin < 3 && isfield(dij,'fixedCurrent')
    I = ones(length(w),1)*dij.fixedCurrent;
    matRad_cfg.dispInfo('Using fixed current  of %g nA for each beamlet!\n',dij.fixedCurrent);
end

%We assume I to be in nA
nA2particles = 1e-9 ./ 1.602176634e-19 ./ 1e6; %
I = nA2particles * I;

scenNum = 1;

% get bixel - beam correspondence  
for i = 1:dij.numOfBeams
    beamInfo(i).suffix = ['_beam', num2str(i)];
    beamInfo(i).logIx  = (dij.beamNum == i);      
end
beamInfo(dij.numOfBeams+1).suffix = '';
beamInfo(dij.numOfBeams+1).logIx  = true(size(w));
%

matRad_cfg.dispWarning('Dose rate computation uses the physical Dose (not RBE-weighted dose!)');

% compute physical dose for all beams individually and together
index=0;
for i = 1:length(beamInfo)
    % Computation of the DADR
    doseRateBeam = reshape(full(dij.physicalDose{scenNum}.^2 * (resultGUI.w .* I .* beamInfo(i).logIx)),dij.doseGrid.dimensions);
    beamDose =  reshape(full(dij.physicalDose{scenNum} * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions);
    doseRateBeam(beamDose > 0) = doseRateBeam(beamDose > 0) ./ beamDose(beamDose > 0);

    % the computed dose rate is then added to the resultGUI struct
    resultGUI.(['DADR', beamInfo(i).suffix]) = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
                                              doseRateBeam, ...
                                             dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);
                                        
end

doseRate = reshape(full(dij.physicalDose{scenNum}.^2 * (resultGUI.w .* I)),dij.doseGrid.dimensions);
dose =  reshape(full(dij.physicalDose{scenNum} * resultGUI.w),dij.doseGrid.dimensions);
doseRate(dose > 0) = doseRate(dose > 0) ./ dose(dose > 0);

% the computed dose rate is then added to the resultGUI struct
resultGUI.DADR = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
    doseRate, ...
    dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);

resultGUI.I = I;

end


           