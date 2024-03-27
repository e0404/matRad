function dij = matRad_calcPhotonDoseOmpMC(ct,stf,pln,cst,calcDoseDirect)
% matRad ompMC monte carlo photon dose calculation wrapper
%
% call
%   dij = matRad_calcPhotonDoseMC(ct,stf,pln,cst)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   visBool:                    binary switch to enable visualization
% output
%   dij:                        matRad dij struct
%
% References
%   -
%
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instance of MatRad_Config class
matRad_cfg =  MatRad_Config.instance();

% handle inputs
if nargin < 6
    calcDoseDirect = false;
end

if calcDoseDirect
    matRad_cfg.dispWarning('ompMC does not support forward dose calculation natively! Will compute full dij and then apply fluence vector.');
end

% load default parameters for matRad_OmpConfig class in case they haven't been set yet
if ~isa(pln.propMC,'matRad_OmpConfig')
    pln = matRad_cfg.getDefaultClass(pln,'propMC','matRad_OmpConfig');
end

% load default parameters for doseCalc in case they haven't been set yet
pln = matRad_cfg.getDefaultProperties(pln,{'propDoseCalc'});

% compile ompMC interface
if ~matRad_checkMexFileExists('omc_matrad') %exist('matRad_ompInterface','file') ~= 3
    matRad_cfg.dispWarning('Compiled mex interface not found. Trying to compile the ompMC interface on the fly!');
    try
        matRad_OmpConfig.compileOmpMCInterface();
    catch MException
        matRad_cfg.dispError('Could not find/generate mex interface for MC dose calculation.\nCause of error:\n%s\n Please compile it yourself (preferably with OpenMP support).',MException.message);
    end
end

%% Initialize dose grid and dij

% load calcDoseInit as usual
matRad_calcDoseInit;

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfBixels,1);
dij.rayNum   = NaN*ones(dij.totalNumOfBixels,1);
dij.beamNum  = NaN*ones(dij.totalNumOfBixels,1);

dij.numHistoriesPerBeamlet = pln.propMC.numHistories;

%% Setup OmpMC options / parameters
ompMCoptions = pln.propMC.getOmpMCoptions(machine);

% conversion from HU to densities & materials
[~,cubeMatIx,cubeRho] = pln.propMC.materialConversion(dij.ctGrid,dij.doseGrid,ct);

% Get ompMC geometry
ompMCgeo = pln.propMC.getOmpMCgeometry(dij.doseGrid);


%% debug visualization
if pln.propMC.visBool
    
    figure
    hold on
    
    axis equal
    
    % ct box
    ctCorner1 = [ompMCgeo.xBounds(1) ompMCgeo.yBounds(1) ompMCgeo.zBounds(1)];
    ctCorner2 = [ompMCgeo.xBounds(end) ompMCgeo.yBounds(end) ompMCgeo.zBounds(end)];
    plot3([ctCorner1(1) ctCorner2(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner1(3) ctCorner1(3)],'k' )
    plot3([ctCorner1(1) ctCorner2(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner1(3) ctCorner1(3)],'k' )
    plot3([ctCorner1(1) ctCorner1(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner1(3) ctCorner1(3)],'k' )
    plot3([ctCorner2(1) ctCorner2(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner1(3) ctCorner1(3)],'k' )
    plot3([ctCorner1(1) ctCorner2(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner2(3) ctCorner2(3)],'k' )
    plot3([ctCorner1(1) ctCorner2(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner2(3) ctCorner2(3)],'k' )
    plot3([ctCorner1(1) ctCorner1(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner2(3) ctCorner2(3)],'k' )
    plot3([ctCorner2(1) ctCorner2(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner2(3) ctCorner2(3)],'k' )
    plot3([ctCorner1(1) ctCorner1(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner1(3) ctCorner2(3)],'k' )
    plot3([ctCorner2(1) ctCorner2(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner1(3) ctCorner2(3)],'k' )
    plot3([ctCorner1(1) ctCorner1(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner1(3) ctCorner2(3)],'k' )
    plot3([ctCorner2(1) ctCorner2(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner1(3) ctCorner2(3)],'k' )
    
    xlabel('x [cm]')
    ylabel('y [cm]')
    zlabel('z [cm]')
    
    rotate3d on
    
end


% Now we have to calibrate to the the beamlet width.
pln.propMC.absCalibrationFactor = pln.propMC.absCalibrationFactor * (pln.propStf.bixelWidth/50)^2;


%% Create beamlet source
scenCount = 0;

for scenarioIx = 1:pln.multScen.totNumScen
    
    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter + pln.multScen.isoShift(scenarioIx,:);
    end

    ctScen = pln.multScen.linearMask(scenarioIx,1);
    shiftScen = pln.multScen.linearMask(scenarioIx,2);
    rangeShiftScen = pln.multScen.linearMask(scenarioIx,3);
    
    if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
        scenCount = scenCount + 1;
        
        % load ompMC source
        ompMCsource = pln.propMC.getOmpMCsource(dij.numOfBeams,dij.totalNumOfBixels,stf);
        
        % Book keeping for dij
        counter = 0;
        for i = 1:dij.numOfBeams
            for j = 1:stf(i).numOfRays
                counter = counter + 1;
                dij.beamNum(counter)  = i;
                dij.rayNum(counter)   = j;
                dij.bixelNum(counter) = j;
            end
        end
            
        %% Call the OmpMC interface
        if pln.multScen.totNumScen > 1
            matRad_cfg.dispInfo('matRad: OmpMC photon dose calculation... \n');
        else
            matRad_cfg.dispInfo('matRad: OmpMC photon dose calculation for scenario %d of %d... \n',scenCount,pln.multScen.totNumScen);
        end
        
        ompMCgeo.isoCenter = [stf(:).isoCenter];
        
        %Call the Monte Carlo simulation and catch  possible mex
        %interface issues       
        try
            %If we ask for variance, a field in the dij will be filled
            if pln.propMC.outputVariance
                [dij.physicalDose{ctScen,shiftScen,rangeShiftScen},dij.physicalDose_MCvar{ctScen,shiftScen,rangeShiftScen}] = omc_matrad(cubeRho{ctScen},cubeMatIx{ctScen},ompMCgeo,ompMCsource,ompMCoptions);
            else
                [dij.physicalDose{ctScen,shiftScen,rangeShiftScen}] = omc_matrad(cubeRho{ctScen},cubeMatIx{ctScen},ompMCgeo,ompMCsource,ompMCoptions);
            end
        catch ME
            errorString = [ME.message '\nThis error was thrown by the MEX-interface of ompMC.\nMex interfaces can raise compatability issues which may be resolved by compiling them by hand directly on your particular system.'];
            matRad_cfg.dispError(errorString);
        end
        
        dij.physicalDose{ctScen,shiftScen,rangeShiftScen} = dij.physicalDose{ctScen,shiftScen,rangeShiftScen} * pln.propMC.absCalibrationFactor;
        if isfield(dij,'physicalDose_MCvar')
            dij.physicalDose_MCvar{s} = dij.physicalDose_MCvar{ctScen,shiftScen,rangeShiftScen} * pln.propMC.absCalibrationFactor^2;
        end
    
        matRad_cfg.dispInfo('matRad: MC photon dose calculation done!\n');
        
        try
            % wait 0.1s for closing all waitbars
            allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar');
            delete(allWaitBarFigures);
            pause(0.1);
        catch
        end
        
        % manipulate isocenter back
        for k = 1:length(stf)
            stf(k).isoCenter = stf(k).isoCenter - pln.multScen.isoShift(scenarioIx,:);
        end
    end
end


end
