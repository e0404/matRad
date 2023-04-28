classdef matRad_OmpConfig < handle
    % matRad_TopasConfig class definition
    %
    %
    % References
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2022 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        % This parameter can be overwritten through MatRad_Config default parameters
        numHistories = 1e6; %Number of histories to compute

        visBool = false; % disable visualiazation by default

        engine = 'ompMC'; %parameter for continuity
        
        %ompMC for matRad returns dose/history * nHistories. 
        % This factor calibrates to 1 Gy in a %(5x5)cm^2 open field (1 bixel) at
        % 5cm depth for SSD = 900 which corresponds to the calibration for the
        % analytical base data.
        absCalibrationFactor = 3.49056 * 1e12; %Approximate!
        
        scale = 10; % to convert to cm
        
        outputVariance;
        useCornersSCD = true; %false -> use ISO corners
    end

    properties% (SetAccess = private)
        omcFolder;

    end

    methods
        function obj = matRad_OmpConfig()
            matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class

            % Default execution paths are set here
            obj.omcFolder = [matRad_cfg.matRadRoot filesep 'ompMC'];
 
            % Set default histories from MatRad_Config
            if isfield(matRad_cfg.propMC,'ompMC_defaultHistories')
                obj.numHistories = matRad_cfg.propMC.ompMC_defaultHistories;
            end
            
            outputVariance = matRad_cfg.propMC.ompMC_defaultOutputVariance;

        end


        function ompMCoptions = getOmpMCoptions(obj,machine,pln)
            matRad_cfg = MatRad_Config.instance(); % Instance of matRad configuration class
            
            %display options
            ompMCoptions.verbose = matRad_cfg.logLevel - 1;
            
            % start MC control
            ompMCoptions.nHistories = obj.numHistories;
            ompMCoptions.nSplit = 20;
            ompMCoptions.nBatches = 10;
            ompMCoptions.randomSeeds = [97 33];
            
            %start source definition
            % TODO make spectrum file a variable that can be loaded in and changed
            ompMCoptions.spectrumFile       = [obj.omcFolder filesep 'spectra' filesep 'mohan6.spectrum'];
            ompMCoptions.monoEnergy         = 0.1;
            ompMCoptions.charge             = 0;
            ompMCoptions.sourceGeometry     = 'gaussian';
            ompMCoptions.sourceGaussianWidth = 0.1 * obj.setGaussianPenumbra(machine);
            
            % start MC transport
            ompMCoptions.dataFolder   = [obj.omcFolder filesep 'data' filesep];
            ompMCoptions.pegsFile     = [obj.omcFolder filesep 'pegs4' filesep '700icru.pegs4dat'];
            ompMCoptions.pgs4formFile = [obj.omcFolder filesep 'pegs4' filesep 'pgs4form.dat'];
            
            ompMCoptions.global_ecut = 0.7;
            ompMCoptions.global_pcut = 0.010;
            
            % Relative Threshold for dose
            ompMCoptions.relDoseThreshold = 1 - pln.propDoseCalc.lateralCutOff;
            
            % Output folders
            ompMCoptions.outputFolder = [obj.omcFolder filesep 'output' filesep];
            
        end
        
        function ompMCgeo = getOmpMCgeometry(obj,doseGrid)
            ompMCgeo.xBounds = (doseGrid.resolution.y * (0.5 + [0:doseGrid.dimensions(1)])) ./ obj.scale;
            ompMCgeo.yBounds = (doseGrid.resolution.x * (0.5 + [0:doseGrid.dimensions(2)])) ./ obj.scale;
            ompMCgeo.zBounds = (doseGrid.resolution.z * (0.5 + [0:doseGrid.dimensions(3)])) ./ obj.scale;
        
            % Create Material Density Cube
            ompMCgeo.material    = obj.setupMaterials();
        end
                
        function ompMCsource = getOmpMCsource(obj,numOfBeams,totalNumOfBixels,stf)          
                numOfBixels = [stf(:).numOfRays];
                beamSource = zeros(numOfBeams, 3);
                
                bixelCorner = zeros(totalNumOfBixels,3);
                bixelSide1 = zeros(totalNumOfBixels,3);
                bixelSide2 = zeros(totalNumOfBixels,3);
                
                beamNum = zeros(1,prod(sum(numOfBixels),numOfBeams));
                counter = 0;
                
                for i = 1:numOfBeams % loop over all beams
                    
                    % define beam source in physical coordinate system in cm
                    beamSource(i,:) = (stf(i).sourcePoint + stf(i).isoCenter)/10;
                    
                    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
                        
                        counter = counter + 1;
                        
                        beamNum(counter)  = i;
                        
                        if obj.useCornersSCD
                            beamletCorners = stf(i).ray(j).rayCorners_SCD;
                        else
                            beamletCorners = stf(i).ray(j).beamletCornersAtIso;
                        end
                        
                        % get bixel corner and delimiting vectors.
                        % a) change coordinate system (Isocenter cs-> physical cs) and units mm -> cm
                        currCorner = (beamletCorners(1,:) + stf(i).isoCenter) ./ obj.scale;
                        bixelCorner(counter,:) = currCorner;
                        bixelSide1(counter,:) = (beamletCorners(2,:) + stf(i).isoCenter) ./ obj.scale - currCorner;
                        bixelSide2(counter,:) = (beamletCorners(4,:) + stf(i).isoCenter) ./ obj.scale - currCorner;
                        
                        if obj.visBool
                            for k = 1:4
                                currCornerVis = (beamletCorners(k,:) + stf(i).isoCenter)/10;
                                % rays connecting source and ray corner
                                plot3([beamSource(i,1) currCornerVis(1)],[beamSource(i,2) currCornerVis(2)],[beamSource(i,3) currCornerVis(3)],'b')
                                % connection between corners
                                lRayCorner = (beamletCorners(mod(k,4) + 1,:) + stf(i).isoCenter)/10;
                                plot3([lRayCorner(1) currCornerVis(1)],[lRayCorner(2) currCornerVis(2)],[lRayCorner(3) currCornerVis(3)],'r')
                            end
                        end
                        
                    end
                    
                end
                
                ompMCsource.nBeams = numOfBeams;
                ompMCsource.iBeam = beamNum(:);
                
                % Switch x and y directions to match ompMC cs.
                ompMCsource.xSource = beamSource(:,2);
                ompMCsource.ySource = beamSource(:,1);
                ompMCsource.zSource = beamSource(:,3);
                
                ompMCsource.nBixels = sum(numOfBixels(:));
                ompMCsource.xCorner = bixelCorner(:,2);
                ompMCsource.yCorner = bixelCorner(:,1);
                ompMCsource.zCorner = bixelCorner(:,3);
                
                ompMCsource.xSide1 = bixelSide1(:,2);
                ompMCsource.ySide1 = bixelSide1(:,1);
                ompMCsource.zSide1 = bixelSide1(:,3);
                
                ompMCsource.xSide2 = bixelSide2(:,2);
                ompMCsource.ySide2 = bixelSide2(:,1);
                ompMCsource.zSide2 = bixelSide2(:,3);
                
                if obj.visBool
                    plot3(ompMCsource.ySource,ompMCsource.xSource,ompMCsource.zSource,'rx')
                end
                
        end
        
        function sigmaGauss = setGaussianPenumbra(~,machine)
            % gaussian filter to model penumbra from (measured) machine output / see diploma thesis siggel 4.1.2
            matRad_cfg = MatRad_Config.instance();
            
            if isfield(machine.data,'penumbraFWHMatIso')
                penumbraFWHM = machine.data.penumbraFWHMatIso;
            else
                penumbraFWHM = 5;
                matRad_cfg.dispWarning('photon machine file does not contain measured penumbra width in machine.data.penumbraFWHMatIso. Assuming 5 mm.');
            end
            
            sourceFWHM = penumbraFWHM * machine.meta.SCD/(machine.meta.SAD - machine.meta.SCD);
            sigmaGauss = sourceFWHM / sqrt(8*log(2)); % [mm]
            
        end
        
        function [HUcube,cubeMatIx,cubeRho] = materialConversion(obj,ctGrid,doseGrid,ct)
            % conversion from HU to densities & materials
            HUcube      = cell(1,ct.numOfCtScen);
            cubeMatIx   = cell(1,ct.numOfCtScen);
            cubeRho     = cell(1,ct.numOfCtScen);
            
            % Create Material Density Cube
            material    = obj.setupMaterials();
            
            for s = 1:ct.numOfCtScen
                
                HUcube{s} =  matRad_interp3(ctGrid.x,ctGrid.y',ctGrid.z,ct.cubeHU{s}, ...
                    doseGrid.x,doseGrid.y',doseGrid.z,'nearest');
                
                % projecting out of bounds HU values where necessary
                if max(HUcube{s}(:)) > material{end,3}
                    matRad_cfg.dispWarning('Projecting out of range HU values');
                    HUcube{s}(HUcube{s}(:) > material{end,3}) = material{end,3};
                end
                if min(HUcube{s}(:)) < material{1,2}
                    matRad_cfg.dispWarning('Projecting out of range HU values');
                    HUcube{s}(HUcube{s}(:) < material{1,2}) = material{1,2};
                end
                
                % find material index
                cubeMatIx{s} = NaN*ones(doseGrid.dimensions,'int32');
                for i = size(material,1):-1:1
                    cubeMatIx{s}(HUcube{s} <= material{i,3}) = i;
                end
                
                % create an artificial HU lookup table
                hlut = [];
                for i = 1:size(material,1)
                    hlut = [hlut;material{i,2} material{i,4};material{i,3}-1e-10 material{i,5}]; % add eps for interpolation
                end
                
                cubeRho{s} = interp1(hlut(:,1),hlut(:,2),HUcube{s});
                
            end
        end
                

    end
    
    methods (Access = private)
        function material = setupMaterials(~)
            material = cell(4,5);
            material{1,1} = 'AIR700ICRU';
            material{1,2} = -1024;
            material{1,3} = -974;
            material{1,4} = 0.001;
            material{1,5} = 0.044;
            material{2,1} = 'LUNG700ICRU';
            material{2,2} = -974;
            material{2,3} = -724;
            material{2,4} = 0.044;
            material{2,5} = 0.302;
            material{3,1} = 'ICRUTISSUE700ICRU';
            material{3,2} = -724;
            material{3,3} = 101;
            material{3,4} = 0.302;
            material{3,5} = 1.101;
            material{4,1} = 'ICRPBONE700ICRU';
            material{4,2} = 101;
            material{4,3} = 1976;
            material{4,4} = 1.101;
            material{4,5} = 2.088;
            
        end
        
    end
    
    methods (Static)
        function compileOmpMCInterface(dest,omcFolder)
        % Compiles the ompMC interface (integrated as submodule)
        %
        % call
        %   matRad_OmpConfig.compileOmpMCInterface()
        %   matRad_OmpConfig.compileOmpMCInterface(dest)
        %   matRad_OmpConfig.compileOmpMCInterface(dest,sourceFolder)
        % if an object is instantiated, matRad_OmpConfig can be replaced by the 
        % object handle
        %
        % input:
        %   dest:           (optional) destination for mex file. Default: location
        %                   of this file
        %   sourceFolder:   (optional) path to ompMC . Default assumes its checked
        %                   out in the submodules folder of matRad
        %
        % References
        %
        %
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Copyright 2020 the matRad development team.
        %
        % This file is part of the matRad project. It is subject to the license
        % terms in the LICENSE file found in the top-level directory of this
        % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
        % of the matRad project, including this file, may be copied, modified,
        % propagated, or distributed except according to the terms contained in the
        % LICENSE file.
        %
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            matRad_cfg = MatRad_Config.instance();
            
            env = matRad_getEnvironment();
            
            if nargin < 1
                dest = fileparts(mfilename('fullpath'));
            end
            
            if nargin < 2
                omcFolder = [matRad_cfg.matRadRoot filesep 'submodules' filesep 'ompMC'];
            end
            
            sourceFolder = [omcFolder filesep 'src'];
            interfaceFolder = [omcFolder filesep 'ucodes' filesep 'omc_matrad'];
            
            mainFile = [interfaceFolder filesep 'omc_matrad.c'];
            
            addFiles = {'ompmc.c','omc_utilities.c','omc_random.c'};
            addFiles = cellfun(@(f) fullfile(sourceFolder,f),addFiles,'UniformOutput',false);
            
            addFiles = strjoin(addFiles,' ');
            
            if exist ('OCTAVE_VERSION','builtin')
                ccName = eval('mkoctfile -p CC');
            else
                myCCompiler = mex.getCompilerConfigurations('C','Selected');
                ccName = myCCompiler.ShortName;
            end
            
            %These settings have only been tested for MSVC and g++. You may need to adapt for other compilers
            if ~isempty(strfind(ccName,'MSVC')) %Not use contains(...) because of octave
                flags{1,1} = 'COMPFLAGS';
                flags{1,2} = '/openmp';
                flags{2,1} = 'OPTIMFLAGS';
                flags{2,2} = '/O2';
            else
                flags{1,1} = 'CFLAGS';
                flags{1,2} = '-std=gnu99 -fopenmp -O3';
                flags{2,1} = 'LDFLAGS';
                flags{2,2} = '-fopenmp';
                
            end
            
            includestring =  ['-I' sourceFolder];
            
            flagstring = '';
            
            %For Octave, the flags will be set in the environment, while they
            %will be parsed as string arguments in MATLAB
            for flag = 1:size(flags,1)
                if strcmp(env,'OCTAVE')
                    preFlagContent = eval(['mkoctfile -p ' flags{flag,1}]);
                    if ~isempty(preFlagContent)
                        preFlagContent = preFlagContent(1:end-1); %Strip newline
                    end
                    newContent = [preFlagContent ' ' flags{flag,2}];
                    setenv(flags{flag,1},newContent);
                    matRad_cfg.dispDebug('Set compiler flag %s to %s\n',flags{flag,1},newContent);
                else
                    flagstring = [flagstring flags{flag,1} '="$' flags{flag,1} ' ' flags{flag,2} '" '];
                end
            end
            
            mexCall = ['mex -largeArrayDims ' flagstring ' ' includestring ' ' mainFile ' ' addFiles];
            matRad_cfg.dispDebug('Compiler call: %s\n',mexCall);
            
            currDir = pwd;
            cd(dest);
            eval(mexCall);
            cd(currDir);
        end
    end
end

