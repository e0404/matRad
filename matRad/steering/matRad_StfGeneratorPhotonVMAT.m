classdef matRad_StfGeneratorPhotonVMAT < matRad_StfGeneratorPhotonRayBixelAbstract

    properties (Constant)
        name = 'Photon VMAT stf Generator';
        shortName = 'PhotonVMAT';
        possibleRadiationModes = {'photons'};
    end 

    properties
        DAOGantryAngles;
        FMOGantryAngles;
        startingAngle;
        finishingAngle;
        continuousAperture;
    end
    
    methods 
        function this = matRad_StfGeneratorPhotonVMAT(pln)
            if nargin < 1
                pln = [];
            end

            this@matRad_StfGeneratorPhotonRayBixelAbstract(pln);
            
            % TODO VMAT: We should not use VMAT options but directly have
            % that as propStf, I think.
            if ~isempty(pln)
                this.startingAngle = pln.propOpt.VMAToptions.startingAngle;
                this.finishingAngle = pln.propOpt.VMAToptions.finishingAngle;
                this.continuousAperture = pln.propOpt.VMAToptions.continuousAperture;
            end

            if isempty(this.radiationMode)
                this.radiationMode = 'photons';
            end
        end            
    end

    methods (Access = protected)        
        function pbMargin = getPbMargin(this)
            pbMargin = this.bixelWidth;
        end

        function stf = generateSourceGeometry(this)
            stf = this.generateSourceGeometry@matRad_StfGeneratorPhotonRayBixelAbstract();
               
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo('Apply VMAT configuration to stf...\n');

            masterRayPosBEV = zeros(0,3);
            for i = 1:numel(stf)
                rayPosBEV = reshape([stf(i).ray(:).rayPos_bev]',3,stf(i).numOfRays)';
                masterRayPosBEV = union(masterRayPosBEV,rayPosBEV,'rows');
            end

            % masterRayPosBEV
            x = masterRayPosBEV(:,1);
            y = masterRayPosBEV(:,2);
            z = masterRayPosBEV(:,3);
            uniZ = unique(z);
            for j = 1:numel(uniZ)
                x_loc = x(z == uniZ(j));
                x_min = min(x_loc);
                x_max = max(x_loc);
                x = [x; (x_min:this.bixelWidth:x_max)'];
                y = [y; zeros((x_max-x_min)/this.bixelWidth+1,1)];
                z = [z; uniZ(j)*ones((x_max-x_min)/this.bixelWidth+1,1)];
            end

            SAD = this.machine.meta.SAD;
            
            masterRayPosBEV = [x,y,z];
            masterRayPosBEV = unique(masterRayPosBEV,'rows');
            masterTargetPointBEV = [2*masterRayPosBEV(:,1) SAD*ones(size(masterRayPosBEV,1),1) 2*masterRayPosBEV(:,3)];
            
            % post-processing function for VMAT
            matRad_cfg.dispInfo('VMAT stf post-processing (1/2)... ');

            numDAO = 1;
            DAODoseAngleBorders = zeros(2*numel(this.DAOGantryAngles),1);
            offset = 1;
            timeFacIndOffset = 1;
            
            for i = 1:length(this.gantryAngles)
                
                % Determine which FMO beam the current beam belongs to
                [~,stf(i).propVMAT.beamParentFMOIndex] = min(abs(this.FMOGantryAngles-this.gantryAngles(i)));
                stf(i).propVMAT.beamParentGantryAngle = this.FMOGantryAngles(stf(i).propVMAT.beamParentFMOIndex);
                stf(i).propVMAT.beamParentIndex = find(abs(this.gantryAngles - stf(i).propVMAT.beamParentGantryAngle) < 1e-6);
                
                % Indicate if this beam is to be included in DOA/FMO or not. All beams
                % are still considered in dose calc for objective function in DAO
                stf(i).propVMAT.FMOBeam = any(abs(this.FMOGantryAngles - this.gantryAngles(i)) < 1e-6);
                stf(i).propVMAT.DAOBeam = any(abs(this.DAOGantryAngles - this.gantryAngles(i)) < 1e-6);
                
                %% Determine different angle borders
                
                % doseAngleBorders are the angular borders over which dose is deposited
                if i == 1                    
                    stf(i).propVMAT.doseAngleBorders = [this.startingAngle (this.gantryAngles(i+1)+this.gantryAngles(i))/2];
                elseif i == length(this.gantryAngles)
                    
                    stf(i).propVMAT.doseAngleBorders = [(this.gantryAngles(i-1)+this.gantryAngles(i))/2 this.finishingAngle];
                else
                    
                    stf(i).propVMAT.doseAngleBorders = ([this.gantryAngles(i-1) this.gantryAngles(i+1)]+this.gantryAngles(i))/2;
                end
                
                stf(i).propVMAT.doseAngleBorderCentreDiff = [stf(i).gantryAngle-stf(i).propVMAT.doseAngleBorders(1) stf(i).propVMAT.doseAngleBorders(2)-stf(i).gantryAngle];
                stf(i).propVMAT.doseAngleBordersDiff = sum(stf(i).propVMAT.doseAngleBorderCentreDiff);
                
                %Assign beam to its Parent, either as child (optimized) or subchild
                %(interpolated)
                if stf(i).propVMAT.DAOBeam
                    DAODoseAngleBorders((offset):(offset+1)) = stf(i).propVMAT.doseAngleBorders;
                    offset= offset+2;
                    
                    if ~isfield(stf(stf(i).propVMAT.beamParentIndex).propVMAT,'beamChildrenGantryAngles') || isempty(stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenGantryAngles)
                        stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren = 0;
                        stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenGantryAngles = nan(1000,1);
                        stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenIndex = nan(1000,1);
                    end
                    
                    stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren = stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren+1;
                    stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenGantryAngles(stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren) = this.gantryAngles(i);
                    stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenIndex(stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren) = i;
                    
                    %optAngleBorders are the angular borders over which an optimized control point
                    %has influence
                    DAOIndex = find(abs(this.DAOGantryAngles - this.gantryAngles(i)) < 1e-8);
                    
                    if DAOIndex == 1
                        stf(i).propVMAT.DAOAngleBorders = [this.startingAngle (this.DAOGantryAngles(DAOIndex+1)+this.DAOGantryAngles(DAOIndex))/2];
                        
                        lastDAOIndex = i;
                        nextDAOIndex = find(abs(this.gantryAngles - this.DAOGantryAngles(DAOIndex+1)) < 1e-8);
                        
                        stf(i).propVMAT.lastDAOIndex = lastDAOIndex;
                        stf(i).propVMAT.nextDAOIndex = nextDAOIndex;
                    elseif DAOIndex == length(this.DAOGantryAngles)
                        stf(i).propVMAT.DAOAngleBorders = [(this.DAOGantryAngles(DAOIndex-1)+this.DAOGantryAngles(DAOIndex))/2 this.finishingAngle];
                        
                        stf(i).propVMAT.lastDAOIndex = find(abs(this.gantryAngles - this.DAOGantryAngles(DAOIndex-1)) < 1e-8);
                        stf(i).propVMAT.nextDAOIndex = i;
                    else
                        stf(i).propVMAT.DAOAngleBorders = ([this.DAOGantryAngles(DAOIndex-1) this.DAOGantryAngles(DAOIndex+1)]+this.DAOGantryAngles(DAOIndex))/2;
                        
                        lastDAOIndex = i;
                        nextDAOIndex = find(abs(this.gantryAngles - this.DAOGantryAngles(DAOIndex+1)) < 1e-8);
                        
                        stf(i).propVMAT.lastDAOIndex = lastDAOIndex;
                        stf(i).propVMAT.nextDAOIndex = nextDAOIndex;
                    end
                    
                    stf(i).propVMAT.DAOIndex = numDAO;
                    numDAO = numDAO+1;
                    
                    stf(i).propVMAT.DAOAngleBorderCentreDiff = [stf(i).gantryAngle-stf(i).propVMAT.DAOAngleBorders(1) stf(i).propVMAT.DAOAngleBorders(2)-stf(i).gantryAngle];
                    stf(i).propVMAT.DAOAngleBordersDiff = sum(stf(i).propVMAT.DAOAngleBorderCentreDiff);
                    
                    %This is the factor that relates the total time in the
                    %optimized arc sector to the total time in the current dose
                    %sector
                    stf(i).propVMAT.timeFacCurr =  stf(i).propVMAT.doseAngleBordersDiff./stf(i).propVMAT.DAOAngleBordersDiff;
                    
                    if this.continuousAperture
                        %These are the factors that relate the total time in the
                        %optimized arc sector to the total time in the previous and
                        %next dose sectors
                        stf(i).propVMAT.timeFac = zeros(1,3);
                        
                        stf(i).propVMAT.timeFac(1) = (stf(i).propVMAT.DAOAngleBorderCentreDiff(1)-stf(i).propVMAT.doseAngleBorderCentreDiff(1))/stf(i).propVMAT.DAOAngleBordersDiff;
                        stf(i).propVMAT.timeFac(2) = stf(i).propVMAT.timeFacCurr;
                        stf(i).propVMAT.timeFac(3) = (stf(i).propVMAT.DAOAngleBorderCentreDiff(2)-stf(i).propVMAT.doseAngleBorderCentreDiff(2))/stf(i).propVMAT.DAOAngleBordersDiff;
                        
                        % keep entries with a non-0 timeFac
                        delInd     = stf(i).propVMAT.timeFac == 0;
                        
                        % write timeFacInd
                        stf(i).propVMAT.timeFacInd          = [timeFacIndOffset-1 timeFacIndOffset timeFacIndOffset+1];
                        stf(i).propVMAT.timeFacInd(delInd)  = 0;
                        
                        % update offset
                        if delInd(3)
                            timeFacIndOffset = timeFacIndOffset+1;
                        else
                            timeFacIndOffset = timeFacIndOffset+2;
                        end
                        
                    else
                        %These are the factors that relate the total time in the
                        %optimized arc sector to the total time in the previous and
                        %next dose sectors
                        stf(i).propVMAT.timeFac = zeros(1,2);
                        
                        stf(i).propVMAT.timeFac(1) = stf(i).propVMAT.DAOAngleBorderCentreDiff(1)/stf(i).propVMAT.DAOAngleBordersDiff;
                        stf(i).propVMAT.timeFac(2) = stf(i).propVMAT.DAOAngleBorderCentreDiff(2)/stf(i).propVMAT.DAOAngleBordersDiff;
                    end
                    
                else
                    if ~isfield(stf(stf(i).propVMAT.beamParentIndex).propVMAT,'beamSubChildrenGantryAngles') || isempty(stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenGantryAngles)
                        stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren = 0;
                        stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenGantryAngles = nan(1000,1);
                        stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenIndex = nan(1000,1);
                    end
                    
                    stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren = stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren+1;
                    stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenGantryAngles(stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren) = this.gantryAngles(i);
                    stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenIndex(stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren) = i;
                    
                    stf(i).propVMAT.fracFromLastDAO = (this.gantryAngles(nextDAOIndex)-this.gantryAngles(i))./(this.gantryAngles(nextDAOIndex)-this.gantryAngles(lastDAOIndex));
                    stf(i).propVMAT.lastDAOIndex = lastDAOIndex;
                    stf(i).propVMAT.nextDAOIndex = nextDAOIndex;
                end
                
                
                if stf(i).propVMAT.FMOBeam
                    % FMOAngleBorders are the angular borders over which an optimized
                    % control point has influence
                    FMOIndex = find(abs(this.FMOGantryAngles - this.gantryAngles(i)) < 1e-8);
                    
                    if FMOIndex == 1
                        
                        stf(i).propVMAT.FMOAngleBorders = [this.startingAngle (this.FMOGantryAngles(FMOIndex+1)+this.FMOGantryAngles(FMOIndex))/2];
                    elseif FMOIndex == length(this.FMOGantryAngles)
                        
                        stf(i).propVMAT.FMOAngleBorders = [(this.FMOGantryAngles(FMOIndex-1)+this.FMOGantryAngles(FMOIndex))/2 this.finishingAngle];
                    else
                        
                        stf(i).propVMAT.FMOAngleBorders = ([this.FMOGantryAngles(FMOIndex-1) this.FMOGantryAngles(FMOIndex+1)]+this.FMOGantryAngles(FMOIndex))/2;
                    end
                    stf(i).propVMAT.FMOAngleBorderCentreDiff = [stf(i).gantryAngle-stf(i).propVMAT.FMOAngleBorders(1) stf(i).propVMAT.FMOAngleBorders(2)-stf(i).gantryAngle];
                    stf(i).propVMAT.FMOAngleBordersDiff = sum(stf(i).propVMAT.FMOAngleBorderCentreDiff);
                end
                
                %% transformation of union of rays
                
                stf(i).numOfRays = size(masterRayPosBEV,1);
                stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
                stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
                
                
                % source position in bev
                stf(i).sourcePoint_bev = [0 -SAD 0];
                
                % get (active) rotation matrix
                % transpose matrix because we are working with row vectors
                rotMat_vectors_T = transpose(matRad_getRotationMatrix(this.gantryAngles(i),this.couchAngles(i)));
                
                stf(i).sourcePoint = stf(i).sourcePoint_bev*rotMat_vectors_T;
                
                % Save ray and target position in lps system.
                for j = 1:stf(i).numOfRays
                    stf(i).ray(j).rayPos_bev = masterRayPosBEV(j,:);
                    stf(i).ray(j).targetPoint_bev = masterTargetPointBEV(j,:);
                    
                    stf(i).ray(j).rayPos      = stf(i).ray(j).rayPos_bev*rotMat_vectors_T;
                    stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev*rotMat_vectors_T;
                    
                    stf(i).ray(j).rayCorners_SCD = (repmat([0, this.machine.meta.SCD - SAD, 0],4,1)+ (this.machine.meta.SCD/SAD) * ...
                        [masterRayPosBEV(j,:) + [+stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                        masterRayPosBEV(j,:) + [-stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                        masterRayPosBEV(j,:) + [-stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2];...
                        masterRayPosBEV(j,:) + [+stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2]])*rotMat_vectors_T;
                end
                
                % loop over all rays to determine meta information for each ray
                stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
                stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
                
                for j = stf(i).numOfRays:-1:1
                    
                    % find appropriate energies for particles
                    if strcmp(stf(i).radiationMode,'photons')
                        
                        % book keeping for photons
                       stf(i).ray(j).energy = this.machine.data.energy;
                    else
                        matRad_cfg.dispError('Error generating stf struct: invalid radiation modality for VMAT.');
                    end
                end
                
                matRad_progress(i,length(this.gantryAngles));
            end
            
            
            %% final cleanup and calculation of factors we couldn't calc before
            matRad_cfg.dispInfo('VMAT post-processing (2/2)... ');
            
            for i = 1:length(this.gantryAngles)
                if stf(i).propVMAT.FMOBeam
                    %remove NaNs from beamChildren and beamSubChildren
                    if isfield(stf(i).propVMAT,'beamChildrenGantryAngles')
                        stf(i).propVMAT.beamChildrenGantryAngles(isnan(stf(i).propVMAT.beamChildrenGantryAngles)) = [];
                        stf(i).propVMAT.beamChildrenIndex(isnan(stf(i).propVMAT.beamChildrenIndex)) = [];
                    else
                        stf(i).propVMAT.numOfBeamChildren = 0;
                    end
                    if isfield(stf(i).propVMAT,'beamSubChildrenGantryAngles')
                        stf(i).propVMAT.beamSubChildrenGantryAngles(isnan(stf(i).propVMAT.beamSubChildrenGantryAngles)) = [];
                        stf(i).propVMAT.beamSubChildrenIndex(isnan(stf(i).propVMAT.beamSubChildrenIndex)) = [];
                    else
                        stf(i).propVMAT.numOfBeamSubChildren = 0;
                    end
                end
                
                if stf(i).propVMAT.DAOBeam
                    if this.continuousAperture
                        
                        stf(i).propVMAT.doseAngleDAO = ones(1,2);
                        
                        if sum(DAODoseAngleBorders == stf(i).propVMAT.doseAngleBorders(2)) > 1
                            %final dose angle is repeated
                            %do not count twice in optimization
                            stf(i).propVMAT.doseAngleDAO(2) = 0;
                        end
                    end
                end
                
                if ~stf(i).propVMAT.FMOBeam && ~stf(i).propVMAT.DAOBeam
                    
                    % for leaf position interpolation
                    stf(i).propVMAT.fracFromLastDAO_I = (stf(stf(i).propVMAT.nextDAOIndex).propVMAT.doseAngleBorders(1)-stf(i).propVMAT.doseAngleBorders(1))./(stf(stf(i).propVMAT.nextDAOIndex).propVMAT.doseAngleBorders(1)-stf(stf(i).propVMAT.lastDAOIndex).propVMAT.doseAngleBorders(2));
                    stf(i).propVMAT.fracFromLastDAO_F = (stf(stf(i).propVMAT.nextDAOIndex).propVMAT.doseAngleBorders(1)-stf(i).propVMAT.doseAngleBorders(2))./(stf(stf(i).propVMAT.nextDAOIndex).propVMAT.doseAngleBorders(1)-stf(stf(i).propVMAT.lastDAOIndex).propVMAT.doseAngleBorders(2));
                    stf(i).propVMAT.fracFromNextDAO_I = (stf(i).propVMAT.doseAngleBorders(1)-stf(stf(i).propVMAT.lastDAOIndex).propVMAT.doseAngleBorders(2))./(stf(stf(i).propVMAT.nextDAOIndex).propVMAT.doseAngleBorders(1)-stf(stf(i).propVMAT.lastDAOIndex).propVMAT.doseAngleBorders(2));
                    stf(i).propVMAT.fracFromNextDAO_F = (stf(i).propVMAT.doseAngleBorders(2)-stf(stf(i).propVMAT.lastDAOIndex).propVMAT.doseAngleBorders(2))./(stf(stf(i).propVMAT.nextDAOIndex).propVMAT.doseAngleBorders(1)-stf(stf(i).propVMAT.lastDAOIndex).propVMAT.doseAngleBorders(2));
                    
                    % for time interpolation
                    stf(i).propVMAT.timeFracFromLastDAO = (stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOAngleBorders(2)-stf(i).propVMAT.doseAngleBorders(1))./stf(i).propVMAT.doseAngleBordersDiff;
                    stf(i).propVMAT.timeFracFromNextDAO = (stf(i).propVMAT.doseAngleBorders(2)-stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOAngleBorders(2))./stf(i).propVMAT.doseAngleBordersDiff;
                    if stf(i).propVMAT.timeFracFromLastDAO > 1
                        stf(i).propVMAT.timeFracFromLastDAO = 1;
                    elseif stf(i).propVMAT.timeFracFromLastDAO < 0
                        stf(i).propVMAT.timeFracFromLastDAO = 0;
                    end
                    if stf(i).propVMAT.timeFracFromNextDAO > 1
                        stf(i).propVMAT.timeFracFromNextDAO = 1;
                    elseif stf(i).propVMAT.timeFracFromNextDAO < 0
                        stf(i).propVMAT.timeFracFromNextDAO = 0;
                    end
                end
                
                matRad_progress(i,length(this.gantryAngles));
            end

        end
    end

    methods (Static)
        function [available,msg] = isAvailable(pln,machine)
            % see superclass for information            
                   
            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            % Check superclass availability
            [available,msg] = matRad_StfGeneratorPhotonRayBixelAbstract.isAvailable(pln,machine);

            if ~available
                return;
            else
                available = false;
                msg = [];
            end
    
            %checkBasic
            try
                checkBasic = isfield(machine,'meta') && isfield(machine,'data');
    
                %check modality
                checkModality = any(strcmp(matRad_StfGeneratorPhotonVMAT.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_StfGeneratorPhotonVMAT.possibleRadiationModes, pln.radiationMode));
                
                %Sanity check compatibility
                if checkModality
                    checkModality = strcmp(machine.meta.radiationMode,pln.radiationMode);
                end
    
                preCheck = checkBasic && checkModality;
    
                if ~preCheck
                    return;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end

            available = preCheck;
        end
    end
end
