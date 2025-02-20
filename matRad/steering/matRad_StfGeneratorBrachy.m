classdef matRad_StfGeneratorBrachy < matRad_StfGeneratorBase

    properties (Constant)
        name = 'Basic Brachytherapy Template';
        shortName = 'SimpleBrachy';
        possibleRadiationModes = {'brachy'};
    end
    properties
        needle
        template
        bixelWidth
    end
    
    methods 
        function this = matRad_StfGeneratorBrachy(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_StfGeneratorBase(pln);
            if isempty(this.radiationMode)
                this.radiationMode = 'brachy';
            end
        end

        function setDefaults(this)
            this.setDefaults@matRad_StfGeneratorBase();
            this.machine = 'HDR';

            this.needle = struct(...
                'seedDistance',10,...
                'seedsNo',6);

            this.template.type = 'checkerboard';

            this.bixelWidth = 5;
        end
    end

    methods (Access = protected)
        function initialize(this)
            this.initialize@matRad_StfGeneratorBase();
            if ~isa(this.multScen,'matRad_NominalScenario') && ~strcmp(this.multScen,'nomScen')
                matRad_cfg.dispError('Brachy Therapy does only work with a single nominal scenario model!');
            end
        end        
        
        function stf = generateSourceGeometry(this)
            matRad_cfg = MatRad_Config.instance();

            if ~isfield(this.template,'root')
                this.template.root = matRad_getTemplateRoot(this.ct,this.cst);
            end

            if ~isfield(this.needle,'seedDistance') || ~isfield(this.needle,'seedsNo')
                matRad_cfg.dispError('Needle information missing!');
            end
            
                        
            stf.targetVolume.Xvox = this.voxTargetWorldCoords(:,1); 
            stf.targetVolume.Yvox = this.voxTargetWorldCoords(:,2);
            stf.targetVolume.Zvox = this.voxTargetWorldCoords(:,3);

            if ~isfield(this.template,'type')
                matRad_cfg.dispError('No template type specified!');
            end

            if strcmp(this.template.type, 'manual')
                %nothing to be done
                if ~isfield(this.template,'activeNeedles')
                    matRad_cfg.dispError('No active needle mask defined for template!');
                end

            elseif strcmp(this.template.type, 'checkerboard')

                %Bounding box
                xMin = min(stf.targetVolume.Xvox);
                xMax = max(stf.targetVolume.Xvox);
                yMin = min(stf.targetVolume.Yvox);
                yMax = max(stf.targetVolume.Yvox);

                %Calculate checkerboard size
                nCheckerboardX = ceil((xMax - xMin) ./ this.bixelWidth);
                nCheckerBoardY = ceil((yMax - yMin) ./ this.bixelWidth);

                %Create checkerboard
                this.template.activeNeedles = this.createCheckerboard(nCheckerboardX, nCheckerBoardY);



            else
                matRad_cfg.dispError('Template type ''%s'' invalid / not implemented!',this.template.type);
            end

            [row, col] = find(this.template.activeNeedles);            
            templX = col * this.bixelWidth + this.template.root(1) - (size(this.template.activeNeedles,1) + 1) / 2 * this.bixelWidth;
            templY = row * this.bixelWidth + this.template.root(2) - (size(this.template.activeNeedles,2) + 1) / 2 * this.bixelWidth;
            templZ = ones(size(col)) + this.template.root(3);


            %% meta info from pln
            stf.radiationMode = this.radiationMode;
            stf.numOfSeedsPerNeedle = this.needle.seedsNo;
            stf.numOfNeedles = sum(this.template.activeNeedles(:));
            stf.totalNumOfBixels = stf.numOfSeedsPerNeedle * stf.numOfNeedles; % means total number of seeds

            %% generate 2D template points
            % the template origin is set at its center. In the image coordinate system,
            % the center will be positioned at the bottom of the volume of interest.
            stf.template = [templX'; templY'; templZ'];
            stf.templateNormal = [0, 0, 1];

            %% generate seed positions
            % seed positions can be generated from needles, template and orientation
            % needles are assumed to go through the template vertically

            % needle position
            d = this.needle.seedDistance;
            seedsNo = this.needle.seedsNo;
            needleDist(1, 1, :) = d .* [0:seedsNo - 1]'; % 1x1xN Array with seed positions on needle
            needleDir = needleDist .* [0; 0; 1];
            seedPos_coord_need_seed = needleDir + stf.template;
            seedPos_need_seed_coord = shiftdim(seedPos_coord_need_seed, 1);
            % the output array has the dimensions (needleNo, seedNo, coordinates)
            X = seedPos_need_seed_coord(:, :, 1);
            Y = seedPos_need_seed_coord(:, :, 2);
            Z = seedPos_need_seed_coord(:, :, 3);

            stf.seedPoints.x = reshape(X, 1, []);
            stf.seedPoints.y = reshape(Y, 1, []);
            stf.seedPoints.z = reshape(Z, 1, []);

            matRad_cfg.dispInfo('Processing completed: %d%%\n', 100);

            %% visualize results of visMode is nonzero
            if this.visMode > 0
                clf
                SeedPoints = plot3(stf.seedPoints.x, stf.seedPoints.y, stf.seedPoints.z, '.', 'DisplayName', 'seed points', 'Color', 'black', 'markersize', 5);
                title('3D Visualization of seed points')
                xlabel('X (left) [mm]')
                ylabel('Y (posterior) [mm]')
                zlabel('Z (superior) [mm]')
                hold on

                % plot 3d VOI points
                TargX = stf.targetVolume.Xvox;
                TargY = stf.targetVolume.Yvox;
                TargZ = stf.targetVolume.Zvox;
                % Prostate = plot3(TargX, TargY, TargZ, '.', 'Color', 'b', 'DisplayName', 'prostate');

                % Prepare points for boundary calculation
                P = [TargX, TargY, TargZ];

                if ~isempty(P)
                    % Determine the environment
                    if matRad_cfg.isOctave
                        % Octave environment
                        [uni, ~] = sort(unique(TargX));
                        n = length(uni);
                        outline = zeros(2 * n, 2);

                        for i = 1:n
                            y_list = TargY(TargX == uni(i));
                            y_max = max(y_list);
                            y_min = min(y_list);
                            outline(i, :) = [uni(i), y_max];
                            outline(2 * n - i + 1, :) = [uni(i), y_min];
                        end

                        % Plot the points and the computed outline
                        figure;
                        plot(TargX, TargY, 'b+', 'DisplayName', 'VOI Points');
                        hold on;
                        plot(outline(:, 1), outline(:, 2), 'g-', 'LineWidth', 3, 'DisplayName', 'Computed Outline');

                        % Calculate the area enclosed by the outline
                        area = polyarea(outline(:, 1), outline(:, 2));
                        disp(['Polygon area: ', num2str(area)]);

                        hold off;
                    else
                        % MATLAB environment
                        k = boundary(P, 1);
                        trisurf(k, P(:, 1), P(:, 2), P(:, 3), 'FaceColor', 'red', 'FaceAlpha', 0.1, 'LineStyle', 'none')
                    end
                else
                    matRad_cfg.dispWarning('Target volume points are empty, boundary cannot be computed.');
                    throw(MException('MATLAB:class:AbstractMember','Target Volumes need to be implemented!'));
                end
            end

            % throw warning if seed points are more than twice the central
            % distance outside the TARGET volume or if no seed points are in the
            % target volume

            seedPointsX = stf.seedPoints.x - this.template.root(1);
            seedPointsY = stf.seedPoints.y - this.template.root(2);
            seedPointsZ = stf.seedPoints.z - this.template.root(3);

            targetVolumeX = stf.targetVolume.Xvox - this.template.root(1);
            targetVolumeY = stf.targetVolume.Yvox - this.template.root(2);
            targetVolumeZ = stf.targetVolume.Zvox - this.template.root(3);

            if any([max(seedPointsX) >= 4 * max(targetVolumeX), ...
                    min(seedPointsX) <= 4 * min(targetVolumeX), ...
                    max(seedPointsY) >= 4 * max(targetVolumeY), ...
                    min(seedPointsY) <= 4 * min(targetVolumeY), ...
                    max(seedPointsZ) >= 4 * max(targetVolumeZ), ...
                    min(seedPointsZ) <= 4 * min(targetVolumeZ)])
                matRad_cfg.dispWarning('Seeds far outside the target volume');
            end

            if all([max(targetVolumeX) <= min(seedPointsX), min(targetVolumeX) >= max(seedPointsX), ...
                    max(targetVolumeY) <= min(seedPointsY), min(targetVolumeY) >= max(seedPointsY), ...
                    max(targetVolumeZ) <= min(seedPointsZ), min(targetVolumeZ) >= max(seedPointsZ)])
                matRad_cfg.dispWarning('No seed points in VOI')
            end
        end
    end
    methods
        function checkerboardMatrix = createCheckerboard(this, numRows, numCols)
            % Initialize checkerboard matrix
            checkerboardMatrix = false(numRows, numCols);

            % Fill checkerboard matrix
            ix = 1:numel(checkerboardMatrix);
            [i,j] = ind2sub(size(checkerboardMatrix), ix);
            checkerboardMatrix(ix) = mod(i + j, 2) == 0;
        end
    end

    methods (Static)
        function [available,msg] = isAvailable(pln,machine)
            % see superclass for information
    
            msg = [];
            available = false;
    
            if nargin < 2
                machine = matRad_loadMachine(pln);
            end
    
            %checkBasic
            try
                checkBasic = isfield(machine,'meta') && isfield(machine,'data');
    
                %check modality
                checkModality = any(strcmp(matRad_StfGeneratorBrachy.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_StfGeneratorBrachy.possibleRadiationModes, pln.radiationMode));
                
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



