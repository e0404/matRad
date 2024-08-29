classdef matRad_brachyStfGenerator < matRad_StfGeneratorBase

    properties (Constant)
        name = 'brachyStfGen';
        shortName = 'brachyStfGen';
    end
    properties
        numRows
        numCols
        Xgrid
        Ygrid
    end
    
    methods 
        function this = matRad_brachyStfGenerator(pln)
            this@matRad_StfGeneratorBase(pln);
            matRad_cfg = MatRad_Config.instance();
            addpath(fullfile(matRad_cfg.matRadRoot));

            

            if ~isfield(pln, 'propStf')
                matRad_cfg.dispError('no applicator information in pln struct');
            end
         end
    end

    methods (Access = protected)

        function initializePatientGeometry(this,ct, cst, visMode)
            initializePatientGeometry@matRad_StfGeneratorBase(this,ct, cst, visMode)
            matRad_cfg = MatRad_Config.instance;



            pln = this.pln; 

           


            if ~isfield(pln.propStf, 'needle') || ~isfield(pln.propStf.needle, 'seedDistance') || ~isfield(pln.propStf.needle, 'seedsNo')
                matRad_cfg.dispError('Needle information missing in pln.propStf.needle');
            end
            if isfield(pln.propStf.template, 'type') && strcmp(pln.propStf.template.type, 'matrix')
                pln.propStf.template.activeNeedles = [0 0 0 1 0 1 0 1 0 1 0 0 0;... % 7.0
                                      0 0 1 0 1 0 0 0 1 0 1 0 0;... % 6.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 6.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 5.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 5.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 4.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 3.0
                                      1 0 1 0 1 0 1 0 1 0 1 0 1;... % 2.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 2.0
                                      1 0 1 0 1 0 0 0 0 0 1 0 1;... % 1.5
                                      0 0 0 0 0 0 0 0 0 0 0 0 0];   % 1.0
                                     %A a B b C c D d E e F f G
            elseif isfield(pln.propStf.template, 'type') && strcmp(pln.propStf.template.type, 'checkerboard')
                pln.propStf.template.activeNeedles = this.createCheckerboard(this.numRows,this.numCols);
            end

            if ~isfield(pln.propStf, 'template') || ~isfield(pln.propStf.template, 'activeNeedles')
                matRad_cfg.dispError('Template information missing in pln.propStf.template!');
            end


            if ~isfield(pln.propStf, 'templateRoot')
                matRad_cfg.dispError('TemplateRoot information missing in pln.propStf.templateRoot!');
            end

            if ~isa(pln.multScen, 'matRad_NominalScenario') && ~strcmp(pln.multScen, 'nomScen')
                matRad_cfg.dispError('Brachy Therapy does only work with a nominal scenario for now!');
            end

            if isempty(this.coordsX_vox) || isempty(this.coordsY_vox) || isempty(this.coordsZ_vox)
                matRad_cfg.dispWarning('coordsXYZ are empty, boundary cannot be computed.');   % they AREN'T EMPTY here problem is in brachyStfGen
                throw(MException('MATLAB:class:AbstractMember','coordsXYZ_vox need to be implemented!'));
            end
        end
        
        
        function stf = generateSourceGeometry(this, ct, cst, visMode)
            matRad_cfg = MatRad_Config.instance;
            pln = this.pln; 
            %translate to geometric coordinates and save in stf

            stf.targetVolume.Xvox = ct.x(this.coordsX_vox); % angabe in mm  
            stf.targetVolume.Yvox = ct.y(this.coordsY_vox);
            stf.targetVolume.Zvox = ct.z(this.coordsZ_vox);

            if strcmp(pln.propStf.template.type, 'matrix')
                
                pln.propStf.template.activeNeedles = [0 0 0 1 0 1 0 1 0 1 0 0 0;... % 7.0
                                      0 0 1 0 1 0 0 0 1 0 1 0 0;... % 6.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 6.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 5.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 5.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 4.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 3.0
                                      1 0 1 0 1 0 1 0 1 0 1 0 1;... % 2.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 2.0
                                      1 0 1 0 1 0 0 0 0 0 1 0 1;... % 1.5
                                      0 0 0 0 0 0 0 0 0 0 0 0 0];   % 1.0
                                     %A a B b C c D d E e F f G

                [row, col] = find(pln.propStf.template.activeNeedles);
                templX = col * pln.propStf.bixelWidth + pln.propStf.templateRoot(1) - (13 + 1)/2 * pln.propStf.bixelWidth;
                templY = row * pln.propStf.bixelWidth + pln.propStf.templateRoot(2) - (13 + 1)/2 * pln.propStf.bixelWidth;
                templZ = ones(size(col)) + pln.propStf.templateRoot(3);

            elseif strcmp(pln.propStf.template.type, 'checkerboard')

                [this.Xgrid, this.Ygrid] = meshgrid(stf.targetVolume.Xvox,stf.targetVolume.Yvox);

                xMin = min(stf.targetVolume.Xvox);
                xMax = max(stf.targetVolume.Xvox);
                yMin = min(stf.targetVolume.Yvox);
                yMax = max(stf.targetVolume.Yvox);

                checkerBoardX = length(xMin:pln.propStf.bixelWidth:xMax);
                checkerBoardY = length(yMin:pln.propStf.bixelWidth:yMax);


                pln.propStf.template.activeNeedles = this.createCheckerboard(checkerBoardX, checkerBoardY);
                [row, col] = find(pln.propStf.template.activeNeedles);
                templX = col * pln.propStf.bixelWidth + pln.propStf.templateRoot(1) - checkerBoardX / 2 * pln.propStf.bixelWidth;
                templY = row * pln.propStf.bixelWidth + pln.propStf.templateRoot(2) - checkerBoardY / 2 * pln.propStf.bixelWidth;
                templZ = ones(size(col)) + pln.propStf.templateRoot(3);

            else
                error('ActiveNeedles of your StfGenerator needs to be implemented!');
            end


            

            %% meta info from pln
            stf.radiationMode = pln.radiationMode;
            stf.numOfSeedsPerNeedle = pln.propStf.needle.seedsNo;
            stf.numOfNeedles = nnz(pln.propStf.template.activeNeedles);
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
            d = pln.propStf.needle.seedDistance;
            seedsNo = pln.propStf.needle.seedsNo;
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
                P = [TargX', TargY', TargZ'];

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

            seedPointsX = stf.seedPoints.x - pln.propStf.templateRoot(1);
            seedPointsY = stf.seedPoints.y - pln.propStf.templateRoot(2);
            seedPointsZ = stf.seedPoints.z - pln.propStf.templateRoot(3);

            targetVolumeX = stf.targetVolume.Xvox - pln.propStf.templateRoot(1);
            targetVolumeY = stf.targetVolume.Yvox - pln.propStf.templateRoot(2);
            targetVolumeZ = stf.targetVolume.Zvox - pln.propStf.templateRoot(3);

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
            checkerboardMatrix = zeros(numRows, numCols);

            for row = 1:numRows
                for col = 1:numCols
                    if mod(row + col, 2) == 0 
                        checkerboardMatrix(row, col) = 1;
                    else
                        checkerboardMatrix(row, col) = 0;
                    end
                end
            end
        end
    end
end



