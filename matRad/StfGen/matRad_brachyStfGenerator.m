classdef matRad_brachyStfGenerator < matRad_StfGeneratorBase

    properties (Constant)
        name = 'brachyStfGen';
        shortName = 'brachyStfGen';
    end 
    
    properties (Access = protected)
        SeedPoints
    end
    
    methods 
        function this = matRad_brachyStfGenerator(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_StfGeneratorBase(pln);
            this.radiationMode = pln.radiationMode;

            matRad_cfg = MatRad_Config.instance();
            addpath(fullfile(matRad_cfg.matRadRoot));
            
            if ~isfield(pln, 'propStf')
                matRad_cfg.dispError('no applicator information in pln struct');
            end
         end
    end

    methods 

        function stf = generate(this, ct, cst)
            stf = generate@matRad_StfGeneratorBase(this, ct, cst);
            this.initializePatientGeometry(ct, cst);
            stf = this.generateSourceGeometry(ct, cst);
        end
    end

    methods (Access = protected)
        function initializePatientGeometry(this, ct, cst)
            % Initialize the patient geometry
            initializePatientGeometry@matRad_StfGeneratorBase(this, ct, cst);

            pln = this.pln; % Ensure `pln` is accessible within this method

            if ~isfield(pln.propStf, 'needle') || ~isfield(pln.propStf.needle, 'seedDistance') || ~isfield(pln.propStf.needle, 'seedsNo')
                matRad_cfg.dispError('Needle information missing in pln.propStf.needle');
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
        end
        
        function pbMargin = getPbMargin(this)
            pbMargin = getPbMargin@matRad_StfGeneratorBase(this);
        end
        
        function stf = generateSourceGeometry(this, ct, cst)
            stf = generateSourceGeometry@matRad_StfGeneratorBase(this);

            pln = this.pln; % Ensure `pln` is accessible within this method

            %translate to geometric coordinates and save in stf

            stf.targetVolume.Xvox = ct.x(this.coordsX_vox); % angabe in mm
            stf.targetVolume.Yvox = ct.y(this.coordsY_vox);
            stf.targetVolume.Zvox = ct.z(this.coordsZ_vox);

            %% meta info from pln
            stf.radiationMode = pln.radiationMode;
            stf.numOfSeedsPerNeedle = pln.propStf.needle.seedsNo;
            stf.numOfNeedles = nnz(pln.propStf.template.activeNeedles);
            stf.totalNumOfBixels = stf.numOfSeedsPerNeedle * stf.numOfNeedles; % means total number of seeds

            %% generate 2D template points
            % the template origin is set at its center. In the image coordinate system,
            % the center will be positioned at the bottom of the volume of interest.
            [row, col] = find(pln.propStf.template.activeNeedles);
            templX = col * pln.propStf.bixelWidth + pln.propStf.templateRoot(1) - (13 + 1) / 2 * pln.propStf.bixelWidth;
            templY = row * pln.propStf.bixelWidth + pln.propStf.templateRoot(2) - (13 + 1) / 2 * pln.propStf.bixelWidth;
            templZ = ones(size(col)) + pln.propStf.templateRoot(3);

            stf.template = [templX'; templY'; templZ'];
            stf.templateNormal = [0, 0, 1];

            %% generate seed positions
            % seed positions can be generated from needles, template and orientation
            % needles are assumed to go through the template vertically

            % needle position
            d = pln.propStf.needle.seedDistance;
            seedsNo = pln.propStf.needle.seedsNo;
            needleDist(1, 1, :) = d .* (0:seedsNo - 1)'; % 1x1xN Array with seed positions on needle
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
            if visMode > 0
                clf
                this.SeedPoints = plot3(stf.seedPoints.x, stf.seedPoints.y, stf.seedPoints.z, '.', 'DisplayName', 'seed points', 'Color', 'black', 'markersize', 5);
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
            end

            % throw warning if seed points are more than twice the central
            % distance outside the TARGET volume or if no seed points are in the
            % target volume

            if (max(stf.seedPoints.x - pln.propStf.templateRoot(1)) >= 4 * max(stf.targetVolume.Xvox - pln.propStf.templateRoot(1)) || ...
                    min(stf.seedPoints.x - pln.propStf.templateRoot(1)) <= 4 * min(stf.targetVolume.Xvox - pln.propStf.templateRoot(1)) || ...
                    max(stf.seedPoints.y - pln.propStf.templateRoot(2)) >= 4 * max(stf.targetVolume.Yvox - pln.propStf.templateRoot(2)) || ...
                    min(stf.seedPoints.y - pln.propStf.templateRoot(2)) <= 4 * min(stf.targetVolume.Yvox - pln.propStf.templateRoot(2)) || ...
                    max(stf.seedPoints.z - pln.propStf.templateRoot(3)) >= 4 * max(stf.targetVolume.Zvox - pln.propStf.templateRoot(3)) || ...
                    min(stf.seedPoints.z - pln.propStf.templateRoot(3)) <= 4 * min(stf.targetVolume.Zvox - pln.propStf.templateRoot(3)))
                matRad_cfg.dispWarning('Seeds far outside the target volume');
            end
            if (max(stf.targetVolume.Xvox) <= min(stf.seedPoints.x) || min(stf.targetVolume.Xvox) >= max(stf.seedPoints.x) || ...
                    max(stf.targetVolume.Yvox) <= min(stf.seedPoints.y) || min(stf.targetVolume.Yvox) >= max(stf.seedPoints.y) || ...
                    max(stf.targetVolume.Zvox) <= min(stf.seedPoints.z) || min(stf.targetVolume.Zvox) >= max(stf.seedPoints.z))
                matRad_cfg.dispWarning('no seed points in VOI')
            end
        end
    end
end


