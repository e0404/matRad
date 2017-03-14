function sb_resultGUI = matRad_singleBeamView(pln, cst, weights, viewBeamNum, dij, ct, stf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  calculate single beam dose distributions
%   this can be done either by providing a dij matrix or by calculating the
%   single beam dose directly from the corresponding spotweights
%  
%
% call
%   resultGUI = matRad_singleBeamView(pln, cst, resultGUI.w, viewBeamNum, dij)
%   or
%   resultGUI = matRad_singleBeamView(pln,cst,resultGUI.w,viewBeamNum,[],ct,stf)
%
% input
%   pln:         matRad pln struct
%   cst:         matRad cst struct
%   weights:     weights for beam spots (matRad resultGUI.w)
%   viewBeamNum: beam to be visualised, either integer for beamNum or 'all'
%                to calculate total dose
%
%   dij:         matRad dij struct (optional, single beam dose calculation
%                will use this, otherwise pass empty matrix)
%               
%                if you provide these parameters dose will be calculated
%                directly from spotweights:
%   ct:          matRad ct struct (optional, only needed if no dij provided)
%   stf:         matRad stf struct (optional, only needed if no dij provided)
%
% output
%  resultGUI:    matRad resultGUI struct containing single beam dose of
%                beam viewBeamNum
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(dij)
    % calculate single beam dose from dij
    if strcmp(viewBeamNum, 'all')
        % recalculate standard resultGUI for all beams
        fprintf('Calculating total dose \n');
        sb_resultGUI = matRad_calcCubes(weights,dij,cst,1);

    elseif ~isnumeric(viewBeamNum)
        error('wrong beam number format in matRad_singleBeamView \n');

    else
        % calculate single beam dose
        
        % adjust cst for single beams
        sb_cst = cst;
        for i=1:size(sb_cst,1)
            for j = 1:size(sb_cst{i,6},1)
                % biological dose splitting
                if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
                    ab = sb_cst{i,5}.alphaX / sb_cst{i,5}.betaX;
                    sb_cst{i,6}(j).dose = -0.5*ab +sqrt( 0.25*ab^2 + ...
                        sb_cst{i,6}(j).dose/pln.numOfBeams + ...
                        (sb_cst{i,6}(j).dose + ab));
                % physical dose splitting
                else
                    sb_cst{i,6}(j).dose = sb_cst{i,6}(j).dose/pln.numOfBeams;
                end
            end
        end

        for i = viewBeamNum
            fprintf(['Calculating Dose Beam ' num2str(viewBeamNum) '\n']);
           % columns in total dij for single beam
            sb_col = find(dij.beamNum == i);
            % construct dij for single beam
            sb_dij.numOfBeams = 1;
            sb_dij.numOfVoxels = dij.numOfVoxels;
            sb_dij.resolution = dij.resolution;
            sb_dij.numOfRaysPerBeam = dij.numOfRaysPerBeam(i);
            sb_dij.totalNumOfRays = sb_dij.numOfRaysPerBeam;
            sb_dij.totalNumOfBixels = size(sb_col, 1);
            sb_dij.dimensions = dij.dimensions;
            sb_dij.numOfScenarios = dij.numOfScenarios;
            sb_dij.bixelNum = dij.bixelNum(sb_col);
            sb_dij.rayNum = dij.rayNum(sb_col);
            sb_dij.beamNum = dij.beamNum(sb_col);
            sb_dij.physicalDose{1} = dij.physicalDose{1}(:, sb_col);
            if isfield(dij, 'RBE')
                sb_dij.RBE = dij.RBE;
            end
            if isfield(dij, 'mLETDose')
                sb_dij.mLETDose = dij.mLETDose(:, sb_col);
            end
            if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
                sb_dij.mAlphaDose{1} = dij.mAlphaDose{1}(:, sb_col);
                sb_dij.mSqrtBetaDose{1} = dij.mSqrtBetaDose{1}(:, sb_col);
            end

            % adjust pln to one beam only
            sb_pln = pln;
            sb_pln.gantryAngles = pln.gantryAngles(i);
            sb_pln.couchAngles = pln.couchAngles(i);

            sb_w = weights(sb_col);
            sb_resultGUI = matRad_calcCubes(sb_w,sb_dij,sb_cst,1);

            % keep full set of weights and for other beams
            sb_resultGUI.w = weights;
        end
    end
 
   
else
    % calculate single beam dose directly from spotweights
    
    if strcmp(viewBeamNum, 'all')
    % recalculate standard resultGUI for all beams
    fprintf('Calculating total dose \n');
    sb_resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,weights);

    elseif ~isnumeric(viewBeamNum)
        fprintf('Error: wrong beam number format in singleBeamView \n');

    else
        % calculate single beam dose

        % adjust cst for single beams
        sb_cst = cst;
        for i=1:size(sb_cst,1)
            for j = 1:size(sb_cst{i,6},1)
                % biological dose splitting
                if strcmp(pln.bioOptimization, 'LEMIV_effect') || ...
                        strcmp(pln.bioOptimization, 'LEMIV_RBExD')
                    ab = sb_cst{i,5}.alphaX / sb_cst{i,5}.betaX;
                    sb_cst{i,6}(j).dose = -0.5*ab +sqrt( 0.25*ab^2 + ...
                        sb_cst{i,6}(j).dose/pln.numOfBeams +...
                        (sb_cst{i,6}(j).dose + ab));
                % physical dose splitting
                else
                    sb_cst{i,6}(j).dose = sb_cst{i,6}(j).dose/pln.numOfBeams;
                end
            end
        end


        fprintf(['Calculating Dose Beam ' num2str(viewBeamNum) '\n']);

        % singlebeam stf
        sb_stf = stf(viewBeamNum);

        % singlebeam weights
        if viewBeamNum==1
            offset = 0;
        else
            offset = sum(stf(1:(viewBeamNum-1)).totalNumOfBixels);
        end
        sb_weights = weights(offset + (1:stf(viewBeamNum).totalNumOfBixels));

        % adjust pln to one beam only
        sb_pln = pln;
        sb_pln.isoCenter = pln.isoCenter(viewBeamNum,:);
        sb_pln.numOfBeams = 1;
        sb_pln.gantryAngles = pln.gantryAngles(viewBeamNum);
        sb_pln.couchAngles = pln.couchAngles(viewBeamNum);

        % calculate single beam dose
        sb_resultGUI = matRad_calcDoseDirect(ct,sb_stf,sb_pln,sb_cst,...
                                                    sb_weights);

        % keep full set of weights and for other beams
        sb_resultGUI.w = weights;

    end

end

end % eof