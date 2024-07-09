classdef (Abstract) matRad_MonteCarloEngineAbstract < DoseEngines.matRad_DoseEngineBase
% matRad_MonteCarloEngineAbstract: abstract superclass for all dose calculation 
%   engines which are based on monte carlo calculation 
%   for more informations see superclass
%   DoseEngines.matRad_DoseEngineBase
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (SetAccess = public, GetAccess = public)
                         
        numHistoriesPerBeamlet; %number of histories per beamlet
        numHistoriesDirect;     %number of histories for a forward dose calculation
        
        outputMCvariance;       %boolean value to decide if variance information for the MC calculation should be scored

        relativeDosimetricCutOff;       %relative dosimetric cut-off (i.e., 1 - minimum relative dose-value to be stored)
    end
        
    methods
        
        function this = matRad_MonteCarloEngineAbstract(pln)   
            if nargin < 1
                pln = [];
            end

            % call superclass constructor
            this = this@DoseEngines.matRad_DoseEngineBase(pln);
        end

        function setDefaults(this)
            setDefaults@DoseEngines.matRad_DoseEngineBase(this);

            % create config instance
            matRad_cfg = MatRad_Config.instance();

            %set number of particles simulated per pencil beam
            this.numHistoriesPerBeamlet     = matRad_cfg.defaults.propDoseCalc.numHistoriesPerBeamlet;
            this.numHistoriesDirect         = matRad_cfg.defaults.propDoseCalc.numHistoriesDirect;
            this.outputMCvariance           = matRad_cfg.defaults.propDoseCalc.outputMCvariance;

            this.relativeDosimetricCutOff   = matRad_cfg.defaults.propDoseCalc.dosimetricLateralCutOff;
        end

        function resultGUI = calcDoseForward(this,ct,cst,stf,w)
            
            matRad_cfg = MatRad_Config.instance();
            if nargin < 5 && ~isfield([stf.ray],'weight')
                matRad_cfg.dispError('No weight vector available. Please provide w or add info to stf')
            end

            % copy bixel weight vector into stf struct
            if nargin == 5
                if sum([stf.totalNumOfBixels]) ~= numel(w) && ~isfield([stf.ray],'shapes')
                    matRad_cfg.dispError('weighting does not match steering information')
                end
                counter = 0;
                for i = 1:size(stf,2)
                    for j = 1:stf(i).numOfRays
                        for k = 1:stf(i).numOfBixelsPerRay(j)
                            counter = counter + 1;
                            stf(i).ray(j).weight(k) = w(counter);
                        end
                    end
                end
            else % weights need to be in stf!
                w = NaN*ones(sum([stf.totalNumOfBixels]),1);
                counter = 0;
                for i = 1:size(stf,2)
                    for j = 1:stf(i).numOfRays
                        for k = 1:stf(i).numOfBixelsPerRay(j)
                            counter = counter + 1;
                            w(counter) = stf(i).ray(j).weight(k);
                        end
                    end
                end
            end            
            
            %Set direct dose calculation and compute "dij"
            this.calcDoseDirect = true;
            dij = this.calcDose(ct,cst,stf);

            % hack dij struct
            dij.numOfBeams = 1;
            dij.beamNum = 1;

            % calculate cubes; use uniform weights here, weighting with actual fluence
            % already performed in dij construction
            resultGUI    = matRad_calcCubes(sum(w),dij);
            
            % remember original fluence weights
            resultGUI.w  = w;
        end
    end


    %% Deprecated properties
    properties (Dependent)
        relDoseCutOff;
    end
    
    methods
        function set.relDoseCutOff(this,relDoseCutOff_)
            this.relativeDosimetricCutOff = 1 - relDoseCutOff_;
            this.warnDeprecatedEngineProperty('relDoseCutOff','','relativeDosimetricCutOff');
        end
        function relDoseCutOff_ = get.relDoseCutOff(this)
            relDoseCutOff_ = 1 - this.relativeDosimetricCutOff;
            this.warnDeprecatedEngineProperty('relDoseCutOff','','relativeDosimetricCutOff');
        end
    end
    
end

