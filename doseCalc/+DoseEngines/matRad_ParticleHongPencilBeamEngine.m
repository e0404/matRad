classdef matRad_ParticleHongPencilBeamEngine < DoseEngines.matRad_ParticlePencilBeamEngineAbstract
% matRad_ParticlePencilBeamEngineAbstractGaussian: 
%   Implements an engine for particle based dose calculation 
%   For detailed information see superclass matRad_DoseEngine
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the
% help edit

% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Constant)
           possibleRadiationModes = {'protons', 'carbon'}
           name = 'Particle Pencil-Beam';
    end
       
    methods 
        
        function this = matRad_ParticleHongPencilBeamEngine(pln)
            % Constructor
            %
            % call
            %   engine = DoseEngines.matRad_ParticleAnalyticalPencilBeamDoseEngine(ct,stf,pln,cst)
            %
            % input
            %   pln:                        matRad plan meta information struct
             
            this = this@DoseEngines.matRad_ParticlePencilBeamEngineAbstract(pln);
        end
        
    end

    methods (Access = protected)

        function bixel = calcParticleBixel(this,bixel)
            kernels = this.interpolateKernelsInDepth(bixel);
            
            %Lateral Component
            switch this.lateralModel
                case 'single'
                    %compute lateral sigma
                    sigmaSq = kernels.sigma.^2 + bixel.sigmaIniSq;
                    L = exp( -bixel.radialDist_sq ./ (2*sigmaSq))./ (2*pi*sigmaSq);
                case 'double'
                    % compute lateral sigmas
                    sigmaSqNarrow = kernels.sigma1.^2 + bixel.sigmaIniSq;
                    sigmaSqBroad  = kernels.sigma2.^2 + bixel.sigmaIniSq;
    
                    % calculate lateral profile
                    L_Narr =  exp( -bixel.radialDist_sq ./ (2*sigmaSqNarrow))./(2*pi*sigmaSqNarrow);
                    L_Bro  =  exp( -bixel.radialDist_sq./ (2*sigmaSqBroad ))./(2*pi*sigmaSqBroad );
                    L = (1-kernels.weight).*L_Narr + kernels.weight.*L_Bro;
                case 'multi'
                    sigmaSq = kernels.sigmaMulti.^2 + bixel.sigmaIniSq;
                    L = sum([1 - sum(kernels.weightMulti,2), kernels.weightMulti] .* exp(-radialDist_sq ./ (2*sigmaSq))./(2*pi*sigmaSq),2);
                otherwise
                    %Sanity check
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid Lateral Model');
            end
                        
            bixel.physicalDose = bixel.baseData.LatCutOff.CompFac * L .* kernels.Z;
            
            % check if we have valid dose values
            if any(isnan(bixel.physicalDose)) || any(bixel.physicalDose<0)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Error in particle dose calculation.');
            end

            if this.calcLET
                bixel.mLETDose = bixel.physicalDose.*kernels.LET;
            end
            
            if this.calcBioDose
                bixel.mAlphaDose = bixel.physicalDose;
                bixel.mSqrtBetaDose = bixel.physicalDose;
                %From matRad_calcLQParameter
                numOfTissueClass = size(bixel.baseData.alpha,2);
                for i = 1:numOfTissueClass
                    mask = bixel.vTissueIndex == i;
                    if any(mask)
                        bixel.mAlphaDose(mask) = bixel.mAlphaDose(mask) .* kernels.alpha(mask);
                        bixel.mSqrtBetaDose(mask)  = bixel.mSqrtBetaDose(mask) .* kernels.beta(mask);
                    end
                end
            end  
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
                checkModality = any(strcmp(DoseEngines.matRad_ParticleHongPencilBeamEngine.possibleRadiationModes, machine.meta.radiationMode));
                
                preCheck = checkBasic && checkModality;

                if ~preCheck
                    return;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end

            checkMeta = all(isfield(machine.meta,{'SAD','BAMStoIsoDist','LUT_bxWidthminFWHM','dataType'}));

            dataType = machine.meta.dataType;
            if strcmp(dataType,'singleGauss')
                checkData = all(isfield(machine.data,{'energy','depths','Z','peakPos','sigma','offset','initFocus'}));
            elseif strcmp(dataType,'doubleGauss')
                checkData = all(isfield(machine.data,{'energy','depths','Z','peakPos','weight','sigma1','sigma2','offset','initFocus'}));
            elseif strcmp(dataType,'multipleGauss')
                checkData = all(isfield(machine.data,{'energy','depths','Z','peakPos','weightMulti','sigmaMulti','offset','initFocus'}));
            else
                checkData = false;
            end
            
            available = checkMeta && checkData;
        end
    end
end

