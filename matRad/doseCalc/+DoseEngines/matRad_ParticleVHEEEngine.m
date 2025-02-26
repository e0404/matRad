classdef matRad_ParticleVHEEEngine < DoseEngines.matRad_ParticlePencilBeamEngineAbstract
% matRad_ParticleVHEEEngine: 
%   Implements an engine for VHEE dose calculation
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2025 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Constant)
           possibleRadiationModes = {'VHEE'}
           name = 'VHEE Particle Pencil-Beam';
           shortName = 'VHEE';
    end
       
    methods 
        
        function this = matRad_ParticleVHEEEngine(pln)
            % Constructor
            %
            % call
            %   engine = DoseEngines.matRad_ParticleAnalyticalPencilBeamDoseEngine(ct,stf,pln,cst)
            %
            % input
            %   pln:                        matRad plan meta information struct

            if nargin < 1
                pln = [];
            end
             
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
                    L = sum([1 - sum(kernels.weightMulti,2), kernels.weightMulti] .* exp(-bixel.radialDist_sq ./ (2*sigmaSq))./(2*pi*sigmaSq),2);
                case 'singleXY'
                    %compute lateral sigma in both directions
                    sigmaSq_x = kernels.sigmaX.^2 + bixel.sigmaIniSq;
                    sigmaSq_y = kernels.sigmaY.^2 + bixel.sigmaIniSq;
                    sigma_x = sqrt(sigmaSq_x);
                    sigma_y = sqrt(sigmaSq_y);

                    L = exp( - (bixel.radialDist_sq.^2)./(2*sigmaSq_x) - (bixel.radialDist_sq.^2)./(2*sigmaSq_y) ) ./(2*pi*sigma_x.*sigma_y);
                    
                otherwise
                    %Sanity check
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid Lateral Model');
            end
                        
            if length(bixel.baseData.LatCutOff.CompFac) > 1
                bixel.baseData.LatCutOff.CompFac = matRad_interp1(bixel.baseData.LatCutOff.depths', bixel.baseData.LatCutOff.CompFac', bixel.radDepths);
            end
            bixel.physicalDose = bixel.baseData.LatCutOff.CompFac .* L .* kernels.Z;
            
            % check if we have valid dose values
            if any(isnan(bixel.physicalDose)) || any(bixel.physicalDose<0)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Error in particle dose calculation.');
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
                checkModality = any(strcmp(DoseEngines.matRad_ParticleVHEEEngine.possibleRadiationModes, machine.meta.radiationMode));
                
                preCheck = checkBasic && checkModality;

                if ~preCheck
                    return;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end

            checkMeta = all(isfield(machine.meta,{'SAD','BAMStoIsoDist','dataType'}));
            
            %Superseded names from older machine file versions
            checkMeta = checkMeta && any(isfield(machine.meta,{'LUTspotSize','LUT_bxWidthminFWHM'}));

            dataType = machine.meta.dataType;
            if strcmp(dataType,'singleGaussXY')
                checkData = all(isfield(machine.data,{'energy','depths','Z','offset','initFocus','sigmaXY'}));
            elseif strcmp(dataType,'doubleGauss')
                checkData = all(isfield(machine.data,{'energy','depths','Z','weight','sigma1','sigma2','offset','initFocus'}));
            else
                checkData = false;
            end
            
            available = checkMeta && checkData;
        end
    end
end

