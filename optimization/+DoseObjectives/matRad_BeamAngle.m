classdef matRad_BeamAngle < DoseObjectives.matRad_DoseObjective
% matRad_BeamAngle minimizes spatial dose averaged LET within OARs (and target volume) 
% (so maximizes dose) by maximizing the angle between treatment fields (beams) to obtain a more
% robust plan with respect to RBE uncertainties
% First BAO then FMO
% (Note that after the optimum beam configuration has been obtained, the final IMRT plan 
% using the optimum beam configuration is obtained by reoptimizing fluences using the full accurate dose calculation.)
%
% References 
% Exploration and application of phenomenological RBE models for proton therapy
%   DOI 10.1088/1361-6560/aad9db
% Development of methods for beam angle optimization for IMRT using an accelerated exhaustive search strategy
%   doi:10.1016/j.ijrobp.2004.06.007
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Constant)
        name = 'Beam Angle';
        parameterNames = {'d^{min}', 'a^{min}', 'V^{max}', 'i^{max}'};
        parameterTypes = {'dose', 'dose', 'numeric', 'numeric'};
    end
    
    properties
        parameters = {60}; % ?
        penaltyT = 1;
        penaltyOAR = 2;
    end
    
    methods
        function obj = matRad_BeamAngle(penaltyT,penaltyOAR,dMin,anglebeamlet,vMaxPercent,index)
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penaltyT)
                inputStructT = penaltyT;
                initFromStructT = true;
            else
                initFromStructT = false;
                inputStructT = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DoseObjectives.matRad_DoseObjective(inputStructT);
            
            %now handle initialization from other parameters
            if ~initFromStructT
                if nargin >= 6 && isscalar(index)
                    obj.parameters{5} = index;
                end

                if nargin >= 5 && isscalar(vMaxPercent)
                    obj.parameters{4} = vMaxPercent;
                end

                if nargin >= 4 && isscalar(anglebeamlet)
                    obj.parameters{3} = anglebeamlet;
                end

                if nargin >= 3 && isscalar(dMin)
                    obj.parameters{2} = dMin;
                end

                if nargin == 2 && isscalar(penaltyOAR)
                    obj.parameters{1} = penaltyOAR;
                end

                if nargin >= 1 && isscalar(penaltyT)
                    obj.penalty = penaltyT;
                end
            end
        end
        
        %% Calculates the Objective Function value
        function fDose = computeDoseObjectiveFunction(obj,dose)
            % get reference Volume
            refVol = obj.parameters{4}/100;

            % calc deviation
            deviation = dose - obj.parameters{1};

            % calc d_ref2: V(d_ref2) = refVol
            d_ref2 = matRad_calcInversDVH(refVol,dose);

            
            deviation(dose < obj.parameters{1} | dose > d_ref2) = 0;

            V1 = (1/numel(dose))*(deviation'*deviation);
            
            i = index;
            Hi2 = V1(i) - V1(2);

            if Hi2 >= 0
                H2 = 1;
            else
                H2 = 0;
            end

            H2i = V1(2) - V1(i);

            if H2i >= 0
                H2i = 1;
            else
                H2i = 0;
            end

            Hi = V1(1) - V1(i);

            if Hi >= 0
                Hi = 1;
            else
                Hi = 0;
            end

            H1 = V1(i) - V1(1);

            if H1 >= 0
                H = 1;
            else
                H = 0;
            end

            % calculate objective function
            % for targets:
            fT = (1/NT)*obj.penalty*sum(H2*Hi*H1^2);
            % for OARs:
            fOARs = (1/NOARs)*obj.parameters{1}*sum(H*H2i*H1^2);

            fDose = fT + fOARs;

        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
            % get reference Volume
            refVol = obj.parameters{4}/100;

            % calc deviation
            deviation = dose - obj.parameters{1};

            % calc d_ref2: V(d_ref2) = refVol
            d_ref2 = matRad_calcInversDVH(refVol,dose);

            
            deviation(dose < obj.parameters{1} | dose > d_ref2) = 0;

            V1 = (1/numel(dose))*(deviation'*deviation);
            
            i = index;
            Hi2 = V1(i) - V1(2);

            if Hi2 >= 0
                H2 = 1;
            else
                H2 = 0;
            end

            H2i = V1(2) - V1(i);

            if H2i >= 0
                H2i = 1;
            else
                H2i = 0;
            end

            Hi = V1(1) - V1(i);

            if Hi >= 0
                Hi = 1;
            else
                Hi = 0;
            end

            H1 = V1(i) - V1(1);

            if H1 >= 0
                H = 1;
            else
                H = 0;
            end

            % calculate objective function
            % for targets:
            fT = (1/NT)*obj.penalty*sum(H2*Hi*H1^2);
            % for OARs:
            fOARs = (1/NOARs)*obj.parameters{1}*sum(H*H2i*H1^2);
            
            % fDose = fT + fOARs

            % calculate delta
            fDoseGrad = 2/numel(dose) * underdose;
        end
    end
    
end