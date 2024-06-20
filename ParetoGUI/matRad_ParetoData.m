classdef matRad_ParetoData < handle
    % matRad_ParetoData implements a class that allows easy storing of
    % variables related to the pareto Navigation
    %
    % References
    %   -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2024 the matRad development team. 
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
        currentWeights
        allWeights
        currentPoint %objective function values of last calculated plan
        allPoints % matrix storing all objective function values
        availablePoints % stores only the points not blocked by boundaries
        initialUpbound
        availableUpbound
        upboundSlider
        lowboundSliderInit
        lowboundSlider
        linestyle
    end

    methods
        function this = matRad_ParetoData(allPoints,allWeights)
            %initialize first set
            this.allWeights = allWeights;
            this.allPoints = allPoints;
            this.availablePoints = allPoints;
            this.initialUpbound = max(allPoints,[],1);
            this.availableUpbound = max(allPoints,[],1);
            this.upboundSlider = max(allPoints,[],1);
            this.lowboundSliderInit = min(allPoints,[],1);
            this.lowboundSlider = min(allPoints,[],1);
            this.linestyle = 2;

            %calculate an initial point based on the shortes distance to the
            %polyhedral center
            center = mean(this.allPoints,1);
            shiftedPoints = this.allPoints-center;
            distToCenter = sum(shiftedPoints.^2,2);
            [~,idx] = min(distToCenter);
            this.currentPoint = this.allPoints(idx,:);
            this.currentWeights = this.allWeights(:,idx);

        end       

   
        function [sliderLowBound,sliderUpBound] = restrictObjective(this,i,bound)
            %restrict the upper value of a specific slider
            this.availableUpbound(i) = bound;
            this.availablePoints = this.availablePoints(this.availablePoints(:,i) <= bound,:);
            this.upboundSlider= max([this.availablePoints;this.currentPoint'],[],1);
            this.upboundSlider(i) = bound;
            this.lowboundSlider = min([this.availablePoints;this.currentPoint'],[],1);
            sliderLowBound = this.lowboundSlider;
            sliderUpBound = this.upboundSlider;
        end

        function releaseObjectiveBounds(this)
            %reset objective bounds
            this.availableUpbound = this.initialUpbound;
            this.availablePoints = this.allPoints;    
            this.upboundSlider = this.initialUpbound;
            this.lowboundSlider = this.lowboundSliderInit;
        end

    end
end