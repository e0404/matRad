classdef matRad_UIData < handle
    % matRad_UIData implements a class that allows easy storing of
    % variables related to the pareto Navigation
    %
    % References
    %   -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team. 
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
        wRef %weight vector of last calculated plan
        fRef %objective function values of last calculated plan
        fIndsAll % stores all objective function values
        fIndsRed % stores only the "reduced points"
        upboundInit
        upboundRed
        upboundSlider
        lowboundSliderInit
        lowboundSlider
        linestyle
    end

    methods
        function obj = matRad_UIData(wRef,fRef,fInds)
            obj.wRef = wRef; 
            obj.fRef = fRef; 
            obj.fIndsAll = fInds;
            obj.fIndsRed = fInds;
            obj.upboundInit = max(fInds,[],1);
            obj.upboundRed = max(fInds,[],1);
            obj.upboundSlider= max(fInds,[],1);
            obj.lowboundSliderInit = min(fInds,[],1);
            obj.lowboundSlider = min(fInds,[],1);
            obj.linestyle = 2;
        end       
   
        function [sliderLowBound,sliderUpBound] = restrictObjective(obj,i,bound)
            obj.upboundRed(i) = bound;
            obj.fIndsRed = obj.fIndsRed(obj.fIndsRed(:,i) <= bound,:);
            obj.upboundSlider= max([obj.fIndsRed;obj.fRef],[],1);
            obj.upboundSlider(i) = bound;
            obj.lowboundSlider = min([obj.fIndsRed;obj.fRef],[],1);
            sliderLowBound = obj.lowboundSlider;
            sliderUpBound = obj.upboundSlider;
        end

        function releaseObjectiveBounds(obj)
            obj.upboundRed = obj.upboundInit;
            obj.fIndsRed = obj.fIndsAll;    
            obj.upboundSlider = obj.upboundInit;
            obj.lowboundSlider = obj.lowboundSliderInit;
        end

    end
end