classdef MatRadScenario < handle

    properties (Access = private)
        dvhQiReady = false
        % during calcQiDVH
        cst
        pln
    end
    
    properties (SetAccess = private)
        weight;
        subIx;
        shift;
        relRangeShift;
        absRangeShift;
        isoShift;
        ctShiftIdentifier;
        radiationQuantity;
        doseCubeDim;
        
        % computed properties which need to be recalculated when dose
        % changes. is not in Dependent properties because calculation is
        % expensive
        dvh;
        qi;
    end

    properties (SetObservable, AbortSet, SetAccess = private)
        % container of scenarios
        dose;
    end
    
    properties (Dependent=true, SetAccess = private)
        doseLin;
    end
    
    methods
        function obj = MatRadScenario(dose, subIx, radiationQuantity, ctShiftIdentifier, shift, ...
                relRangeShift, absRangeShift, doseCubeDim, weight)
            
            obj.dose = dose;
            if ndims(dose) == 3
                obj.subIx = subIx;
                obj.doseCubeDim = size(obj.dose);
            else
                obj.doseCubeDim = doseCubeDim;
                obj.subIx = subIx;
            end
            
            if exist('weight','var') && ~isempty(weight)
                obj.weight = weight;
            else
                obj.weight = NaN;
            end
            obj.radiationQuantity = radiationQuantity;
            
            obj.ctShiftIdentifier = ctShiftIdentifier;
            obj.shift = shift;
            obj.relRangeShift = relRangeShift;
            obj.absRangeShift = absRangeShift;
            
            addlistener(obj,'dose','PostSet',@obj.handleChangeOfDose);
        end % eof constructor
        
        function obj = cumulateDose(obj, dose, ctShiftIdentifier, shift, absRangeShift, relRangeShift)
            obj.dose = obj.dose + dose;
            obj.ctShiftIdentifier = [obj.ctShiftIdentifier ctShiftIdentifier];
            obj.shift = [obj.shift; shift];
            obj.absRangeShift = [obj.absRangeShift absRangeShift];
            obj.relRangeShift = [obj.relRangeShift relRangeShift];
        end
        
        function obj = calcQiDVH(obj, cst, pln, dvhType, doseGrid, refGy, refVol)
            obj.cst = cst;
            obj.pln = pln;
            obj.dvh = matRad_calcDVH(cst, obj.dose, dvhType, doseGrid);
            obj.qi = matRad_calcQualityIndicators(cst, pln, obj.dose, refGy, refVol);
            obj.dvhQiReady = true;
        end
        
        function plotQiDVH(obj)
            if obj.dvhQiReady
                subplot(2,1,1)
                matRad_showDVH(obj.dvh,obj.cst,obj.pln);
                subplot(2,1,2)
                matRad_showQualityIndicators(obj.qi);
            else
                error('You need to compute first.');
            end
        end
        
        function plotDoseSlice(obj, ax, ct, plane, slice, doseWindow, legendOn)
            if ~exist('doseWindow', 'var') || isempty(doseWindow)
                doseWindow = [0 max(obj.dose(:))];
            end
            if ~exist('legendOn', 'var') || isempty(legendOn)
                legendOn = false;
            end
            colorMapLabel = obj.radiationQuantity;
            matRad_plotSliceWrapper(ax,ct,obj.cst,1,obj.dose,plane,slice,[],[],jet,[],doseWindow,[],[],colorMapLabel,legendOn);
        end
        
        function singleStructDVH = getSingleStructDVH(obj, voi)
            ix = obj.getIndexOfStruct(obj.dvh, voi, 'name');
            singleStructDVH = obj.dvh(ix);
        end
        
        function singleStructQi = getSingleStructQi(obj, voi)
            ix = obj.getIndexOfStruct(obj.qi, voi, 'name');
            singleStructQi = obj.qi(ix);
        end
        
        function dose = get.dose(obj)
            if isempty(obj.dose)
                dose = [];
            elseif ndims(obj.dose) ~= 3
                dose = obj.get3dDoseFromLin(obj.dose, obj.subIx, obj.doseCubeDim);
            else
                dose = obj.dose;
            end
        end
        
        function set.dose(obj, v)
            obj.dose = v;
        end
        
        function doseLin = get.doseLin(obj)
            % is a derived quantity to ensure subIx like linear output of
            % dose.
            if isempty(obj.dose)
                doseLin = [];
            else
                doseLin = obj.dose(obj.subIx);
            end  
        end
        
        function dvh = get.dvh(obj)
            if isempty(obj.dvh)
                dvh = NaN;
                warning('DVH is not yet calculated.');
            else
                dvh = obj.dvh;
            end
        end
        
        function qi = get.qi(obj)
            if isempty(obj.qi)
                qi = NaN;
                warning('QI are not yet calculated.');
            else
                qi = obj.qi;
            end
        end
        
        function shift = set.shift(obj, v)
            if isempty(v)
                obj.shift = [0 0 0];
            else
                obj.shift = v;
            end
        end
        
    end % eof methods
    
    methods (Access = private)
        function obj = handleChangeOfDose(obj,src,data)
            obj.dvhQiReady = false;
            obj.resetComputedProperties;
        end
        
        function obj = resetComputedProperties(obj)
            obj.dvhQiReady = false;
            obj.dvh = [];
            obj.qi = [];
            obj.cst = [];
            obj.pln = [];
        end
        
    end
    
    methods (Static, Access = private)
        
        function ix = getIndexOfStruct(voiInStruct, voi, fieldname)
            if isnumeric(voi)
                ix = voi;
            elseif ischar(voi)
                for i = 1:numel(voiInStruct)
                    if strcmp(voiInStruct(i).(fieldname), voi)
                        ix = i;
                        break;
                    end
                end
            else
                ix = NaN;
                error('VOI needs to be either string or numeric');
            end
        end
        
        function dose3D = get3dDoseFromLin(doseLin, subIx, doseCubeDim)
            dose3D = zeros(doseCubeDim(1), doseCubeDim(2), doseCubeDim(3));
            dose3D(subIx) = doseLin;
        end
        
    end
    
    methods (Static)
        function cubes = computeDerivedCubes(physicalDose, alpha, beta, bioParam)
          RBExD = physicalDose * dij.RBE;
          %%% WRITE HERE %%%
          
        end
    end
    
end % eof classdef
