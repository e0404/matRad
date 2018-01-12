classdef MatRadScenario < handle

    properties (Access = private)
        dvhQiReady = false
        % during calcQiDVH
        cst
        pln
    end
    
    properties (SetAccess = private)
        shift;
        relRangeShift;
        absRangeShift;
        ctShiftIdentifier;
        weight;
        radiationQuantity;
        
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
    end
    
    methods
        function obj = MatRadScenario(dose, weight, radiationQuantity, ctShiftIdentifier, shift, ...
                relRangeShift, absRangeShift)
            
            obj.dose = dose;
            obj.weight = weight;
            obj.radiationQuantity = radiationQuantity;
            
            obj.ctShiftIdentifier = ctShiftIdentifier;
            obj.shift = shift;
            obj.relRangeShift = relRangeShift;
            obj.absRangeShift = absRangeShift;
            
            addlistener(obj,'dose','PostSet',@obj.handleChangeOfDose);
        end % eof constructor
        
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
        
        function plotDoseSlice(obj, ax, ct, plane, slice)
            matRad_plotSliceWrapper(ax,ct,obj.cst,1,obj.dose,plane,slice)
        end
        
        function singleStructDVH = getSingleStructDVH(obj, voi)
            ix = obj.getIndexOfStruct(obj.dvh, voi, 'name');
            singleStructDVH = obj.dvh(ix);
        end
        
        function singleStructQi = getSingleStructQi(obj, voi)
            ix = obj.getIndexOfStruct(obj.qi, voi, 'name');
            singleStructQi = obj.qi(ix);
        end
        
        function dvh = get.dvh(obj)
            if isempty(obj.dvh)
                error('DVH is not yet calculated.');
            else
                dvh = obj.dvh;
            end
        end
        
        function qi = get.qi(obj)
            if isempty(obj.qi)
                error('QI are not yet calculated.');
            else
                qi = obj.qi;
            end
        end
        
        
    end % eof methods
    
    methods (Access = private)
        function obj = handleChangeOfDose(obj,src,data)
            obj.dvhQiReady = false;
            obj.resetComputedProperties;
            fprintf('Changed.');
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
                    if voiInStruct(i).(fieldname) == voi
                        ix = i;
                        break;
                    end
                end
            else
                ix = NaN;
                error('VOI needs to be either string or numeric');
            end
        end
        
    end
    
end % eof classdef
