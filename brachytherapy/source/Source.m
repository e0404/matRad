classdef Source < matlab.mixin.Copyable
    %SOURCE Summary of this class goes here
    %   Detailed explanation goes here
    
    %      properties
    %          distance
    %          angle
    %          dictionary
    %      end
    
    properties (SetAccess = private)
        
        % lambda(Lambda): Dose-rate constant in water
        lambda
        
        SourceIsotopeName
        SourceIsotopeHalfLife
        
        % ActiveSourceLength(L): Active length of the source
        ActiveSourceLength
        ActiveSourceDiameter
        SourceLength
        SourceDiameter
        
        SourceStrengthImplanted
        SourceStrengthUnit
        SourceStrengthConversionFactor %[U/mCi]
        
        
        
        %%% Radial Dose Function
        % RadialDoseDistance: Vector of radial distance
        % RadialDoseValue: Vector of tabulated values of Radial Dose Function
        RadialDoseDistance
        RadialDoseValue
        
        %%% 1D anisotropy function
        % AnisotropyFactorRadialDistance: Vector of radial distance
        % AnisotropyFactorValue: Vector of tabulated values of the 1D anisotropy function
        AnisotropyFactorRadialDistance
        AnisotropyFactorValue
        
        %%% 2D anisotropy function
        % AnisotropyRadialDistances: Vector of radial distance
        % AnisotropyPolarAngles: Vector of polar angles
        % AnisotropyFunctionValue: Anisotropy matric on 1D array form
        % AnisotropyMatrix: Matrix of tabulated values of the 2D anisotropy function
        AnisotropyRadialDistances
        AnisotropyPolarAngles
        AnisotropyFunctionValue
        AnisotropyMatrix
        
        % Treatment Type indicates whether LDR or HDR is used for
        % calculation
        SourceType
        
        %% for qa
        QA
        
    end
    
    methods
        %construction function
        function obj = Source(input,var1,var2)
            
            if isa(input,'Source')
                sourceStruct = [];
                obj = copyElement(input);
                switch nargin
                    case 2
                        obj.SourceStrengthImplanted = var1;
                    case 3
                        obj.setSourceStrength(var1,var2)
                end
            elseif isstruct(input)
                sourceStruct = input;
            elseif ischar(input)
                global d;
                sourceStruct = obj.loadSource(input);
                obj.SourceIsotopeHalfLife = obj.SourceIsotopeHalfLife*d;
            else
                error('Input is not valid')
            end
            
            if ~isempty(sourceStruct)
                fNames = fieldnames(obj);
                for i=1:length(fNames)
                    if isfield(sourceStruct,fNames{i})
                        obj.(fNames{i}) = sourceStruct.(fNames{i});
                    else
                        obj.(fNames{i}) = [];
                    end
                end
            end
            
            if isempty(obj.AnisotropyMatrix)
                obj.AnisotropyPolarAngles = obj.AnisotropyPolarAngles;
                columns = length(obj.AnisotropyRadialDistances);
                rows    = length(obj.AnisotropyPolarAngles);
                obj.AnisotropyMatrix = zeros(rows,columns);
                for rr=1:rows
                    obj.AnisotropyMatrix(rr,:) = obj.AnisotropyFunctionValue(((rr-1)*columns+1):rr*columns);
                end
            end
            
            if isempty(obj.AnisotropyFunctionValue)
                obj.AnisotropyFunctionValue = reshape(obj.AnisotropyMatrix',[],1);
            end
            
            
        end
        
        
        DoseRate = getDoseRate1D(obj,r)
        DoseRate = getDoseRate2D(obj,r,theta)
        
        Dose = getDose1D(obj,r,time)
        Dose = getDose2D(obj,r,theta,time)
        
        function setSourceStrength(obj,strength,unit)
            switch nargin
                case 2
                    obj.SourceStrengthImplanted = strength;
                case 3
                    if strcmp(unit,'U')
                        obj.SourceStrengthImplanted = strength;
                    elseif strcmp(unit,'mCi')
                        obj.SourceStrengthImplanted = strength*obj.SourceStrengthConversionFactor;
                    elseif strcmp(unit,'Ci')
                        obj.SourceStrengthImplanted = strength*obj.SourceStrengthConversionFactor*1e3;
                    else
                        error('Unit not recognized')
                    end
            end
        end
        
        function strength = getSourceStrength(obj,unit)
            if strcmp(unit,'U')
                strength = obj.SourceStrengthImplanted;
            elseif strcmp(unit,'mCi')
                strength = obj.SourceStrengthImplanted/obj.SourceStrengthConversionFactor;
            end
        end
        
        function setSourceLength(obj,length)
            obj.SourceLength = length;
        end
        
        function setSourceDiameter(obj,diameter)
            obj.SourceDiameter = diameter;
        end
        
    end
    
    
    methods (Access = private)
        
        Phi_an   = get1DAnisotropyFunction(obj,r);
        Phi_an   = get2DAnisotropyFunction(obj,r,theta);
        
        gL       = getRadialDoseFunction(obj,r);
        GL       = getTransverseGeometryFunction(obj,r);
        
        DoseRate = getDoseRateAtPoint(obj,r,theta,time)
        doseRateTable = getAbsoluteDoseRateTable(obj,x,y,direction)
        plot1DDoseRate(obj,r)
        
    end
    
    methods (Access = protected)
        function cpObj = copyElement(obj)
            cpObj = copyElement@matlab.mixin.Copyable(obj);
        end
    end
    
    methods (Static)
        data = loadSource(input)
    end
    
    
    methods
        doseStruct = getDoseStruct(obj,time,distance,theta);
        
        function time = getTime(obj,distance,dose)
            time = dose/obj.getDoseRate1D(distance);
        end
    end
    
end

