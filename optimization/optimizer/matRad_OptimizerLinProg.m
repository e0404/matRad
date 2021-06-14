classdef matRad_OptimizerLinProg < matRad_Optimizer
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    %
    %   subclass of matRad_Optimizer
    %   Inherated properties
    %       options
    %       wResult
    %       resultInfo
    
    properties
        InitialWeights
        
    end
    
    methods
        function obj = matRad_OptimizerLinProg(wInit,optiProb,dij,cst)
            %CONSTRUCTOR Construct an instance of this class
            %   Call obj.optimize(wInit,optiProb,dij,cst) for optimization.
            
            obj.wResult     = [];
            obj.resultInfo  = [];
            obj.wInit = wInit
            % snack dose objectives from cst struct
            opj.
            
        end
        
        function obj = optimize();
            %OPTIMIZATION Carries out the linear optimization
            %   Detailed explanation goes here (fill out at the end.)
            
            % obtain lower and upper weight bounds
            if isfield(pln.propOpt,'lowerWeightBounds')
                lb = pln.propOpt.lowerWeightBounds*ones(size(wInit));
                ub = pln.propOpt.upperWeightBounds*ones(size(wInit));
            else
                lb = optiProb.lowerBounds(w0);
                ub = optiProb.upperBounds(w0);
            end
            % snack dose objectives from cst struct
            
            
            % optimization
            
        end
            
            
        function 
                
            
            
        end
    end
end

