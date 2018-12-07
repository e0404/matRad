classdef (Abstract) matRad_Optimizer < handle
    %OPTIMIZER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract)
        options %options struct
        wResult
        resultInfo
    end
    
   
    %These should be abstract methods, however Octave can't parse them. As soon 
    %as Octave is able to do this, they should be made abstract again 
    methods %(Abstract)        
        function obj = optimize(obj,w0,optiProb,dij,cst)
          error('Function needs to be implemented!');
        end
        
        function [msg,statusflag] = GetStatus(obj)
          error('Function needs to be implemented!');
        end
    end
end

