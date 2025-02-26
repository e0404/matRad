classdef matRad_OptimizerSimulannealbnd < matRad_Optimizer
    % matRad_OptimizerSimulannealbnd implements the interface for the Simulated Annealing optimizer 
    % of the MATLAB Global Optimization toolbox
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
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        options     % the optimoptions for simulannealbnd
    end

    properties (SetAccess = protected)
        wResult     % last optimization result
        resultInfo  % info struct about last results
    end
    
    methods
        function obj = matRad_OptimizerSimulannealbnd
            % matRad_OptimizerSimulannealbnd Constructor
            matRad_cfg = MatRad_Config.instance();

            if ~matRad_OptimizerSimulannealbnd.IsAvailable()
                matRad_cfg.dipsError('matRad_OptimizerSimulannealbnd cannot be constructed as simulannealbnd is not available!');
            end
            
            obj.wResult = [];
            obj.resultInfo = [];

            if matRad_cfg.logLevel >= 4
                optDisplay = 'diagnose';
            elseif matRad_cfg.logLevel == 3
                optDisplay = 'iter';
            elseif matRad_cfg.logLevel == 2
                optDisplay = 'final';
            else
                optDisplay = 'off';
            end

            if matRad_cfg.disableGUI
                pltFcns = {[]};                
            else
                pltFcns = {@saplotbestf,@saplotbestx, @saplotf,@saplotx,@saplotstopping,@saplottemperature};
            end
            
            % Create default optimizer options
            obj.options = optimoptions('simulannealbnd', ...
                'InitialTemperature', 0.7, ...
                'TemperatureFcn',@temperatureboltz, ... 
                'Display', optDisplay, ...
                'MaxIterations', matRad_cfg.defaults.propOpt.maxIter, ...
                'MaxFunctionEvaluations', 120000, ...
                'PlotFcn',pltFcns);
        end
                
        function obj = optimize(obj, w0, optiProb, dij, cst)
            matRad_cfg = MatRad_Config.instance();
            % optimize Carries out the optimization
            
            % Obtain lower and upper variable bounds
            lb = optiProb.lowerBounds(w0);
            ub = optiProb.upperBounds(w0);
                        
            % Informing user to press q to terminate optimization
            matRad_cfg.dispInfo('Optimization initiating...\n');
            matRad_cfg.dispInfo('Press q to terminate the optimization...\n');
                
            % Define the objective function
            objectiveFunction = @(x) optiProb.matRad_objectiveFunction(x, dij, cst);
            
            % Run simulated annealing optimization
            [obj.wResult, fVal, exitflag, info] = simulannealbnd(objectiveFunction, w0, lb, ub, obj.options);
            
            obj.resultInfo = info;
            obj.resultInfo.fVal = fVal;
            obj.resultInfo.exitflag = exitflag;
        end

        function [statusmsg, statusflag] = GetStatus(obj)
            try 
                statusmsg = obj.resultInfo.message;
                if obj.resultInfo.exitflag == 0
                    statusflag = 0;
                elseif obj.resultInfo.exitflag > 0
                    statusflag = 1;
                else 
                    statusflag = -1;
                end
            catch
                statusmsg = 'No Last Optimizer Status Available!';
                statusflag = -1;
            end
        end
    end
    
    methods (Static)    
        function available = IsAvailable()
            available = exist('simulannealbnd', 'file') ~= 0;
        end
    end
end
