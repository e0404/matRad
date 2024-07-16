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
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
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
            
            % Create default optimizer options
            obj.options = optimoptions('simulannealbnd', ...
                'InitialTemperature', 0.7, ...
                'TemperatureFcn',@temperatureboltz, ... 
                'Display', 'iter', ...
                'MaxIterations', matRad_cfg.defaults.propOpt.maxIter, ...
                'MaxFunctionEvaluations', 120000, ...
                'PlotFcn', {@saplotbestf,@saplotbestx, @saplotf,@saplotx,@saplotstopping,@saplottemperature});
        end
                
        function obj = optimize(obj, w0, optiProb, dij, cst)
            % optimize Carries out the optimization
            
            % Obtain lower and upper variable bounds
            lb = optiProb.lowerBounds(w0);
            ub = optiProb.upperBounds(w0);
                        
            % Informing user to press q to terminate optimization
            fprintf('\nOptimization initiating...\n');
            fprintf('Press q to terminate the optimization...\n');

            matRad_cfg = MatRad_Config.instance();
            if matRad_cfg.isMatlab && str2double(matRad_cfg.envVersion) <= 9.13 && strcmp(obj.options.Diagnostics, 'on')
                matRad_cfg.dispWarning('Diagnostics in simulannealbnd will be turned off due to a bug when using lbfgs with specified number of histories!');
                obj.options.Diagnostics = 'off';
            end
                
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
