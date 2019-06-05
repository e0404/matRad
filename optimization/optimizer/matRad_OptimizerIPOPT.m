classdef matRad_OptimizerIPOPT < matRad_Optimizer
% matRad_OptimizerIPOPT implements the interface for ipopt
%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
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
        options
        wResult
        resultInfo
        env
    end
    
    properties (Access = private)
        allObjectiveFunctionValues
        axesHandle
        abortRequested
    end
    
    methods
        function obj = matRad_OptimizerIPOPT
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            %obj = createDefaultOptimizerOptions;
            %obj.Property1 = inputArg1 + inputArg2;
            obj.wResult = [];
            obj.resultInfo = [];
            obj.axesHandle = [];
            obj.allObjectiveFunctionValues = [];
            obj.abortRequested = false;
            
            %Set Default Options
            obj.options.print_level                   = 5;
            obj.options.print_user_options            = 'no';
            obj.options.print_options_documentation   = 'no';
            
            % Termination (C.2)
            obj.options.tol                           = 1e-8; % (Opt1)
            obj.options.dual_inf_tol                  = 1;    % (Opt2)
            obj.options.constr_viol_tol               = 1e-4; % (Opt3)
            obj.options.compl_inf_tol                 = 1e-4; % (Opt4), Optimal Solution Found if (Opt1),...,(Opt4) fullfiled
            
            obj.options.acceptable_iter               = 3;    % (Acc1)
            obj.options.acceptable_tol                = 1e10; % (Acc2)
            obj.options.acceptable_constr_viol_tol    = 1e10; % (Acc3)
            obj.options.acceptable_dual_inf_tol       = 1e10; % (Acc4)
            obj.options.acceptable_compl_inf_tol      = 1e10; % (Acc5)
            obj.options.acceptable_obj_change_tol     = 1e-3; % (Acc6), Solved To Acceptable Level if (Acc1),...,(Acc6) fullfiled
            
            obj.options.max_iter                      = 1000;
            obj.options.max_cpu_time                  = 3000;
            
            % Barrier Parameter (C.6)
            obj.options.mu_strategy = 'adaptive';
            
            % Line Sarch (C.8)
            %obj.options.accept_every_trial_step = 'yes';
            %obj.options.line_search_method = 'cg-penalty';
            
            % Restoration Phase (C.10)
            %obj.options.soft_resto_pderror_reduction_factor = 100;
            %obj.options.required_infeasibility_reduction    = 0.9999;
            
            % Quasi-Newton (C.13)
            obj.options.hessian_approximation         = 'limited-memory';
            obj.options.limited_memory_max_history    = 6;
            obj.options.limited_memory_initialization = 'scalar2';
            
            % determine once if Matlab or Octave
            obj.env = matRad_getEnvironment();
            
            % for derivate checking
            % obj.options.derivative_test              = 'first-order'; % none / first-order / second-order / only-second-order
            % obj.options.derivative_test_perturbation = 1e-6; % default 1e-8
            % obj.options.derivative_test_tol          = 1e-6;

        end
        
        function obj = optimize(obj,w0,optiProb,dij,cst)
            % set optimization options
            %options.radMod          = pln.radiationMode;
            %options.bioOpt          = pln.propOpt.bioOptimization;
            %options.ID              = [pln.radiationMode '_' pln.propOpt.bioOptimization];
            %options.numOfScenarios  = dij.numOfScenarios;
            
            %Set up ipopt structure
            ipoptStruct = struct;
            
            %optimizer options
            ipoptStruct.ipopt = obj.options;
            
            %variable bounds
            ipoptStruct.lb = optiProb.lowerBounds(w0);
            ipoptStruct.ub = optiProb.upperBounds(w0);
            
            %constraint bounds;
            [ipoptStruct.cl,ipoptStruct.cu] = optiProb.matRad_getConstraintBounds(cst);
            
            % set callback functions.
            
            funcs.objective         = @(x) optiProb.matRad_objectiveFunction(x,dij,cst);
            funcs.constraints       = @(x) optiProb.matRad_constraintFunctions(x,dij,cst);
            funcs.gradient          = @(x) optiProb.matRad_objectiveGradient(x,dij,cst);
            funcs.jacobian          = @(x) optiProb.matRad_constraintJacobian(x,dij,cst);
            funcs.jacobianstructure = @( ) optiProb.matRad_getJacobianStructure(w0,dij,cst);
            funcs.iterfunc          = @(iter,objective,paramter) obj.iterFunc(iter,objective,paramter,ipoptStruct.ipopt.max_iter);
            
            % Informing user to press q to terminate optimization
            fprintf('\nOptimzation initiating...\n');
            fprintf('Press q to terminate the optimization...\n');
            
            % set Callback
            
            if ~isdeployed % only if _not_ running as standalone                               
                switch obj.env
                    case 'MATLAB'
                        % get handle to Matlab command window
                        mde         = com.mathworks.mde.desk.MLDesktop.getInstance;
                        cw          = mde.getClient('Command Window');
                        xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);
                        h_cw        = handle(xCmdWndView,'CallbackProperties');
                        
                        % set Key Pressed Callback of Matlab command window
                        set(h_cw, 'KeyPressedCallback', @(h,event) obj.abortCallbackKey(h,event));
                end                
            end
            
            %ipoptStruct.options = obj.options;
            obj.abortRequested = false;
            % Run IPOPT.
            [obj.wResult, obj.resultInfo] = ipopt(w0,funcs,ipoptStruct);
            
            % unset Key Pressed Callback of Matlab command window
            if ~isdeployed && strcmp(obj.env,'MATLAB')
                set(h_cw, 'KeyPressedCallback',' ');
            end
            
            obj.abortRequested = false;
            % Empty the array of stored function values
            obj.allObjectiveFunctionValues = [];
        end
        
        function [statusmsg,statusflag] = GetStatus(obj)
            try
                switch obj.resultInfo.status
                    case 0
                        statusmsg = 'solved';
                    case 1
                        statusmsg = 'solved to acceptable level';
                    case 2
                        statusmsg = 'infeasible problem detected';
                    case 3
                        statusmsg = 'search direction too small';
                    case 4
                        statusmsg = 'diverging iterates';
                    case 5
                        statusmsg = 'user requested stop';
                    case -1
                        statusmsg = 'maximum number of iterations';
                    case -2
                        statusmsg = 'restoration phase failed';
                    case -3
                        statusmsg = 'error in step computation';
                    case -4
                        statusmsg = 'maximum CPU time exceeded';
                    case -10
                        statusmsg = 'not enough degrees of freedom';
                    case -11
                        statusmsg = 'invalid problem definition';
                    case -12
                        statusmsg = 'invalid option';
                    case -13
                        statusmsg = 'invalid number detected';
                    case -100
                        statusmsg = 'unrecoverable exception';
                    case -101
                        statusmsg = 'non-IPOPT exception thrown';
                    case -102
                        statusmsg = 'insufficient memory';
                    case -199
                        statusmsg = 'IPOPT internal error';
                    otherwise
                        statusmsg = 'IPOPT returned unknown status';
                end
                
                if obj.resultInfo.status == 0
                    statusflag = 0;
                elseif obj.resultInfo.status > 0
                    statusflag = 1;
                else
                    statusflag = -1;
                end
            catch
                statusmsg = 'No Last IPOPT Status Available!';
                statusflag = -1;
            end
        end
        
        function flag = iterFunc(obj,iter,objective,~,~)
            obj.allObjectiveFunctionValues(iter + 1) = objective;
            if strcmp(obj.env,'MATLAB')
                obj.plotFunction();
            end
            flag = ~obj.abortRequested;
        end
        
        function plotFunction(obj)
            % plot objective function output 
            if isempty(obj.axesHandle) || ~isgraphics(obj.axesHandle,'axes')
                hFig = figure('Name','Progress of IPOPT Optimization','NumberTitle','off','Color',[.5 .5 .5]);
                obj.axesHandle = axes(hFig);
                hold(obj.axesHandle,'on');
                grid(obj.axesHandle,'on');
                grid(obj.axesHandle,'minor');
                c = uicontrol;
                c.String = 'Stop';
                c.Position(1) = 5;
                c.Position(2) = 5;
                c.Callback = @obj.abortCallbackButton;                
            else
                hFig = obj.axesHandle.Parent;
            end
            
            defaultFontSize = 14;
            set(obj.axesHandle,'YScale','log');
            title(obj.axesHandle,'Progress of Optimization','LineWidth',defaultFontSize),
            xlabel(obj.axesHandle,'# iterations','Fontsize',defaultFontSize),ylabel(obj.axesHandle,'objective function value','Fontsize',defaultFontSize)
            
            % draw updated axes
            plot(obj.axesHandle,1:numel(obj.allObjectiveFunctionValues),obj.allObjectiveFunctionValues,'xb','LineWidth',1.5);
            drawnow;
            
            % ensure to bring optimization window to front also for a re-optimization
            figure(hFig);
            
        end
        
        function abortCallbackKey(obj,~,KeyEvent)
            
            % check if user pressed q
            if  get(KeyEvent,'keyCode') == 81
                
                obj.abortRequested = true;
                
            end
            
        end
        
        function abortCallbackButton(obj,~,~,~)
            obj.abortRequested = true;
        end       
        
    end
end
