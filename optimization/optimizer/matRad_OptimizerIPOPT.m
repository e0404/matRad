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
% References
%   -
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
        plotHandle
        abortRequested
        plotFailed
    end
    
    methods
        function obj = matRad_OptimizerIPOPT
            %matRad_OptimizerIPOPT
            %   Construct an instance of the IPOPT optimizer (mex
            %   interface)
            
            matRad_cfg = MatRad_Config.instance();
            
            obj.wResult = [];
            obj.resultInfo = [];
            obj.axesHandle = [];
            obj.allObjectiveFunctionValues = [];
            obj.abortRequested = false;
            
            %Set Default Options
            if matRad_cfg.logLevel <= 1
                lvl = 0;
            elseif matRad_cfg.logLevel <= 2
                lvl = 2;
            elseif matRad_cfg.logLevel <= 3
                lvl = 5;
            else 
                %There seems to be a problem with higher log levels in
                %IPOPT!
                lvl = 5;
            end
                
            obj.options.print_level                   = lvl;
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
            
            obj.options.max_iter                      = matRad_cfg.propOpt.defaultMaxIter;
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
            
            if ~matRad_checkMexFileExists('ipopt')
                matRad_cfg.dispError('IPOPT mex interface not available for %s!',obj.env);
            end

        end
        
        function obj = optimize(obj,w0,optiProb,dij,cst)
            matRad_cfg = MatRad_Config.instance();
            
            % set optimization options            
            
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
            matRad_cfg.dispInfo('\nOptimzation initiating...\n');
            
            % set Callback
            qCallbackSet = false;
            if ~isdeployed % only if _not_ running as standalone                               
                switch obj.env
                    case 'MATLAB'
                        try
                            % get handle to Matlab command window
                            mde         = com.mathworks.mde.desk.MLDesktop.getInstance;
                            cw          = mde.getClient('Command Window');
                            xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);
                            h_cw        = handle(xCmdWndView,'CallbackProperties');
                        
                            % set Key Pressed Callback of Matlab command window
                            set(h_cw, 'KeyPressedCallback', @(h,event) obj.abortCallbackKey(h,event));
                            fprintf('Press q to terminate the optimization...\n');
                            qCallbackSet = true;
                        catch
                            matRad_cfg.dispInfo('Manual termination with q not possible due to failing callback setup.\n');
                        end
                end                
            end
            
            %ipoptStruct.options = obj.options;
            obj.abortRequested = false;
            obj.plotFailed = false;
            
            % Run IPOPT.
            try
                [obj.wResult, obj.resultInfo] = ipopt(w0,funcs,ipoptStruct);
            catch ME
                errorString = [ME.message '\nThis error was thrown by the MEX-interface of IPOPT.\nMex interfaces can raise compatability issues which may be resolved by compiling them by hand directly on your particular system.'];
                matRad_cfg.dispError(errorString);
            end
            
            % unset Key Pressed Callback of Matlab command window
            if qCallbackSet
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
            %We don't want the optimization to crash because of drawing
            %errors
            if ~obj.plotFailed
                try            
                    obj.plotFunction();
                catch ME
                    matRad_cfg = MatRad_Config.instance();
                    %Put a warning at iteration 1 that plotting failed
                    matRad_cfg.dispWarning('Objective Function plotting failed and thus disabled. Message:\n%s',ME.message);
                    obj.plotFailed = true;
                end                
            end
            flag = ~obj.abortRequested;
        end
        
        function plotFunction(obj)
            % plot objective function output
            y = obj.allObjectiveFunctionValues;
            x = 1:numel(y);
            
            if isempty(obj.axesHandle) || ~isgraphics(obj.axesHandle,'axes')
                %Create new Fiure and store axes handle
                hFig = figure('Name','Progress of IPOPT Optimization','NumberTitle','off','Color',[.5 .5 .5]);
                hAx = axes(hFig);
                hold(hAx,'on');
                grid(hAx,'on');
                grid(hAx,'minor');
                set(hAx,'YScale','log');
                 
                %Add a Stop button with callback to change abort flag
                c = uicontrol;
                cPos = get(c,'Position');
                cPos(1) = 5;
                cPos(2) = 5;
                set(c,  'String','Stop',...
                        'Position',cPos,...
                        'Callback',@(~,~) abortCallbackButton(obj));                
                
                %Set up the axes scaling & labels
                defaultFontSize = 14;
                set(hAx,'YScale','log');
                title(hAx,'Progress of Optimization','LineWidth',defaultFontSize);
                xlabel(hAx,'# iterations','Fontsize',defaultFontSize),ylabel(hAx,'objective function value','Fontsize',defaultFontSize);
                
                %Create plot handle and link to data for faster update
                hPlot = plot(hAx,x,y,'xb','LineWidth',1.5,'XDataSource','x','YDataSource','y');
                obj.plotHandle = hPlot;
                obj.axesHandle = hAx;
                                
            else %Figure already exists, retreive from axes handle
                hFig = get(obj.axesHandle,'Parent');
                hAx = obj.axesHandle;
                hPlot = obj.plotHandle;
            end
            
            % draw updated axes by refreshing data of the plot handle (which is linked to y and y) 
            % in the caller workspace. Octave needs and works on figure handles, which
            % is substantially (factor 10) slower, thus we check explicitly
            switch obj.env
                case 'OCTAVE'
                    refreshdata(hFig,'caller');
                otherwise
                    refreshdata(hPlot,'caller');
            end
            drawnow;
            
            % ensure to bring optimization window to front
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
    
    methods (Static)
        function available = IsAvailable()
            available = matRad_checkMexFileExists('ipopt');                   
        end
    end
end
