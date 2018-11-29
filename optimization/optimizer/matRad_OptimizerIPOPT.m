classdef matRad_OptimizerIPOPT < matRad_Optimizer
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        options
        wResult
        resultInfo
    end
    
    methods
        function obj = matRad_OptimizerIPOPT
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            %obj = createDefaultOptimizerOptions;
            %obj.Property1 = inputArg1 + inputArg2;
            obj.wResult = [];
            obj.resultInfo = [];
        end
        
        function obj  = createDefaultOptimizerOptions(obj)
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
            
            obj.options.max_iter                      = 500;
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
            
            %constraint bounds
            %ipoptStruct.cl = optiProb.getLowerConstraintBounds();
            %ipoptStruct.cu = optiProb.getUpperConstraintBounds();
            [ipoptStruct.cl,ipoptStruct.cu] = optiProb.matRad_getConstraintBounds(cst);
            
            % set callback functions.
            %[options.cl,options.cu] = matRad_getConstBoundsWrapper(cst,options);
            funcs.objective         = @(x) optiProb.matRad_objectiveFunction(x,dij,cst);
            funcs.constraints       = @(x) optiProb.matRad_constraintFunctions(x,dij,cst);
            funcs.gradient          = @(x) optiProb.matRad_objectiveGradient(x,dij,cst);
            funcs.jacobian          = @(x) optiProb.matRad_constraintJacobian(x,dij,cst);
            funcs.jacobianstructure = @( ) matRad_getJacobStruct(w,dij,cst);
            
            % Informing user to press q to terminate optimization
            fprintf('\nOptimzation initiating...\n');
            fprintf('Press q to terminate the optimization...\n');
            
            %ipoptStruct.options = obj.options;
            
            % Run IPOPT.
            [obj.wResult, obj.resultInfo] = ipopt(w0,funcs,ipoptStruct);
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
    end
    
    
end
