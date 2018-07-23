classdef matRad_OptimizerIPOPT < matRad_Optimizer
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        options
    end
    
    methods
        function obj = matRad_OptimizerIPOPT
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            %obj = createDefaultOptimizerOptions;
            %obj.Property1 = inputArg1 + inputArg2;
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
        
        function obj = optimize(obj,optiProb)
            % set optimization options
            %options.radMod          = pln.radiationMode;
            %options.bioOpt          = pln.propOpt.bioOptimization;
            %options.ID              = [pln.radiationMode '_' pln.propOpt.bioOptimization];
            %options.numOfScenarios  = dij.numOfScenarios;
            
            %Set up ipopt structure
            ipoptStruct = struct;
            
            %optimizer options
            ipoptStruct.ipopt.options = obj.options;
            
            %variable bounds
            ipoptStruct.lb = optiProb.lowerBounds();
            ipoptStruct.ub = optiProb.upperBounds();
            
            %constraint bounds
            %ipoptStruct.cl = optiProb.getLowerConstraintBounds();
            %ipoptStruct.cu = optiProb.getUpperConstraintBounds();
            [ipoptStruct.cl,ipopstStruct.cu] = optiProb.matRad_getConstraintBounds(optiProb.cst);
            
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
            
            ipoptStruct.options = obj.options;
            
            % Run IPOPT.
            [wOpt, info]            = ipopt(wInit,funcs,ipoptStruct);
        end
    end
end
