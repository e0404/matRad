classdef matRad_OptimizerFmincon < matRad_Optimizer
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        options
        wResult
        resultInfo
    end
    
    methods
        function obj = matRad_OptimizerFmincon
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            %obj = createDefaultOptimizerOptions;
            %obj.Property1 = inputArg1 + inputArg2;
            obj.wResult = [];
            obj.resultInfo = [];
        end
        
        function obj  = createDefaultOptimizerOptions(obj)
            obj.options = optimoptions('fmincon',...
                'Algorithm','interior-point',...
                'Display','iter-detailed',...
                'SpecifyObjectiveGradient',true,...
                'SpecifyConstraintGradient',true,...
                'AlwaysHonorConstraints', 'bounds',...
                'MaxIterations',500,...
                'MaxFunctionEvaluations',3000,...
                'CheckGradients',false,...
                'HessianApproximation',{'lbfgs',6},...
                'UseParallel',true,...
                'Diagnostics','on',...
                'ScaleProblem',true,...
                'PlotFcn',{@optimplotfval,@optimplotx,@optimplotfunccount,@optimplotconstrviolation,@optimplotstepsize,@optimplotfirstorderopt});
            
            %obj.options
        end
        
        function obj = optimize(obj,w0,optiProb,dij,cst)
            
            % set callback functions.
            lb = optiProb.lowerBounds(w0);
            ub = optiProb.upperBounds(w0);
                        
            % Informing user to press q to terminate optimization
            fprintf('\nOptimzation initiating...\n');
            fprintf('Press q to terminate the optimization...\n');
            
            % Run fmincon.
            [obj.wResult,~,obj.resultInfo] = fmincon(@(x) obj.fmincon_objAndGradWrapper(x,optiProb,dij,cst),...
                w0,... % Starting Point
                [],[],... % Linear Constraints we do not explicitly use
                [],[],... % Also no linear inequality constraints
                [],[],... lb,ub,... % Lower and upper bounds for optimization variable
                @(x) obj.fmincon_nonlconWrapper(x,optiProb,dij,cst),...
                obj.options); % Non linear constraint structure);
        end
        
        function [f, fGrad] = fmincon_objAndGradWrapper(obj,x,optiProb,dij,cst)
            f = optiProb.matRad_objectiveFunction(x,dij,cst);
            fGrad = optiProb.matRad_objectiveGradient(x,dij,cst);
        end
        
        function [c,cEq,cJacob,cEqJacob] = fmincon_nonlconWrapper(obj,x,optiProb,dij,cst)
            [cl,cu] = optiProb.matRad_getConstraintBounds(cst);
            
            % Some checks
            assert(isequal(size(cl),size(cu)));
            assert(all(cl <= cu));
            
            %For fmincon we need to separate into equalty and inequality
            %constraints
            isEqConstr = cl == cu;
            eqIx = find(isEqConstr);
            ineqIx = find(~isEqConstr);
            
            cVals = optiProb.matRad_constraintFunctions(x,dij,cst);
            cJacob = optiProb.matRad_constraintJacobian(x,dij,cst);
            
            cL = cl(ineqIx) - cVals(ineqIx);
            cU = cVals(ineqIx) - cu(ineqIx);
            cJacobL = -cJacob(ineqIx,:);
            cJacobU = cJacob(ineqIx,:);
            
            c = [cL; cU];
            cJacob = transpose([cJacobL; cJacobU]);
            
            cEq = cVals(eqIx);
            cEqJacob = cJacob(eqIx,:)';
        end
    end
end