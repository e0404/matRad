classdef matRad_OptimizerFmincon < matRad_Optimizer
% matRad_OptimizerFmincon implements the interface for the fmincon optimizer 
% of the MATLAB Optiization toolbox
%    
% References
%   -
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
        options     %the optimoptions for fmincon
        wResult     %last optimization result
        resultInfo  %info struct about last results
    end
    
    methods
        function obj = matRad_OptimizerFmincon
            %matRad_OptimizerFmincon 
            %   Construct an instance of the fmincon optimizer from the Optimization Toolbox           
            matRad_cfg = MatRad_Config.instance();
            
            obj.wResult = [];
            obj.resultInfo = [];
            
            %createDefaultOptimizerOptions Constructs a set of default
            %options for the optimizer to use
            obj.options = optimoptions('fmincon',...
                'Algorithm','interior-point',...
                'Display','iter-detailed',...
                'SpecifyObjectiveGradient',true,...
                'SpecifyConstraintGradient',true,...
                'AlwaysHonorConstraints', 'bounds',...
                'MaxIterations',matRad_cfg.propOpt.defaultMaxIter,...
                'MaxFunctionEvaluations',3000,...
                'CheckGradients',false,...
                'HessianApproximation',{'lbfgs',6},...
                'UseParallel',true,...
                'Diagnostics','on',...
                'ScaleProblem',true,...
                'PlotFcn',{@optimplotfval,@optimplotx,@optimplotfunccount,@optimplotconstrviolation,@optimplotstepsize,@optimplotfirstorderopt});                    
        end
                
        function obj = optimize(obj,w0,optiProb,dij,cst)
            %optimize Carries Out the optimization
            
            % obtain lower and upper variable bounds
            lb = optiProb.lowerBounds(w0);
            ub = optiProb.upperBounds(w0);
                        
            % Informing user to press q to terminate optimization
            %fprintf('\nOptimzation initiating...\n');
            %fprintf('Press q to terminate the optimization...\n');
            
            % Run fmincon.
            [obj.wResult,fVal,exitflag,info] = fmincon(@(x) obj.fmincon_objAndGradWrapper(x,optiProb,dij,cst),...
                w0,... % Starting Point
                [],[],... % Linear Constraints we do not explicitly use
                [],[],... % Also no linear inequality constraints
                lb,ub,... % Lower and upper bounds for optimization variable
                @(x) obj.fmincon_nonlconWrapper(x,optiProb,dij,cst),...
                obj.options); % Non linear constraint structure);
            
            obj.resultInfo = info;
            obj.resultInfo.fVal = fVal;
            obj.resultInfo.exitflag = exitflag;
        end
        
        function [f, fGrad] = fmincon_objAndGradWrapper(obj,x,optiProb,dij,cst)
            f = optiProb.matRad_objectiveFunction(x,dij,cst);
            fGrad = optiProb.matRad_objectiveGradient(x,dij,cst);
        end
        
        function [c,cEq,cJacob,cEqJacob] = fmincon_nonlconWrapper(obj,x,optiProb,dij,cst)
            %Get the bounds of the constraint
            [cl,cu] = optiProb.matRad_getConstraintBounds(cst);
                    
            %Get finite bounds
            clFinIx = isfinite(cl);
            cuFinIx = isfinite(cu);
            
            % Some checks
            assert(isequal(size(cl),size(cu)));
            assert(all(cl <= cu));
            
            %For fmincon we need to separate into equalty and inequality
            %constraints
            isEqConstr = (cl == cu);
            eqIx = isEqConstr;
            ineqIx = ~isEqConstr;
            
            %Obtain all constraint functions and derivatives
            cVals = optiProb.matRad_constraintFunctions(x,dij,cst);
            cJacob = optiProb.matRad_constraintJacobian(x,dij,cst);
            
            %Subselection of equality constraints
            cEq = cVals(eqIx & clFinIx); %We can only rely on cl indices here due to the equality index
            cEqJacob = cJacob(eqIx & clFinIx,:)';
            
            %Prepare inequality constraints:
            %We need to separate upper and lower bound constraints for
            %fmincon
            cL = cl(ineqIx & clFinIx) - cVals(ineqIx & clFinIx);
            cU = cVals(ineqIx & cuFinIx) - cu(ineqIx & cuFinIx);
            cJacobL = -cJacob(ineqIx & clFinIx,:);
            cJacobU = cJacob(ineqIx & cuFinIx,:);
            
            %build the inequality jacobian
            c = [cL; cU];
            cJacob = transpose([cJacobL; cJacobU]);
        end
        
        function [statusmsg,statusflag] = GetStatus(obj)
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
            %'fmincon' is a p-code file in the optimization toolbox
            available = exist('fmincon') == 6;
        end
    end
end