classdef matRad_OptimizerSuperization < matRad_Optimizer
    properties
        options     
        wResult     %last optimization result
        resultInfo  %info struct about last results
        
        %Set default options
        max_iter {mustBeInteger} = 500;
        max_time {mustBeNumeric, mustBePositive} = 3000; %second
        alpha {mustBeNumeric, mustBePositive, mustBeLessThan(alpha, 1)} = 0.99;
        lambda {mustBeNumeric, mustBePositive, mustBeLessThan(lambda, 200)} = 1.9;
        feasibility_seeker {mustBeMember(feasibility_seeker, {'AMS_sim', 'AMS_sequential'})} = 'AMS_sim';
        tol_obj {mustBeNumeric, mustBePositive}            = 1e-20;
        tol_violation {mustBeNumeric, mustBePositive}      = 1e-20;
        accepted_violation {mustBeNumeric, mustBePositive} = 1e-20;
        num_reductions {mustBeInteger}  = 2;
        
    end
    
    properties (Access = private)
    A_norm; 
    M; % temp variable for the feasibility seeker
    end
    
    
    methods
        function obj = matRad_OptimizerSuperization(pln)
            %matRad_OptimizerFmincon 
            %   Construct an instance of the superization optimizer
            
            obj.wResult = [];
            obj.resultInfo = struct();
            obj.A_norm = [];
            obj.M = [];
            
            % Overwrite default values with values provided in pln.propOpt,
            % IF pln is given as input.
            if nargin > 0
                  if isfield(pln, 'propOpt')
                        fields = fieldnames(pln.propOpt);
                        props = properties(obj);
                        for field=fields'
                            for i=1:numel(props)
                                if strcmp(props{i}, field{1})
                                	obj.(field{1})= pln.propOpt.(props{i});
                                end
                            end
                        end
                  end
            end
            
            % Init performance measures
            obj.resultInfo.obj_values = zeros(obj.max_iter+1,1);
            obj.resultInfo.norm2_violations = zeros(obj.max_iter+1, 1);
            obj.resultInfo.max_violations = zeros(obj.max_iter+1, 1);
            obj.resultInfo.betas = zeros(obj.max_iter+1, 1);
            
         end
                
        function obj = optimize(obj, x_0,  optiProb, dij, cst)
            global A b c f
            
            [A, b, c] = obj.getlinearinequalities(dij ,cst);
            
            % Check if inequalities are present
            
            if isempty(A)
                error("No linear inequalities found. Please, consider using IPOPT.")
            elseif strcmp(obj.feasibility_seeker, 'AMS_sequential') && size(A, 2)>50000
                warning("sequential AMS does not work well for many voxels. Consider unsing simultaneous AMS.")
            end
            
            
            fprintf("Starting optimization ... \n")
            
            startSuper = tic; %Start timing the optimization

            %Set objectiv function and gradient of objectiv function.    
            f = @(x) (optiProb.matRad_objectiveFunction(x,dij,cst));
            fGrad = @(x) (optiProb.matRad_objectiveGradient(x,dij,cst));
            
            % Choose feasibility seeker
            if strcmp(obj.feasibility_seeker, 'AMS_sim')
                seeker = @obj.AMS_sim;
            elseif strcmp(obj.feasibility_seeker, 'AMS_sequential')
                seeker = @obj.AMS_sequential;
            else
                error('The choosen feasibility seeker (%s) is not implemented.', ...
                    obj.feasibility_seeker)
            end

            % Superizaztion algorithm
            
            x=x_0;
            l=-1;
            n_A = size(A, 2); % number of voxels with constraints
            terminated = false;
            beta = 1;

            for i=0:obj.max_iter+1
                
                %Calculation of performance measures
                res = obj.residual(x, A, b, c);
                obj.resultInfo.obj_values(i+1) = f(x);
                obj.resultInfo.norm2_violations(i+1) = (1/n_A)*norm(res);
                obj.resultInfo.max_violations(i+1) = norm(res, 'inf');
                obj.resultInfo.betas(i+1) = beta;
                
                if mod(i ,10)==0
                    fprintf('iter \t objective \t ||residual|| \t max(residual) \t beta\n');
                end
                
                fprintf(sprintf('%i \t %0.4f \t %0.6f \t %0.4f \t %0.4f\n', ...
                [i, obj.resultInfo.obj_values(i+1), ...
                 obj.resultInfo.norm2_violations(i+1), ...
                 obj.resultInfo.max_violations(i+1), ...
                 obj.resultInfo.betas(i+1)]));

                %Stopping rules
                if (i ~= 0 ... 
                && (obj.resultInfo.norm2_violations(i+1) < obj.accepted_violation) ...
                && ((obj.resultInfo.obj_values(i)-obj.resultInfo.obj_values(i+1))/max(1, obj.resultInfo.obj_values(i))< obj.tol_obj))
                    terminated = true;
                    obj.resultInfo.stoppedby = "Acceptable_solution";
                    fprintf('Acceptable solution found. Terminating now... \n')
                elseif (i ~= 0  ... 
                && ((obj.resultInfo.norm2_violations(i)-obj.resultInfo.norm2_violations(i+1))/max(1, obj.resultInfo.norm2_violations(i)) < obj.tol_violation) ...
                && ((obj.resultInfo.obj_values(i)-obj.resultInfo.obj_values(i+1))/max(1, obj.resultInfo.obj_values(i))< obj.tol_obj))
                    terminated = true;
                    obj.resultInfo.stoppedby = "No_change";
                    fprintf('No significant change in residual or objectiv function. Terminating now... \n')
                elseif toc(startSuper) > obj.max_time
                    terminated = true;
                    obj.resultInfo.stoppedby = "Time_exceeded";
                    fprintf('Time exceeded. Terminating now...\n')
                elseif i == obj.max_iter
                    terminated = true;
                    obj.resultInfo.stoppedby = "Max_iter";
                    fprintf('Max number of iterations exceeded. Terminating now...\n')
                end
               
                % Stopping and returns
                if terminated
                    obj.resultInfo.obj_values =  obj.resultInfo.obj_values(1:i+1);
                    obj.resultInfo.norm2_violations = obj.resultInfo.norm2_violations(1:i+1);
                    obj.resultInfo.max_violations = obj.resultInfo.max_violations(1:i+1);
                    obj.resultInfo.betas = obj.resultInfo.betas(1:i+1);
                    obj.resultInfo.time = toc(startSuper);
                    
                    obj.wResult = x;
                    
                    % Deconstruct persistent variables
                    clear obj.m obj.A_norm
                    return
                end
                
                n=0;
                while n<obj.num_reductions
                   
                    v=fGrad(x);
                                        
                    loop = true;
                    while loop
                        l=l+1;
                        beta=(obj.alpha)^l;
                        z=x-beta*v;
                        if f(z) <= obj.resultInfo.obj_values(i+1)
                            n=n+1;
                            x=z;
                            loop=false;
                        end
                    end
                end
                
                x = seeker(x, A, b, c);
            end

        end
        
        function [A, b, c] = getlinearinequalities(obj, dij,cst)
            fprintf("Obtaining linear inequalities from cst \n")
            b=zeros(dij.ctGrid.numOfVoxels, 1); 
            c=zeros(dij.ctGrid.numOfVoxels, 1);
            idxs=[];
            
            for i=1:size(cst, 1)
               for j=1:numel(cst{i, 6})
                   if strcmp(cst{i, 6}{j}.name, 'Min/Max dose constraint')
                        cl=cst{i, 6}{j}.parameters{1, 1};
                        cu=cst{i, 6}{j}.parameters{1, 2};
                        idx = cst{i, 4}{1, 1};
                        idxs=[idxs; idx];
                        if ~isempty(idx)
                            c(idx)=cl;
                            b(idx)=cu;
                       end
                   end
               end
            end
            A=dij.physicalDose{1,1}';
            
            %Only include voxels that have a constraint.
            A=A(:, idxs);
            b=b(idxs);
            c=c(idxs);
            
            %Remove all voxels that are not effected by a ray.
            idxs = sum(A.^2, 1)>0;
            A=A(:, idxs);
            b=b(idxs);
            c=c(idxs);
            if isempty(b)
                error("No Min/Max dose constraints found. Please use IPOPT!")
            end
        end
        
        function [x] = AMS_sim(obj, x, A, b, c)
            
            L = obj.lambda; 
            
            if isempty(obj.M)  % Set M if not already set
                obj.M = A*obj.A_norm; 
            end
            
            A_x= x'*A;
            res_b=b-A_x';
            res_c=A_x'-c;
            x=x+(L/sum(res_b<0))*obj.M(:, res_b<0)*res_b(res_b<0)-(L/sum(res_c<0))*obj.M(:, res_c<0)*res_c(res_c<0);
            x(x<0)=0;   
            
        end
        
        
        function [x] = AMS_sequential(obj, x, A, b, c)
            
            L = obj.lambda; 
            
            % Precalculat persistent variables
            if isempty(obj.M)
                obj.M = A*obj.A_norm; 
            end
                                   
            for i=1:n
                
                A_x = x'*A(:, i);
                res_b=b(i)-A_x;
                res_c=A_x-c(i);
                
                if res_b < 0
                    x=x+L*(obj.M(:, i)*res_b);
                elseif res_c <0
                    x=x-L*(obj.M(:, i)*res_c);
                end
                
            end
            
            x(x<0) = 0;
        end
        
        
        function [res] = residual(obj, x, A, b, c)
            
            if isempty(obj.A_norm) % Set M if not already set
                obj.A_norm = diag(sum(A.^2, 1).^-1);
            end  
            res=max(0, obj.A_norm*(x'*A-b')') + max(0, obj.A_norm*(c'-x'*A)');
        end
    end
        
    methods (Static)
        function available = IsAvailable()
           available = true; % if you made it to here the function exists.
        end
    end
end