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
        tol_obj {mustBeNumeric}            = 1e-20;
        tol_violation {mustBeNumeric}      = 1e-20;
        accepted_violation {mustBeNumeric} = 1e-20;
        num_reductions {mustBeInteger}  = 2;
        weighted = false;
        
    end
    
    properties (Access = private)
    M; % temp variables for the feasibility seeker
    end
    
    
    methods
        function obj = matRad_OptimizerSuperization(pln)
            %matRad_OptimizerFmincon 
            %Construct an instance of the superization optimizer
            
            obj.wResult = [];
            obj.resultInfo = struct();
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
            
            if obj.weighted
                fprintf("Getting linear inequalities and weights from cst. \n")
                [A, b, c, A_norm, weights] = obj.getvariables(dij ,cst);
            else
                fprintf("Getting linear inequalities from cst. \n")
                [A, b, c, A_norm, ~] = obj.getvariables(dij ,cst);
            end
            
            % Check if inequalities are present
            
            if isempty(A)
                error("No linear inequalities found. Please, consider using IPOPT.")
            elseif strcmp(obj.feasibility_seeker, 'AMS_sequential') && size(A, 2)>50000
                warning("sequential AMS does not work well for many voxels. Consider unsing simultaneous AMS.")
            end
            
            fprintf("Starting optimization ... \n")
            
            startSuper = tic; %Start timing the optimization

            %Set objective function and gradient of objective function.    
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
                violation = (max(0, x'*A-b')+ max(0, c'-x'*A))'./A_norm; 
                obj.resultInfo.obj_values(i+1) = f(x);
                obj.resultInfo.norm2_violations(i+1) = (1/n_A)*norm(violation);
                obj.resultInfo.max_violations(i+1) = norm(violation, 'inf');
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
                %if (i ~= 0 ... 
                %&& (obj.resultInfo.norm2_violations(i+1) < obj.accepted_violation) ...
                %&& (abs((obj.resultInfo.obj_values(i)-obj.resultInfo.obj_values(i+1))/max(1, obj.resultInfo.obj_values(i)))< obj.tol_obj))
                 %   terminated = true;
                 %   obj.resultInfo.stoppedby = "Acceptable_solution";
                 %   fprintf('Acceptable solution found. Terminating now... \n')
                %elseif (i ~= 0  ... 
                %&& (abs((obj.resultInfo.norm2_violations(i)-obj.resultInfo.norm2_violations(i+1))/max(1, obj.resultInfo.norm2_violations(i))) < obj.tol_violation) ...
                %&& (abs((obj.resultInfo.obj_values(i)-obj.resultInfo.obj_values(i+1))/max(1, obj.resultInfo.obj_values(i))< obj.tol_obj))
                %    terminated = true;
                %    obj.resultInfo.stoppedby = "No_change";
                %    fprintf('No significant change in residual or objectiv function. Terminating now... \n')
                if toc(startSuper) > obj.max_time
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
                    
                    % Deconstruct temp variables
                    clear obj.M obj.A_norm obj.weights 
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
                if obj.weighted
                    x = seeker(x, A, b, c, A_norm, weights);
                else
                    x = seeker(x, A, b, c, A_norm);
                end
            end

        end
        
        function [A, b, c, A_norm, weights] = getvariables(obj, dij,cst)
           
            b=zeros(dij.ctGrid.numOfVoxels, 1); 
            c=zeros(dij.ctGrid.numOfVoxels, 1);
            weights = zeros(dij.ctGrid.numOfVoxels, 1);
            idxs=[];
            
            for i=1:size(cst, 1)
               for j=1:numel(cst{i, 6})
                   if strcmp(cst{i, 6}{j}.name, 'Min/Max dose constraint')
                        
                        cl=cst{i, 6}{j}.parameters{1, 1};
                        cu=cst{i, 6}{j}.parameters{1, 2};
                        if obj.weighted
                            if ( strcmp(cst{i, 6}{j}.name, 'Min/Max dose constraint') && ...
                            numel(cst{i, 6}{j}.parameters) == 4 )
                                weight = cst{i, 6}{j}.parameters{1, 4};  
                            elseif ( strcmp(cst{i, 6}{j}.name, 'Min/Max dose constraint') && ...
                            numel(cst{i, 6}{j}.parameters) ~= 4 )
                                error("Either you forgot to set weights or mixed weighted and unweighted constraints")
                            end
                            
                        end
                        
                        idx = cst{i, 4}{1, 1};
                        idxs=[idxs; idx];
                        if ~isempty(idx)
                            c(idx)=cl;
                            b(idx)=cu;
                            if obj.weighted
                                weights(idx) = weight;
                            end
                       end
                   end
               end
            end
            A=dij.physicalDose{1,1}';
            
            %Only include voxels that have a constraint.
            A=A(:, idxs);
            b=b(idxs);
            c=c(idxs);
            weights = weights(idxs);
            
            %Remove all voxels that are not effected by a ray.
            A_norm = sum(A.^2, 1);
            idxs = sum(A.^2, 1)>0;
            A=A(:, idxs);
            b=b(idxs);
            c=c(idxs);
            weights = weights(idxs);
            weights = weights./max(weights);
            A_norm = A_norm(idxs)';
            
        end
        
  
        
        
        function [x] = AMS_sim(obj, x, A, b, c, A_norm, weights)
            
            % Precalculat persistent variables
            if isempty(obj.M) && nargin == 6
                obj.M = obj.lambda*A./A_norm'; 
            elseif isempty(obj.M) && nargin == 7
                weights = obj.lambda*weights;
                obj.M = A.*(weights./A_norm)';
            end
            
            A_x= x'*A;
            res_b=b-A_x';
            res_c=A_x'-c;
            x=x+(1/max(1, sum(res_b<0)+sum(res_c<0)))*(obj.M(:, res_b<0)*res_b(res_b<0)-obj.M(:, res_c<0)*res_c(res_c<0));
            x(x<0)=0;   
        end
        
        
        function [x] = AMS_sequential(obj, x, A, b, c, A_norm, weights)
            
            % Precalculat persistent variables
            if isempty(obj.M) && nargin == 6
                obj.M = obj.lambda*A./A_norm'; 
            elseif isempty(obj.M) && nargin == 7
                weights = obj.lambda*weights;
                obj.M = A.*(weights./A_norm)';
            end
                               
            for i=1:size(A, 2)
                
                A_x = x'*A(:, i);
                res_b=b(i)-A_x;
                res_c=A_x-c(i);
               
                if res_b < 0
                    x=x+obj.M(:, i)*res_b;
                elseif res_c <0
                    x=x-obj.M(:, i)*res_c;
                end
             end
            
            x(x<0) = 0;
        end
    end
        
    methods (Static)
        function available = IsAvailable()
           available = true; % if you made it to here the function exists.
        end
    end
end
