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
        control_sequence {mustBeMember(control_sequence, {'sequential', 'random', 'weight', 'weight_inv'})} = 'random';
        weight_decay = 0.99;
        temp_weight_decay = 1;
        
        tol_obj {mustBeNumeric}            = 1e-6;
        tol_violation {mustBeNumeric}      = 1e-6;
        tol_max_violation {mustBeNumeric}   = 1e-3;
        
        accepted_tol_change {mustBeNumeric} = 1e-3;
        accepted_violation {mustBeNumeric} = 1e-5;
        accepted_max_violation {mustBeNumeric} = 1e-2;
        accepted_iter = 3;
        
        num_reductions {mustBeInteger}  = 2;
        weighted = false;
        
        ignoreObjective = false;
    end
    
    properties (SetAccess = private)
        allObjectiveFunctionValues;
        allConstraintViolations;
        allOptVars;
        timeInit;
        timeStart;
        timeIter;
        
        abortRequested = false;
    end
    
    properties (Access = private)
        M; % temp variables for the feasibility seeker
        axesHandles;
        plotHandles;
        figHandle;
    end
    
    
    methods
        function obj = matRad_OptimizerSuperization(pln)
            %matRad_OptimizerFmincon
            %Construct an instance of the superization optimizer
            
            obj.wResult = [];
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
        end
        
        function obj = optimize(obj, x_0,  optiProb, dij, cst)
            
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo('Starting superiorization ... \n');
            
            obj.resultInfo = struct();
            
            obj.timeInit = tic;
            
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
            elseif strcmp(obj.feasibility_seeker, 'AMS_sequential') && size(A, 2)> 50000
                warning("sequential AMS does not work well for many voxels. Consider unsing simultaneous AMS.")
            end
            
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
            
            obj.timeStart = tic; %Start timing the optimization
            % Superizaztion algorithm
            
            x=x_0;
            l=-1;
            n_A = size(A, 2); % number of voxels with constraints
            terminated = false;
            beta = 1;
            
            %Initial values
            violation = (max(0, x'*A-b')+ max(0, c'-x'*A))'./A_norm;
            normViolation = (1/n_A)*norm(violation);
            maxViolation = norm(violation,'inf');
            
            fVal = f(x);
            obj.allObjectiveFunctionValues(1) = fVal;
            obj.allConstraintViolations(1) = maxViolation;
            obj.timeIter(1) = toc(obj.timeStart);
            
            obj.allOptVars.obj_values(1) = fVal;
            obj.allOptVars.norm2_violations(1) = normViolation;
            obj.allOptVars.max_violations(1) = maxViolation;
            obj.allOptVars.betas(1) = beta;
            obj.allOptVars.fEvals(1) = 1;
            obj.allOptVars.numRedus(1) = 0;
            obj.allOptVars.seekTime(1) = 0;
            obj.allOptVars.supTime(1) = 0;
            obj.allOptVars.relObjChange(1) = Inf;
            
                    fprintf('iter \t objective \t ||res|| \t max(res) \t beta \t #f_red \t #f_evals \t t_Seek \t t_Sup \n');
                
                fprintf('%i \t %0.4f \t %0.6f \t %0.4f \t %0.4f \t %d \t %d \t %0.4f \t %0.4f \n', ...
                0, ...
                obj.allOptVars.obj_values(1), ...
                obj.allOptVars.norm2_violations(1), ...
                obj.allOptVars.max_violations(1), ...
                obj.allOptVars.betas(1), ...
                obj.allOptVars.fEvals(1),...
                obj.allOptVars.numRedus(1), ...
                obj.allOptVars.seekTime(1), ...
                obj.allOptVars.supTime(1));
            
            obj.plotFunction();
            
            for i=1:obj.max_iter+1
                
                %Start Superiorization Step
                n=0;
                supStart = tic;
                f_evals = 0;
                
                
                
                if ~obj.ignoreObjective
                    while n<obj.num_reductions
                        
                        v=fGrad(x);
                        
                        loop = true;
                        while loop
                            l=l+1;
                            beta=(obj.alpha)^l;
                            z=x-beta*v;
                            
                            if f(z) <= fVal
                                n=n+1;
                                x=z;
                                loop=false;
                            end
                            f_evals = f_evals + 1;
                        end
                    end
                end
                supStop = toc(supStart);
                
                
                
                seekStart = tic;
                if obj.weighted
                    x = seeker(x, A, b, c, A_norm, weights);
                else
                    x = seeker(x, A, b, c, A_norm);
                end
                seekStop = toc(seekStart);
                
                fVal = f(x);
                f_evals = f_evals + 1;
                
                %Calculation of performance measures
                violation = (max(0, x'*A-b')+ max(0, c'-x'*A))'./A_norm;
                normViolation = (1/n_A)*norm(violation);
                maxViolation = norm(violation,'inf');
                
                obj.allObjectiveFunctionValues(i+1) = fVal;
                obj.allConstraintViolations(i+1) = maxViolation;
                obj.timeIter(i+1) = toc(obj.timeStart);
                
                obj.allOptVars.obj_values(i+1) = fVal;
                obj.allOptVars.norm2_violations(i+1) = normViolation;
                obj.allOptVars.max_violations(i+1) = maxViolation;
                obj.allOptVars.betas(i+1) = beta;
                obj.allOptVars.numRedus(i+1) = n;
                obj.allOptVars.fEvals(i+1) = f_evals;
                obj.allOptVars.seekTime(i+1) = seekStop;
                obj.allOptVars.supTime(i+1) = supStop;
                obj.allOptVars.relObjChange(i+1) = (fVal - obj.allOptVars.obj_values(i))/obj.allOptVars.obj_values(i);
                
                
                if mod(i ,10)==0
                    fprintf('iter \t objective \t ||res|| \t max(res) \t beta \t #f_red \t #f_evals \t t_Seek \t t_Sup \n');
                end
                
                fprintf('%i \t %0.4f \t %0.6f \t %0.4f \t %0.4f \t %d \t %d \t %0.4f \t %0.4f \n', ...
                    i, ...
                    obj.allOptVars.obj_values(i+1), ...
                    obj.allOptVars.norm2_violations(i+1), ...
                    obj.allOptVars.max_violations(i+1), ...
                    obj.allOptVars.betas(i+1), ...
                    obj.allOptVars.numRedus(i+1), ...
                    obj.allOptVars.fEvals(i+1), ...
                    obj.allOptVars.seekTime(i+1), ...
                    obj.allOptVars.supTime(i+1));
                
                obj.plotFunction();
                
                %Stopping rules
                if i >= 1 ...
                    && obj.allConstraintViolations(end) < obj.tol_max_violation ...
                    && obj.allOptVars.norm2_violations(end) < obj.tol_violation ...
                    && obj.allObjectiveFunctionValues < obj.tol_obj
                
                  terminated = true;
                  obj.resultInfo.stoppedby = "solution_found";
                  fprintf('Solution found. Terminating now... \n')
                  
                elseif i >= obj.accepted_iter  ...
                        &&  all(obj.allConstraintViolations(end-obj.accepted_iter:end) < obj.accepted_max_violation) ...
                        &&  all(obj.allOptVars.norm2_violations(end-obj.accepted_iter:end) < obj.accepted_violation) ...
                        &&  all(abs(obj.allOptVars.relObjChange(end-obj.accepted_iter:end)) < obj.accepted_tol_change)
                    
                   terminated = true;
                   obj.resultInfo.stoppedby = "acceptable_solution_found";
                   fprintf('Acceptable solution found. Terminating now... \n');
                   
                elseif toc(obj.timeStart) > obj.max_time
                    terminated = true;
                    obj.resultInfo.stoppedby = "time_exceeded";
                    fprintf('Time exceeded. Terminating now...\n');
                elseif i == obj.max_iter
                    terminated = true;
                    obj.resultInfo.stoppedby = "max_iter";
                    fprintf('Max number of iterations exceeded. Terminating now...\n');
                elseif obj.abortRequested
                    terminated = true;
                    obj.resultInfo.stoppedby = "user";
                    fprintf('Abort requested by user. Terminating now...\n');
                end
                
                % Stopping and returns
                if terminated
                    
                    obj.resultInfo.obj_value =  obj.allOptVars.obj_values(end);
                    obj.resultInfo.norm2_violations = obj.allOptVars.norm2_violations(end);
                    obj.resultInfo.max_violations = obj.allOptVars.max_violations(end);
                    obj.resultInfo.betas = obj.allOptVars.betas(end);
                    obj.resultInfo.time = toc(obj.timeStart);
                    
                    obj.wResult = x;
                    
                    % Deconstruct temp variables
                    clear obj.M obj.A_norm obj.weights
                    return
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
                            weight = weight ./ numel(cst{i,4}{1});
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
            
            % Choose control sequence.
            if strcmp(obj.control_sequence, 'random')
                cycle = randperm(size(A, 2));
            elseif strcmp(obj.control_sequence, 'sequential')
                cycle = 1:size(A, 2);
            elseif strcmp(obj.control_sequence, 'weight')
                [~,cycle] = sort(weights);
                cycle = cycle';
            elseif strcmp(obj.control_sequence, 'weight_inv')
                [~,cycle] = sort(weights,'descend');
                cycle = cycle';
            end
            
            for i=cycle
                
                A_x = x'*A(:, i);
                res_b=b(i)-A_x;
                res_c=A_x-c(i);
                
                if res_b < 0
                    x=x+obj.temp_weight_decay*obj.M(:, i)*res_b;
                elseif res_c <0
                    x=x-obj.temp_weight_decay*obj.M(:, i)*res_c;
                end
            end
            
            x(x<0) = 0;
            obj.temp_weight_decay = obj.temp_weight_decay* obj.weight_decay;
        end
        
        function plotFunction(obj)
            % plot objective function output
            yF = obj.allObjectiveFunctionValues;
            yViol = obj.allConstraintViolations;
            yProx = obj.allOptVars.norm2_violations;
            x = 1:numel(yF);
            
            if isempty(obj.figHandle) || ~isgraphics(obj.figHandle,'figure')
                %Create new Fiure and store axes handle
                obj.figHandle = figure('Name','Progress of Superiorization','NumberTitle','off','Color',[.5 .5 .5]);
                hFig = obj.figHandle;
                               
                titles = {'Objective Function','Max. constr. violation','Proximity'};
                yscale = {'log','log','log'};
                ylabels = {'Objective Function','Max. constr. violation','Proximity'};
                
                defaultFontSize = 14;
                               
                for p = 1:3
                    hAx(p) = subplot(1,3,p,'Parent',hFig);
                    hold(hAx(p),'on');
                    grid(hAx(p),'on');
                    grid(hAx(p),'minor');
                    set(hAx(p),'YScale',yscale{p});
                    title(hAx(p),titles{p},'LineWidth',defaultFontSize);
                    xlabel(hAx(p),'# iterations','Fontsize',defaultFontSize);
                    ylabel(hAx(p),ylabels{p},'Fontsize',defaultFontSize);
                   
                end
                   
                %Create plot handle and link to data for faster update
                hPlot(1) = plot(hAx(1),x,yF,'xb','LineWidth',1.5,'XDataSource','x','YDataSource','yF');
                hPlot(2) = plot(hAx(2),x,yViol,'xb','LineWidth',1.5,'XDataSource','x','YDataSource','yViol');
                hPlot(3) = plot(hAx(3),x,yProx,'xb','LineWidth',1.5,'XDataSource','x','YDataSource','yProx');
                               
                %Add a Stop button with callback to change abort flag
                c = uicontrol;
                cPos = get(c,'Position');
                cPos(1) = 5;
                cPos(2) = 5;
                set(c,  'String','Stop',...
                    'Position',cPos,...
                    'Callback',@(~,~) abortCallbackButton(obj));

                obj.plotHandles = hPlot;
                obj.axesHandles = hAx;
                
            else %Figure already exists, retreive from axes handle
                hFig = obj.figHandle;
                hAx = obj.axesHandles;
                hPlot = obj.plotHandles;
            end
            
            % draw updated axes by refreshing data of the plot handle (which is linked to y and y)
            % in the caller workspace. Octave needs and works on figure handles, which
            % is substantially (factor 10) slower, thus we check explicitly
            matRad_cfg = MatRad_Config.instance();
            
            if matRad_cfg.isOctave
                refreshdata(hFig,'caller');
            else
                for p = 1:3
                    refreshdata(hPlot(p),'caller');
                end
            end
            
            drawnow;
        end
        
        function abortCallbackButton(obj,~,~,~)
            obj.abortRequested = true;
        end
    end
    
    methods (Static)
        function available = IsAvailable()
            available = true; % if you made it to here the function exists.
        end
    end
end
