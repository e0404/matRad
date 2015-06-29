function optResult = matRad_projectedLBFGS(objFunc,wInit,visBool,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projected L-BFGS optimizer including a positivity constraints on the
% optimization variable
% 
% call
%   optResult = matRad_projectedLBFGS(objFunc,wInit)
%
% input
%   objFunc:    objective function to be optimized
%   wInit:      start solution for optimizer
%   visBool:    plots the objective function value in dependence of the
%               number of iterations
%   varargin:   optional: number of iterations and precision
%
% output
%   optResult:  struct containing the optimized fluence vector
%
% References
%   [1] Kelley: Iterative methods for optimization 1999
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% variables for termination criteria
iter      = 0;
if isempty(varargin{1,1})
    numOfIter = 1000;
    prec      = 1e-3;
else
    optParam = varargin{1,1};
    numOfIter = optParam{1,1}.numOfIter;
    prec      = optParam{1,1}.prec;
end

numOfParameters = numel(wInit);
% plot objective function output
if visBool
    figOpt=figure('Name','Progress of Optimization','NumberTitle','off');
    hold on, grid on, grid minor, xlabel('# iterations','Fontsize',14),ylabel('objective value','Fontsize',14)
    set(gca,'YScale','log');
    title('Progress of Optimization','LineWidth',14),
end
% initialize LBFGS optimizer
historyCounter = 0;
mem            = 10;        % number of past gradients and function values used for inverse hessian contruction
x              = NaN*ones(numOfParameters,mem);
x(:,1)         = wInit;

objFuncValue   = NaN*ones(1,mem);
dx             = NaN*ones(numOfParameters,mem);

s_k            = ones(numOfParameters,mem-1);
y_k            = ones(numOfParameters,mem-1);
r_k            = ones(mem-1,1);
a_k            = ones(1,mem-1);

% 1st calculation of objective function and gradient
[objFuncValue(1),dx(:,1)] = objFunc(wInit);
objFuncValue(2) = 2*objFuncValue(1);

% convergence if change in objective function smaller than prec or maximum
% number of iteration reached. no convergence if lbfgs has just been rest
continueOpt = 1;
while continueOpt == 1
    % implementation of L-BFGS according to
    % http://en.wikipedia.org/wiki/L-BFGS
    
    % inverse hessian update
    q = dx(:,1);
    for i = 1:historyCounter
        a_k(i) = r_k(i)*s_k(:,i)'*q;
        q      = q - a_k(i)*y_k(:,i);
    end
    z = s_k(:,1)'*y_k(:,1)/(y_k(:,1)'*y_k(:,1))*q; % this corresponds to H*q where H is approximated as described in Nocedal 9.1
    for i = historyCounter:-1:1
        b = r_k(i)*y_k(:,i)'*z;
        z = z + s_k(:,i)*(a_k(i)-b);
    end
    
    % obtain search direction
    dir = -z; % this is the lbfgs direction!
    
    % 2.1 armijo linesearch to to find acceptable stepsize alpha
    alpha        = 10;
    fac          = 1/10; % < 1 reduction factor of alpha
    c_1          = 1e-10;

    continueLineSearch = true;
    
    expectedDescend = ((x(:,1)>0).*dir)'*dx(:,1);
    
    fprintf('Starting line search ')
    while continueLineSearch
        
        alpha = fac*alpha;
        
        candidateX = x(:,1) + alpha*dir;

        % project candidate to feasible set
        candidateX(candidateX<0) = 0;
        
        lineSearchObjFuncValue = objFunc(candidateX);
        
        continueLineSearch = lineSearchObjFuncValue > objFuncValue(1) + c_1*alpha*expectedDescend;
    
        if alpha < 1e-10;
            %fprintf('Error in Line search - alpha close to working precision...\n');
            fprintf(1,'Performed 10 line searches - Resetting LBFGS update\n');
        
            s_k = ones(numOfParameters,mem-1);
            y_k = ones(numOfParameters,mem-1);
            r_k = ones(mem-1,1);
            a_k = ones(1,mem-1);

            historyCounter = 0;
            break;
        end
    
        fprintf('.')
    end
    fprintf('\n')
        
    % 2.2 update x
    x(:,2:end)      = x(:,1:end-1);
    dx(:,2:end)     = dx(:,1:end-1);
    x(:,1)          = candidateX;

    objFuncValue(2) = objFuncValue(1);
    
    [objFuncValue(1),dx(:,1)] = objFunc(x(:,1));
        
    s_k = -diff(x,[],2);
    y_k = -diff(dx,[],2);

    s_k(x(:,1)<=0,1) = 0;
    y_k(x(:,1)<=0,1) = 0;
    
    r_k = 1./diag(y_k'*s_k);
        
    % increment iteration counter
    iter = iter + 1;
    historyCounter = min(historyCounter+1,mem-1);
    
    if (s_k(:,1)'*y_k(:,1)) <= 0
        
        fprintf(1,'Resetting LBFGS update\n');
        
        s_k = ones(numOfParameters,mem-1);
        y_k = ones(numOfParameters,mem-1);
        r_k = ones(mem-1,1);
        a_k = ones(1,mem-1);

        historyCounter = 0;
        
    end
    
    fprintf(1,'Iteration %d: alpha = %f, Obj func = %f\n',iter,alpha,objFuncValue(1));

    continueOpt = (iter < numOfIter && abs((objFuncValue(2)-objFuncValue(1))/objFuncValue(1))>prec) || historyCounter < 2;
    if  objFuncValue(2)== 0 && objFuncValue(1) == 0 
        continueOpt = 0;
        disp('!!! please review your constraints !!!')
    end
    
    if visBool
        objFuncValues{iter}=objFuncValue(1);
        plot(1:1:iter,cell2mat(objFuncValues),'b','Linewidth',3); drawnow
    end
    
end

fprintf(['\n' num2str(iter) ' iteration(s) performed to converge\n'])

optResult.w = x(:,1);
