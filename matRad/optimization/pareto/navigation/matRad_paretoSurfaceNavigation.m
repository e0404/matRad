function v = matRad_paretoSurfaceNavigation(Y,yr,tau,j,b)
    %matRad_navigation problem allows the exploration of plans that are not in the set
    %of precalculated plans. 
    %
    % input
    %   Y:              Set of precalculated Pareto optimal plans
    %   yr:             Reference point that is used as a starting point to calculate the new plna
    %   tau:            Value that objective j is fixed for calculation of next point
    %   j:              Index of the objectives
    %   b:              Array containing the upper bounds for each objective
    %
    % output
    %   facetPoints:    vertices of the hyperplane in objective space
    %   refPoint:       Reference Point, used to construct the hyperplane
    %   normal:         Normal of the facet
    %
    %
    % References
    %   [1] DOI: 10.1088/0031-9155/53/4/011
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [k,m] = size(Y);
     

    f = zeros(m+k+1,1);
    f(end) = 1;

    %Inequality constraints
    A = zeros(k,m+k+1);
    A(:,1:m) = Y;



    %Equality constraints
    Aeq = zeros(k+1,m+k+1);
    Aeq(1:k,1:m) = Y;
    Aeq(k+1,1:m) = 1;
    Aeq(1:k,m+1:m+k) = 1;
    Aeq(1:k,end) = -1;
    Aeq(j,m+1:end) = 0;



    beq = zeros(k+1,1);
    beq(1:k,1) = yr;
    beq(j,1) = tau;
    beq(end) = 1;

    %Bounds
    lb = zeros(m+k,1);
    ub = ones(1,m);
    
    options = optimoptions('linprog','Display','iter');
    %options = optimoptions('linprog','Algorithm','interior-point','Display','iter');
    v = linprog(f,A,b,Aeq,beq,lb,ub,options); 
    if numel(v) > 0
        v = v(1:m);
    end
end