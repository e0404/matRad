function [resultGUI,returnStruct] = matRad_paretoOptimization(dij,cst,pln,nIter,wInit)
% matRad inverse pareto planning wrapper function
% 
% call
%   [returnStruct] = matRad_paretoOptimization(dij,cst,pln,nIter)
%   [returnstruct] = matRad_paretoOptimization(dij,cst,pln,nIter,wInit)
%
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct
%   pln:        matRad pln struct
%   nIter:      number of iterations to run/points to calculate after anchor point generation
%   wInit:      (optional) custom weights to initialize problems
%
% output
%   retStruct:  Structure containing the weights of the final plans and other important information
%
% References
%   - https://dx.doi.org/10.2139/ssrn.1427721 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();
if exist('wInit','var')
    [dij,cst,pln,wInit,optiProb,FLAG_ROB_OPT] = matRad_initOptimization(dij,cst,pln,wInit);
else
    [dij,cst,pln,wInit,optiProb,FLAG_ROB_OPT] = matRad_initOptimization(dij,cst,pln);
end

if ~isfield(pln.propOpt,'optimizer')
    pln.propOpt.optimizer = 'IPOPT';
end

switch pln.propOpt.optimizer
    case 'IPOPT'
        optimizer = matRad_OptimizerIPOPT;
    case 'fmincon'
        optimizer = matRad_OptimizerFmincon;
    otherwise
        warning(['Optimizer ''' pln.propOpt.optimizer ''' not known! Fallback to IPOPT!']);
        optimizer = matRad_OptimizerIPOPT;
end
        
if ~optimizer.IsAvailable()
    matRad_cfg.dispError(['Optimizer ''' pln.propOpt.optimizer ''' not available!']);
end


% PARETO PART
if optimizer.options.acceptable_obj_change_tol > 1e-5
    warning(['Pareto Optimization requires more accurate results and therefore small objective change tolerance!']);
end
%get number of objectives at the start

objcount = numel(optiProb.objectives);

%% generaete Anchor Points for optimization
penGrid = matRad_generateParetoAnchorPoints(objcount);

%initialize matrices for weight vectors and associated objective function values
weights = zeros(numel(wInit),size(penGrid,1));
fInd = zeros(size(penGrid,1),objcount);
objectiveFunctionVals = {};
% loop over all anchor points

for i = 1:size(penGrid,1)
    
    optiProb.updatePenalties(penGrid(i,:));
    
    
    optimizer = optimizer.optimize(wInit,optiProb,dij,cst);
    wOpt = optimizer.wResult;
    info = optimizer.resultInfo;
    weights(:,i) = wOpt;

    %set values for warm start
    %optimizer.optionsWarmStart.use      = true;
    %optimizer.optionsWarmStart.zl       = info.zl;
    %optimizer.optionsWarmStart.zu       = info.zu;
    %optimizer.optionsWarmStart.lambda   = info.lambda;
    fInd(i,:) = matRad_objectiveFunctions(optiProb,wOpt,dij,cst)';
    objectiveFunctionVals(end + 1) = {optimizer.allObjectiveFunctionValues};
    
end



%%normalize objectives

optiProb.normalizationScheme.name = 'UL';
optiProb.normalizationScheme.U = max(fInd,[],1);
optiProb.normalizationScheme.L = min(fInd,[],1);

fInd = (fInd-optiProb.normalizationScheme.L)./(optiProb.normalizationScheme.U-optiProb.normalizationScheme.L);


%% Maybe: Grouping based on correlation of objectives
%% Generaete further points

%initialize array storing OPS boundaries
OPSA = [];
OPSb = [];
np = size(penGrid,1);

% first on: "Balanced point" after anchor points (slightly different calculation for normal)
%

[a,b,firstNormal] = matRad_normalFromFacet(fInd,1:np,1);

penGrid(np+1,:) = abs(firstNormal);
newPen = abs(firstNormal);

%update cst values
optiProb.updatePenalties(newPen); % change for groupings!

%calculate first new point

optimizer = optimizer.optimize(wInit,optiProb,dij,cst);
wOpt = optimizer.wResult;
weights = [weights,wOpt];
objectiveFunctionVals(end + 1) = {optimizer.allObjectiveFunctionValues};
%wInit = wOpt;
fIndv = matRad_objectiveFunctions(optiProb,wOpt,dij,cst)';
fIndv = optiProb.normalizeObjectives(fIndv);

fInd = [fInd;fIndv];

%initialize some variables that might be used if points are removed later on (should only happen with too high stopping criteria)
removedfInd = [];
removedpenGrid = [];
removedweights = [];
removedOPSA = [];
removedOPSb = [];


OPSA = [OPSA;-penGrid(np+1,:)];
OPSb = [OPSb;-1*fInd(np+1,:)*penGrid(np+1,:)'];
errors = [];
allErrors = {};

initSize = size(penGrid,1);

%calculating for distance measure
L = min(fInd,[],1);
U = max(fInd,[],1);
eps = U - L;

%% calculate remaining facets

for i = 1:nIter
    fprintf('Now in iteration %i',i)
    %Step 1 calculate convex Hull -> Inner approximation (IPS) and gives facets
    %Rennen Algorithm
    fVals = fInd;

    fValsMod = matRad_generateParetoDummyPoints(fVals,U); %generate dummy points
    %
    [kmod,vol] = convhulln(fValsMod);
    %[kred,vol] = convhulln(fVals);

    %check for relevant facets (those that contain points of the original
    %fVals set)
    IPSidxs = 1:size(fVals,1);
    relFacetidxs = [];
    for j = 1:size(kmod,1)
        if any(ismember(kmod(j,:),IPSidxs))
            relFacetidxs = [relFacetidxs,j];
        end
    end
    facetMods = kmod(relFacetidxs,:);
    facetErrors = zeros(size(facetMods,1),1);
    normals = zeros(size(facetMods));

    %Loop over facets and calculate error
    for j = 1:size(relFacetidxs,2)
        [facetPoints,refPoint,normal] = matRad_normalFromFacet(fValsMod,facetMods,j);

        
        %check for sign of normals (might be redundant)
        if all(normal<0)
            continue
        end
        
        %now check for OPS point for facet
        lb = min(fVals,[],1);
        ub = max(fVals,[],1);
        z = linprog(normal,OPSA,OPSb,[],[],lb,ub); 
        
        %hyperplane distance
        b = refPoint*normal;

        %calculate error for each facet
        
        facetErrors(j) = (b-z'*normal)/(eps*normal); 
        normals(j,:) = normal;

    end
    allErrors(end+1) = {facetErrors};
    [A,I] = sort(facetErrors,'descend');

    %%check for next facet to run
    found = false;
    facetNum= 1;
    accuracy = 4;

    while ~found && facetNum <= numel(I) %loop over facets
        idx = I(facetNum);
        norm = normals(idx,:);

        %check if facet has already been run
        if ~any(ismember(round(penGrid,accuracy),round(norm,accuracy),'rows'))

            %update weights in cst
            optiProb.updatePenalties(norm);
            paretoOptimal = false;
            calcPointsForFacet = 0;

            while ~paretoOptimal && calcPointsForFacet < 4 %more than 4 points?
                if calcPointsForFacet < 3
                    optimizer = optimizer.optimize(weights(:,end-calcPointsForFacet),optiProb,dij,cst);
                else
                    optimizer = optimizer.optimize(wInit,optiProb,dij,cst);
                end

                %really necessary?
                if optimizer.resultInfo.iter < 15 %sometimes optimizer will terminate fast (points tend to be not optimal!)
                    calcPointsForFacet = calcPointsForFacet + 1;
                    continue
                end
                %check number of iterations needed for optimzation
                wOpt = optimizer.wResult;
                
                fIndv = matRad_objectiveFunctions(optiProb,wOpt,dij,cst)';
                fIndv = optiProb.normalizeObjectives(fIndv);
                %how does the newly generated point influence the pareto
                %surface?
    
                dominatedPoints = matRad_paretoCheckPoint(fInd,fIndv); 
                % if no point is dominated -> returns []
                % if the newly generated point is not optimal -> Returns array with a 0 inside 
                % if (a) previously generated point(s) is/are dominated by the
                % newly generated point -> returns index in fInd -> Need to
                % update penGrid, weights, fInd,OPSa,OPSb

                if isempty(dominatedPoints)
                    paretoOptimal = true;   
                    
                elseif dominatedPoints == 0
                    calcPointsForFacet = calcPointsForFacet + 1;
                    continue %% found point is dominated by previously calculated point -> has to be recalculated

                else %newly calculated point dominates points on the previous pareto surface
                    paretoOptimal = true;
                    %remove dominated points and related equations %might
                    %have to check if this fully works
                    removedfInd  = [removedfInd;fInd(dominatedPoints,:)];
                    removedpenGrid = [removedpenGrid;penGrid(dominatedPoints,:)];
                    removedweights = [removedweights,weights(:,dominatedPoints)];


                    fInd(dominatedPoints,:) = []; 
                    penGrid(dominatedPoints,:) = []; 
                    weights(:,dominatedPoints) = [];
                    %Issue if 
                    if ~(dominatedPoints-objcount <=0)
                        removedOPSA = [removedOPSA;OPSA(dominatedPoints-objcount,:)];
                        removedOPSb = [removedOPSb;OPSb(dominatedPoints-objcount,:)];
                        OPSA(dominatedPoints-objcount,:) = []; % shifted index
                        OPSb(dominatedPoints-objcount,:) = []; 

                    end
                    %remove dominated points
                    
                end
            end
            objectiveFunctionVals(end + 1) = {optimizer.allObjectiveFunctionValues};
            errors = [errors,facetErrors(idx)];
            newPen = norm;
            %awInit = wOpt;
            info = optimizer.resultInfo;
            weights = [weights,wOpt];
            penGrid = [penGrid;newPen];

            
            %fIndv = cellfun(@(group)sum(fIndv(:,group)),groups)
            fInd = [fInd;fIndv];
                      
            found = true;
        end
        facetNum = facetNum +1;
    end    

    % when final point is found: Update OPSA and OPSb
    OPSA = [OPSA;-newPen]; %add normal vector of facet that was run 
    OPSb = [OPSb;-fIndv*newPen'];

end





%defining a structure that contains all relevant information
returnStruct.OPSA = OPSA;
returnStruct.OPSb = OPSb;
returnStruct.weights = weights;
returnStruct.finds = fInd;
returnStruct.penGrid = penGrid;
returnStruct.errors = errors;
returnStruct.allErrors = allErrors;
returnStruct.removed = {removedfInd,removedpenGrid,removedweights,removedOPSA,removedOPSb};

returnStruct.weights = weights;
returnStruct.finds = fInd;
returnStruct.penGrid = penGrid;

returnStruct.optiProb = optiProb;
returnStruct.allObj = objectiveFunctionVals;
returnStruct.modcst = cst;

%calculate a single plan
resultGUI = matRad_calcCubes(wOpt,dij);
resultGUI.wUnsequenced = wOpt;
resultGUI.usedOptimizer = optimizer;
resultGUI.info = info;
resultGUI.optiProb = optiProb;



% unblock mex files
clear mex
