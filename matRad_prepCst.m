function cst = matRad_prepCst(cst,sparecst)
% function to prepare cst for mixed mod optimization
% sets alpha beta parameters in cst{x,6}
% explicit order of objectives : underdosing and THEN overdosing
% also empty obj priorities are set to 99 
if nargin < 2
    sparecst = 0;
end
if sparecst ==0
    for i= 1:size(cst,1)
        
        if ~isempty(cst{i,6})
            for j=1:numel(cst{i,6})
                cst{i,6}{j}.alphaX = cst{i,5}.alphaX;
                cst{i,6}{j}.betaX = cst{i,5}.betaX;
            end
        else 
            cst{i,5}.Priority = 99;
        end

        
    end
else
    
%     sparecst = voiPTV;
    %
    %
    cst{sparecst,6}(2).alphaX = 0.1;
    cst{sparecst,6}(1).alphaX = 0.5;
    cst{sparecst,6}(1).betaX = 0.05;
    cst{sparecst,6}(2).betaX = 0.05;
    for i= 1:size(cst,1)
        
        if ~isempty(cst{i,6})
            if ~(i==sparecst)
                for j=1:numel(cst{i,6})
                    cst{i,6}{j}.alphaX = cst{i,5}.alphaX;
                    cst{i,6}{j}.betaX = cst{i,5}.betaX;
                end
            end
        else
            cst{i,5}.Priority = 99;
        end
        
    end
end
end