function cst = matRad_computeAllVoiSurfaces(ct,cst)
% matRad function that computes all VOI surfaces
%
% call
%   cst = matRad_computeAllVoiSurfaces(ct,cst)
%
% input
%   ct  matRad ct struct
%   cst matRad cst struct
%
% output
%   cst the new cst with the column containing the precomputed surface
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Computing 3D Surfaces...');

    % initialize waitbar
    figureWait = waitbar(0,'Computing 3D Surfaces...');
    % prevent closure of waitbar and show busy state
    set(figureWait,'pointer','watch');
    
    xCoord = ct.resolution.x * double(1:ct.cubeDim(2));
    yCoord = ct.resolution.y * double(1:ct.cubeDim(1));
    zCoord = ct.resolution.z * double(1:ct.cubeDim(3));
    
    [xMesh,yMesh,zMesh] = meshgrid(xCoord,yCoord,zCoord);

    numVois = size(cst,1);
    

    
    for s = 1:numVois
        mask = zeros(ct.cubeDim); % create zero cube with same dimeonsions like dose cube
        mask(cst{s,4}{1}) = 1;
        
        %Smooth the VOI
        v = smooth3(mask,'gaussian',[5 5 5],2);
        isoSurface = isosurface(xMesh,yMesh,zMesh,v,0.5);
        
        %reduce the complexity
        isoSurface = reducepatch(isoSurface,0.05);
    
        %compute isonormals
        isoNormals = isonormals(xMesh,yMesh,zMesh,v,isoSurface.vertices);
        
        cst{s,8}{1} = isoSurface;
        cst{s,8}{2} = isoNormals;
        matRad_progress(s,numVois);
        waitbar(s/numVois,figureWait);
    end    
    try
    % wait 0.1s for closing all waitbars
    allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar'); 
    delete(allWaitBarFigures);
    pause(0.1); 
    catch
        
    end
end

