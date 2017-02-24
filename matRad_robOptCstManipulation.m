function cst = matRad_robOptCstManipulation(ct,cst,pln,multScen,targetExpansion,normalTissueRingVOIName,varargin)

covFlag = 0;
Counter = 0;


% determine scenarios used for optimization
for i = 1:size(cst,1)
    for j = 1:numel(cst{i,6})
       if ~isfield(cst{i,6}(j),'coverage')
           cst{i,6}(j).coverage = 0;    % ensure backwards compability
       end
    end
end


for  i = 1:size(cst,1)
    
    if ~isempty(cst{i,6})
        
        if strcmp(pln.robOpt,'coverage') % sum(strcmp({cst{i,6}(:).robustness},'coverage'))
                
            covFlag = 1;
            Counter = Counter + 1;  

            % create VOI scenario union structure
            cstVOIScenUnion{Counter,1}          = size(cst,1) - 1 + Counter;
            cstVOIScenUnion{Counter,2}          = [cst{i,2},' ScenUnion'];
            cstVOIScenUnion{Counter,3}          = cst{i,3};
            cstVOIScenUnion{Counter,5}          = cst{i,5};
            cstVOIScenUnion{Counter,5}.Priority = NaN;
            cstVOIScenUnion{Counter,6}          = [];


            if multScen.numOfShiftScen > 1
                % create Scenario Union via dij shift scenarios

                % shift VOI according to dij shifts
                [yCoordsVOI_vox, xCoordsVOI_vox, zCoordsVOI_vox] = ind2sub(ct.cubeDim,cst{i,4}{1});
                voxelProbCube                                    = zeros(ct.cubeDim);

                for k = 1:multScen.numOfShiftScen

                    % round shifts to voxel dimensions
                    xShift_vox = round(multScen.shifts(1,k)/ct.resolution.x);
                    yShift_vox = round(multScen.shifts(2,k)/ct.resolution.y);
                    zShift_vox = round(multScen.shifts(3,k)/ct.resolution.z);

                    shiftedVOIidx = sub2ind(ct.cubeDim, yCoordsVOI_vox - yShift_vox,...
                                                        xCoordsVOI_vox - xShift_vox,...
                                                        zCoordsVOI_vox - zShift_vox);

                    % fill prob cube                                   
                    voxelProbCube(shiftedVOIidx) = voxelProbCube(shiftedVOIidx) + 1/multScen.numOfShiftScen;                                   

                end

                % save voxelIDs and voxel probabilities
                cstVOIScenUnion{Counter,4}{1}        = find(voxelProbCube > 0);
                cstVOIScenUnion{Counter,5}.voxelID   = cstVOIScenUnion{Counter,4}{1};
                cstVOIScenUnion{Counter,5}.voxelProb = voxelProbCube(cstVOIScenUnion{Counter,4}{1})';

            elseif multScen.numOfShiftScen == 1
                % create Scenario Union via VOI shift scenarios

                % sample VOI shifts
                cstVOIScenUnion{Counter,5}.VOIShift = matRad_sampleVOIShift(cst,ct,multScen.shiftSD,cst{i,2},multScen.numOfIntSegShiftScen);
                cst{i,5}.VOIShift                   = cstVOIScenUnion{Counter,5}.VOIShift; 

                % save voxelIDs and voxel probabilities
                cstVOIScenUnion{Counter,4}{1}        = find(cstVOIScenUnion{Counter,5}.VOIShift.voxelProbCube > 0);
                cstVOIScenUnion{Counter,5}.voxelID   = cstVOIScenUnion{Counter,4}{1};
                cstVOIScenUnion{Counter,5}.voxelProb = cstVOIScenUnion{Counter,5}.VOIShift.voxelProbCube(cstVOIScenUnion{Counter,4}{1});

            end                                                   
        end 
    end
end

if covFlag
    cst = [cst(:,1:6);cstVOIScenUnion];
end

if targetExpansion > 0
    for  i = 1:size(cst,1)
        if ~isempty(cst{i,6})
            if strcmp(pln.robOpt,'coverage') && isequal(cst{i,3},'TARGET')
                
                cstidx = find(strcmp([cst(:,2)],cst{i,2}));
                
                cstTmp{1,1} = cst{end,1} + 1;
                cstTmp{1,2} = [cst{cstidx,2},' initial expansion'];
                cstTmp{1,3} = cst{cstidx,3};
                cstTmp{1,5} = cst{cstidx,5};
                cstTmp{1,6} = cst{cstidx,6};
                
                % set penalties to zero (only need for initialization)
                [cstTmp{1,6}(:).penalty]    = deal(0);
                [cstTmp{1,6}(:).type]       = deal('square deviation');
                [cstTmp{1,6}(:).robustness] = deal('none');
                
                voiCube                   = zeros(ct.cubeDim);
                voiCube(cst{cstidx,4}{1}) = 1;
                
                % expand VOI
                [margin.x, margin.y, margin.z] = deal(targetExpansion);
                VOIExpansion   = matRad_addMargin(voiCube,cst,ct.resolution,margin,true);
                cstTmp{1,4}{1} = find(VOIExpansion>0); 
                
                cst = [cst;cstTmp];
                
            end
        end
    end
end

if ischar(normalTissueRingVOIName)  
    cstidx = find(strcmp([cst(:,2)],'BODY'));
    
    cstTmp{1,1} = cst{end,1} + 1;
    cstTmp{1,2} = 'normal tissue Ring';
    cstTmp{1,3} = cst{cstidx,3};
    cstTmp{1,5} = cst{cstidx,5};
    cstTmp{1,6} = cst{cstidx,6};
    cst{cstidx,6} = [];
    
    voiCube                                                             = zeros(ct.cubeDim);
    voiCube(cst{find(strcmp([cst(:,2)],normalTissueRingVOIName)),4}{1}) = 1;
    
    % inner margin
    [margin.x, margin.y, margin.z] = deal(varargin{1});
    voiCubeInner = matRad_addMargin(voiCube,cst,ct.resolution,margin,true);
    VInner    = find(voiCubeInner>0);

    % outer margin
    [margin.x, margin.y, margin.z] = deal(varargin{2});
    voiCubeOuter = matRad_addMargin(voiCube,cst,ct.resolution,margin,true);
    VOuter    = find(voiCubeOuter>0);
    
    cstTmp{1,4}{1} = setdiff(VOuter,VInner);
    
    cst = [cst;cstTmp];
end

end