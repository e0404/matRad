function cst = matRad_coverageBasedCstManipulation(cst,ct,multScen,ringCreationType,voxelWeightingType,normalTissueRingFlag)

covFlag = 0;
Counter = 0;

for  i = 1:size(cst,1)
    if ~isempty(cst{i,6})
        if sum(strcmp({cst{i,6}(:).robustness},'coverage')) > 0 & isequal(cst{i,6}.type,'min DCH objective') |...
           sum(strcmp({cst{i,6}(:).robustness},'coverage')) > 0 & isequal(cst{i,6}.type,'max DCH objective')
       
           if isempty(strfind(cst{i,2},'Ring'))
                
                covFlag = 1;
                Counter = Counter + 1;  

                % create ring structure around VOI      
                cstRing{Counter,1} = size(cst,1) - 1 + Counter;
                cstRing{Counter,2} = [cst{i,2},' Ring'];
                cstRing{Counter,3} = cst{i,3};
                cstRing{Counter,5} = cst{i,5};
                
                if isequal(ringCreationType, 'exact')
                    if multScen.numOfShiftScen > 1
                        % use dij scenarios
                        
                        % get VOI voxel coordinates
                        [yCoordsVOI_vox, xCoordsVOI_vox, zCoordsVOI_vox] = ind2sub(ct.cubeDim,cst{i,4}{1});

                        % create empty Cube
                        voxelProbCube = zeros(ct.cubeDim);

                        VOIRingidx = cst{i,4}{1};

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

                            % build ring structure
                            VOIRingidx = union(VOIRingidx,shiftedVOIidx);

                        end

                        % create cst with ring structure
                        cstRing{Counter,4}{1}        = setdiff(VOIRingidx,cst{i,4}{1});
                        cstRing{Counter,5}.voxelProb = voxelProbCube(cstRing{Counter,4}{1})';
                    
                    elseif multScen.numOfShiftScen == 1
                        % create scnearios with shifts
                        
                        % sample voxel probabilities
                        [voxelProbCube,voxelShift,idxShift] = matRad_sampleVoxelProb(cst,ct,multScen.shiftSize,cst{i,2},1000);

                        % create cst with ring structure
                        cstRing{Counter,4}{1}         = setdiff(find(voxelProbCube > 0),cst{i,4}{1});
                        cstRing{Counter,5}.voxelProb  = voxelProbCube(cstRing{Counter,4}{1})';
                        cstRing{Counter,5}.voxelShift = voxelShift;
                        cstRing{Counter,5}.idxShift   = idxShift;

                        for k  = 1:size(cst,1)
                            cst{k,5}.voxelShift = voxelShift;
                            cst{k,5}.idxShift   = idxShift;
                        end                        
                    end                    

                elseif isequal(ringCreationType, 'predefinedMargin') % only for heuristic weighting

                    voiCube                      = zeros(ct.cubeDim);
                    voiCube(cst{i,4}{1})         = 1;
                    [margin.x,margin.y,margin.z] = deal(14); 
                    voiCube                      = matRad_addMargin(voiCube,cst,ct.resolution,margin,true);
                    cstRing{Counter,4}{1}        = find(voiCube>0);
                    
                    % create shifts
                    [~,voxelShift,idxShift]       = matRad_sampleVoxelProb(cst,ct,multScen.shiftSize,cst{i,2},1000);
                    cstRing{Counter,5}.voxelShift = voxelShift;
                    cstRing{Counter,5}.idxShift   = idxShift;

                    for k  = 1:size(cst,1)
                        cst{k,5}.voxelShift = voxelShift;
                        cst{k,5}.idxShift   = idxShift;
                    end 
                    
                else
                    error('no valid ringCreationType')
                end
                
                % set voxel weighting type (can be changed afterwards)
                if isequal(voxelWeightingType,'heurWeighting')
                    cstRing{Counter,5}.voxelWeightingType = 'heurWeighting';
                elseif isequal(voxelWeightingType,'probWeighting')
                    if isequal(ringCreationType, 'predefinedMargin')
                        cstRing{Counter,5}.voxelWeightingType = 'heurWeighting';
                        warning('in case of ringCreationType=predefinedMargin, only heurWeighting possible. set voxelWeightingType = heurWeighting')
                    else
                        cstRing{Counter,5}.voxelWeightingType = 'probWeighting';
                    end
                end
                    
                % calculate mininum distance to VOI (heurWeighting)
                [cstRing{Counter,5}.minDistToVOI,~] = matRad_calcMinDist(ct,cstRing{Counter,4}{1},cst{i,4}{1});

                % generate coverage based objectives
                % pass coverage based objective specification to VOI ring structure
                logidx             = strcmp({cst{i,6}(:).robustness},'coverage');
                cstRing{Counter,6} = cst{i,6}(logidx);
                cst{i,6}           = cst{i,6}(~logidx);
                
                % create EXP structure (VOI + VOIRing)
                Counter = Counter + 1;
                
                cstRing{Counter,1}    = size(cst,1) - 1 + Counter;
                cstRing{Counter,2}    = [cst{i,2},' EXP'];
                cstRing{Counter,3}    = cst{i,3};
                cstRing{Counter,4}{1} = union(cstRing{Counter-1,4}{1},cst{i,4}{1}); 
                cstRing{Counter,5}    = cst{i,5};

                if isequal(cstRing{Counter-1,6}.type,'min DCH objective')

                    % add square underdosing objective to VOI
                    [tmp.type,tmp.dose,tmp.penalty,tmp.EUD,tmp.volume,tmp.coverage,tmp.robustness] = deal('square underdosing',cstRing{Counter-1,6}(1).dose,100,NaN,NaN,NaN,'none');
                    cst{i,6}(end + 1,1) = tmp;

                    % add square overdosing objective to EXP structure
                    [tmp.type,tmp.dose,tmp.penalty,tmp.EUD,tmp.volume,tmp.coverage,tmp.robustness] = deal('square overdosing',round(1.07*cstRing{Counter-1,6}(1).dose*10)/10,2250,NaN,NaN,NaN,'none');
                    cstRing{Counter,6}(end + 1,1) = tmp;

                elseif isequal(cstRing{Counter,6}.type,'min DCH objective')

                    error('min DCH objective specification not defined yet')

                end
               
           end
        end 
    end
end

if covFlag
    cst = [cst;cstRing];
end

if normalTissueRingFlag
    cstidx = find(strcmp([cst(:,2)],'BODY'));
    
    cstTmp{1,1} = cst{cstidx,1};
    cstTmp{1,2} = 'nomral tissue Ring';
    cstTmp{1,3} = cst{cstidx,3};
    cstTmp{1,5} = cst{cstidx,5};
    cstTmp{1,6} = cst{cstidx,6};
    cst{cstidx,6} = [];
    
    voiCube                                                    = zeros(ct.cubeDim);
    voiCube(cst{find(strcmp([cst(:,2)],'prostate bed')),4}{1}) = 1;
    
    % inner margin
    [margin.x margin.y margin.z] = deal(20);
    voiCubeInner = matRad_addMargin(voiCube,cst,ct.resolution,margin,true);
    VInner    = find(voiCubeInner>0);

    % outer margin
    [margin.x margin.y margin.z] = deal(40);
    voiCubeOuter = matRad_addMargin(voiCube,cst,ct.resolution,margin,true);
    VOuter    = find(voiCubeOuter>0);
    
    cstTmp{1,4}{1} = setdiff(VOuter,VInner);
    
    cst = [cst;cstTmp];
end

end