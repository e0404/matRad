function cst = matRad_coverageBasedCstManipulation(cst,ct,multScen,voxelWeightingType,ringCreationType)

covFlag = 0;
Counter = 0;

for  i = 1:size(cst,1)
    if ~isempty(cst{i,6})
        if sum(strcmp({cst{i,6}(:).robustness},'coverage')) > 0 
           if isempty(strfind(cst{i,2},'Ring'))
                
                covFlag = 1;
                Counter = Counter + 1;  

                % create ring structure around VOI
                    if isequal(ringCreationType, 'sampling')
                        % sample voxel probabilities
                        [voxelProbCube,shift_vox] = matRad_sampleVoxelProb(cst,ct,multScen.shiftSize,cst{i,2},1000);

                        % create cst with ring structure
                        probTreshold = 1e-4;
                        cstRing{Counter,1}           = size(cst,1) - 1 + Counter;
                        cstRing{Counter,2}           = [cst{i,2},'Ring'];
                        cstRing{Counter,3}           = cst{i,3};
                        cstRing{Counter,4}{1}        = setdiff(find(voxelProbCube > probTreshold*max(voxelProbCube(:))),cst{i,4}{1});
                        cstRing{Counter,5}           = cst{i,5};
                        cstRing{Counter,5}.voxelProb = voxelProbCube(cstRing{Counter,4}{1})';
                        cstRing{Counter,5}.shift_vox = shift_vox;
                    elseif isequal(ringCreationType, 'exact')
                        if multScen.numOfShiftScen > 1
                            
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
                            cstRing{Counter,1}           = size(cst,1) - 1 + Counter;
                            cstRing{Counter,2}           = [cst{i,2},'Ring'];
                            cstRing{Counter,3}           = cst{i,3};
                            cstRing{Counter,4}{1}        = setdiff(VOIRingidx,cst{i,4}{1});
                            cstRing{Counter,5}           = cst{i,5};
                            cstRing{Counter,5}.Priority  = cst{i,5}.Priority + 1;
                            cstRing{Counter,5}.voxelProb = voxelProbCube(cstRing{Counter,4}{1})';
                                                  
                        end
                        
                    end
                    
                    if isequal(voxelWeightingType,'heurWeighting')
                        cstRing{Counter,5}.voxelWeightingType = 'heurWeighting';
                    elseif isequal(voxelWeightingType,'probWeighting')
                        cstRing{Counter,5}.voxelWeightingType = 'probWeighting';
                    end
                    
                    % calculate mininum distance to VOI (heurWeighting)
                    [cstRing{Counter,5}.minDistToVOI,~] = matRad_calcMinDist(ct,cstRing{Counter,4}{1},cst{i,4}{1});

                % generate coverage based objectives/constraints
                    % pass coverage based objective/constraint specification to VOI ring structure
                    logidx             = strcmp({cst{i,6}(:).robustness},'coverage');
                    cstRing{Counter,6} = cst{i,6}(logidx);
                    cst{i,6}           = cst{i,6}(~logidx);

                    if isequal(cstRing{Counter,6}.type,'min DCH objective')

                        % add square underdosing objective to VOI
                        [tmp.type,tmp.dose,tmp.penalty,tmp.EUD,tmp.volume,tmp.coverage,tmp.robustness] = deal('square underdosing',cstRing{Counter,6}(1).dose,100,NaN,NaN,NaN,'none');
                        cst{i,6}(end + 1,1) = tmp;

                        % add square overdosing objective to VOI
                        [tmp.type,tmp.dose,tmp.penalty,tmp.EUD,tmp.volume,tmp.coverage,tmp.robustness] = deal('square overdosing',1.07*cstRing{Counter,6}(1).dose,2250,NaN,NaN,NaN,'none');
                        cst{i,6}(end + 1,1) = tmp;

                        % add square overdosing objective to VOI ring structure
                        [tmp.type,tmp.dose,tmp.penalty,tmp.EUD,tmp.volume,tmp.coverage,tmp.robustness] = deal('square overdosing',1.07*cstRing{Counter,6}(1).dose,2250,NaN,NaN,NaN,'none');
                        cstRing{Counter,6}(end + 1,1) = tmp;

                    end
               
           end
        end 
    end
end

if covFlag
    cst = [cst;cstRing];
end

end