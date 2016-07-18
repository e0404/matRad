function c = matRad_constFuncWrapper(w,dij,cst,type)

global matRad_DCH_ScenarioFlag;
% % A1 % 
% global matRad_DVH_Scaling;
% global matRad_DCH_Scaling;
% % A1 % 

% get current dose / effect / RBExDose vector
d = matRad_backProjection(w,dij,type);

% Initializes constraints
c = [];

% compute objective function for every VOI.
for  i = 1:size(cst,1)

    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            % only perform computations for constraints
            if ~isempty(strfind(cst{i,6}(j).type,'constraint'))
                
                % compute reference
                if (~isequal(cst{i,6}(j).type, 'max dose constraint') && ~isequal(cst{i,6}(j).type, 'min dose constraint') &&...
                    ~isequal(cst{i,6}(j).type, 'min mean dose constraint') && ~isequal(cst{i,6}(j).type, 'max mean dose constraint') &&...
                    ~isequal(cst{i,6}(j).type, 'min max mean dose constraint') && ~isequal(cst{i,6}(j).type, 'min EUD constraint') &&...
                    ~isequal(cst{i,6}(j).type, 'max EUD constraint') && ~isequal(cst{i,6}(j).type, 'min max EUD constraint')) &&...
                    isequal(type,'effect')
                     
                    d_ref = cst{i,5}.alphaX*cst{i,6}(j).dose + cst{i,5}.betaX*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end

                % if conventional opt: just add constraints of nominal dose
                if strcmp(cst{i,6}(j).robustness,'none')

                    d_i = d{1}(cst{i,4}{1});

                    c = [c; matRad_constFunc(d_i,cst{i,6}(j),d_ref)];

                % if prob opt or voxel-wise worst case: add constraints of all dose scenarios
                elseif strcmp(cst{i,6}(j).robustness,'probabilistic') || strcmp(cst{i,6}(j).robustness,'voxel-wise worst case')
                    
                    for k = 1:dij.numOfScenarios
                        
                        d_i = d{k}(cst{i,4}{1});
                        
                        c = [c; matRad_constFunc(d_i,cst{i,6}(j),d_ref)];
                        
                    end
                    
                % if coveraged based opt   
                elseif strcmp(cst{i,6}(j).robustness,'coverage')
                    
                    if isequal(cst{i,6}(j).type, 'max DCH constraint') || ...
                       isequal(cst{i,6}(j).type, 'min DCH constraint')
                    
                        for k = 1:dij.numOfScenarios

                            % get current dose
                            d_i = d{k}(cst{i,4}{1});

                            % inverse DVH calculation
                            d_pi(k) = matRad_calcInversDVH(cst{i,6}(j).volume/100,d_i);

                        end

                        c = [c; matRad_constFunc(d_i,cst{i,6}(j),d_ref,d_pi)];
                    
                    elseif isequal(cst{i,6}(j).type, 'max DCH constraint2') || ...
                           isequal(cst{i,6}(j).type, 'min DCH constraint2')                       
                        
                        d_i = [];
                        
                        % get cst index of VOI that corresponds to VOI ring
                        cstidx = find(strcmp(cst(:,2),cst{i,2}(1:end-5)));
                       
                        % get dose of VOI that corresponds to VOI ring
                        for k = 1:dij.numOfScenarios
                            d_i{k} = d{k}(cst{cstidx,4}{1});
                        end

                        % calc invers DCH of VOI
                        refQ   = cst{i,6}(j).coverage/100;
                        refVol = cst{i,6}(j).volume/100;
                        d_ref2 = matRad_calcInversDCH(refVol,refQ,d_i,dij.numOfScenarios);

                        % get dose of VOI ring
                        d_i    = d{1}(cst{i,4}{1});

                        % calc voxel dependent weighting
                        voxelWeighting = 5*cst{i,5}.voxelProb;

                        c = [c; matRad_constFunc(d_i,cst{i,6}(j),d_ref,1,d_ref2,voxelWeighting)];
                        
                    elseif isequal(cst{i,6}(j).type, 'max DCH constraint3') || ...
                           isequal(cst{i,6}(j).type, 'min DCH constraint3')
                       
                        % calculate scenario approximation scaling
                        if dij.numOfScenarios > 1
                            for k = 1:dij.numOfScenarios

                                % get current dose
                                d_i = d{k}(cst{i,4}{1});

                                % calculate volume
                                volume_pi(k) = sum(d_i >= d_ref)/numel(d_i);
                            end
                            
                            % calculate coverage probabilty
                            scenProb = 1/dij.numOfScenarios;  % assume scenarios with equal probabilities
                            
                        else
                            
%                             % A2 % 
%                             % calculate logistic function scaling and volumes
%                             DVHScaling = matRad_DVH_Scaling(j);
%                             % A2 %
                            
                            for k = 1:cst{i,5}.VOIShift.ncase
                                
                                % get current dose
                                if isequal(cst{i,5}.VOIShift.shiftType,'rounded')
                                    d_i = d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.roundedShift.idxShift(k));
                                elseif isequal(cst{i,5}.VOIShift.shiftType,'linInterp')
                                    error('linInterp for constraints not implemented yet')
                                end
                                
%                                 % A3 %                                
%                                 volume_pi(k) = sum(1./(1+exp(-2*DVHScaling*(d_i-d_ref))))/numel(d_i);
%                                 % A3 %
                                    
                                % B1 %
                                % calculate volumes
                                volume_pi(k) = sum(d_i >= d_ref)/numel(d_i);
                                % B1 %  
                            end
                            
                            % calculate coverage probabilty
                            scenProb = 1/cst{i,5}.VOIShift.ncase;  % assume scenarios with equal probabilities 
                            
                        end
                        
%                         % A4 % 
%                         % calculate logistic function scaling and coverage probability
%                         DCHScaling = matRad_DCH_Scaling(j);
%                         c          = [c; sum(scenProb*(1./(1+exp(-2*DCHScaling*(volume_pi - cst{i,6}(j).volume/100)))))];
%                         % A4 % 
                        
                        % B2 %
                        % calculate coverage probability
                        c = [c; sum(scenProb*(volume_pi >= cst{i,6}(j).volume/100))];
                        % B2 %

                    elseif isequal(cst{i,6}(j).type, 'max DCH constraint4') || ...
                           isequal(cst{i,6}(j).type, 'min DCH constraint4')
                       
                       % update scenario flag
                       for k = 1:dij.numOfScenarios
                             
                            % get current dose
                            d_i = d{k}(cst{i,4}{1});
                            
                            % calculate volume and dose of scenario k that
                            % correspomd to ref values
                            %volume(k) = matRad_constFunc(d_i,cst{i,6}(j),d_ref);
                            dose(k)   = matRad_calcInversDVH(cst{i,6}(j).volume/100,d_i);
                            
                       end
                       
%                        volume_sorted = sort(volume(2:end),'descend'); 
%                        idx           = ceil(round((cst{i,6}(j).coverage/100 - 1/dij.numOfScenarios)*numel(volume)*10)/10);
                       
                       dose_sorted = sort(dose(2:end),'descend'); 
                       idx         = ceil(round((cst{i,6}(j).coverage/100 - 1/dij.numOfScenarios)*numel(dose)*10)/10);
                       
                       if idx == 0
                           
                            matRad_DCH_ScenarioFlag = [true, false(1,dij.numOfScenarios-1)];
                           
                       else
                           
%                             volumeTmp               = [volume(1),volume_sorted(1:idx)];            
%                             matRad_DCH_ScenarioFlag = ismember(volume,volumeTmp);
                            doseTmp                 = [dose(1),dose_sorted(1:idx)];            
                            matRad_DCH_ScenarioFlag = ismember(dose,doseTmp);
                       
                       end
                       
                       for k = 1:dij.numOfScenarios
                           if matRad_DCH_ScenarioFlag(k)
                                d_i = d{k}(cst{i,4}{1});
                                c   = [c;matRad_constFunc(d_i,cst{i,6}(j),d_ref)];
                           else
                                c = [c;cst{i,6}(j).volume/100];
                           end
                            
                       end
                       
                    elseif isequal(cst{i,6}(j).type, 'max DCH constraint5') || ...
                           isequal(cst{i,6}(j).type, 'min DCH constraint5')
                       
                       % update scenario flag
                       for k = 1:cst{i,5}.VOIShift.ncase
                            if isequal(cst{i,5}.VOIShift.shiftType,'rounded')
                               doseVec = d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.roundedShift.idxShift(k));

                            elseif isequal(cst{i,5}.VOIShift.shiftType,'linInterp')
                                % lin interpolation in x
                                c00 = d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X0Y0Z0(k)).*(1-cst{i,5}.VOIShift.linInterpShift.idxShift.x(k)) +...
                                      d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X1Y0Z0(k)).*cst{i,5}.VOIShift.linInterpShift.idxShift.x(k);
                                c01 = d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X0Y0Z1(k)).*(1-cst{i,5}.VOIShift.linInterpShift.idxShift.x(k)) +...
                                      d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X1Y0Z1(k)).*cst{i,5}.VOIShift.linInterpShift.idxShift.x(k);
                                c10 = d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X0Y1Z0(k)).*(1-cst{i,5}.VOIShift.linInterpShift.idxShift.x(k)) +...
                                      d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X1Y1Z0(k)).*cst{i,5}.VOIShift.linInterpShift.idxShift.x(k);
                                c11 = d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X0Y1Z1(k)).*(1-cst{i,5}.VOIShift.linInterpShift.idxShift.x(k)) +...
                                      d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X1Y1Z1(k)).*cst{i,5}.VOIShift.linInterpShift.idxShift.x(k);

                                % lin interpolation in y  
                                c0  = c00.*(1-cst{i,5}.VOIShift.linInterpShift.idxShift.y(k))+c10.*cst{i,5}.VOIShift.linInterpShift.idxShift.y(k);
                                c1  = c01.*(1-cst{i,5}.VOIShift.linInterpShift.idxShift.y(k))+c11.*cst{i,5}.VOIShift.linInterpShift.idxShift.y(k);

                                % lin interpolation in z
                                doseVec = c0.*(1-cst{i,5}.VOIShift.linInterpShift.idxShift.z(k))+c1.*cst{i,5}.VOIShift.linInterpShift.idxShift.z(k);
                                
                            end
                           
                            volume(k)     = matRad_constFunc(doseVec,cst{i,6}(j),d_ref);
                            dose(k)       = matRad_calcInversDVH(cst{i,6}(j).volume/100,doseVec);
%                             DVHdevArea(k) = (abs((volume(k)*100-cst{i,6}(j).volume)/cst{i,6}(j).volume)) * (abs((dose(k)-d_ref)/d_ref));
                            DVHdevArea(k) = volume(k) - volume(1);
%                             DVHdevArea(k) = dose(k) - dose(1);
                       end
                       
                       if isequal(cst{i,6}(j).type, 'max DCH constraint5')
                           [DVHdevAreaSorted,DVHdevAreaSortedidx] = sort(DVHdevArea(2:end),'descend');
                       elseif isequal(cst{i,6}(j).type, 'min DCH constraint5')
                           [DVHdevAreaSorted,DVHdevAreaSortedidx] = sort(DVHdevArea(2:end),'ascend');
                       end
                        
                        idx                                    = ceil((cst{i,6}(j).coverage/100 - 1/cst{i,5}.VOIShift.ncase)*numel(DVHdevArea));
                        
                       if idx == 0
                           
                            matRad_DCH_ScenarioFlag = [true, false(1,cst{i,5}.VOIShift.ncase-1)];
                           
                       else
                            matRad_DCH_ScenarioFlag = [true, false(1,cst{i,5}.VOIShift.ncase-1)];
%                             doseTmp                 = [dose(1),dose_sorted(1:idx)];            
%                             matRad_DCH_ScenarioFlag = ismember(dose,doseTmp);
%                             DVHdevAreaTmp           = [DVHdevArea(1),DVHdevAreaSorted(1:idx)];
%                             for m = 1:length(DVHdevAreaTmp)
%                                 idx = find(ismember(DVHdevArea,DVHdevAreaTmp(m)));
%                                 matRad_DCH_ScenarioFlag(idx(1)) = true;
%                             end
%                             matRad_DCH_ScenarioFlag = ismember(DVHdevArea,DVHdevAreaTmp);
%                             matRad_DCH_ScenarioFlag = [true false(1,size(cst{1,5}.voxelShift,2)-1)];
                            matRad_DCH_ScenarioFlag(DVHdevAreaSortedidx(1:idx)+1) = true;
                            
                       end                       

                       for k = 1:cst{i,5}.VOIShift.ncase
                           if matRad_DCH_ScenarioFlag(k)
                                if isequal(cst{i,5}.VOIShift.shiftType,'rounded')
                                   doseVec = d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.roundedShift.idxShift(k));

                                elseif isequal(cst{i,5}.VOIShift.shiftType,'linInterp')
                                    % lin interpolation in x
                                    c00 = d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X0Y0Z0(k)).*(1-cst{i,5}.VOIShift.linInterpShift.idxShift.x(k)) +...
                                          d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X1Y0Z0(k)).*cst{i,5}.VOIShift.linInterpShift.idxShift.x(k);
                                    c01 = d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X0Y0Z1(k)).*(1-cst{i,5}.VOIShift.linInterpShift.idxShift.x(k)) +...
                                          d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X1Y0Z1(k)).*cst{i,5}.VOIShift.linInterpShift.idxShift.x(k);
                                    c10 = d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X0Y1Z0(k)).*(1-cst{i,5}.VOIShift.linInterpShift.idxShift.x(k)) +...
                                          d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X1Y1Z0(k)).*cst{i,5}.VOIShift.linInterpShift.idxShift.x(k);
                                    c11 = d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X0Y1Z1(k)).*(1-cst{i,5}.VOIShift.linInterpShift.idxShift.x(k)) +...
                                          d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.linInterpShift.idxShift.X1Y1Z1(k)).*cst{i,5}.VOIShift.linInterpShift.idxShift.x(k);

                                    % lin interpolation in y  
                                    c0  = c00.*(1-cst{i,5}.VOIShift.linInterpShift.idxShift.y(k))+c10.*cst{i,5}.VOIShift.linInterpShift.idxShift.y(k);
                                    c1  = c01.*(1-cst{i,5}.VOIShift.linInterpShift.idxShift.y(k))+c11.*cst{i,5}.VOIShift.linInterpShift.idxShift.y(k);

                                    % lin interpolation in z
                                    doseVec = c0.*(1-cst{i,5}.VOIShift.linInterpShift.idxShift.z(k))+c1.*cst{i,5}.VOIShift.linInterpShift.idxShift.z(k);

                                end                            
                                c = [c;matRad_constFunc(doseVec,cst{i,6}(j),d_ref)];
                           else
                                c = [c;cst{i,6}(j).volume/100];
                           end
                            
                       end
                       
                    end

                end % if we are in the nominal sceario or rob opt
            
            end

        end % over all defined constraints & objectives

    end % if structure not empty and oar or target

end % over all structures

global cScaling
 c = cScaling.*c;