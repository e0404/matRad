function [f, g] = matRad_DAOobjFunc(shapeInfoVect,shapeInfo,addInfoVect,dij,cst)

% update shapeInfoVect und indVect
[shapeInfo] = tk_updateShapeInfo(shapeInfo,shapeInfoVect);
indVect = tk_createIndVect(shapeInfoVect,addInfoVect,shapeInfo);

% calculate w (bixel weight vector)
totalNumOfBixels = shapeInfo.totalNumOfBixels;
beamNumVect = dij.beamNum;
w = zeros(totalNumOfBixels,1);
for i = 1:totalNumOfBixels    
    % find the bixelposition corresponding to the bixel i    
        % 1. find the corresponding beam
        beamNum = beamNumVect(i);        
        % 2. find index in the MLC map
        MLCPosInd = shapeInfo.beam(beamNum).bixelIndMap == i;    
    % add up the fluence from every shape of this beam        
        for j=1:shapeInfo.beam(beamNum).numOfShapes
            w(i) = w(i) + shapeInfo.beam(beamNum).shape(j).weight * ...
                    shapeInfo.beam(beamNum).shape(j).shapeMap(MLCPosInd);
        end    
end     


% Calculate dose
d = dij.physicalDose*w;

% Numbers of voxels
numVoxels = size(dij.physicalDose,1);

% Initializes f
f = 0;

% Initializes delta
delta_underdose = zeros(numVoxels,1);
delta_overdose  = zeros(numVoxels,1);
delta_deviation = zeros(numVoxels,1);
delta_mean      = zeros(numVoxels,1);
delta_EUD       = zeros(numVoxels,1);

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
    
        % get dose vector in current VOI
        d_i = d(cst{i,4});
                
        % loop over the number of constraints for the current VOI
        for j = 1:size(cst{i,6},1)
            
            % get Penalty
            rho = cst{i,6}(j).parameter(1);
            
            if isequal(cst{i,6}(j).type, 'square underdosing') 
               
                if ~isequal(cst{i,3},'OAR')
                    % underdose : Dose minus prefered dose
                    underdose = d_i - cst{i,6}(j).parameter(2);

                    % apply positive operator
                    underdose(underdose>0) = 0;

                    % calculate objective function
                    f = f + (rho/size(cst{i,4},1))*(underdose'*underdose);

                    % calculate delta
                    delta_underdose(cst{i,4}) = delta_underdose(cst{i,4}) +...
                        (rho/size(cst{i,4},1))*underdose;
                else
                    disp(['square underdosing constraint for ' cst{i,2} ' will be skipped'])
                end
                
            elseif isequal(cst{i,6}(j).type, 'square overdosing')
                
                    % overdose : Dose minus prefered dose
                    overdose = d_i - cst{i,6}(j).parameter(2);

                    % apply positive operator
                    overdose(overdose<0) = 0;

                    % calculate objective function
                    f = f + (rho/size(cst{i,4},1))*(overdose'*overdose);

                    %calculate delta
                    delta_overdose(cst{i,4}) = delta_overdose(cst{i,4}) + ...
                        (rho/size(cst{i,4},1))*overdose;
                
           elseif isequal(cst{i,6}(j).type, 'square deviation')
               
               if ~isequal(cst{i,3},'OAR')
                    % deviation : Dose minus prefered dose
                    deviation = d_i - cst{i,6}(j).parameter(2);

                    % claculate objective function
                    f = f + (rho/size(cst{i,4},1))*(deviation'*deviation);

                    % calculate delta
                    delta_deviation(cst{i,4}) = delta_deviation(cst{i,4}) +...
                        (rho/size(cst{i,4},1))*deviation;
                else
                    disp(['square deviation constraint for ' cst{i,2} ' will be skipped'])
                end
                
            elseif isequal(cst{i,6}(j).type, 'mean')              
                
                if ~isequal(cst{i,3},'TARGET')
                    % calculate objective function
                    f = f + (rho/size(cst{i,4},1))*sum(d_i);

                    % calculate delta
                    delta_mean(cst{i,4}) = delta_mean(cst{i,4}) + ...
                        (rho/size(cst{i,4},1))*ones(size(cst{i,4},1),1);
                else
                    disp(['mean constraint for ' cst{i,2} ' will be skipped'])
                end
                
            elseif isequal(cst{i,6}(j).type, 'EUD') 
               
               if ~isequal(cst{i,3},'TARGET')
                    % get exponent for EUD
                    exponent = cst{i,6}(j).parameter(2);

                    % calculate objective function and delta
                    if sum(d_i.^exponent)>0

                        f = f + rho*nthroot((1/size(cst{i,4},1))*sum(d_i.^exponent),exponent);

                        delta_EUD(cst{i,4}) = delta_EUD(cst{i,4}) + ...
                            rho*nthroot(1/size(cst{i,4},1),exponent) * sum(d_i.^exponent)^((1-exponent)/exponent) * (d_i.^(exponent-1));

                    end

                    if sum(~isfinite(delta_EUD)) > 0 % check for inf and nan for numerical stability
                        error(['EUD computation for ' cst{i,2} ' failed. Reduce exponent to resolve numerical problems.']);
                    end
               else
                    disp(['EUD constraint for ' cst{i,2} ' will be skipped'])
                end
               
            else
                
                error('undefined objective in cst struct');
                
            end
      
        end
        
    end
end
   
if nargout > 1
    % Calculate gradient of bixel weights
    bixelGrad = (( 2*(delta_underdose + delta_overdose + delta_deviation) + delta_mean + delta_EUD )' * dij.physicalDose)';

    g = tk_getGradients(bixelGrad,addInfoVect,indVect,shapeInfo);
end

end