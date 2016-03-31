function f = matRad_objFuncWrapper(w,dij,cst,type)

% get current dose / effect / RBExDose vector
d = matRad_backProjection(w,dij,type);

% Initializes f
f = 0;

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
        
            % reference dose/effect/RBExDose
            if isempty(strfind(cst{i,6}(j).type,'constraint'))

                if (~isequal(cst{i,6}(j).type, 'mean') && ~isequal(cst{i,6}(j).type, 'EUD')) &&...
                    isequal(type,'effect') 

                    d_ref = dij.ax(cst{i,4}).*cst{i,6}(j).dose + dij.bx(cst{i,4})*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
            else
                d_ref = [];
            end
       
            if strcmp(cst{i,6}(j).robustness,'none')
                
                d_i = d{1}(cst{i,4});
                
                f = f + matRad_objFunc(d_i,cst{i,6}(j),d_ref);
                
            elseif strcmp(cst{i,6}(j).robustness,'probabilistic')
            elseif strcmp(cst{i,6}(j).robustness,'objective-wise worst case')
            elseif strcmp(cst{i,6}(j).robustness,'voxel-wise worst case')
                [d_max,maxDoseIx] = max([d{:}],[],2);
                [d_min,minDoseIx] = min([d{:}],[],2);
            end
        end
            
    end
    
end

        